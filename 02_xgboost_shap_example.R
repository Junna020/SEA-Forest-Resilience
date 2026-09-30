#!/usr/bin/env Rscript

# Reproducible XGBoost + SHAP example using simulated pixel-year data.
# The example uses a pixel-level split, so observations from one pixel never
# appear in more than one of the training, validation, and test sets.

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[1] else "outputs/xgboost_shap_example"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required <- c("data.table", "xgboost", "ggplot2", "svglite", "ragg")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install required packages first: ", paste(missing, collapse = ", "))
}

set.seed(20260930)
n_pixels <- 1000L
years <- 2005:2020

# -----------------------------------------------------------------------------
# 1. Simulate an example environmental panel
# -----------------------------------------------------------------------------
pixel <- data.table::data.table(
  pixel_id = seq_len(n_pixels),
  MAT_base = stats::rnorm(n_pixels, 24, 2.2),
  MAP_base = stats::rgamma(n_pixels, shape = 8, scale = 280),
  HFP_base = pmax(0, stats::rgamma(n_pixels, shape = 1.8, scale = 3)),
  SM_base = pmin(0.42, pmax(0.12, stats::rnorm(n_pixels, 0.27, 0.045)))
)
dat <- data.table::CJ(pixel_id = pixel$pixel_id, year = years)
dat <- pixel[dat, on = "pixel_id"]
dat[, year_index := year - min(year)]
dat[, `:=`(
  MAT = MAT_base + 0.025 * year_index + rnorm(.N, 0, 0.45),
  MAP = pmax(400, MAP_base + rnorm(.N, 0, 220)),
  HFP = pmax(0, HFP_base + 0.08 * year_index + rnorm(.N, 0, 0.35)),
  SM = pmin(0.48, pmax(0.08, SM_base + rnorm(.N, 0, 0.018)))
)]
dat[, VPD := pmax(0.08, 0.16 * MAT - 4.5 * SM - 0.00008 * MAP + rnorm(.N, 0, 0.12))]
dat[, SolarRad := 185 + 3.2 * MAT - 0.006 * MAP + rnorm(.N, 0, 8)]
dat[, MAT_variability := abs(rnorm(.N, 1.1 + 0.015 * year_index, 0.22))]
dat[, VPD_variability := abs(rnorm(.N, 0.16 + 0.018 * VPD, 0.035))]

z <- function(x) as.numeric(scale(x))
dat[, AR1 :=
  0.34 +
  0.050 * z(VPD) +
  0.040 * z(HFP) -
  0.035 * z(MAP) -
  0.025 * z(SM) +
  0.030 * z(VPD) * z(HFP) +
  0.018 * pmax(z(MAT), 0)^2 +
  0.015 * z(VPD_variability) +
  rnorm(.N, 0, 0.055)
]
dat[, AR1 := pmin(0.90, pmax(-0.20, AR1))]

feature_names <- c(
  "MAT", "MAP", "VPD", "SM", "HFP", "SolarRad",
  "MAT_variability", "VPD_variability"
)

# -----------------------------------------------------------------------------
# 2. Split by pixel: 70% train, 15% validation, 15% independent test
# -----------------------------------------------------------------------------
pixel_order <- sample(pixel$pixel_id)
n_train <- floor(0.70 * n_pixels)
n_validation <- floor(0.15 * n_pixels)
train_pixels <- pixel_order[seq_len(n_train)]
validation_pixels <- pixel_order[n_train + seq_len(n_validation)]
test_pixels <- setdiff(pixel_order, c(train_pixels, validation_pixels))

train_idx <- dat$pixel_id %in% train_pixels
validation_idx <- dat$pixel_id %in% validation_pixels
test_idx <- dat$pixel_id %in% test_pixels

x_train <- as.matrix(dat[train_idx, ..feature_names])
x_validation <- as.matrix(dat[validation_idx, ..feature_names])
x_test <- as.matrix(dat[test_idx, ..feature_names])
y_train <- dat[train_idx, AR1]
y_validation <- dat[validation_idx, AR1]
y_test <- dat[test_idx, AR1]

dtrain <- xgboost::xgb.DMatrix(x_train, label = y_train, feature_names = feature_names)
dvalidation <- xgboost::xgb.DMatrix(x_validation, label = y_validation, feature_names = feature_names)
dtest <- xgboost::xgb.DMatrix(x_test, label = y_test, feature_names = feature_names)

# -----------------------------------------------------------------------------
# 3. Fit XGBoost with early stopping
# -----------------------------------------------------------------------------
params <- list(
  objective = "reg:squarederror",
  eval_metric = "rmse",
  eta = 0.03,
  max_depth = 5L,
  min_child_weight = 5,
  subsample = 0.80,
  colsample_bytree = 0.80,
  lambda = 1,
  alpha = 0
)

model <- xgboost::xgb.train(
  params = params,
  data = dtrain,
  nrounds = 1500L,
  evals = list(train = dtrain, validation = dvalidation),
  early_stopping_rounds = 40L,
  print_every_n = 50L,
  verbose = 1L
)

prediction <- predict(model, dtest)
rmse <- sqrt(mean((y_test - prediction)^2))
r2 <- 1 - sum((y_test - prediction)^2) / sum((y_test - mean(y_test))^2)
best_iteration <- suppressWarnings(as.integer(xgboost::xgb.attr(model, "best_iteration")))
if (!length(best_iteration) || !is.finite(best_iteration)) best_iteration <- NA_integer_
metrics <- data.frame(
  split = "held-out pixels",
  n_pixels = length(test_pixels),
  n_observations = length(y_test),
  best_iteration = best_iteration,
  RMSE = rmse,
  R_squared = r2
)
utils::write.csv(metrics, file.path(output_dir, "model_performance.csv"), row.names = FALSE)

test_predictions <- dat[test_idx, .(pixel_id, year, observed_AR1 = AR1)]
test_predictions[, predicted_AR1 := prediction]
data.table::fwrite(test_predictions, file.path(output_dir, "held_out_predictions.csv"))
xgboost::xgb.save(model, file.path(output_dir, "xgboost_model.ubj"))

# -----------------------------------------------------------------------------
# 4. Calculate exact TreeSHAP values for held-out observations
# -----------------------------------------------------------------------------
shap_with_bias <- predict(model, dtest, predcontrib = TRUE, approxcontrib = FALSE)
shap_values <- shap_with_bias[, feature_names, drop = FALSE]
importance <- data.frame(
  feature = feature_names,
  mean_abs_SHAP = colMeans(abs(shap_values))
)
importance <- importance[order(importance$mean_abs_SHAP, decreasing = TRUE), ]
utils::write.csv(importance, file.path(output_dir, "shap_importance.csv"), row.names = FALSE)

top_features <- head(importance$feature, 8L)
plot_rows <- sample(seq_len(nrow(x_test)), min(4000L, nrow(x_test)))
shap_long <- data.table::rbindlist(lapply(top_features, function(v) {
  raw_value <- x_test[plot_rows, v]
  value_01 <- if (diff(range(raw_value)) == 0) rep(0.5, length(raw_value)) else {
    (raw_value - min(raw_value)) / diff(range(raw_value))
  }
  data.table::data.table(
    feature = v,
    SHAP = shap_values[plot_rows, v],
    feature_value = raw_value,
    feature_value_01 = value_01
  )
}))
shap_long[, feature := factor(feature, levels = rev(top_features))]

theme_pub <- ggplot2::theme_classic(base_size = 8, base_family = "Arial") +
  ggplot2::theme(
    axis.line = ggplot2::element_line(linewidth = 0.35, colour = "#222222"),
    axis.ticks = ggplot2::element_line(linewidth = 0.35),
    panel.grid.major.y = ggplot2::element_line(colour = "#ECECEC", linewidth = 0.25),
    legend.position = "right"
  )

p_beeswarm <- ggplot2::ggplot(
  shap_long,
  ggplot2::aes(x = SHAP, y = feature, colour = feature_value_01)
) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "#888888") +
  ggplot2::geom_jitter(height = 0.18, width = 0, size = 0.55, alpha = 0.45) +
  ggplot2::scale_colour_gradientn(
    colours = c("#2166AC", "#F7F7F7", "#B2182B"),
    limits = c(0, 1), breaks = c(0, 1), labels = c("Low", "High"),
    name = "Feature value"
  ) +
  ggplot2::labs(x = "SHAP value (effect on predicted AR1)", y = NULL) +
  theme_pub

importance_plot <- importance[seq_len(min(8L, nrow(importance))), ]
importance_plot$feature <- factor(importance_plot$feature, levels = rev(importance_plot$feature))
p_importance <- ggplot2::ggplot(
  importance_plot,
  ggplot2::aes(mean_abs_SHAP, feature)
) +
  ggplot2::geom_col(width = 0.68, fill = "#4C8DAE") +
  ggplot2::labs(x = "Mean |SHAP|", y = NULL) +
  theme_pub

dependence_features <- head(importance$feature, 4L)
dependence_data <- data.table::rbindlist(lapply(dependence_features, function(v) {
  data.table::data.table(
    feature = v,
    feature_value = x_test[, v],
    SHAP = shap_values[, v]
  )
}))
dependence_data[, feature := factor(feature, levels = dependence_features)]
p_dependence <- ggplot2::ggplot(
  dependence_data,
  ggplot2::aes(feature_value, SHAP)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "#999999", linewidth = 0.3) +
  ggplot2::geom_point(size = 0.45, alpha = 0.22, colour = "#326D88") +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                       colour = "#C84C36", fill = "#EFB9AA", linewidth = 0.7) +
  ggplot2::facet_wrap(~feature, scales = "free_x", ncol = 2) +
  ggplot2::labs(x = "Feature value", y = "SHAP value") +
  theme_pub

save_plot <- function(plot, stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  print(plot)
  grDevices::dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Arial")
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in,
                 units = "in", res = 600, compression = "lzw")
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(paste0(stem, "_preview.png"), width = width_in, height = height_in,
                units = "in", res = 220, background = "white")
  print(plot)
  grDevices::dev.off()
}

save_plot(p_beeswarm, file.path(output_dir, "SHAP_beeswarm"), 120, 92)
save_plot(p_importance, file.path(output_dir, "SHAP_importance"), 90, 80)
save_plot(p_dependence, file.path(output_dir, "SHAP_dependence"), 150, 120)

data.table::fwrite(dat, file.path(output_dir, "simulated_pixel_year_data.csv.gz"), compress = "gzip")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))

message(sprintf("Held-out performance: R2 = %.3f; RMSE = %.3f", r2, rmse))
message("Outputs written to: ", normalizePath(output_dir, winslash = "/"))
