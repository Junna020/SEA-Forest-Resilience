#!/usr/bin/env Rscript

# Simplified calculation of SIF-based ecosystem resilience (AR1).
#
# Workflow for each pixel:
#   1. Decode monthly GOSIF values and remove fill values.
#   2. Set negative SIF values to zero and remove very-low-SIF pixels.
#   3. Remove the mean seasonal cycle.
#   4. Remove a linear long-term trend.
#   5. Calculate lag-1 autocorrelation in a moving window.
#
# Usage:
# Rscript R/01_calculate_sif_ar1.R INPUT_DIR OUTPUT_DIR \
#   [WINDOW_MONTHS=60] [CORES=1] [SCALE_FACTOR=0.0001] [MASK_FILE=NONE]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    paste(
      "Usage:",
      "Rscript R/01_calculate_sif_ar1.R INPUT_DIR OUTPUT_DIR",
      "[WINDOW_MONTHS=60] [CORES=1] [SCALE_FACTOR=0.0001] [MASK_FILE=NONE]"
    )
  )
}

input_dir <- normalizePath(args[1], mustWork = TRUE)
output_dir <- normalizePath(args[2], mustWork = FALSE)
window_months <- if (length(args) >= 3L) as.integer(args[3]) else 60L
cores <- if (length(args) >= 4L) as.integer(args[4]) else 1L
scale_factor <- if (length(args) >= 5L) as.numeric(args[5]) else 0.0001
mask_file <- if (length(args) >= 6L && toupper(args[6]) != "NONE") {
  normalizePath(args[6], mustWork = TRUE)
} else {
  NA_character_
}

if (!requireNamespace("terra", quietly = TRUE)) {
  stop("Package 'terra' is required. Install it with install.packages('terra').")
}
if (!is.finite(window_months) || window_months < 12L) stop("WINDOW_MONTHS must be >= 12.")
if (!is.finite(cores) || cores < 1L) stop("CORES must be >= 1.")
if (!is.finite(scale_factor) || scale_factor <= 0) stop("SCALE_FACTOR must be > 0.")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
terra::terraOptions(progress = 1, memfrac = 0.65)

file_pattern <- "^GOSIF_[0-9]{4}\\.M[0-9]{2}\\.tif$"
files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
if (!length(files)) stop("No monthly GOSIF GeoTIFFs were found in: ", input_dir)

file_names <- basename(files)
years <- as.integer(sub("^GOSIF_([0-9]{4})\\.M[0-9]{2}\\.tif$", "\\1", file_names))
months <- as.integer(sub("^GOSIF_[0-9]{4}\\.M([0-9]{2})\\.tif$", "\\1", file_names))
dates <- as.Date(sprintf("%04d-%02d-01", years, months))
ord <- order(dates)
files <- files[ord]
dates <- dates[ord]
years <- years[ord]
months <- months[ord]

if (anyDuplicated(dates)) stop("Duplicate monthly files were found.")
expected_dates <- seq(min(dates), max(dates), by = "month")
if (!identical(dates, expected_dates)) {
  missing_dates <- setdiff(expected_dates, dates)
  stop("The monthly series is incomplete. Missing: ", paste(missing_dates, collapse = ", "))
}
if (length(files) < window_months) stop("The time series is shorter than the moving window.")

message("Reading ", length(files), " monthly layers from ", min(dates), " to ", max(dates), ".")
sif <- terra::rast(files)

if (!is.na(mask_file)) {
  message("Applying analysis mask: ", mask_file)
  analysis_mask <- terra::rast(mask_file)[[1]]
  if (!terra::same.crs(analysis_mask, sif)) {
    analysis_mask <- terra::project(analysis_mask, sif[[1]], method = "near")
  } else if (!isTRUE(terra::compareGeom(analysis_mask, sif[[1]], stopOnError = FALSE))) {
    analysis_mask <- terra::resample(analysis_mask, sif[[1]], method = "near")
  }
  analysis_mask <- terra::ifel(is.na(analysis_mask) | analysis_mask <= 0, NA, 1)
  sif <- terra::mask(sif, analysis_mask)
}

# These defaults follow the GOSIF documentation and the analysis used in this
# study. Change SCALE_FACTOR to 1 only if the input files are already decoded.
fill_values <- c(32766, 32767)
minimum_mean_sif <- 0.01
minimum_valid_fraction <- 0.90
n_output <- length(files) - window_months + 1L

calculate_pixel_ar1 <- function(
    x, month_index, window_size, value_scale, fill_codes,
    mean_sif_threshold, min_valid_fraction) {
  result_na <- rep(NA_real_, length(x) - window_size + 1L)
  if (all(is.na(x))) return(result_na)

  x[x %in% fill_codes] <- NA_real_
  x <- x * value_scale
  x[x < 0] <- 0

  if (sum(is.finite(x)) < ceiling(length(x) * min_valid_fraction)) return(result_na)
  if (mean(x, na.rm = TRUE) < mean_sif_threshold) return(result_na)

  seasonal_mean <- vapply(1:12, function(m) mean(x[month_index == m], na.rm = TRUE), numeric(1))
  if (any(!is.finite(seasonal_mean))) return(result_na)
  deseasonalized <- x - seasonal_mean[month_index]

  time_index <- seq_along(x)
  trend_fit <- stats::lm(deseasonalized ~ time_index, na.action = stats::na.exclude)
  anomalies <- deseasonalized - stats::predict(trend_fit, newdata = data.frame(time_index = time_index))

  ar1 <- rep(NA_real_, length(x) - window_size + 1L)
  for (i in seq_along(ar1)) {
    z <- anomalies[i:(i + window_size - 1L)]
    if (all(is.finite(z)) && stats::sd(z) > 0) {
      ar1[i] <- as.numeric(stats::acf(z, lag.max = 1L, plot = FALSE)$acf[2])
    }
  }
  ar1
}

message("Calculating ", window_months, "-month moving-window AR1 using ", cores, " core(s).")
ar1_stack <- terra::app(
  sif,
  fun = calculate_pixel_ar1,
  month_index = months,
  window_size = window_months,
  value_scale = scale_factor,
  fill_codes = fill_values,
  mean_sif_threshold = minimum_mean_sif,
  min_valid_fraction = minimum_valid_fraction,
  cores = cores
)

endpoint_dates <- dates[window_months:length(dates)]
stopifnot(terra::nlyr(ar1_stack) == n_output)
names(ar1_stack) <- paste0("AR1_", format(endpoint_dates, "%Y.M%m"))
terra::time(ar1_stack) <- endpoint_dates

stack_file <- file.path(output_dir, sprintf("SIF_AR1_sliding_%02d_month.tif", window_months))
terra::writeRaster(
  ar1_stack, stack_file, overwrite = TRUE,
  wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "BIGTIFF=YES"))
)

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
mean_ar1 <- terra::app(ar1_stack, safe_mean)
names(mean_ar1) <- "mean_AR1"
mean_file <- file.path(output_dir, "SIF_AR1_temporal_mean.tif")
terra::writeRaster(
  mean_ar1, mean_file, overwrite = TRUE,
  wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
)

settings <- data.frame(
  setting = c(
    "input_layers", "first_month", "last_month", "window_months", "output_layers",
    "scale_factor", "fill_values", "negative_values", "minimum_mean_sif",
    "minimum_valid_fraction", "mask_file"
  ),
  value = c(
    length(files), as.character(min(dates)), as.character(max(dates)), window_months,
    n_output, scale_factor, paste(fill_values, collapse = ";"), "set to zero",
    minimum_mean_sif, minimum_valid_fraction, ifelse(is.na(mask_file), "none", mask_file)
  )
)
utils::write.csv(settings, file.path(output_dir, "SIF_AR1_processing_settings.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))

message("Completed.")
message("AR1 time series: ", normalizePath(stack_file, winslash = "/"))
message("Temporal mean AR1: ", normalizePath(mean_file, winslash = "/"))
