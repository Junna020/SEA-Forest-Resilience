# #Resilience Analysis

This repository provides compact R workflows for calculating SIF-based ecosystem resilience and demonstrating XGBoost-SHAP analysis. The included datasets are synthetic and are intended only for software testing.

## Contents

```text
R/
  01_calculate_sif_ar1.R
  02_xgboost_shap_example.R
data/
  sif_demo/                 72 synthetic monthly GeoTIFFs
  xgboost_demo_data.csv     synthetic pixel-year table
expected_output/
  sif_ar1_processing_settings.csv
  xgboost_model_performance.csv
  xgboost_shap_importance.csv
.gitignore
README.md
```

## System requirements

Tested on Windows 11 x64 with R 4.6.1. No GPU or non-standard hardware is required. A standard desktop with at least 8 GB RAM is recommended. Full-resolution GOSIF processing requires additional memory and temporary disk space.

Tested package versions:

| Package    | Version |
| ---------- | -------:|
| terra      | 1.9-34  |
| data.table | 1.18.4  |
| xgboost    | 3.2.1.1 |
| ggplot2    | 4.0.3   |
| svglite    | 2.2.2   |
| ragg       | 1.5.2   |

## Installation

Install R from <https://cran.r-project.org/> and then run:

```r
install.packages(c("terra", "data.table", "xgboost", "ggplot2", "svglite", "ragg"))
```

Installation typically takes 5–15 minutes on a standard desktop with precompiled packages and a stable internet connection. Compilation from source may take longer.

## Demo 1: SIF-based AR1

The workflow removes the monthly seasonal cycle and a linear temporal trend, then calculates lag-1 autocorrelation in a 60-month moving window.

From the repository root, run:

```bash
Rscript R/01_calculate_sif_ar1.R data/sif_demo outputs/sif_demo 60 1 0.0001 NONE
```

The included 72-layer dataset is fully synthetic, follows the expected GOSIF filename convention, and is not used in the manuscript. The demo takes approximately 4 seconds on the tested Windows desktop.

Expected files:

- `SIF_AR1_sliding_60_month.tif`
- `SIF_AR1_temporal_mean.tif`
- `SIF_AR1_processing_settings.csv`
- `R_sessionInfo.txt`

The settings file should match `expected_output/sif_ar1_processing_settings.csv`. Minor differences in compressed raster file size are harmless.

### Running on GOSIF data

Download monthly GOSIF rasters separately:

- Product page: <https://globalecology.unh.edu/data/GOSIF.html>
- Monthly GOSIF v2 files: <https://data.globalecology.unh.edu/data/GOSIF_v2/Monthly/>
- Documentation and data-use policy: <https://data.globalecology.unh.edu/data/GOSIF_v2/Fair_Data_Use_Policy_and_Readme_GOSIF_v2.pdf>

The script expects sequential files named `GOSIF_YYYY.MMM.tif`. The default scale factor is `0.0001`; fill values `32766` and `32767` are converted to missing values.

```bash
Rscript R/01_calculate_sif_ar1.R INPUT_DIRECTORY OUTPUT_DIRECTORY 60 4 0.0001 NONE
```

The final argument may be replaced with the path to a positive-valued raster mask. Use the same GOSIF version, temporal range, mask, and preprocessing settings reported in the manuscript.

Reference: Li, X. & Xiao, J. *Remote Sensing* **11**, 517 (2019). <https://doi.org/10.3390/rs11050517>

## Demo 2: XGBoost and SHAP

The workflow reads a pixel-year table, splits observations by pixel to limit spatial leakage, trains an XGBoost model with early stopping, evaluates held-out pixels, and calculates exact TreeSHAP values.

```bash
Rscript R/02_xgboost_shap_example.R data/xgboost_demo_data.csv outputs/xgboost_demo
```

The demo takes approximately 7 seconds on the tested Windows desktop. With the supplied dataset, expected held-out performance is approximately RMSE = 0.060 and R² = 0.773. Small numerical differences between platforms are possible.

Expected files include:

- `model_performance.csv`
- `held_out_predictions.csv`
- `xgboost_model.ubj`
- `shap_importance.csv`
- SHAP beeswarm, importance, and dependence plots in PDF, SVG, TIFF, and PNG formats
- `R_sessionInfo.txt`

Reference results are provided in `expected_output/`.

### Running on another dataset

Supply a CSV with one row per pixel-year observation and these columns:

```text
pixel_id, year, AR1, MAT, MAP, VPD, SM, HFP,
SolarRad, MAT_variability, VPD_variability
```

`pixel_id` identifies the spatial unit, `year` identifies time, `AR1` is the response, and all remaining columns are numeric predictors. At least 20 unique pixels are required. Edit `feature_names` in the script if a different predictor set is used.

The supplied synthetic data demonstrate the software only and must not be interpreted as empirical evidence or as a reproduction of the manuscript results.


