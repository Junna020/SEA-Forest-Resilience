# SIF Resilience Analysis

This repository contains reproducible R examples for:

1. calculating SIF-based ecosystem resilience using moving-window lag-1 autocorrelation (AR1); and
2. fitting an XGBoost model and interpreting it with SHAP values.

Large raster datasets are not included in this repository.

## Repository structure

```text
R/
  01_calculate_sif_ar1.R
  02_xgboost_shap_example.R
.gitignore
README.md
```

## Requirements

- R >= 4.3
- `terra`
- `data.table`
- `xgboost`
- `ggplot2`
- `svglite`
- `ragg`

```r
install.packages(c("terra", "data.table", "xgboost", "ggplot2", "svglite", "ragg"))
```

## GOSIF data

Monthly GOSIF rasters must be downloaded separately:

- Product page: <https://globalecology.unh.edu/data/GOSIF.html>
- Monthly GOSIF v2 files: <https://data.globalecology.unh.edu/data/GOSIF_v2/Monthly/>
- Documentation and data-use policy: <https://data.globalecology.unh.edu/data/GOSIF_v2/Fair_Data_Use_Policy_and_Readme_GOSIF_v2.pdf>

Expected filenames:

```text
GOSIF_2000.M03.tif
GOSIF_2000.M04.tif
...
GOSIF_2024.M12.tif
```

The default scale factor is `0.0001`. Fill values `32766` and `32767` are treated as missing data. Use the same GOSIF version and study period reported in the manuscript.

Reference:

> Li, X. & Xiao, J. A global, 0.05-degree product of solar-induced chlorophyll fluorescence derived from OCO-2, MODIS, and reanalysis data. *Remote Sensing* **11**, 517 (2019). <https://doi.org/10.3390/rs11050517>

## SIF-based AR1

The script removes the monthly seasonal cycle and linear trend, then calculates AR1 in a 60-month moving window.

```bash
Rscript R/01_calculate_sif_ar1.R INPUT_DIR OUTPUT_DIR 60 1 0.0001 NONE
```

Arguments after the input and output directories are window length, CPU cores, scale factor, and an optional raster mask. Use `NONE` when no mask is required.

Main outputs:

- `SIF_AR1_sliding_60_month.tif`
- `SIF_AR1_temporal_mean.tif`
- `SIF_AR1_processing_settings.csv`
- `R_sessionInfo.txt`

## XGBoost and SHAP

The example creates a simulated pixel-year dataset, splits it by pixel, trains an XGBoost model, evaluates held-out predictions, and calculates exact TreeSHAP values.

```bash
Rscript R/02_xgboost_shap_example.R outputs/xgboost_shap_example
```

The simulated data demonstrate the workflow only and must not be presented as empirical results. Manuscript reproduction should use the processed study data while retaining spatially independent data splitting.

## Upload to GitHub

Create an empty repository at <https://github.com/new> using these settings:

- Repository name: `sif-resilience-analysis`
- Visibility: `Private` during preparation or `Public` after release
- Add README: unchecked
- Add `.gitignore`: none
- License: MIT is recommended

Then run:

```bash
git init
git add README.md .gitignore R
git commit -m "Initial release"
git branch -M main
git remote add origin https://github.com/YOUR_ACCOUNT/sif-resilience-analysis.git
git push -u origin main
```

Do not upload raw rasters, full-resolution outputs, credentials, local package libraries, or restricted third-party data. Large research outputs can be deposited in Zenodo, Dryad, or Figshare.

## Code availability

The R code for calculating moving-window SIF-based AR1 and demonstrating the XGBoost-SHAP workflow is available at `https://github.com/ACCOUNT/REPOSITORY`. Original GOSIF data are available from the University of New Hampshire Global Ecology Data Repository.

Replace the placeholder repository address before publication.
