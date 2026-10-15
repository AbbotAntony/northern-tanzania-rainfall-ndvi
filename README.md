# Rainfall and rangeland vegetation in Northern Tanzania (2003–2023)

Code for analysing seasonal and annual rainfall (CHIRPS) and vegetation greenness (MODIS NDVI) across Arusha, Kilimanjaro, Manyara and Tanga, Tanzania, over 2003–2023. It covers the analyses in the MSc dissertation of Abbot Josephat Antony (GIS and Remote Sensing, Ardhi University): validation against Tanzania Meteorological Authority (TMA) station records, spatial and temporal variability, Mann–Kendall and Sen's slope trend tests, regression, and pixel-wise correlation between rainfall and NDVI.

## About this code

The analysis pipeline was rebuilt after the dissertation was completed, from the original Google Earth Engine scripts and the methods described in the dissertation, so that the results can be reproduced end to end. It extends the original scripts with annual composites, station-level validation and data-quality checks. Re-running it on the full data is how the reported results are being verified, so small differences from the dissertation tables are possible.

## Files

| File | Purpose |
|---|---|
| `GEE_rainfall_ndvi_composites.js` | Google Earth Engine: seasonal (MAM, OND) and annual rainfall and NDVI stacks, monthly CHIRPS at the TMA stations, and study-area mean series |
| `R_thesis_analysis.R` | R: quality control, validation, maps, coefficient of variation, trends, regression, correlation, pixel-wise correlation |
| `TMA_template.csv` | Long-format layout for station rainfall, if not using the TMA Excel layout |
| `data/` | Input folder (data are not included - see `data/README.md`) |

## Requirements

- A Google Earth Engine account, and the study-area boundary uploaded as an Earth Engine asset (update the asset path at the top of the GEE script).
- R (4.1 or later) with: `terra`, `ggplot2`, `patchwork`, `dplyr`, `tidyr`, `zoo`, `readxl`. Optional: `trend`, `modifiedmk` (the script falls back to built-in Mann–Kendall and Sen's slope code if these are missing).

```r
install.packages(c("terra", "ggplot2", "patchwork", "dplyr", "tidyr", "zoo",
                   "readxl", "trend", "modifiedmk"))
```

## How to run

1. Paste `GEE_rainfall_ndvi_composites.js` into the Earth Engine Code Editor and click **Run**. In the **Tasks** tab, run all eight exports (six GeoTIFFs, two CSVs).
2. Download the exports from the Drive folder `Thesis_GEE_exports` into `data/`.
3. Put the TMA monthly rainfall workbook in `data/` (see `data/README.md`).
4. Open `R_thesis_analysis.R`, check the settings block at the top, and run the whole script.
5. Figures are written to `outputs/figures/`, tables to `outputs/tables/`, and the key results are printed at the end.

## Outputs

| Dissertation output | File |
|---|---|
| Figure 4.1 | `Figure_4_1_CHIRPS_vs_TMA.png` (+ `Figure_4_1b_validation_by_station.png`, `validation_metrics.csv`) |
| Figures 4.2–4.7 | `Figure_4_2_...` to `Figure_4_7_...` (map grids) |
| Section 4.4.3 | `CV_table.csv` |
| Figures 4.8–4.10 | `Figure_4_8_MAM_trend.png`, `Figure_4_9_OND_trend.png`, `Figure_4_10_annual_trend.png` |
| Tables 4.1–4.2 | `Table_4_1_seasonal_trends.csv`, `Table_4_2_annual_trends.csv` |
| Tables 4.3–4.5 | `Table_4_3_regression_MAM.csv`, `Table_4_4_regression_OND.csv`, `Table_4_5_regression_annual.csv` |
| Figures 4.11–4.12 | `Figure_4_11_seasonal_correlation.png`, `Figure_4_12_annual_correlation.png` |
| Figure 4.13 | `Figure_4_13_pixelwise_correlation.png` (+ `pixelwise_correlation.tif`) |

## Settings worth checking

- `RANGELAND_MASK`, `MASK_METHOD`, `NDVI_RANGE` (R): the 0.25–0.55 NDVI rangeland definition.
- `PIXEL_GRID` (R): grid used for the pixel-wise correlation.
- `STATIONS` (R) and `stations` (GEE): TMA station coordinates - confirm against TMA station metadata.
- `TMA_FILE`, `STATION_NAME_FIX` (R): input workbook name and station-name spelling fixes.
- `APPLY_NDVI_QA` (GEE): MODIS quality filter, off by default to match the original exports.

## Data sources

- Funk, C. et al. (2015). The climate hazards infrared precipitation with stations - a new environmental record for monitoring extremes. *Scientific Data*, 2, 150066. CHIRPS Daily via Earth Engine: `UCSB-CHG/CHIRPS/DAILY`.
- Didan, K. (2021). MODIS/Terra Vegetation Indices 16-Day L3 Global 250m SIN Grid V061. NASA EOSDIS LP DAAC. Earth Engine: `MODIS/061/MOD13Q1`.
- Tanzania Meteorological Authority (TMA): monthly station rainfall for Arusha, Kilimanjaro, Manyara and Tanga (not redistributed here).

## Citation

If you use this code, please cite the dissertation:

Antony, A. J. (2025). *Spatio-temporal analysis of rainfall trends in Northern Tanzania's rangelands*. MSc dissertation, Department of Geospatial Sciences and Technology, Ardhi University, Dar es Salaam.

## Licence

Code released under the MIT Licence - see `LICENSE`.
