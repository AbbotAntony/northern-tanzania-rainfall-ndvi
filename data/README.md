# Input data (not included)

Place these files here before running `R_thesis_analysis.R`.

**From Google Earth Engine** (run `GEE_rainfall_ndvi_composites.js`, then download from the Drive folder `Thesis_GEE_exports`):

- `MAM_Rainfall_2003_2023.tif`, `OND_Rainfall_2003_2023.tif`, `ANN_Rainfall_2003_2023.tif`
- `MAM_NDVI_2003_2023.tif`, `OND_NDVI_2003_2023.tif`, `ANN_NDVI_2003_2023.tif`
- `CHIRPS_monthly_at_TMA_stations.csv`
- `Area_mean_series_2003_2023.csv`

**From the Tanzania Meteorological Authority**:

- Monthly station rainfall for Arusha, Kilimanjaro, Manyara and Tanga, 2003–2023. The script reads the workbook in TMA's layout (one block per station: a title row, a `Year, JAN ... DEC` header row, then one row per year), or a long CSV in the `TMA_template.csv` layout. Set the file name in `TMA_FILE`.

TMA station data are obtained from TMA under its own terms and are not redistributed in this repository. Request them from TMA headquarters, Dar es Salaam.
