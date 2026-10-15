/**************************************************************************
 * MSc thesis - Rainfall and rangeland vegetation, Northern Tanzania 2003-2023
 * STEP 1 of 2: Google Earth Engine
 *
 * Improved version of the two original thesis scripts (CHIRPS rainfall and
 * MODIS NDVI). Same datasets, seasons, clipping and export scales, with
 * these changes:
 *
 *  1. ONE script for rainfall and NDVI instead of two.
 *  2. ANNUAL composites added (Jan-Dec). The originals built only MAM and
 *     OND, but thesis Figures 4.5-4.7, 4.10, 4.12-4.13 and Tables 4.2 and 4.5
 *     need annual rainfall and NDVI.
 *  3. Six multi-band exports (one band per year) instead of 126 single-year
 *     tasks, and no getInfo() calls inside loops (slow, and can time out).
 *     Set EXPORT_PER_YEAR = true to also get the original file names
 *     (MAM_Rainfall_2003.tif ...).
 *  4. NEW: monthly CHIRPS rainfall at the four TMA stations, exported as CSV
 *     for validation against the TMA Excel data (month by month, not only
 *     annual totals).
 *  5. NEW: study-area mean series for all six variables, exported as CSV,
 *     so the thesis numbers (CV, Mann-Kendall, regression) can be checked
 *     directly, before any raster work in R.
 *  6. Optional MODIS quality filter (off by default, as in the original).
 *
 * Exports (Drive folder below):
 *   MAM/OND/ANN_Rainfall_2003_2023.tif   21 bands, mm, 5000 m
 *   MAM/OND/ANN_NDVI_2003_2023.tif       21 bands, raw NDVI x 10000, 500 m
 *   CHIRPS_monthly_at_TMA_stations.csv   station, year, month, chirps_mm
 *   Area_mean_series_2003_2023.csv       year + six study-area means
 **************************************************************************/

// ================================================================ SETTINGS
var START_YEAR = 2003;
var END_YEAR   = 2023;

var studyArea = ee.FeatureCollection('projects/ee-abbotantony37/assets/C_boundary');

var RAIN_SCALE    = 5000;   // as in the original rainfall export
var NDVI_SCALE    = 500;    // as in the original NDVI export
var CHIRPS_NATIVE = 5566;   // CHIRPS pixel size (0.05 deg) for station extraction
var EXPORT_CRS    = 'EPSG:4326';
var DRIVE_FOLDER  = 'Thesis_GEE_exports';

var APPLY_NDVI_QA   = false;  // true = keep only good/marginal MODIS pixels
var EXPORT_PER_YEAR = false;  // true = also export the 126 single-year files

// TMA stations (thesis Table 3.1). CHECK these against TMA station metadata:
// the CHIRPS values used for validation are taken at exactly these points.
var stations = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([35.59, -2.29]), {station: 'Arusha'}),
  ee.Feature(ee.Geometry.Point([37.68, -3.55]), {station: 'Kilimanjaro'}),
  ee.Feature(ee.Geometry.Point([39.03, -5.08]), {station: 'Tanga'}),
  ee.Feature(ee.Geometry.Point([36.89, -4.19]), {station: 'Manyara'})
]);

// ============================================================= COLLECTIONS
var startDate = ee.Date.fromYMD(START_YEAR, 1, 1);
var endDate   = ee.Date.fromYMD(END_YEAR + 1, 1, 1);
var years     = ee.List.sequence(START_YEAR, END_YEAR);
var bandNames = years.map(function (y) {
  return ee.String('y').cat(ee.Number(y).format('%d'));
});

var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .filterDate(startDate, endDate)
  .select('precipitation');

var modis = ee.ImageCollection('MODIS/061/MOD13Q1')
  .filterDate(startDate, endDate);
if (APPLY_NDVI_QA) {
  modis = modis.map(function (img) {
    return img.updateMask(img.select('SummaryQA').lte(1));
  });
}
var modisNDVI = modis.select('NDVI');

// ============================================================== COMPOSITES
// Rainfall total (mm) for months m1-m2 of one year - same logic as the
// original getSeasonalRainfall(): daily CHIRPS summed.
function rainTotal(year, m1, m2) {
  year = ee.Number(year);
  return chirps
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(m1, m2, 'month'))
    .sum()
    .rename('b');
}

// Mean NDVI (raw MODIS integers) for months m1-m2 of one year - same logic
// as the original getMAMNDVI() / getONDNDVI(): 16-day composites averaged.
function ndviMean(year, m1, m2) {
  year = ee.Number(year);
  return modisNDVI
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(m1, m2, 'month'))
    .mean()
    .rename('b');
}

// 21-band image (y2003 ... y2023), clipped to the study area.
function yearlyStack(fn, m1, m2) {
  var col = ee.ImageCollection.fromImages(
    years.map(function (y) { return fn(y, m1, m2); })
  );
  return col.toBands().rename(bandNames).clip(studyArea);
}

var stacks = {
  MAM_Rainfall: yearlyStack(rainTotal, 3, 5).toFloat(),
  OND_Rainfall: yearlyStack(rainTotal, 10, 12).toFloat(),
  ANN_Rainfall: yearlyStack(rainTotal, 1, 12).toFloat(),
  MAM_NDVI:     yearlyStack(ndviMean, 3, 5).toInt16(),
  OND_NDVI:     yearlyStack(ndviMean, 10, 12).toInt16(),
  ANN_NDVI:     yearlyStack(ndviMean, 1, 12).toInt16()
};

// ================================================================ PREVIEW
// One layer per variable (the originals added 84 layers via getInfo, which
// is very slow). Change PREVIEW_YEAR to look at another year.
var PREVIEW_YEAR = 'y2018';
Map.centerObject(studyArea, 7);
Map.addLayer(stacks.MAM_Rainfall.select(PREVIEW_YEAR),
  {min: 0, max: 800, palette: ['ffffd9', '7fcdbb', '1d91c0', '081d58']},
  'MAM rainfall ' + PREVIEW_YEAR);
Map.addLayer(stacks.OND_Rainfall.select(PREVIEW_YEAR),
  {min: 0, max: 500, palette: ['ffffd9', '7fcdbb', '1d91c0', '081d58']},
  'OND rainfall ' + PREVIEW_YEAR, false);
Map.addLayer(stacks.MAM_NDVI.select(PREVIEW_YEAR),
  {min: 0, max: 9000, palette: ['d73027', 'fee08b', '1a9850']},
  'MAM NDVI ' + PREVIEW_YEAR + ' (raw)', false);
Map.addLayer(studyArea.style({color: 'black', fillColor: '00000000'}), {}, 'Boundary');
Map.addLayer(stations, {color: 'red'}, 'TMA stations');

// ===================================================== RASTER EXPORTS
function exportImage(img, name, scale) {
  Export.image.toDrive({
    image: img,
    description: name,
    fileNamePrefix: name,
    folder: DRIVE_FOLDER,
    region: studyArea.geometry(),
    scale: scale,
    crs: EXPORT_CRS,
    fileFormat: 'GeoTIFF',
    maxPixels: 1e10
  });
}

exportImage(stacks.MAM_Rainfall, 'MAM_Rainfall_2003_2023', RAIN_SCALE);
exportImage(stacks.OND_Rainfall, 'OND_Rainfall_2003_2023', RAIN_SCALE);
exportImage(stacks.ANN_Rainfall, 'ANN_Rainfall_2003_2023', RAIN_SCALE);
exportImage(stacks.MAM_NDVI,     'MAM_NDVI_2003_2023',     NDVI_SCALE);
exportImage(stacks.OND_NDVI,     'OND_NDVI_2003_2023',     NDVI_SCALE);
exportImage(stacks.ANN_NDVI,     'ANN_NDVI_2003_2023',     NDVI_SCALE);

// Optional: the original one-file-per-year exports, same names as before.
if (EXPORT_PER_YEAR) {
  var yearList = [];
  for (var y = START_YEAR; y <= END_YEAR; y++) yearList.push(y);
  Object.keys(stacks).forEach(function (key) {
    var scale = key.indexOf('Rainfall') > -1 ? RAIN_SCALE : NDVI_SCALE;
    var prefix = key.replace('ANN_', 'Annual_');
    yearList.forEach(function (yr) {
      exportImage(stacks[key].select('y' + yr), prefix + '_' + yr, scale);
    });
  });
}

// ============================== NEW: monthly CHIRPS at TMA stations (CSV)
var nMonths = (END_YEAR - START_YEAR + 1) * 12;
var monthlyRain = ee.ImageCollection.fromImages(
  ee.List.sequence(0, nMonths - 1).map(function (m) {
    var start = startDate.advance(m, 'month');
    return chirps.filterDate(start, start.advance(1, 'month')).sum()
      .set('year', start.get('year'))
      .set('month', start.get('month'));
  })
);

var stationMonthly = monthlyRain.map(function (img) {
  return img.reduceRegions({
    collection: stations,
    reducer: ee.Reducer.first().setOutputs(['chirps_mm']),
    scale: CHIRPS_NATIVE
  }).map(function (f) {
    return f.set('year', img.get('year'), 'month', img.get('month'))
            .setGeometry(null);
  });
}).flatten();

print('Station extraction preview:', stationMonthly.limit(8));

Export.table.toDrive({
  collection: stationMonthly,
  description: 'CHIRPS_monthly_at_TMA_stations',
  folder: DRIVE_FOLDER,
  fileFormat: 'CSV',
  selectors: ['station', 'year', 'month', 'chirps_mm']
});

// ========================= NEW: study-area mean series for checking (CSV)
// Means over ALL pixels in the boundary (no rangeland mask - that is done
// in R). NDVI is rescaled to -1..1 here.
function meanDict(img, scale) {
  return img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: studyArea.geometry(),
    scale: scale,
    maxPixels: 1e13,
    tileScale: 4
  });
}
var means = {
  MAM_rain: meanDict(stacks.MAM_Rainfall, RAIN_SCALE),
  OND_rain: meanDict(stacks.OND_Rainfall, RAIN_SCALE),
  ANN_rain: meanDict(stacks.ANN_Rainfall, RAIN_SCALE),
  MAM_ndvi: meanDict(stacks.MAM_NDVI, NDVI_SCALE),
  OND_ndvi: meanDict(stacks.OND_NDVI, NDVI_SCALE),
  ANN_ndvi: meanDict(stacks.ANN_NDVI, NDVI_SCALE)
};

var seriesTable = ee.FeatureCollection(years.map(function (y) {
  var b = ee.String('y').cat(ee.Number(y).format('%d'));
  return ee.Feature(null, {
    year: y,
    MAM_rain: means.MAM_rain.get(b),
    OND_rain: means.OND_rain.get(b),
    ANN_rain: means.ANN_rain.get(b),
    MAM_ndvi: ee.Number(means.MAM_ndvi.get(b)).multiply(0.0001),
    OND_ndvi: ee.Number(means.OND_ndvi.get(b)).multiply(0.0001),
    ANN_ndvi: ee.Number(means.ANN_ndvi.get(b)).multiply(0.0001)
  });
}));

Export.table.toDrive({
  collection: seriesTable,
  description: 'Area_mean_series_2003_2023',
  folder: DRIVE_FOLDER,
  fileFormat: 'CSV',
  selectors: ['year', 'MAM_rain', 'OND_rain', 'ANN_rain',
              'MAM_ndvi', 'OND_ndvi', 'ANN_ndvi']
});

/**************************************************************************
 * NEXT: open the Tasks tab and run all 8 exports. Download them from the
 * Drive folder into the R project's data/ folder, together with
 * ABBOT_DATA_2003-2023.xlsx, then run R_thesis_analysis.R.
 **************************************************************************/
