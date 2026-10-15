###########################################################################
# MSc thesis - Spatio-temporal analysis of rainfall and rangeland vegetation
# in Northern Tanzania, 2003-2023
#
# STEP 2 of 2: R analysis
#
# Reproduces every Chapter 4 output:
#   Sec 4.2  data quality control            -> tables/QC_report.csv
#   Sec 4.3  CHIRPS vs TMA validation         -> Figure 4.1, r
#   Sec 4.4  seasonal & annual variability    -> Figures 4.2-4.7, CV table
#   Sec 4.5  trends (graphical, MK, Sen)      -> Figures 4.8-4.10, Tables 4.1-4.2
#   Sec 4.6  regression & correlation         -> Tables 4.3-4.5, Figures 4.11-4.13
#
# Inputs (put them in DATA_DIR):
#   six GeoTIFF stacks from GEE_rainfall_ndvi_composites.js
#     (or the older per-year files, e.g. MAM_NDVI_2003.tif - both work)
#   ABBOT_DATA_2003-2023.xlsx           TMA monthly rainfall, as received
#     (a long CSV in the TMA_template.csv layout also works)
#   CHIRPS_monthly_at_TMA_stations.csv  from the GEE script (optional but
#     recommended - gives exact monthly CHIRPS at the stations)
#   Area_mean_series_2003_2023.csv      from the GEE script (optional
#     cross-check of the R area means)
#
# Packages: terra, ggplot2, patchwork, dplyr, tidyr, zoo, readxl
#           trend and modifiedmk are used if installed; otherwise the script
#           falls back to its own Mann-Kendall / Sen's slope code (identical
#           formulas, Equations 1-2 of the thesis).
###########################################################################

# ------------------------------------------------------------------ 0. setup
# install.packages(c("terra","ggplot2","patchwork","dplyr","tidyr","zoo",
#                    "readxl","trend","modifiedmk"))
suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(zoo)
  library(readxl)
})

DATA_DIR <- "data"       # GeoTIFFs, TMA Excel and GEE CSVs live here
OUT_DIR  <- "outputs"    # figures/ and tables/ are created here
YEARS    <- 2003:2023

# Rangeland definition from thesis Sec 3.2.2: pixels with moderate NDVI
# (0.25-0.55); bare / sparse land below 0.20 excluded.
# "static" builds ONE mask from each pixel's long-term mean NDVI and applies
# it to every year, so the same pixels are followed through time (needed for
# valid trend tests). "per_layer" masks each year separately.
# Set RANGELAND_MASK <- FALSE to analyse every vegetated pixel.
RANGELAND_MASK <- TRUE
MASK_METHOD    <- "static"   # "static" or "per_layer"
NDVI_RANGE     <- c(0.25, 0.55)

# Pixel-wise correlation grid (thesis Sec 3.3.5.2 aggregates NDVI to CHIRPS).
# "chirps": NDVI aggregated to the ~5 km rainfall grid (as described in text)
# "ndvi"  : rainfall resampled to the 500 m NDVI grid (finer-looking map)
PIXEL_GRID <- "chirps"

# TMA stations, thesis Table 3.1 (south latitudes entered as negative).
STATIONS <- data.frame(
  station = c("Arusha", "Kilimanjaro", "Tanga", "Manyara"),
  lat     = c(-2.29, -3.55, -5.08, -4.19),
  lon     = c(35.59, 37.68, 39.03, 36.89)
)

# TMA observations and GEE station extraction
TMA_FILE         <- "ABBOT_DATA_2003-2023.xlsx"            # or a long CSV
CHIRPS_STN_FILE  <- "CHIRPS_monthly_at_TMA_stations.csv"   # from GEE
AREA_MEAN_FILE   <- "Area_mean_series_2003_2023.csv"       # from GEE

# Spelling fixes for station headers in the Excel file.
STATION_NAME_FIX <- c(Kilmanjaro = "Kilimanjaro")

FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------- 1. read rasters
# Reads a 21-band stack if present, otherwise per-year files.
read_series <- function(season, var) {
  stack_file <- file.path(DATA_DIR, sprintf("%s_%s_2003_2023.tif", season, var))
  if (file.exists(stack_file)) {
    r <- rast(stack_file)
  } else {
    prefixes <- if (season == "ANN") c("ANN", "Annual", "annual") else season
    pat   <- sprintf("^(%s)_%s_(\\d{4})\\.tif$", paste(prefixes, collapse = "|"), var)
    files <- list.files(DATA_DIR, pattern = pat, full.names = TRUE)
    if (length(files) == 0)
      stop("No input found for ", season, " ", var, " in ", DATA_DIR)
    yrs   <- as.integer(sub(pat, "\\2", basename(files)))
    files <- files[order(yrs)]
    r     <- rast(files)
  }
  if (nlyr(r) != length(YEARS))
    stop(season, " ", var, ": expected ", length(YEARS), " layers, found ", nlyr(r))
  names(r) <- paste0("y", YEARS)
  r
}

rain <- list(MAM = read_series("MAM", "Rainfall"),
             OND = read_series("OND", "Rainfall"),
             ANN = read_series("ANN", "Rainfall"))

ndvi_raw <- list(MAM = read_series("MAM", "NDVI"),
                 OND = read_series("OND", "NDVI"),
                 ANN = read_series("ANN", "NDVI"))

# MODIS NDVI is stored x 10000; rescale to -1..1 when needed.
rescale_ndvi <- function(r) {
  mx <- max(global(r, "max", na.rm = TRUE)$max, na.rm = TRUE)
  if (mx > 1.5) r * 1e-4 else r
}
ndvi_raw <- lapply(ndvi_raw, rescale_ndvi)

# ------------------------------------------------------ 2. rangeland mask
apply_rangeland_mask <- function(r) {
  if (!RANGELAND_MASK) return(r)
  if (MASK_METHOD == "static") {
    m <- mean(ndvi_raw$ANN, na.rm = TRUE)          # long-term mean annual NDVI
    keep <- (m >= NDVI_RANGE[1]) & (m <= NDVI_RANGE[2])
    keep <- resample(keep, r, method = "near")
    mask(r, keep, maskvalues = c(0, NA))
  } else {
    r[(r < NDVI_RANGE[1]) | (r > NDVI_RANGE[2])] <- NA
    r
  }
}
ndvi <- lapply(ndvi_raw, apply_rangeland_mask)

# -------------------------------------------------- 3. area-mean series
area_mean <- function(r) global(r, "mean", na.rm = TRUE)$mean

series <- data.frame(
  year     = YEARS,
  idx      = seq_along(YEARS),          # x used in the trend-line equations
  MAM_rain = area_mean(rain$MAM),
  OND_rain = area_mean(rain$OND),
  ANN_rain = area_mean(rain$ANN),
  MAM_ndvi = area_mean(ndvi$MAM),
  OND_ndvi = area_mean(ndvi$OND),
  ANN_ndvi = area_mean(ndvi$ANN)
)
write.csv(series, file.path(TAB_DIR, "area_mean_series.csv"), row.names = FALSE)

# -------------------------------------------- 4. data quality control (4.2)
qc_series <- function(x, name) {
  q  <- quantile(x, c(.25, .75), na.rm = TRUE)
  iq <- diff(q)
  data.frame(
    variable     = name,
    n            = length(x),
    missing      = sum(is.na(x)),
    mean         = mean(x, na.rm = TRUE),
    sd           = sd(x, na.rm = TRUE),
    min          = min(x, na.rm = TRUE),
    max          = max(x, na.rm = TRUE),
    outliers_3IQR = sum(x < q[1] - 3 * iq | x > q[2] + 3 * iq, na.rm = TRUE)
  )
}

qc <- bind_rows(lapply(names(series)[3:8], function(v) qc_series(series[[v]], v)))

# ---- TMA station data
# Reads the Excel exactly as received from TMA: one sheet with a block per
# station - a title row ("ARUSHA MONTHLY TOTAL RAINFALL (mm)"), a header row
# (Year, JAN ... DEC), then one row per year. A long CSV with columns
# station, year, month, rainfall_mm is also accepted.
read_tma <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    d <- read.csv(path, stringsAsFactors = FALSE)
  } else {
    raw <- suppressMessages(as.data.frame(
      read_excel(path, sheet = 1, col_names = FALSE, col_types = "text")))
    rows <- list(); station <- NA
    for (i in seq_len(nrow(raw))) {
      a <- raw[i, 1]
      if (!is.na(a) && grepl("RAINFALL", toupper(a))) {
        nm <- tolower(strsplit(trimws(a), "\\s+")[[1]][1])
        station <- paste0(toupper(substring(nm, 1, 1)), substring(nm, 2))
        next
      }
      yr <- suppressWarnings(as.integer(a))
      if (!is.na(yr) && !is.na(station)) {
        vals <- suppressWarnings(as.numeric(unlist(raw[i, 2:13])))
        rows[[length(rows) + 1]] <- data.frame(station = station, year = yr,
                                               month = 1:12, rainfall_mm = vals)
      }
    }
    d <- bind_rows(rows)
  }
  fix <- d$station %in% names(STATION_NAME_FIX)
  d$station[fix] <- STATION_NAME_FIX[d$station[fix]]
  d
}

tma_path <- file.path(DATA_DIR, TMA_FILE)
have_tma <- file.exists(tma_path)

if (have_tma) {
  tma <- read_tma(tma_path) %>%
    filter(year %in% YEARS) %>%
    complete(station, year = YEARS, month = 1:12) %>%
    arrange(station, year, month) %>%
    group_by(station) %>%
    mutate(missing_before_fill = is.na(rainfall_mm),
           rainfall_mm = na.approx(rainfall_mm, maxgap = 2, na.rm = FALSE)) %>%
    ungroup()

  # Rain gauges are read to 0.1 mm. Values with many more decimals usually
  # come from gridded / satellite products, not a gauge - flag them, because
  # validating CHIRPS against satellite-derived "observations" inflates r.
  tma_qc <- tma %>%
    group_by(station) %>%
    summarise(
      n_months             = n(),
      months_missing       = sum(missing_before_fill),
      months_still_missing = sum(is.na(rainfall_mm)),
      share_not_0.1mm      = round(mean(abs(rainfall_mm - round(rainfall_mm, 1)) > 1e-6,
                                        na.rm = TRUE), 3),
      mean_annual_mm       = round(sum(rainfall_mm, na.rm = TRUE) / length(YEARS), 1),
      negative_values      = sum(rainfall_mm < 0, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(gauge_precision = ifelse(share_not_0.1mm > 0.10, "SUSPECT", "OK"),
           variable = paste0("TMA_", station))
  qc <- bind_rows(qc, tma_qc)

  suspect_stations <- tma_qc$station[tma_qc$gauge_precision == "SUSPECT"]
  if (length(suspect_stations) > 0)
    warning("TMA series not recorded at gauge precision (0.1 mm): ",
            paste(suspect_stations, collapse = ", "),
            ". Confirm with TMA that these are gauge observations.")

  tma_annual <- tma %>%
    group_by(station, year) %>%
    summarise(tma_mm = if (any(is.na(rainfall_mm))) NA_real_ else sum(rainfall_mm),
              .groups = "drop")
} else {
  suspect_stations <- character(0)
}
write.csv(qc, file.path(TAB_DIR, "QC_report.csv"), row.names = FALSE)

# ------------------------------------------- 5. validation, Figure 4.1 (4.3)
# CHIRPS at the stations: the GEE monthly CSV if present (exact monthly
# totals from daily CHIRPS); otherwise annual values read from the annual
# rainfall raster.
val_metrics <- function(obs, est) {
  ok <- !is.na(obs) & !is.na(est)
  data.frame(n = sum(ok),
             r = cor(obs[ok], est[ok]),
             bias_mm = mean(est[ok] - obs[ok]),
             pbias_pct = 100 * sum(est[ok] - obs[ok]) / sum(obs[ok]),
             rmse_mm = sqrt(mean((est[ok] - obs[ok])^2)))
}

r_val <- NA
if (have_tma) {
  stn_file <- file.path(DATA_DIR, CHIRPS_STN_FILE)
  if (file.exists(stn_file)) {
    chirps_m <- read.csv(stn_file, stringsAsFactors = FALSE) %>%
      mutate(year = as.integer(year), month = as.integer(month))
    val_month <- inner_join(tma, chirps_m, by = c("station", "year", "month")) %>%
      filter(!is.na(rainfall_mm), !is.na(chirps_mm))
    chirps_annual <- chirps_m %>% group_by(station, year) %>%
      summarise(chirps_mm = sum(chirps_mm), .groups = "drop")
    source_note <- "GEE monthly extraction"
  } else {
    pts <- vect(STATIONS, geom = c("lon", "lat"), crs = "EPSG:4326")
    pts <- project(pts, crs(rain$ANN))
    ex  <- terra::extract(rain$ANN, pts)
    chirps_annual <- data.frame(station = STATIONS$station, ex[, -1]) %>%
      pivot_longer(-station, names_to = "year", values_to = "chirps_mm") %>%
      mutate(year = as.integer(sub("y", "", year)))
    val_month <- NULL
    source_note <- "annual raster extraction"
  }

  val <- inner_join(tma_annual, chirps_annual, by = c("station", "year")) %>%
    filter(!is.na(tma_mm), !is.na(chirps_mm))
  if (nrow(val) == 0) stop("No matching station-years: check station names in ",
                           TMA_FILE, " against STATIONS / the GEE CSV.")
  write.csv(val, file.path(TAB_DIR, "validation_station_years.csv"), row.names = FALSE)

  # thesis statistic: pooled annual r across all stations
  r_val <- cor(val$tma_mm, val$chirps_mm)

  val_tab <- bind_rows(
    cbind(scale = "annual", station = "All (thesis)", val_metrics(val$tma_mm, val$chirps_mm)),
    cbind(scale = "annual", station = "All except suspect",
          val_metrics(val$tma_mm[!val$station %in% suspect_stations],
                      val$chirps_mm[!val$station %in% suspect_stations])),
    bind_rows(lapply(split(val, val$station), function(d)
      cbind(scale = "annual", station = d$station[1], val_metrics(d$tma_mm, d$chirps_mm))))
  )
  if (!is.null(val_month)) {
    val_tab <- bind_rows(val_tab,
      cbind(scale = "monthly", station = "All", val_metrics(val_month$rainfall_mm, val_month$chirps_mm)),
      bind_rows(lapply(split(val_month, val_month$station), function(d)
        cbind(scale = "monthly", station = d$station[1],
              val_metrics(d$rainfall_mm, d$chirps_mm)))))
  }
  val_tab <- val_tab %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
  write.csv(val_tab, file.path(TAB_DIR, "validation_metrics.csv"), row.names = FALSE)

  # Figure 4.1 - pooled annual scatter (as in the thesis)
  fig41 <- ggplot(val, aes(tma_mm, chirps_mm)) +
    geom_point(colour = "blue", alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, colour = "red", formula = y ~ x) +
    labs(x = "TMA Rainfall (mm)", y = "CHIRPS Rainfall (mm)") +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, colour = "grey30"))
  ggsave(file.path(FIG_DIR, "Figure_4_1_CHIRPS_vs_TMA.png"), fig41,
         width = 7, height = 4.5, dpi = 300)

  # Supplementary - one panel per station, with the 1:1 line
  fig41b <- ggplot(val, aes(tma_mm, chirps_mm)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(colour = "blue", alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, colour = "red", formula = y ~ x) +
    facet_wrap(~ station, scales = "free") +
    labs(x = "TMA annual rainfall (mm)", y = "CHIRPS annual rainfall (mm)") +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, colour = "grey30"))
  ggsave(file.path(FIG_DIR, "Figure_4_1b_validation_by_station.png"), fig41b,
         width = 8, height = 6, dpi = 300)

  message(sprintf("Validation (%s): pooled annual r = %.4f (n = %d station-years)",
                  source_note, r_val, nrow(val)))
} else {
  message(TMA_FILE, " not found in ", DATA_DIR, " - skipping validation (Figure 4.1).")
}

# Optional cross-check of the R area means against the GEE CSV. Only
# comparable when RANGELAND_MASK is FALSE (GEE averages all pixels).
am_file <- file.path(DATA_DIR, AREA_MEAN_FILE)
if (file.exists(am_file)) {
  gee_series <- read.csv(am_file)
  cmp <- merge(series, gee_series, by = "year", suffixes = c("_R", "_GEE"))
  write.csv(cmp, file.path(TAB_DIR, "area_mean_R_vs_GEE.csv"), row.names = FALSE)
}

# ------------------------------------- 6. variability maps, Figures 4.2-4.7
map_panel <- function(r, title, legend, palette) {
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(df)[3] <- "value"
  ggplot(df, aes(x, y, fill = value)) +
    geom_raster() +
    scale_fill_distiller(palette = palette, direction = 1,
                         name = legend, na.value = "grey50") +
    coord_equal(expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 7) +
    theme(axis.text = element_blank(), panel.grid = element_blank(),
          legend.key.width = unit(0.25, "cm"),
          legend.key.height = unit(0.45, "cm"))
}
rain_panel <- function(r, t) map_panel(r, t, "Rainfall (mm)", "YlGnBu")
ndvi_panel <- function(r, t) map_panel(r, t, "NDVI", "RdYlGn")

# Maps show all vegetated pixels, as in the thesis figures.
save_grid <- function(panels, file, ncol) {
  nrow <- ceiling(length(panels) / ncol)
  g <- wrap_plots(panels, ncol = ncol)
  ggsave(file.path(FIG_DIR, file), g,
         width = 2.2 * ncol, height = 1.9 * nrow, dpi = 300, limitsize = FALSE)
}

panels_42 <- panels_43 <- panels_44 <- panels_45 <- panels_46 <- panels_47 <- list()
for (i in seq_along(YEARS)) {
  y <- YEARS[i]
  mr <- rain_panel(rain$MAM[[i]], paste("MAM Rainfall", y))
  or <- rain_panel(rain$OND[[i]], paste("OND Rainfall", y))
  mn <- ndvi_panel(ndvi_raw$MAM[[i]], paste("MAM NDVI", y))
  on <- ndvi_panel(ndvi_raw$OND[[i]], paste("OND NDVI", y))
  ar <- rain_panel(rain$ANN[[i]], paste("Annual Rainfall", y))
  an <- ndvi_panel(ndvi_raw$ANN[[i]], paste("Annual NDVI", y))

  panels_42 <- c(panels_42, list(mr, or))
  panels_43 <- c(panels_43, list(mn, on))
  panels_44 <- c(panels_44, list(mr, mn, or, on))
  panels_45 <- c(panels_45, list(ar))
  panels_46 <- c(panels_46, list(an))
  panels_47 <- c(panels_47, list(ar, an))
}
save_grid(panels_42, "Figure_4_2_seasonal_rainfall_maps.png", 8)
save_grid(panels_43, "Figure_4_3_seasonal_NDVI_maps.png", 8)
save_grid(panels_44, "Figure_4_4_seasonal_rainfall_NDVI_maps.png", 8)
save_grid(panels_45, "Figure_4_5_annual_rainfall_maps.png", 7)
save_grid(panels_46, "Figure_4_6_annual_NDVI_maps.png", 7)
save_grid(panels_47, "Figure_4_7_annual_rainfall_NDVI_maps.png", 8)

# --------------------------------------- 7. coefficient of variation (4.4.3)
cv <- function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100
cv_class <- function(v) ifelse(v < 20, "Low", ifelse(v <= 30, "Moderate", "High"))
cv_tab <- data.frame(
  variable = c("Annual rainfall", "MAM rainfall", "OND rainfall",
               "Annual vegetation", "MAM vegetation", "OND vegetation"),
  CV_percent = c(cv(series$ANN_rain), cv(series$MAM_rain), cv(series$OND_rain),
                 cv(series$ANN_ndvi), cv(series$MAM_ndvi), cv(series$OND_ndvi))
) %>% mutate(CV_percent = round(CV_percent, 2), class = cv_class(CV_percent))
write.csv(cv_tab, file.path(TAB_DIR, "CV_table.csv"), row.names = FALSE)

# --------------------------------- 8. trends: graphical, Figures 4.8-4.10
trend_panel <- function(var, ylab, line_col, pt_col, digits) {
  fit <- lm(series[[var]] ~ series$idx)
  b   <- coef(fit)
  eq  <- sprintf("y = %sx + %s", signif(b[2], 5), round(b[1], digits))
  ggplot(series, aes(year, .data[[var]])) +
    geom_line(colour = line_col, linewidth = 1.2) +
    geom_point(colour = pt_col, size = 2.5) +
    geom_smooth(method = "lm", se = FALSE, colour = "black",
                linewidth = 0.6, formula = y ~ x) +
    annotate("text", x = -Inf, y = Inf, label = eq, hjust = -0.1, vjust = 1.4,
             fontface = "bold", size = 5) +
    scale_x_continuous(breaks = YEARS) +
    labs(x = "Years", y = ylab) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, face = "bold"),
          axis.title = element_text(face = "bold"))
}
two_panel <- function(rv, nv, file) {
  p <- trend_panel(rv, "Rainfall (mm)", "blue", "red", 2) |
       trend_panel(nv, "Rangelands", "green4", "orange", 2)
  ggsave(file.path(FIG_DIR, file), p, width = 11, height = 5.5, dpi = 300)
}
two_panel("MAM_rain", "MAM_ndvi", "Figure_4_8_MAM_trend.png")
two_panel("OND_rain", "OND_ndvi", "Figure_4_9_OND_trend.png")
two_panel("ANN_rain", "ANN_ndvi", "Figure_4_10_annual_trend.png")

# ------------------------- 9. trends: Mann-Kendall & Sen's slope (4.5.2)
# Fallback implementation (thesis Equations 1-2, normal approximation with
# continuity correction and tie correction - the same as trend::mk.test).
mk_own <- function(x) {
  x <- x[!is.na(x)]; n <- length(x)
  s <- 0
  for (k in 1:(n - 1)) s <- s + sum(sign(x[(k + 1):n] - x[k]))
  ties <- table(x); ties <- ties[ties > 1]
  var_s <- (n * (n - 1) * (2 * n + 5) - sum(ties * (ties - 1) * (2 * ties + 5))) / 18
  z <- if (s > 0) (s - 1) / sqrt(var_s) else if (s < 0) (s + 1) / sqrt(var_s) else 0
  list(tau = s / (n * (n - 1) / 2), p = 2 * (1 - pnorm(abs(z))), S = s, z = z)
}
sen_own <- function(x) {
  n <- length(x); sl <- c()
  for (k in 1:(n - 1)) for (j in (k + 1):n) sl <- c(sl, (x[j] - x[k]) / (j - k))
  median(sl, na.rm = TRUE)
}
has_trend <- requireNamespace("trend", quietly = TRUE)
has_mmk   <- requireNamespace("modifiedmk", quietly = TRUE)

trend_row <- function(var, label) {
  x <- series[[var]]
  if (has_trend) {
    mk  <- trend::mk.test(x)
    tau <- unname(mk$estimates["tau"]); p <- mk$p.value
    sen <- unname(trend::sens.slope(x)$estimates)
  } else {
    mk  <- mk_own(x); tau <- mk$tau; p <- mk$p; sen <- sen_own(x)
  }
  ols  <- unname(coef(lm(x ~ seq_along(x)))[2])
  acf1 <- acf(x, lag.max = 1, plot = FALSE)$acf[2]
  p_mmk <- if (has_mmk) unname(modifiedmk::mmkh(x)["new P-value"]) else NA
  data.frame(
    Variable = label,
    Tau = round(tau, 4),
    P_value = signif(p, 4),
    Sens_slope = signif(sen, 4),
    OLS_slope = signif(ols, 4),
    Lag1_autocorr = round(acf1, 3),
    Autocorr_significant = abs(acf1) > 1.96 / sqrt(length(x)),
    P_value_modified_MK = signif(p_mmk, 4),
    Trend = ifelse(p < 0.05, "Significant", "Not significant")
  )
}
table_41 <- bind_rows(trend_row("MAM_rain", "MAM Rainfall"),
                      trend_row("OND_rain", "OND Rainfall"),
                      trend_row("MAM_ndvi", "MAM vegetation"),
                      trend_row("OND_ndvi", "OND vegetation"))
table_42 <- bind_rows(trend_row("ANN_rain", "Annual Rainfall"),
                      trend_row("ANN_ndvi", "Annual vegetation"))
write.csv(table_41, file.path(TAB_DIR, "Table_4_1_seasonal_trends.csv"), row.names = FALSE)
write.csv(table_42, file.path(TAB_DIR, "Table_4_2_annual_trends.csv"), row.names = FALSE)

# -------------------------------- 10. regression, Tables 4.3-4.5 (4.6)
reg_table <- function(yv, xv) {
  fit <- lm(series[[yv]] ~ series[[xv]])
  s   <- summary(fit)
  data.frame(
    Statistic = c("Intercept", "Coefficient (Rainfall)", "p-value (Rainfall)",
                  "R-squared", "Adjusted R-squared", "n"),
    Value = c(signif(coef(fit)[1], 5), signif(coef(fit)[2], 5),
              signif(s$coefficients[2, 4], 4), round(s$r.squared, 4),
              round(s$adj.r.squared, 4), nobs(fit))
  )
}
write.csv(reg_table("MAM_ndvi", "MAM_rain"),
          file.path(TAB_DIR, "Table_4_3_regression_MAM.csv"), row.names = FALSE)
write.csv(reg_table("OND_ndvi", "OND_rain"),
          file.path(TAB_DIR, "Table_4_4_regression_OND.csv"), row.names = FALSE)
write.csv(reg_table("ANN_ndvi", "ANN_rain"),
          file.path(TAB_DIR, "Table_4_5_regression_annual.csv"), row.names = FALSE)

# --------------------------- 11. correlation, Figures 4.11-4.12 (4.6)
r_mam <- cor(series$MAM_rain, series$MAM_ndvi)
r_ond <- cor(series$OND_rain, series$OND_ndvi)
r_ann <- cor(series$ANN_rain, series$ANN_ndvi)

scatter <- function(xv, yv, xl, yl, pt, ln, se) {
  ggplot(series, aes(.data[[xv]], .data[[yv]])) +
    geom_point(colour = pt, size = 2.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = se, colour = ln, formula = y ~ x) +
    labs(x = xl, y = yl) +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, colour = "grey30"))
}
fig411 <- scatter("MAM_rain", "MAM_ndvi", "MAM Rainfall (mm)", "MAM Vegetation",
                  "#E74C3C", "blue", FALSE) /
          scatter("OND_rain", "OND_ndvi", "OND Rainfall (mm)", "OND Vegetation",
                  "#E74C3C", "red", FALSE)
ggsave(file.path(FIG_DIR, "Figure_4_11_seasonal_correlation.png"), fig411,
       width = 7, height = 6, dpi = 300)

fig412 <- scatter("ANN_rain", "ANN_ndvi", "Rainfall (mm)", "Vegetation",
                  "blue", "red", TRUE)
ggsave(file.path(FIG_DIR, "Figure_4_12_annual_correlation.png"), fig412,
       width = 7, height = 4, dpi = 300)

cor_tab <- data.frame(
  Scale = c("MAM", "OND", "Annual"),
  Pearson_r = round(c(r_mam, r_ond, r_ann), 4),
  p_value = signif(c(cor.test(series$MAM_rain, series$MAM_ndvi)$p.value,
                     cor.test(series$OND_rain, series$OND_ndvi)$p.value,
                     cor.test(series$ANN_rain, series$ANN_ndvi)$p.value), 4)
)
write.csv(cor_tab, file.path(TAB_DIR, "correlation_table.csv"), row.names = FALSE)

# -------------------------- 12. pixel-wise correlation, Figure 4.13 (4.6.2.2.2)
# Pearson r at every pixel between the 21-year annual rainfall series and the
# 21-year annual NDVI series (all vegetated pixels, as in the thesis map).
if (PIXEL_GRID == "chirps") {
  fact <- round(res(rain$ANN)[1] / res(ndvi_raw$ANN)[1])
  nd   <- aggregate(ndvi_raw$ANN, fact = max(fact, 1), fun = "mean", na.rm = TRUE)
  nd   <- resample(nd, rain$ANN, method = "bilinear")
  rr   <- rain$ANN
} else {
  nd <- ndvi_raw$ANN
  rr <- resample(rain$ANN, nd, method = "bilinear")
}
pix_cor <- app(c(rr, nd), fun = function(v) {
  n <- length(v) / 2
  a <- v[1:n]; b <- v[(n + 1):(2 * n)]
  ok <- !is.na(a) & !is.na(b)
  if (sum(ok) < 10 || sd(a[ok]) == 0 || sd(b[ok]) == 0) return(NA)
  cor(a[ok], b[ok])
})
names(pix_cor) <- "correlation"
writeRaster(pix_cor, file.path(OUT_DIR, "pixelwise_correlation.tif"), overwrite = TRUE)

pc_df <- as.data.frame(pix_cor, xy = TRUE, na.rm = TRUE)
fig413 <- ggplot(pc_df, aes(x, y, fill = correlation)) +
  geom_raster() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1), name = "Correlation") +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, colour = "grey30"))
ggsave(file.path(FIG_DIR, "Figure_4_13_pixelwise_correlation.png"), fig413,
       width = 7, height = 6, dpi = 300)

pix_summary <- data.frame(
  share_positive = mean(pc_df$correlation > 0),
  share_significant_positive = mean(pc_df$correlation > 0.433),  # r crit, n=21, p<.05
  share_negative = mean(pc_df$correlation < 0),
  median_r = median(pc_df$correlation)
)
write.csv(pix_summary, file.path(TAB_DIR, "pixelwise_correlation_summary.csv"), row.names = FALSE)

# ------------------------------------------------------------- 13. summary
cat("\n================ THESIS NUMBERS ================\n")
cat(sprintf("Rangeland mask: %s (%s, NDVI %.2f-%.2f)\n",
            RANGELAND_MASK, MASK_METHOD, NDVI_RANGE[1], NDVI_RANGE[2]))
cat(sprintf("Validation r (CHIRPS vs TMA): %s\n",
            ifelse(is.na(r_val), "not run", sprintf("%.4f", r_val))))
if (exists("val_tab")) { cat("\nValidation metrics:\n"); print(val_tab, row.names = FALSE) }
if (length(suspect_stations) > 0)
  cat("\n!! Not at gauge precision (check with TMA):", paste(suspect_stations, collapse = ", "), "\n")
cat("\nCoefficient of variation (%):\n"); print(cv_tab, row.names = FALSE)
cat("\nTable 4.1 - seasonal trends:\n");   print(table_41, row.names = FALSE)
cat("\nTable 4.2 - annual trends:\n");     print(table_42, row.names = FALSE)
cat(sprintf("\nPearson r: MAM %.3f | OND %.3f | Annual %.3f\n", r_mam, r_ond, r_ann))
cat("\nPixel-wise correlation summary:\n"); print(pix_summary, row.names = FALSE)
cat(sprintf("\nMK/Sen computed with: %s\n",
            ifelse(has_trend, "trend package", "built-in fallback")))
cat("Outputs written to:", normalizePath(OUT_DIR), "\n")
