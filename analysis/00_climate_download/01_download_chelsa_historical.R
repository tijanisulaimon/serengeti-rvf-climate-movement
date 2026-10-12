#!/usr/bin/env Rscript
# 01 - Download CHELSA-monthly historical climate (V2.1) for the Serengeti
#' study domain.
#'
# This is the observed historical baseline paired with the CHELSA-CMIP6 future
# projections. Files are pulled directly from CHELSA's Cloud-Optimized GeoTIFF
# S3 bucket through /vsicurl, so only data inside the study bounding box are
# transferred.
#'
#' Source: https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/monthly/{var}/
#' Citation: Karger, D.N., Brun, P., Zilker, F. (2025). CHELSA-monthly climate
#' data at high resolution. EnviDat. https://doi.org/10.16904/envidat.686
#'
#' Period: 1981-2018 (38 complete years shared by all four variables)
#' Variables: tas, tasmax, tasmin, pr
#' Bounding box: Serengeti study domain + 0.2 degree safety margin
#'
#' Output: data/processed/chelsa_historical/
#'         CHELSA_monthly_{var}_1981_2018.nc
#'         (one NetCDF per variable, 38*12 = 456 monthly layers each)
#'
#' Idempotent: skips per-month TIFs already cached locally; rebuilds the
#' final stitched NetCDF only if a variable's expected layer count is short
#' or the file is missing.

suppressPackageStartupMessages({
  library(terra)
  library(readr)
})

# ----------------------------------------------------------------------
# CONFIG
# ----------------------------------------------------------------------

# Existing future-scenario crop is 34.21652-35.26652 lon, -2.50014--1.41681 lat.
# Pad by 0.2 deg on each side to keep extraction stable for the village +
# resource points and any future bbox tweak in module 01 / 02.
study_bbox <- terra::ext(34.02, 35.47, -2.70, -1.22)

# Although the CHELSA-monthly product page says coverage runs to 2021-12,
# the actual S3 bucket only hosts `pr` through 2019-06. To keep all four
# variables on a clean whole-year schedule we use 1981-2018 (38 years).
historical_years <- 1981:2018
variables <- c("tas", "tasmax", "tasmin", "pr")

base_url <- "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/monthly"

# Paths (resolved relative to project root)
# Layout convention:
#   data/raw/chelsa_historical/per_month_tifs/   per-month cropped GeoTIFF cache
#   data/processed/chelsa_historical/            final stitched per-variable NCs
project_root <- here::here()
cache_dir <- file.path(project_root, "data", "raw", "chelsa_historical",
                       "per_month_tifs")
out_dir <- file.path(project_root, "data", "processed", "chelsa_historical")
log_path <- file.path(project_root, "data", "download_records",
                      "chelsa_historical_download_log.csv")

dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------
# PER-MONTH CACHE
# ----------------------------------------------------------------------

# We cache cropped per-month rasters as small GeoTIFFs so a partial run
# can be resumed cheaply. Layer order in the final stack is enforced by
# the year/month ordering, NOT by file mtime.
cache_path <- function(var, year, month) {
  file.path(cache_dir, sprintf("CHELSA_%s_%02d_%d_V.2.1_crop.tif",
                               var, month, year))
}

remote_url <- function(var, year, month) {
  sprintf("%s/%s/CHELSA_%s_%02d_%d_V.2.1.tif",
          base_url, var, var, month, year)
}

# Fetch one month's cropped raster. Individual range requests occasionally
# time out, so each month receives up to three attempts before the run fails.
fetch_one_month <- function(var, year, month) {
  out_path <- cache_path(var, year, month)
  if (file.exists(out_path)) {
    # Verify the cached file is readable; if it's a corrupt zero-byte file,
    # delete and re-fetch.
    file_ok <- tryCatch({
      r <- terra::rast(out_path)
      ok <- terra::ncell(r) > 0 && all(is.finite(terra::ext(r)[1:4]))
      ok
    }, error = function(e) FALSE)
    if (file_ok) return(invisible(out_path))
    file.remove(out_path)
  }

  url <- remote_url(var, year, month)
  vsi_url <- paste0("/vsicurl/", url)

# CHELSA V2.1 GeoTIFFs store integer-scaled values. Temperature is multiplied
# by 0.1 to obtain kelvin; precipitation is multiplied by 0.01 to obtain
# kg m-2 month-1, numerically equivalent to mm/month. This matches the units
# of the processed CHELSA-CMIP6 inputs.
  scale_factor <- if (var %in% c("tas", "tasmax", "tasmin")) {
    0.1
  } else if (var == "pr") {
    0.01
  } else {
    1
  }
  attempts <- 0
  max_attempts <- 3
  last_err <- NULL
  while (attempts < max_attempts) {
    attempts <- attempts + 1
    res <- tryCatch({
      r <- terra::rast(vsi_url)
      r_crop <- terra::crop(r, study_bbox, snap = "out")
      r_scaled <- r_crop * scale_factor
      terra::writeRaster(r_scaled, out_path, overwrite = TRUE,
                         filetype = "GTiff",
                         gdal = c("COMPRESS=DEFLATE", "ZLEVEL=6"))
      list(ok = TRUE, ncell = terra::ncell(r_scaled))
    }, error = function(e) {
      list(ok = FALSE, err = conditionMessage(e))
    })
    if (isTRUE(res$ok)) return(invisible(out_path))
    last_err <- res$err
    Sys.sleep(2)
  }
  stop(sprintf("Failed to fetch %s after %d attempts: %s\n  url = %s",
               basename(out_path), max_attempts, last_err, url))
}

# ----------------------------------------------------------------------
# STITCH PER-VARIABLE NETCDF
# ----------------------------------------------------------------------

stitched_path <- function(var) {
  file.path(out_dir,
            sprintf("CHELSA_monthly_%s_%d_%d.nc", var,
                    min(historical_years), max(historical_years)))
}

# Stack all months for one variable in chronological order and write a
# multi-band NetCDF with monthly midpoint dates.
stitch_variable <- function(var) {
  out_path <- stitched_path(var)

  # Build time-ordered file list. Use lubridate-free arithmetic to avoid
  # depending on additional packages just for date construction.
  grid <- expand.grid(month = 1:12, year = historical_years)
  grid <- grid[order(grid$year, grid$month), ]
  files <- mapply(cache_path, var, grid$year, grid$month, USE.NAMES = FALSE)

  if (any(!file.exists(files))) {
    stop(sprintf("Cannot stitch %s: %d of %d per-month files missing.",
                 var, sum(!file.exists(files)), length(files)))
  }

  message(sprintf("[%s] Stacking %d monthly layers into %s",
                  var, length(files), basename(out_path)))
  stack <- terra::rast(files)

  # Time encoding: midpoint of each month (15th), seconds since 1970-01-01
  # converted to days. terra encodes time when we set it explicitly.
  dates <- as.Date(sprintf("%d-%02d-15", grid$year, grid$month))
  terra::time(stack) <- dates
  names(stack) <- format(dates, "%Y_%m")

  terra::writeCDF(
    stack, out_path, varname = var, longname = sprintf("CHELSA-monthly %s", var),
    unit = if (var == "pr") "kg m-2 month-1" else "K",
    overwrite = TRUE, compression = 6
  )
  out_path
}

# ----------------------------------------------------------------------
# MAIN LOOP
# ----------------------------------------------------------------------

start_time <- Sys.time()
log_rows <- list()

for (var in variables) {
  message(sprintf("\n=== Variable: %s ===", var))
  fetched <- 0
  cached <- 0
  failed <- 0

  for (year in historical_years) {
    for (month in 1:12) {
      out_path <- cache_path(var, year, month)
      pre_existed <- file.exists(out_path)
      tryCatch({
        fetch_one_month(var, year, month)
        if (pre_existed) cached <- cached + 1L else fetched <- fetched + 1L
      }, error = function(e) {
        failed <<- failed + 1L
        message(sprintf("  FAILED %s %d-%02d: %s", var, year, month, conditionMessage(e)))
      })
    }
    if (year %% 5 == 0 || year == max(historical_years)) {
      message(sprintf("  ... %s %d done (cached=%d, fetched=%d, failed=%d)",
                      var, year, cached, fetched, failed))
    }
  }

  if (failed > 0) {
    stop(sprintf("Variable %s had %d failed months — refusing to stitch a partial NetCDF.",
                 var, failed))
  }

  out_path <- stitch_variable(var)
  fi <- file.info(out_path)
  log_rows[[length(log_rows) + 1L]] <- data.frame(
    variable = var,
    year_min = min(historical_years),
    year_max = max(historical_years),
    n_months = length(historical_years) * 12L,
    n_cached_pre = cached,
    n_fetched = fetched,
    output_path = out_path,
    output_size_mb = round(fi$size / 1024^2, 2),
    stringsAsFactors = FALSE
  )
}

end_time <- Sys.time()
elapsed_min <- as.numeric(difftime(end_time, start_time, units = "mins"))

download_log <- do.call(rbind, log_rows)
download_log$run_finished_at <- format(end_time, "%Y-%m-%d %H:%M:%S")
download_log$elapsed_minutes <- round(elapsed_min, 2)
write_csv(download_log, log_path)

message(sprintf("\nDone. Elapsed %.1f min. Download log: %s", elapsed_min, log_path))
print(download_log)
