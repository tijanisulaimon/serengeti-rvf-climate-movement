#!/usr/bin/env Rscript
#' 03 — Process raw chelsa_cmip6 per-year NetCDFs into per-(GCM, variable)
#' monthly time-series NetCDFs that match the structure of the historical
#' CHELSA-monthly products in `data/processed/chelsa_historical/`.
#'
#' Input layout (written by 02_generate_chelsa_cmip6_future.py):
#'   data/raw/chelsa_cmip6_future/{scenario}/{GCM}/{year}/
#'     CHELSA_{INST}_{GCM}_{var}_{scenario}_{m}_{Y}-01-01_{Y}-12-{31|30}.nc
#'   plus a 1981-2010 reference-climatology NC per variable that we ignore
#'   (these are identical across years).
#'
#' Output:
#'   data/processed/chelsa_cmip6_future/{scenario}/{GCM}_{var}.nc
#'   1032 monthly layers (86 years * 12 months), 2015-01..2100-12, time axis
#'   set to month midpoints (the 15th).
#'
#' Behaviour:
#'   - Processes every scenario listed in analysis/_future_climate_config.R.
#'   - For each (scenario, GCM, variable) combination, only writes the
#'     stitched NetCDF when all expected years are available on disk.
#'     Incomplete combinations are reported and skipped.
#'   - Re-uses an existing output file if it is newer than every per-year
#'     input (idempotent across re-runs).

suppressPackageStartupMessages({
  library(terra)
})

# ---- CONFIG -----------------------------------------------------------------
source(here::here("analysis", "_future_climate_config.R"))

cfg <- list(
  project_root = here::here(),
  scenarios  = future_climate$scenarios,
  gcms       = future_climate$download_gcms,
  variables  = c("tas", "tasmax", "tasmin", "pr"),
  years      = future_climate$climate_years,
  raw_root   = file.path(here::here(), "data", "raw",
                         "chelsa_cmip6_future"),
  out_root   = file.path(here::here(), "data", "processed",
                         "chelsa_cmip6_future")
)

# ---- HELPERS ----------------------------------------------------------------

# Shared configuration uses underscore-separated identifiers in filenames and
# analysis tables. Raw CHELSA directories retain the official hyphenated GCM
# names used by the downloader.
gcm_to_raw_name <- function(gcm) gsub("_", "-", gcm, fixed = TRUE)

# Locate the per-year NetCDF for one (gcm, var, year).
# Filenames follow chelsa_cmip6's convention. The final day is 31 for normal
# calendars and 30 for UKESM's 360-day calendar, so we wildcard the day.
year_nc_path <- function(gcm, var, year, raw_base, experiment) {
  dir <- file.path(raw_base, gcm_to_raw_name(gcm), year)
  if (!dir.exists(dir)) return(character(0))
  patt <- sprintf("^CHELSA_.*_%s_%s_.*_%d-01-01_%d-12-..[.]nc$",
                  var, experiment, year, year)
  list.files(dir, pattern = patt, full.names = TRUE)
}

# Build a year-ordered list of per-year input paths. Returns NA where missing.
collect_year_files <- function(gcm, var, years, raw_base, experiment) {
  vapply(years, function(y) {
    fs <- year_nc_path(gcm, var, y, raw_base, experiment)
    if (length(fs) > 1L) {
      stop(
        "Expected one ", experiment, " ", gcm, " ", var, " file for ", y,
        " but found ", length(fs), ":\n",
        paste(fs, collapse = "\n")
      )
    }
    if (length(fs) == 0L) NA_character_ else fs[[1]]
  }, character(1))
}

# Output filenames use underscores instead of dashes (e.g. GFDL_ESM4_tas.nc),
# matching module 03's `climate_models` identifier convention. Raw input dirs
# keep the chelsa_cmip6 dashed form (GFDL-ESM4).
gcm_to_id <- function(gcm) gsub("-", "_", gcm, fixed = TRUE)

# Process one (gcm, var) pair: skip if incomplete OR if output is fresh.
process_gcm_var <- function(gcm, var, experiment, raw_base, out_dir,
                            years = cfg$years) {
  out_path <- file.path(out_dir,
                        sprintf("%s_%s.nc", gcm_to_id(gcm), var))

  per_year_files <- collect_year_files(
    gcm,
    var,
    years,
    raw_base,
    experiment
  )
  n_present <- sum(!is.na(per_year_files))
  n_total <- length(years)
  expected_n <- n_total * 12L

  result_row <- function(status, n_layers = NA_integer_) {
    data.frame(
      scenario = experiment,
      gcm = gcm,
      variable = var,
      years_present = n_present,
      years_expected = n_total,
      layers = n_layers,
      status = status,
      output_file = out_path,
      stringsAsFactors = FALSE
    )
  }

  if (n_present < n_total) {
    message(sprintf(
      "[%s | %s] %d/%d years on disk — skipping (download in progress)",
      gcm, var, n_present, n_total))
    return(result_row("incomplete_download"))
  }

  # Reuse only a current output with the expected monthly layer count.
  if (file.exists(out_path)) {
    out_m <- file.info(out_path)$mtime
    in_m <- file.info(per_year_files)$mtime
    if (all(out_m >= in_m, na.rm = TRUE)) {
      existing_n <- tryCatch(
        terra::nlyr(terra::rast(out_path)),
        error = function(e) NA_integer_
      )
      if (identical(existing_n, expected_n)) {
        message(sprintf("[%s | %s] output already up-to-date, skipping",
                        gcm, var))
        return(result_row("up_to_date", existing_n))
      }
      message(sprintf(
        "[%s | %s] existing output has %s layers; rebuilding",
        gcm,
        var,
        ifelse(is.na(existing_n), "an unreadable number of", existing_n)
      ))
    }
  }

  message(sprintf("[%s | %s] stacking %d years (%d-%d) -> %s",
                  gcm, var, n_total, min(years), max(years),
                  basename(out_path)))

  # Stack all per-year rasters. Each year file has 12 monthly layers; terra's
  # SpatRaster `c()` concatenates layers, preserving order.
  # We rebuild the stack from raw values so that conflicting time/calendar
  # metadata from the individual per-year NetCDFs (UKESM uses a 360-day
  # calendar; others use proleptic Gregorian) does not survive into the
  # concatenated object — writeCDF otherwise rejects it.
  stacks <- lapply(per_year_files, terra::rast)
  full_with_meta <- do.call(c, stacks)
  template <- stacks[[1]][[1]]                # one layer for geometry only
  full <- terra::rast(template, nlyrs = terra::nlyr(full_with_meta))
  terra::values(full) <- terra::values(full_with_meta)

  if (terra::nlyr(full) != expected_n) {
    stop(sprintf(
      "[%s | %s] expected %d layers, got %d (unexpected month count)",
      gcm, var, expected_n, terra::nlyr(full)))
  }

  # Anchor the time axis to month midpoints (matches the historical
  # CHELSA-monthly product structure). The fresh-rebuild above already
  # detached the per-year metadata, so plain time<- works here.
  dates <- as.Date(unlist(lapply(years, function(y) {
    sprintf("%d-%02d-15", y, 1:12)
  })))
  terra::time(full) <- dates
  names(full) <- format(dates, "%Y_%m")

  unit_str <- if (var == "pr") "kg m-2 month-1" else "K"
  longname <- sprintf("CHELSA-CMIP6 monthly %s, %s, %s",
                      var, gcm, experiment)

  terra::writeCDF(
    full, out_path, varname = var, longname = longname, unit = unit_str,
    overwrite = TRUE, compression = 6
  )

  fi <- file.info(out_path)
  message(sprintf("  written: %s  (%.1f MB)",
                  basename(out_path), fi$size / 1024^2))
  result_row("written", terra::nlyr(full))
}

# ---- MAIN -------------------------------------------------------------------

start_t <- Sys.time()
n_done <- 0L
n_skipped <- 0L
status_rows <- list()
for (experiment in cfg$scenarios) {
  raw_base <- file.path(cfg$raw_root, experiment)
  out_dir <- file.path(cfg$out_root, experiment)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  message("\nProcessing CHELSA-CMIP6 scenario: ", experiment)
  message("Raw directory: ", raw_base)
  message("Output directory: ", out_dir)

  for (gcm in cfg$gcms) {
    message(sprintf("\n=== %s | %s ===", experiment, gcm))
    for (var in cfg$variables) {
      result <- process_gcm_var(
        gcm,
        var,
        experiment,
        raw_base,
        out_dir
      )
      status_rows[[length(status_rows) + 1L]] <- result
      if (identical(result$status, "incomplete_download")) {
        n_skipped <- n_skipped + 1L
      } else {
        n_done <- n_done + 1L
      }
    }
  }
}
processing_status <- do.call(rbind, status_rows)
status_path <- file.path(cfg$out_root, "processing_status.csv")
utils::write.csv(processing_status, status_path, row.names = FALSE)

elapsed <- as.numeric(difftime(Sys.time(), start_t, units = "mins"))
message(sprintf("\nDone. available=%d incomplete=%d elapsed=%.1f min",
                n_done, n_skipped, elapsed))
message("Processing status: ", status_path)
