#' ####################################################################### #
#' 01c — Prepare Historical Feature Tables
#'
#' Combine GIMMS-3G+ NDVI with CHELSA-monthly historical climate (pr, tasmax,
#' tasmin) and CHELSA-CMIP6 future projections onto a common monthly grid.
#' Produces the feather files used by all downstream modeling scripts.
#'
#' Historical climate switched from CHIRPS+CHIRTS to CHELSA-monthly
#' (Karger et al. 2025) in 2026-05 so that historical and future come from
#' the same V2.1 backbone, removing a dataset-mismatch source of bias.
#'
#' Historical overlap: 1981-01 to 2018-12 (CHELSA-monthly coverage).
#' Working grid: GIMMS cells that also have CHELSA coverage for all GCMs.
#'
#' INPUTS:  data/raw/GIMMS_3Gplus/,
#'          data/processed/chelsa_historical/CHELSA_monthly_{var}_1981_2018.nc
#'          data/processed/chelsa_cmip6_future/{scenario}/{GCM}_{var}.nc
#' OUTPUTS: derived/gimms_support/all_historical_observed.feather
#'          derived/gimms_support/future/{scenario}/all_future_{GCM}.feather
#'          (only when the corresponding processed climate files are available)
#'          derived/gimms_support/prep_observed_gimms_support_summary.csv
#'
#' Behaviour: future feathers are skipped (with a clear message) when the
#' CHELSA-CMIP6 processed files are not yet on disk, so this script can be
#' run as soon as historical is ready and re-run later for the future side.
#' ####################################################################### #

suppressPackageStartupMessages({
  library(terra)
  library(tidyverse)
  library(lubridate)
  library(arrow)
  library(readr)
})

source(here::here("analysis", "01_ndvi_prediction", "_config.R"))

cfg <- ndvi_config
Dir.Data <- cfg$data_dir
Dir.Out <- cfg$gimms_support_dir

model_names <- cfg$chelsa_models
future_scenarios <- cfg$future_scenarios

read_climate_nc <- function(filepath, kelvin_to_celsius = FALSE) {
  r <- rast(filepath)
  if (kelvin_to_celsius) {
    r <- r - 273.15
  }
  if (all(is.na(time(r)))) {
    stop("Missing time dimension in: ", filepath)
  }
  names(r) <- format(as.Date(time(r)), "%b_%Y")
  crs(r) <- crs("EPSG:4326")
  r
}

crop_mask <- function(r, boundary) {
  mask(crop(r, boundary), boundary)
}

to_df <- function(r, varname) {
  # Floor `date` to the first of each month so the join across data sources
  # works regardless of which day each NetCDF anchored its time axis to.
  # GIMMS uses month-start (-01); CHELSA-monthly uses month-midpoint (-15);
  # CHELSA-CMIP6 future uses (-15) as well. Without this floor the join by
  # date silently produces all-NA climate columns.
  as.data.frame(r, xy = TRUE, time = TRUE, wide = FALSE) |>
    rename(!!varname := values, date = time) |>
    mutate(date = lubridate::floor_date(as.Date(date), "month")) |>
    select(x, y, date, all_of(varname))
}

load_gimms_monthly <- function(gimms_dir, boundary) {
  gimms_files <- list.files(gimms_dir, pattern = "\\.nc4$", full.names = TRUE) |>
    sort()

  if (!length(gimms_files)) {
    stop("No GIMMS granules found in: ", gimms_dir)
  }

  message("Loading ", length(gimms_files), " GIMMS granules...")

  ndvi_parts <- lapply(gimms_files, function(path) {
    r <- rast(path)
    ndvi <- r[[grepl("^ndvi_", names(r))]]
    crs(ndvi) <- crs("EPSG:4326")
    ndvi <- crop_mask(ndvi, boundary)
    ndvi <- ndvi / 10000
    ndvi <- ifel(ndvi < -0.3 | ndvi > 1, NA, ndvi)
    ndvi
  })

  ndvi_halfmonthly <- do.call(c, ndvi_parts)
  ndvi_times <- as.Date(time(ndvi_halfmonthly))
  ord <- order(ndvi_times)
  ndvi_halfmonthly <- ndvi_halfmonthly[[ord]]
  ndvi_times <- ndvi_times[ord]

  month_keys <- format(ndvi_times, "%Y-%m")
  month_index <- factor(month_keys, levels = unique(month_keys))
  monthly_ndvi <- tapp(ndvi_halfmonthly, month_index, fun = mean, na.rm = TRUE)

  monthly_dates <- as.Date(paste0(levels(month_index), "-01"))
  time(monthly_ndvi) <- monthly_dates
  names(monthly_ndvi) <- format(monthly_dates, "%b_%Y")
  monthly_ndvi
}

compute_future_common_extent <- function(future_dirs, models) {
  # CHELSA-CMIP6 outputs are on a fixed 30-arcsec grid for a given bbox, so
  # the per-GCM extents are identical when present. We still take the
  # intersection of whichever files do exist — and fall back to the
  # downloader's bbox (analysis/00_climate_download/02_*.py) so 01c can
  # produce historical feathers even when downloads are mid-flight.
  files <- unlist(
    lapply(
      future_dirs,
      function(fut_dir) file.path(fut_dir, paste0(models, "_pr.nc"))
    ),
    use.names = FALSE
  )
  files <- files[file.exists(files)]

  if (length(files) == 0) {
    message(
      "No CHELSA-CMIP6 future _pr.nc files found yet (downloads in progress?); ",
      "falling back to the configured download bbox for the working grid."
    )
    return(ext(34.02, 35.47, -2.70, -1.22))
  }

  common_ext <- NULL
  for (f in files) {
    r <- rast(f)
    crs(r) <- crs("EPSG:4326")
    current_ext <- ext(r)
    common_ext <- if (is.null(common_ext)) {
      current_ext
    } else {
      intersect(common_ext, current_ext)
    }
  }
  common_ext
}

boundary <- read_serengeti_boundary(cfg)

hist_obs_dir <- cfg$observed_dir
gimms_monthly <- load_gimms_monthly(cfg$gimms_raw_dir, boundary)

# CHELSA-monthly historical (one product, K-encoded temps).
# 01c keeps the feather schema {pr, tasmax, tasmin} (downstream computes
# tasmean from tmin/tmax) so we only need three of the four CHELSA vars.
hist_pr <- rast(cfg$path_chelsa_pr)
hist_tmax <- rast(cfg$path_chelsa_tasmax) - 273.15
hist_tmin <- rast(cfg$path_chelsa_tasmin) - 273.15

# CHELSA-monthly covers 1981-2018; GIMMS-3G+ covers 1981-2018 too. The
# shared start is whatever each product happens to begin reporting on.
common_start <- max(
  min(as.Date(time(gimms_monthly))),
  min(as.Date(time(hist_pr))),
  min(as.Date(time(hist_tmax))),
  min(as.Date(time(hist_tmin)))
)
common_end <- min(
  max(as.Date(time(gimms_monthly))),
  max(as.Date(time(hist_pr))),
  max(as.Date(time(hist_tmax))),
  max(as.Date(time(hist_tmin))),
  cfg$end_hist
)

gimms_monthly <- gimms_monthly[[which(as.Date(time(gimms_monthly)) >= common_start &
                                      as.Date(time(gimms_monthly)) <= common_end)]]
hist_pr <- hist_pr[[which(as.Date(time(hist_pr)) >= common_start &
                          as.Date(time(hist_pr)) <= common_end)]]
hist_tmax <- hist_tmax[[which(as.Date(time(hist_tmax)) >= common_start &
                              as.Date(time(hist_tmax)) <= common_end)]]
hist_tmin <- hist_tmin[[which(as.Date(time(hist_tmin)) >= common_start &
                              as.Date(time(hist_tmin)) <= common_end)]]

future_climate_dirs <- vapply(
  future_scenarios,
  ndvi_future_climate_dir,
  character(1),
  cfg = cfg
)
future_common_ext <- compute_future_common_extent(
  future_climate_dirs,
  model_names
)

template <- crop(gimms_monthly[[1]], future_common_ext)
template <- mask(template, boundary)
gimms_monthly <- crop(gimms_monthly, template)
gimms_monthly <- mask(gimms_monthly, template)

hist_pr_g <- resample(hist_pr, template, method = "average")
hist_tmax_g <- resample(hist_tmax, template, method = "average")
hist_tmin_g <- resample(hist_tmin, template, method = "average")

ele <- rast(cfg$elevation_path)
crs(ele) <- crs("EPSG:4326")
ele <- crop_mask(ele, boundary)
ele_g <- resample(ele, template, method = "average")
names(ele_g) <- "ele"

ndvi_df <- to_df(gimms_monthly, "ndvi")
pr_df <- to_df(hist_pr_g, "pr")
tasmax_df <- to_df(hist_tmax_g, "tasmax")
tasmin_df <- to_df(hist_tmin_g, "tasmin")
ele_df <- as.data.frame(ele_g, xy = TRUE, wide = FALSE) |>
  rename(ele = values) |>
  select(x, y, ele)

all_hist <- ndvi_df |>
  left_join(pr_df, by = c("x", "y", "date")) |>
  left_join(tasmin_df, by = c("x", "y", "date")) |>
  left_join(tasmax_df, by = c("x", "y", "date")) |>
  left_join(ele_df, by = c("x", "y")) |>
  drop_na() |>
  arrange(x, y, date)

write_feather(all_hist, file.path(Dir.Out, "all_historical_observed.feather"))

prepare_future_model <- function(model, scenario, template, boundary, ele_df,
                                 out_dir, fut_dir) {
  pr_file <- file.path(fut_dir, paste0(model, "_pr.nc"))
  tmin_file <- file.path(fut_dir, paste0(model, "_tasmin.nc"))
  tmax_file <- file.path(fut_dir, paste0(model, "_tasmax.nc"))

  if (!all(file.exists(c(pr_file, tmin_file, tmax_file)))) {
    message("[", scenario, " | ", model,
            "] CHELSA-CMIP6 processed files not yet on disk; ",
            "skipping (re-run 01c after 00_climate_download/03_process_*.R).")
    return(NULL)
  }

  message("Preparing CHELSA future forcing for: ", scenario, " | ", model)
  fut_pr <- read_climate_nc(pr_file)
  fut_tmin <- read_climate_nc(tmin_file, kelvin_to_celsius = TRUE)
  fut_tmax <- read_climate_nc(tmax_file, kelvin_to_celsius = TRUE)

  fut_pr <- crop_mask(fut_pr, boundary)
  fut_tmin <- crop_mask(fut_tmin, boundary)
  fut_tmax <- crop_mask(fut_tmax, boundary)

  fut_pr_g <- resample(fut_pr, template, method = "average")
  fut_tmin_g <- resample(fut_tmin, template, method = "average")
  fut_tmax_g <- resample(fut_tmax, template, method = "average")

  fut_pr_df <- to_df(fut_pr_g, "pr")
  fut_tmin_df <- to_df(fut_tmin_g, "tasmin")
  fut_tmax_df <- to_df(fut_tmax_g, "tasmax")

  fut_df <- fut_pr_df |>
    left_join(fut_tmin_df, by = c("x", "y", "date")) |>
    left_join(fut_tmax_df, by = c("x", "y", "date")) |>
    left_join(ele_df, by = c("x", "y")) |>
    drop_na() |>
    arrange(x, y, date)

  write_feather(
    fut_df,
    file.path(out_dir, paste0("all_future_", model, ".feather"))
  )

  tibble(
    scenario = scenario,
    climate_model = model,
    future_rows = nrow(fut_df),
    future_cells = n_distinct(fut_df$x, fut_df$y),
    future_start = as.character(min(fut_df$date)),
    future_end = as.character(max(fut_df$date))
  )
}

future_summary <- purrr::map_dfr(future_scenarios, function(scenario) {
  purrr::map_dfr(
    model_names,
    prepare_future_model,
    scenario = scenario,
    template = template,
    boundary = boundary,
    ele_df = ele_df,
    out_dir = ndvi_future_features_dir(scenario, cfg, create = TRUE),
    fut_dir = ndvi_future_climate_dir(scenario, cfg)
  )
})

hist_summary <- tibble(
  common_start = as.character(common_start),
  common_end = as.character(common_end),
  historical_rows = nrow(all_hist),
  historical_cells = n_distinct(all_hist$x, all_hist$y),
  template_xres = res(template)[1],
  template_yres = res(template)[2],
  buffer_km = cfg$serengeti_buffer_km
)

summary_df <- if (nrow(future_summary) > 0L) {
  crossing(hist_summary, future_summary)
} else {
  hist_summary
}

write_csv(summary_df,
          file.path(Dir.Out, "prep_observed_gimms_support_summary.csv"))

message("Observed-history preparation complete.")
message("Output directory: ", Dir.Out)
message("Historical overlap: ", common_start, " to ", common_end)
if (nrow(future_summary) == 0L) {
  message("No future GCM feathers written — re-run 01c after ",
          "00_climate_download/03_process_chelsa_cmip6_future.R completes.")
}
