#' ####################################################################### #
#' MODULE 01 — NDVI Prediction Configuration
#' PURPOSE: Paths, parameters, and helper functions for the CHELSA-only
#'          NDVI prediction pipeline (GIMMS-3G+ + CHELSA-monthly/CHELSA-CMIP6).
#' ####################################################################### #

source(here::here("_project_paths.R"))
source(here::here("analysis", "_future_climate_config.R"))

make_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  path
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

ndvi_config <- local({
  data_dir      <- proj$data_dir
  raw_dir       <- make_dir(file.path(data_dir, "raw"))
  observed_dir  <- make_dir(file.path(data_dir, "observed"))
  records_dir   <- make_dir(file.path(data_dir, "download_records"))

  serengeti_buffer_km <- 5

  chelsa_future_root <- file.path(data_dir, "processed", "chelsa_cmip6_future")
  gimms_support_dir <- make_dir(file.path(proj$ndvi_derived,
                                          "gimms_support"))
  future_features_root <- make_dir(file.path(gimms_support_dir, "future"))
  scenario_outputs_root <- make_dir(file.path(proj$ndvi_outputs,
                                              "scenarios"))

  list(
    data_dir = data_dir,
    raw_dir = raw_dir,
    observed_dir = observed_dir,
    records_dir = records_dir,
    serengeti_buffer_km = serengeti_buffer_km,

    # module directories
    derived_dir = make_dir(proj$ndvi_derived),
    outputs_dir = make_dir(proj$ndvi_outputs),
    gimms_support_dir = gimms_support_dir,
    future_features_root = future_features_root,
    scenario_outputs_root = scenario_outputs_root,

    # time range
    start_hist = as.Date("1982-01-01"),
    end_hist = as.Date("2022-12-01"),
    temporal_holdout_start_year = 2015L,

    # Future climate
    # Use the same four quality-screened GCMs for both SSPs. MPI-ESM1-2-LR was
    # downloaded as a technical substitute for MPI-ESM1-2-HR but is retained
    # only for climate QA because its SSP2-4.5 `tas` channel is inconsistent
    # with `tasmin` and `tasmax`.
    future_scenarios = future_climate$scenarios,
    primary_future_scenario = future_climate$primary_scenario,
    future_scenario_labels = future_climate$scenario_labels,
    projection_start_year = future_climate$projection_start_year,
    analysis_windows = future_climate$analysis_windows,
    chelsa_models = future_climate$reporting_gcms,
    chelsa_quality_models = future_climate$download_gcms,
    chelsa_reporting_models = future_climate$reporting_gcms,
    chelsa_future_root = chelsa_future_root,

    # CHELSA-monthly historical (Karger et al. 2025) — single observational
    # source for all four climate variables, replacing the previous
    # CHIRPS+CHIRTS pipeline. Pairs cleanly with CHELSA-CMIP6 future
    # (delta-change downscaling on the same V2.1 backbone).
    chelsa_historical_dir = file.path(data_dir, "processed",
                                      "chelsa_historical"),
    path_chelsa_pr = file.path(data_dir, "processed", "chelsa_historical",
                               "CHELSA_monthly_pr_1981_2018.nc"),
    path_chelsa_tas = file.path(data_dir, "processed", "chelsa_historical",
                                "CHELSA_monthly_tas_1981_2018.nc"),
    path_chelsa_tasmax = file.path(data_dir, "processed", "chelsa_historical",
                                   "CHELSA_monthly_tasmax_1981_2018.nc"),
    path_chelsa_tasmin = file.path(data_dir, "processed", "chelsa_historical",
                                   "CHELSA_monthly_tasmin_1981_2018.nc"),

    gimms_collection_concept_id = "C2759076389-ORNL_CLOUD",
    gimms_cmr_granules_url = "https://cmr.earthdata.nasa.gov/search/granules.json",

    # GIMMS
    gimms_raw_dir = make_dir(file.path(raw_dir, "GIMMS_3Gplus")),
    gimms_granule_list = file.path(records_dir, "gimms_3gplus_granules.csv"),

    # shapefiles & elevation
    shapefile_path = file.path(data_dir, "shapefiles", "gadm41_TZA_2_pk.rds"),
    elevation_path = file.path(data_dir, "elevation", "serengeti_ele.tif")
  )
})

ndvi_future_climate_dir <- function(scenario, cfg = ndvi_config) {
  file.path(cfg$chelsa_future_root, scenario)
}

ndvi_future_features_dir <- function(scenario, cfg = ndvi_config,
                                     create = FALSE) {
  path <- file.path(cfg$future_features_root, scenario)
  if (create) make_dir(path) else path
}

ndvi_scenario_outputs_dir <- function(scenario, cfg = ndvi_config,
                                      create = FALSE) {
  path <- file.path(cfg$scenario_outputs_root, scenario)
  if (create) make_dir(path) else path
}

read_serengeti_boundary <- function(cfg = ndvi_config) {
  tz_shape <- terra::vect(cfg$shapefile_path)
  serengeti <- tz_shape[tz_shape$NAME_2 == "Serengeti"]
  terra::crs(serengeti) <- terra::crs("EPSG:4326")

  if (cfg$serengeti_buffer_km > 0) {
    serengeti_utm <- terra::project(serengeti, "EPSG:32736")
    serengeti_utm <- terra::buffer(serengeti_utm, width = cfg$serengeti_buffer_km * 1000)
    serengeti <- terra::project(serengeti_utm, "EPSG:4326")
  }

  serengeti
}
