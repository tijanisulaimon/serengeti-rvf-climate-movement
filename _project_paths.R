#' Shared project paths — source this from any script.
#' Uses here::here() anchored to the .Rproj file at the project root.

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Install the 'here' package: install.packages('here')")
}

proj <- list(
  root = here::here(),

  # --- raw external data (downloaded by 01 scripts) -------------------------
  data_dir        = here::here("data"),
  raw_dir         = here::here("data", "raw"),
  observed_dir    = here::here("data", "observed"),
  chelsa_dir      = here::here("data", "chelsa"),
  shapefiles_dir  = here::here("data", "shapefiles"),
  elevation_dir   = here::here("data", "elevation"),
  movement_dir    = here::here("data", "observed_movements"),

  # --- module 01: NDVI prediction -------------------------------------------
  ndvi_dir     = here::here("analysis", "01_ndvi_prediction"),
  ndvi_derived = here::here("analysis", "01_ndvi_prediction", "derived"),
  ndvi_outputs = here::here("analysis", "01_ndvi_prediction", "outputs"),

  # --- module 02: movement networks ----------------------------------------
  mvt_dir     = here::here("analysis", "02_movement_networks"),
  mvt_outputs = here::here("analysis", "02_movement_networks", "outputs"),

  # --- module 03: RVF transmission ------------------------------------------
  rvf_dir      = here::here("analysis", "03_rvf_transmission"),
  rvf_data     = here::here("analysis", "03_rvf_transmission", "data"),
  rvf_outputs  = here::here("analysis", "03_rvf_transmission", "outputs")
)
