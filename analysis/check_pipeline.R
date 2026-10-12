#!/usr/bin/env Rscript
# Check that the completed outputs passed between Modules 00-03 are present and
# consistent with the shared SSP, GCM, and reporting-window configuration.
# This is a light check: it reads the saved index and summary tables but does
# not rerun models.

project_root <- normalizePath(getwd(), mustWork = TRUE)
shared_config <- file.path(project_root, "analysis", "_future_climate_config.R")

if (!file.exists(shared_config)) {
  stop(
    "Run this script from the repository root:\n",
    "  Rscript analysis/check_pipeline.R",
    call. = FALSE
  )
}

source(shared_config)
module03_config <- file.path(
  project_root, "analysis", "03_rvf_transmission", "R", "config.R"
)
source(module03_config)

check <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

read_required_table <- function(path, label) {
  check(file.exists(path), paste0(label, " is missing: ", path))
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

resolve_output_path <- function(path) {
  if (file.exists(path)) return(path)
  file.path(project_root, path)
}

check_files <- function(paths, label) {
  resolved <- vapply(paths, resolve_output_path, character(1))
  missing <- paths[!file.exists(resolved)]
  check(
    length(missing) == 0L,
    paste0(label, " missing ", length(missing), " file(s), including: ", missing[1])
  )
}

expected_pairs <- as.vector(outer(
  future_climate$scenarios,
  future_climate$reporting_gcms,
  paste,
  sep = "::"
))

# Module 00: every downloaded GCM has four complete monthly climate channels.
climate_status <- read_required_table(
  file.path(
    project_root,
    "data", "processed", "chelsa_cmip6_future", "processing_status.csv"
  ),
  "Module 00 processing status"
)
expected_climate_rows <- length(future_climate$scenarios) *
  length(future_climate$download_gcms) * 4L
check(nrow(climate_status) == expected_climate_rows,
      "Module 00 processing status has an unexpected number of rows.")
check(all(climate_status$status %in% c("written", "up_to_date")),
      "Module 00 contains incomplete climate products.")
check(all(climate_status$layers == length(future_climate$climate_years) * 12L),
      "Module 00 climate products do not all have the expected monthly layers.")
check_files(climate_status$output_file, "Module 00")

quality_checks <- read_required_table(
  file.path(
    project_root, "analysis", "00_climate_download", "outputs",
    "climate_quality", "future_climate_quality_checks.csv"
  ),
  "Module 00 climate-quality table"
)
quality_reporting <- quality_checks[
  quality_checks$climate_model %in% future_climate$reporting_gcms,
]
quality_pairs <- paste(
  quality_reporting$scenario,
  quality_reporting$climate_model,
  sep = "::"
)
check(setequal(quality_pairs, expected_pairs),
      "Module 00 quality checks do not cover the reporting ensemble.")
check(all(quality_reporting$review_status == "pass"),
      "At least one reporting GCM failed the Module 00 quality screen.")
cat("PASS Module 00: climate products and reporting ensemble\n")

# Module 01: one complete NDVI projection for every SSP-GCM pair.
ndvi_index <- read_required_table(
  file.path(
    project_root, "analysis", "01_ndvi_prediction", "outputs",
    "scenarios", "projection_index.csv"
  ),
  "Module 01 projection index"
)
ndvi_pairs <- paste(ndvi_index$scenario, ndvi_index$climate_model, sep = "::")
check(setequal(ndvi_pairs, expected_pairs),
      "Module 01 projection index does not match the reporting ensemble.")
check(all(ndvi_index$status %in% c("written", "reused")),
      "Module 01 contains incomplete NDVI projections.")
check(all(ndvi_index$start == "2019-01-01"),
      "Module 01 projections do not all begin in January 2019.")
check(all(ndvi_index$end == "2100-12-01"),
      "Module 01 projections do not all end in December 2100.")
check_files(ndvi_index$output, "Module 01")
check(
  file.exists(file.path(
    project_root, "analysis", "01_ndvi_prediction", "outputs", "models",
    "final", "model_fit.rds"
  )),
  "Module 01 final model object is missing."
)
cat("PASS Module 01: fitted model and eight future projections\n")

# Module 02: ten observed years and every year in the three future windows.
bundle_index <- read_required_table(
  file.path(
    project_root, "analysis", "02_movement_networks", "outputs",
    "projections", "summaries", "network_bundle_index.csv"
  ),
  "Module 02 network-bundle index"
)
check(all(bundle_index$status %in% c("written", "existing")),
      "Module 02 contains incomplete network bundles.")

observed_networks <- bundle_index[bundle_index$dataset == "observed",]
future_networks <- bundle_index[bundle_index$dataset == "future",]
check(
  nrow(observed_networks) == 10L &&
    !anyDuplicated(observed_networks$year) &&
    setequal(observed_networks$year, 2009:2018),
      "Module 02 observed bundles do not cover 2009-2018 exactly once.")

future_years <- sort(unique(unlist(future_climate$analysis_windows)))
expected_future_rows <- length(future_climate$scenarios) *
  length(future_climate$reporting_gcms) * length(future_years)
check(nrow(future_networks) == expected_future_rows,
      "Module 02 future bundle count does not match the configured design.")
network_pairs <- unique(paste(
  future_networks$scenario,
  future_networks$climate_model,
  sep = "::"
))
check(setequal(network_pairs, expected_pairs),
      "Module 02 future bundles do not match the reporting ensemble.")
check(all(vapply(
  split(future_networks$year, paste(
    future_networks$scenario,
    future_networks$climate_model,
    sep = "::"
  )),
  function(year) setequal(year, future_years),
  logical(1)
)), "At least one Module 02 SSP-GCM pair has incomplete yearly coverage.")
check_files(bundle_index$output_file, "Module 02 bundles")
check_files(bundle_index$summary_file, "Module 02 job summaries")
cat("PASS Module 02: 490 observed and future network bundles\n")

# Module 03 inputs: the materialised bundle grid must match the current
# publication configuration before any expensive transmission run starts.
publication_replicates <- sort(unique(as.integer(
  rvf_config$analysis_network_replicates
)))
network_index <- read_required_table(
  file.path(rvf_config$networks_dir, "network_scenario_index.csv"),
  "Module 03 network-scenario index"
)
scenario_index <- read_required_table(
  file.path(rvf_config$scenarios_dir, "scenario_bundle_index.csv"),
  "Module 03 scenario-bundle index"
)

expected_observed_scenarios <-
  length(rvf_config$current_years) * length(publication_replicates)
expected_future_scenarios <-
  length(rvf_config$future_scenarios) *
  length(rvf_config$climate_models) *
  length(rvf_config$future_years) *
  length(publication_replicates)
expected_scenarios <- expected_observed_scenarios + expected_future_scenarios

check(nrow(network_index) == expected_scenarios,
      "Module 03 network index does not match the current replicate design.")
check(nrow(scenario_index) == expected_scenarios,
      "Module 03 scenario index does not match the current replicate design.")
check(!anyDuplicated(network_index$scenario_id),
      "Module 03 network index contains duplicated scenario IDs.")
check(!anyDuplicated(scenario_index$scenario_id),
      "Module 03 scenario index contains duplicated scenario IDs.")
check(
  identical(
    sort(unique(as.integer(scenario_index$network_replicate))),
    publication_replicates
  ),
  "Module 03 scenario bundles use a stale network-replicate set."
)
check_files(scenario_index$file, "Module 03 scenario bundles")
cat(
  "PASS Module 03 inputs: ",
  expected_observed_scenarios, " observed and ",
  expected_future_scenarios, " future scenario-replicate bundles\n",
  sep = ""
)

# Module 03 results: every publication output must use the current replicate
# count. This prevents an earlier 5- or 10-replicate run from satisfying the
# checker after R/config.R has been changed to 20 replicates.
expected_n_replicates <- length(publication_replicates)
calibration <- read_required_table(
  file.path(
    rvf_config$results_dir, "calibrate_seroprevalence",
    "calibrate_recommended.csv"
  ),
  "Module 03 calibration recommendations"
)
check(
  "n_network_replicates" %in% names(calibration) &&
    all(calibration$n_network_replicates == expected_n_replicates),
  "Module 03 calibration recommendations use a stale replicate count."
)

result_tables <- c(
  file.path(
    rvf_config$results_dir, "mobility_capacity_profile_comparison",
    "mobility_profile_summary.csv"
  ),
  file.path(
    rvf_config$results_dir, "future_cold_endemic",
    "future_cold_endemic_window_burden_aggregate.csv"
  ),
  file.path(
    rvf_config$results_dir,
    "future_cold_endemic_profile_peridomestic_dominant",
    "future_cold_endemic_window_burden_aggregate.csv"
  ),
  file.path(
    rvf_config$results_dir,
    "future_cold_endemic_quadratic_rainfall",
    "future_cold_endemic_window_burden_aggregate.csv"
  ),
  file.path(
    rvf_config$results_dir, "climate_movement_attribution_cold_start",
    "climate_movement_attribution_decomposition_aggregate.csv"
  ),
  file.path(
    rvf_config$results_dir,
    "climate_movement_attribution_cold_start_profile_peridomestic_dominant",
    "climate_movement_attribution_decomposition_aggregate.csv"
  )
)
result_labels <- c(
  "mobility profiles",
  "future grazing/floodplain profile",
  "future village-associated profile",
  "future quadratic-rainfall sensitivity",
  "attribution grazing/floodplain profile",
  "attribution village-associated profile"
)
for (i in seq_along(result_tables)) {
  result <- read_required_table(
    result_tables[[i]], paste("Module 03", result_labels[[i]])
  )
  check(
    "n_network_replicates" %in% names(result) &&
      all(result$n_network_replicates == expected_n_replicates),
    paste0(
      "Module 03 ", result_labels[[i]],
      " output uses a stale replicate count."
    )
  )
  if ("window_length_years" %in% names(result)) {
    check(
      all(result$window_length_years == 40L),
      paste0(
        "Module 03 ", result_labels[[i]],
        " output does not use the common 40-year reporting period."
      )
    )
  }
}

# Observed baselines are stored in profile-specific folders. Require the three
# configured profile/rainfall combinations and the complete replicate set.
baseline_files <- list.files(
  file.path(rvf_config$results_dir, "observed_cycle_endemic"),
  pattern = "_annual_system\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
check(length(baseline_files) >= length(rvf_config$execution$baseline_runs),
      "Module 03 observed-baseline outputs are incomplete.")
baseline_metadata <- do.call(rbind, lapply(baseline_files, function(path) {
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  unique(x[c(
    "site_profile", "rainfall_model", "capacity_scale",
    "network_replicate"
  )])
}))
for (run in rvf_config$execution$baseline_runs) {
  recommended_scale <- calibration$capacity_scale[
    calibration$site_profile == run$site_profile &
      calibration$rainfall_model == run$rainfall_model
  ]
  check(
    length(recommended_scale) == 1L && is.finite(recommended_scale),
    paste0(
      "Calibration does not define one capacity scale for ",
      run$site_profile, " / ", run$rainfall_model, "."
    )
  )
  run_rows <- baseline_metadata[
    baseline_metadata$site_profile == run$site_profile &
      baseline_metadata$rainfall_model == run$rainfall_model &
      abs(baseline_metadata$capacity_scale - recommended_scale) < 1e-12,
  ]
  check(
    setequal(
      sort(unique(as.integer(run_rows$network_replicate))),
      sort(unique(as.integer(run$network_replicates)))
    ),
    paste0(
      "Observed baseline is incomplete for ",
      run$site_profile, " / ", run$rainfall_model, "."
    )
  )
}

# Spatial and seasonal products combine four GCMs with the configured network
# replicates in future windows, and only network replicates under current
# conditions. Their n_runs fields provide a compact stale-output check.
expected_future_runs_per_cell <-
  length(future_climate$reporting_gcms) * expected_n_replicates
spatial_summary <- read_required_table(
  file.path(
    rvf_config$results_dir, "spatial_patterns",
    "spatial_concentration_summary.csv"
  ),
  "Module 03 spatial summary"
)
check(
  all(spatial_summary$n_runs[spatial_summary$scenario == "current"] ==
        expected_n_replicates) &&
    all(spatial_summary$n_runs[spatial_summary$scenario != "current"] ==
          expected_future_runs_per_cell),
  "Module 03 spatial summary uses a stale or incomplete replicate grid."
)
seasonal_summary <- read_required_table(
  file.path(
    rvf_config$results_dir, "seasonal_cycles",
    "seasonal_characteristics_summary.csv"
  ),
  "Module 03 seasonal summary"
)
check(
  all(seasonal_summary$n_runs[seasonal_summary$scenario == "current"] ==
        expected_n_replicates) &&
    all(seasonal_summary$n_runs[seasonal_summary$scenario != "current"] ==
          expected_future_runs_per_cell),
  "Module 03 seasonal summary uses a stale or incomplete replicate grid."
)

# Attribution is additive at the paired replicate level. Checking that table,
# rather than sums of medians, respects the matched factorial design.
attribution_replicate_tables <- c(
  file.path(
    rvf_config$results_dir, "climate_movement_attribution_cold_start",
    "climate_movement_attribution_decomposition_replicates.csv"
  ),
  file.path(
    rvf_config$results_dir,
    "climate_movement_attribution_cold_start_profile_peridomestic_dominant",
    "climate_movement_attribution_decomposition_replicates.csv"
  )
)
for (path in attribution_replicate_tables) {
  attribution_replicates <- read_required_table(
    path,
    "Module 03 replicate-level attribution"
  )
  check(
    all(abs(attribution_replicates$decomposition_identity_error) < 1e-8),
    "Module 03 replicate-level attribution is not additive."
  )
}

cat("PASS Module 03: calibration, scenario results and summaries\n")

cat("\nWhole-pipeline preflight passed.\n")
