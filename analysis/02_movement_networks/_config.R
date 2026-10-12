# Module 02 configuration.
#
# These are scientific workflow settings, not command-line options. Edit this
# file when changing the analysis design so that every run leaves an auditable
# record of the assumptions used.

source(here::here("_project_paths.R"))
source(here::here("analysis", "_future_climate_config.R"))

make_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  path
}

module_dir <- proj$mvt_dir
outputs_dir <- make_dir(proj$mvt_outputs)

mvt_config <- list(
  input_schema_version = 3L,

  # Input data
  historical_ndvi_file = file.path(
    proj$ndvi_derived, "gimms_support", "all_historical_observed.feather"
  ),
  future_ndvi_dir = file.path(proj$ndvi_outputs, "scenarios"),
  shape_dir = file.path(proj$shapefiles_dir, "resource_areas"),
  # Sites whose land designation restricts livestock access are identified from
  # their recorded names. The name patterns are supplied with the participatory
  # data and are not part of the public code.
  restricted_site_patterns_file = file.path(
    proj$shapefiles_dir, "resource_areas", "restricted_site_patterns.txt"
  ),
  movement_files = c(
    grazing = file.path(proj$movement_dir, "Local_movt_Gza_tj.csv"),
    watering = file.path(proj$movement_dir, "Local_movt_wp_tj.csv")
  ),

  # A common historical reference and common GCM composition are used for
  # both pathways. Movement parameters are calibrated once to observed data.
  observed_years = 2009:2018,
  scenarios = future_climate$scenarios,
  scenario_labels = future_climate$scenario_labels,
  climate_models = future_climate$reporting_gcms,
  future_windows = future_climate$analysis_windows,
  future_years = sort(unique(unlist(future_climate$analysis_windows))),

  # Output contract
  inputs_dir = make_dir(file.path(outputs_dir, "inputs")),
  scenario_inputs_dir = make_dir(file.path(outputs_dir, "inputs", "scenarios")),
  observed_inputs_file = file.path(outputs_dir, "inputs", "observed_inputs.rds"),
  calibration_dir = make_dir(file.path(outputs_dir, "calibration")),
  validation_dir = make_dir(file.path(outputs_dir, "validation")),
  projection_dir = make_dir(file.path(outputs_dir, "projections")),
  bundles_dir = make_dir(file.path(outputs_dir, "projections", "bundles")),
  projection_job_summaries_dir = make_dir(
    file.path(outputs_dir, "projections", "job_summaries")
  ),
  projection_summaries_dir = make_dir(file.path(outputs_dir, "projections", "summaries")),
  projection_figures_dir = make_dir(file.path(outputs_dir, "projections", "figures")),
  sensitivity_dir = make_dir(file.path(outputs_dir, "sensitivity")),

  # Regime definition. The lower threshold isolates the known August-
  # September peak drought; the median separates routine from intermediate
  # resource conditions. Both thresholds remain fixed in future scenarios.
  drought_percentile = 0.12,
  transition_percentile = 0.50,
  availability_reference_months = c("Aug", "Sep"),
  availability_target_unavailable = list(grazing = 0.15, watering = 0.225),
  availability_floor = 0.02,
  availability_slope = NA_real_,

  # Sensitivity to the uncertain watering-site response. The grazing deficit
  # is held at the current working value, whereas watering spans an equal
  # response, the current assumption, and a stronger dry-season response.
  # These are mean availability deficits, not literal proportions of sites
  # that dry completely.
  availability_sensitivity = list(
    grazing_deficit = 0.15,
    watering_deficits = c(0.15, 0.225, 0.30),
    labels = c("Equal response", "Current assumption", "Stronger water loss"),
    replicates = 100L
  ),

  # Mirror screen for the grazing response, which the watering screen above
  # holds fixed. Participatory-mapping work in Longido district reports that
  # about 40% of rangeland is in use during the extreme dry season and under
  # 20% in a drought year, against the 15% mean deficit assumed here. Longido
  # is a different district, so these are bounds for a structural screen, not
  # values to adopt. Watering is held at the current working value.
  grazing_sensitivity = list(
    watering_deficit = 0.225,
    grazing_deficits = c(0.15, 0.40, 0.60),
    labels = c("Current assumption", "Moderate forage loss", "Severe forage loss"),
    replicates = 100L
  ),

  # Within a regime, lower resource availability increases the expected
  # number of additional choices. The previous month may carry scarcity into
  # the current month's movement response, but destination weights use the
  # current month's site conditions.
  availability_severity_strength = 3,
  availability_severity_cap = 2.5,
  availability_severity_memory = "lag_min",

  # Model assumptions
  use_restriction = TRUE,
  # A small common turnover fraction prevents routine and transition routes
  # from being identical in every month. The value was selected from a bounded
  # screen against annual-union size and is shared across resource types to
  # avoid fitting poorly identified resource-specific turnover parameters.
  monthly_route_turnover_fraction = 0.02,
  drought_tail_turnover_fraction = 0.25,
  fit_year_union_weight = 0.5,

  # Reproducible simulation controls. Six cores is the project-wide ceiling.
  fit_replicates = 20L,
  fit_starts = c(year_round = 6L, transition = 7L, drought = 7L),
  fit_max_iterations = 500L,
  validation_replicates = 100L,
  projection_replicates = 100L,
  n_cores = 6L,
  seeds = list(fit_routine = 42L, fit_transition = 43L,
               fit_drought = 44L, validation = 1000L,
               projection = 2000L)
)

scenario_input_file <- function(scenario, climate_model) {
  file.path(mvt_config$scenario_inputs_dir, scenario,
            paste0(climate_model, ".rds"))
}

network_bundle_file <- function(dataset, year, scenario = NULL,
                                climate_model = NULL) {
  if (identical(dataset, "observed")) {
    return(file.path(mvt_config$bundles_dir, "observed",
                     paste0(year, ".rds")))
  }

  file.path(mvt_config$bundles_dir, scenario, climate_model,
            paste0(year, ".rds"))
}

network_job_summary_file <- function(dataset, year, scenario = NULL,
                                     climate_model = NULL) {
  if (identical(dataset, "observed")) {
    return(file.path(mvt_config$projection_job_summaries_dir, "observed",
                     paste0(year, ".rds")))
  }

  file.path(mvt_config$projection_job_summaries_dir, scenario, climate_model,
            paste0(year, ".rds"))
}
