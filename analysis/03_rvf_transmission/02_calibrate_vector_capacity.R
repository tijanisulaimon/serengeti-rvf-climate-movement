#!/usr/bin/env Rscript

#' 02_calibrate_vector_capacity.R
#' -------------------------------------------------------------------------
#' Literature-anchored endemic plausibility tuning.
#'
#' Compare fixed rainfall-response assumptions across a single tunable
#' vector-abundance scale and identify settings that place endemic cattle
#' seroprevalence in a literature-informed 4-10 % inter-epidemic target band
#' for Tanzania. Lower anchor:
#' de Glanville et al. (2022, PLoS Negl Trop Dis 16:e0010871) reported
#' 4.4 % cattle seroprevalence in northern Tanzania. Upper anchor:
#' Salekwa et al. (2019, Am J Trop Med Hyg 101:459-466) reported 7.7 %
#' overall cattle seroprevalence across selected Tanzanian regions.
#'
#' The host-vector transmission probabilities are held fixed at the
#' literature values in `R/model_parameters.R`
#' (alpha_hv = 0.6, alpha_vh = 0.4). Rainfall-shape parameters and the fixed
#' site-capacity profiles are also assumptions rather than fitted quantities.
#' The only tuned quantity within each rainfall/profile combination is the
#' overall vector-contact scale (`capacity_scale`), because the available
#' evidence is too weak to identify transmission, rainfall shape, site profile,
#' and vector abundance simultaneously.
#'
#' OUTPUTS
#'   - compact tables: `outputs/results/calibrate_seroprevalence/`
#'   - replicate object: `outputs/intermediate/calibrate_seroprevalence/`
#'   - manuscript figure: `outputs/figures/`
#'   - calibrate_grid.csv                    settings evaluated in this run
#'   - calibrate_results.csv                 aggregate over network replicates
#'   - calibrate_results_by_replicate.rds    compact replicate-level rows when complete
#'   - calibrate_results_by_replicate.csv    optional readable copy of those rows
#'   - calibrate_recommended.csv             per-rainfall-family anchor settings
#'   - calibrate_regime_brackets.csv         low / anchor / high endemic settings
#'   - calibrate_effective_profile_shares.csv effective FOI-share targets after residence time
#'   - outputs/figures/01_calibration_profile_gate.png
#' -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

source(here::here("analysis", "03_rvf_transmission", "R", "transmission_model.R"))

cfg <- load_rvf_config()
sim <- cfg$sim
calibration_cfg <- cfg$execution$calibration
sim$years_endemic <- calibration_cfg$years_endemic %||% sim$years_endemic
sim$years_burn <- calibration_cfg$years_burn %||% sim$years_burn
default_long_window_years <- max(1L, sim$years_endemic - (sim$summary_burnin_years %||% 90L))
calibrate_long_window_years <- calibration_cfg$terminal_window_years %||%
  cfg$execution$reporting_years %||%
  default_long_window_years
calibrate_long_window_years <- min(calibrate_long_window_years, sim$years_endemic)
calibration_metric <- "mean_seroprevalence_terminal_window"
calibrate_seed_mode <- calibration_cfg$seed_mode
if (!calibrate_seed_mode %in% c("single_village", "all_villages")) {
  stop("execution$calibration$seed_mode must be 'single_village' or 'all_villages'.")
}
sim$n_seed <- as.integer(calibration_cfg$n_seed)
if (identical(calibrate_seed_mode, "all_villages")) {
  sim$seed_village <- "all"
} else {
  sim$seed_village <- as.integer(calibration_cfg$seed_village)
}
if (!is.finite(sim$years_endemic) || sim$years_endemic < 1L ||
  !is.finite(sim$years_burn) || sim$years_burn < 0L ||
  !is.finite(sim$n_seed) || sim$n_seed < 0L) {
  stop("Invalid calibration horizon or seeding settings in R/config.R.")
}
TARGET_SERO_MIN <- sim$sero_target_min
TARGET_SERO_MAX <- sim$sero_target_max
TARGET_SERO_MID <- sim$sero_target_mid

out_subdir <- "calibrate_seroprevalence"
out_suffix <- trimws(calibration_cfg$output_suffix %||% "")
if (nzchar(out_suffix)) {
  out_suffix <- gsub("[^A-Za-z0-9_.-]+", "_", out_suffix)
  out_subdir <- paste0(out_subdir, "_", out_suffix)
}
out_dir <- file.path(cfg$results_dir, out_subdir)
intermediate_dir <- file.path(cfg$intermediate_dir, out_subdir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

# Reuse a valid aggregate table only when explicitly requested. This supports
# narrow grid refinement without rerunning completed cells. The scientific
# dependency check below prevents reuse after model or scenario changes.
resume_existing_results <- isTRUE(calibration_cfg$resume_existing)
existing_results_path <- file.path(out_dir, "calibrate_results.csv")
existing_results <- if (resume_existing_results && file.exists(existing_results_path)) {
  read_csv(existing_results_path, show_col_types = FALSE)
} else {
  tibble()
}

# Remove only files owned by this script. This prevents a disabled optional
# output, such as the replicate-level table, from surviving from an older run.
calibration_results <- c(
  "calibrate_results.csv",
  "calibrate_results_by_replicate.csv",
  "calibrate_grid.csv",
  "calibrate_recommended.csv",
  "calibrate_regime_brackets.csv",
  "calibrate_effective_profile_shares.csv"
)
if (!resume_existing_results) {
  unlink(file.path(out_dir, calibration_results))
  unlink(file.path(intermediate_dir, "calibrate_results_by_replicate.rds"))
  # Do not leave a completed-run figure visible while a replacement
  # calibration is still running or has failed partway through.
  unlink(file.path(cfg$figures_dir, "01_calibration_profile_gate.png"))
}

baseline_years <- sim$baseline_years %||% sim$baseline_year
scenario_base_ids <- sprintf("OBSERVED_CURRENT_%d", baseline_years)
scenario_index <- rvf_select_scenarios(
  dataset = "observed",
  year = baseline_years,
  network_replicate = cfg$analysis_network_replicates
)
if (nrow(scenario_index) == 0) {
  stop(
    "No observed scenario bundles found for baseline year(s) ", paste(baseline_years, collapse = ", "),
    ". Rebuild replicate-aware scenario bundles first."
  )
}
network_replicates <- rvf_complete_replicates_for_years(
  index = scenario_index,
  dataset = "observed",
  years = baseline_years,
  network_replicate = cfg$analysis_network_replicates
)
if (length(network_replicates) == 0) {
  stop(
    "No complete observed network replicate set found for baseline years: ",
    paste(baseline_years, collapse = ", ")
  )
}
requested_calibration_replicates <- calibration_cfg$replicates
network_replicates <- intersect(network_replicates, sort(unique(requested_calibration_replicates)))
if (length(network_replicates) == 0) {
  stop(
    "Configured calibration replicates did not match any complete observed replicate. ",
    "Available complete replicates: ",
    paste(
      rvf_complete_replicates_for_years(
        index = scenario_index,
        dataset = "observed",
        years = baseline_years,
        network_replicate = cfg$analysis_network_replicates
      ),
      collapse = ", "
    )
  )
}
cat("Calibration network replicates:", paste(network_replicates, collapse = ", "), "\n")

# Calibration checkpoints must track the prepared scenario contents, not only
# their stable IDs and paths. The scenario index can remain byte-identical when
# Module 02 regenerates different network matrices under the same scenario IDs.
calibration_scenario_ids <- scenario_index |>
  filter(network_replicate %in% network_replicates) |>
  pull(scenario_id) |>
  unique() |>
  sort()
calibration_scenario_files <- file.path(
  cfg$scenarios_dir,
  paste0(calibration_scenario_ids, ".rds")
)
missing_calibration_scenarios <- calibration_scenario_files[
  !file.exists(calibration_scenario_files)
]
if (length(missing_calibration_scenarios) > 0L) {
  stop(
    "Calibration scenario files are missing:\n - ",
    paste(missing_calibration_scenarios, collapse = "\n - "),
    "\nRun 01_prepare_inputs.R before calibration."
  )
}

# These files define either the transmission model or the exact climate and
# network inputs used by the calibration. Their hashes form the cache protocol.
scientific_dependency_files <- c(
  here::here(
    "analysis", "03_rvf_transmission",
    "02_calibrate_vector_capacity.R"
  ),
  here::here("analysis", "03_rvf_transmission", "R", "config.R"),
  here::here("analysis", "03_rvf_transmission", "R", "scenario_io.R"),
  here::here("analysis", "03_rvf_transmission", "R", "model_parameters.R"),
  here::here("analysis", "03_rvf_transmission", "R", "transmission_model.R"),
  here::here(
    "analysis", "03_rvf_transmission", "R", "workflows", "preparation",
    "prepare_networks.R"
  ),
  here::here(
    "analysis", "03_rvf_transmission", "R", "workflows", "preparation",
    "extract_site_climate.R"
  ),
  here::here(
    "analysis", "03_rvf_transmission", "R", "workflows", "preparation",
    "assemble_scenario_bundles.R"
  ),
  cfg$mvt_bundle_index,
  file.path(cfg$scenarios_dir, "scenario_bundle_index.csv"),
  calibration_scenario_files
)
scientific_dependency_files <- unique(scientific_dependency_files)
missing_scientific_dependencies <- scientific_dependency_files[
  !file.exists(scientific_dependency_files)
]
if (length(missing_scientific_dependencies) > 0L) {
  stop(
    "Calibration scientific dependencies are missing:\n - ",
    paste(missing_scientific_dependencies, collapse = "\n - ")
  )
}

if (resume_existing_results) {
  if (nrow(existing_results) == 0L) {
    message("No existing aggregate calibration table found; all requested cells will run.")
  } else {
    results_mtime <- file.info(existing_results_path)$mtime
    newer_dependencies <- scientific_dependency_files[
      file.info(scientific_dependency_files)$mtime > results_mtime
    ]
    if (length(newer_dependencies) > 0L) {
      stop(
        "Cannot resume from the aggregate calibration table because scientific ",
        "dependencies are newer than it:\n - ",
        paste(newer_dependencies, collapse = "\n - "),
        "\nSet execution$calibration$resume_existing to FALSE and rebuild the full grid."
      )
    }

    required_resume_columns <- c(
      "years_burn", "years_endemic", "calibration_metric",
      "calibration_window_years", "seed_mode", "n_seed", "p_mobile",
      "n_network_replicates", "transmission_scale", "target_sero_min",
      "target_sero_max"
    )
    missing_resume_columns <- setdiff(required_resume_columns, names(existing_results))
    if (length(missing_resume_columns) > 0L) {
      stop(
        "Existing calibration results lack resume metadata: ",
        paste(missing_resume_columns, collapse = ", "),
        ". Rebuild the full grid."
      )
    }

    expected_protocol <- list(
      years_burn = sim$years_burn,
      years_endemic = sim$years_endemic,
      calibration_metric = calibration_metric,
      calibration_window_years = calibrate_long_window_years,
      seed_mode = calibrate_seed_mode,
      n_seed = sim$n_seed,
      p_mobile = sim$p_mobile,
      n_network_replicates = length(network_replicates),
      transmission_scale = sim$transmission_scale,
      target_sero_min = TARGET_SERO_MIN,
      target_sero_max = TARGET_SERO_MAX
    )
    mismatched_protocol <- names(expected_protocol)[vapply(
      names(expected_protocol),
      function(field) {
        observed <- unique(existing_results[[field]])
        expected <- expected_protocol[[field]]
        if (is.numeric(expected)) {
          length(observed) != 1L || !is.finite(observed) || abs(observed - expected) > 1e-12
        } else {
          length(observed) != 1L || !identical(as.character(observed), as.character(expected))
        }
      },
      logical(1)
    )]
    if (length(mismatched_protocol) > 0L) {
      stop(
        "Existing calibration results use a different protocol: ",
        paste(mismatched_protocol, collapse = ", "),
        ". Rebuild the full grid rather than mixing results."
      )
    }
    message(
      "Validated ", nrow(existing_results),
      " existing aggregate calibration cells for reuse."
    )
  }
}

# Screening ranges differ by rainfall response because the Hill curve has a
# higher mean capacity fraction than the quadratic curve. The centres come
# from the pre-repair exploratory calibration, but the ranges are deliberately
# wider because every value must now be re-evaluated with the corrected
# 2009--2018 climate cycle.
default_capacity_grids <- list(
  saturating_rainfall = list(
    peridomestic_dominant = c(0.12, 0.15, 0.17, 0.20, 0.23),
    balanced = c(0.62, 0.78, 0.90, 1.03, 1.17),
    grazing_floodplain_dominant = c(0.40, 0.46, 0.48, 0.49, 0.50, 0.52, 0.58, 0.64),
    resource_edge_dominant = c(0.82, 1.03, 1.12, 1.13, 1.14, 1.19, 1.35, 1.54),
    waterpoint_extreme = c(2.50, 3.10, 3.60, 4.10, 4.65)
  ),
  quadratic_rainfall = list(
    peridomestic_dominant = c(0.22, 0.27, 0.30, 0.34, 0.38),
    balanced = c(1.40, 1.56, 1.72, 1.90, 2.10),
    grazing_floodplain_dominant = c(0.70, 0.82, 0.88, 0.96, 1.08),
    resource_edge_dominant = c(1.70, 1.90, 2.05, 2.30, 2.60),
    waterpoint_extreme = c(5.00, 5.80, 6.30, 7.10, 8.00)
  )
)

# Return the prespecified capacity-scale grid for one ecological profile.
capacity_grid_for_profile <- function(site_profile, rainfall_model) {
  default_grid <- default_capacity_grids[[rainfall_model]][[site_profile]]
  if (is.null(default_grid)) {
    stop(
      "No default capacity grid for rainfall_model=", rainfall_model,
      " and site_profile=", site_profile, "."
    )
  }
  default_grid
}

# Keep the transmission block and rainfall-shape parameters fixed. Compare
# rainfall assumptions while tuning only the single vector-scale parameter.
rainfall_models_to_run <- calibration_cfg$rainfall_models
unknown_rainfall <- setdiff(rainfall_models_to_run, names(cfg$rainfall_models))
if (length(unknown_rainfall) > 0) {
  stop("Unknown rainfall model(s) in calibration request: ", paste(unknown_rainfall, collapse = ", "))
}

site_profiles_to_run <- calibration_cfg$profiles
site_profiles_to_run <- rvf_resolve_site_profile_ids(site_profiles_to_run, cfg = cfg)
unknown_profiles <- setdiff(site_profiles_to_run, names(cfg$site_profiles))
if (length(unknown_profiles) > 0) {
  stop("Unknown site profile(s) in calibration request: ", paste(unknown_profiles, collapse = ", "))
}

# Stable key used to match requested cells to aggregate results and internal
# checkpoints. Twelve significant digits are ample for the explicit grids.
calibration_cell_id <- function(rainfall_model, site_profile, capacity_scale) {
  paste(
    rainfall_model,
    site_profile,
    formatC(capacity_scale, digits = 12, format = "fg", flag = "#"),
    sep = "__"
  )
}

grid <- bind_rows(lapply(rainfall_models_to_run, function(rainfall_model) {
  bind_rows(lapply(site_profiles_to_run, function(site_profile) {
    tibble(
      site_profile = site_profile,
      rainfall_model = rainfall_model,
      capacity_scale = capacity_grid_for_profile(site_profile, rainfall_model)
    )
  }))
})) |>
  mutate(rainfall_order = match(.data$rainfall_model, rainfall_models_to_run)) |>
  arrange(.data$rainfall_order, .data$site_profile, .data$capacity_scale) |>
  select(-all_of("rainfall_order")) |>
  mutate(
    cell_id = calibration_cell_id(
      .data$rainfall_model,
      .data$site_profile,
      .data$capacity_scale
    )
  )

N_CELLS <- nrow(grid)
if (nrow(existing_results) > 0L) {
  existing_results <- existing_results |>
    mutate(
      cell_id = calibration_cell_id(
        .data$rainfall_model,
        .data$site_profile,
        .data$capacity_scale
      )
    )
}
grid_to_run <- if (nrow(existing_results) > 0L) {
  grid |>
    anti_join(existing_results |> distinct(.data$cell_id), by = "cell_id")
} else {
  grid
}
N_RUN_CELLS <- nrow(grid_to_run)

N_CORES <- rvf_max_cores(
  default = cfg$execution$max_cores,
  n_tasks = max(1L, N_RUN_CELLS)
)

cat(sprintf(
  "Calibration grid: %d cells (%d reused, %d to run), %d site profile(s), %d core(s)\n",
  N_CELLS, N_CELLS - N_RUN_CELLS, N_RUN_CELLS,
  length(site_profiles_to_run), N_CORES
))
grid_summary <- grid |>
  group_by(rainfall_model, site_profile) |>
  summarise(
    capacity_values = paste(format(sort(unique(capacity_scale)), trim = TRUE), collapse = ", "),
    .groups = "drop"
  )
print(grid_summary, n = Inf, width = Inf)
write_csv(grid |> select(-all_of("cell_id")), file.path(out_dir, "calibrate_grid.csv"))

effective_profile_shares <- bind_rows(lapply(site_profiles_to_run, function(site_profile) {
  rvf_effective_profile_shares(
    site_profile = site_profile,
    p_mobile = sim$p_mobile,
    cfg = cfg
  )
}))
write_csv(effective_profile_shares, file.path(out_dir, "calibrate_effective_profile_shares.csv"))

# Cell checkpoints are internal computation caches, not manuscript outputs.
# Each stores its scientific protocol and is ignored automatically if model
# code, scenario inputs, or core assumptions change.
use_checkpoints <- isTRUE(calibration_cfg$use_checkpoints)
checkpoint_dir <- file.path(cfg$cache_dir, "calibration_cells", out_subdir)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
checkpoint_protocol <- list(
  version = 2L,
  scientific_dependency_md5 = as.list(
    tools::md5sum(scientific_dependency_files)
  ),
  baseline_years = as.integer(baseline_years),
  network_replicates = as.integer(network_replicates),
  years_burn = sim$years_burn,
  years_endemic = sim$years_endemic,
  calibration_metric = calibration_metric,
  calibration_window_years = calibrate_long_window_years,
  seed_mode = calibrate_seed_mode,
  seed_village = as.character(sim$seed_village),
  n_seed = sim$n_seed,
  p_mobile = sim$p_mobile,
  # Store the actual residence-time assumptions in the checkpoint protocol.
  time_budget = unlist(rvf_get_time_budget(), use.names = TRUE),
  transmission_scale = sim$transmission_scale,
  external_introduction_per_year = sim$external_introduction_per_year,
  serology_target = c(min = TARGET_SERO_MIN, max = TARGET_SERO_MAX),
  ode_controls = rvf_ode_controls(),
  rainfall_models = cfg$rainfall_models[rainfall_models_to_run],
  site_profiles = cfg$site_profiles[site_profiles_to_run]
)

# Map a stable calibration-cell identifier to its cache file.
checkpoint_path <- function(cell_id) {
  safe_id <- gsub("[^A-Za-z0-9_.-]+", "_", cell_id)
  file.path(checkpoint_dir, paste0(safe_id, ".rds"))
}

# Optional audit mode: materialise the planned grid and profile shares without
# starting any ODE integrations.
if (isTRUE(calibration_cfg$plan_only)) {
  message("Calibration plan written without running the ODE: ", out_dir)
  quit(save = "no", status = 0L, runLast = FALSE)
}

#' Run one calibration cell for one replicate-specific observed-year sequence.
#'
#' @param rep_idx            network replicate used across all baseline years.
#' @param site_profile       named village/grazing/watering m_max profile.
#' @param rainfall_model     rainfall-response family to use.
#' @param capacity_scale     multiplier on m_max_{v,g,w} (K0_* recomputed).
run_one_replicate <- function(rep_idx, site_profile, rainfall_model, capacity_scale) {
  # Calibration is replicate-aware: the same parameter pair is evaluated
  # against the configured observed climate/network cycle for each retained
  # network realisation, then aggregated below.
  scenario_ids <- rvf_with_replicate(scenario_base_ids, rep_idx)
  scenario <- build_rvf_sequence_scenario(
    scenario_ids,
    cfg = cfg,
    validate = TRUE,
    scenario_type = "observed_current_cycle",
    sequence_label = sprintf("OBSERVED_%d_%d_CYCLE", min(baseline_years), max(baseline_years))
  )
  scenario <- rvf_apply_network_mode(scenario, network_mode = sim$network_mode)
  rainfall_spec <- rvf_get_rainfall_model_spec(rainfall_model = rainfall_model, cfg = cfg)
  calibrated <- rvf_adjust_calibrated_params(
    transmission_scale = sim$transmission_scale,
    capacity_scale = capacity_scale,
    site_profile = site_profile,
    rainfall_model = rainfall_model
  )
  share_cols <- rvf_effective_profile_share_wide(
    site_profile = site_profile,
    p_mobile = sim$p_mobile,
    cfg = cfg
  )

  params <- rvf_assemble_params(
    scenario = scenario,
    calibrated = calibrated,
    use_temp_dependence = TRUE,
    p_mobile = sim$p_mobile,
    waning_years = sim$waning_years,
    external_introduction_per_year = sim$external_introduction_per_year,
    verbose = FALSE
  )

  res <- rvf_run_simulation(
    params = params,
    years_burn = sim$years_burn,
    years_endemic = sim$years_endemic,
    output_step = sim$output_step,
    n_seed = sim$n_seed,
    seed_village = sim$seed_village,
    init_state = NULL,
    verbose = FALSE
  )

  ts <- res$time_series
  annual <- rvf_annual_metrics(res$output_endemic, params, t_seed = res$t_seed)
  late <- ts |> filter(year_since_seed >= (sim$years_endemic - 5))
  terminal_window_start <- max(
    0,
    sim$years_endemic - calibrate_long_window_years
  )
  # The ODE output contains regular weekly samples plus exact forcing
  # boundaries. Duration weighting prevents those extra boundaries from
  # receiving more influence in the calibration average.
  late_window <- ts |>
    arrange(.data$year_since_seed) |>
    mutate(
      dt_years = lead(.data$year_since_seed) - .data$year_since_seed
    ) |>
    filter(
      .data$year_since_seed >= terminal_window_start,
      .data$year_since_seed < sim$years_endemic,
      is.finite(.data$dt_years),
      .data$dt_years > 0
    )
  if (nrow(late_window) == 0L) {
    stop("The calibration terminal window contains no weighted time steps.")
  }
  seroprevalence_terminal_mean <- weighted.mean(
    late_window$seroprevalence,
    late_window$dt_years,
    na.rm = TRUE
  )
  seroprevalence_terminal_sd <- sqrt(weighted.mean(
    (late_window$seroprevalence - seroprevalence_terminal_mean)^2,
    late_window$dt_years,
    na.rm = TRUE
  ))
  prevalence_terminal_mean <- weighted.mean(
    late_window$prevalence,
    late_window$dt_years,
    na.rm = TRUE
  )
  final_year <- max(annual$system_year$analysis_year, na.rm = TRUE)
  annual_last5 <- annual$system_year |>
    filter(.data$analysis_year >= (final_year - 4))
  terminal_window_years <- min(calibrate_long_window_years, final_year)
  annual_terminal_window <- annual$system_year |>
    filter(.data$analysis_year >= (final_year - terminal_window_years + 1))
  total_seeded_cattle <- if (identical(sim$seed_village, "all")) {
    sim$n_seed * params$n_vil
  } else {
    sim$n_seed
  }

  tibble(
    scenario_id = paste(scenario_ids, collapse = ";"),
    scenario_base_id = paste(scenario_base_ids, collapse = ";"),
    network_replicate = rep_idx,
    baseline_years = paste(baseline_years, collapse = ","),
    years_burn = sim$years_burn,
    years_endemic = sim$years_endemic,
    calibration_metric = calibration_metric,
    calibration_window_years = terminal_window_years,
    seed_mode = calibrate_seed_mode,
    seed_village = as.character(sim$seed_village),
    n_seed = sim$n_seed,
    total_seeded_cattle = total_seeded_cattle,
    run_status = "ok",
    error_message = NA_character_,
    site_profile = calibrated$site_profile,
    site_profile_label = calibrated$site_profile_label,
    site_profile_ratio_v = calibrated$site_profile_ratio_v,
    site_profile_ratio_g = calibrated$site_profile_ratio_g,
    site_profile_ratio_w = calibrated$site_profile_ratio_w,
    effective_share_village_routine = share_cols$effective_share_village_routine,
    effective_share_grazing_routine = share_cols$effective_share_grazing_routine,
    effective_share_watering_routine = share_cols$effective_share_watering_routine,
    effective_share_village_transition = share_cols$effective_share_village_transition,
    effective_share_grazing_transition = share_cols$effective_share_grazing_transition,
    effective_share_watering_transition = share_cols$effective_share_watering_transition,
    effective_share_village_drought = share_cols$effective_share_village_drought,
    effective_share_grazing_drought = share_cols$effective_share_grazing_drought,
    effective_share_watering_drought = share_cols$effective_share_watering_drought,
    m_max_v = calibrated$m_max_v,
    m_max_g = calibrated$m_max_g,
    m_max_w = calibrated$m_max_w,
    rainfall_model = rainfall_model,
    rainfall_label = rvf_rainfall_model_label(rainfall_model, cfg = cfg),
    p_mobile = sim$p_mobile,
    R_half = rainfall_spec$R_half %||% NA_real_,
    n_hill = rainfall_spec$n_hill %||% NA_real_,
    min_K = rainfall_spec$min_K %||% NA_real_,
    rainfall_lower = rainfall_spec$rainfall_lower %||% NA_real_,
    rainfall_upper = rainfall_spec$rainfall_upper %||% NA_real_,
    transmission_scale = sim$transmission_scale,
    capacity_scale = capacity_scale,
    external_introduction_per_year = sim$external_introduction_per_year,
    final_seroprevalence = tail(ts$seroprevalence, 1),
    mean_seroprevalence_last5y = mean(late$seroprevalence, na.rm = TRUE),
    sd_seroprevalence_last5y = sd(late$seroprevalence, na.rm = TRUE),
    cv_seroprevalence_last5y = sd(late$seroprevalence, na.rm = TRUE) /
      mean(late$seroprevalence, na.rm = TRUE),
    mean_seroprevalence_terminal_window = seroprevalence_terminal_mean,
    sd_seroprevalence_terminal_window = seroprevalence_terminal_sd,
    cv_seroprevalence_terminal_window =
      seroprevalence_terminal_sd / seroprevalence_terminal_mean,
    mean_seroprevalence_calibration = seroprevalence_terminal_mean,
    final_prevalence = tail(ts$prevalence, 1),
    mean_prevalence_last5y = mean(late$prevalence, na.rm = TRUE),
    mean_prevalence_terminal_window = prevalence_terminal_mean,
    final_year = final_year,
    annual_mean_daily_foi_last5y = mean(annual_last5$annual_mean_daily_foi, na.rm = TRUE),
    mean_annual_cum_hazard_last5y = mean(annual_last5$mean_annual_cum_hazard, na.rm = TRUE),
    mean_annual_foi_probability_last5y = mean(annual_last5$mean_annual_foi_probability, na.rm = TRUE),
    annual_new_infections_last5y = mean(annual_last5$annual_new_infections, na.rm = TRUE),
    annual_mean_daily_foi_terminal_window =
      mean(annual_terminal_window$annual_mean_daily_foi, na.rm = TRUE),
    mean_annual_cum_hazard_terminal_window =
      mean(annual_terminal_window$mean_annual_cum_hazard, na.rm = TRUE),
    mean_annual_foi_probability_terminal_window =
      mean(annual_terminal_window$mean_annual_foi_probability, na.rm = TRUE),
    annual_new_infections_terminal_window =
      mean(annual_terminal_window$annual_new_infections, na.rm = TRUE)
  )
}

# Preserve failed runs in the replicate table so failure rates remain visible.
failed_replicate_row <- function(rep_idx, site_profile, rainfall_model, capacity_scale, error_message) {
  scenario_ids <- rvf_with_replicate(scenario_base_ids, rep_idx)
  rainfall_spec <- rvf_get_rainfall_model_spec(rainfall_model = rainfall_model, cfg = cfg)
  calibrated <- rvf_adjust_calibrated_params(
    transmission_scale = sim$transmission_scale,
    capacity_scale = capacity_scale,
    site_profile = site_profile,
    rainfall_model = rainfall_model
  )
  share_cols <- rvf_effective_profile_share_wide(
    site_profile = site_profile,
    p_mobile = sim$p_mobile,
    cfg = cfg
  )

  tibble(
    scenario_id = paste(scenario_ids, collapse = ";"),
    scenario_base_id = paste(scenario_base_ids, collapse = ";"),
    network_replicate = rep_idx,
    baseline_years = paste(baseline_years, collapse = ","),
    years_burn = sim$years_burn,
    years_endemic = sim$years_endemic,
    calibration_metric = calibration_metric,
    calibration_window_years = calibrate_long_window_years,
    seed_mode = calibrate_seed_mode,
    seed_village = as.character(sim$seed_village),
    n_seed = sim$n_seed,
    total_seeded_cattle = NA_real_,
    run_status = "failed",
    error_message = error_message,
    site_profile = calibrated$site_profile,
    site_profile_label = calibrated$site_profile_label,
    site_profile_ratio_v = calibrated$site_profile_ratio_v,
    site_profile_ratio_g = calibrated$site_profile_ratio_g,
    site_profile_ratio_w = calibrated$site_profile_ratio_w,
    effective_share_village_routine = share_cols$effective_share_village_routine,
    effective_share_grazing_routine = share_cols$effective_share_grazing_routine,
    effective_share_watering_routine = share_cols$effective_share_watering_routine,
    effective_share_village_transition = share_cols$effective_share_village_transition,
    effective_share_grazing_transition = share_cols$effective_share_grazing_transition,
    effective_share_watering_transition = share_cols$effective_share_watering_transition,
    effective_share_village_drought = share_cols$effective_share_village_drought,
    effective_share_grazing_drought = share_cols$effective_share_grazing_drought,
    effective_share_watering_drought = share_cols$effective_share_watering_drought,
    m_max_v = calibrated$m_max_v,
    m_max_g = calibrated$m_max_g,
    m_max_w = calibrated$m_max_w,
    rainfall_model = rainfall_model,
    rainfall_label = rvf_rainfall_model_label(rainfall_model, cfg = cfg),
    p_mobile = sim$p_mobile,
    R_half = rainfall_spec$R_half %||% NA_real_,
    n_hill = rainfall_spec$n_hill %||% NA_real_,
    min_K = rainfall_spec$min_K %||% NA_real_,
    rainfall_lower = rainfall_spec$rainfall_lower %||% NA_real_,
    rainfall_upper = rainfall_spec$rainfall_upper %||% NA_real_,
    transmission_scale = sim$transmission_scale,
    capacity_scale = capacity_scale,
    external_introduction_per_year = sim$external_introduction_per_year,
    final_seroprevalence = NA_real_,
    mean_seroprevalence_last5y = NA_real_,
    sd_seroprevalence_last5y = NA_real_,
    cv_seroprevalence_last5y = NA_real_,
    mean_seroprevalence_terminal_window = NA_real_,
    sd_seroprevalence_terminal_window = NA_real_,
    cv_seroprevalence_terminal_window = NA_real_,
    mean_seroprevalence_calibration = NA_real_,
    final_prevalence = NA_real_,
    mean_prevalence_last5y = NA_real_,
    mean_prevalence_terminal_window = NA_real_,
    final_year = NA_integer_,
    annual_mean_daily_foi_last5y = NA_real_,
    mean_annual_cum_hazard_last5y = NA_real_,
    mean_annual_foi_probability_last5y = NA_real_,
    annual_new_infections_last5y = NA_real_,
    annual_mean_daily_foi_terminal_window = NA_real_,
    mean_annual_cum_hazard_terminal_window = NA_real_,
    mean_annual_foi_probability_terminal_window = NA_real_,
    annual_new_infections_terminal_window = NA_real_
  )
}

# Summary helpers return NA, rather than NaN, when every replicate failed.
safe_mean <- function(x) {
  out <- mean(x, na.rm = TRUE)
  if (is.nan(out)) NA_real_ else out
}

safe_sd <- function(x) {
  out <- sd(x, na.rm = TRUE)
  if (is.nan(out)) NA_real_ else out
}

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  as.numeric(quantile(x, probs = probability, names = FALSE, type = 8))
}

# Evaluate one rainfall-profile-scale cell across all network replicates.
run_one_cell <- function(site_profile, rainfall_model, capacity_scale) {
  # Raw replicate-level results are retained within each checkpoint so we can
  # inspect whether a parameter pair depends on an unusual network draw.
  raw <- bind_rows(lapply(network_replicates, function(rep_idx) {
    tryCatch(
      run_one_replicate(rep_idx, site_profile, rainfall_model, capacity_scale),
      error = function(e) {
        msg <- conditionMessage(e)
        message(sprintf(
          "Calibration replicate failed | profile=%s rainfall=%s cap=%.2f rep=%s | %s",
          site_profile, rainfall_model, capacity_scale, rep_idx, msg
        ))
        failed_replicate_row(rep_idx, site_profile, rainfall_model, capacity_scale, msg)
      }
    )
  }))

  # The profile gate and recommendation table use these summaries, which
  # average over network uncertainty at a fixed grid cell.
  summary <- raw |>
    summarise(
      years_burn = dplyr::first(years_burn),
      years_endemic = dplyr::first(years_endemic),
      calibration_metric = dplyr::first(calibration_metric),
      calibration_window_years = dplyr::first(calibration_window_years),
      seed_mode = dplyr::first(seed_mode),
      seed_village = dplyr::first(seed_village),
      n_seed = dplyr::first(n_seed),
      total_seeded_cattle = dplyr::first(total_seeded_cattle),
      site_profile = dplyr::first(site_profile),
      site_profile_label = dplyr::first(site_profile_label),
      site_profile_ratio_v = dplyr::first(site_profile_ratio_v),
      site_profile_ratio_g = dplyr::first(site_profile_ratio_g),
      site_profile_ratio_w = dplyr::first(site_profile_ratio_w),
      effective_share_village_routine = dplyr::first(effective_share_village_routine),
      effective_share_grazing_routine = dplyr::first(effective_share_grazing_routine),
      effective_share_watering_routine = dplyr::first(effective_share_watering_routine),
      effective_share_village_transition = dplyr::first(effective_share_village_transition),
      effective_share_grazing_transition = dplyr::first(effective_share_grazing_transition),
      effective_share_watering_transition = dplyr::first(effective_share_watering_transition),
      effective_share_village_drought = dplyr::first(effective_share_village_drought),
      effective_share_grazing_drought = dplyr::first(effective_share_grazing_drought),
      effective_share_watering_drought = dplyr::first(effective_share_watering_drought),
      m_max_v = dplyr::first(m_max_v),
      m_max_g = dplyr::first(m_max_g),
      m_max_w = dplyr::first(m_max_w),
      rainfall_model = dplyr::first(rainfall_model),
      rainfall_label = dplyr::first(rainfall_label),
      p_mobile = dplyr::first(p_mobile),
      R_half = dplyr::first(R_half),
      n_hill = dplyr::first(n_hill),
      min_K = dplyr::first(min_K),
      rainfall_lower = dplyr::first(rainfall_lower),
      rainfall_upper = dplyr::first(rainfall_upper),
      transmission_scale = sim$transmission_scale,
      capacity_scale = dplyr::first(capacity_scale),
      external_introduction_per_year = sim$external_introduction_per_year,
      n_network_replicates = n(),
      n_success = sum(run_status == "ok", na.rm = TRUE),
      n_failed = sum(run_status != "ok", na.rm = TRUE),
      success_fraction = n_success / n_network_replicates,
      failure_messages = paste(
        unique(error_message[run_status != "ok" & !is.na(error_message)]),
        collapse = " | "
      ),
      between_replicate_sd_seroprevalence = safe_sd(mean_seroprevalence_calibration),
      mean_seroprevalence_calibration_p10 = safe_quantile(mean_seroprevalence_calibration, 0.10),
      mean_seroprevalence_calibration_p90 = safe_quantile(mean_seroprevalence_calibration, 0.90),
      final_seroprevalence = safe_mean(final_seroprevalence),
      mean_seroprevalence_last5y = safe_mean(mean_seroprevalence_last5y),
      sd_seroprevalence_last5y = safe_mean(sd_seroprevalence_last5y),
      cv_seroprevalence_last5y = safe_mean(cv_seroprevalence_last5y),
      mean_seroprevalence_terminal_window = safe_mean(mean_seroprevalence_terminal_window),
      sd_seroprevalence_terminal_window = safe_mean(sd_seroprevalence_terminal_window),
      cv_seroprevalence_terminal_window = safe_mean(cv_seroprevalence_terminal_window),
      mean_seroprevalence_calibration = safe_mean(mean_seroprevalence_calibration),
      final_prevalence = safe_mean(final_prevalence),
      mean_prevalence_last5y = safe_mean(mean_prevalence_last5y),
      mean_prevalence_terminal_window = safe_mean(mean_prevalence_terminal_window),
      annual_mean_daily_foi_last5y = safe_mean(annual_mean_daily_foi_last5y),
      mean_annual_cum_hazard_last5y = safe_mean(mean_annual_cum_hazard_last5y),
      mean_annual_foi_probability_last5y = safe_mean(mean_annual_foi_probability_last5y),
      annual_new_infections_last5y = safe_mean(annual_new_infections_last5y),
      annual_mean_daily_foi_terminal_window = safe_mean(annual_mean_daily_foi_terminal_window),
      mean_annual_cum_hazard_terminal_window = safe_mean(mean_annual_cum_hazard_terminal_window),
      mean_annual_foi_probability_terminal_window =
        safe_mean(mean_annual_foi_probability_terminal_window),
      annual_new_infections_terminal_window = safe_mean(annual_new_infections_terminal_window)
    )

  list(raw = raw, summary = summary)
}

# Load one matching checkpoint or run and atomically save one new cell.
# Reuse only checkpoints generated under the identical scientific protocol.
run_or_load_cell <- function(i) {
  cell <- grid_to_run[i, , drop = FALSE]
  cell_file <- checkpoint_path(cell$cell_id[[1]])
  if (use_checkpoints && file.exists(cell_file)) {
    cached <- tryCatch(readRDS(cell_file), error = function(e) NULL)
    valid_cache <- is.list(cached) &&
      identical(cached$cell_id, cell$cell_id[[1]]) &&
      isTRUE(all.equal(cached$protocol, checkpoint_protocol, check.attributes = TRUE)) &&
      is.list(cached$result) && all(c("raw", "summary") %in% names(cached$result))
    if (valid_cache) {
      cat(sprintf(
        "[%d/%d cached] profile=%s rainfall=%s cap=%.2f\n",
        i, N_RUN_CELLS, cell$site_profile, cell$rainfall_model, cell$capacity_scale
      ))
      flush.console()
      return(cached$result)
    }
  }

  started <- proc.time()[["elapsed"]]
  cat(sprintf(
    "[%d/%d start] profile=%s rainfall=%s cap=%.2f\n",
    i, N_RUN_CELLS, cell$site_profile, cell$rainfall_model, cell$capacity_scale
  ))
  flush.console()
  result <- run_one_cell(
    site_profile = cell$site_profile,
    rainfall_model = cell$rainfall_model,
    capacity_scale = cell$capacity_scale
  )

  if (use_checkpoints) {
    tmp_file <- tempfile("calibration_cell_", tmpdir = checkpoint_dir, fileext = ".rds")
    saveRDS(
      list(
        protocol = checkpoint_protocol,
        cell_id = cell$cell_id[[1]],
        result = result
      ),
      tmp_file,
      compress = "gzip"
    )
    if (file.exists(cell_file)) unlink(cell_file)
    if (!file.rename(tmp_file, cell_file)) {
      unlink(tmp_file)
      stop("Could not write calibration checkpoint: ", cell_file)
    }
  }

  elapsed_minutes <- (proc.time()[["elapsed"]] - started) / 60
  cat(sprintf(
    "[%d/%d done] profile=%s rainfall=%s cap=%.2f elapsed=%.1f min\n",
    i, N_RUN_CELLS, cell$site_profile, cell$rainfall_model,
    cell$capacity_scale, elapsed_minutes
  ))
  flush.console()
  result
}

# Modest parallelism is safe here because grid cells are independent.
cell_results <- if (N_RUN_CELLS > 0L) {
  rvf_parallel_lapply(seq_len(N_RUN_CELLS), run_or_load_cell, n_cores = N_CORES)
} else {
  message("All requested calibration cells were reused; no ODE integrations are needed.")
  list()
}

results_raw <- bind_rows(lapply(cell_results, `[[`, "raw"))
results_new <- bind_rows(lapply(cell_results, `[[`, "summary"))
results <- bind_rows(
  results_new,
  existing_results |> select(-any_of("cell_id"))
) |>
  mutate(
    cell_id = calibration_cell_id(
      .data$rainfall_model,
      .data$site_profile,
      .data$capacity_scale
    )
  ) |>
  distinct(.data$cell_id, .keep_all = TRUE) |>
  filter(.data$cell_id %in% grid$cell_id) |>
  arrange(.data$rainfall_model, .data$site_profile, .data$capacity_scale) |>
  select(-all_of("cell_id"))
if (nrow(results) != N_CELLS) {
  stop(
    "Calibration assembly produced ", nrow(results), " aggregate cells for a ",
    N_CELLS, "-cell grid. Check resume and checkpoint inputs before proceeding."
  )
}

# Tag cells by where they sit relative to the broad literature target band.
results <- results |>
  mutate(
    target_sero_min = TARGET_SERO_MIN,
    target_sero_max = TARGET_SERO_MAX,
    target_sero_mid = TARGET_SERO_MID,
    target_band = !is.na(mean_seroprevalence_calibration) &
      mean_seroprevalence_calibration >= TARGET_SERO_MIN &
      mean_seroprevalence_calibration <= TARGET_SERO_MAX,
    target_distance = if_else(
      is.na(mean_seroprevalence_calibration),
      NA_real_,
      pmax(
        0,
        TARGET_SERO_MIN - mean_seroprevalence_calibration,
        mean_seroprevalence_calibration - TARGET_SERO_MAX
      )
    ),
    target_mid_distance = if_else(
      is.na(mean_seroprevalence_calibration),
      NA_real_,
      abs(mean_seroprevalence_calibration - TARGET_SERO_MID)
    ),
    endemic_regime = case_when(
      is.na(mean_seroprevalence_calibration) ~ "failed",
      mean_seroprevalence_calibration < TARGET_SERO_MIN ~ "low",
      mean_seroprevalence_calibration > TARGET_SERO_MAX ~ "high",
      TRUE ~ "anchor"
    ),
    profile_decision = case_when(
      endemic_regime == "failed" ~ "failed_solver_or_input",
      target_band & success_fraction < 0.8 ~ "unstable_anchor_too_many_failed_replicates",
      target_band ~ "main_or_sensitivity_anchor",
      endemic_regime == "low" ~ "rejected_low_under_current_time_budget",
      endemic_regime == "high" ~ "rejected_high_under_current_time_budget",
      TRUE ~ "unclassified"
    ),
    regime_reference_seroprevalence = case_when(
      endemic_regime == "failed" ~ NA_real_,
      endemic_regime == "low" ~ TARGET_SERO_MIN,
      endemic_regime == "high" ~ TARGET_SERO_MAX,
      TRUE ~ TARGET_SERO_MID
    ),
    regime_distance = abs(mean_seroprevalence_calibration - regime_reference_seroprevalence)
  ) |>
  arrange(site_profile, rainfall_model, target_distance, target_mid_distance, mean_seroprevalence_calibration)

write_csv(results, file.path(out_dir, "calibrate_results.csv"))
raw_results_complete <- nrow(results_raw) == N_CELLS * length(network_replicates)
if (raw_results_complete) {
  saveRDS(
    results_raw,
    file.path(intermediate_dir, "calibrate_results_by_replicate.rds"),
    compress = "gzip"
  )
} else {
  unlink(file.path(intermediate_dir, "calibrate_results_by_replicate.rds"))
  message(
    "Replicate-level RDS not written because aggregate resume supplied ",
    N_CELLS - N_RUN_CELLS, " cell(s) without raw rows."
  )
}
if (rvf_write_derived_csvs() && raw_results_complete) {
  write_csv(results_raw, file.path(out_dir, "calibrate_results_by_replicate.csv"))
} else {
  unlink(file.path(out_dir, "calibrate_results_by_replicate.csv"))
  message(
    "Skipping replicate-level calibration CSV; a complete compact RDS is ",
    "written only for full/checkpoint-complete runs."
  )
}

# Recommended cells are selected within each rainfall family. A candidate must
# have finite output and at least 80% successful network-replicate runs before
# it can be used downstream.
eligible_recommendations <- results |>
  filter(
    is.finite(.data$mean_seroprevalence_calibration),
    .data$success_fraction >= 0.80,
    !isTRUE(cfg$execution$require_calibration_target) | .data$target_band
  )
recommendation_coverage <- eligible_recommendations |>
  count(.data$site_profile, .data$rainfall_model)
expected_recommendation_groups <- tidyr::crossing(
  site_profile = site_profiles_to_run,
  rainfall_model = rainfall_models_to_run
)
missing_recommendation_groups <- anti_join(
  expected_recommendation_groups,
  recommendation_coverage,
  by = c("site_profile", "rainfall_model")
)
if (nrow(missing_recommendation_groups) > 0L) {
  stop(
    "No eligible calibration cell for: ",
    paste(
      paste(
        missing_recommendation_groups$site_profile,
        missing_recommendation_groups$rainfall_model,
        sep = " / "
      ),
      collapse = ", "
    )
  )
}

# Selecting separately avoids treating the broad serology range as evidence
# that one rainfall response family is more likely than the other.
recommended <- eligible_recommendations |>
  group_by(site_profile, site_profile_label, rainfall_model, rainfall_label) |>
  arrange(
    target_distance,
    target_mid_distance,
    between_replicate_sd_seroprevalence,
    cv_seroprevalence_terminal_window,
    .by_group = TRUE
  ) |>
  slice_head(n = 1) |>
  ungroup()

regime_brackets <- results |>
  group_by(site_profile, site_profile_label, rainfall_model, rainfall_label, endemic_regime) |>
  arrange(regime_distance, between_replicate_sd_seroprevalence, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  arrange(site_profile, rainfall_model, factor(endemic_regime, levels = c("low", "anchor", "high")))

write_csv(recommended, file.path(out_dir, "calibrate_recommended.csv"))
write_csv(regime_brackets, file.path(out_dir, "calibrate_regime_brackets.csv"))

# Profile gate plot: the weak serology band is an anchor, not an objective
# function. Showing raw points against the band is clearer than a heatmap when
# each profile uses its own narrow capacity-scale bracket.
profile_levels <- recommended |>
  distinct(.data$site_profile_label) |>
  arrange(.data$site_profile_label) |>
  pull(.data$site_profile_label)
rainfall_levels <- vapply(
  rainfall_models_to_run,
  function(model) cfg$rainfall_models[[model]]$label,
  character(1)
)

p <- ggplot(
  results,
  aes(
    x = .data$capacity_scale,
    y = .data$mean_seroprevalence_calibration
  )
) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = TARGET_SERO_MIN, ymax = TARGET_SERO_MAX,
    fill = "grey75", alpha = 0.28
  ) +
  geom_line(colour = "grey55", linewidth = 0.4) +
  geom_errorbar(
    aes(
      ymin = .data$mean_seroprevalence_calibration_p10,
      ymax = .data$mean_seroprevalence_calibration_p90
    ),
    width = 0,
    linewidth = 0.35,
    colour = "grey35"
  ) +
  geom_point(aes(fill = .data$endemic_regime), shape = 21, size = 2.6, colour = "grey20") +
  geom_point(
    data = recommended,
    shape = 21, size = 4.2, stroke = 1.1,
    fill = NA, colour = "black"
  ) +
  facet_grid(
    rows = vars(factor(site_profile_label, levels = profile_levels)),
    cols = vars(factor(rainfall_label, levels = rainfall_levels)),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    breaks = scales::breaks_pretty(n = 6),
    expand = expansion(mult = c(0, 0.06))
  ) +
  scale_fill_manual(
    values = c(anchor = "#74add1", high = "#d73027", low = "#4575b4", failed = "grey70"),
    breaks = c("low", "anchor", "high", "failed"),
    labels = c(low = "Below anchor", anchor = "Within 4-10%", high = "Above anchor", failed = "Failed")
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", colour = "grey50"),
    legend.position = "bottom",
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    legend.text = element_text(colour = "black")
  ) +
  labs(
    x = "Capacity scale",
    y = "Mean cattle seroprevalence (%)",
    fill = NULL
  )
figures_dir <- cfg$figures_dir
write_headline_figures <- isTRUE(calibration_cfg$write_headline_figures) &&
  !nzchar(out_suffix)
if (write_headline_figures) {
  ggsave(
    file.path(figures_dir, "01_calibration_profile_gate.png"),
    p,
    width = 10.5,
    height = 10.5,
    dpi = 300
  )
} else {
  message(
    "Not updating the shared manuscript calibration figure for this run."
  )
}

print(results)
message("Calibration complete. Results written to: ", out_dir)
