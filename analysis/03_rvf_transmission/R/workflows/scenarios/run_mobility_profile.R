#!/usr/bin/env Rscript

#' 04_run_mobility_scenarios.R
#' -------------------------------------------------------------------------
#' Test how the mobile-herd fraction changes endemic potential under one
#' calibrated exposure profile. Every mobility value is initialised
#' independently, seeded with one infectious cow in every village, and run
#' under the same repeated observed climate and movement sequence.
#'
#' The experiment changes residence-time exposure, not vector habitat
#' capacity: carrying capacity is evaluated at the p_mobile = 1 reference for
#' every run.
#'
#' A short workflow test can be enabled under `execution$smoke` in R/config.R.
#' Smoke outputs are isolated and are not suitable for scientific inference.
#' -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source(here::here("analysis", "03_rvf_transmission", "R", "transmission_model.R"))
source(here::here("analysis", "03_rvf_transmission", "R", "cycle_diagnostics.R"))

cfg <- load_rvf_config()

# -------- Analysis specification -----------------------------------------
sim <- cfg$sim
if (exists("analysis_profile", inherits = FALSE)) {
  sim$site_profile <- analysis_profile
}
if (isTRUE(cfg$execution$use_calibrated_profiles)) {
  sim <- rvf_apply_profile_calibration(
    sim,
    cfg = cfg,
    require_target_band = isTRUE(cfg$execution$require_calibration_target)
  )
}

smoke_mode <- isTRUE(cfg$execution$smoke$mobility)
forcing_years <- sim$baseline_years %||% sim$baseline_year
forcing_cycle_years <- length(forcing_years)

if (!is.finite(forcing_cycle_years) || forcing_cycle_years < 1L) {
  stop("The observed forcing cycle must contain at least one year.")
}
if (!isTRUE(all.equal(as.numeric(sim$external_introduction_per_year), 0))) {
  stop("The canonical mobility experiment requires zero external introductions.")
}

# Full analysis: five positive mobility fractions and the configured network
# replicates. Smoke mode reduces both dimensions without changing the main
# publication settings.
# Smoke test: two mobility fractions, two replicates, and four forcing cycles.
p_mobile_values <- if (smoke_mode) c(0.1, 1) else c(0.1, 0.25, 0.5, 0.75, 1)
reference_p_mobile <- 1
# Four observed 10-year cycles give the same 40-year reporting duration used
# for future scenarios, whose forcing cycles are 20 years long.
reporting_years <- as.integer(cfg$execution$reporting_years)
if (reporting_years %% forcing_cycle_years != 0L) {
  stop("The reporting duration must contain complete observed forcing cycles.")
}
terminal_cycles <- if (smoke_mode) {
  2L
} else {
  reporting_years %/% forcing_cycle_years
}
window_cycles <- c(2L, 4L, 6L)
years_burn <- sim$years_burn
years_endemic <- if (smoke_mode) 4L * forcing_cycle_years else sim$years_endemic
output_step <- sim$output_step
seed_per_village <- 1
seed_village <- "all"
required_replicate_fraction <- 0.80
relative_tolerance <- 0.10
seroprevalence_tolerance_pp <- 0.50

if (years_endemic %% forcing_cycle_years != 0L) {
  stop("years_endemic must be a whole number of observed forcing cycles.")
}
if (years_endemic < 2L * terminal_cycles * forcing_cycle_years) {
  stop(
    "The run must contain two non-overlapping reporting blocks of ",
    terminal_cycles, " complete forcing cycles each."
  )
}
if (any(!is.finite(p_mobile_values)) ||
    any(p_mobile_values <= 0) ||
    any(p_mobile_values > 1)) {
  stop("Headline mobility values must be greater than zero and no greater than one.")
}
if (!any(abs(p_mobile_values - reference_p_mobile) < 1e-12)) {
  stop("The mobility grid must include p_mobile = 1 as the paired reference.")
}

scenario_index <- rvf_select_scenarios(
  dataset = "observed",
  year = forcing_years,
  network_replicate = cfg$analysis_network_replicates
)
if (nrow(scenario_index) == 0L) {
  stop(
    "No observed scenario bundles were found for ",
    paste(forcing_years, collapse = ", "),
    ". Run 01_prepare_inputs.R first."
  )
}

network_replicates <- rvf_complete_replicates_for_years(
  index = scenario_index,
  dataset = "observed",
  years = forcing_years,
  network_replicate = cfg$analysis_network_replicates
)
if (length(network_replicates) == 0L) {
  stop("No complete network replicate spans the observed forcing cycle.")
}
network_replicates <- sort(unique(network_replicates))
if (smoke_mode) network_replicates <- head(network_replicates, 2L)
expected_replicates <- if (smoke_mode) {
  head(sort(unique(cfg$analysis_network_replicates)), 2L)
} else {
  sort(unique(cfg$analysis_network_replicates))
}
if (!identical(network_replicates, expected_replicates)) {
  stop(
    "The mobility input contains an incomplete or stale network-replicate set."
  )
}

calibrated <- rvf_adjust_calibrated_params(
  transmission_scale = sim$transmission_scale,
  capacity_scale = sim$capacity_scale,
  site_profile = sim$site_profile,
  rainfall_model = sim$rainfall_model
)

out_name <- if (smoke_mode) {
  paste0("mobility_capacity_smoke_", sim$site_profile)
} else if (identical(sim$network_mode, "monthly")) {
  paste0("mobility_capacity_cold_start_", sim$site_profile, "_pmobile_capacity")
} else {
  paste0(
    "mobility_capacity_cold_start_annualized_",
    sim$site_profile,
    "_pmobile_capacity"
  )
}
out_dir <- if (smoke_mode) {
  file.path(cfg$cache_dir, "smoke", out_name)
} else {
  file.path(cfg$results_dir, out_name)
}
intermediate_dir <- if (smoke_mode) {
  out_dir
} else {
  file.path(cfg$intermediate_dir, out_name)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

scenario_base_ids <- sprintf("OBSERVED_CURRENT_%d", forcing_years)
scenario_label <- sprintf(
  "OBSERVED_%d_%d_CYCLE",
  min(forcing_years),
  max(forcing_years)
)

safe_median <- function(x) {
  if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}

safe_p10 <- function(x) {
  if (all(is.na(x))) {
    NA_real_
  } else {
    as.numeric(stats::quantile(x, 0.10, na.rm = TRUE))
  }
}

safe_p90 <- function(x) {
  if (all(is.na(x))) {
    NA_real_
  } else {
    as.numeric(stats::quantile(x, 0.90, na.rm = TRUE))
  }
}

safe_ratio <- function(value, reference) {
  ifelse(
    is.finite(reference) & abs(reference) > .Machine$double.eps,
    value / reference,
    NA_real_
  )
}

# Function: attach the same run metadata to annual, within-year, and spatial
# summaries without repeating a long mutate block.
append_run_metadata <- function(data, metadata) {
  if (nrow(data) == 0L) return(data)
  metadata_columns <- setdiff(names(metadata), names(data))
  dplyr::bind_cols(
    data,
    metadata[rep(1L, nrow(data)), metadata_columns, drop = FALSE]
  )
}

# Function: define identifiers and scientific assumptions for one run.
run_metadata <- function(p_mobile, network_replicate) {
  shares <- rvf_effective_profile_share_wide(
    site_profile = calibrated$site_profile,
    p_mobile = p_mobile,
    cfg = cfg
  )
  scenario_ids <- rvf_with_replicate(scenario_base_ids, network_replicate)

  tibble(
    source = "observed_current_cycle",
    scenario_design = "independent_endemic_potential_scenario",
    scenario = scenario_label,
    scenario_id = paste(scenario_ids, collapse = ";"),
    scenario_base_id = paste(scenario_base_ids, collapse = ";"),
    network_replicate = network_replicate,
    network_mode = sim$network_mode,
    p_mobile = p_mobile,
    reference_p_mobile = reference_p_mobile,
    years_burn = years_burn,
    years_endemic = years_endemic,
    seed_mode = "all_villages",
    seed_per_village = seed_per_village,
    external_introduction_per_year = sim$external_introduction_per_year,
    rainfall_model = sim$rainfall_model,
    rainfall_label = rvf_rainfall_model_label(sim$rainfall_model, cfg = cfg),
    transmission_scale = sim$transmission_scale,
    capacity_scale = sim$capacity_scale,
    site_profile = calibrated$site_profile,
    site_profile_label = calibrated$site_profile_label,
    m_max_v = calibrated$m_max_v,
    m_max_g = calibrated$m_max_g,
    m_max_w = calibrated$m_max_w,
    effective_share_village_routine =
      shares[["effective_share_village_routine"]],
    effective_share_grazing_routine =
      shares[["effective_share_grazing_routine"]],
    effective_share_watering_routine =
      shares[["effective_share_watering_routine"]]
  )
}

# Function: simulate one paired mobility x network-replicate scenario.
run_one <- function(p_mobile, network_replicate) {
  metadata <- run_metadata(p_mobile, network_replicate)
  scenario_ids <- rvf_with_replicate(scenario_base_ids, network_replicate)
  scenario <- build_rvf_sequence_scenario(
    scenario_ids,
    cfg = cfg,
    validate = TRUE,
    scenario_type = "observed_current_cycle",
    sequence_label = scenario_label
  )
  scenario <- rvf_apply_network_mode(scenario, network_mode = sim$network_mode)

  params <- rvf_assemble_params(
    scenario = scenario,
    calibrated = calibrated,
    use_temp_dependence = TRUE,
    p_mobile = p_mobile,
    # Freeze habitat capacity at the p = 1 host denominator. This prevents a
    # change in host residence from silently changing vector carrying capacity.
    capacity_reference_p_mobile = reference_p_mobile,
    waning_years = sim$waning_years,
    external_introduction_per_year = sim$external_introduction_per_year,
    verbose = FALSE
  )

  result <- rvf_run_simulation(
    params = params,
    years_burn = years_burn,
    years_endemic = years_endemic,
    output_step = output_step,
    n_seed = seed_per_village,
    seed_village = seed_village,
    init_state = NULL,
    verbose = FALSE
  )
  annual <- rvf_annual_metrics(
    result$output_endemic,
    params,
    t_seed = result$t_seed
  )

  annual_system <- append_run_metadata(annual$system_year, metadata)
  time_series <- result$time_series |>
    select(
      year_since_seed,
      step_label,
      seroprevalence,
      prevalence,
      infectious_hosts
    ) |>
    append_run_metadata(metadata)

  # Keep only a compact village summary for the same final reporting block used
  # in the system-level comparison. Full weekly trajectories are not saved.
  terminal_start_year <- years_endemic - terminal_cycles * forcing_cycle_years + 1L
  village_terminal <- annual$village_year |>
    filter(.data$analysis_year >= terminal_start_year) |>
    group_by(.data$village_id) |>
    summarise(
      annual_new_infections = mean(.data$annual_new_infections, na.rm = TRUE),
      annual_mean_prevalence = mean(.data$annual_mean_prevalence, na.rm = TRUE),
      annual_mean_daily_foi = mean(.data$annual_mean_daily_foi, na.rm = TRUE),
      .groups = "drop"
    ) |>
    append_run_metadata(metadata)

  list(
    status = metadata |>
      mutate(run_status = "ok", error_message = NA_character_),
    annual_system = annual_system,
    time_series = time_series,
    village_terminal = village_terminal
  )
}

jobs <- tidyr::crossing(
  p_mobile = p_mobile_values,
  network_replicate = network_replicates
)
n_cores <- rvf_max_cores(
  default = if (smoke_mode) 2L else cfg$execution$max_cores,
  n_tasks = nrow(jobs)
)

message(sprintf(
  "Running %d mobility scenario(s) for %s on %d core(s).",
  nrow(jobs),
  sim$site_profile,
  n_cores
))

results <- rvf_parallel_lapply(seq_len(nrow(jobs)), function(i) {
  job <- jobs[i, ]
  message(sprintf(
    "Mobility scenario %d/%d: p=%.2f, replicate %03d",
    i,
    nrow(jobs),
    job$p_mobile,
    job$network_replicate
  ))

  tryCatch(
    run_one(job$p_mobile, job$network_replicate),
    error = function(error) {
      metadata <- run_metadata(job$p_mobile, job$network_replicate)
      list(
        status = metadata |>
          mutate(run_status = "failed", error_message = conditionMessage(error)),
        annual_system = tibble(),
        time_series = tibble(),
        village_terminal = tibble()
      )
    }
  )
}, n_cores = n_cores)

bind_result <- function(name) {
  bind_rows(lapply(results, function(result) result[[name]]))
}

run_status <- bind_result("status")
annual_system <- bind_result("annual_system")
time_series <- bind_result("time_series")
village_terminal <- bind_result("village_terminal")

rds_path <- file.path(intermediate_dir, "mobility_capacity_runs.rds")
if (nrow(annual_system) == 0L || nrow(time_series) == 0L) {
  saveRDS(
    list(
      status = run_status,
      settings = list(sim = sim, smoke_mode = smoke_mode, jobs = jobs)
    ),
    rds_path
  )
  stop("All mobility scenarios failed. Inspect status in ", rds_path)
}

# -------- Phase-aware recurrent-state diagnostics -----------------------
diagnostic_groups <- c(
  "p_mobile",
  "network_replicate",
  "site_profile",
  "rainfall_model",
  "capacity_scale"
)
cycle_diagnostics <- rvf_build_cycle_diagnostics(
  annual_system = annual_system,
  time_series = time_series,
  group_cols = diagnostic_groups,
  forcing_cycle_years = forcing_cycle_years,
  years_endemic = years_endemic,
  terminal_cycles = terminal_cycles,
  relative_tolerance = relative_tolerance,
  seroprevalence_tolerance_pp = seroprevalence_tolerance_pp,
  required_replicate_fraction = required_replicate_fraction,
  window_cycles = window_cycles
)

join_keys <- diagnostic_groups
terminal_by_run <- cycle_diagnostics$window_sensitivity |>
  filter(.data$summary_cycles == terminal_cycles) |>
  transmute(
    across(all_of(join_keys)),
    terminal_summary_cycles = .data$summary_cycles,
    terminal_summary_years = .data$summary_years,
    annual_new_infections_terminal_window = .data$annual_new_infections,
    annual_mean_prevalence_terminal_window = .data$annual_mean_prevalence,
    annual_mean_daily_foi_terminal_window = .data$annual_mean_daily_foi,
    mean_seroprevalence_terminal_window = .data$mean_seroprevalence_end
  )

final_cycle_by_run <- cycle_diagnostics$cycle_summary |>
  filter(.data$forcing_cycle == cycle_diagnostics$protocol$complete_cycles) |>
  transmute(
    across(all_of(join_keys)),
    annual_new_infections_terminal_cycle = .data$annual_new_infections
  )

final_state_by_run <- time_series |>
  group_by(across(all_of(join_keys))) |>
  slice_max(.data$year_since_seed, n = 1L, with_ties = FALSE) |>
  ungroup() |>
  transmute(
    across(all_of(join_keys)),
    final_seroprevalence = .data$seroprevalence,
    final_prevalence = .data$prevalence,
    final_infectious_hosts = .data$infectious_hosts
  )

summary_by_run <- run_status |>
  filter(.data$run_status == "ok") |>
  left_join(terminal_by_run, by = join_keys) |>
  left_join(final_cycle_by_run, by = join_keys) |>
  left_join(final_state_by_run, by = join_keys) |>
  left_join(cycle_diagnostics$convergence_by_run, by = join_keys) |>
  left_join(cycle_diagnostics$periodicity_by_run, by = join_keys) |>
  left_join(cycle_diagnostics$persistence_by_run, by = join_keys) |>
  mutate(
    persistence_class = case_when(
      .data$fraction_time_below_1e_6_infectious_hosts > 0 ~
        "near_deterministic_extinction",
      .data$fraction_time_below_one_infectious_host > 0 ~
        "fractional_deterministic_persistence",
      TRUE ~ "persistent_above_one_host"
    )
  )

# Pair every mobility result with p = 1 from the same network replicate.
pair_keys <- c(
  "network_replicate",
  "site_profile",
  "rainfall_model",
  "capacity_scale"
)
reference_by_run <- summary_by_run |>
  filter(abs(.data$p_mobile - reference_p_mobile) < 1e-12) |>
  select(
    all_of(pair_keys),
    reference_annual_new_infections = annual_new_infections_terminal_window,
    reference_mean_seroprevalence = mean_seroprevalence_terminal_window,
    reference_final_seroprevalence = final_seroprevalence
  )

summary_by_run <- summary_by_run |>
  left_join(reference_by_run, by = pair_keys) |>
  mutate(
    relative_annual_new_infections_terminal_window = safe_ratio(
      .data$annual_new_infections_terminal_window,
      .data$reference_annual_new_infections
    ),
    relative_mean_seroprevalence_terminal_window = safe_ratio(
      .data$mean_seroprevalence_terminal_window,
      .data$reference_mean_seroprevalence
    ),
    relative_final_seroprevalence = safe_ratio(
      .data$final_seroprevalence,
      .data$reference_final_seroprevalence
    )
  )

# Compare the final 2-cycle estimate with longer 4- and 6-cycle averages.
window_reference <- cycle_diagnostics$window_sensitivity |>
  filter(.data$summary_cycles == terminal_cycles) |>
  select(
    all_of(join_keys),
    reference_window_infections = annual_new_infections,
    reference_window_seroprevalence = mean_seroprevalence_end
  )
window_sensitivity_by_run <- cycle_diagnostics$window_sensitivity |>
  filter(.data$summary_cycles > terminal_cycles) |>
  left_join(window_reference, by = join_keys) |>
  mutate(
    incidence_relative_difference = abs(
      rvf_symmetric_relative_change(
        .data$annual_new_infections,
        .data$reference_window_infections
      )
    ),
    seroprevalence_difference_pp = abs(
      .data$mean_seroprevalence_end -
        .data$reference_window_seroprevalence
    )
  ) |>
  group_by(across(all_of(join_keys))) |>
  summarise(
    window_sensitivity_incidence_max_relative_difference =
      max(.data$incidence_relative_difference, na.rm = TRUE),
    window_sensitivity_seroprevalence_max_difference_pp =
      max(.data$seroprevalence_difference_pp, na.rm = TRUE),
    .groups = "drop"
  )

# -------- Replicate aggregation and interpretation gates ----------------
status_by_cell <- run_status |>
  group_by(
    .data$p_mobile,
    .data$site_profile,
    .data$rainfall_model,
    .data$capacity_scale
  ) |>
  summarise(
    n_planned = n(),
    n_successful = sum(.data$run_status == "ok"),
    success_fraction = .data$n_successful / .data$n_planned,
    .groups = "drop"
  )

aggregate_groups <- c(
  "source",
  "scenario_design",
  "network_mode",
  "p_mobile",
  "reference_p_mobile",
  "seed_mode",
  "seed_per_village",
  "years_burn",
  "years_endemic",
  "terminal_summary_years",
  "rainfall_model",
  "rainfall_label",
  "transmission_scale",
  "capacity_scale",
  "site_profile",
  "site_profile_label",
  "m_max_v",
  "m_max_g",
  "m_max_w",
  "effective_share_village_routine",
  "effective_share_grazing_routine",
  "effective_share_watering_routine"
)

aggregate <- summary_by_run |>
  group_by(across(all_of(aggregate_groups))) |>
  summarise(
    n_network_replicates = n(),
    mean_seroprevalence_terminal_window_mid =
      safe_median(.data$mean_seroprevalence_terminal_window),
    mean_seroprevalence_terminal_window_lo =
      safe_p10(.data$mean_seroprevalence_terminal_window),
    mean_seroprevalence_terminal_window_hi =
      safe_p90(.data$mean_seroprevalence_terminal_window),
    annual_new_infections_terminal_window_mid =
      safe_median(.data$annual_new_infections_terminal_window),
    annual_new_infections_terminal_window_lo =
      safe_p10(.data$annual_new_infections_terminal_window),
    annual_new_infections_terminal_window_hi =
      safe_p90(.data$annual_new_infections_terminal_window),
    annual_new_infections_terminal_cycle_mid =
      safe_median(.data$annual_new_infections_terminal_cycle),
    annual_new_infections_terminal_cycle_lo =
      safe_p10(.data$annual_new_infections_terminal_cycle),
    annual_new_infections_terminal_cycle_hi =
      safe_p90(.data$annual_new_infections_terminal_cycle),
    final_seroprevalence_mid = safe_median(.data$final_seroprevalence),
    final_seroprevalence_lo = safe_p10(.data$final_seroprevalence),
    final_seroprevalence_hi = safe_p90(.data$final_seroprevalence),
    relative_mean_seroprevalence_terminal_window_mid =
      safe_median(.data$relative_mean_seroprevalence_terminal_window),
    relative_mean_seroprevalence_terminal_window_lo =
      safe_p10(.data$relative_mean_seroprevalence_terminal_window),
    relative_mean_seroprevalence_terminal_window_hi =
      safe_p90(.data$relative_mean_seroprevalence_terminal_window),
    relative_annual_new_infections_terminal_window_mid =
      safe_median(.data$relative_annual_new_infections_terminal_window),
    relative_annual_new_infections_terminal_window_lo =
      safe_p10(.data$relative_annual_new_infections_terminal_window),
    relative_annual_new_infections_terminal_window_hi =
      safe_p90(.data$relative_annual_new_infections_terminal_window),
    minimum_infectious_hosts_mid =
      safe_median(.data$minimum_infectious_hosts),
    persistent_above_one_host_fraction = mean(
      .data$persistence_class == "persistent_above_one_host",
      na.rm = TRUE
    ),
    deterministic_persistence_supported_fraction = mean(
      .data$persistence_class != "near_deterministic_extinction",
      na.rm = TRUE
    ),
    period_one_fraction = mean(
      .data$recurrent_pattern == "period_one",
      na.rm = TRUE
    ),
    period_two_fraction = mean(
      .data$recurrent_pattern == "period_two",
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  left_join(
    status_by_cell,
    by = c("p_mobile", "site_profile", "rainfall_model", "capacity_scale")
  )

change_summary <- cycle_diagnostics$convergence_by_run |>
  group_by(
    .data$p_mobile,
    .data$site_profile,
    .data$rainfall_model,
    .data$capacity_scale
  ) |>
  summarise(
    annual_new_infections_window_abs_sym_rel_change_mid =
      safe_median(.data$incidence_abs_symmetric_relative_change),
    annual_new_infections_window_abs_sym_rel_change_lo =
      safe_p10(.data$incidence_abs_symmetric_relative_change),
    annual_new_infections_window_abs_sym_rel_change_hi =
      safe_p90(.data$incidence_abs_symmetric_relative_change),
    .groups = "drop"
  )

window_sensitivity_summary <- window_sensitivity_by_run |>
  group_by(
    .data$p_mobile,
    .data$site_profile,
    .data$rainfall_model,
    .data$capacity_scale
  ) |>
  summarise(
    window_sensitivity_incidence_max_relative_difference_mid = safe_median(
      .data$window_sensitivity_incidence_max_relative_difference
    ),
    window_sensitivity_incidence_max_relative_difference_hi = safe_p90(
      .data$window_sensitivity_incidence_max_relative_difference
    ),
    window_sensitivity_seroprevalence_max_difference_pp_mid = safe_median(
      .data$window_sensitivity_seroprevalence_max_difference_pp
    ),
    window_sensitivity_seroprevalence_max_difference_pp_hi = safe_p90(
      .data$window_sensitivity_seroprevalence_max_difference_pp
    ),
    .groups = "drop"
  )

convergence <- cycle_diagnostics$convergence_summary |>
  left_join(
    change_summary,
    by = c("p_mobile", "site_profile", "rainfall_model", "capacity_scale")
  ) |>
  left_join(
    window_sensitivity_summary,
    by = c("p_mobile", "site_profile", "rainfall_model", "capacity_scale")
  ) |>
  left_join(
    aggregate |>
      select(
        p_mobile,
        site_profile,
        rainfall_model,
        capacity_scale,
        source,
        scenario_design,
        network_mode,
        reference_p_mobile,
        seed_mode,
        seed_per_village,
        transmission_scale,
        site_profile_label,
        n_planned,
        n_successful,
        success_fraction,
        persistent_above_one_host_fraction,
        deterministic_persistence_supported_fraction
      ),
    by = c("p_mobile", "site_profile", "rainfall_model", "capacity_scale")
  ) |>
  mutate(
    convergence_required_fraction = .data$required_replicate_fraction,
    incidence_quasi_equilibrium_supported =
      .data$incidence_stable_fraction >= .data$required_replicate_fraction,
    prevalence_quasi_equilibrium_supported =
      .data$prevalence_stable_fraction >= .data$required_replicate_fraction,
    foi_quasi_equilibrium_supported =
      .data$foi_stable_fraction >= .data$required_replicate_fraction,
    quasi_equilibrium_supported = .data$recurrent_state_supported,
    deterministic_persistence_supported =
      .data$deterministic_persistence_supported_fraction >=
        .data$required_replicate_fraction,
    run_complete = .data$success_fraction >= .data$required_replicate_fraction,
    headline_supported =
      .data$run_complete &
        .data$quasi_equilibrium_supported &
        .data$deterministic_persistence_supported,
    outcome_class = case_when(
      !.data$run_complete ~ "incomplete_run_cell",
      # Relative convergence measures are ill-conditioned when incidence in
      # both comparison windows is effectively zero. Classify those cells by
      # their persistence diagnostic before applying the recurrent-state gate.
      !.data$deterministic_persistence_supported ~
        "near_deterministic_extinction",
      !.data$quasi_equilibrium_supported ~ "finite_horizon_not_settled",
      .data$persistent_above_one_host_fraction <
        .data$required_replicate_fraction ~
        "recurrent_fractional_persistence",
      TRUE ~ "recurrent_persistent_above_one_host"
    ),
    # Compact aliases used by the mobility-summary script.
    annual_new_infections_abs_sym_rel_change_mid =
      .data$annual_new_infections_window_abs_sym_rel_change_mid,
    annual_new_infections_abs_sym_rel_change_lo =
      .data$annual_new_infections_window_abs_sym_rel_change_lo,
    annual_new_infections_abs_sym_rel_change_hi =
      .data$annual_new_infections_window_abs_sym_rel_change_hi
  ) |>
  arrange(.data$p_mobile)

aggregate <- aggregate |>
  left_join(
    convergence |>
      select(
        p_mobile,
        site_profile,
        rainfall_model,
        capacity_scale,
        quasi_equilibrium_supported,
        deterministic_persistence_supported,
        headline_supported,
        outcome_class
      ),
    by = c("p_mobile", "site_profile", "rainfall_model", "capacity_scale")
  ) |>
  arrange(.data$p_mobile)

# Compact cycle trajectories support script 05's diagnostic plot without
# retaining weekly ODE output.
cycle_trajectory_by_run <- annual_system |>
  mutate(
    cycle_index = (.data$analysis_year - 1L) %/% forcing_cycle_years + 1L
  ) |>
  group_by(
    .data$p_mobile,
    .data$network_replicate,
    .data$cycle_index,
    .data$source,
    .data$scenario_design,
    .data$network_mode,
    .data$reference_p_mobile,
    .data$seed_mode,
    .data$seed_per_village,
    .data$rainfall_model,
    .data$transmission_scale,
    .data$capacity_scale,
    .data$site_profile,
    .data$site_profile_label
  ) |>
  summarise(
    first_analysis_year = min(.data$analysis_year),
    final_analysis_year = max(.data$analysis_year),
    annual_new_infections = mean(.data$annual_new_infections, na.rm = TRUE),
    annual_infectious_host_days = mean(
      .data$annual_infectious_host_days,
      na.rm = TRUE
    ),
    annual_mean_prevalence = mean(.data$annual_mean_prevalence, na.rm = TRUE),
    annual_mean_daily_foi = mean(.data$annual_mean_daily_foi, na.rm = TRUE),
    .groups = "drop"
  )

cycle_trajectory <- cycle_trajectory_by_run |>
  group_by(
    .data$source,
    .data$scenario_design,
    .data$network_mode,
    .data$p_mobile,
    .data$reference_p_mobile,
    .data$seed_mode,
    .data$seed_per_village,
    .data$rainfall_model,
    .data$transmission_scale,
    .data$capacity_scale,
    .data$site_profile,
    .data$site_profile_label,
    .data$cycle_index,
    .data$first_analysis_year,
    .data$final_analysis_year
  ) |>
  summarise(
    n_network_replicates = n(),
    annual_new_infections_mid = safe_median(.data$annual_new_infections),
    annual_new_infections_lo = safe_p10(.data$annual_new_infections),
    annual_new_infections_hi = safe_p90(.data$annual_new_infections),
    annual_infectious_host_days_mid =
      safe_median(.data$annual_infectious_host_days),
    annual_mean_prevalence_mid = safe_median(.data$annual_mean_prevalence),
    annual_mean_daily_foi_mid = safe_median(.data$annual_mean_daily_foi),
    .groups = "drop"
  ) |>
  arrange(.data$p_mobile, .data$cycle_index)

settings <- list(
  sim = sim,
  smoke_mode = smoke_mode,
  p_mobile_values = p_mobile_values,
  reference_p_mobile = reference_p_mobile,
  network_replicates = network_replicates,
  forcing_years = forcing_years,
  years_burn = years_burn,
  years_endemic = years_endemic,
  terminal_cycles = terminal_cycles,
  terminal_window_years = terminal_cycles * forcing_cycle_years,
  seeding = list(mode = "all_villages", infectious_cattle_per_village = 1)
)

saveRDS(
  list(
    settings = settings,
    status = run_status,
    summary = summary_by_run,
    aggregate = aggregate,
    convergence = convergence,
    cycle_trajectory = cycle_trajectory,
    annual_system = annual_system,
    village_terminal = village_terminal,
    cycle_diagnostics = cycle_diagnostics
  ),
  rds_path
)
write_csv(aggregate, file.path(out_dir, "mobility_capacity_aggregate.csv"))
write_csv(convergence, file.path(out_dir, "mobility_capacity_convergence.csv"))

# Hard execution gate: comparisons are unsafe when a cell has too few
# successful runs or when successful non-reference runs lack their paired p=1
# network replicate. Scientific non-establishment remains a reported result.
incomplete_cells <- status_by_cell |>
  filter(.data$success_fraction < required_replicate_fraction)
successful_replicates <- run_status |>
  filter(.data$run_status == "ok") |>
  distinct(.data$network_replicate)
reference_replicates <- run_status |>
  filter(
    .data$run_status == "ok",
    abs(.data$p_mobile - reference_p_mobile) < 1e-12
  ) |>
  distinct(.data$network_replicate)
missing_reference_replicates <- anti_join(
  successful_replicates,
  reference_replicates,
  by = "network_replicate"
)

if (nrow(incomplete_cells) > 0L || nrow(missing_reference_replicates) > 0L) {
  stop(
    "Mobility run-completeness gate failed. Outputs were saved for audit in ",
    out_dir
  )
}

if (!smoke_mode && any(!convergence$headline_supported)) {
  warning(
    "At least one mobility cell is not supported for recurrent endemic-state ",
    "interpretation. Inspect mobility_capacity_convergence.csv before plotting.",
    call. = FALSE
  )
}

message("Mobility counterfactual outputs written to: ", out_dir)
if (smoke_mode) {
  message("Smoke test complete. Do not interpret its finite-horizon estimates.")
}
