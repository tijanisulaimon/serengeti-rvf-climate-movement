#!/usr/bin/env Rscript

#' 07_summarise_future_scenarios.R
#' -------------------------------------------------------------------------
#' Compare the observed recurrent baseline with the future endemic-potential
#' scenarios produced by script 06. Observed and future runs are each averaged
#' over a common 40-year reporting period (four observed cycles and two future
#' cycles), then annualised before comparison.
#'
#' Absolute burden is summarised within each SSP x GCM x future window. Observed
#' and future movement draws were generated independently, so descriptive fold
#' changes use the median observed burden as a common denominator. Raw burdens
#' and their network-replicate quantiles remain the primary results.
#' -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

source(here::here("analysis", "03_rvf_transmission", "R", "transmission_model.R"))
source(here::here("analysis", "03_rvf_transmission", "R", "plot_helpers.R"))

cfg <- load_rvf_config()
sim <- cfg$sim
if (exists("analysis_profile", inherits = FALSE)) {
  sim$site_profile <- analysis_profile
}
if (exists("analysis_rainfall_model", inherits = FALSE)) {
  sim$rainfall_model <- analysis_rainfall_model
}
primary_future_profile <- rvf_resolve_site_profile_id(
  sim$future_site_profile %||% "grazing_floodplain_dominant"
)
sim$site_profile <- rvf_resolve_site_profile_id(
  sim$site_profile %||% primary_future_profile
)
is_primary_profile <- identical(sim$site_profile, primary_future_profile)
is_primary_rainfall <- identical(
  sim$rainfall_model,
  cfg$primary_rainfall_model
)
sim <- rvf_apply_profile_calibration(
  sim,
  cfg = cfg,
  require_target_band = TRUE
)

smoke_mode <- isTRUE(cfg$execution$smoke$future)
future_windows <- cfg$future_windows
future_scenarios <- cfg$future_scenarios
climate_models <- cfg$climate_models
network_replicates <- if (exists("analysis_network_replicates", inherits = FALSE)) {
  sort(unique(analysis_network_replicates))
} else {
  sort(unique(cfg$analysis_network_replicates))
}
if (smoke_mode) {
  future_windows <- future_windows[1]
  future_scenarios <- future_scenarios[1]
  climate_models <- climate_models[1]
  network_replicates <- head(network_replicates, 2L)
}

baseline_years <- sim$baseline_years %||% sim$baseline_year
baseline_cycle_years <- length(baseline_years)
future_cycle_years <- unique(vapply(future_windows, length, integer(1)))
if (length(future_cycle_years) != 1L) {
  stop("All future forcing windows must contain the same number of years.")
}

# Compare mean annual burden over the same 40-year duration. Both periods
# contain complete forcing cycles: four observed 10-year cycles and two future
# 20-year cycles. Using complete cycles avoids dependence on starting phase.
reporting_years <- as.integer(cfg$execution$reporting_years)
if (reporting_years %% baseline_cycle_years != 0L ||
    reporting_years %% future_cycle_years != 0L) {
  stop("The reporting period must contain complete observed and future cycles.")
}
baseline_reporting_cycles <- reporting_years %/% baseline_cycle_years
future_reporting_cycles <- reporting_years %/% future_cycle_years

future_name <- if (smoke_mode) {
  "future_endemic_smoke"
} else if (identical(sim$network_mode, "monthly")) {
  "future_cold_endemic"
} else {
  "future_cold_endemic_annualized"
}
if (!is_primary_profile) {
  profile_tag <- gsub("[^A-Za-z0-9]+", "_", sim$site_profile)
  future_name <- paste0(future_name, "_profile_", profile_tag)
}
future_name <- rvf_rainfall_output_subdir(
  future_name,
  sim$rainfall_model,
  cfg
)
future_dir <- if (smoke_mode) {
  file.path(cfg$cache_dir, "smoke", future_name)
} else {
  file.path(cfg$results_dir, future_name)
}
future_intermediate_dir <- if (smoke_mode) {
  future_dir
} else {
  file.path(cfg$intermediate_dir, future_name)
}
figures_dir <- if (smoke_mode) future_dir else cfg$figures_dir
dir.create(future_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(future_intermediate_dir, recursive = TRUE, showWarnings = FALSE)

future_paths <- c(
  status = file.path(future_dir, "future_cold_endemic_run_status.csv"),
  summary = file.path(future_dir, "future_cold_endemic_summary.csv"),
  convergence = file.path(future_dir, "future_cold_endemic_convergence.csv"),
  annual = file.path(future_intermediate_dir, "future_cold_endemic_annual_system.rds"),
  cycles = file.path(future_intermediate_dir, "future_cold_endemic_cycle_diagnostics.rds")
)
missing_inputs <- future_paths[!file.exists(future_paths)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing script-06 output(s):\n  ",
    paste(missing_inputs, collapse = "\n  "),
    "\nRun 06_run_future_scenarios.R first."
  )
}

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

# Function: annualise epidemiological outcomes over complete terminal cycles.
# annual_df: annual system output for one network replicate.
# window_years: final complete years retained from each run.
summarise_annual_window <- function(annual_df, window_years, total_hosts) {
  final_year <- max(annual_df$analysis_year, na.rm = TRUE)
  window_start <- final_year - window_years + 1L
  window_df <- annual_df |>
    filter(
      .data$analysis_year >= window_start,
      .data$analysis_year <= final_year
    )

  observed_years <- n_distinct(window_df$analysis_year)
  if (observed_years != window_years) {
    stop(
      "Terminal reporting window contains ", observed_years,
      " years; expected ", window_years, "."
    )
  }

  tibble(
    terminal_window_start_analysis_year = window_start,
    terminal_window_end_analysis_year = final_year,
    window_length_years = observed_years,
    window_total_new_infections = sum(
      window_df$annual_new_infections,
      na.rm = TRUE
    ),
    annual_new_infections = mean(
      window_df$annual_new_infections,
      na.rm = TRUE
    ),
    annual_new_infections_per_1000 =
      mean(window_df$annual_new_infections, na.rm = TRUE) /
        total_hosts * 1000,
    annual_mean_infectious_hosts = mean(
      window_df$annual_mean_infectious_hosts,
      na.rm = TRUE
    ),
    annual_mean_prevalence = mean(
      window_df$annual_mean_prevalence,
      na.rm = TRUE
    ),
    mean_annual_foi_probability = mean(
      window_df$mean_annual_foi_probability,
      na.rm = TRUE
    ),
    mean_annual_cum_hazard = mean(
      window_df$mean_annual_cum_hazard,
      na.rm = TRUE
    ),
    annual_mean_daily_foi = mean(
      window_df$annual_mean_daily_foi,
      na.rm = TRUE
    )
  )
}

# -------- Validate the run grid before summarising -----------------------
expected_jobs <- tidyr::crossing(
  scenario = future_scenarios,
  climate_model = climate_models,
  future_window = names(future_windows),
  network_replicate = network_replicates
)
future_status <- read_csv(future_paths[["status"]], show_col_types = FALSE)
future_status <- future_status |>
  mutate(future_window = as.character(.data$future_window))
status_keys <- c(
  "scenario",
  "climate_model",
  "future_window",
  "network_replicate"
)

duplicate_status <- future_status |>
  count(across(all_of(status_keys)), name = "n") |>
  filter(.data$n != 1L)
missing_status <- anti_join(
  expected_jobs,
  future_status |> distinct(across(all_of(status_keys))),
  by = status_keys
)
unexpected_status <- anti_join(
  future_status |> distinct(across(all_of(status_keys))),
  expected_jobs,
  by = status_keys
)
failed_status <- future_status |>
  filter(.data$run_status != "ok")
if (nrow(duplicate_status) > 0L || nrow(missing_status) > 0L ||
    nrow(unexpected_status) > 0L || nrow(failed_status) > 0L) {
  stop(
    "Future run-completeness gate failed: duplicates=", nrow(duplicate_status),
    ", missing=", nrow(missing_status),
    ", unexpected=", nrow(unexpected_status),
    ", failed=", nrow(failed_status),
    ". Inspect ", future_paths[["status"]]
  )
}

future_convergence <- read_csv(
  future_paths[["convergence"]],
  show_col_types = FALSE
) |>
  mutate(future_window = as.character(.data$future_window))
future_summary <- read_csv(
  future_paths[["summary"]],
  show_col_types = FALSE
) |>
  mutate(future_window = as.character(.data$future_window))

expected_cells <- expected_jobs |>
  distinct(.data$scenario, .data$climate_model, .data$future_window)
convergence_cells <- future_convergence |>
  distinct(.data$scenario, .data$climate_model, .data$future_window)
missing_convergence <- anti_join(
  expected_cells,
  convergence_cells,
  by = c("scenario", "climate_model", "future_window")
)
if (nrow(missing_convergence) > 0L) {
  stop("Future convergence diagnostics do not cover the complete run grid.")
}

# The strict diagnostic requires year-by-year phase matching for incidence,
# prevalence, FOI, and seroprevalence. A bounded multi-cycle pattern can fail
# that test even when the cycle-averaged quantities reported here are stable.
# Reporting support therefore also requires agreement between the final 2-,
# 4-, and 6-cycle summaries; this does not change either tolerance.
cell_keys <- c("scenario", "climate_model", "future_window")
window_support <- future_summary |>
  left_join(
    future_convergence |>
      select(
        all_of(cell_keys),
        relative_tolerance,
        seroprevalence_tolerance_pp
      ),
    by = cell_keys,
    relationship = "many-to-one"
  ) |>
  group_by(across(all_of(cell_keys))) |>
  summarise(
    incidence_window_stable_fraction = mean(
      .data$window_sensitivity_incidence_max_relative_difference <=
        .data$relative_tolerance,
      na.rm = TRUE
    ),
    seroprevalence_window_stable_fraction = mean(
      .data$window_sensitivity_seroprevalence_max_difference_pp <=
        .data$seroprevalence_tolerance_pp,
      na.rm = TRUE
    ),
    incidence_window_max_relative_difference = max(
      .data$window_sensitivity_incidence_max_relative_difference,
      na.rm = TRUE
    ),
    seroprevalence_window_max_difference_pp = max(
      .data$window_sensitivity_seroprevalence_max_difference_pp,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

reporting_support <- future_convergence |>
  left_join(window_support, by = cell_keys, relationship = "one-to-one") |>
  mutate(
    phase_matched_recurrence_supported = .data$recurrent_state_supported,
    cycle_averaged_burden_supported =
      .data$run_complete &
      .data$incidence_stable_fraction >=
        .data$required_replicate_fraction &
      .data$prevalence_stable_fraction >=
        .data$required_replicate_fraction &
      .data$foi_stable_fraction >=
        .data$required_replicate_fraction &
      .data$incidence_window_stable_fraction >=
        .data$required_replicate_fraction &
      .data$seroprevalence_window_stable_fraction >=
        .data$required_replicate_fraction,
    reporting_support_class = case_when(
      .data$phase_matched_recurrence_supported ~
        "phase_matched_recurrent",
      .data$cycle_averaged_burden_supported ~
        "phase_sensitive_bounded",
      TRUE ~ "finite_horizon_not_settled"
    )
  )
write_csv(
  reporting_support,
  file.path(future_dir, "future_cold_endemic_reporting_support.csv")
)

unsupported_reporting <- reporting_support |>
  filter(!.data$cycle_averaged_burden_supported)
if (!smoke_mode && nrow(unsupported_reporting) > 0L) {
  stop(
    nrow(unsupported_reporting),
    " future cell(s) lack stable cycle-averaged burden. Extend or diagnose ",
    "script 06 before plotting."
  )
}
phase_sensitive_cells <- reporting_support |>
  filter(.data$reporting_support_class == "phase_sensitive_bounded")
if (nrow(phase_sensitive_cells) > 0L) {
  message(
    nrow(phase_sensitive_cells),
    " future cell(s) have stable cycle-averaged burden but not strict ",
    "year-by-year phase matching; see future_cold_endemic_reporting_support.csv."
  )
}

if (!all(future_summary$site_profile == sim$site_profile)) {
  stop("Script-06 outputs use a different exposure profile from script 07.")
}
if (!all(future_summary$rainfall_model == sim$rainfall_model)) {
  stop("Script-06 outputs use a different rainfall response from script 07.")
}
if (!all(abs(future_summary$p_mobile - 1) < 1e-12)) {
  stop("Script-06 outputs do not use the canonical p_mobile = 1 setting.")
}

# -------- Observed recurrent baseline ------------------------------------
total_hosts <- rvf_baseline_total_hosts(sim, cfg)
baseline_annual_path <- rvf_baseline_annual_system_path(sim, cfg)
baseline_cycle_path <- file.path(
  rvf_baseline_intermediate_dir(sim, cfg),
  paste0(rvf_baseline_run_tag(sim), "_cycle_diagnostics.rds")
)
if (!file.exists(baseline_annual_path) || !file.exists(baseline_cycle_path)) {
  stop(
    "Missing observed recurrent-baseline outputs. Run ",
    "03_run_observed_baseline.R first."
  )
}

baseline_annual_source <- read_csv(
  baseline_annual_path,
  show_col_types = FALSE
) |>
  filter(.data$network_replicate %in% network_replicates)
baseline_diagnostics <- readRDS(baseline_cycle_path)

missing_baseline_replicates <- setdiff(
  network_replicates,
  unique(baseline_annual_source$network_replicate)
)
if (length(missing_baseline_replicates) > 0L) {
  stop(
    "Observed baseline is missing network replicate(s): ",
    paste(missing_baseline_replicates, collapse = ", ")
  )
}
if (!isTRUE(all(baseline_diagnostics$convergence_summary$recurrent_state_supported))) {
  stop("Observed baseline lacks recurrent-state support.")
}

baseline_serology <- baseline_diagnostics$window_sensitivity |>
  filter(
    .data$summary_cycles == baseline_reporting_cycles,
    .data$network_replicate %in% network_replicates
  ) |>
  transmute(
    network_replicate = .data$network_replicate,
    mean_seroprevalence = .data$mean_seroprevalence_end
  )
baseline_recurrence <- baseline_diagnostics$convergence_by_run |>
  filter(.data$network_replicate %in% network_replicates) |>
  select(network_replicate, recurrent_state_supported)

baseline_metrics <- baseline_annual_source |>
  group_by(.data$network_replicate) |>
  group_modify(function(.x, .y) {
    summarise_annual_window(
      annual_df = .x,
      window_years = reporting_years,
      total_hosts = total_hosts
    )
  }) |>
  ungroup() |>
  left_join(baseline_serology, by = "network_replicate") |>
  left_join(baseline_recurrence, by = "network_replicate") |>
  mutate(
    scenario_type = "current_baseline",
    scenario = "observed",
    scenario_label = "Observed baseline",
    climate_model = "observed_current",
    future_year = NA_integer_,
    future_window = "current",
    window_start = min(baseline_years),
    window_end = max(baseline_years),
    window_years = paste(baseline_years, collapse = ","),
    network_mode = sim$network_mode,
    external_introduction_per_year = 0,
    rainfall_model = sim$rainfall_model,
    rainfall_label = rvf_rainfall_model_label(sim$rainfall_model, cfg = cfg),
    p_mobile = sim$p_mobile,
    transmission_scale = sim$transmission_scale,
    capacity_scale = sim$capacity_scale,
    site_profile = sim$site_profile,
    site_profile_label = rvf_site_profile_label(sim$site_profile, cfg = cfg),
    annual_new_infections_ratio_vs_baseline = 1,
    annual_mean_infectious_hosts_ratio_vs_baseline = 1,
    mean_annual_foi_probability_ratio_vs_baseline = 1
  )

# -------- Future terminal-cycle summaries --------------------------------
future_annual <- readRDS(future_paths[["annual"]])
future_group_keys <- c(
  "scenario",
  "scenario_label",
  "climate_model",
  "future_year",
  "future_window",
  "window_start",
  "window_end",
  "window_years",
  "network_replicate",
  "network_mode",
  "external_introduction_per_year",
  "rainfall_model",
  "rainfall_label",
  "p_mobile",
  "transmission_scale",
  "capacity_scale",
  "site_profile",
  "site_profile_label"
)
future_metrics <- future_annual |>
  filter(
    .data$scenario %in% future_scenarios,
    .data$climate_model %in% climate_models,
    .data$future_window %in% names(future_windows),
    .data$network_replicate %in% network_replicates
  ) |>
  group_by(across(all_of(future_group_keys))) |>
  group_modify(function(.x, .y) {
    summarise_annual_window(
      annual_df = .x,
      window_years = reporting_years,
      total_hosts = total_hosts
    )
  }) |>
  ungroup() |>
  mutate(scenario_type = "future_window")

future_diagnostics <- readRDS(future_paths[["cycles"]])
diagnostic_keys <- c(
  "scenario",
  "climate_model",
  "future_window",
  "network_replicate",
  "network_mode",
  "external_introduction_per_year",
  "rainfall_model",
  "site_profile"
)
future_serology <- future_diagnostics$window_sensitivity |>
  filter(.data$summary_cycles == future_reporting_cycles) |>
  transmute(
    across(all_of(diagnostic_keys)),
    mean_seroprevalence = .data$mean_seroprevalence_end
  )
future_recurrence <- future_diagnostics$convergence_by_run |>
  select(all_of(diagnostic_keys), recurrent_state_supported)
future_metrics <- future_metrics |>
  left_join(future_serology, by = diagnostic_keys) |>
  left_join(future_recurrence, by = diagnostic_keys)

if (nrow(future_metrics) != nrow(expected_jobs)) {
  stop(
    "Future terminal summaries contain ", nrow(future_metrics),
    " runs; expected ", nrow(expected_jobs), "."
  )
}

# Observed and future movement draws do not share a stochastic seed, so equal
# replicate numbers are labels rather than matched pairs. Use the median across
# observed network replicates as the common denominator for descriptive ratios.
# Raw burdens and their network-replicate quantiles remain the primary results.
baseline_reference <- baseline_metrics |>
  summarise(
    baseline_annual_new_infections = safe_median(
      .data$annual_new_infections
    ),
    baseline_annual_mean_infectious_hosts = safe_median(
      .data$annual_mean_infectious_hosts
    ),
    baseline_mean_annual_foi_probability = safe_median(
      .data$mean_annual_foi_probability
    )
  )
future_metrics <- future_metrics |>
  mutate(
    baseline_annual_new_infections =
      baseline_reference$baseline_annual_new_infections[[1]],
    baseline_annual_mean_infectious_hosts =
      baseline_reference$baseline_annual_mean_infectious_hosts[[1]],
    baseline_mean_annual_foi_probability =
      baseline_reference$baseline_mean_annual_foi_probability[[1]]
  ) |>
  mutate(
    annual_new_infections_ratio_vs_baseline = safe_ratio(
      .data$annual_new_infections,
      .data$baseline_annual_new_infections
    ),
    annual_mean_infectious_hosts_ratio_vs_baseline = safe_ratio(
      .data$annual_mean_infectious_hosts,
      .data$baseline_annual_mean_infectious_hosts
    ),
    mean_annual_foi_probability_ratio_vs_baseline = safe_ratio(
      .data$mean_annual_foi_probability,
      .data$baseline_mean_annual_foi_probability
    )
  )

replicate_output_columns <- c(
  "scenario_type", "scenario", "scenario_label", "climate_model",
  "future_year", "future_window", "window_start", "window_end",
  "window_years", "network_replicate", "network_mode",
  "external_introduction_per_year", "rainfall_model", "rainfall_label",
  "p_mobile", "transmission_scale", "capacity_scale", "site_profile",
  "site_profile_label", "terminal_window_start_analysis_year",
  "terminal_window_end_analysis_year", "window_length_years",
  "window_total_new_infections", "annual_new_infections",
  "annual_new_infections_per_1000", "annual_mean_infectious_hosts",
  "annual_mean_prevalence", "mean_annual_foi_probability",
  "mean_annual_cum_hazard", "annual_mean_daily_foi",
  "mean_seroprevalence", "annual_new_infections_ratio_vs_baseline",
  "annual_mean_infectious_hosts_ratio_vs_baseline",
  "mean_annual_foi_probability_ratio_vs_baseline",
  "recurrent_state_supported"
)
replicate_metrics <- bind_rows(baseline_metrics, future_metrics) |>
  select(all_of(replicate_output_columns))

# -------- Network-replicate and GCM summaries ----------------------------
aggregate_groups <- c(
  "scenario_type",
  "scenario",
  "scenario_label",
  "climate_model",
  "future_year",
  "future_window",
  "window_start",
  "window_end",
  "window_years",
  "network_mode",
  "external_introduction_per_year",
  "rainfall_model",
  "rainfall_label",
  "p_mobile",
  "transmission_scale",
  "capacity_scale",
  "site_profile",
  "site_profile_label",
  "window_length_years"
)
aggregate_metrics <- replicate_metrics |>
  group_by(across(all_of(aggregate_groups))) |>
  summarise(
    n_network_replicates = n_distinct(.data$network_replicate),
    recurrent_state_supported_fraction = mean(
      .data$recurrent_state_supported,
      na.rm = TRUE
    ),
    annual_new_infections_mid = safe_median(.data$annual_new_infections),
    annual_new_infections_lo = safe_p10(.data$annual_new_infections),
    annual_new_infections_hi = safe_p90(.data$annual_new_infections),
    annual_new_infections_per_1000_mid = safe_median(
      .data$annual_new_infections_per_1000
    ),
    annual_mean_infectious_hosts_mid = safe_median(
      .data$annual_mean_infectious_hosts
    ),
    annual_mean_infectious_hosts_lo = safe_p10(
      .data$annual_mean_infectious_hosts
    ),
    annual_mean_infectious_hosts_hi = safe_p90(
      .data$annual_mean_infectious_hosts
    ),
    mean_seroprevalence_mid = safe_median(.data$mean_seroprevalence),
    mean_seroprevalence_lo = safe_p10(.data$mean_seroprevalence),
    mean_seroprevalence_hi = safe_p90(.data$mean_seroprevalence),
    mean_annual_foi_probability_mid = safe_median(
      .data$mean_annual_foi_probability
    ),
    annual_new_infections_ratio_vs_baseline_mid = safe_median(
      .data$annual_new_infections_ratio_vs_baseline
    ),
    annual_new_infections_ratio_vs_baseline_lo = safe_p10(
      .data$annual_new_infections_ratio_vs_baseline
    ),
    annual_new_infections_ratio_vs_baseline_hi = safe_p90(
      .data$annual_new_infections_ratio_vs_baseline
    ),
    .groups = "drop"
  )

baseline_anchor <- aggregate_metrics |>
  filter(.data$scenario_type == "current_baseline")
if (nrow(baseline_anchor) != 1L) {
  stop("Expected exactly one aggregated observed-baseline row.")
}

ensemble_metrics <- aggregate_metrics |>
  filter(.data$scenario_type == "future_window") |>
  group_by(
    .data$scenario,
    .data$scenario_label,
    .data$future_year,
    .data$future_window,
    .data$window_start,
    .data$window_end,
    .data$network_mode,
    .data$rainfall_model,
    .data$rainfall_label,
    .data$p_mobile,
    .data$capacity_scale,
    .data$site_profile,
    .data$site_profile_label
  ) |>
  summarise(
    n_climate_models = n_distinct(.data$climate_model),
    n_network_replicates_per_model = min(.data$n_network_replicates),
    annual_new_infections_ensemble_median = median(
      .data$annual_new_infections_mid,
      na.rm = TRUE
    ),
    annual_new_infections_ensemble_lo = min(
      .data$annual_new_infections_mid,
      na.rm = TRUE
    ),
    annual_new_infections_ensemble_hi = max(
      .data$annual_new_infections_mid,
      na.rm = TRUE
    ),
    annual_new_infections_ratio_ensemble_median = median(
      .data$annual_new_infections_ratio_vs_baseline_mid,
      na.rm = TRUE
    ),
    annual_new_infections_ratio_ensemble_lo = min(
      .data$annual_new_infections_ratio_vs_baseline_mid,
      na.rm = TRUE
    ),
    annual_new_infections_ratio_ensemble_hi = max(
      .data$annual_new_infections_ratio_vs_baseline_mid,
      na.rm = TRUE
    ),
    mean_seroprevalence_ensemble_median = median(
      .data$mean_seroprevalence_mid,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

if (any(ensemble_metrics$n_climate_models != length(climate_models))) {
  stop("At least one SSP x window summary is missing a reporting GCM.")
}

write_csv(
  replicate_metrics,
  file.path(future_dir, "future_cold_endemic_window_burden_replicates.csv")
)
write_csv(
  aggregate_metrics,
  file.path(future_dir, "future_cold_endemic_window_burden_aggregate.csv")
)
write_csv(
  ensemble_metrics,
  file.path(future_dir, "future_cold_endemic_window_burden_ensemble.csv")
)

# -------- Publication figures --------------------------------------------
window_levels <- vapply(
  future_windows,
  function(years) sprintf("%d-%d", min(years), max(years)),
  character(1)
)
plot_gcm <- aggregate_metrics |>
  filter(.data$scenario_type == "future_window") |>
  mutate(
    ssp_label = rvf_ssp_label(.data$scenario),
    climate_model = factor(.data$climate_model, levels = climate_models),
    window_label = factor(
      sprintf("%d-%d", .data$window_start, .data$window_end),
      levels = window_levels
    )
  )
plot_ensemble <- ensemble_metrics |>
  mutate(
    ssp_label = rvf_ssp_label(.data$scenario),
    window_label = factor(
      sprintf("%d-%d", .data$window_start, .data$window_end),
      levels = window_levels
    )
  )
baseline_band <- tibble(
  ssp_label = rvf_ssp_label(future_scenarios),
  baseline_mid = baseline_anchor$annual_new_infections_mid[[1]],
  baseline_lo = baseline_anchor$annual_new_infections_lo[[1]],
  baseline_hi = baseline_anchor$annual_new_infections_hi[[1]]
)

dodge_width <- 0.52
p_absolute <- ggplot() +
  geom_rect(
    data = baseline_band,
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = .data$baseline_lo,
      ymax = .data$baseline_hi
    ),
    inherit.aes = FALSE,
    fill = "grey88",
    colour = NA
  ) +
  geom_hline(
    data = baseline_band,
    aes(yintercept = .data$baseline_mid),
    colour = "black",
    linetype = "dashed",
    linewidth = 0.55
  ) +
  geom_linerange(
    data = plot_ensemble,
    aes(
      x = .data$window_label,
      ymin = .data$annual_new_infections_ensemble_lo,
      ymax = .data$annual_new_infections_ensemble_hi
    ),
    inherit.aes = FALSE,
    colour = "grey55",
    linewidth = 1.1,
    alpha = 1
  ) +
  geom_point(
    data = plot_ensemble,
    aes(
      x = .data$window_label,
      y = .data$annual_new_infections_ensemble_median
    ),
    inherit.aes = FALSE,
    colour = "black",
    shape = 18,
    size = 3.2
  ) +
  geom_pointrange(
    data = plot_gcm,
    aes(
      x = .data$window_label,
      y = .data$annual_new_infections_mid,
      ymin = .data$annual_new_infections_lo,
      ymax = .data$annual_new_infections_hi,
      colour = .data$climate_model
    ),
    position = position_dodge(width = dodge_width),
    linewidth = 0.55
  ) +
  facet_ssp(scales = "fixed") +
  scale_colour_gcm(name = "GCM", drop = FALSE) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x = "Future climate window",
    y = "Mean annual new infections",
    colour = "GCM"
  ) +
  theme_rvf(base_size = 13) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(colour = "black", size = 13),
    strip.text = element_text(colour = "black", size = 12),
    legend.text = element_text(colour = "black", size = 10)
  ) +
  guides(colour = guide_legend(nrow = 1))

# Clamp only graphical intervals at a small positive value when a replicate
# interval touches zero. The unchanged numerical values remain in the tables.
positive_ratios <- c(
  plot_gcm$annual_new_infections_ratio_vs_baseline_mid,
  plot_gcm$annual_new_infections_ratio_vs_baseline_lo,
  plot_ensemble$annual_new_infections_ratio_ensemble_median,
  plot_ensemble$annual_new_infections_ratio_ensemble_lo
)
positive_ratios <- positive_ratios[is.finite(positive_ratios) & positive_ratios > 0]
ratio_floor <- if (length(positive_ratios) == 0L) 1 / 32 else min(positive_ratios) / 2
plot_gcm <- plot_gcm |>
  mutate(
    ratio_mid_plot = pmax(.data$annual_new_infections_ratio_vs_baseline_mid, ratio_floor),
    ratio_lo_plot = pmax(.data$annual_new_infections_ratio_vs_baseline_lo, ratio_floor),
    ratio_hi_plot = pmax(.data$annual_new_infections_ratio_vs_baseline_hi, ratio_floor)
  )
plot_ensemble <- plot_ensemble |>
  mutate(
    ratio_mid_plot = pmax(.data$annual_new_infections_ratio_ensemble_median, ratio_floor),
    ratio_lo_plot = pmax(.data$annual_new_infections_ratio_ensemble_lo, ratio_floor),
    ratio_hi_plot = pmax(.data$annual_new_infections_ratio_ensemble_hi, ratio_floor)
  )

p_fold <- ggplot() +
  geom_hline(
    yintercept = 1,
    colour = "black",
    linetype = "dashed",
    linewidth = 0.55
  ) +
  geom_linerange(
    data = plot_ensemble,
    aes(
      x = .data$window_label,
      ymin = .data$ratio_lo_plot,
      ymax = .data$ratio_hi_plot
    ),
    inherit.aes = FALSE,
    colour = "grey55",
    linewidth = 1.1,
    alpha = 1
  ) +
  geom_point(
    data = plot_ensemble,
    aes(x = .data$window_label, y = .data$ratio_mid_plot),
    inherit.aes = FALSE,
    colour = "black",
    shape = 18,
    size = 3.2
  ) +
  geom_pointrange(
    data = plot_gcm,
    aes(
      x = .data$window_label,
      y = .data$ratio_mid_plot,
      ymin = .data$ratio_lo_plot,
      ymax = .data$ratio_hi_plot,
      colour = .data$climate_model
    ),
    position = position_dodge(width = dodge_width),
    linewidth = 0.55
  ) +
  facet_ssp(scales = "fixed") +
  scale_colour_gcm(name = "GCM", drop = FALSE) +
  scale_y_continuous(
    trans = "log2",
    breaks = c(1 / 8, 1 / 4, 1 / 2, 1, 2, 4, 8, 16, 32),
    labels = function(x) paste0(formatC(x, format = "g"), "x")
  ) +
  labs(
    x = "Future climate window",
    y = "Annual infections relative to observed-climate baseline",
    colour = "GCM"
  ) +
  theme_rvf(base_size = 13) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(colour = "black", size = 13),
    strip.text = element_text(colour = "black", size = 12),
    legend.text = element_text(colour = "black", size = 10)
  ) +
  guides(colour = guide_legend(nrow = 1))

# -------- Two-profile publication comparison ----------------------------
# Build this figure only after both calibrated profile analyses exist. The
# profile facets share one y-scale so absolute burdens can be compared without
# conflating profile differences with panel scaling.
primary_future_dir <- file.path(
  cfg$results_dir,
  rvf_rainfall_output_subdir(
    "future_cold_endemic",
    sim$rainfall_model,
    cfg
  )
)
contrast_future_dir <- file.path(
  cfg$results_dir,
  rvf_rainfall_output_subdir(
    "future_cold_endemic_profile_peridomestic_dominant",
    sim$rainfall_model,
    cfg
  )
)
profile_aggregate_paths <- file.path(
  c(primary_future_dir, contrast_future_dir),
  "future_cold_endemic_window_burden_aggregate.csv"
)
profile_ensemble_paths <- file.path(
  c(primary_future_dir, contrast_future_dir),
  "future_cold_endemic_window_burden_ensemble.csv"
)

if (all(file.exists(c(profile_aggregate_paths, profile_ensemble_paths)))) {
  profile_levels <- c(
    "grazing_floodplain_dominant",
    "peridomestic_dominant"
  )
  profile_labels <- c(
    grazing_floodplain_dominant = "Shared grazing/floodplain",
    peridomestic_dominant = "Village-associated"
  )

  combined_aggregate <- bind_rows(lapply(
    profile_aggregate_paths,
    read_csv,
    show_col_types = FALSE
  ))
  combined_ensemble <- bind_rows(lapply(
    profile_ensemble_paths,
    read_csv,
    show_col_types = FALSE
  ))

  combined_gcm <- combined_aggregate |>
    filter(.data$scenario_type == "future_window") |>
    mutate(
      ssp_label = rvf_ssp_label(.data$scenario),
      profile_label = factor(
        unname(profile_labels[.data$site_profile]),
        levels = unname(profile_labels[profile_levels])
      ),
      climate_model = factor(.data$climate_model, levels = climate_models),
      window_label = factor(
        sprintf("%d-%d", .data$window_start, .data$window_end),
        levels = window_levels
      )
    )
  combined_ensemble_plot <- combined_ensemble |>
    mutate(
      ssp_label = rvf_ssp_label(.data$scenario),
      profile_label = factor(
        unname(profile_labels[.data$site_profile]),
        levels = unname(profile_labels[profile_levels])
      ),
      window_label = factor(
        sprintf("%d-%d", .data$window_start, .data$window_end),
        levels = window_levels
      )
    )
  combined_baseline <- combined_aggregate |>
    filter(.data$scenario_type == "current_baseline") |>
    transmute(
      profile_label = factor(
        unname(profile_labels[.data$site_profile]),
        levels = unname(profile_labels[profile_levels])
      ),
      baseline_mid = .data$annual_new_infections_mid,
      baseline_lo = .data$annual_new_infections_lo,
      baseline_hi = .data$annual_new_infections_hi
    ) |>
    tidyr::crossing(
      ssp_label = factor(rvf_ssp_levels, levels = rvf_ssp_levels)
    )

  p_profiles_absolute <- ggplot() +
    geom_rect(
      data = combined_baseline,
      aes(
        xmin = -Inf,
        xmax = Inf,
        ymin = .data$baseline_lo,
        ymax = .data$baseline_hi
      ),
      inherit.aes = FALSE,
      fill = "grey88",
      colour = NA
    ) +
    geom_hline(
      data = combined_baseline,
      aes(yintercept = .data$baseline_mid),
      colour = "black",
      linetype = "dashed",
      linewidth = 0.55
    ) +
    geom_linerange(
      data = combined_ensemble_plot,
      aes(
        x = .data$window_label,
        ymin = .data$annual_new_infections_ensemble_lo,
        ymax = .data$annual_new_infections_ensemble_hi
      ),
      inherit.aes = FALSE,
      colour = "grey55",
      linewidth = 1.0
    ) +
    geom_point(
      data = combined_ensemble_plot,
      aes(
        x = .data$window_label,
        y = .data$annual_new_infections_ensemble_median
      ),
      inherit.aes = FALSE,
      colour = "black",
      shape = 18,
      size = 3.0
    ) +
    geom_pointrange(
      data = combined_gcm,
      aes(
        x = .data$window_label,
        y = .data$annual_new_infections_mid,
        ymin = .data$annual_new_infections_lo,
        ymax = .data$annual_new_infections_hi,
        colour = .data$climate_model
      ),
      position = position_dodge(width = dodge_width),
      linewidth = 0.50
    ) +
    facet_grid(
      rows = vars(.data$profile_label),
      cols = vars(.data$ssp_label),
      scales = "fixed"
    ) +
    scale_colour_gcm(name = "GCM", drop = FALSE) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(
      x = "Future climate window",
      y = "Mean annual new infections",
      colour = "GCM"
    ) +
    theme_rvf(base_size = 13) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = 10.5),
      axis.title = element_text(colour = "black", size = 13),
      strip.text = element_text(colour = "black", size = 12),
      legend.text = element_text(colour = "black", size = 10.5),
      # Separate profile rows so adjacent tick labels cannot appear merged.
      panel.spacing.y = grid::unit(0.65, "cm")
    ) +
    guides(colour = guide_legend(nrow = 1))

  combined_gcm <- combined_gcm |>
    mutate(
      ratio_mid_plot = pmax(
        .data$annual_new_infections_ratio_vs_baseline_mid,
        ratio_floor
      ),
      ratio_lo_plot = pmax(
        .data$annual_new_infections_ratio_vs_baseline_lo,
        ratio_floor
      ),
      ratio_hi_plot = pmax(
        .data$annual_new_infections_ratio_vs_baseline_hi,
        ratio_floor
      )
    )
  combined_ensemble_plot <- combined_ensemble_plot |>
    mutate(
      ratio_mid_plot = pmax(
        .data$annual_new_infections_ratio_ensemble_median,
        ratio_floor
      ),
      ratio_lo_plot = pmax(
        .data$annual_new_infections_ratio_ensemble_lo,
        ratio_floor
      ),
      ratio_hi_plot = pmax(
        .data$annual_new_infections_ratio_ensemble_hi,
        ratio_floor
      )
    )

  p_profiles_fold <- ggplot() +
    geom_hline(
      yintercept = 1,
      colour = "black",
      linetype = "dashed",
      linewidth = 0.55
    ) +
    geom_linerange(
      data = combined_ensemble_plot,
      aes(
        x = .data$window_label,
        ymin = .data$ratio_lo_plot,
        ymax = .data$ratio_hi_plot
      ),
      inherit.aes = FALSE,
      colour = "grey55",
      linewidth = 1.0
    ) +
    geom_point(
      data = combined_ensemble_plot,
      aes(
        x = .data$window_label,
        y = .data$ratio_mid_plot,
        shape = "Four-GCM median"
      ),
      inherit.aes = FALSE,
      colour = "black",
      size = 3.7
    ) +
    geom_pointrange(
      data = combined_gcm,
      aes(
        x = .data$window_label,
        y = .data$ratio_mid_plot,
        ymin = .data$ratio_lo_plot,
        ymax = .data$ratio_hi_plot,
        colour = .data$climate_model
      ),
      position = position_dodge(width = dodge_width),
      linewidth = 0.50
    ) +
    facet_grid(
      rows = vars(.data$profile_label),
      cols = vars(.data$ssp_label),
      scales = "fixed"
    ) +
    scale_colour_gcm(name = "GCM", drop = FALSE) +
    scale_shape_manual(
      values = c("Four-GCM median" = 18),
      name = NULL
    ) +
    scale_y_continuous(
      limits = c(0.75, NA),
      breaks = 1:8,
      labels = function(x) paste0(x, "x"),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x = "Future climate window",
      y = "Annual infections relative to observed-climate baseline",
      colour = "GCM"
    ) +
    theme_rvf(base_size = 13) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = 10.5),
      axis.title = element_text(colour = "black", size = 13),
      strip.text = element_text(colour = "black", size = 12),
      legend.text = element_text(colour = "black", size = 10.5)
    ) +
    guides(
      colour = guide_legend(nrow = 1, order = 1),
      shape = guide_legend(order = 2)
    )

  comparison_dir <- file.path(cfg$results_dir, "future_profile_comparison")
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(
    combined_ensemble,
    file.path(comparison_dir, "future_profile_comparison_ensemble.csv")
  )
  rvf_ggsave(
    file.path(
      cfg$figures_dir,
      "04_future_window_annual_infections_profiles.png"
    ),
    p_profiles_absolute,
    width = 10.5,
    height = 8.2
  )
  rvf_ggsave(
    file.path(
      cfg$figures_dir,
      "04_future_window_annual_infections_profiles_fold.png"
    ),
    p_profiles_fold,
    width = 10.5,
    height = 8.7
  )
}

# -------- Rainfall-response sensitivity ---------------------------------
# When this script summarises a non-primary rainfall response, compare it with
# the primary Hill run using only network replicates available under both
# formulations. This keeps the robustness comparison paired and prevents a
# replicate-set difference from being mistaken for a rainfall-curve effect.
primary_replicate_path <- file.path(
  cfg$results_dir,
  rvf_rainfall_output_subdir(
    "future_cold_endemic",
    cfg$primary_rainfall_model,
    cfg
  ),
  "future_cold_endemic_window_burden_replicates.csv"
)

if (!is_primary_rainfall && is_primary_profile &&
    file.exists(primary_replicate_path)) {
  primary_replicates <- read_csv(
    primary_replicate_path,
    show_col_types = FALSE
  )
  common_replicates <- intersect(
    unique(primary_replicates$network_replicate),
    unique(replicate_metrics$network_replicate)
  )
  if (length(common_replicates) == 0L) {
    stop("No common network replicates for the rainfall-response comparison.")
  }

  rainfall_gcm_summary <- bind_rows(
    primary_replicates |>
      select(
        scenario_type,
        scenario,
        scenario_label,
        climate_model,
        future_window,
        window_start,
        window_end,
        network_replicate,
        annual_new_infections,
        annual_new_infections_ratio_vs_baseline
      ) |>
      mutate(
        rainfall_short = "hill",
        rainfall_comparison_label = "Saturating Hill"
      ),
    replicate_metrics |>
      select(
        scenario_type,
        scenario,
        scenario_label,
        climate_model,
        future_window,
        window_start,
        window_end,
        network_replicate,
        annual_new_infections,
        annual_new_infections_ratio_vs_baseline
      ) |>
      mutate(
        rainfall_short = "quadratic",
        rainfall_comparison_label = "Quadratic"
      )
  ) |>
    filter(
      .data$scenario_type == "future_window",
      .data$network_replicate %in% common_replicates
    ) |>
    group_by(
      .data$rainfall_short,
      .data$rainfall_comparison_label,
      .data$scenario,
      .data$scenario_label,
      .data$future_window,
      .data$window_start,
      .data$window_end,
      .data$climate_model
    ) |>
    summarise(
      annual_new_infections = median(
        .data$annual_new_infections,
        na.rm = TRUE
      ),
      annual_infections_relative_to_baseline = median(
        .data$annual_new_infections_ratio_vs_baseline,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  rainfall_ensemble_summary <- rainfall_gcm_summary |>
    group_by(
      .data$rainfall_short,
      .data$rainfall_comparison_label,
      .data$scenario,
      .data$scenario_label,
      .data$future_window,
      .data$window_start,
      .data$window_end
    ) |>
    summarise(
      annual_new_infections_min = min(.data$annual_new_infections),
      annual_new_infections_max = max(.data$annual_new_infections),
      annual_new_infections = median(.data$annual_new_infections),
      relative_to_baseline_min = min(
        .data$annual_infections_relative_to_baseline
      ),
      relative_to_baseline_max = max(
        .data$annual_infections_relative_to_baseline
      ),
      annual_infections_relative_to_baseline = median(
        .data$annual_infections_relative_to_baseline
      ),
      .groups = "drop"
    )

  rainfall_comparison_table <- rainfall_ensemble_summary |>
    select(-rainfall_comparison_label) |>
    tidyr::pivot_wider(
      names_from = rainfall_short,
      values_from = c(
        annual_new_infections,
        annual_new_infections_min,
        annual_new_infections_max,
        annual_infections_relative_to_baseline,
        relative_to_baseline_min,
        relative_to_baseline_max
      ),
      names_glue = "{.value}_{rainfall_short}"
    ) |>
    mutate(
      absolute_burden_difference_percent = 100 * (
        .data$annual_new_infections_quadratic /
          .data$annual_new_infections_hill - 1
      ),
      fold_change_difference =
        .data$annual_infections_relative_to_baseline_quadratic -
        .data$annual_infections_relative_to_baseline_hill
    )

  write_csv(
    rainfall_comparison_table,
    file.path(future_dir, "rainfall_response_sensitivity.csv")
  )

  rainfall_plot_data <- rainfall_ensemble_summary |>
    mutate(
      ssp_label = factor(
        rvf_ssp_label(.data$scenario),
        levels = rvf_ssp_levels
      ),
      window_label = factor(
        sprintf("%d-%d", .data$window_start, .data$window_end),
        levels = window_levels
      )
    )
  p_rainfall_sensitivity <- ggplot(
    rainfall_plot_data,
    aes(
      x = .data$window_label,
      y = .data$annual_infections_relative_to_baseline,
      ymin = .data$relative_to_baseline_min,
      ymax = .data$relative_to_baseline_max,
      colour = .data$rainfall_comparison_label,
      shape = .data$rainfall_comparison_label
    )
  ) +
    geom_hline(yintercept = 1, colour = "black", linewidth = 0.5) +
    geom_pointrange(
      position = position_dodge(width = 0.35),
      linewidth = 0.65
    ) +
    facet_wrap(vars(.data$ssp_label), nrow = 1) +
    scale_colour_manual(
      values = c("Saturating Hill" = "#0072B2", "Quadratic" = "#D55E00"),
      breaks = c("Saturating Hill", "Quadratic"),
      name = "Rainfall response"
    ) +
    scale_shape_manual(
      values = c("Saturating Hill" = 16, "Quadratic" = 17),
      breaks = c("Saturating Hill", "Quadratic"),
      name = "Rainfall response"
    ) +
    scale_y_continuous(breaks = 1:8) +
    labs(
      x = "Future climate window",
      y = "Annual infections relative to\nobserved-climate baseline"
    ) +
    theme_rvf(base_size = 13) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black"),
      strip.text = element_text(colour = "black")
    ) +
    guides(
      colour = guide_legend(nrow = 1),
      shape = guide_legend(nrow = 1)
    )

  rvf_ggsave(
    file.path(figures_dir, "04_rainfall_response_sensitivity.png"),
    p_rainfall_sensitivity,
    width = 9.2,
    height = 4.8
  )
}

cat("\nObserved annual burden over the common 40-year reporting period:\n")
baseline_anchor |>
  select(
    n_network_replicates,
    annual_new_infections_mid,
    annual_new_infections_lo,
    annual_new_infections_hi,
    mean_seroprevalence_mid
  ) |>
  print(width = Inf)

cat("\nFuture ensemble summaries:\n")
ensemble_metrics |>
  select(
    scenario_label,
    future_window,
    n_climate_models,
    annual_new_infections_ensemble_median,
    annual_new_infections_ensemble_lo,
    annual_new_infections_ensemble_hi,
    annual_new_infections_ratio_ensemble_median
  ) |>
  arrange(.data$scenario_label, .data$future_window) |>
  print(n = Inf, width = Inf)

message("Future burden summaries and figures written to: ", future_dir)
