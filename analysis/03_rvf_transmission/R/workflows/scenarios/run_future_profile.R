#!/usr/bin/env Rscript

#' 06_run_future_scenarios.R
#' -------------------------------------------------------------------------
#' Compare RVFV endemic potential under repeated future climate and movement
#' windows. Each SSP x GCM x 20-year window x network-replicate run starts
#' independently and receives one infectious cow in every village. The
#' publication run uses the profile and rainfall response declared in
#' R/config.R. Contrasting profile and rainfall-sensitivity runs write to
#' separate output directories.
#'
#' Each 20-year window starts from the same initial state and is simulated
#' separately, rather than as one continuous run from the present to 2100.
#' Script 06 writes the simulation outputs and recurrent-cycle diagnostics,
#' and script 07 performs the burden aggregation and plotting.
#'
#' A four-cycle workflow test can be enabled under `execution$smoke` in
#' R/config.R. Smoke outputs are isolated and are not suitable for inference.
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
sim <- rvf_apply_profile_calibration(
  sim,
  cfg = cfg,
  require_target_band = isTRUE(cfg$execution$require_calibration_target)
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

cycle_lengths <- vapply(future_windows, length, integer(1))
if (length(unique(cycle_lengths)) != 1L) {
  stop("All future forcing windows must contain the same number of years.")
}
forcing_cycle_years <- unique(cycle_lengths)
reporting_years <- as.integer(cfg$execution$reporting_years)
if (reporting_years %% forcing_cycle_years != 0L) {
  stop("The reporting duration must contain complete future forcing cycles.")
}
terminal_cycles <- if (smoke_mode) {
  2L
} else {
  reporting_years %/% forcing_cycle_years
}
window_cycles <- c(2L, 4L, 6L)
required_replicate_fraction <- 0.80
relative_tolerance <- 0.10
seroprevalence_tolerance_pp <- 0.50
seed_per_village <- 1L
seed_village <- "all"
years_burn <- sim$years_burn
years_endemic <- if (smoke_mode) {
  4L * forcing_cycle_years
} else {
  ceiling(sim$years_endemic / forcing_cycle_years) * forcing_cycle_years
}

if (years_endemic < 2L * terminal_cycles * forcing_cycle_years) {
  stop("Future runs require at least four complete forcing cycles.")
}
if (!isTRUE(all.equal(as.numeric(sim$p_mobile), 1))) {
  stop("The canonical future comparison requires p_mobile = 1.")
}
if (!isTRUE(all.equal(as.numeric(sim$external_introduction_per_year), 0))) {
  stop("The canonical future comparison requires zero external introductions.")
}
if (length(network_replicates) == 0L) {
  stop("No analysis network replicates are configured.")
}

calibrated <- rvf_adjust_calibrated_params(
  transmission_scale = sim$transmission_scale,
  capacity_scale = sim$capacity_scale,
  site_profile = sim$site_profile,
  rainfall_model = sim$rainfall_model
)

out_name <- if (smoke_mode) {
  "future_endemic_smoke"
} else if (identical(sim$network_mode, "monthly")) {
  "future_cold_endemic"
} else {
  "future_cold_endemic_annualized"
}
# Keep the canonical profile at the established path used by the manuscript.
# Contrasting profiles receive isolated folders and cannot overwrite it.
if (!is_primary_profile) {
  profile_tag <- gsub("[^A-Za-z0-9]+", "_", sim$site_profile)
  out_name <- paste0(out_name, "_profile_", profile_tag)
}
out_name <- rvf_rainfall_output_subdir(out_name, sim$rainfall_model, cfg)
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

output_paths <- c(
  status = file.path(out_dir, "future_cold_endemic_run_status.csv"),
  summary = file.path(out_dir, "future_cold_endemic_summary.csv"),
  convergence = file.path(out_dir, "future_cold_endemic_convergence.csv"),
  annual = file.path(intermediate_dir, "future_cold_endemic_annual_system.rds"),
  village = file.path(intermediate_dir, "future_cold_endemic_village_terminal.rds"),
  cycles = file.path(intermediate_dir, "future_cold_endemic_cycle_diagnostics.rds")
)

# Remove the previous output contract before a new run. If the simulation is
# interrupted, script 07 will then fail clearly instead of reading stale files.
unlink(output_paths[file.exists(output_paths)])

# -------- Complete SSP x GCM x window x replicate grid -------------------
scenario_index <- list_rvf_scenarios(cfg)
future_cells <- tidyr::crossing(
  scenario = future_scenarios,
  climate_model = climate_models,
  future_window = names(future_windows)
)

future_index <- bind_rows(lapply(seq_len(nrow(future_cells)), function(i) {
  cell <- future_cells[i, ]
  window_years <- future_windows[[cell$future_window]]
  complete_replicates <- rvf_complete_replicates_for_years(
    index = scenario_index,
    dataset = "future",
    years = window_years,
    climate_model = cell$climate_model,
    network_replicate = network_replicates,
    scenario = cell$scenario
  )

  tibble(
    scenario = cell$scenario,
    scenario_label = unname(cfg$scenario_labels[[cell$scenario]]),
    climate_model = cell$climate_model,
    future_year = as.integer(names(future_windows)[
      match(cell$future_window, names(future_windows))
    ]),
    future_window = cell$future_window,
    window_start = min(window_years),
    window_end = max(window_years),
    window_years = paste(window_years, collapse = ","),
    network_replicate = network_replicates,
    forcing_complete = network_replicates %in% complete_replicates
  )
})) |>
  arrange(.data$scenario, .data$climate_model, .data$future_year,
          .data$network_replicate)

missing_forcing <- future_index |>
  filter(!.data$forcing_complete)
if (nrow(missing_forcing) > 0L) {
  missing_preview <- missing_forcing |>
    transmute(
      cell = sprintf(
        "%s / %s / %s / replicate %03d",
        .data$scenario,
        .data$climate_model,
        .data$future_window,
        .data$network_replicate
      )
    ) |>
    pull("cell")
  stop(
    "Incomplete future forcing grid. Missing:\n  ",
    paste(head(missing_preview, 20L), collapse = "\n  "),
    if (length(missing_preview) > 20L) "\n  ..." else "",
    "\nRerun Module 02 projections and Module 03 script 01."
  )
}
future_index <- future_index |>
  select(-forcing_complete)

# Function: add one run's metadata without repeating large mutate blocks.
append_run_metadata <- function(data, metadata) {
  if (nrow(data) == 0L) return(data)
  metadata_columns <- setdiff(names(metadata), names(data))
  bind_cols(
    data,
    metadata[rep(1L, nrow(data)), metadata_columns, drop = FALSE]
  )
}

# Function: record the assumptions and identifiers for one future run.
run_metadata <- function(job) {
  tibble(
    source = "future_repeated_window",
    scenario_design = "independent_endemic_potential_scenario",
    scenario = job$scenario[[1]],
    scenario_label = job$scenario_label[[1]],
    climate_model = job$climate_model[[1]],
    future_year = job$future_year[[1]],
    future_window = job$future_window[[1]],
    window_start = job$window_start[[1]],
    window_end = job$window_end[[1]],
    window_years = job$window_years[[1]],
    network_replicate = job$network_replicate[[1]],
    network_mode = sim$network_mode,
    p_mobile = sim$p_mobile,
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
    m_max_w = calibrated$m_max_w
  )
}

# Function: run one SSP x GCM x window x network-replicate scenario.
run_one <- function(job) {
  metadata <- run_metadata(job)
  window_years <- as.integer(strsplit(
    job$window_years[[1]],
    ",",
    fixed = TRUE
  )[[1]])
  scenario_base_ids <- sprintf(
    "%s_%s_%d",
    toupper(job$scenario[[1]]),
    job$climate_model[[1]],
    window_years
  )
  scenario_ids <- rvf_with_replicate(
    scenario_base_ids,
    job$network_replicate[[1]]
  )
  scenario <- build_rvf_sequence_scenario(
    scenario_ids,
    cfg = cfg,
    validate = TRUE,
    scenario_type = "future_window",
    sequence_label = sprintf(
      "%s_%s_%s_WINDOW",
      toupper(job$scenario[[1]]),
      job$climate_model[[1]],
      job$future_window[[1]]
    )
  )
  scenario <- rvf_apply_network_mode(
    scenario,
    network_mode = sim$network_mode
  )

  params <- rvf_assemble_params(
    scenario = scenario,
    calibrated = calibrated,
    use_temp_dependence = TRUE,
    p_mobile = sim$p_mobile,
    # The canonical future design uses full mobility, so the host denominator
    # defining site capacity is fixed explicitly at that reference.
    capacity_reference_p_mobile = 1,
    waning_years = sim$waning_years,
    external_introduction_per_year = sim$external_introduction_per_year,
    verbose = FALSE
  )
  result <- rvf_run_simulation(
    params = params,
    years_burn = years_burn,
    years_endemic = years_endemic,
    output_step = sim$output_step,
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

  # Spatial summaries use the same final reporting block as the system
  # comparison.
  terminal_start_year <-
    years_endemic - terminal_cycles * forcing_cycle_years + 1L
  village_terminal <- annual$village_year |>
    filter(.data$analysis_year >= terminal_start_year) |>
    group_by(.data$village_id) |>
    summarise(
      terminal_window_years = n_distinct(.data$analysis_year),
      annual_new_infections = mean(.data$annual_new_infections, na.rm = TRUE),
      annual_mean_infectious_hosts = mean(
        .data$annual_mean_infectious_hosts,
        na.rm = TRUE
      ),
      annual_mean_prevalence = mean(
        .data$annual_mean_prevalence,
        na.rm = TRUE
      ),
      annual_mean_daily_foi = mean(
        .data$annual_mean_daily_foi,
        na.rm = TRUE
      ),
      annual_cum_hazard = mean(.data$annual_cum_hazard, na.rm = TRUE),
      annual_foi_probability = mean(
        .data$annual_foi_probability,
        na.rm = TRUE
      ),
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

# -------- Run and validate the complete future grid ----------------------
n_cores <- rvf_max_cores(
  default = if (smoke_mode) 2L else cfg$execution$max_cores,
  n_tasks = nrow(future_index)
)
message(sprintf(
  "Running %d future endemic-potential scenario(s) on %d core(s).",
  nrow(future_index),
  n_cores
))
message(sprintf(
  "Protocol: %d years (%d x %d-year cycles), profile=%s, rainfall=%s.",
  years_endemic,
  years_endemic / forcing_cycle_years,
  forcing_cycle_years,
  sim$site_profile,
  sim$rainfall_model
))

results <- rvf_parallel_lapply(seq_len(nrow(future_index)), function(i) {
  job <- future_index[i, ]
  message(sprintf(
    "Future scenario %d/%d: %s %s %s replicate %03d",
    i,
    nrow(future_index),
    job$scenario,
    job$climate_model,
    job$future_window,
    job$network_replicate
  ))

  tryCatch(
    run_one(job),
    error = function(error) {
      list(
        status = run_metadata(job) |>
          mutate(
            run_status = "failed",
            error_message = conditionMessage(error)
          ),
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
write_csv(run_status, output_paths[["status"]])

failed_runs <- run_status |>
  filter(.data$run_status != "ok")
if (nrow(failed_runs) > 0L) {
  stop(
    nrow(failed_runs),
    " future run(s) failed. No scientific outputs were written; inspect ",
    output_paths[["status"]]
  )
}

annual_system <- bind_result("annual_system")
time_series <- bind_result("time_series")
village_terminal <- bind_result("village_terminal")
if (nrow(annual_system) == 0L || nrow(time_series) == 0L) {
  stop("Future runs returned no annual or within-year output.")
}

# -------- Phase-matched recurrent-state diagnostics ---------------------
diagnostic_groups <- c(
  "scenario",
  "climate_model",
  "future_window",
  "network_replicate",
  "network_mode",
  "external_introduction_per_year",
  "rainfall_model",
  "site_profile"
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

terminal_reference <- cycle_diagnostics$window_sensitivity |>
  filter(.data$summary_cycles == terminal_cycles) |>
  transmute(
    across(all_of(join_keys)),
    reference_window_infections = .data$annual_new_infections,
    reference_window_seroprevalence = .data$mean_seroprevalence_end
  )
window_sensitivity_by_run <- cycle_diagnostics$window_sensitivity |>
  filter(.data$summary_cycles > terminal_cycles) |>
  left_join(terminal_reference, by = join_keys) |>
  mutate(
    incidence_relative_difference = rvf_symmetric_relative_change(
      .data$annual_new_infections,
      .data$reference_window_infections
    ),
    seroprevalence_difference_pp = abs(
      .data$mean_seroprevalence_end -
        .data$reference_window_seroprevalence
    )
  ) |>
  group_by(across(all_of(join_keys))) |>
  summarise(
    window_sensitivity_incidence_max_relative_difference = max(
      .data$incidence_relative_difference,
      na.rm = TRUE
    ),
    window_sensitivity_seroprevalence_max_difference_pp = max(
      .data$seroprevalence_difference_pp,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

summary_by_run <- run_status |>
  left_join(terminal_by_run, by = join_keys) |>
  left_join(final_state_by_run, by = join_keys) |>
  left_join(cycle_diagnostics$convergence_by_run, by = join_keys) |>
  left_join(cycle_diagnostics$periodicity_by_run, by = join_keys) |>
  left_join(cycle_diagnostics$persistence_by_run, by = join_keys) |>
  left_join(window_sensitivity_by_run, by = join_keys) |>
  mutate(
    persistence_class = case_when(
      .data$fraction_time_below_1e_6_infectious_hosts > 0 ~
        "near_deterministic_extinction",
      .data$fraction_time_below_one_infectious_host > 0 ~
        "fractional_deterministic_persistence",
      TRUE ~ "persistent_above_one_host"
    )
  )

status_by_cell <- run_status |>
  group_by(
    .data$scenario,
    .data$climate_model,
    .data$future_window,
    .data$rainfall_model,
    .data$site_profile
  ) |>
  summarise(
    n_planned = n(),
    n_successful = sum(.data$run_status == "ok"),
    success_fraction = .data$n_successful / .data$n_planned,
    .groups = "drop"
  )

convergence <- cycle_diagnostics$convergence_summary |>
  left_join(
    status_by_cell,
    by = c(
      "scenario",
      "climate_model",
      "future_window",
      "rainfall_model",
      "site_profile"
    )
  ) |>
  mutate(
    run_complete = .data$n_successful == .data$n_planned,
    quasi_equilibrium_supported = .data$recurrent_state_supported,
    headline_supported =
      .data$run_complete & .data$quasi_equilibrium_supported,
    outcome_class = case_when(
      !.data$run_complete ~ "incomplete_run_cell",
      !.data$quasi_equilibrium_supported ~ "finite_horizon_not_settled",
      TRUE ~ "recurrent_state_supported"
    )
  ) |>
  arrange(.data$scenario, .data$climate_model, .data$future_window)

cycle_diagnostics$analysis_settings <- list(
  smoke_mode = smoke_mode,
  future_scenarios = future_scenarios,
  climate_models = climate_models,
  network_replicates = network_replicates,
  forcing_cycle_years = forcing_cycle_years,
  terminal_cycles = terminal_cycles,
  years_burn = years_burn,
  years_endemic = years_endemic,
  site_profile = sim$site_profile,
  rainfall_model = sim$rainfall_model,
  seeding = list(mode = "all_villages", infectious_cattle_per_village = 1L)
)

# Weekly trajectories have served their diagnostic purpose and are not saved.
write_csv(summary_by_run, output_paths[["summary"]])
write_csv(convergence, output_paths[["convergence"]])
saveRDS(annual_system, output_paths[["annual"]])
saveRDS(village_terminal, output_paths[["village"]])
saveRDS(cycle_diagnostics, output_paths[["cycles"]])

if (!smoke_mode && any(!convergence$headline_supported)) {
  warning(
    "At least one future cell lacks recurrent-state support. Inspect ",
    output_paths[["convergence"]],
    " before running script 07.",
    call. = FALSE
  )
}

message("Future endemic-potential outputs written to: ", out_dir)
if (smoke_mode) {
  message("Smoke test complete. Do not interpret its finite-horizon estimates.")
}
