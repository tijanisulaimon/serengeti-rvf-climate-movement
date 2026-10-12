#!/usr/bin/env Rscript

#' 08_run_climate_movement_attribution.R
#' -------------------------------------------------------------------------
#' Attribute future endemic-potential burden to climate, movement-network
#' change, and their interaction using a matched 2 x 2 factorial:
#'
#'   A = current climate + current movement
#'   B = future climate  + current movement
#'   C = current climate + future movement
#'   D = future climate  + future movement
#'
#' Every cell is initialised independently, seeded with one infectious cow in
#' every village, and run under repeated forcing. The publication analysis uses
#' both SSPs, four GCMs, three future windows, and the paired network
#' replicates configured in `R/config.R`.
#'
#' A two-replicate, four-cycle workflow test can be enabled under
#' `execution$smoke` in R/config.R. Smoke outputs are isolated and must not be
#' interpreted scientifically.
#'
#' The publication analysis is run separately for the primary
#' grazing/floodplain profile and the contrasting peridomestic profile. The
#' configured profile determines which calibrated run is produced; separate
#' output directories keep their expensive caches and summaries independent.
#' -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

source(here::here("analysis", "03_rvf_transmission", "R", "transmission_model.R"))
source(here::here("analysis", "03_rvf_transmission", "R", "cycle_diagnostics.R"))
source(here::here("analysis", "03_rvf_transmission", "R", "plot_helpers.R"))

cfg <- load_rvf_config()

# -------- Analysis specification -----------------------------------------
sim <- cfg$sim
if (exists("analysis_profile", inherits = FALSE)) {
  sim$site_profile <- analysis_profile
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

if (!identical(sim$network_mode, "monthly")) {
  stop("The publication attribution requires monthly movement networks.")
}
if (!isTRUE(all.equal(as.numeric(sim$p_mobile), 1))) {
  stop("The publication attribution requires p_mobile = 1.")
}
if (!isTRUE(all.equal(as.numeric(sim$external_introduction_per_year), 0))) {
  stop("The publication attribution requires zero external introductions.")
}

smoke_mode <- isTRUE(cfg$execution$smoke$attribution)
reuse_cache <- isTRUE(cfg$execution$reuse$attribution)
future_scenarios <- cfg$future_scenarios
climate_models <- cfg$climate_models
future_windows <- cfg$future_windows

# The factorial is substantially larger than the direct future comparison.
# Keep the paired replicate set explicit because each A/B/C/D contrast must use
# the same network realisation.
n_attribution_replicates <- length(cfg$execution$attribution_replicates)
network_replicates <- head(
  intersect(
    sort(unique(cfg$analysis_network_replicates)),
    sort(unique(cfg$execution$attribution_replicates))
  ),
  n_attribution_replicates
)

if (smoke_mode) {
  future_scenarios <- future_scenarios[1]
  climate_models <- climate_models[1]
  future_windows <- future_windows[1]
  network_replicates <- head(network_replicates, 2L)
}
if (length(network_replicates) < if (smoke_mode) 2L else n_attribution_replicates) {
  stop("Insufficient analysis network replicates for the attribution design.")
}

cycle_lengths <- vapply(future_windows, length, integer(1))
if (length(unique(cycle_lengths)) != 1L) {
  stop("All attribution forcing windows must contain the same number of years.")
}
forcing_cycle_years <- unique(cycle_lengths) # 20 years in the current design.
reporting_years <- as.integer(cfg$execution$reporting_years)
if (reporting_years %% forcing_cycle_years != 0L) {
  stop("The reporting duration must contain complete attribution cycles.")
}
terminal_cycles <- if (smoke_mode) {
  2L
} else {
  reporting_years %/% forcing_cycle_years
}
window_cycles <- if (smoke_mode) c(2L, 4L) else c(2L, 4L, 6L)
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
  stop("Attribution runs require at least four complete forcing cycles.")
}

calibrated <- rvf_adjust_calibrated_params(
  transmission_scale = sim$transmission_scale,
  capacity_scale = sim$capacity_scale,
  site_profile = sim$site_profile,
  rainfall_model = sim$rainfall_model
)
total_hosts <- rvf_baseline_total_hosts(sim, cfg)

out_name <- if (smoke_mode) {
  "attribution_smoke"
} else {
  "climate_movement_attribution_cold_start"
}
# Keep the primary profile at the paths used by the manuscript. Contrasting
# profiles receive isolated folders so they cannot overwrite primary results.
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
  status = file.path(out_dir, "climate_movement_attribution_run_status.csv"),
  cells = file.path(out_dir, "climate_movement_attribution_aggregate.csv"),
  convergence = file.path(out_dir, "climate_movement_attribution_convergence.csv"),
  decomposition = file.path(
    out_dir,
    "climate_movement_attribution_decomposition_aggregate.csv"
  ),
  decomposition_replicates = file.path(
    out_dir,
    "climate_movement_attribution_decomposition_replicates.csv"
  ),
  window_summary = file.path(
    out_dir,
    "climate_movement_attribution_decomposition_window_summary.csv"
  ),
  cache = file.path(intermediate_dir, "climate_movement_attribution_runs.rds")
)

if (!reuse_cache) {
  # Remove only files in the active output contract. An interrupted run cannot
  # then leave apparently complete summaries from an older protocol.
  unlink(output_paths[file.exists(output_paths)])
}

safe_median <- function(x) {
  if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}
safe_p10 <- function(x) {
  if (all(is.na(x))) NA_real_ else as.numeric(quantile(x, 0.10, na.rm = TRUE))
}
safe_p90 <- function(x) {
  if (all(is.na(x))) NA_real_ else as.numeric(quantile(x, 0.90, na.rm = TRUE))
}
safe_min <- function(x) {
  if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}
safe_max <- function(x) {
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}
safe_ratio <- function(x, denominator) {
  ifelse(
    is.finite(denominator) & abs(denominator) > .Machine$double.eps,
    x / denominator,
    NA_real_
  )
}

# Function: construct scenario identifiers for one climate or movement source.
# period: "current" uses the observed cycle; "future" uses one SSP/GCM window.
# years: ordered years forming the 20-year forcing sequence.
base_ids_for_period <- function(period, years, climate_model, scenario) {
  if (identical(period, "current")) {
    sprintf("OBSERVED_CURRENT_%d", years)
  } else if (identical(period, "future")) {
    sprintf("%s_%s_%d", toupper(scenario), climate_model, years)
  } else {
    stop("Unknown forcing period: ", period)
  }
}

# Function: concatenate annual scenario bundles without losing repeated years.
# The shared sequence helper uses year-based names, which are not unique when
# the observed 10-year cycle is repeated to make a 20-year sequence. Unique
# month-step names preserve every repeated month and retain the month prefix.
build_unique_sequence <- function(scenario_ids, scenario_type, sequence_label) {
  bundles <- lapply(
    scenario_ids,
    function(id) load_rvf_scenario(id, cfg = cfg, validate = TRUE)
  )
  reference <- bundles[[1]]

  for (bundle in bundles[-1]) {
    stopifnot(bundle$n_vil == reference$n_vil)
    stopifnot(bundle$n_gra == reference$n_gra)
    stopifnot(bundle$n_wat == reference$n_wat)
    stopifnot(identical(bundle$vil_names, reference$vil_names))
    stopifnot(identical(bundle$gra_names, reference$gra_names))
    stopifnot(identical(bundle$wat_names, reference$wat_names))
  }

  G_list <- list()
  W_list <- list()
  climate_vil <- list()
  climate_gra <- list()
  climate_wat <- list()
  sequence_year_index <- integer(0)
  sequence_month <- character(0)

  for (year_index in seq_along(bundles)) {
    bundle <- bundles[[year_index]]
    for (month in month.abb) {
      step_number <- length(G_list) + 1L
      key <- sprintf("%s_%03d", month, step_number)
      G_list[[key]] <- bundle$G_list[[month]]
      W_list[[key]] <- bundle$W_list[[month]]
      climate_vil[[key]] <- bundle$climate$vil[[month]]
      climate_gra[[key]] <- bundle$climate$gra[[month]]
      climate_wat[[key]] <- bundle$climate$wat[[month]]
      sequence_year_index <- c(sequence_year_index, year_index)
      sequence_month <- c(sequence_month, month)
    }
  }

  sequence <- list(
    G_list = G_list,
    W_list = W_list,
    climate = list(
      vil = climate_vil,
      gra = climate_gra,
      wat = climate_wat
    ),
    n_vil = reference$n_vil,
    n_gra = reference$n_gra,
    n_wat = reference$n_wat,
    vil_names = reference$vil_names,
    gra_names = reference$gra_names,
    wat_names = reference$wat_names,
    metadata = list(
      scenario_type = scenario_type,
      sequence_label = sequence_label,
      component_scenarios = scenario_ids,
      sequence_years = vapply(
        bundles,
        function(x) x$metadata$scenario_year,
        numeric(1)
      ),
      sequence_year_index = sequence_year_index,
      sequence_month = sequence_month,
      sequence_length_months = length(G_list),
      network_mode = "monthly"
    )
  )
  validate_rvf_scenario(sequence)
  sequence
}

# Function: combine climate from one sequence with movement from another.
# Both source sequences refer to the same sites and contain the same number of
# monthly steps; only their climate and G/W components differ.
build_mixed_forcing_scenario <- function(climate_ids,
                                         movement_ids,
                                         climate_period,
                                         movement_period,
                                         sequence_label) {
  climate_sequence <- build_unique_sequence(
    climate_ids,
    scenario_type = paste0(climate_period, "_climate_source"),
    sequence_label = paste0(sequence_label, "_climate")
  )
  movement_sequence <- build_unique_sequence(
    movement_ids,
    scenario_type = paste0(movement_period, "_movement_source"),
    sequence_label = paste0(sequence_label, "_movement")
  )

  if (length(climate_sequence$climate$vil) !=
      length(movement_sequence$G_list)) {
    stop("Climate and movement sequences contain different month counts.")
  }
  stopifnot(climate_sequence$n_vil == movement_sequence$n_vil)
  stopifnot(climate_sequence$n_gra == movement_sequence$n_gra)
  stopifnot(climate_sequence$n_wat == movement_sequence$n_wat)
  stopifnot(identical(climate_sequence$vil_names, movement_sequence$vil_names))
  stopifnot(identical(climate_sequence$gra_names, movement_sequence$gra_names))
  stopifnot(identical(climate_sequence$wat_names, movement_sequence$wat_names))

  mixed <- climate_sequence
  mixed$G_list <- movement_sequence$G_list
  mixed$W_list <- movement_sequence$W_list
  mixed$metadata <- list(
    scenario_type = "climate_movement_attribution",
    sequence_label = sequence_label,
    climate_period = climate_period,
    movement_period = movement_period,
    climate_component_scenarios = climate_ids,
    movement_component_scenarios = movement_ids,
    climate_sequence_years = climate_sequence$metadata$sequence_years,
    movement_sequence_years = movement_sequence$metadata$sequence_years,
    sequence_length_months = length(mixed$G_list),
    network_mode = "monthly"
  )
  validate_rvf_scenario(mixed)
  mixed
}

# -------- Complete SSP x GCM x window x replicate grid -------------------
scenario_index <- list_rvf_scenarios(cfg)
current_complete_replicates <- rvf_complete_replicates_for_years(
  index = scenario_index,
  dataset = "observed",
  years = sim$baseline_years,
  climate_model = "OBSERVED_CURRENT",
  network_replicate = network_replicates
)

future_cells <- crossing(
  scenario = future_scenarios,
  climate_model = climate_models,
  future_window = names(future_windows)
)
context_index <- bind_rows(lapply(seq_len(nrow(future_cells)), function(i) {
  cell <- future_cells[i, ]
  window_years <- future_windows[[cell$future_window]]
  future_complete_replicates <- rvf_complete_replicates_for_years(
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
    future_year = as.integer(cell$future_window),
    future_window = cell$future_window,
    window_start = min(window_years),
    window_end = max(window_years),
    window_years = paste(window_years, collapse = ","),
    network_replicate = network_replicates,
    forcing_complete =
      network_replicates %in% current_complete_replicates &
        network_replicates %in% future_complete_replicates
  )
})) |>
  arrange(
    .data$scenario,
    .data$climate_model,
    .data$future_year,
    .data$network_replicate
  )

missing_contexts <- context_index |>
  filter(!.data$forcing_complete)
if (nrow(missing_contexts) > 0L) {
  stop(
    "The attribution forcing grid is incomplete. Rerun Module 03 script 01."
  )
}
context_index <- context_index |>
  select(-forcing_complete) |>
  mutate(context_row = row_number())

expected_contexts <-
  length(future_scenarios) * length(climate_models) *
  length(future_windows) * length(network_replicates)
if (nrow(context_index) != expected_contexts) {
  stop("Unexpected number of complete attribution contexts.")
}

cell_grid <- tibble(
  attribution_cell = c("A", "B", "C", "D"),
  climate_period = c("current", "future", "current", "future"),
  movement_period = c("current", "current", "future", "future"),
  cell_label = c(
    "Current climate and movement",
    "Future climate only",
    "Future movement only",
    "Future climate and movement"
  )
)

# Cell A is independent of SSP, GCM, and future window. Run it once for each
# network replicate, then broadcast its summary to matched future contexts.
a_context_rows <- context_index |>
  group_by(.data$network_replicate) |>
  slice_min(.data$context_row, n = 1L, with_ties = FALSE) |>
  ungroup() |>
  pull("context_row")

attribution_jobs <- bind_rows(
  tibble(
    context_row = a_context_rows,
    attribution_cell = "A"
  ),
  crossing(
    context_row = context_index$context_row,
    attribution_cell = c("B", "C", "D")
  )
) |>
  mutate(job_id = row_number())

# Function: run and diagnose one factorial cell.
# job: one row of attribution_jobs, identifying a context and A/B/C/D cell.
run_attribution_job <- function(job) {
  context <- context_index |>
    filter(.data$context_row == job$context_row[[1]])
  cell <- cell_grid |>
    filter(.data$attribution_cell == job$attribution_cell[[1]])

  future_year_values <- as.integer(strsplit(
    context$window_years[[1]],
    ",",
    fixed = TRUE
  )[[1]])
  current_year_values <- rep(
    as.integer(sim$baseline_years),
    length.out = length(future_year_values)
  )
  climate_years <- if (cell$climate_period[[1]] == "current") {
    current_year_values
  } else {
    future_year_values
  }
  movement_years <- if (cell$movement_period[[1]] == "current") {
    current_year_values
  } else {
    future_year_values
  }

  replicate_id <- context$network_replicate[[1]]
  climate_base_ids <- base_ids_for_period(
    cell$climate_period[[1]],
    climate_years,
    context$climate_model[[1]],
    context$scenario[[1]]
  )
  movement_base_ids <- base_ids_for_period(
    cell$movement_period[[1]],
    movement_years,
    context$climate_model[[1]],
    context$scenario[[1]]
  )
  climate_ids <- rvf_with_replicate(climate_base_ids, replicate_id)
  movement_ids <- rvf_with_replicate(movement_base_ids, replicate_id)
  run_key <- sprintf(
    "%s_%s_%s_%s_REP%03d",
    toupper(context$scenario[[1]]),
    context$climate_model[[1]],
    context$future_window[[1]],
    cell$attribution_cell[[1]],
    replicate_id
  )

  scenario <- build_mixed_forcing_scenario(
    climate_ids = climate_ids,
    movement_ids = movement_ids,
    climate_period = cell$climate_period[[1]],
    movement_period = cell$movement_period[[1]],
    sequence_label = run_key
  )
  params <- rvf_assemble_params(
    scenario = scenario,
    calibrated = calibrated,
    use_temp_dependence = TRUE,
    p_mobile = sim$p_mobile,
    capacity_reference_p_mobile = sim$p_mobile,
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
  annual_for_diagnostics <- annual$system_year |>
    mutate(run_key = .env$run_key)
  time_for_diagnostics <- result$time_series |>
    select(
      year_since_seed,
      step_label,
      seroprevalence,
      prevalence,
      infectious_hosts
    ) |>
    mutate(run_key = .env$run_key)

  diagnostics <- rvf_build_cycle_diagnostics(
    annual_system = annual_for_diagnostics,
    time_series = time_for_diagnostics,
    group_cols = "run_key",
    forcing_cycle_years = forcing_cycle_years,
    years_endemic = years_endemic,
    terminal_cycles = terminal_cycles,
    relative_tolerance = relative_tolerance,
    seroprevalence_tolerance_pp = seroprevalence_tolerance_pp,
    required_replicate_fraction = 1,
    window_cycles = window_cycles
  )

  terminal <- diagnostics$window_sensitivity |>
    filter(.data$summary_cycles == terminal_cycles)
  if (nrow(terminal) != 1L) {
    stop("Expected one terminal summary for attribution job ", run_key)
  }

  longer_windows <- diagnostics$window_sensitivity |>
    filter(.data$summary_cycles > terminal_cycles)
  incidence_window_difference <- if (nrow(longer_windows) == 0L) {
    NA_real_
  } else {
    max(
      rvf_symmetric_relative_change(
        longer_windows$annual_new_infections,
        terminal$annual_new_infections
      ),
      na.rm = TRUE
    )
  }
  serology_window_difference <- if (nrow(longer_windows) == 0L) {
    NA_real_
  } else {
    max(
      abs(
        longer_windows$mean_seroprevalence_end -
          terminal$mean_seroprevalence_end
      ),
      na.rm = TRUE
    )
  }

  diagnostic_row <- diagnostics$convergence_by_run |>
    left_join(diagnostics$periodicity_by_run, by = "run_key") |>
    left_join(diagnostics$persistence_by_run, by = "run_key")

  metadata <- tibble(
    job_id = job$job_id[[1]],
    context_row = context$context_row[[1]],
    run_key = run_key,
    attribution_cell = cell$attribution_cell[[1]],
    cell_label = cell$cell_label[[1]],
    climate_period = cell$climate_period[[1]],
    movement_period = cell$movement_period[[1]],
    scenario = context$scenario[[1]],
    scenario_label = context$scenario_label[[1]],
    climate_model = context$climate_model[[1]],
    future_year = context$future_year[[1]],
    future_window = context$future_window[[1]],
    window_start = context$window_start[[1]],
    window_end = context$window_end[[1]],
    window_years = context$window_years[[1]],
    network_replicate = replicate_id,
    years_burn = years_burn,
    years_endemic = years_endemic,
    terminal_cycles = terminal_cycles,
    terminal_summary_years = terminal$summary_years,
    seed_mode = "all_villages",
    seed_per_village = seed_per_village,
    p_mobile = sim$p_mobile,
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

  persistence_class <- case_when(
    diagnostic_row$fraction_time_below_1e_6_infectious_hosts > 0 ~
      "near_deterministic_extinction",
    diagnostic_row$fraction_time_below_one_infectious_host > 0 ~
      "fractional_deterministic_persistence",
    TRUE ~ "persistent_above_one_host"
  )

  list(
    status = metadata |>
      select(
        job_id,
        context_row,
        attribution_cell,
        scenario,
        climate_model,
        future_window,
        network_replicate
      ) |>
      mutate(run_status = "ok", error_message = NA_character_),
    summary = metadata |>
      mutate(
        annual_new_infections_mean = terminal$annual_new_infections,
        annual_new_infections_per_1000 =
          terminal$annual_new_infections / .env$total_hosts * 1000,
        annual_mean_prevalence = terminal$annual_mean_prevalence,
        annual_mean_daily_foi = terminal$annual_mean_daily_foi,
        mean_seroprevalence = terminal$mean_seroprevalence_end,
        recurrent_state_supported = diagnostic_row$recurrent_state_supported,
        recurrent_pattern = diagnostic_row$recurrent_pattern,
        incidence_stable = diagnostic_row$incidence_stable,
        prevalence_stable = diagnostic_row$prevalence_stable,
        foi_stable = diagnostic_row$foi_stable,
        seroprevalence_stable = diagnostic_row$seroprevalence_stable,
        persistence_class = persistence_class,
        incidence_window_max_relative_difference =
          incidence_window_difference,
        seroprevalence_window_max_difference_pp =
          serology_window_difference
      )
  )
}

# -------- Run or reload the compact simulation cache ---------------------
n_cores <- rvf_max_cores(
  default = if (smoke_mode) 2L else cfg$execution$max_cores,
  n_tasks = nrow(attribution_jobs)
)
message(sprintf(
  paste0(
    "Attribution design: %d context(s), %d job(s), %d replicate(s), ",
    "%d years, %d core(s)."
  ),
  nrow(context_index),
  nrow(attribution_jobs),
  length(network_replicates),
  years_endemic,
  n_cores
))

if (reuse_cache) {
  if (!file.exists(output_paths[["cache"]])) {
    stop("Attribution-cache reuse is enabled in R/config.R, but no cache exists.")
  }
  cache <- readRDS(output_paths[["cache"]])
  if (is.null(cache$jobs) || is.null(cache$results) ||
      !isTRUE(all.equal(
        cache$jobs[, c("context_row", "attribution_cell", "job_id")],
        attribution_jobs[, c("context_row", "attribution_cell", "job_id")]
      ))) {
    stop("The saved attribution cache does not match the current job grid.")
  }
  results <- cache$results
} else {
  results <- rvf_parallel_lapply(seq_len(nrow(attribution_jobs)), function(i) {
    job <- attribution_jobs[i, ]
    context <- context_index |>
      filter(.data$context_row == job$context_row[[1]])
    message(sprintf(
      "Attribution %d/%d: %s %s %s cell %s replicate %03d",
      i,
      nrow(attribution_jobs),
      context$scenario[[1]],
      context$climate_model[[1]],
      context$future_window[[1]],
      job$attribution_cell[[1]],
      context$network_replicate[[1]]
    ))

    tryCatch(
      run_attribution_job(job),
      error = function(error) {
        list(
          status = tibble(
            job_id = job$job_id[[1]],
            context_row = context$context_row[[1]],
            attribution_cell = job$attribution_cell[[1]],
            scenario = context$scenario[[1]],
            climate_model = context$climate_model[[1]],
            future_window = context$future_window[[1]],
            network_replicate = context$network_replicate[[1]],
            run_status = "failed",
            error_message = conditionMessage(error)
          ),
          summary = tibble()
        )
      }
    )
  }, n_cores = n_cores)

  # The expensive output is compact because diagnostics are calculated within
  # each worker; weekly ODE trajectories are not retained in the cache.
  saveRDS(
    list(
      protocol = list(
        smoke_mode = smoke_mode,
        years_endemic = years_endemic,
        forcing_cycle_years = forcing_cycle_years,
        terminal_cycles = terminal_cycles,
        reporting_years = terminal_cycles * forcing_cycle_years,
        network_replicates = network_replicates,
        future_scenarios = future_scenarios,
        climate_models = climate_models,
        future_windows = names(future_windows),
        site_profile = sim$site_profile,
        rainfall_model = sim$rainfall_model
      ),
      contexts = context_index,
      jobs = attribution_jobs,
      results = results
    ),
    output_paths[["cache"]]
  )
}

run_status <- bind_rows(lapply(results, function(x) x$status))
write_csv(run_status, output_paths[["status"]])
failed_runs <- run_status |>
  filter(.data$run_status != "ok")
if (nrow(failed_runs) > 0L) {
  stop(
    nrow(failed_runs),
    " attribution job(s) failed. Inspect ",
    output_paths[["status"]]
  )
}

run_summaries <- bind_rows(lapply(results, function(x) x$summary))
if (nrow(run_summaries) != nrow(attribution_jobs)) {
  stop("Attribution results do not match the planned job grid.")
}

# -------- Broadcast cell A and validate the complete factorial -----------
context_columns <- c(
  "context_row",
  "scenario",
  "scenario_label",
  "climate_model",
  "future_year",
  "future_window",
  "window_start",
  "window_end",
  "window_years"
)
a_summaries <- run_summaries |>
  filter(.data$attribution_cell == "A") |>
  select(-job_id, -all_of(context_columns))
if (nrow(a_summaries) != length(network_replicates) ||
    anyDuplicated(a_summaries$network_replicate)) {
  stop("Expected one cell-A summary per network replicate.")
}
a_broadcast <- context_index |>
  left_join(a_summaries, by = "network_replicate", relationship = "many-to-one")
non_a_summaries <- run_summaries |>
  filter(.data$attribution_cell != "A") |>
  select(-job_id)

cell_replicates <- bind_rows(a_broadcast, non_a_summaries) |>
  arrange(
    .data$scenario,
    .data$climate_model,
    .data$future_year,
    .data$network_replicate,
    .data$attribution_cell
  )

factorial_check <- cell_replicates |>
  count(
    .data$scenario,
    .data$climate_model,
    .data$future_window,
    .data$network_replicate,
    name = "n_cells"
  )
if (nrow(cell_replicates) != 4L * nrow(context_index) ||
    any(factorial_check$n_cells != 4L)) {
  stop("At least one context lacks a complete A/B/C/D factorial.")
}

# Recurrent-state support is assessed across the configured paired network
# replicates within each SSP/GCM/window/cell. Near-extinction is a scientific
# outcome and is not treated as a convergence failure.
convergence_aggregate <- cell_replicates |>
  group_by(
    .data$scenario,
    .data$scenario_label,
    .data$climate_model,
    .data$future_year,
    .data$future_window,
    .data$window_start,
    .data$window_end,
    .data$attribution_cell,
    .data$cell_label
  ) |>
  summarise(
    n_network_replicates = n(),
    phase_matched_recurrence_supported_fraction = mean(
      .data$recurrent_state_supported,
      na.rm = TRUE
    ),
    incidence_stable_fraction = mean(.data$incidence_stable, na.rm = TRUE),
    prevalence_stable_fraction = mean(.data$prevalence_stable, na.rm = TRUE),
    foi_stable_fraction = mean(.data$foi_stable, na.rm = TRUE),
    seroprevalence_stable_fraction = mean(
      .data$seroprevalence_stable,
      na.rm = TRUE
    ),
    incidence_window_stable_fraction = mean(
      .data$incidence_window_max_relative_difference <=
        .env$relative_tolerance,
      na.rm = TRUE
    ),
    seroprevalence_window_stable_fraction = mean(
      .data$seroprevalence_window_max_difference_pp <=
        .env$seroprevalence_tolerance_pp,
      na.rm = TRUE
    ),
    incidence_window_max_relative_difference_mid = safe_median(
      .data$incidence_window_max_relative_difference
    ),
    incidence_window_max_relative_difference_hi = safe_p90(
      .data$incidence_window_max_relative_difference
    ),
    seroprevalence_window_max_difference_pp_mid = safe_median(
      .data$seroprevalence_window_max_difference_pp
    ),
    seroprevalence_window_max_difference_pp_hi = safe_p90(
      .data$seroprevalence_window_max_difference_pp
    ),
    near_extinction_fraction = mean(
      .data$persistence_class == "near_deterministic_extinction",
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  mutate(
    phase_matched_recurrence_supported =
      .data$phase_matched_recurrence_supported_fraction >=
        .env$required_replicate_fraction,
    cycle_averaged_burden_supported =
      .data$incidence_stable_fraction >=
        .env$required_replicate_fraction &
      .data$prevalence_stable_fraction >=
        .env$required_replicate_fraction &
      .data$foi_stable_fraction >=
        .env$required_replicate_fraction &
      .data$incidence_window_stable_fraction >=
        .env$required_replicate_fraction &
      .data$seroprevalence_window_stable_fraction >=
        .env$required_replicate_fraction,
    reporting_support_class = case_when(
      .data$phase_matched_recurrence_supported ~
        "phase_matched_recurrent",
      .data$cycle_averaged_burden_supported ~
        "phase_sensitive_bounded",
      TRUE ~ "finite_horizon_not_settled"
    )
  )
write_csv(convergence_aggregate, output_paths[["convergence"]])

if (!smoke_mode && any(!convergence_aggregate$cycle_averaged_burden_supported)) {
  stop(
    "At least one attribution cell lacks stable cycle-averaged burden. Inspect ",
    output_paths[["convergence"]],
    " before interpreting the factorial."
  )
}
phase_sensitive_cells <- convergence_aggregate |>
  filter(.data$reporting_support_class == "phase_sensitive_bounded")
if (nrow(phase_sensitive_cells) > 0L) {
  message(
    nrow(phase_sensitive_cells),
    " attribution cell(s) have stable cycle-averaged burden but not strict ",
    "year-by-year phase matching; see the convergence table."
  )
}

# -------- Cell estimates and matched decomposition -----------------------
cell_aggregate <- cell_replicates |>
  group_by(
    .data$scenario,
    .data$scenario_label,
    .data$climate_model,
    .data$future_year,
    .data$future_window,
    .data$window_start,
    .data$window_end,
    .data$window_years,
    .data$attribution_cell,
    .data$cell_label,
    .data$climate_period,
    .data$movement_period,
    .data$rainfall_model,
    .data$rainfall_label,
    .data$site_profile,
    .data$site_profile_label,
    .data$capacity_scale,
    .data$p_mobile,
    .data$terminal_summary_years
  ) |>
  summarise(
    n_network_replicates = n(),
    annual_new_infections_mean_mid = safe_median(
      .data$annual_new_infections_mean
    ),
    annual_new_infections_mean_lo = safe_p10(
      .data$annual_new_infections_mean
    ),
    annual_new_infections_mean_hi = safe_p90(
      .data$annual_new_infections_mean
    ),
    mean_seroprevalence_mid = safe_median(.data$mean_seroprevalence),
    mean_seroprevalence_lo = safe_p10(.data$mean_seroprevalence),
    mean_seroprevalence_hi = safe_p90(.data$mean_seroprevalence),
    .groups = "drop"
  )

decomposition_replicates <- cell_replicates |>
  select(
    scenario,
    scenario_label,
    climate_model,
    future_year,
    future_window,
    window_start,
    window_end,
    window_years,
    network_replicate,
    rainfall_model,
    rainfall_label,
    site_profile,
    site_profile_label,
    capacity_scale,
    p_mobile,
    terminal_summary_years,
    attribution_cell,
    annual_new_infections_mean
  ) |>
  pivot_wider(
    names_from = attribution_cell,
    values_from = annual_new_infections_mean
  ) |>
  mutate(
    total_change = .data$D - .data$A,
    climate_contribution = .data$B - .data$A,
    movement_contribution = .data$C - .data$A,
    interaction_contribution = .data$D - .data$B - .data$C + .data$A,
    decomposition_identity_error =
      .data$total_change -
      .data$climate_contribution -
      .data$movement_contribution -
      .data$interaction_contribution
  )

if (anyNA(decomposition_replicates[c("A", "B", "C", "D")])) {
  stop("At least one matched replicate lacks an A/B/C/D cell estimate.")
}
if (max(abs(decomposition_replicates$decomposition_identity_error)) > 1e-8) {
  stop("The replicate-level attribution decomposition is not additive.")
}

decomposition_aggregate <- decomposition_replicates |>
  group_by(
    .data$scenario,
    .data$scenario_label,
    .data$climate_model,
    .data$future_year,
    .data$future_window,
    .data$window_start,
    .data$window_end,
    .data$window_years,
    .data$rainfall_model,
    .data$rainfall_label,
    .data$site_profile,
    .data$site_profile_label,
    .data$capacity_scale,
    .data$p_mobile,
    .data$terminal_summary_years
  ) |>
  summarise(
    n_network_replicates = n(),
    A_mid = safe_median(.data$A),
    B_mid = safe_median(.data$B),
    C_mid = safe_median(.data$C),
    D_mid = safe_median(.data$D),
    total_change_mid = safe_median(.data$total_change),
    total_change_lo = safe_p10(.data$total_change),
    total_change_hi = safe_p90(.data$total_change),
    climate_contribution_mid = safe_median(.data$climate_contribution),
    climate_contribution_lo = safe_p10(.data$climate_contribution),
    climate_contribution_hi = safe_p90(.data$climate_contribution),
    movement_contribution_mid = safe_median(.data$movement_contribution),
    movement_contribution_lo = safe_p10(.data$movement_contribution),
    movement_contribution_hi = safe_p90(.data$movement_contribution),
    interaction_contribution_mid = safe_median(
      .data$interaction_contribution
    ),
    interaction_contribution_lo = safe_p10(.data$interaction_contribution),
    interaction_contribution_hi = safe_p90(.data$interaction_contribution),
    decomposition_identity_error_max = max(
      abs(.data$decomposition_identity_error)
    ),
    .groups = "drop"
  ) |>
  mutate(
    # Central estimates and intervals both summarise paired replicate-level
    # contrasts. Their medians need not add exactly because the median is
    # nonlinear; exact additivity is checked for every replicate above.
    central_additivity_gap =
      .data$total_change_mid -
      .data$climate_contribution_mid -
      .data$movement_contribution_mid -
      .data$interaction_contribution_mid,
    climate_share_of_total_mid = safe_ratio(
      .data$climate_contribution_mid,
      .data$total_change_mid
    ),
    movement_share_of_total_mid = safe_ratio(
      .data$movement_contribution_mid,
      .data$total_change_mid
    ),
    interaction_share_of_total_mid = safe_ratio(
      .data$interaction_contribution_mid,
      .data$total_change_mid
    ),
    dominant_component_mid = case_when(
      abs(.data$climate_contribution_mid) >=
        abs(.data$movement_contribution_mid) &
        abs(.data$climate_contribution_mid) >=
        abs(.data$interaction_contribution_mid) ~ "climate",
      abs(.data$movement_contribution_mid) >=
        abs(.data$interaction_contribution_mid) ~ "movement",
      TRUE ~ "interaction"
    ),
    metric = "annual_new_infections_mean"
  )

decomposition_window_summary <- decomposition_aggregate |>
  group_by(
    .data$scenario,
    .data$scenario_label,
    .data$future_year,
    .data$future_window,
    .data$window_start,
    .data$window_end,
    .data$rainfall_model,
    .data$rainfall_label,
    .data$site_profile,
    .data$site_profile_label,
    .data$capacity_scale,
    .data$p_mobile,
    .data$terminal_summary_years
  ) |>
  summarise(
    n_contexts = n_distinct(.data$climate_model),
    n_network_replicates = safe_median(.data$n_network_replicates),
    baseline_A_mid = safe_median(.data$A_mid),
    both_D_mid = safe_median(.data$D_mid),
    total_change_mid = safe_median(.data$total_change_mid),
    total_change_min = safe_min(.data$total_change_mid),
    total_change_max = safe_max(.data$total_change_mid),
    climate_contribution_mid = safe_median(.data$climate_contribution_mid),
    climate_contribution_min = safe_min(.data$climate_contribution_mid),
    climate_contribution_max = safe_max(.data$climate_contribution_mid),
    movement_contribution_mid = safe_median(.data$movement_contribution_mid),
    movement_contribution_min = safe_min(.data$movement_contribution_mid),
    movement_contribution_max = safe_max(.data$movement_contribution_mid),
    interaction_contribution_mid = safe_median(
      .data$interaction_contribution_mid
    ),
    interaction_contribution_min = safe_min(
      .data$interaction_contribution_mid
    ),
    interaction_contribution_max = safe_max(
      .data$interaction_contribution_mid
    ),
    climate_share_of_total_mid = safe_median(.data$climate_share_of_total_mid),
    movement_share_of_total_mid = safe_median(
      .data$movement_share_of_total_mid
    ),
    interaction_share_of_total_mid = safe_median(
      .data$interaction_share_of_total_mid
    ),
    .groups = "drop"
  ) |>
  arrange(.data$scenario, .data$future_year)

write_csv(cell_aggregate, output_paths[["cells"]])
write_csv(
  decomposition_replicates,
  output_paths[["decomposition_replicates"]]
)
write_csv(decomposition_aggregate, output_paths[["decomposition"]])
write_csv(decomposition_window_summary, output_paths[["window_summary"]])

# -------- Publication figures --------------------------------------------
cell_plot_data <- cell_aggregate |>
  mutate(
    scenario_label = factor(
      .data$scenario_label,
      levels = unname(cfg$scenario_labels[future_scenarios])
    ),
    future_window = factor(
      .data$future_window,
      levels = names(future_windows)
    ),
    cell_label = factor(
      .data$cell_label,
      levels = cell_grid$cell_label
    )
  )

p_cells <- ggplot(
  cell_plot_data,
  aes(
    x = cell_label,
    y = annual_new_infections_mean_mid,
    colour = climate_model
  )
) +
  geom_linerange(
    aes(
      ymin = annual_new_infections_mean_lo,
      ymax = annual_new_infections_mean_hi
    ),
    position = position_dodge(width = 0.60),
    linewidth = 0.60
  ) +
  geom_point(
    position = position_dodge(width = 0.60),
    size = 2.5
  ) +
  facet_grid(scenario_label ~ future_window) +
  coord_flip() +
  scale_colour_gcm() +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(
    x = NULL,
    y = "Mean annual cattle infections"
  ) +
  theme_rvf(base_size = 12.5) +
  theme(
    text = element_text(colour = "black"),
    axis.title = element_text(size = 13.5, colour = "black"),
    axis.text = element_text(size = 10.5, colour = "black"),
    strip.text = element_text(size = 11, face = "bold", colour = "black"),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10.5, colour = "black")
  )

# Build the component table explicitly so the plotted quantities and their
# paired-replicate intervals are easy to audit.
decomposition_plot_data <- bind_rows(
  decomposition_aggregate |>
    transmute(
      across(c(
        scenario,
        scenario_label,
        climate_model,
        future_window
      )),
      component = "Direct climate",
      mid = climate_contribution_mid,
      lo = climate_contribution_lo,
      hi = climate_contribution_hi
    ),
  decomposition_aggregate |>
    transmute(
      across(c(
        scenario,
        scenario_label,
        climate_model,
        future_window
      )),
      component = "Movement network",
      mid = movement_contribution_mid,
      lo = movement_contribution_lo,
      hi = movement_contribution_hi
    ),
  decomposition_aggregate |>
    transmute(
      across(c(
        scenario,
        scenario_label,
        climate_model,
        future_window
      )),
      component = "Interaction",
      mid = interaction_contribution_mid,
      lo = interaction_contribution_lo,
      hi = interaction_contribution_hi
    )
) |>
  mutate(
    scenario_label = factor(
      .data$scenario_label,
      levels = unname(cfg$scenario_labels[future_scenarios])
    ),
    future_window = factor(
      .data$future_window,
      levels = names(future_windows)
    ),
    component = factor(
      .data$component,
      levels = c("Direct climate", "Movement network", "Interaction")
    )
  )

p_decomposition <- ggplot(
  decomposition_plot_data,
  aes(x = component, y = mid, colour = climate_model)
) +
  geom_hline(yintercept = 0, linewidth = 0.65, colour = "black") +
  geom_linerange(
    aes(ymin = lo, ymax = hi),
    position = position_dodge(width = 0.60),
    linewidth = 0.60
  ) +
  geom_point(
    position = position_dodge(width = 0.60),
    size = 2.5
  ) +
  facet_grid(scenario_label ~ future_window) +
  scale_colour_gcm() +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(
    x = NULL,
    y = "Contribution to annual cattle infections"
  ) +
  theme_rvf(base_size = 12.5) +
  theme(
    text = element_text(colour = "black"),
    axis.title = element_text(size = 13.5, colour = "black"),
    axis.text = element_text(size = 10.5, colour = "black"),
    strip.text = element_text(size = 11, face = "bold", colour = "black"),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10.5, colour = "black")
  )

# The raw A/B/C/D cell figure is retained for the primary exposure profile as
# a supplementary check. Script 10 builds the cross-profile decomposition.
if (!smoke_mode && is_primary_profile) {
  rvf_ggsave(
    file.path(cfg$figures_dir, "05_attribution_cells.png"),
    p_cells,
    width = 11.5,
    height = 7.4
  )
}

message("Climate-movement attribution outputs written to: ", out_dir)
if (smoke_mode) {
  message("Smoke test complete. Do not interpret its finite-horizon estimates.")
}
