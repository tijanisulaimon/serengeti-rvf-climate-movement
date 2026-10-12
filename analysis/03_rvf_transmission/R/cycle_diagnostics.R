# Shared diagnostics for periodically forced endemic simulations.
#
# The observed sequence repeats every 10 years and each future window every
# 20 years. These systems need not approach a constant equilibrium. The
# functions below therefore compare complete forcing cycles and save compact
# summaries rather than relying on a terminal year or full weekly trajectories.

rvf_symmetric_relative_change <- function(current, previous) {
  denominator <- (abs(current) + abs(previous)) / 2
  ifelse(
    denominator <= sqrt(.Machine$double.eps),
    0,
    abs(current - previous) / denominator
  )
}

rvf_max_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

# Function: summarise convergence and seasonality over complete forcing cycles.
# annual_system: one row per simulated year and analysis unit.
# time_series: within-year host trajectory for the same analysis units.
# group_cols: columns that uniquely identify one simulation run.
# terminal_cycles: complete cycles in each comparison block and reporting window.
rvf_build_cycle_diagnostics <- function(annual_system,
                                        time_series,
                                        group_cols,
                                        forcing_cycle_years,
                                        years_endemic,
                                        terminal_cycles = 2L,
                                        relative_tolerance = 0.10,
                                        seroprevalence_tolerance_pp = 0.50,
                                        required_replicate_fraction = 0.80,
                                        window_cycles = c(2L, 4L, 6L)) {
  forcing_cycle_years <- as.integer(forcing_cycle_years)
  years_endemic <- as.integer(years_endemic)
  terminal_cycles <- as.integer(terminal_cycles)

  if (!is.finite(forcing_cycle_years) || forcing_cycle_years < 1L) {
    stop("forcing_cycle_years must be a positive integer.")
  }
  if (!is.finite(years_endemic) || years_endemic < forcing_cycle_years) {
    stop("years_endemic must contain at least one complete forcing cycle.")
  }
  if (!is.finite(terminal_cycles) || terminal_cycles < 1L) {
    stop("terminal_cycles must be a positive integer.")
  }

  required_annual <- c(
    group_cols, "analysis_year", "annual_new_infections",
    "annual_mean_prevalence", "annual_mean_daily_foi"
  )
  required_time <- c(
    group_cols, "year_since_seed", "step_label", "seroprevalence",
    "prevalence", "infectious_hosts"
  )
  missing_annual <- setdiff(required_annual, names(annual_system))
  missing_time <- setdiff(required_time, names(time_series))
  if (length(missing_annual) > 0L) {
    stop("annual_system is missing: ", paste(missing_annual, collapse = ", "))
  }
  if (length(missing_time) > 0L) {
    stop("time_series is missing: ", paste(missing_time, collapse = ", "))
  }

  complete_cycles <- years_endemic %/% forcing_cycle_years
  complete_years <- complete_cycles * forcing_cycle_years
  if (complete_cycles < 2L * terminal_cycles) {
    stop(
      "Need at least ", 2L * terminal_cycles,
      " complete forcing cycles for the block convergence check."
    )
  }

  # Annual boundaries are included explicitly in every ODE output. Using the
  # boundary values avoids weighting artefacts from the additional monthly
  # forcing boundaries in the weekly output grid.
  serology_year <- time_series |>
    dplyr::filter(
      .data$year_since_seed > 0,
      abs(.data$year_since_seed - round(.data$year_since_seed)) < 1e-7
    ) |>
    dplyr::mutate(analysis_year = as.integer(round(.data$year_since_seed))) |>
    dplyr::filter(.data$analysis_year <= complete_years) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, "analysis_year")))) |>
    dplyr::summarise(
      seroprevalence_end = mean(.data$seroprevalence, na.rm = TRUE),
      infectious_hosts_end = mean(.data$infectious_hosts, na.rm = TRUE),
      .groups = "drop"
    )

  annual_with_cycle <- annual_system |>
    dplyr::filter(
      .data$analysis_year >= 1L,
      .data$analysis_year <= complete_years
    ) |>
    dplyr::left_join(
      serology_year,
      by = c(stats::setNames(group_cols, group_cols), "analysis_year")
    ) |>
    dplyr::mutate(
      forcing_cycle = (.data$analysis_year - 1L) %/% forcing_cycle_years + 1L,
      year_in_cycle = (.data$analysis_year - 1L) %% forcing_cycle_years + 1L
    )

  cycle_summary <- annual_with_cycle |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, "forcing_cycle")))) |>
    dplyr::summarise(
      forcing_cycle_years = .env$forcing_cycle_years,
      cycle_start_analysis_year = min(.data$analysis_year),
      cycle_end_analysis_year = max(.data$analysis_year),
      annual_new_infections = mean(.data$annual_new_infections, na.rm = TRUE),
      annual_mean_prevalence = mean(.data$annual_mean_prevalence, na.rm = TRUE),
      annual_mean_daily_foi = mean(.data$annual_mean_daily_foi, na.rm = TRUE),
      mean_seroprevalence_end = mean(.data$seroprevalence_end, na.rm = TRUE),
      minimum_infectious_hosts_end = min(.data$infectious_hosts_end, na.rm = TRUE),
      .groups = "drop"
    )

  previous_cycles <- seq.int(
    complete_cycles - 2L * terminal_cycles + 1L,
    complete_cycles - terminal_cycles
  )
  terminal_cycle_ids <- seq.int(
    complete_cycles - terminal_cycles + 1L,
    complete_cycles
  )

  # Average replicate cycles within each comparison block at matching forcing
  # phases. Taking a single mean over the entire block could conceal a rise in
  # one climate year that is offset by a decline in another.
  block_phase_summary <- annual_with_cycle |>
    dplyr::mutate(
      comparison_block = dplyr::case_when(
        .data$forcing_cycle %in% previous_cycles ~ "previous",
        .data$forcing_cycle %in% terminal_cycle_ids ~ "terminal",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(.data$comparison_block)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(
      group_cols, "comparison_block", "year_in_cycle"
    )))) |>
    dplyr::summarise(
      annual_new_infections = mean(.data$annual_new_infections, na.rm = TRUE),
      annual_mean_prevalence = mean(.data$annual_mean_prevalence, na.rm = TRUE),
      annual_mean_daily_foi = mean(.data$annual_mean_daily_foi, na.rm = TRUE),
      seroprevalence_end = mean(.data$seroprevalence_end, na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      names_from = "comparison_block",
      values_from = c(
        "annual_new_infections", "annual_mean_prevalence",
        "annual_mean_daily_foi", "seroprevalence_end"
      ),
      names_glue = "{.value}_{comparison_block}"
    )

  convergence_by_run <- block_phase_summary |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n_phase_years = dplyr::n(),
      annual_new_infections_previous = mean(
        .data$annual_new_infections_previous,
        na.rm = TRUE
      ),
      annual_new_infections_terminal = mean(
        .data$annual_new_infections_terminal,
        na.rm = TRUE
      ),
      annual_mean_prevalence_previous = mean(
        .data$annual_mean_prevalence_previous,
        na.rm = TRUE
      ),
      annual_mean_prevalence_terminal = mean(
        .data$annual_mean_prevalence_terminal,
        na.rm = TRUE
      ),
      annual_mean_daily_foi_previous = mean(
        .data$annual_mean_daily_foi_previous,
        na.rm = TRUE
      ),
      annual_mean_daily_foi_terminal = mean(
        .data$annual_mean_daily_foi_terminal,
        na.rm = TRUE
      ),
      mean_seroprevalence_end_previous = mean(
        .data$seroprevalence_end_previous,
        na.rm = TRUE
      ),
      mean_seroprevalence_end_terminal = mean(
        .data$seroprevalence_end_terminal,
        na.rm = TRUE
      ),
      incidence_abs_symmetric_relative_change = rvf_max_or_na(
        rvf_symmetric_relative_change(
          .data$annual_new_infections_terminal,
          .data$annual_new_infections_previous
        )
      ),
      prevalence_abs_symmetric_relative_change = rvf_max_or_na(
        rvf_symmetric_relative_change(
          .data$annual_mean_prevalence_terminal,
          .data$annual_mean_prevalence_previous
        )
      ),
      foi_abs_symmetric_relative_change = rvf_max_or_na(
        rvf_symmetric_relative_change(
          .data$annual_mean_daily_foi_terminal,
          .data$annual_mean_daily_foi_previous
        )
      ),
      seroprevalence_abs_difference_pp = rvf_max_or_na(abs(
        .data$seroprevalence_end_terminal -
          .data$seroprevalence_end_previous
      )),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      forcing_cycle_years = forcing_cycle_years,
      terminal_cycles = terminal_cycles,
      comparison_window_years = terminal_cycles * forcing_cycle_years,
      previous_cycle_start = min(previous_cycles),
      previous_cycle_end = max(previous_cycles),
      terminal_cycle_start = min(terminal_cycle_ids),
      terminal_cycle_end = max(terminal_cycle_ids),
      incidence_stable = .data$incidence_abs_symmetric_relative_change <= relative_tolerance,
      prevalence_stable = .data$prevalence_abs_symmetric_relative_change <= relative_tolerance,
      foi_stable = .data$foi_abs_symmetric_relative_change <= relative_tolerance,
      seroprevalence_stable =
        .data$seroprevalence_abs_difference_pp <= seroprevalence_tolerance_pp,
      recurrent_state_supported =
        .data$incidence_stable & .data$prevalence_stable &
        .data$foi_stable & .data$seroprevalence_stable
    )

  aggregate_cols <- setdiff(group_cols, "network_replicate")
  convergence_summary <- convergence_by_run |>
    dplyr::group_by(dplyr::across(dplyr::all_of(aggregate_cols))) |>
    dplyr::summarise(
      n_network_replicates = dplyr::n(),
      relative_tolerance = relative_tolerance,
      seroprevalence_tolerance_pp = seroprevalence_tolerance_pp,
      required_replicate_fraction = required_replicate_fraction,
      incidence_stable_fraction = mean(.data$incidence_stable, na.rm = TRUE),
      prevalence_stable_fraction = mean(.data$prevalence_stable, na.rm = TRUE),
      foi_stable_fraction = mean(.data$foi_stable, na.rm = TRUE),
      seroprevalence_stable_fraction = mean(.data$seroprevalence_stable, na.rm = TRUE),
      recurrent_state_supported_fraction =
        mean(.data$recurrent_state_supported, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      recurrent_state_supported =
        .data$incidence_stable_fraction >= required_replicate_fraction &
        .data$prevalence_stable_fraction >= required_replicate_fraction &
        .data$foi_stable_fraction >= required_replicate_fraction &
        .data$seroprevalence_stable_fraction >= required_replicate_fraction
    )

  compare_cycles <- function(current_cycle, previous_cycle, label) {
    metric_cols <- c(
      "annual_new_infections", "annual_mean_prevalence",
      "annual_mean_daily_foi"
    )
    current <- annual_with_cycle |>
      dplyr::filter(.data$forcing_cycle == current_cycle) |>
      dplyr::select(dplyr::all_of(c(
        group_cols, "year_in_cycle", metric_cols, "seroprevalence_end"
      ))) |>
      dplyr::rename_with(
        ~ paste0(.x, "_current"),
        dplyr::all_of(c(metric_cols, "seroprevalence_end"))
      )
    previous <- annual_with_cycle |>
      dplyr::filter(.data$forcing_cycle == previous_cycle) |>
      dplyr::select(dplyr::all_of(c(
        group_cols, "year_in_cycle", metric_cols, "seroprevalence_end"
      ))) |>
      dplyr::rename_with(
        ~ paste0(.x, "_previous"),
        dplyr::all_of(c(metric_cols, "seroprevalence_end"))
      )
    dplyr::inner_join(
      current,
      previous,
      by = c(group_cols, "year_in_cycle")
    ) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
      dplyr::summarise(
        !!paste0(label, "_max_relative_change") := pmax(
          rvf_max_or_na(rvf_symmetric_relative_change(
              .data$annual_new_infections_current,
              .data$annual_new_infections_previous
            )),
          rvf_max_or_na(rvf_symmetric_relative_change(
              .data$annual_mean_prevalence_current,
              .data$annual_mean_prevalence_previous
            )),
          rvf_max_or_na(rvf_symmetric_relative_change(
              .data$annual_mean_daily_foi_current,
              .data$annual_mean_daily_foi_previous
            )),
          na.rm = TRUE
        ),
        !!paste0(label, "_seroprevalence_difference_pp") := rvf_max_or_na(abs(
          .data$seroprevalence_end_current -
            .data$seroprevalence_end_previous
        )),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        !!paste0(label, "_supported") :=
          .data[[paste0(label, "_max_relative_change")]] <= relative_tolerance &
          .data[[paste0(label, "_seroprevalence_difference_pp")]] <=
            seroprevalence_tolerance_pp
      )
  }

  # A periodically forced system may repeat after one forcing cycle or alternate
  # between two cycles. Diagnose both possibilities before describing a run as
  # having reached a recurrent state.
  period_one <- compare_cycles(
    complete_cycles,
    complete_cycles - 1L,
    "period_one"
  )
  period_two_a <- compare_cycles(
    complete_cycles,
    complete_cycles - 2L,
    "period_two_a"
  )
  period_two_b <- compare_cycles(
    complete_cycles - 1L,
    complete_cycles - 3L,
    "period_two_b"
  )
  periodicity_by_run <- period_one |>
    dplyr::left_join(period_two_a, by = group_cols) |>
    dplyr::left_join(period_two_b, by = group_cols) |>
    dplyr::mutate(
      period_two_supported =
        .data$period_two_a_supported & .data$period_two_b_supported,
      recurrent_pattern = dplyr::case_when(
        .data$period_one_supported ~ "period_one",
        .data$period_two_supported ~ "period_two",
        TRUE ~ "not_resolved"
      )
    )

  # Recalculate headline metrics over several complete terminal windows. This
  # exposes conclusions that depend on an arbitrary reporting-window length.
  window_cycles <- sort(unique(as.integer(window_cycles)))
  window_cycles <- window_cycles[
    is.finite(window_cycles) & window_cycles >= 1L & window_cycles <= complete_cycles
  ]
  window_sensitivity <- dplyr::bind_rows(lapply(window_cycles, function(n_cycles) {
    cycle_summary |>
      dplyr::filter(.data$forcing_cycle > complete_cycles - n_cycles) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
      dplyr::summarise(
        summary_cycles = .env$n_cycles,
        summary_years = .env$n_cycles * .env$forcing_cycle_years,
        annual_new_infections = mean(.data$annual_new_infections, na.rm = TRUE),
        annual_mean_prevalence = mean(.data$annual_mean_prevalence, na.rm = TRUE),
        annual_mean_daily_foi = mean(.data$annual_mean_daily_foi, na.rm = TRUE),
        mean_seroprevalence_end = mean(.data$mean_seroprevalence_end, na.rm = TRUE),
        .groups = "drop"
      )
  }))

  terminal_window_years <- terminal_cycles * forcing_cycle_years
  terminal_window_start <- years_endemic - terminal_window_years
  # ODE output includes both regular reporting times and monthly forcing
  # boundaries. Duration weights prevent densely sampled months from receiving
  # greater influence in the settled seasonal summaries.
  weighted_time_series <- time_series |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::arrange(.data$year_since_seed, .by_group = TRUE) |>
    dplyr::mutate(
      dt_days = 365 * (dplyr::lead(.data$year_since_seed) - .data$year_since_seed)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(
      .data$year_since_seed >= terminal_window_start,
      .data$year_since_seed < years_endemic,
      is.finite(.data$dt_days),
      .data$dt_days > 0
    ) |>
    dplyr::mutate(month = substr(.data$step_label, 1L, 3L))

  seasonal_by_run <- weighted_time_series |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, "month")))) |>
    dplyr::summarise(
      comparison_window_years = .env$terminal_window_years,
      mean_infectious_hosts = stats::weighted.mean(
        .data$infectious_hosts, .data$dt_days, na.rm = TRUE
      ),
      mean_prevalence = stats::weighted.mean(
        .data$prevalence, .data$dt_days, na.rm = TRUE
      ),
      mean_seroprevalence = stats::weighted.mean(
        .data$seroprevalence, .data$dt_days, na.rm = TRUE
      ),
      .groups = "drop"
    )

  # These thresholds distinguish continuous deterministic persistence from
  # trajectories that approach practical extinction during part of the cycle.
  # They are diagnostics only and do not alter the ODE state.
  persistence_by_run <- weighted_time_series |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      comparison_window_years = .env$terminal_window_years,
      minimum_infectious_hosts = min(.data$infectious_hosts, na.rm = TRUE),
      fraction_time_below_one_infectious_host =
        sum(.data$dt_days[.data$infectious_hosts < 1], na.rm = TRUE) /
        sum(.data$dt_days, na.rm = TRUE),
      fraction_time_below_1e_6_infectious_hosts =
        sum(.data$dt_days[.data$infectious_hosts < 1e-6], na.rm = TRUE) /
        sum(.data$dt_days, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    protocol = list(
      forcing_cycle_years = forcing_cycle_years,
      years_endemic = years_endemic,
      complete_cycles = complete_cycles,
      terminal_cycles = terminal_cycles,
      terminal_window_years = terminal_window_years,
      relative_tolerance = relative_tolerance,
      seroprevalence_tolerance_pp = seroprevalence_tolerance_pp,
      required_replicate_fraction = required_replicate_fraction
    ),
    cycle_summary = cycle_summary,
    convergence_by_run = convergence_by_run,
    convergence_summary = convergence_summary,
    periodicity_by_run = periodicity_by_run,
    window_sensitivity = window_sensitivity,
    seasonal_by_run = seasonal_by_run,
    persistence_by_run = persistence_by_run
  )
}
