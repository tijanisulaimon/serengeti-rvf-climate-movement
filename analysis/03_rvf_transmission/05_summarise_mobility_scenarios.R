#!/usr/bin/env Rscript

#' 05_summarise_mobility_scenarios.R
#' -------------------------------------------------------------------------
#' Combine the five independently calibrated mobility x exposure-profile runs
#' produced by script 04. This script validates the complete result grid before
#' writing any figures, so incomplete or unsettled runs cannot silently enter
#' the manuscript outputs.
#'
#' Outputs are deliberately limited to:
#'   1. one profile-by-mobility analysis table;
#'   2. one convergence and persistence audit table;
#'   3. one five-profile interpretation table;
#'   4. a two-profile headline figure showing the direction reversal;
#'   5. a five-profile supplementary figure and cycle-settling diagnostic.
#' -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

source(here::here("analysis", "03_rvf_transmission", "R", "transmission_model.R"))
source(here::here("analysis", "03_rvf_transmission", "R", "plot_helpers.R"))

cfg <- load_rvf_config()
out_dir <- file.path(cfg$results_dir, "mobility_capacity_profile_comparison")
figures_dir <- cfg$figures_dir

# Order profiles from predominantly village exposure towards increasingly
# resource-associated exposure. This order is used in tables and legends.
profile_order <- c(
  "peridomestic_dominant",
  "balanced",
  "grazing_floodplain_dominant",
  "resource_edge_dominant",
  "waterpoint_extreme"
)
profile_labels <- c(
  peridomestic_dominant = "Village-associated",
  balanced = "Balanced",
  grazing_floodplain_dominant = "Grazing/floodplain",
  resource_edge_dominant = "Resource-edge",
  waterpoint_extreme = "Waterpoint-extreme"
)
profile_colours <- c(
  peridomestic_dominant = "#9C443C",
  balanced = "#666666",
  grazing_floodplain_dominant = "#177E7A",
  resource_edge_dominant = "#3B5BA5",
  waterpoint_extreme = "#B17800"
)
profile_shapes <- c(
  peridomestic_dominant = 16,
  balanced = 17,
  grazing_floodplain_dominant = 15,
  resource_edge_dominant = 18,
  waterpoint_extreme = 8
)

expected_p_mobile <- c(0.1, 0.25, 0.5, 0.75, 1) # Values run in script 04.
expected_result_dirs <- setNames(
  file.path(
    cfg$results_dir,
    paste0("mobility_capacity_cold_start_", profile_order, "_pmobile_capacity")
  ),
  profile_order
)
expected_intermediate_dirs <- setNames(
  file.path(
    cfg$intermediate_dir,
    paste0("mobility_capacity_cold_start_", profile_order, "_pmobile_capacity")
  ),
  profile_order
)

# Function: fail early when a required table field is unavailable.
assert_columns <- function(data, required, object_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      object_name,
      " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

# Read the three compact products from one profile-specific mobility run.
read_profile_outputs <- function(profile, result_directory, intermediate_directory) {
  paths <- c(
    aggregate = file.path(result_directory, "mobility_capacity_aggregate.csv"),
    convergence = file.path(result_directory, "mobility_capacity_convergence.csv"),
    runs = file.path(intermediate_directory, "mobility_capacity_runs.rds")
  )
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop(
      "Profile ", profile, " is incomplete. Missing: ",
      paste(basename(missing), collapse = ", "),
      ". Allow script 04 to finish before running script 05."
    )
  }

  runs <- readRDS(paths[["runs"]])
  cycle_trajectory <- runs$cycle_trajectory
  if (is.null(cycle_trajectory) || nrow(cycle_trajectory) == 0L) {
    stop("Profile ", profile, " has no saved cycle trajectory.")
  }

  list(
    aggregate = read_csv(paths[["aggregate"]], show_col_types = FALSE),
    convergence = read_csv(paths[["convergence"]], show_col_types = FALSE),
    cycle_trajectory = cycle_trajectory
  )
}

profile_results <- lapply(
  profile_order,
  function(profile) {
    read_profile_outputs(
      profile,
      expected_result_dirs[[profile]],
      expected_intermediate_dirs[[profile]]
    )
  }
)
names(profile_results) <- profile_order

# Stack the same product from each independently calibrated exposure profile.
aggregate_all <- bind_rows(lapply(profile_results, function(x) x$aggregate))
convergence_all <- bind_rows(lapply(profile_results, function(x) x$convergence))
cycle_trajectory_all <- bind_rows(
  lapply(profile_results, function(x) x$cycle_trajectory)
)

# Profile outputs generated before the extinction-boundary fix may label a
# near-zero cell as unsettled because relative changes are unstable around
# zero. Preserve the failed recurrence diagnostic, but classify the biological
# outcome from the independent persistence check.
aggregate_all <- aggregate_all |>
  mutate(
    outcome_class = if_else(
      !is.na(.data$deterministic_persistence_supported) &
        !.data$deterministic_persistence_supported,
      "near_deterministic_extinction",
      .data$outcome_class
    )
  )
convergence_all <- convergence_all |>
  mutate(
    outcome_class = if_else(
      .data$run_complete & !.data$deterministic_persistence_supported,
      "near_deterministic_extinction",
      .data$outcome_class
    )
  )

aggregate_required <- c(
  "site_profile", "site_profile_label", "p_mobile", "capacity_scale",
  "terminal_summary_years",
  "effective_share_village_routine", "effective_share_grazing_routine",
  "effective_share_watering_routine", "n_network_replicates", "n_planned",
  "n_successful", "success_fraction",
  "annual_new_infections_terminal_window_mid",
  "annual_new_infections_terminal_window_lo",
  "annual_new_infections_terminal_window_hi",
  "relative_annual_new_infections_terminal_window_mid",
  "relative_annual_new_infections_terminal_window_lo",
  "relative_annual_new_infections_terminal_window_hi",
  "mean_seroprevalence_terminal_window_mid",
  "mean_seroprevalence_terminal_window_lo",
  "mean_seroprevalence_terminal_window_hi",
  "quasi_equilibrium_supported", "deterministic_persistence_supported",
  "headline_supported", "outcome_class"
)
convergence_required <- c(
  "site_profile", "p_mobile", "n_planned", "n_successful",
  "success_fraction", "incidence_stable_fraction",
  "prevalence_stable_fraction", "foi_stable_fraction",
  "seroprevalence_stable_fraction", "quasi_equilibrium_supported",
  "deterministic_persistence_supported", "persistent_above_one_host_fraction",
  "window_sensitivity_incidence_max_relative_difference_mid",
  "window_sensitivity_incidence_max_relative_difference_hi",
  "window_sensitivity_seroprevalence_max_difference_pp_mid",
  "window_sensitivity_seroprevalence_max_difference_pp_hi",
  "headline_supported", "outcome_class"
)
assert_columns(aggregate_all, aggregate_required, "Mobility aggregate")
assert_columns(convergence_all, convergence_required, "Mobility diagnostics")

# -------- Check that every planned scenario is available -----------------
expected_grid <- tidyr::crossing(
  site_profile = profile_order,
  p_mobile = expected_p_mobile
)
observed_grid <- aggregate_all |>
  distinct(.data$site_profile, .data$p_mobile)
missing_grid <- anti_join(
  expected_grid,
  observed_grid,
  by = c("site_profile", "p_mobile")
)
extra_grid <- anti_join(
  observed_grid,
  expected_grid,
  by = c("site_profile", "p_mobile")
)
duplicate_grid <- aggregate_all |>
  count(.data$site_profile, .data$p_mobile, name = "n") |>
  filter(.data$n != 1L)

if (nrow(missing_grid) > 0L ||
    nrow(extra_grid) > 0L ||
    nrow(duplicate_grid) > 0L) {
  stop("The profile-by-mobility grid is incomplete or duplicated.")
}
expected_replicates <- length(cfg$analysis_network_replicates)
if (anyNA(aggregate_all$n_planned) ||
    anyNA(aggregate_all$n_successful) ||
    anyNA(aggregate_all$success_fraction) ||
    any(aggregate_all$n_planned != expected_replicates) ||
    any(aggregate_all$n_successful != expected_replicates) ||
    any(aggregate_all$success_fraction < 1)) {
  stop(
    "All ", expected_replicates,
    " configured network replicates must succeed in every publication cell."
  )
}
unsettled_persistent_cells <- aggregate_all |>
  filter(
    is.na(.data$quasi_equilibrium_supported) |
      (!.data$quasi_equilibrium_supported &
        .data$outcome_class != "near_deterministic_extinction")
  )
if (nrow(unsettled_persistent_cells) > 0L) {
  stop(
    "At least one non-extinct cell has not reached a supported recurrent state. ",
    "Inspect the diagnostics and extend that run before plotting."
  )
}

# A low-mobility cell may genuinely approach deterministic extinction. This is
# retained as a categorical result, but other unsettled outcomes are rejected.
invalid_outcomes <- setdiff(
  unique(aggregate_all$outcome_class),
  c(
    "recurrent_persistent_above_one_host",
    "recurrent_fractional_persistence",
    "near_deterministic_extinction"
  )
)
if (length(invalid_outcomes) > 0L) {
  stop("Unsupported outcome class(es): ", paste(invalid_outcomes, collapse = ", "))
}

reference_check <- aggregate_all |>
  filter(abs(.data$p_mobile - 1) < 1e-12)
reference_ratios <-
  reference_check$relative_annual_new_infections_terminal_window_mid
if (anyNA(reference_ratios) || any(abs(reference_ratios - 1) > 1e-8)) {
  stop("The paired p_mobile = 1 reference does not equal one for every profile.")
}

if (any(
  convergence_all$window_sensitivity_incidence_max_relative_difference_mid >
    0.10,
  na.rm = TRUE
)) {
  warning(
    "At least one 2-cycle estimate differs by more than 10% from its longer ",
    "4/6-cycle averages. Inspect the diagnostic table.",
    call. = FALSE
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

plot_data <- aggregate_all |>
  mutate(
    site_profile = factor(.data$site_profile, levels = profile_order),
    profile_label = factor(
      unname(profile_labels[as.character(.data$site_profile)]),
      levels = unname(profile_labels[profile_order])
    ),
    recurrent_endemic_state =
      .data$outcome_class != "near_deterministic_extinction"
  ) |>
  arrange(.data$site_profile, .data$p_mobile)

# Keep one compact table containing the estimates that support manuscript
# statements. Full diagnostic objects remain in each profile's RDS file.
comparison_table <- plot_data |>
  select(
    site_profile,
    profile_label,
    p_mobile,
    capacity_scale,
    terminal_summary_years,
    effective_share_village_routine,
    effective_share_grazing_routine,
    effective_share_watering_routine,
    n_network_replicates,
    annual_new_infections_terminal_window_mid,
    annual_new_infections_terminal_window_lo,
    annual_new_infections_terminal_window_hi,
    relative_annual_new_infections_terminal_window_mid,
    relative_annual_new_infections_terminal_window_lo,
    relative_annual_new_infections_terminal_window_hi,
    mean_seroprevalence_terminal_window_mid,
    mean_seroprevalence_terminal_window_lo,
    mean_seroprevalence_terminal_window_hi,
    quasi_equilibrium_supported,
    deterministic_persistence_supported,
    outcome_class
  )

diagnostic_table <- convergence_all |>
  mutate(
    site_profile = factor(.data$site_profile, levels = profile_order),
    profile_label = unname(profile_labels[as.character(.data$site_profile)])
  ) |>
  select(
    site_profile,
    profile_label,
    p_mobile,
    n_planned,
    n_successful,
    success_fraction,
    incidence_stable_fraction,
    prevalence_stable_fraction,
    foi_stable_fraction,
    seroprevalence_stable_fraction,
    quasi_equilibrium_supported,
    deterministic_persistence_supported,
    persistent_above_one_host_fraction,
    window_sensitivity_incidence_max_relative_difference_mid,
    window_sensitivity_incidence_max_relative_difference_hi,
    window_sensitivity_seroprevalence_max_difference_pp_mid,
    window_sensitivity_seroprevalence_max_difference_pp_hi,
    headline_supported,
    outcome_class
  ) |>
  arrange(.data$site_profile, .data$p_mobile)

# Describe the direction of the median response. The tolerance below only
# prevents solver-scale numerical noise from being counted as a real step;
# it is not an epidemiological effect threshold.
direction_table <- plot_data |>
  group_by(.data$site_profile, .data$profile_label) |>
  arrange(.data$p_mobile, .by_group = TRUE) |>
  summarise(
    capacity_scale = first(.data$capacity_scale),
    effective_share_village_at_p1 =
      last(.data$effective_share_village_routine),
    effective_share_grazing_at_p1 =
      last(.data$effective_share_grazing_routine),
    effective_share_watering_at_p1 =
      last(.data$effective_share_watering_routine),
    annual_infections_p01 = first(
      .data$annual_new_infections_terminal_window_mid
    ),
    annual_infections_p1 = last(
      .data$annual_new_infections_terminal_window_mid
    ),
    relative_infections_p01 = first(
      .data$relative_annual_new_infections_terminal_window_mid
    ),
    n_increasing_steps = sum(
      diff(.data$annual_new_infections_terminal_window_mid) >
        1e-8 * max(.data$annual_new_infections_terminal_window_mid),
      na.rm = TRUE
    ),
    n_decreasing_steps = sum(
      diff(.data$annual_new_infections_terminal_window_mid) <
        -1e-8 * max(.data$annual_new_infections_terminal_window_mid),
      na.rm = TRUE
    ),
    n_near_extinction_cells = sum(!.data$recurrent_endemic_state),
    .groups = "drop"
  ) |>
  mutate(
    response_direction = case_when(
      .data$n_increasing_steps > 0 & .data$n_decreasing_steps == 0 ~
        "increases_with_mobility",
      .data$n_decreasing_steps > 0 & .data$n_increasing_steps == 0 ~
        "decreases_with_mobility",
      .data$n_increasing_steps == 0 & .data$n_decreasing_steps == 0 ~
        "little_change",
      TRUE ~ "non_monotonic"
    )
  ) |>
  arrange(.data$site_profile)

write_csv(
  comparison_table,
  file.path(out_dir, "mobility_profile_summary.csv")
)
write_csv(
  diagnostic_table,
  file.path(out_dir, "mobility_profile_diagnostics.csv")
)
write_csv(
  direction_table,
  file.path(out_dir, "mobility_profile_interpretation.csv")
)

# -------- Mobility figures: burden relative to full mobility -------------
# Each estimate is paired with p_mobile = 1 from the same network replicate
# in script 04. A linear axis keeps zero visible and uses regular one-fold tick
# intervals. Lines join only the five evaluated values; they do not represent
# a fitted continuous response.
active_profiles <- profile_order[
  profile_order %in% as.character(unique(plot_data$site_profile))
]
relative_upper <- ceiling(max(
  plot_data$relative_annual_new_infections_terminal_window_hi,
  na.rm = TRUE
))

# The main figure isolates the two exposure geographies needed to demonstrate
# the direction reversal. The five-profile version below shows that the result
# is not an artefact of selecting only the two endpoints of the profile set.
headline_profiles <- c(
  "peridomestic_dominant",
  "grazing_floodplain_dominant"
)
headline_labels <- c(
  peridomestic_dominant = "Village-associated exposure",
  grazing_floodplain_dominant = "Shared grazing/floodplain exposure"
)
headline_data <- plot_data |>
  filter(as.character(.data$site_profile) %in% headline_profiles) |>
  mutate(
    site_profile = factor(
      as.character(.data$site_profile),
      levels = headline_profiles
    )
  )

p_relative_headline <- ggplot(
  headline_data,
  aes(
    x = p_mobile,
    y = relative_annual_new_infections_terminal_window_mid,
    colour = site_profile,
    shape = site_profile,
    group = site_profile
  )
) +
  geom_hline(yintercept = 1, linewidth = 0.7, colour = "black") +
  geom_errorbar(
    aes(
      ymin = relative_annual_new_infections_terminal_window_lo,
      ymax = relative_annual_new_infections_terminal_window_hi
    ),
    width = 0.018,
    linewidth = 0.65,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 3.2, stroke = 0.8) +
  scale_x_continuous(
    breaks = expected_p_mobile,
    limits = c(0.08, 1.02)
  ) +
  scale_y_continuous(
    breaks = seq.int(0, relative_upper, by = 1),
    labels = paste0(seq.int(0, relative_upper, by = 1), "x"),
    limits = c(0, relative_upper),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_colour_manual(
    name = NULL,
    values = profile_colours[headline_profiles],
    breaks = headline_profiles,
    labels = headline_labels[headline_profiles]
  ) +
  scale_shape_manual(
    name = NULL,
    values = profile_shapes[headline_profiles],
    breaks = headline_profiles,
    labels = headline_labels[headline_profiles]
  ) +
  guides(
    colour = guide_legend(nrow = 1, byrow = TRUE),
    shape = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  labs(
    x = expression("Mobile-herd fraction (" * p[mobile] * ")"),
    y = expression("Annual infections relative to " * p[mobile] * " = 1")
  ) +
  theme_rvf(base_size = 13) +
  theme(
    text = element_text(colour = "black"),
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 11, colour = "black"),
    legend.key.width = grid::unit(1.15, "cm")
  )

rvf_ggsave(
  file.path(
    figures_dir,
    "03_pmobile_contrasting_profiles_relative_incidence.png"
  ),
  p_relative_headline,
  width = 8.8,
  height = 5.6
)

# All five calibrated profiles are retained as a supplementary sensitivity
# analysis. This version preserves intermediate and extreme exposure settings.
p_relative_all_profiles <- ggplot(
  plot_data,
  aes(
    x = p_mobile,
    y = relative_annual_new_infections_terminal_window_mid,
    colour = site_profile,
    shape = site_profile,
    group = site_profile
  )
) +
  geom_hline(yintercept = 1, linewidth = 0.7, colour = "black") +
  geom_errorbar(
    aes(
      ymin = relative_annual_new_infections_terminal_window_lo,
      ymax = relative_annual_new_infections_terminal_window_hi
    ),
    width = 0.018,
    linewidth = 0.65,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.9, stroke = 0.8) +
  scale_x_continuous(
    breaks = expected_p_mobile,
    limits = c(0.08, 1.02)
  ) +
  scale_y_continuous(
    breaks = seq.int(0, relative_upper, by = 1),
    labels = paste0(seq.int(0, relative_upper, by = 1), "x"),
    limits = c(0, relative_upper),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_colour_manual(
    name = NULL,
    values = profile_colours,
    breaks = active_profiles,
    labels = profile_labels[active_profiles]
  ) +
  scale_shape_manual(
    name = NULL,
    values = profile_shapes,
    breaks = active_profiles,
    labels = profile_labels[active_profiles]
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  labs(
    x = expression(p[mobile]),
    y = expression("Annual infections relative to " * p[mobile] * " = 1")
  ) +
  theme_rvf(base_size = 13) +
  theme(
    text = element_text(colour = "black"),
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 11, colour = "black"),
    legend.key.width = grid::unit(1.15, "cm")
  )

rvf_ggsave(
  file.path(figures_dir, "03_pmobile_all_profiles_relative_incidence.png"),
  p_relative_all_profiles,
  width = 8.8,
  height = 5.8
)

# -------- Diagnostic figure: settling across repeated forcing cycles ------
# Only the two endpoint mobility values are shown. Intermediate values add
# clutter without improving the assessment of whether the contrasting systems
# approach recurrent trajectories.
cycle_plot_data <- cycle_trajectory_all |>
  filter(.data$p_mobile %in% c(0.1, 1)) |>
  mutate(
    site_profile = factor(.data$site_profile, levels = profile_order),
    profile_label = factor(
      unname(profile_labels[as.character(.data$site_profile)]),
      levels = unname(profile_labels[profile_order])
    ),
    mobility_label = factor(
      if_else(
        abs(.data$p_mobile - 0.1) < 1e-12,
        "p = 0.1",
        "p = 1.0"
      ),
      levels = c("p = 0.1", "p = 1.0")
    )
  )

endpoint_colours <- c("p = 0.1" = "#6B6B6B", "p = 1.0" = "#1F6E9C")
p_cycles <- ggplot(
  cycle_plot_data,
  aes(
    x = cycle_index,
    y = annual_new_infections_mid,
    colour = mobility_label,
    fill = mobility_label,
    group = mobility_label
  )
) +
  geom_ribbon(
    aes(ymin = annual_new_infections_lo, ymax = annual_new_infections_hi),
    colour = NA,
    alpha = 0.15,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.85) +
  facet_wrap(vars(profile_label), ncol = 3, scales = "free_y") +
  scale_colour_manual(values = endpoint_colours) +
  scale_fill_manual(values = endpoint_colours) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 6)) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(
    x = "Repeated 10-year forcing cycle",
    y = "Mean annual infections"
  ) +
  theme_rvf(base_size = 11.5) +
  theme(
    text = element_text(colour = "black"),
    axis.title = element_text(size = 12.5, colour = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    strip.text = element_text(size = 11, face = "bold", colour = "black"),
    legend.text = element_text(size = 11, colour = "black")
  )

if (rvf_write_diagnostic_plots()) {
  diagnostic_dir <- file.path(figures_dir, "diagnostics")
  rvf_ggsave(
    file.path(diagnostic_dir, "mobility_profile_cycle_diagnostic.png"),
    p_cycles,
    width = 9,
    height = 6.5
  )
}

message("Validated all five mobility-profile runs.")
message("Wrote comparison tables to: ", out_dir)
message("Wrote manuscript figures to: ", figures_dir)
