#!/usr/bin/env Rscript
#' Collate Module 03 headline evidence for manuscript integration.
#'
#' This script does not run the ODE. It reads the current Q1--Q4 aggregate
#' outputs and writes compact summary tables plus a short Markdown status note.
#' The aim is to make the academic claims easy to audit without digging through
#' every long-run output folder.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tibble)
})

source(here::here("analysis", "03_rvf_transmission", "R", "transmission_model.R"))
source(here::here("analysis", "03_rvf_transmission", "R", "plot_helpers.R"))

cfg <- load_rvf_config()
out_dir <- file.path(cfg$results_dir, "headline_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

q1_path <- rvf_calibration_recommended_path(cfg)
q2_path <- file.path(
  cfg$results_dir,
  "future_cold_endemic",
  "future_cold_endemic_window_burden_aggregate.csv"
)
q3_path <- file.path(
  cfg$results_dir,
  "mobility_capacity_profile_comparison",
  "mobility_profile_interpretation.csv"
)
q3_summary_path <- file.path(
  cfg$results_dir,
  "mobility_capacity_profile_comparison",
  "mobility_profile_summary.csv"
)
q4_primary_path <- file.path(
  cfg$results_dir,
  "climate_movement_attribution_cold_start",
  "climate_movement_attribution_decomposition_aggregate.csv"
)
q4_primary_window_path <- file.path(
  cfg$results_dir,
  "climate_movement_attribution_cold_start",
  "climate_movement_attribution_decomposition_window_summary.csv"
)
q4_contrast_path <- file.path(
  cfg$results_dir,
  "climate_movement_attribution_cold_start_profile_peridomestic_dominant",
  "climate_movement_attribution_decomposition_aggregate.csv"
)
q4_contrast_window_path <- file.path(
  cfg$results_dir,
  "climate_movement_attribution_cold_start_profile_peridomestic_dominant",
  "climate_movement_attribution_decomposition_window_summary.csv"
)

input_list <- tibble(
  question = c(
    "Q1", "Q2", "Q3", "Q3-summary", "Q4-primary", "Q4-primary-window",
    "Q4-contrast", "Q4-contrast-window"
  ),
  role = c(
    "calibration profile gate",
    "matched-window current/future burden",
    "mobility x vector-pressure profile interpretation",
    "mobility replicate and reporting-period check",
    "primary-profile climate/movement/interaction attribution",
    "primary-profile window-level attribution summary",
    "contrasting-profile climate/movement/interaction attribution",
    "contrasting-profile window-level attribution summary"
  ),
  path = c(
    q1_path, q2_path, q3_path, q3_summary_path,
    q4_primary_path, q4_primary_window_path,
    q4_contrast_path, q4_contrast_window_path
  ),
  exists = file.exists(path)
)
write_csv(input_list, file.path(out_dir, "module03_headline_inputs.csv"))

missing_inputs <- input_list |> filter(!.data$exists)
if (nrow(missing_inputs) > 0L) {
  stop(
    "Missing required headline input(s):\n",
    paste0(" - ", missing_inputs$question, ": ", missing_inputs$path, collapse = "\n")
  )
}

fmt_num <- function(x, digits = 1) {
  ifelse(
    is.na(x),
    "NA",
    formatC(x, format = "f", digits = digits, big.mark = ",")
  )
}

fmt_ratio <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", paste0(formatC(x, format = "f", digits = digits), "x"))
}

fmt_range <- function(x, digits = 0) {
  rng <- range(x, na.rm = TRUE)
  if (length(unique(rng)) == 1L) {
    fmt_num(rng[1], digits)
  } else {
    paste0(fmt_num(rng[1], digits), "--", fmt_num(rng[2], digits))
  }
}

q1 <- read_csv(q1_path, show_col_types = FALSE)
q2 <- read_csv(q2_path, show_col_types = FALSE)
q3 <- read_csv(q3_path, show_col_types = FALSE)
q3_runs <- read_csv(q3_summary_path, show_col_types = FALSE)
q4 <- bind_rows(
  read_csv(q4_primary_path, show_col_types = FALSE),
  read_csv(q4_contrast_path, show_col_types = FALSE)
)
q4_window <- bind_rows(
  read_csv(q4_primary_window_path, show_col_types = FALSE),
  read_csv(q4_contrast_window_path, show_col_types = FALSE)
)

# Refuse to turn exploratory 5- or 10-replicate outputs into headline evidence.
# The manuscript builder should only consume the configured publication grid and
# the common phase-matched reporting period.
expected_replicates <- length(unique(cfg$analysis_network_replicates))
expected_reporting_years <- as.integer(cfg$execution$reporting_years)

assert_protocol <- function(data, label, require_window = TRUE) {
  if (!"n_network_replicates" %in% names(data) ||
      any(data$n_network_replicates != expected_replicates)) {
    stop(label, " does not use ", expected_replicates, " network replicates.")
  }
  if (require_window &&
      (!"terminal_summary_years" %in% names(data) ||
       any(data$terminal_summary_years != expected_reporting_years))) {
    stop(label, " does not use the common ", expected_reporting_years,
         "-year reporting period.")
  }
}

assert_protocol(q1, "Q1 calibration", require_window = FALSE)
if (!"calibration_window_years" %in% names(q1) ||
    any(q1$calibration_window_years != expected_reporting_years)) {
  stop("Q1 calibration does not use the common reporting period.")
}
assert_protocol(q2, "Q2 future burden", require_window = FALSE)
if (!"window_length_years" %in% names(q2) ||
    any(q2$window_length_years != expected_reporting_years)) {
  stop("Q2 future burden does not use the common reporting period.")
}
assert_protocol(q3_runs, "Q3 mobility")
assert_protocol(q4, "Q4 attribution")

q1_selected <- q1
if ("rainfall_model" %in% names(q1_selected)) {
  configured_rainfall <- cfg$sim$rainfall_model %||% NA_character_
  rainfall_models <- sort(unique(q1_selected$rainfall_model))
  if (configured_rainfall %in% rainfall_models) {
    q1_selected <- q1_selected |>
      filter(.data$rainfall_model == configured_rainfall)
  } else if (length(rainfall_models) > 1L || length(rainfall_models) == 1L) {
    stop(
      "Q1 calibration table does not contain the configured rainfall model '",
      configured_rainfall,
      "'. Available rainfall_model values: ",
      paste(rainfall_models, collapse = ", "),
      ". Rerun script 02 for the configured rainfall model."
    )
  }
} else {
  q1_selected$rainfall_model <- NA_character_
}
if (!"rainfall_label" %in% names(q1_selected)) {
  q1_selected$rainfall_label <- NA_character_
}

q1_summary <- q1_selected |>
  transmute(
    site_profile,
    site_profile_label,
    rainfall_model,
    rainfall_label,
    manuscript_role = case_when(
      .data$site_profile == "grazing_floodplain_dominant" ~ "primary biological profile",
      .data$site_profile == "peridomestic_dominant" ~ "opposing biological contrast",
      .data$site_profile %in% c("balanced", "resource_edge_dominant") ~
        "intermediate sensitivity profile",
      .data$site_profile == "waterpoint_extreme" ~ "stress-test sensitivity profile",
      TRUE ~ "supporting profile"
    ),
    capacity_scale,
    target_band,
    mean_seroprevalence_calibration,
    annual_new_infections_terminal_window,
    n_network_replicates,
    n_success,
    profile_decision
  ) |>
  arrange(match(
    .data$site_profile,
    c(
      "grazing_floodplain_dominant",
      "peridomestic_dominant",
      "balanced",
      "resource_edge_dominant",
      "waterpoint_extreme"
    )
  ))

q2_summary <- q2 |>
  transmute(
    scenario_type,
    scenario,
    scenario_label,
    climate_model,
    future_window,
    window_start,
    window_end,
    window_length_years,
    site_profile,
    site_profile_label,
    annual_new_infections_mean_mid = annual_new_infections_mid,
    annual_new_infections_mean_lo = annual_new_infections_lo,
    annual_new_infections_mean_hi = annual_new_infections_hi,
    annual_new_infections_per_1000_mean_mid = annual_new_infections_per_1000_mid,
    annual_mean_infectious_hosts_mean_mid = annual_mean_infectious_hosts_mid,
    mean_seroprevalence_window_mid = mean_seroprevalence_mid,
    annual_new_infections_mean_ratio_vs_baseline =
      annual_new_infections_ratio_vs_baseline_mid
  ) |>
  arrange(
    desc(.data$scenario_type == "current_baseline"),
    desc(.data$annual_new_infections_mean_mid)
  )

q2_baseline <- q2_summary |> filter(.data$scenario_type == "current_baseline") |> slice(1)
q2_future <- q2_summary |> filter(.data$scenario_type != "current_baseline")
q2_top <- q2_future |> slice_max(.data$annual_new_infections_mean_mid, n = 1, with_ties = FALSE)

q3_summary <- q3 |>
  transmute(
    site_profile,
    profile_label,
    manuscript_role = case_when(
      .data$site_profile == "grazing_floodplain_dominant" ~
        "primary biological profile",
      .data$site_profile == "peridomestic_dominant" ~
        "opposing biological contrast",
      TRUE ~ "supporting sensitivity profile"
    ),
    capacity_scale,
    qualitative_effect = response_direction,
    annual_new_infections_p01 = annual_infections_p01,
    annual_new_infections_p1 = annual_infections_p1,
    relative_annual_infections_p01 = relative_infections_p01,
    n_increasing_steps,
    n_decreasing_steps,
    n_near_extinction_cells
  ) |>
  arrange(match(
    .data$site_profile,
    c(
      "grazing_floodplain_dominant",
      "peridomestic_dominant",
      "balanced",
      "resource_edge_dominant",
      "waterpoint_extreme"
    )
  ))

q4_summary <- q4 |>
  filter(.data$metric %in% c(
    "annual_new_infections_mean",
    "annual_mean_infectious_hosts_mean",
    "mean_annual_foi_probability_mean"
  )) |>
  transmute(
    metric,
    scenario,
    scenario_label,
    climate_model,
    future_window,
    window_start,
    window_end,
    site_profile,
    site_profile_label,
    capacity_scale,
    A_current_climate_current_movement = A_mid,
    B_future_climate_current_movement = B_mid,
    C_current_climate_future_movement = C_mid,
    D_future_climate_future_movement = D_mid,
    total_change_mid,
    climate_contribution_mid,
    movement_contribution_mid,
    interaction_contribution_mid,
    climate_share_of_total_mid,
    movement_share_of_total_mid,
    interaction_share_of_total_mid,
    dominant_component_mid,
    n_network_replicates
  ) |>
  arrange(.data$metric, .data$scenario, .data$future_window, .data$climate_model)

q4_incidence_all <- q4_summary |>
  filter(.data$metric == "annual_new_infections_mean")

q4_window_summary <- q4_window |>
  transmute(
    scenario,
    scenario_label,
    future_window,
    window_start,
    window_end,
    site_profile,
    site_profile_label,
    capacity_scale,
    n_contexts,
    n_network_replicates,
    baseline_A_mid,
    both_D_mid,
    total_change_mid,
    climate_contribution_mid,
    movement_contribution_mid,
    interaction_contribution_mid,
    climate_share_of_total_mid,
    movement_share_of_total_mid,
    interaction_share_of_total_mid
  ) |>
  arrange(.data$site_profile, .data$scenario, .data$future_window)

# One attribution context is one profile x SSP x GCM x future-window
# comparison. Each profile contributes 24 contexts (2 x 4 x 3).
q4_n_contexts <- nrow(dplyr::distinct(
  q4_incidence_all,
  site_profile,
  scenario,
  climate_model,
  future_window
))
q4_n_profiles <- n_distinct(q4_incidence_all$site_profile)
q4_n_scenarios <- n_distinct(q4_incidence_all$scenario)
q4_n_models <- n_distinct(q4_incidence_all$climate_model)
q4_n_windows <- n_distinct(q4_incidence_all$future_window)
q4_n_climate_dominant <- sum(q4_incidence_all$dominant_component_mid == "climate", na.rm = TRUE)

# Summarise the four GCMs without treating them as a random sample. The point
# is the ensemble median and the line spans the full four-model range.
q4_component_shares <- bind_rows(
  q4_incidence_all |>
    transmute(
      scenario,
      scenario_label,
      future_window,
      climate_model,
      site_profile,
      site_profile_label,
      component = "Direct climate",
      contribution_share = climate_share_of_total_mid
    ),
  q4_incidence_all |>
    transmute(
      scenario,
      scenario_label,
      future_window,
      climate_model,
      site_profile,
      site_profile_label,
      component = "Movement network",
      contribution_share = movement_share_of_total_mid
    ),
  q4_incidence_all |>
    transmute(
      scenario,
      scenario_label,
      future_window,
      climate_model,
      site_profile,
      site_profile_label,
      component = "Interaction",
      contribution_share = interaction_share_of_total_mid
    )
)

q4_profile_plot_summary <- q4_component_shares |>
  group_by(
    .data$scenario,
    .data$scenario_label,
    .data$future_window,
    .data$site_profile,
    .data$site_profile_label,
    .data$component
  ) |>
  summarise(
    n_gcms = n_distinct(.data$climate_model),
    contribution_percent_mid = 100 * median(
      .data$contribution_share,
      na.rm = TRUE
    ),
    contribution_percent_min = 100 * min(
      .data$contribution_share,
      na.rm = TRUE
    ),
    contribution_percent_max = 100 * max(
      .data$contribution_share,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  mutate(
    ssp_label = rvf_ssp_label(.data$scenario),
    future_window = factor(
      .data$future_window,
      levels = names(cfg$future_windows)
    ),
    profile_panel = case_when(
      .data$site_profile == "grazing_floodplain_dominant" ~
        "Grazing/floodplain-dominant",
      .data$site_profile == "peridomestic_dominant" ~
        "Peridomestic-dominant",
      TRUE ~ .data$site_profile_label
    ),
    profile_panel = factor(
      .data$profile_panel,
      levels = c(
        "Grazing/floodplain-dominant",
        "Peridomestic-dominant"
      )
    ),
    component = factor(
      .data$component,
      levels = c("Direct climate", "Movement network", "Interaction")
    )
  )

# Plot absolute contributions because percentages cluster near 100% for direct
# climate and zero for the other terms. A shared-scale panel shows dominance;
# a second panel reveals the sign and GCM range of the smaller components.
q4_component_absolute <- bind_rows(
  q4_incidence_all |>
    transmute(
      scenario,
      future_window,
      climate_model,
      site_profile,
      site_profile_label,
      component = "Direct climate",
      contribution = climate_contribution_mid
    ),
  q4_incidence_all |>
    transmute(
      scenario,
      future_window,
      climate_model,
      site_profile,
      site_profile_label,
      component = "Movement network",
      contribution = movement_contribution_mid
    ),
  q4_incidence_all |>
    transmute(
      scenario,
      future_window,
      climate_model,
      site_profile,
      site_profile_label,
      component = "Interaction",
      contribution = interaction_contribution_mid
    )
)

q4_absolute_plot_summary <- q4_component_absolute |>
  group_by(
    .data$scenario,
    .data$future_window,
    .data$site_profile,
    .data$site_profile_label,
    .data$component
  ) |>
  summarise(
    n_gcms = n_distinct(.data$climate_model),
    contribution_mid = median(.data$contribution, na.rm = TRUE),
    contribution_min = min(.data$contribution, na.rm = TRUE),
    contribution_max = max(.data$contribution, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    ssp_label = rvf_ssp_label(.data$scenario),
    future_window = factor(
      .data$future_window,
      levels = names(cfg$future_windows)
    ),
    profile_panel = case_when(
      .data$site_profile == "grazing_floodplain_dominant" ~
        "Shared grazing/floodplain",
      .data$site_profile == "peridomestic_dominant" ~
        "Village-associated",
      TRUE ~ .data$site_profile_label
    ),
    profile_panel = factor(
      .data$profile_panel,
      levels = c(
        "Shared grazing/floodplain",
        "Village-associated"
      )
    ),
    component = factor(
      .data$component,
      levels = c("Direct climate", "Movement network", "Interaction")
    )
  )

attribution_colours <- c(
  "Direct climate" = "#D55E00",
  "Movement network" = "#0072B2",
  "Interaction" = "#009E73"
)

attribution_theme <- theme_rvf(base_size = 12) +
  theme(
    text = element_text(colour = "black"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 9.5, colour = "black"),
    strip.text = element_text(size = 10, face = "bold", colour = "black"),
    legend.text = element_text(size = 10, colour = "black"),
    legend.position = "bottom",
    plot.margin = margin(5.5, 16, 5.5, 5.5)
  )

p_attribution_full <- q4_absolute_plot_summary |>
  filter(.data$component == "Direct climate") |>
  ggplot(
  aes(
    x = future_window,
    y = contribution_mid,
    colour = component,
    group = component
  )
) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.55) +
  geom_linerange(
    aes(
      ymin = contribution_min,
      ymax = contribution_max
    ),
    position = position_dodge(width = 0.55),
    linewidth = 0.70
  ) +
  geom_point(
    position = position_dodge(width = 0.55),
    size = 2.6
  ) +
  facet_grid(profile_panel ~ ssp_label) +
  scale_colour_manual(
    values = attribution_colours,
    limits = names(attribution_colours),
    drop = FALSE
  ) +
  scale_x_discrete(labels = c(
    "2050" = "2041-2060",
    "2080" = "2061-2080",
    "2100" = "2081-2100"
  )) +
  labs(
    x = NULL,
    y = "Change in annual cattle infections",
    tag = "a"
  ) +
  attribution_theme +
  guides(colour = "none")

p_attribution_zoom <- q4_absolute_plot_summary |>
  filter(.data$component != "Direct climate") |>
  ggplot(
    aes(
      x = future_window,
      y = contribution_mid,
      colour = component,
      group = component
    )
  ) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.55) +
  geom_linerange(
    aes(ymin = contribution_min, ymax = contribution_max),
    position = position_dodge(width = 0.45),
    linewidth = 0.70
  ) +
  geom_point(position = position_dodge(width = 0.45), size = 2.5) +
  facet_grid(profile_panel ~ ssp_label) +
  scale_colour_manual(
    values = attribution_colours,
    limits = c("Movement network", "Interaction"),
    drop = TRUE
  ) +
  scale_x_discrete(labels = c(
    "2050" = "2041-2060",
    "2080" = "2061-2080",
    "2100" = "2081-2100"
  )) +
  labs(
    x = "Future climate window",
    y = "Change in annual cattle infections",
    tag = "b"
  ) +
  attribution_theme

# Panel a isolates the dominant direct-climate contribution. Panel b uses an
# expanded scale for the smaller movement-network and interaction terms.
p_attribution_profiles <- p_attribution_full / p_attribution_zoom +
  plot_layout(heights = c(1, 1.05), guides = "collect")

question_status <- tibble(
  question = c("Q1", "Q2", "Q3", "Q4"),
  headline_evidence = c(
    sprintf(
      "%d/%d calibrated exposure profiles pass the 4--10%% serology gate.",
      sum(q1_summary$target_band %in% TRUE, na.rm = TRUE),
      nrow(q1_summary)
    ),
    sprintf(
      "Highest future annual infection burden is %s %s at %s infections/year (%s baseline).",
      q2_top$climate_model,
      q2_top$future_window,
      fmt_num(q2_top$annual_new_infections_mean_mid, 0),
      fmt_ratio(q2_top$annual_new_infections_mean_ratio_vs_baseline, 2)
    ),
    paste(
      "Profile effects span",
      paste(sort(unique(q3_summary$qualitative_effect)), collapse = ", "),
      "under the positive p_mobile headline grid."
    ),
    sprintf(
      "Across %d attribution contexts (%d profiles x %d SSP x %d GCM x %d windows), climate dominates annual-infection change in %d/%d (climate %s vs movement %s and interaction %s).",
      q4_n_contexts,
      q4_n_profiles,
      q4_n_scenarios,
      q4_n_models,
      q4_n_windows,
      q4_n_climate_dominant,
      q4_n_contexts,
      fmt_range(q4_incidence_all$climate_contribution_mid, 0),
      fmt_range(q4_incidence_all$movement_contribution_mid, 0),
      fmt_range(q4_incidence_all$interaction_contribution_mid, 0)
    )
  ),
  primary_table = c(
    "module03_q1_calibration_summary.csv",
    "module03_q2_window_burden_summary.csv",
    "module03_q3_mobility_profile_summary.csv",
    "module03_q4_attribution_summary.csv"
  )
)

write_csv(q1_summary, file.path(out_dir, "module03_q1_calibration_summary.csv"))
write_csv(q2_summary, file.path(out_dir, "module03_q2_window_burden_summary.csv"))
write_csv(q3_summary, file.path(out_dir, "module03_q3_mobility_profile_summary.csv"))
write_csv(q4_summary, file.path(out_dir, "module03_q4_attribution_summary.csv"))
write_csv(q4_window_summary, file.path(out_dir, "module03_q4_attribution_window_summary.csv"))
write_csv(
  q4_profile_plot_summary,
  file.path(out_dir, "module03_q4_attribution_profile_window_summary.csv")
)
write_csv(
  q4_absolute_plot_summary,
  file.path(out_dir, "module03_q4_attribution_absolute_window_summary.csv")
)
write_csv(question_status, file.path(out_dir, "module03_question_status.csv"))
rvf_ggsave(
  file.path(cfg$figures_dir, "05_attribution_decomposition.png"),
  p_attribution_profiles,
  width = 10.5,
  height = 10.5
)

md <- c(
  "# Module 03 Headline Summary",
  "",
  "## Inputs",
  "",
  paste0("- ", input_list$question, ": ", input_list$path),
  "",
  "## Q1 Calibration Gate",
  "",
  sprintf(
    "- %d of %d effective exposure profiles pass the 4--10%% serology anchor after tuning only `capacity_scale`.",
    sum(q1_summary$target_band %in% TRUE, na.rm = TRUE),
    nrow(q1_summary)
  ),
  paste0(
    "- ",
    q1_summary$site_profile,
    ": c_K = ",
    fmt_num(q1_summary$capacity_scale, 2),
    ", mean serology = ",
    fmt_num(q1_summary$mean_seroprevalence_calibration, 2),
    "%, target band = ",
    q1_summary$target_band
  ),
  "",
  "## Q2 Current/Future Burden",
  "",
  sprintf(
    "- Current observed-cycle baseline: %s annual infections/year, %s infections per 1,000 cattle/year.",
    fmt_num(q2_baseline$annual_new_infections_mean_mid, 0),
    fmt_num(q2_baseline$annual_new_infections_per_1000_mean_mid, 1)
  ),
  sprintf(
    "- Highest future window: %s %s, %s annual infections/year, %s relative to the current baseline.",
    q2_top$climate_model,
    q2_top$future_window,
    fmt_num(q2_top$annual_new_infections_mean_mid, 0),
    fmt_ratio(q2_top$annual_new_infections_mean_ratio_vs_baseline, 2)
  ),
  "",
  "## Q3 Mobility x Vector-Pressure Geography",
  "",
  paste0(
    "- ",
    q3_summary$site_profile,
    ": ",
    q3_summary$qualitative_effect,
    "; annual infections at p_mobile = 0.1 and 1 are ",
    fmt_num(q3_summary$annual_new_infections_p01, 0),
    " and ",
    fmt_num(q3_summary$annual_new_infections_p1, 0),
    ", respectively."
  ),
  "",
  "## Q4 Attribution",
  "",
  sprintf(
    "- Attribution contexts: %d profile-SSP-GCM-window combinations (%d profiles x %d SSPs x %d GCMs x %d windows), with %d network replicates per context.",
    q4_n_contexts,
    q4_n_profiles,
    q4_n_scenarios,
    q4_n_models,
    q4_n_windows,
    unique(q4_incidence_all$n_network_replicates)[1]
  ),
  sprintf(
    "- Annual infections: A = %s; D ranges from %s across future climate and movement contexts.",
    fmt_range(q4_incidence_all$A_current_climate_current_movement, 0),
    fmt_range(q4_incidence_all$D_future_climate_future_movement, 0)
  ),
  sprintf(
    "- Decomposition across contexts: climate = %s, movement = %s, interaction = %s; climate is dominant in %d/%d contexts.",
    fmt_range(q4_incidence_all$climate_contribution_mid, 0),
    fmt_range(q4_incidence_all$movement_contribution_mid, 0),
    fmt_range(q4_incidence_all$interaction_contribution_mid, 0),
    q4_n_climate_dominant,
    q4_n_contexts
  ),
  paste0(
    "- ",
    q4_window_summary$site_profile, ", ",
    q4_window_summary$scenario_label, " ",
    q4_window_summary$future_window,
    " window median: D = ",
    fmt_num(q4_window_summary$both_D_mid, 0),
    ", climate = ",
    fmt_num(q4_window_summary$climate_contribution_mid, 0),
    ", movement = ",
    fmt_num(q4_window_summary$movement_contribution_mid, 0),
    ", interaction = ",
    fmt_num(q4_window_summary$interaction_contribution_mid, 0),
    "."
  )
)

writeLines(md, file.path(out_dir, "module03_headline_summary.md"))

cat("\nModule 03 headline summary written to:\n")
cat(" - ", file.path(out_dir, "module03_headline_summary.md"), "\n", sep = "")
cat(" - ", file.path(out_dir, "module03_question_status.csv"), "\n", sep = "")
