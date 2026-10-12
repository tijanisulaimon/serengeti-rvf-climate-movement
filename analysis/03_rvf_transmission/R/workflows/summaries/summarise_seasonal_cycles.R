#!/usr/bin/env Rscript

#' 09_summarise_spatial_seasonal_patterns.R
#' -------------------------------------------------------------------------
#' Compare settled within-year infectious-cattle cycles under the observed
#' climate cycle and the three future windows. The script reads complete-cycle
#' diagnostics from scripts 03 and 06; it does not rerun transmission.
#'
#' Absolute trajectories describe seasonal burden. Monthly shares of each
#' run's annual infectious-host total separate seasonal timing from the much
#' larger change in overall burden. Both SSPs and the two contrasting exposure
#' profiles are shown.
#' -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tidyr)
})

source(here::here("analysis", "03_rvf_transmission", "R", "transmission_model.R"))
source(here::here("analysis", "03_rvf_transmission", "R", "plot_helpers.R"))

cfg <- load_rvf_config()
sim <- cfg$sim

if (!identical(sim$network_mode, "monthly")) {
  stop("The publication seasonal summary requires monthly movement networks.")
}

primary_profile <- rvf_resolve_site_profile_id(
  sim$future_site_profile %||% "grazing_floodplain_dominant"
)
profiles <- unique(c(primary_profile, "peridomestic_dominant"))
profile_plot_labels <- c(
  grazing_floodplain_dominant = "Shared grazing/floodplain exposure",
  peridomestic_dominant = "Village-associated exposure"
)

out_dir <- file.path(cfg$results_dir, "seasonal_cycles")
figures_dir <- cfg$figures_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Function: read settled monthly diagnostics for one calibrated profile.
# Profile-specific baseline tags ensure that current and future trajectories
# use the same exposure allocation and capacity calibration.
read_profile_seasonal_outputs <- function(profile) {
  profile_sim <- sim
  profile_sim$site_profile <- profile
  profile_sim <- rvf_apply_profile_calibration(
    profile_sim,
    cfg = cfg,
    require_target_band = TRUE
  )

  baseline_path <- file.path(
    rvf_baseline_intermediate_dir(profile_sim, cfg),
    paste0(rvf_baseline_run_tag(profile_sim), "_cycle_diagnostics.rds")
  )

  future_name <- if (identical(profile, primary_profile)) {
    "future_cold_endemic"
  } else {
    paste0("future_cold_endemic_profile_", profile)
  }
  future_name <- rvf_rainfall_output_subdir(
    future_name,
    profile_sim$rainfall_model,
    cfg
  )
  future_path <- file.path(
    cfg$intermediate_dir,
    future_name,
    "future_cold_endemic_cycle_diagnostics.rds"
  )

  missing_paths <- c(baseline_path, future_path)[
    !file.exists(c(baseline_path, future_path))
  ]
  if (length(missing_paths) > 0L) {
    stop(
      "Missing seasonal input(s) for ", profile, ":\n  ",
      paste(missing_paths, collapse = "\n  "),
      "\nRun scripts 03 and 06 for this profile before script 09."
    )
  }

  baseline_diagnostics <- readRDS(baseline_path)
  future_diagnostics <- readRDS(future_path)
  if (is.null(baseline_diagnostics$seasonal_by_run) ||
      is.null(future_diagnostics$seasonal_by_run)) {
    stop("Cycle diagnostics lack seasonal_by_run for profile: ", profile)
  }

  current <- baseline_diagnostics$seasonal_by_run |>
    transmute(
      site_profile = profile,
      profile_label = unname(profile_plot_labels[[profile]]),
      scenario = "current",
      scenario_label = "Current",
      future_window = "Current",
      future_year = NA_integer_,
      climate_model = "Observed",
      .data$network_replicate,
      .data$month,
      .data$mean_infectious_hosts
    )

  future <- future_diagnostics$seasonal_by_run |>
    filter(
      .data$scenario %in% cfg$future_scenarios,
      .data$climate_model %in% cfg$climate_models,
      abs(
        .data$external_introduction_per_year -
          profile_sim$external_introduction_per_year
      ) < 1e-12
    ) |>
    transmute(
      site_profile = profile,
      profile_label = unname(profile_plot_labels[[profile]]),
      .data$scenario,
      scenario_label = unname(cfg$scenario_labels[.data$scenario]),
      .data$future_window,
      future_year = as.integer(.data$future_window),
      .data$climate_model,
      .data$network_replicate,
      .data$month,
      .data$mean_infectious_hosts
    )

  bind_rows(current, future)
}

seasonal_by_run <- bind_rows(lapply(profiles, read_profile_seasonal_outputs)) |>
  mutate(month_index = match(.data$month, month.abb))

# Every run must contribute one estimate for each calendar month.
seasonal_counts <- seasonal_by_run |>
  count(
    .data$site_profile,
    .data$scenario,
    .data$future_window,
    .data$climate_model,
    .data$network_replicate,
    name = "n_months"
  )
if (any(seasonal_counts$n_months != 12L) ||
    anyNA(seasonal_by_run$month_index)) {
  stop("At least one seasonal run does not contain all 12 calendar months.")
}
if (any(!is.finite(seasonal_by_run$mean_infectious_hosts)) ||
    any(seasonal_by_run$mean_infectious_hosts < 0)) {
  stop("Seasonal infectious-host summaries must be finite and non-negative.")
}

# Require the complete 20-replicate publication grid. The 12-month check above
# verifies each run internally but cannot detect an older, smaller run set.
observed_seasonal_runs <- seasonal_by_run |>
  filter(.data$scenario == "current") |>
  distinct(.data$site_profile, .data$network_replicate)
expected_observed_seasonal_runs <- crossing(
  site_profile = profiles,
  network_replicate = sort(unique(cfg$analysis_network_replicates))
)
future_seasonal_runs <- seasonal_by_run |>
  filter(.data$scenario %in% cfg$future_scenarios) |>
  distinct(
    .data$site_profile,
    .data$scenario,
    .data$future_window,
    .data$climate_model,
    .data$network_replicate
  )
expected_future_seasonal_runs <- crossing(
  site_profile = profiles,
  scenario = cfg$future_scenarios,
  future_window = names(cfg$future_windows),
  climate_model = cfg$climate_models,
  network_replicate = sort(unique(cfg$analysis_network_replicates))
)
if (!setequal(
  observed_seasonal_runs,
  expected_observed_seasonal_runs
) || !setequal(
  future_seasonal_runs,
  expected_future_seasonal_runs
)) {
  stop(
    "Seasonal summaries contain an incomplete, duplicated, or stale ",
    "profile-by-scenario-by-replicate grid."
  )
}

# Monthly shares sum to one within each run. This standardisation preserves
# timing and seasonal breadth while removing differences in annual magnitude.
seasonal_by_run <- seasonal_by_run |>
  group_by(
    .data$site_profile,
    .data$scenario,
    .data$future_window,
    .data$climate_model,
    .data$network_replicate
  ) |>
  mutate(
    annual_infectious_host_sum = sum(.data$mean_infectious_hosts),
    monthly_burden_share =
      .data$mean_infectious_hosts / .data$annual_infectious_host_sum
  ) |>
  ungroup()
if (any(seasonal_by_run$annual_infectious_host_sum <= 0)) {
  stop("At least one run has no infectious-host burden to standardise.")
}

# Seasonal breadth is the inverse Simpson concentration of monthly burden.
# It ranges from 1, when all annual burden occurs in one month, to 12, when
# burden is distributed evenly through the year. The metric complements peak
# month by detecting a broader transmission season even when the peak is fixed.
seasonal_breadth_by_run <- seasonal_by_run |>
  group_by(
    .data$site_profile,
    .data$profile_label,
    .data$scenario,
    .data$scenario_label,
    .data$future_window,
    .data$future_year,
    .data$climate_model,
    .data$network_replicate
  ) |>
  summarise(
    effective_months = 1 / sum(.data$monthly_burden_share^2),
    peak_month = .data$month[which.max(.data$mean_infectious_hosts)],
    peak_month_share = max(.data$monthly_burden_share),
    .groups = "drop"
  )

# Translate each settled monthly trajectory into three directly interpretable
# characteristics. Values are expressed per 1,000 cattle across the 94-village
# district; peak month is retained on the calendar-month scale.
district_cattle <- 94 * as.numeric(sim$host_population_per_village)
seasonal_characteristics_by_run <- seasonal_by_run |>
  group_by(
    .data$site_profile,
    .data$profile_label,
    .data$scenario,
    .data$scenario_label,
    .data$future_window,
    .data$future_year,
    .data$climate_model,
    .data$network_replicate
  ) |>
  summarise(
    peak_infectious_per_1000 =
      1000 * max(.data$mean_infectious_hosts) / .env$district_cattle,
    trough_infectious_per_1000 =
      1000 * min(.data$mean_infectious_hosts) / .env$district_cattle,
    amplitude_per_1000 =
      1000 * (max(.data$mean_infectious_hosts) -
        min(.data$mean_infectious_hosts)) / .env$district_cattle,
    peak_month_index = .data$month_index[which.max(.data$mean_infectious_hosts)],
    peak_month = .data$month[which.max(.data$mean_infectious_hosts)],
    .groups = "drop"
  )

seasonal_characteristics_summary <- seasonal_characteristics_by_run |>
  group_by(
    .data$site_profile,
    .data$profile_label,
    .data$scenario,
    .data$scenario_label,
    .data$future_window,
    .data$future_year
  ) |>
  summarise(
    n_runs = n(),
    across(
      c(
        "peak_infectious_per_1000",
        "trough_infectious_per_1000",
        "amplitude_per_1000",
        "peak_month_index"
      ),
      list(
        mid = median,
        lo = ~ as.numeric(quantile(.x, 0.10)),
        hi = ~ as.numeric(quantile(.x, 0.90))
      ),
      .names = "{.col}_{.fn}"
    ),
    july_peak_fraction = mean(.data$peak_month == "Jul"),
    .groups = "drop"
  ) |>
  arrange(.data$site_profile, .data$scenario, .data$future_year)

if (any(!is.finite(seasonal_breadth_by_run$effective_months)) ||
    any(seasonal_breadth_by_run$effective_months < 1 - 1e-8) ||
    any(seasonal_breadth_by_run$effective_months > 12 + 1e-8)) {
  stop("Seasonal effective-month estimates must lie between 1 and 12.")
}

# Current intervals span network replicates; future intervals combine the four
# GCMs and network replicates within each SSP and future window.
seasonal_breadth_summary <- seasonal_breadth_by_run |>
  group_by(
    .data$site_profile,
    .data$profile_label,
    .data$scenario,
    .data$scenario_label,
    .data$future_window,
    .data$future_year
  ) |>
  summarise(
    n_runs = n(),
    effective_months_mid = median(.data$effective_months),
    effective_months_lo = as.numeric(
      quantile(.data$effective_months, 0.10)
    ),
    effective_months_hi = as.numeric(
      quantile(.data$effective_months, 0.90)
    ),
    peak_month_share_mid = median(.data$peak_month_share),
    peak_month_share_lo = as.numeric(
      quantile(.data$peak_month_share, 0.10)
    ),
    peak_month_share_hi = as.numeric(
      quantile(.data$peak_month_share, 0.90)
    ),
    .groups = "drop"
  ) |>
  arrange(.data$site_profile, .data$scenario, .data$future_year)

# Current intervals span network replicates; future intervals combine the four
# GCMs and network replicates within each SSP and future window.
seasonal_summary <- seasonal_by_run |>
  group_by(
    .data$site_profile,
    .data$profile_label,
    .data$scenario,
    .data$scenario_label,
    .data$future_window,
    .data$future_year,
    .data$month,
    .data$month_index
  ) |>
  summarise(
    n_runs = n(),
    infectious_hosts_mid = median(.data$mean_infectious_hosts),
    infectious_hosts_lo = as.numeric(
      quantile(.data$mean_infectious_hosts, 0.10)
    ),
    infectious_hosts_hi = as.numeric(
      quantile(.data$mean_infectious_hosts, 0.90)
    ),
    monthly_share_mid = median(.data$monthly_burden_share),
    monthly_share_lo = as.numeric(
      quantile(.data$monthly_burden_share, 0.10)
    ),
    monthly_share_hi = as.numeric(
      quantile(.data$monthly_burden_share, 0.90)
    ),
    .groups = "drop"
  ) |>
  arrange(.data$site_profile, .data$scenario, .data$future_year,
          .data$month_index)

# Retain the full peak-month distribution rather than reporting only a mode;
# close June-July peaks can otherwise look more certain than they are.
peak_month_distribution <- seasonal_by_run |>
  group_by(
    .data$site_profile,
    .data$profile_label,
    .data$scenario,
    .data$scenario_label,
    .data$future_window,
    .data$future_year,
    .data$climate_model,
    .data$network_replicate
  ) |>
  slice_max(
    order_by = .data$mean_infectious_hosts,
    n = 1L,
    with_ties = FALSE
  ) |>
  ungroup() |>
  count(
    .data$site_profile,
    .data$profile_label,
    .data$scenario,
    .data$scenario_label,
    .data$future_window,
    .data$future_year,
    peak_month = .data$month,
    name = "n_runs"
  ) |>
  group_by(
    .data$site_profile,
    .data$scenario,
    .data$future_window
  ) |>
  mutate(peak_month_fraction = .data$n_runs / sum(.data$n_runs)) |>
  ungroup() |>
  arrange(.data$site_profile, .data$scenario, .data$future_year,
          match(.data$peak_month, month.abb))

write_csv(
  seasonal_summary,
  file.path(out_dir, "seasonal_cycle_summary.csv")
)
write_csv(
  peak_month_distribution,
  file.path(out_dir, "seasonal_peak_month_distribution.csv")
)
write_csv(
  seasonal_breadth_by_run,
  file.path(out_dir, "seasonal_breadth_by_run.csv")
)
write_csv(
  seasonal_breadth_summary,
  file.path(out_dir, "seasonal_breadth_summary.csv")
)
write_csv(
  seasonal_characteristics_by_run,
  file.path(out_dir, "seasonal_characteristics_by_run.csv")
)
write_csv(
  seasonal_characteristics_summary,
  file.path(out_dir, "seasonal_characteristics_summary.csv")
)

# Repeat the observed baseline in each SSP column so both pathways share an
# identical reference. Future lines pool GCM and network-replicate uncertainty.
current_summary <- seasonal_summary |>
  filter(.data$scenario == "current") |>
  select(-any_of(c("scenario", "scenario_label")))
current_plot <- crossing(
  current_summary,
  scenario = cfg$future_scenarios
) |>
  mutate(
    scenario_label = unname(cfg$scenario_labels[.data$scenario]),
    series = "Current"
  )

future_plot <- seasonal_summary |>
  filter(.data$scenario != "current") |>
  mutate(
    series = vapply(
      .data$future_window,
      function(window) {
        years <- cfg$future_windows[[window]]
        sprintf("%d-%d", min(years), max(years))
      },
      character(1)
    )
  )

series_levels <- c(
  "Current",
  vapply(
    cfg$future_windows,
    function(years) sprintf("%d-%d", min(years), max(years)),
    character(1)
  )
)
series_colours <- c(
  Current = "#1A1A1A",
  `2041-2060` = "#0072B2",
  `2061-2080` = "#D55E00",
  `2081-2100` = "#7A5195"
)

plot_data <- bind_rows(current_plot, future_plot) |>
  mutate(
    profile_label = factor(
      .data$profile_label,
      levels = unname(profile_plot_labels[profiles])
    ),
    scenario_label = factor(
      .data$scenario_label,
      levels = unname(cfg$scenario_labels[cfg$future_scenarios])
    ),
    series = factor(.data$series, levels = series_levels)
  )

# Function: draw a settled seasonal trajectory for one burden scale.
# y_mid/y_lo/y_hi identify either absolute infectious hosts or monthly shares.
plot_seasonal_cycle <- function(data, y_mid, y_lo, y_hi, y_label) {
  ggplot(
    data,
    aes(
      x = .data$month_index,
      y = .data[[y_mid]],
      colour = .data$series,
      fill = .data$series,
      group = .data$series
    )
  ) +
    geom_ribbon(
      aes(ymin = .data[[y_lo]], ymax = .data[[y_hi]]),
      alpha = 0.12,
      colour = NA
    ) +
    geom_line(linewidth = 0.9) +
    facet_grid(profile_label ~ scenario_label) +
    scale_x_continuous(breaks = seq_along(month.abb), labels = month.abb) +
    scale_colour_manual(values = series_colours, name = NULL, drop = FALSE) +
    scale_fill_manual(values = series_colours, guide = "none", drop = FALSE) +
    labs(x = NULL, y = y_label) +
    theme_classic(base_size = 12.5) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = 10),
      axis.title = element_text(colour = "black", size = 12.5),
      strip.background = element_blank(),
      strip.text = element_text(colour = "black", size = 11, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(colour = "black", size = 10.5)
    )
}

p_absolute <- plot_seasonal_cycle(
  plot_data,
  y_mid = "infectious_hosts_mid",
  y_lo = "infectious_hosts_lo",
  y_hi = "infectious_hosts_hi",
  y_label = "Mean infectious cattle"
)
p_timing <- plot_seasonal_cycle(
  plot_data,
  y_mid = "monthly_share_mid",
  y_lo = "monthly_share_lo",
  y_hi = "monthly_share_hi",
  y_label = "Monthly share of annual infectious-host burden"
) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1))

rvf_ggsave(
  file.path(figures_dir, "07_seasonal_infectious_cycle.png"),
  p_absolute,
  width = 11.2,
  height = 7.2
)
rvf_ggsave(
  file.path(figures_dir, "07_seasonal_timing.png"),
  p_timing,
  width = 11.2,
  height = 7.2
)

# A characteristic-based summary makes the substantive result easier to read
# than another set of overlaid monthly curves: future peak-to-trough amplitude
# and trough burden rise, whereas the peak remains concentrated in June-July.
time_levels <- c(
  "Current",
  vapply(
    cfg$future_windows,
    function(years) sprintf("%d-%d", min(years), max(years)),
    character(1)
  )
)
future_window_labels <- setNames(time_levels[-1], names(cfg$future_windows))
characteristic_plot_data <- seasonal_characteristics_summary |>
  mutate(
    time_label = if_else(
      .data$scenario == "current",
      "Current",
      unname(.env$future_window_labels[.data$future_window])
    ),
    time_label = factor(.data$time_label, levels = time_levels),
    scenario_label = factor(
      .data$scenario_label,
      levels = c("Current", unname(cfg$scenario_labels[cfg$future_scenarios]))
    ),
    profile_label = factor(
      .data$profile_label,
      levels = unname(profile_plot_labels[profiles])
    )
  )

characteristic_colours <- setNames(
  c("#1A1A1A", "#0072B2", "#D55E00"),
  c("Current", "SSP2-4.5", "SSP3-7.0")
)
characteristic_shapes <- setNames(
  c(18, 16, 17),
  c("Current", "SSP2-4.5", "SSP3-7.0")
)
characteristic_dodge <- position_dodge(width = 0.46)

# Function: draw one median and 10th-90th percentile characteristic panel.
# Called for amplitude and trough abundance, which share units and structure.
plot_seasonal_characteristic <- function(data, mid, lo, hi, y_label) {
  ggplot(
    data,
    aes(
      x = .data$time_label,
      y = .data[[mid]],
      colour = .data$scenario_label,
      shape = .data$scenario_label,
      group = .data$scenario_label
    )
  ) +
    geom_linerange(
      aes(ymin = .data[[lo]], ymax = .data[[hi]]),
      position = characteristic_dodge,
      linewidth = 0.65
    ) +
    geom_point(position = characteristic_dodge, size = 2.6) +
    facet_wrap(vars(.data$profile_label), ncol = 1) +
    scale_colour_manual(values = characteristic_colours, name = NULL) +
    scale_shape_manual(values = characteristic_shapes, name = NULL) +
    labs(x = NULL, y = y_label) +
    theme_classic(base_size = 11.5) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = 9.2),
      axis.text.x = element_text(angle = 18, hjust = 1),
      axis.title = element_text(colour = "black", size = 10.5),
      strip.background = element_blank(),
      strip.text = element_text(colour = "black", face = "bold", size = 10),
      legend.position = "top",
      legend.text = element_text(colour = "black", size = 9.5)
    )
}

p_amplitude <- plot_seasonal_characteristic(
  characteristic_plot_data,
  "amplitude_per_1000_mid",
  "amplitude_per_1000_lo",
  "amplitude_per_1000_hi",
  "Peak-to-trough amplitude\nper 1,000 cattle"
)
p_trough <- plot_seasonal_characteristic(
  characteristic_plot_data,
  "trough_infectious_per_1000_mid",
  "trough_infectious_per_1000_lo",
  "trough_infectious_per_1000_hi",
  "Trough infectious cattle\nper 1,000"
)
# Peak month is a discrete outcome. A frequency heatmap retains the June-July
# split directly; quantile error bars on month numbers would imply a continuous
# scale and make near-identical distributions look erratic.
peak_column_levels <- c(
  "Observed",
  unlist(lapply(
    unname(cfg$scenario_labels[cfg$future_scenarios]),
    function(ssp_label) paste(ssp_label, time_levels[-1], sep = "\n")
  ))
)
peak_month_plot_data <- peak_month_distribution |>
  mutate(
    comparison_label = if_else(
      .data$scenario == "current",
      "Observed",
      paste(
        .data$scenario_label,
        unname(.env$future_window_labels[.data$future_window]),
        sep = "\n"
      )
    )
  ) |>
  select(
    .data$profile_label,
    .data$comparison_label,
    .data$peak_month,
    .data$peak_month_fraction
  )

# Include every peak month present in the simulations. Restricting the plot to
# June and July would silently hide an unexpected shift after a future rerun.
peak_month_levels <- month.abb[
  month.abb %in% unique(peak_month_plot_data$peak_month)
]
peak_month_plot_data <- crossing(
  profile_label = unname(profile_plot_labels[profiles]),
  comparison_label = peak_column_levels,
  peak_month = peak_month_levels
) |>
  left_join(
    peak_month_plot_data,
    by = c("profile_label", "comparison_label", "peak_month")
  ) |>
  mutate(
    peak_month_fraction = replace_na(.data$peak_month_fraction, 0),
    profile_label = factor(
      .data$profile_label,
      levels = unname(profile_plot_labels[profiles])
    ),
    comparison_label = factor(
      .data$comparison_label,
      levels = peak_column_levels
    ),
    peak_month = factor(
      .data$peak_month,
      levels = rev(peak_month_levels)
    )
  )

p_peak_month <- ggplot(
  peak_month_plot_data,
  aes(
    x = .data$comparison_label,
    y = .data$peak_month,
    fill = .data$peak_month_fraction
  )
) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(
    aes(
      label = scales::percent(.data$peak_month_fraction, accuracy = 1),
      colour = .data$peak_month_fraction >= 0.5
    ),
    size = 3
  ) +
  facet_wrap(vars(.data$profile_label), ncol = 1) +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1),
    labels = scales::label_percent(accuracy = 1),
    name = "Runs"
  ) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = "none") +
  labs(x = NULL, y = "Peak month") +
  theme_classic(base_size = 11.5) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black", size = 9.2),
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.title = element_text(colour = "black", size = 10.5),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", face = "bold", size = 10),
    legend.position = "top",
    legend.text = element_text(colour = "black", size = 9.5),
    legend.title = element_text(colour = "black", size = 9.5)
  )

# The continuous characteristics share the top row; the categorical peak-month
# summary spans the full lower row so its seven scenario columns remain legible.
seasonal_characteristics_plot <- (p_amplitude + p_trough) / p_peak_month +
  plot_layout(heights = c(1, 0.82), guides = "collect") +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      legend.position = "top",
      plot.tag = element_text(face = "bold", size = 13)
    )
  ) &
  theme(legend.position = "top")

rvf_ggsave(
  file.path(figures_dir, "07_seasonal_characteristics.png"),
  seasonal_characteristics_plot,
  width = 11.4,
  height = 7.4
)

message("Seasonal-cycle summaries written to: ", out_dir)
