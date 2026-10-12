#!/usr/bin/env Rscript
#' ####################################################################### #
#' 06 - Summarise and Plot Future NDVI Projections
#'
#' Summarise the two SSP ensembles over the common 20-year analysis windows.
#' Save only the tables needed for downstream interpretation and four figure
#' products: annual, seasonal, and covariate-support cross-SSP comparisons,
#' plus one spatial-change map per SSP.
#'
#' INPUTS:  derived/gimms_support/all_historical_observed.feather
#'          derived/gimms_support/future/{scenario}/all_future_{GCM}.feather
#'          outputs/scenarios/{scenario}/projections/{GCM}/monthly_ndvi.feather
#' OUTPUTS: outputs/scenarios/{scenario}/summaries/
#'          outputs/scenarios/{scenario}/figures/
#'          outputs/comparison/summaries/
#'          outputs/comparison/figures/
#' ####################################################################### #

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(purrr)
  library(readr)
  library(slider)
  library(tidyr)
})

source(here::here("analysis", "01_ndvi_prediction", "_config.R"))

theme_projection <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      legend.text = element_text(color = "black"),
      strip.text = element_text(face = "bold", color = "black"),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
}

sample_rows <- function(df, n) {
  if (nrow(df) <= n) df else slice_sample(df, n = n)
}

# --- Shared historical summaries -------------------------------------------

historical_path <- file.path(
  ndvi_config$gimms_support_dir,
  "all_historical_observed.feather"
)
if (!file.exists(historical_path)) {
  stop("Missing historical NDVI table: ", historical_path)
}

historical <- read_feather(historical_path) |>
  mutate(
    date = as.Date(date),
    year = year(date),
    month = month(date)
  )

historical_annual <- historical |>
  group_by(year) |>
  summarise(mean_ndvi = mean(ndvi, na.rm = TRUE), .groups = "drop")

# The movement and transmission baselines use the 2009-2018 observed cycle,
# so use the same period as the seasonal and spatial reference here.
baseline_years <- 2009:2018
historical_seasonal <- historical |>
  filter(year %in% baseline_years) |>
  group_by(month) |>
  summarise(mean_ndvi = mean(ndvi, na.rm = TRUE), .groups = "drop")

historical_spatial <- historical |>
  filter(year %in% baseline_years) |>
  group_by(x, y) |>
  summarise(baseline_ndvi = mean(ndvi, na.rm = TRUE), .groups = "drop")

historical_support <- historical |>
  arrange(x, y, date) |>
  group_by(x, y) |>
  mutate(
    pr_roll3 = slide_dbl(pr, mean, .before = 2, .complete = TRUE),
    tasmax_roll3 = slide_dbl(tasmax, mean, .before = 2, .complete = TRUE)
  ) |>
  ungroup() |>
  filter(!is.na(pr_roll3), !is.na(tasmax_roll3))

support_bounds <- historical_support |>
  summarise(
    pr_q01 = quantile(pr_roll3, 0.01, na.rm = TRUE),
    pr_q99 = quantile(pr_roll3, 0.99, na.rm = TRUE),
    tasmax_q01 = quantile(tasmax_roll3, 0.01, na.rm = TRUE),
    tasmax_q99 = quantile(tasmax_roll3, 0.99, na.rm = TRUE)
  )

window_lookup <- tibble::enframe(
  ndvi_config$analysis_windows,
  name = "window",
  value = "year"
) |>
  unnest_longer(year) |>
  mutate(
    year = as.integer(year),
    window = factor(window, levels = names(ndvi_config$analysis_windows))
  )

gcm_models <- ndvi_config$chelsa_reporting_models
scenario_labels <- ndvi_config$future_scenario_labels

scenario_colors <- c(
  ssp245 = "#0072B2",
  ssp370 = "#D55E00"
)
window_colors <- c(
  `Observed 2009-2018` = "black",
  `2041-2060` = "#0072B2",
  `2061-2080` = "#009E73",
  `2081-2100` = "#D55E00"
)
gcm_colors <- c(
  GFDL_ESM4 = "#1B9E77",
  IPSL_CM6A_LR = "#D95F02",
  MRI_ESM2_0 = "#7570B3",
  UKESM1_0_LL = "#E7298A"
)

# --- Scenario summaries -----------------------------------------------------

summarise_scenario <- function(scenario) {
  scenario_root <- ndvi_scenario_outputs_dir(scenario, create = TRUE)
  summary_dir <- file.path(scenario_root, "summaries")
  figure_dir <- file.path(scenario_root, "figures")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  future <- map_dfr(gcm_models, function(model_name) {
    path <- file.path(
      scenario_root,
      "projections",
      model_name,
      "monthly_ndvi.feather"
    )
    if (!file.exists(path)) stop("Missing NDVI projection: ", path)

    read_feather(path) |>
      mutate(
        date = as.Date(date),
        climate_model = model_name,
        scenario = scenario
      )
  })

  annual_by_gcm <- future |>
    group_by(scenario, climate_model, year) |>
    summarise(mean_ndvi = mean(ndvi_pred, na.rm = TRUE), .groups = "drop")

  annual_ensemble <- annual_by_gcm |>
    group_by(scenario, year) |>
    summarise(
      median_ndvi = median(mean_ndvi),
      min_ndvi = min(mean_ndvi),
      max_ndvi = max(mean_ndvi),
      .groups = "drop"
    )

  seasonal_by_gcm <- future |>
    left_join(window_lookup, by = "year") |>
    filter(!is.na(window)) |>
    group_by(scenario, climate_model, window, month) |>
    summarise(mean_ndvi = mean(ndvi_pred, na.rm = TRUE), .groups = "drop")

  seasonal_ensemble <- seasonal_by_gcm |>
    group_by(scenario, window, month) |>
    summarise(
      median_ndvi = median(mean_ndvi),
      min_ndvi = min(mean_ndvi),
      max_ndvi = max(mean_ndvi),
      .groups = "drop"
    )

  window_by_gcm <- annual_by_gcm |>
    left_join(window_lookup, by = "year") |>
    filter(!is.na(window)) |>
    group_by(scenario, climate_model, window) |>
    summarise(mean_ndvi = mean(mean_ndvi), .groups = "drop")

  future_support <- map_dfr(gcm_models, function(model_name) {
    feature_path <- file.path(
      ndvi_future_features_dir(scenario),
      paste0("all_future_", model_name, ".feather")
    )
    if (!file.exists(feature_path)) stop("Missing future features: ", feature_path)

    read_feather(feature_path) |>
      mutate(
        date = as.Date(date),
        year = year(date),
        climate_model = model_name
      ) |>
      arrange(x, y, date) |>
      group_by(x, y) |>
      mutate(
        pr_roll3 = slide_dbl(pr, mean, .before = 2, .complete = TRUE),
        tasmax_roll3 = slide_dbl(tasmax, mean, .before = 2, .complete = TRUE)
      ) |>
      ungroup() |>
      left_join(window_lookup, by = "year") |>
      filter(!is.na(window), !is.na(pr_roll3), !is.na(tasmax_roll3)) |>
      mutate(
        outside_training_range =
          pr_roll3 < support_bounds$pr_q01 |
          pr_roll3 > support_bounds$pr_q99 |
          tasmax_roll3 < support_bounds$tasmax_q01 |
          tasmax_roll3 > support_bounds$tasmax_q99
      )
  })

  support_summary <- future_support |>
    group_by(climate_model, window) |>
    summarise(
      cell_months = n(),
      pct_outside_training_range = 100 * mean(outside_training_range),
      pr_roll3_min = min(pr_roll3),
      pr_roll3_max = max(pr_roll3),
      tasmax_roll3_min = min(tasmax_roll3),
      tasmax_roll3_max = max(tasmax_roll3),
      .groups = "drop"
    )

  spatial_by_gcm <- future |>
    left_join(window_lookup, by = "year") |>
    filter(!is.na(window)) |>
    group_by(scenario, climate_model, window, x, y) |>
    summarise(future_ndvi = mean(ndvi_pred, na.rm = TRUE), .groups = "drop") |>
    left_join(historical_spatial, by = c("x", "y")) |>
    mutate(delta_ndvi = future_ndvi - baseline_ndvi)

  spatial_ensemble <- spatial_by_gcm |>
    group_by(scenario, window, x, y) |>
    summarise(
      median_delta_ndvi = median(delta_ndvi),
      min_delta_ndvi = min(delta_ndvi),
      max_delta_ndvi = max(delta_ndvi),
      .groups = "drop"
    )

  write_csv(annual_by_gcm, file.path(summary_dir, "annual_ndvi_by_gcm.csv"))
  write_csv(seasonal_by_gcm, file.path(summary_dir, "seasonal_ndvi_by_gcm.csv"))
  write_csv(window_by_gcm, file.path(summary_dir, "window_ndvi_by_gcm.csv"))
  write_csv(support_summary, file.path(summary_dir, "climate_support.csv"))
  write_csv(spatial_by_gcm, file.path(summary_dir, "spatial_ndvi_change.csv"))

  future_sample <- future_support |>
    group_by(climate_model, window) |>
    group_modify(~ sample_rows(.x, 2500)) |>
    ungroup()

  spatial_limit <- max(abs(spatial_ensemble$median_delta_ndvi), na.rm = TRUE)
  spatial_plot <- spatial_ensemble |>
    ggplot(aes(x, y, fill = median_delta_ndvi)) +
    geom_tile() +
    coord_equal() +
    facet_wrap(~ window, nrow = 1) +
    scale_fill_gradient2(
      low = "#8C510A",
      mid = "white",
      high = "#1B7837",
      midpoint = 0,
      limits = c(-spatial_limit, spatial_limit),
      name = "Change in NDVI"
    ) +
    theme_projection() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    ) +
    labs(x = NULL, y = NULL)

  ggsave(
    file.path(figure_dir, "spatial_ndvi_change.png"),
    spatial_plot,
    width = 10,
    height = 4.2,
    dpi = 300
  )

  list(
    annual_by_gcm = annual_by_gcm,
    annual_ensemble = annual_ensemble,
    seasonal_by_gcm = seasonal_by_gcm,
    seasonal_ensemble = seasonal_ensemble,
    window_by_gcm = window_by_gcm,
    support_sample = future_sample,
    support_summary = support_summary
  )
}

scenario_results <- setNames(
  map(ndvi_config$future_scenarios, summarise_scenario),
  ndvi_config$future_scenarios
)

# --- Cross-scenario figures and tables --------------------------------------

comparison_root <- file.path(ndvi_config$outputs_dir, "comparison")
comparison_summary_dir <- file.path(comparison_root, "summaries")
comparison_figure_dir <- file.path(comparison_root, "figures")
dir.create(comparison_summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(comparison_figure_dir, recursive = TRUE, showWarnings = FALSE)

annual_ensemble <- map_dfr(scenario_results, "annual_ensemble") |>
  mutate(scenario_label = unname(scenario_labels[scenario]))
seasonal_ensemble <- map_dfr(scenario_results, "seasonal_ensemble") |>
  mutate(scenario_label = unname(scenario_labels[scenario]))
window_by_gcm <- map_dfr(scenario_results, "window_by_gcm")
future_support_sample <- imap_dfr(
  scenario_results,
  ~ mutate(.x$support_sample, scenario = .y)
)
future_support_summary <- imap_dfr(
  scenario_results,
  ~ mutate(.x$support_summary, scenario = .y)
)

window_ensemble <- window_by_gcm |>
  group_by(scenario, window) |>
  summarise(
    median_ndvi = median(mean_ndvi),
    min_ndvi = min(mean_ndvi),
    max_ndvi = max(mean_ndvi),
    .groups = "drop"
  ) |>
  mutate(scenario_label = unname(scenario_labels[scenario]))

annual_plot <- ggplot() +
  geom_line(
    data = historical_annual,
    aes(year, mean_ndvi, color = "Observed 1982-2018"),
    linewidth = 0.75
  ) +
  geom_ribbon(
    data = annual_ensemble,
    aes(year, ymin = min_ndvi, ymax = max_ndvi, fill = scenario),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(
    data = annual_ensemble,
    aes(year, median_ndvi, color = scenario),
    linewidth = 0.9
  ) +
  geom_vline(
    xintercept = ndvi_config$projection_start_year,
    linetype = 2,
    color = "grey35"
  ) +
  scale_color_manual(
    values = c("Observed 1982-2018" = "black", scenario_colors),
    breaks = c("Observed 1982-2018", ndvi_config$future_scenarios),
    labels = c(
      "Observed 1982-2018",
      unname(scenario_labels[ndvi_config$future_scenarios])
    ),
    name = NULL
  ) +
  scale_fill_manual(values = scenario_colors, guide = "none") +
  theme_projection() +
  theme(legend.position = "top") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x = NULL, y = "Annual mean NDVI")

historical_seasonal_plot <- tidyr::crossing(
  historical_seasonal,
  scenario = ndvi_config$future_scenarios
) |>
  mutate(period = "Observed 2009-2018") |>
  mutate(scenario_label = unname(scenario_labels[scenario]))

seasonal_plot <- ggplot() +
  geom_ribbon(
    data = seasonal_ensemble,
    aes(month, ymin = min_ndvi, ymax = max_ndvi, fill = window),
    alpha = 0.14,
    color = NA
  ) +
  geom_line(
    data = seasonal_ensemble,
    aes(month, median_ndvi, color = window),
    linewidth = 0.9
  ) +
  geom_line(
    data = historical_seasonal_plot,
    aes(month, mean_ndvi, color = period),
    linewidth = 0.9
  ) +
  facet_wrap(~ scenario_label, ncol = 2) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_manual(
    values = window_colors,
    breaks = names(window_colors),
    name = NULL
  ) +
  scale_fill_manual(
    values = window_colors[names(ndvi_config$analysis_windows)],
    guide = "none"
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_projection() +
  labs(x = NULL, y = "Mean NDVI")

# Display the two pathways in the same climate-support panels. GCMs and future
# windows already define the facets, so colour is reserved for SSP. This makes
# extrapolation under SSP2-4.5 and SSP3-7.0 directly comparable without asking
# readers to compare two separate figures.
support_box <- tibble(
  xmin = support_bounds$pr_q01,
  xmax = support_bounds$pr_q99,
  ymin = support_bounds$tasmax_q01,
  ymax = support_bounds$tasmax_q99
)
support_labels <- future_support_summary |>
  mutate(
    scenario_label = unname(scenario_labels[.data$scenario]),
    label_line = sprintf(
      "%s: %.0f%% outside",
      .data$scenario_label,
      .data$pct_outside_training_range
    )
  ) |>
  arrange(.data$scenario) |>
  group_by(.data$climate_model, .data$window) |>
  summarise(label = paste(.data$label_line, collapse = "\n"), .groups = "drop")

set.seed(101)
historical_support_sample <- sample_rows(historical_support, 20000)
climate_support_plot <- ggplot() +
  geom_point(
    data = historical_support_sample,
    aes(.data$pr_roll3, .data$tasmax_roll3),
    colour = "grey70",
    alpha = 0.08,
    size = 0.25
  ) +
  geom_rect(
    data = support_box,
    aes(xmin = .data$xmin, xmax = .data$xmax,
        ymin = .data$ymin, ymax = .data$ymax),
    fill = NA,
    colour = "black",
    linetype = 2,
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = future_support_sample,
    aes(
      .data$pr_roll3,
      .data$tasmax_roll3,
      colour = .data$scenario
    ),
    alpha = 0.16,
    size = 0.3
  ) +
  geom_text(
    data = support_labels,
    aes(x = -Inf, y = Inf, label = .data$label),
    hjust = -0.08,
    vjust = 1.15,
    size = 3.2,
    lineheight = 0.95,
    inherit.aes = FALSE
  ) +
  facet_grid(climate_model ~ window) +
  scale_colour_manual(
    values = scenario_colors,
    labels = unname(scenario_labels[names(scenario_colors)]),
    name = NULL
  ) +
  theme_projection(base_size = 12) +
  theme(legend.position = "top") +
  # Points are drawn faint to cope with overplotting, which leaves the legend
  # keys almost invisible unless they are overridden.
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.6))) +
  labs(
    x = "Three-month mean rainfall (mm)",
    y = "Three-month mean maximum temperature (deg C)"
  )

write_csv(
  annual_ensemble,
  file.path(comparison_summary_dir, "annual_ndvi_ensemble.csv")
)
write_csv(
  seasonal_ensemble,
  file.path(comparison_summary_dir, "seasonal_ndvi_ensemble.csv")
)
write_csv(
  window_ensemble,
  file.path(comparison_summary_dir, "window_ndvi_ensemble.csv")
)

ggsave(
  file.path(comparison_figure_dir, "annual_ndvi_by_ssp.png"),
  annual_plot,
  width = 9,
  height = 4.8,
  dpi = 300
)
ggsave(
  file.path(comparison_figure_dir, "seasonal_ndvi_by_ssp.png"),
  seasonal_plot,
  width = 11,
  height = 5.2,
  dpi = 300
)
ggsave(
  file.path(comparison_figure_dir, "climate_support_by_ssp.png"),
  climate_support_plot,
  width = 11,
  height = 10,
  dpi = 300
)

cat("Saved two-SSP NDVI summaries to:", comparison_root, "\n")
