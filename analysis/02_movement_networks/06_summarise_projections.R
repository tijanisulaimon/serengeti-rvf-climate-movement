#!/usr/bin/env Rscript
# Summarise projected movement regimes and one-mode network change.
#
# Publication figures focus on the two quantities needed downstream:
# the frequency of resource-condition regimes and changes in monthly contact
# density and geographic reach. The large bipartite bundles are not read.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

source(here::here("analysis", "02_movement_networks", "_config.R"))

cfg <- mvt_config
fit_file <- file.path(cfg$calibration_dir, "fit.rds")
if (!file.exists(fit_file)) {
  stop("Run 03_fit_model.R before summarising projections.")
}
current_fit_digest <- unname(tools::md5sum(fit_file))

index_file <- file.path(
  cfg$projection_summaries_dir, "network_bundle_index.csv"
)
if (!file.exists(index_file)) {
  stop("Run 05_project_networks.R before summarising projections.")
}

bundle_index <- read_csv(index_file, show_col_types = FALSE)

# A complete index lists the ten observed years plus every future
# SSP-GCM-year combination. Rejecting an incomplete index prevents a
# partially interrupted projection run from yielding publication summaries.
expected_jobs <- length(cfg$observed_years) +
  length(cfg$scenarios) * length(cfg$climate_models) * length(cfg$future_years)
if (nrow(bundle_index) != expected_jobs || anyDuplicated(bundle_index$job_id)) {
  stop(
    "The projection index is incomplete or duplicated (",
    nrow(bundle_index), " rows; expected ", expected_jobs,
    "). Re-run 05_project_networks.R."
  )
}
required_columns <- c("fit_digest", "n_replicates")
if (!all(required_columns %in% names(bundle_index))) {
  stop(
    "The projection index predates the current fit checks. ",
    "Re-run 05_project_networks.R."
  )
}
if (
  any(bundle_index$fit_digest != current_fit_digest) ||
    any(as.integer(bundle_index$n_replicates) != cfg$projection_replicates)
) {
  stop(
    "The projection index was generated from an earlier movement fit or ",
    "replicate setting. Re-run 05_project_networks.R."
  )
}

missing_summaries <- bundle_index$summary_file[!file.exists(bundle_index$summary_file)]
if (length(missing_summaries) > 0L) {
  stop(length(missing_summaries), " compact projection summaries are missing. ",
       "Re-run 05_project_networks.R; completed jobs will be skipped.")
}

# Read one compact job summary and attach scenario identifiers. Keeping this
# cache separate from the large bundles makes interrupted projection runs cheap
# to resume and avoids loading Module 03 inputs merely to draw figures.
read_job_summary <- function(index) {
  object <- readRDS(bundle_index$summary_file[[index]])
  metadata <- object$metadata
  if (!identical(metadata$fit_digest, current_fit_digest)) {
    stop(
      "Projection summary was generated from an earlier movement fit: ",
      bundle_index$summary_file[[index]]
    )
  }
  if (!identical(
    as.integer(metadata$n_replicates),
    as.integer(cfg$projection_replicates)
  )) {
    stop(
      "Projection summary has ", metadata$n_replicates,
      " replicates; expected ", cfg$projection_replicates, ": ",
      bundle_index$summary_file[[index]]
    )
  }

  identifying_columns <- list(
    dataset = metadata$dataset,
    scenario = metadata$scenario,
    scenario_label = metadata$scenario_label,
    climate_model = metadata$climate_model,
    future_window = metadata$future_window,
    year = metadata$year
  )

  list(
    monthly = object$monthly |>
      mutate(!!!identifying_columns, .before = 1),
    annual = object$annual |>
      mutate(!!!identifying_columns, .before = 1),
    regimes = tibble(
      dataset = metadata$dataset,
      scenario = metadata$scenario,
      scenario_label = metadata$scenario_label,
      climate_model = metadata$climate_model,
      future_window = metadata$future_window,
      year = metadata$year,
      month = names(metadata$regimes),
      regime = unname(metadata$regimes)
    )
  )
}

job_summaries <- if (.Platform$OS.type == "unix" && cfg$n_cores > 1L) {
  parallel::mclapply(
    seq_len(nrow(bundle_index)), read_job_summary,
    mc.cores = min(cfg$n_cores, nrow(bundle_index)), mc.preschedule = FALSE
  )
} else {
  lapply(seq_len(nrow(bundle_index)), read_job_summary)
}

# A try-error marks a failed read; NULL marks a worker process that stopped
# without returning, which mclapply reports only as a warning.
failed <- vapply(job_summaries, function(result) {
  is.null(result) || inherits(result, "try-error")
}, logical(1))
if (any(failed)) {
  first_error <- job_summaries[[which(failed)[[1]]]]
  error_detail <- if (inherits(first_error, "try-error")) {
    conditionMessage(attr(first_error, "condition"))
  } else {
    "a worker process stopped without returning a result"
  }
  stop(
    sum(failed), " compact summary file(s) could not be read.\n",
    "First error: ", error_detail
  )
}

monthly <- bind_rows(lapply(job_summaries, `[[`, "monthly"))
regime_months <- bind_rows(lapply(job_summaries, `[[`, "regimes"))

metrics <- c(
  "mean_degree", "mean_edge_distance_km", "p95_edge_distance_km",
  "long_edges_20"
)

# For each replicate, average monthly structure over the complete observed
# cycle or future 20-year window. Quantiles then describe variation among
# simulated networks, not pseudo-replication among months.
baseline_replicates <- monthly |>
  filter(dataset == "observed") |>
  group_by(resource, replicate) |>
  summarise(across(all_of(metrics), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  pivot_longer(all_of(metrics), names_to = "metric", values_to = "value")

baseline_reference <- baseline_replicates |>
  group_by(resource, metric) |>
  summarise(
    q10 = quantile(value, 0.10, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q90 = quantile(value, 0.90, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(
  baseline_reference,
  file.path(cfg$projection_summaries_dir, "observed_network_reference.csv")
)

future_replicates <- monthly |>
  filter(dataset == "future") |>
  group_by(scenario, scenario_label, climate_model, future_window,
           resource, replicate) |>
  summarise(across(all_of(metrics), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  pivot_longer(all_of(metrics), names_to = "metric", values_to = "value")

future_network_change <- future_replicates |>
  group_by(scenario, scenario_label, climate_model, future_window,
           resource, metric) |>
  summarise(
    q10 = quantile(value, 0.10, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q90 = quantile(value, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    baseline_reference |>
      select(resource, metric, observed_median = median),
    by = c("resource", "metric")
  ) |>
  mutate(
    q10_change_percent = 100 * (q10 / observed_median - 1),
    median_change_percent = 100 * (median / observed_median - 1),
    q90_change_percent = 100 * (q90 / observed_median - 1)
  )
write_csv(
  future_network_change,
  file.path(cfg$projection_summaries_dir, "future_network_change_by_gcm.csv")
)

ensemble_summary <- future_network_change |>
  group_by(scenario, scenario_label, future_window, resource, metric) |>
  summarise(
    ensemble_median_change_percent = median(median_change_percent),
    gcm_minimum_change_percent = min(median_change_percent),
    gcm_maximum_change_percent = max(median_change_percent),
    .groups = "drop"
  )
write_csv(
  ensemble_summary,
  file.path(cfg$projection_summaries_dir, "future_network_ensemble_summary.csv")
)

regime_composition <- regime_months |>
  mutate(
    comparison_period = if_else(dataset == "observed", "2009-2018", future_window),
    climate_model = if_else(dataset == "observed", "Observed", climate_model)
  ) |>
  count(dataset, scenario, scenario_label, climate_model,
        comparison_period, regime, name = "months") |>
  group_by(dataset, scenario, scenario_label, climate_model, comparison_period) |>
  mutate(share = months / sum(months)) |>
  ungroup()
write_csv(
  regime_composition,
  file.path(cfg$projection_summaries_dir, "regime_composition.csv")
)

scenario_colours <- c("ssp245" = "#0072B2", "ssp370" = "#D55E00")
gcm_shapes <- setNames(c(16, 17, 15, 18), cfg$climate_models)
window_levels <- names(cfg$future_windows)

# Use the same quiet styling for both projection figures. Horizontal guides
# support comparisons with zero, while compact side-by-side legends preserve
# space for the data panels.
projection_dodge <- position_dodge(width = 0.32)
movement_projection_theme <- theme_classic(base_size = 11.5) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.margin = margin(t = 3, r = 0, b = 0, l = 0),
    legend.spacing.x = grid::unit(0.25, "lines"),
    legend.key.width = grid::unit(1.0, "lines"),
    legend.text = element_text(size = 9.5),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", face = "bold"),
    panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = 0.35),
    panel.spacing = grid::unit(0.9, "lines"),
    plot.margin = margin(t = 4, r = 6, b = 2, l = 6)
  )

network_plot_data <- future_network_change |>
  filter(metric %in% c("mean_degree", "p95_edge_distance_km")) |>
  mutate(
    future_window = factor(future_window, levels = window_levels),
    scenario_label = factor(
      scenario_label,
      levels = unname(cfg$scenario_labels[cfg$scenarios])
    ),
    climate_model = factor(climate_model, levels = cfg$climate_models),
    resource = factor(
      resource, levels = c("grazing", "watering"),
      labels = c("Grazing network", "Watering network")
    ),
    metric = factor(
      metric,
      levels = c("mean_degree", "p95_edge_distance_km"),
      labels = c("Mean degree", "95th-percentile edge distance")
    )
  )

network_plot <- ggplot(
  network_plot_data,
  aes(
    x = future_window,
    y = median_change_percent,
    colour = scenario_label,
    shape = climate_model,
    group = interaction(scenario_label, climate_model)
  )
) +
  geom_hline(yintercept = 0, colour = "#666666", linewidth = 0.45) +
  geom_errorbar(
    aes(ymin = q10_change_percent, ymax = q90_change_percent),
    position = projection_dodge, width = 0.07, linewidth = 0.4, alpha = 0.8
  ) +
  geom_point(position = projection_dodge, size = 2.65) +
  facet_grid(metric ~ resource, scales = "free_y") +
  scale_colour_manual(
    values = setNames(scenario_colours[cfg$scenarios],
                      unname(cfg$scenario_labels[cfg$scenarios]))
  ) +
  scale_shape_manual(
    values = gcm_shapes,
    labels = gsub("_", "-", cfg$climate_models)
  ) +
  labs(
    x = "Future window",
    y = "Change from observed-cycle median (%)",
    colour = NULL,
    shape = NULL
  ) +
  guides(
    colour = guide_legend(
      order = 1, nrow = 1,
      override.aes = list(shape = 16, size = 2.65)
    ),
    shape = guide_legend(
      order = 2, nrow = 1,
      override.aes = list(colour = "#222222", size = 2.65)
    )
  ) +
  movement_projection_theme

ggsave(
  file.path(cfg$projection_figures_dir, "future_network_change.pdf"),
  network_plot, width = 9.2, height = 6.15, device = cairo_pdf
)

observed_regimes <- regime_composition |>
  filter(dataset == "observed") |>
  transmute(
    regime,
    observed_share = share
  )

regime_change_plot_data <- regime_composition |>
  filter(dataset == "future") |>
  left_join(observed_regimes, by = "regime") |>
  mutate(
    comparison_period = factor(comparison_period, levels = window_levels),
    scenario_label = factor(
      scenario_label,
      levels = unname(cfg$scenario_labels[cfg$scenarios])
    ),
    climate_model = factor(climate_model, levels = cfg$climate_models),
    regime = factor(
      regime,
      levels = c("routine", "transition", "drought"),
      labels = c("Routine", "Transition", "Drought")
    ),
    change_percentage_points = 100 * (share - observed_share)
  )

# Plot all three regime shares as changes from the observed cycle. A common
# zero reference is more direct than separate, unlabeled baseline lines.
regime_change_plot <- ggplot(
  regime_change_plot_data,
  aes(
    comparison_period, change_percentage_points,
    colour = scenario_label,
    shape = climate_model
  )
) +
  geom_hline(yintercept = 0, colour = "#666666", linewidth = 0.45) +
  geom_point(position = projection_dodge, size = 2.7) +
  facet_wrap(~regime, nrow = 1) +
  scale_colour_manual(
    values = setNames(scenario_colours[cfg$scenarios],
                      unname(cfg$scenario_labels[cfg$scenarios]))
  ) +
  scale_shape_manual(
    values = gcm_shapes,
    labels = gsub("_", "-", cfg$climate_models)
  ) +
  labs(
    x = "Future window",
    y = "Change from observed share\n(percentage points)",
    colour = NULL,
    shape = NULL
  ) +
  guides(
    colour = guide_legend(
      order = 1, nrow = 1,
      override.aes = list(shape = 16, size = 2.7)
    ),
    shape = guide_legend(
      order = 2, nrow = 1,
      override.aes = list(colour = "#222222", size = 2.7)
    )
  ) +
  movement_projection_theme

ggsave(
  file.path(cfg$projection_figures_dir, "future_regime_change.pdf"),
  regime_change_plot, width = 9.8, height = 4.35, device = cairo_pdf
)

message("Projection summaries written to: ", cfg$projection_summaries_dir)
message("Projection figures written to: ", cfg$projection_figures_dir)
