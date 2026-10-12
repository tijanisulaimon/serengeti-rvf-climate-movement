#!/usr/bin/env Rscript
# Validate the fitted simulator against the observed one-mode networks.
#
# The validation asks whether replicate simulations reproduce monthly contact
# density, geographic reach, annual union structure, and within-regime edge
# persistence. It does not compare latent bipartite edges with observations,
# because the participatory data do not identify resource-site choices.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tidyr)
})

source(here::here("analysis", "02_movement_networks", "_config.R"))
source(here::here("analysis", "02_movement_networks", "R", "utils.R"))
source(here::here("analysis", "02_movement_networks", "R", "simulator.R"))

cfg <- mvt_config
options(movement_n_cores = cfg$n_cores)

observed <- readRDS(cfg$observed_inputs_file)
targets <- readRDS(file.path(cfg$calibration_dir, "targets.rds"))
fit <- readRDS(file.path(cfg$calibration_dir, "fit.rds"))
current_targets_digest <- unname(tools::md5sum(
  file.path(cfg$calibration_dir, "targets.rds")
))
if (!identical(fit$targets_digest, current_targets_digest)) {
  stop("The fitted model does not match the current calibration targets. ",
       "Run 03_fit_model.R before validation.")
}
if (!identical(
  fit$monthly_route_turnover_fraction,
  cfg$monthly_route_turnover_fraction
)) {
  stop("The fitted model predates the current monthly route-turnover setting. ",
       "Run 03_fit_model.R before validation.")
}
distance_village <- observed$distance_village_village
metric_names <- c(
  "edges", "mean_dist", "p95_dist", "long_edges_15", "long_edges_20"
)

# Summarise the one-mode contact features used for fitting and validation.
# Long-edge thresholds are retained here so monthly and annual summaries use
# exactly the same definitions.
graph_summary <- function(graph) {
  distance <- edge_distance_vector(graph, distance_village)
  tibble(
    edges = igraph::ecount(graph),
    mean_dist = if (length(distance) == 0L) NA_real_ else mean(distance),
    p95_dist = if (length(distance) == 0L) NA_real_ else
      unname(quantile(distance, 0.95)),
    long_edges_15 = sum(distance >= 15),
    long_edges_20 = sum(distance >= 20),
    long_edges_25 = sum(distance >= 25),
    long_edges_30 = sum(distance >= 30)
  )
}

# Extract the fitted parameters in the layer structure expected by the simulator.
selected_parameters <- function(resource) {
  lapply(fit$fits[[resource]], `[[`, "params")
}

# Generate new stochastic realizations under the observed climate cycle.
# These use different random seeds from fitting, but the same 2017 network data,
# so they assess model adequacy rather than out-of-sample predictive validity.
simulate_resource <- function(resource) {
  inputs <- build_resource_inputs(
    observed = observed,
    resource = resource,
    availability_matrix = targets$observed_availability_matrices[[resource]],
    regimes = targets$regimes,
    config = cfg
  )

  simulate_replicates(
    parameters = selected_parameters(resource),
    inputs = inputs,
    regimes = targets$regimes,
    n_reps = cfg$validation_replicates,
    seed = cfg$seeds$validation + match(resource, c("grazing", "watering")) * 10000L,
    n_cores = cfg$n_cores
  )
}

message("Simulating observed-cycle validation replicates...")
simulations <- list(
  grazing = simulate_resource("grazing"),
  watering = simulate_resource("watering")
)

simulated_monthly <- bind_rows(lapply(names(simulations), function(resource) {
  bind_rows(lapply(seq_along(simulations[[resource]]), function(replicate) {
    simulation <- simulations[[resource]][[replicate]]
    bind_rows(lapply(month.abb, function(month) {
      graph_summary(simulation$one_mode[[month]]) |>
        mutate(
          resource = resource,
          replicate = replicate,
          month = month,
          regime = unname(targets$regimes[[month]]),
          .before = 1
        )
    }))
  }))
}))

observed_monthly <- bind_rows(lapply(c("grazing", "watering"), function(resource) {
  bind_rows(lapply(month.abb, function(month) {
    graph_summary(observed$observed_networks[[resource]][[month]]) |>
      mutate(
        resource = resource,
        month = month,
        regime = unname(targets$regimes[[month]]),
        .before = 1
      )
  }))
}))

simulated_monthly_summary <- simulated_monthly |>
  pivot_longer(all_of(metric_names), names_to = "metric", values_to = "value") |>
  group_by(resource, month, regime, metric) |>
  summarise(
    q10 = quantile(value, 0.10, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q90 = quantile(value, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

monthly_validation <- observed_monthly |>
  pivot_longer(all_of(metric_names), names_to = "metric", values_to = "observed") |>
  left_join(
    simulated_monthly_summary,
    by = c("resource", "month", "regime", "metric")
  ) |>
  mutate(within_q10_q90 = observed >= q10 & observed <= q90)

write_csv(
  monthly_validation,
  file.path(cfg$validation_dir, "monthly_network_validation.csv")
)

# Regime summaries retain variation among replicate years rather than pooling
# all simulated months into a single artificial sample.
simulated_regime <- simulated_monthly |>
  group_by(resource, replicate, regime) |>
  summarise(across(all_of(metric_names), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  pivot_longer(all_of(metric_names), names_to = "metric", values_to = "value") |>
  group_by(resource, regime, metric) |>
  summarise(
    q10 = quantile(value, 0.10, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q90 = quantile(value, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

observed_regime <- observed_monthly |>
  group_by(resource, regime) |>
  summarise(across(all_of(metric_names), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  pivot_longer(all_of(metric_names), names_to = "metric", values_to = "observed")

regime_validation <- left_join(
  observed_regime, simulated_regime,
  by = c("resource", "regime", "metric")
)
write_csv(
  regime_validation,
  file.path(cfg$validation_dir, "regime_network_validation.csv")
)

simulated_annual <- bind_rows(lapply(names(simulations), function(resource) {
  bind_rows(lapply(seq_along(simulations[[resource]]), function(replicate) {
    graph_summary(simulations[[resource]][[replicate]]$year_net) |>
      mutate(resource = resource, replicate = replicate, .before = 1)
  }))
}))

observed_annual <- bind_rows(lapply(c("grazing", "watering"), function(resource) {
  graph_summary(union_one_mode_graphs(
    observed$observed_networks[[resource]], observed$village_names
  )) |>
    mutate(resource = resource, .before = 1)
}))

annual_validation <- simulated_annual |>
  pivot_longer(-c(resource, replicate), names_to = "metric", values_to = "value") |>
  group_by(resource, metric) |>
  summarise(
    q10 = quantile(value, 0.10, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q90 = quantile(value, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    observed_annual |>
      pivot_longer(-resource, names_to = "metric", values_to = "observed"),
    by = c("resource", "metric")
  )
write_csv(
  annual_validation,
  file.path(cfg$validation_dir, "annual_network_validation.csv")
)

# Identify all edges, or only long edges, before comparing within-regime overlap.
edge_keys_above <- function(graph, minimum_distance) {
  keys <- edge_keys(graph)
  if (minimum_distance <= 0 || length(keys) == 0L) return(keys)
  keys[edge_distance_vector(graph, distance_village) >= minimum_distance]
}

regime_months <- split(names(targets$regimes), unname(targets$regimes))

# Average pairwise Jaccard similarity across months in one regime. Reporting
# every regime reveals excessive persistence as well as insufficient overlap.
mean_pair_jaccard <- function(networks, months, minimum_distance) {
  if (length(months) < 2L) return(NA_real_)
  month_pairs <- utils::combn(months, 2, simplify = FALSE)
  values <- vapply(month_pairs, function(pair) {
    jaccard(
      edge_keys_above(networks[[pair[1]]], minimum_distance),
      edge_keys_above(networks[[pair[2]]], minimum_distance)
    )
  }, numeric(1))
  mean(values)
}

simulated_persistence <- bind_rows(lapply(names(simulations), function(resource) {
  bind_rows(lapply(seq_along(simulations[[resource]]), function(replicate) {
    current_networks <- simulations[[resource]][[replicate]]$one_mode
    bind_rows(lapply(names(regime_months), function(regime) {
      jaccard_values <- vapply(c(0, 15, 20), function(distance) {
        mean_pair_jaccard(
          current_networks,
          regime_months[[regime]],
          distance
        )
      }, numeric(1))
      tibble(
        resource = resource,
        replicate = replicate,
        regime = regime,
        minimum_edge_distance_km = c(0, 15, 20),
        jaccard = jaccard_values
      )
    }))
  }))
}))

observed_persistence <- bind_rows(lapply(c("grazing", "watering"), function(resource) {
  current_networks <- observed$observed_networks[[resource]]
  bind_rows(lapply(names(regime_months), function(regime) {
    jaccard_values <- vapply(c(0, 15, 20), function(distance) {
      mean_pair_jaccard(
        current_networks,
        regime_months[[regime]],
        distance
      )
    }, numeric(1))
    tibble(
      resource = resource,
      regime = regime,
      minimum_edge_distance_km = c(0, 15, 20),
      observed = jaccard_values
    )
  }))
}))

persistence_validation <- simulated_persistence |>
  group_by(resource, regime, minimum_edge_distance_km) |>
  summarise(
    q10 = quantile(jaccard, 0.10, na.rm = TRUE),
    median = median(jaccard, na.rm = TRUE),
    q90 = quantile(jaccard, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(observed_persistence,
            by = c("resource", "regime", "minimum_edge_distance_km")) |>
  mutate(within_q10_q90 = observed >= q10 & observed <= q90)
write_csv(
  persistence_validation,
  file.path(cfg$validation_dir, "edge_persistence_validation.csv")
)

# Retain the established drought-only filename for downstream manuscript code
# while the complete table provides the broader model diagnostic.
write_csv(
  filter(persistence_validation, regime == "drought"),
  file.path(cfg$validation_dir, "drought_edge_persistence_validation.csv")
)

resource_labels <- c(grazing = "Grazing", watering = "Watering")
metric_labels <- c(
  edges = "Monthly edges",
  p95_dist = "95th-percentile distance (km)"
)

# Monthly panels retain the two quantities that most directly describe contact
# density and spatial reach. The complete five-metric comparison remains in
# monthly_network_validation.csv.
monthly_plot <- monthly_validation |>
  filter(metric %in% names(metric_labels)) |>
  mutate(
    month = factor(month, levels = month.abb),
    metric = factor(metric, levels = names(metric_labels), labels = metric_labels),
    resource = factor(resource, levels = names(resource_labels), labels = resource_labels)
  ) |>
  ggplot(aes(x = month)) +
  geom_ribbon(aes(ymin = q10, ymax = q90, group = 1),
              fill = "#56B4E9", alpha = 0.22) +
  geom_line(aes(y = median, colour = "Simulated", group = 1), linewidth = 0.75) +
  geom_line(aes(y = observed, colour = "Observed", group = 1),
            linewidth = 0.55) +
  geom_point(aes(y = observed, colour = "Observed"), size = 1.65) +
  facet_grid(metric ~ resource, scales = "free_y") +
  scale_colour_manual(
    values = c("Observed" = "#111111", "Simulated" = "#0072B2"),
    breaks = c("Observed", "Simulated")
  ) +
  labs(x = NULL, y = NULL, colour = NULL, tag = "a") +
  theme_classic(base_size = 10.5) +
  theme(
    legend.position = "bottom",
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(size = 8.5),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", face = "bold"),
    panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = 0.3)
  )

# Annual unions test whether monthly route persistence produces a realistic
# number of distinct village contacts over the full year.
annual_plot <- annual_validation |>
  filter(metric == "edges") |>
  mutate(
    resource = factor(resource, levels = names(resource_labels),
                      labels = resource_labels)
  ) |>
  ggplot(aes(resource, median)) +
  geom_linerange(aes(ymin = q10, ymax = q90),
                 colour = "#0072B2", linewidth = 1.0) +
  geom_point(colour = "#0072B2", size = 2.4) +
  geom_point(aes(y = observed), shape = 23, fill = "white",
             colour = "#111111", stroke = 0.9, size = 2.7) +
  labs(x = NULL, y = "Annual-union edges", tag = "b") +
  theme_classic(base_size = 10.5) +
  theme(
    axis.text = element_text(colour = "black"),
    panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = 0.3)
  )

# Jaccard similarity is shown as a diagnostic rather than a fitted target.
# The contrast between observed diamonds and simulated intervals makes any
# remaining over-persistence visible instead of hiding it in aggregate counts.
persistence_plot <- persistence_validation |>
  filter(minimum_edge_distance_km == 0) |>
  mutate(
    resource = factor(resource, levels = names(resource_labels),
                      labels = resource_labels),
    regime = factor(regime, levels = c("routine", "transition", "drought"),
                    labels = c("Routine", "Transition", "Drought"))
  ) |>
  ggplot(aes(regime, median)) +
  geom_linerange(aes(ymin = q10, ymax = q90),
                 colour = "#0072B2", linewidth = 0.9) +
  geom_point(colour = "#0072B2", size = 2.2) +
  geom_point(aes(y = observed), shape = 23, fill = "white",
             colour = "#111111", stroke = 0.85, size = 2.5) +
  facet_wrap(~resource, nrow = 1) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Pairwise Jaccard similarity", tag = "c") +
  theme_classic(base_size = 10.5) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(size = 8.5, angle = 30, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = 0.3)
  )

structure_plot <- annual_plot + persistence_plot +
  plot_layout(widths = c(0.85, 1.15))

validation_plot <- monthly_plot / structure_plot +
  plot_layout(heights = c(1.65, 1), guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(face = "bold", size = 11)
  )

ggsave(
  file.path(cfg$validation_dir, "movement_model_validation.pdf"),
  validation_plot, width = 8.2, height = 7.8, device = cairo_pdf
)

regime_reference <- read_csv(
  file.path(cfg$calibration_dir, "regime_reference.csv"),
  show_col_types = FALSE
)
regime_plot_data <- bind_rows(
  regime_reference |>
    transmute(
      month, regime, panel = "Combined resource availability",
      value = combined_availability
    ),
  observed_monthly |>
    transmute(
      month, regime,
      panel = if_else(resource == "grazing",
                      "Observed grazing contacts", "Observed watering contacts"),
      value = edges
    )
) |>
  mutate(
    month = factor(month, levels = month.abb),
    regime = factor(regime, levels = c("routine", "transition", "drought")),
    panel = factor(
      panel,
      levels = c(
        "Combined resource availability",
        "Observed grazing contacts",
        "Observed watering contacts"
      )
    )
  )

threshold_lines <- tibble(
  panel = factor(
    rep("Combined resource availability", 2),
    levels = levels(regime_plot_data$panel)
  ),
  threshold = unname(targets$thresholds)
)

# Show that regimes are assigned from resource availability, then compare that
# classification with the observed seasonal contact pattern.
regime_definition_plot <- ggplot(regime_plot_data, aes(month, value, group = 1)) +
  geom_hline(
    data = threshold_lines,
    aes(yintercept = threshold),
    colour = "#555555", linetype = 3, linewidth = 0.55
  ) +
  geom_line(colour = "#222222", linewidth = 0.65) +
  geom_point(aes(fill = regime), shape = 21, colour = "white", stroke = 0.4,
             size = 2.8) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  scale_fill_manual(
    values = c(routine = "#009E73", transition = "#E69F00", drought = "#D55E00"),
    labels = c(routine = "Routine", transition = "Transition", drought = "Drought")
  ) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic(base_size = 10.5) +
  theme(
    legend.position = "bottom",
    axis.text = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", face = "bold"),
    panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = 0.3)
  )

ggsave(
  file.path(cfg$validation_dir, "movement_regime_definition.pdf"),
  regime_definition_plot, width = 7.5, height = 6.8, device = cairo_pdf
)

message("Validation outputs written to: ", cfg$validation_dir)
