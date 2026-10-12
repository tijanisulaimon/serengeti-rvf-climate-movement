#!/usr/bin/env Rscript
# Define the three movement regimes and the one-mode calibration targets.
#
# Regime thresholds are derived once from the 2009-2018 resource-condition
# cycle and are held fixed for both SSPs. The 2017 participatory networks are
# observed only as village-to-village contacts, so all fit targets are computed
# at that same one-mode scale.

suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(readr)
  library(tidyr)
})

source(here::here("analysis", "02_movement_networks", "_config.R"))
source(here::here("analysis", "02_movement_networks", "R", "utils.R"))

cfg <- mvt_config
observed <- readRDS(cfg$observed_inputs_file)
distance_village <- observed$distance_village_village
resource_ndvi <- observed$ndvi[c("grazing", "watering")]

availability_parameters <- calibrate_resource_availability_params(
  resource_ndvi,
  years = cfg$observed_years,
  reference_months = cfg$availability_reference_months,
  target_unavailable = cfg$availability_target_unavailable,
  minimum = cfg$availability_floor,
  slope = cfg$availability_slope
)

availability_matrices <- lapply(names(resource_ndvi), function(resource) {
  resource_availability_cycle_mean_matrix(
    resource_ndvi[[resource]], cfg$observed_years,
    availability_parameters[[resource]]
  )
})
names(availability_matrices) <- names(resource_ndvi)

reference_availability <- combined_resource_availability_cycle_mean(
  resource_ndvi, cfg$observed_years, availability_parameters
)
thresholds <- calibrate_availability_regime_thresholds_cycle_current(
  resource_ndvi,
  availability_parameters,
  cfg$observed_years,
  drought_percentile = cfg$drought_percentile,
  transition_percentile = cfg$transition_percentile
)
regimes <- assign_regime_signal_three(reference_availability$combined, thresholds)

graph_summary <- function(graph) {
  distance <- edge_distance_vector(graph, distance_village)
  c(
    edges = igraph::ecount(graph),
    mean_dist = if (length(distance) == 0L) NA_real_ else mean(distance),
    p95_dist = if (length(distance) == 0L) NA_real_ else
      unname(stats::quantile(distance, 0.95)),
    long_edges_15 = sum(distance >= 15),
    long_edges_20 = sum(distance >= 20),
    long_edges_25 = sum(distance >= 25),
    long_edges_30 = sum(distance >= 30)
  )
}

summary_template <- stats::setNames(
  numeric(7),
  c(
    "edges", "mean_dist", "p95_dist", "long_edges_15",
    "long_edges_20", "long_edges_25", "long_edges_30"
  )
)

summarise_resource <- function(resource) {
  networks <- observed$observed_networks[[resource]]
  regime_summaries <- lapply(c("routine", "transition", "drought"), function(regime) {
    months <- names(regimes)[regimes == regime]
    monthly <- vapply(months, function(month) {
      graph_summary(networks[[month]])
    }, summary_template)
    rowMeans(monthly, na.rm = TRUE)
  })
  names(regime_summaries) <- c("routine", "transition", "drought")

  list(
    summaries = regime_summaries,
    monthly = networks,
    annual = graph_summary(union_one_mode_graphs(networks, observed$village_names))
  )
}

# `union_one_mode_graphs()` is defined in the simulator because simulated and
# observed annual contacts must use exactly the same temporal operation.
source(here::here("analysis", "02_movement_networks", "R", "simulator.R"))

targets <- list(
  grazing = summarise_resource("grazing"),
  watering = summarise_resource("watering"),
  regimes = regimes,
  thresholds = thresholds,
  availability_parameters = availability_parameters,
  observed_availability_matrices = availability_matrices,
  observed_years = cfg$observed_years,
  observation_scale = "monthly one-mode village networks"
)
saveRDS(targets, file.path(cfg$calibration_dir, "targets.rds"))

regime_reference <- tibble::tibble(
  month = month.abb,
  month_number = seq_along(month.abb),
  regime = unname(regimes[month.abb]),
  combined_availability = as.numeric(reference_availability$combined[month.abb]),
  grazing_availability = as.numeric(
    reference_availability$by_resource$grazing[month.abb]
  ),
  watering_availability = as.numeric(
    reference_availability$by_resource$watering[month.abb]
  ),
  drought_threshold = unname(thresholds[["drought"]]),
  transition_threshold = unname(thresholds[["transition"]])
)
readr::write_csv(
  regime_reference, file.path(cfg$calibration_dir, "regime_reference.csv")
)

target_table <- dplyr::bind_rows(lapply(c("grazing", "watering"), function(resource) {
  dplyr::bind_rows(lapply(names(targets[[resource]]$summaries), function(regime) {
    tibble::as_tibble_row(targets[[resource]]$summaries[[regime]]) |>
      dplyr::mutate(
        resource = resource,
        regime = regime,
        n_months = sum(regimes == regime),
        .before = 1
      )
  }))
}))
readr::write_csv(
  target_table, file.path(cfg$calibration_dir, "observed_network_targets.csv")
)

availability_table <- dplyr::bind_rows(lapply(names(availability_parameters), function(resource) {
  parameter <- availability_parameters[[resource]]
  tibble::tibble(
    resource = resource,
    midpoint = parameter$midpoint,
    slope = parameter$slope,
    minimum_availability = parameter$minimum,
    target_unavailable = parameter$target_unavailable,
    realised_unavailable = parameter$realised_reference_unavailable,
    reference_months = paste(parameter$reference_months, collapse = ", ")
  )
}))
readr::write_csv(
  availability_table,
  file.path(cfg$calibration_dir, "availability_parameters.csv")
)

within_regime_persistence <- function(resource) {
  networks <- observed$observed_networks[[resource]]
  dplyr::bind_rows(lapply(c("routine", "transition", "drought"), function(regime) {
    months <- names(regimes)[regimes == regime]
    pairs <- utils::combn(months, 2, simplify = FALSE)
    values <- vapply(pairs, function(pair) {
      jaccard(edge_keys(networks[[pair[1]]]), edge_keys(networks[[pair[2]]]))
    }, numeric(1))
    tibble::tibble(
      resource = resource,
      regime = regime,
      n_month_pairs = length(values),
      median_jaccard = median(values),
      minimum_jaccard = min(values),
      maximum_jaccard = max(values)
    )
  }))
}

persistence_table <- dplyr::bind_rows(
  within_regime_persistence("grazing"),
  within_regime_persistence("watering")
)
readr::write_csv(
  persistence_table,
  file.path(cfg$calibration_dir, "observed_edge_persistence.csv")
)

threshold_sensitivity <- dplyr::bind_rows(lapply(
  c(0.10, 0.12, 0.15, 0.18, 0.20),
  function(percentile) {
    candidate <- calibrate_availability_regime_thresholds_cycle_current(
      resource_ndvi,
      availability_parameters,
      cfg$observed_years,
      drought_percentile = percentile,
      transition_percentile = cfg$transition_percentile
    )
    assigned <- assign_regime_signal_three(
      reference_availability$combined, candidate
    )
    tibble::tibble(
      drought_percentile = percentile,
      drought_months = paste(names(assigned)[assigned == "drought"], collapse = ", "),
      drought_threshold = unname(candidate[["drought"]])
    )
  }
))
readr::write_csv(
  threshold_sensitivity,
  file.path(cfg$calibration_dir, "drought_threshold_sensitivity.csv")
)

message("Regimes: ", paste(month.abb, regimes, sep = "=", collapse = ", "))
message("Calibration targets written to: ", cfg$calibration_dir)
