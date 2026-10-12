#!/usr/bin/env Rscript
#' Index the Module 02 movement bundles used by the transmission model.
#'
#' Module 02 already stores every monthly bipartite network draw. This script
#' therefore writes only a compact index for the selected replicates; it does
#' not duplicate the matrices. Script 03 reads the indexed source bundle,
#' normalises the selected matrices, and joins them to site-level climate.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

source(here::here("analysis", "03_rvf_transmission", "R", "config.R"))
ensure_rvf_dirs()

cfg <- rvf_config

if (!file.exists(cfg$mvt_bundle_index)) {
  stop("Module 02 bundle index not found: ", cfg$mvt_bundle_index)
}

# A resumed projection run marks the bundles it skipped as "existing", so both
# statuses mean the bundle is complete.
bundle_index <- read_csv(cfg$mvt_bundle_index, show_col_types = FALSE) |>
  filter(.data$status %in% c("written", "existing")) |>
  filter(
    .data$dataset == "observed" |
      (.data$scenario %in% cfg$future_scenarios &
         .data$climate_model %in% cfg$climate_models)
  ) |>
  arrange(.data$dataset, .data$scenario, .data$climate_model, .data$year)

if (nrow(bundle_index) == 0L) {
  stop("No eligible Module 02 movement bundles were found in the bundle index.")
}

missing_files <- bundle_index$output_file[!file.exists(bundle_index$output_file)]
if (length(missing_files) > 0L) {
  stop(
    "The Module 02 bundle index lists missing bundle(s), for example: ",
    missing_files[[1]]
  )
}

normalise_observed_id <- function(scenario_id) {
  if (grepl("^OBSERVED_[0-9]{4}$", scenario_id)) {
    sub("^OBSERVED_", "OBSERVED_CURRENT_", scenario_id)
  } else {
    scenario_id
  }
}

sample_bundle <- readRDS(bundle_index$output_file[[1]])
required <- c("metadata", "vill_names", "grazing", "watering")
if (!all(required %in% names(sample_bundle))) {
  stop("Unexpected Module 02 bundle structure: ", bundle_index$output_file[[1]])
}

n_available <- min(
  length(sample_bundle$grazing$replicates),
  length(sample_bundle$watering$replicates)
)
selected <- cfg$network_replicates
if (is.null(selected)) selected <- seq_len(n_available)
selected <- sort(unique(as.integer(selected)))
if (any(!selected %in% seq_len(n_available))) {
  stop("Requested network replicate is outside the available range 1-", n_available)
}

n_vil <- length(sample_bundle$vill_names)
n_gra <- length(sample_bundle$grazing$cent_names)
n_wat <- length(sample_bundle$watering$cent_names)

index_rows <- lapply(seq_len(nrow(bundle_index)), function(i) {
  row <- bundle_index[i, ]
  base_id <- normalise_observed_id(row$scenario_id)
  scenario <- if (identical(row$dataset, "observed")) NA_character_ else row$scenario
  scenario_label <- if (identical(row$dataset, "observed")) {
    "Observed"
  } else {
    row$scenario_label
  }
  climate_model <- if (identical(row$dataset, "observed")) {
    "OBSERVED_CURRENT"
  } else {
    row$climate_model
  }

  tibble(
    scenario_id = paste0(base_id, "_REP", sprintf("%03d", selected)),
    scenario_base_id = base_id,
    dataset = row$dataset,
    scenario = scenario,
    scenario_label = scenario_label,
    climate_model = climate_model,
    future_window = row$future_window,
    year = as.integer(row$year),
    network_replicate = selected,
    n_replicates_available = n_available,
    n_vil = n_vil,
    n_gra = n_gra,
    n_wat = n_wat,
    source_bundle = row$output_file
  )
})

network_index <- bind_rows(index_rows) |>
  arrange(
    .data$dataset, .data$scenario, .data$climate_model,
    .data$year, .data$network_replicate
  )

index_file <- file.path(cfg$networks_dir, "network_scenario_index.csv")
write_csv(network_index, index_file)

message("Indexed ", nrow(network_index), " scenario-replicate combinations.")
message("Network index: ", index_file)
