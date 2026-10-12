#!/usr/bin/env Rscript

# Build the remaining figures and the headline tables from the saved
# simulation outputs. This script does not rerun the transmission model.

module_dir <- here::here("analysis", "03_rvf_transmission")

run_stage <- function(path, label = basename(path)) {
  message("\nBuilding manuscript output: ", label)
  sys.source(file.path(module_dir, path), envir = new.env(parent = globalenv()))
}

run_stage(
  "R/figures/plot_climate_vector_functions.R",
  "temperature-trait and rainfall-response figures"
)
run_stage(
  "R/figures/summarise_vector_abundance.R",
  "model-implied vector-seasonality diagnostic"
)
run_stage(
  "R/figures/plot_attribution_scenario_change.R",
  "climate and movement attribution figure"
)
run_stage("R/workflows/manuscript/summarise_headline_results.R")

message("\nModule 03 manuscript outputs complete.")
