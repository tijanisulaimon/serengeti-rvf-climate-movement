#!/usr/bin/env Rscript

# Prepare all movement and climate inputs used by the transmission analyses.
# Each stage writes a documented intermediate product that is consumed by the
# next stage. Separate evaluation environments prevent temporary objects from
# one stage leaking into another.

module_dir <- here::here("analysis", "03_rvf_transmission")
stages <- c(
  "R/workflows/preparation/prepare_networks.R",
  "R/workflows/preparation/extract_site_climate.R",
  "R/workflows/preparation/assemble_scenario_bundles.R"
)

for (stage in stages) {
  message("\nPreparing Module 03 inputs: ", basename(stage))
  sys.source(file.path(module_dir, stage), envir = new.env(parent = globalenv()))
}

message("\nModule 03 input preparation complete.")
