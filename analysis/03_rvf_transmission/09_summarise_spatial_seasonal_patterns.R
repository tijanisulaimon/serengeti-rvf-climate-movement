#!/usr/bin/env Rscript

# Summarise where and when the model supports endemic transmission. These
# analyses read existing baseline and future outputs; they do not rerun the
# ODE system.

module_dir <- here::here("analysis", "03_rvf_transmission")
# The seasonal summary runs first because the spatial summary ends with a
# village map, which needs the village coordinates from the field data. Without
# them the run stops at the map, after both summaries have written their tables.
# The combined figure reads the spatial tables, so it runs last.
stages <- c(
  "R/workflows/summaries/summarise_seasonal_cycles.R",
  "R/workflows/summaries/summarise_spatial_patterns.R",
  "R/workflows/summaries/build_spatial_seasonal_figure.R"
)

for (stage in stages) {
  message("\nBuilding spatial/seasonal summary: ", basename(stage))
  sys.source(file.path(module_dir, stage), envir = new.env(parent = globalenv()))
}

message("\nSpatial and seasonal summaries complete.")
