#!/usr/bin/env Rscript

# Run the independently initialised mobility experiment for every calibrated
# exposure profile. The rainfall response remains the primary Hill formulation.

source(here::here("analysis", "03_rvf_transmission", "R", "config.R"))
stage <- here::here(
  "analysis", "03_rvf_transmission", "R", "workflows", "scenarios",
  "run_mobility_profile.R"
)

for (profile in rvf_config$execution$mobility_profiles) {
  message("\nMobility scenarios: ", profile)
  run_env <- new.env(parent = globalenv())
  run_env$analysis_profile <- profile
  sys.source(stage, envir = run_env)
}

message("\nMobility scenarios are complete for all configured profiles.")
