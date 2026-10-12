#!/usr/bin/env Rscript

# Run the observed-climate baseline required by the main future profiles and
# the quadratic rainfall-response sensitivity. Each run uses its own calibrated
# capacity scale and writes to a profile/response-specific location.

source(here::here("analysis", "03_rvf_transmission", "R", "config.R"))
stage <- here::here(
  "analysis", "03_rvf_transmission", "R", "workflows", "scenarios",
  "run_observed_baseline.R"
)

for (run in rvf_config$execution$baseline_runs) {
  message(
    "\nObserved baseline: ", run$site_profile,
    " / ", run$rainfall_model
  )
  run_env <- new.env(parent = globalenv())
  run_env$analysis_profile <- run$site_profile
  run_env$analysis_rainfall_model <- run$rainfall_model
  run_env$analysis_network_replicates <- run$network_replicates
  sys.source(stage, envir = run_env)
}

message("\nAll configured observed baselines are complete.")
