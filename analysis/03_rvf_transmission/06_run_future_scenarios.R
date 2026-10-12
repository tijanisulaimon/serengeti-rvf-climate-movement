#!/usr/bin/env Rscript

# Run the two contrasting Hill-response exposure profiles and the quadratic
# rainfall-response sensitivity specified in R/config.R.

source(here::here("analysis", "03_rvf_transmission", "R", "config.R"))
stage <- here::here(
  "analysis", "03_rvf_transmission", "R", "workflows", "scenarios",
  "run_future_profile.R"
)

for (run in rvf_config$execution$future_runs) {
  message(
    "\nFuture scenarios: ", run$site_profile,
    " / ", run$rainfall_model
  )
  run_env <- new.env(parent = globalenv())
  run_env$analysis_profile <- run$site_profile
  run_env$analysis_rainfall_model <- run$rainfall_model
  run_env$analysis_network_replicates <- run$network_replicates
  sys.source(stage, envir = run_env)
}

message("\nAll configured future scenarios are complete.")
