#!/usr/bin/env Rscript

# Run the climate x movement factorial for the two contrasting Hill-response
# exposure profiles specified in R/config.R.

source(here::here("analysis", "03_rvf_transmission", "R", "config.R"))
stage <- here::here(
  "analysis", "03_rvf_transmission", "R", "workflows", "scenarios",
  "run_attribution_profile.R"
)

for (profile in rvf_config$execution$attribution_profiles) {
  message("\nClimate-movement attribution: ", profile)
  run_env <- new.env(parent = globalenv())
  run_env$analysis_profile <- profile
  sys.source(stage, envir = run_env)
}

message("\nAttribution is complete for all configured profiles.")
