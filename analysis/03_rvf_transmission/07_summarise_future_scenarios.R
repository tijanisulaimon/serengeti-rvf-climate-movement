#!/usr/bin/env Rscript

# Summarise each configured future run over the common 40-year reporting
# period, then build the two-profile comparison and rainfall sensitivity.

source(here::here("analysis", "03_rvf_transmission", "R", "config.R"))
stage <- here::here(
  "analysis", "03_rvf_transmission", "R", "workflows", "summaries",
  "summarise_future_profile.R"
)

for (run in rvf_config$execution$future_runs) {
  message(
    "\nFuture summary: ", run$site_profile,
    " / ", run$rainfall_model
  )
  run_env <- new.env(parent = globalenv())
  run_env$analysis_profile <- run$site_profile
  run_env$analysis_rainfall_model <- run$rainfall_model
  run_env$analysis_network_replicates <- run$network_replicates
  sys.source(stage, envir = run_env)
}

message("\nAll configured future summaries are complete.")
