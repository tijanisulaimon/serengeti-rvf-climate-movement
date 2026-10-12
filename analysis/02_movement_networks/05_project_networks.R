#!/usr/bin/env Rscript
# Generate observed and future bipartite movement-network bundles.
#
# Each output contains the monthly binary village-resource matrices required
# by Module 03, plus compact one-mode summaries used for Module 02 figures.
# Completed files are checked and skipped, so the 490 year-level jobs can be
# resumed safely after an interruption.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

source(here::here("analysis", "02_movement_networks", "_config.R"))
source(here::here("analysis", "02_movement_networks", "R", "utils.R"))
source(here::here("analysis", "02_movement_networks", "R", "simulator.R"))

cfg <- mvt_config
observed <- readRDS(cfg$observed_inputs_file)
fit_file <- file.path(cfg$calibration_dir, "fit.rds")
fit <- readRDS(fit_file)
fit_digest <- unname(tools::md5sum(fit_file))
current_targets_digest <- unname(tools::md5sum(
  file.path(cfg$calibration_dir, "targets.rds")
))
if (!identical(fit$targets_digest, current_targets_digest)) {
  stop("The fitted model does not match the current calibration targets. ",
       "Run 03_fit_model.R before projection.")
}
if (!identical(
  fit$monthly_route_turnover_fraction,
  cfg$monthly_route_turnover_fraction
)) {
  stop("The fitted model predates the current monthly route-turnover setting. ",
       "Run 03_fit_model.R before projection.")
}
distance_village <- observed$distance_village_village

jobs <- bind_rows(
  tibble(
    dataset = "observed",
    scenario = NA_character_,
    scenario_label = "Observed",
    climate_model = NA_character_,
    year = cfg$observed_years
  ),
  # igraph also exports a function named crossing(); qualify the tidyr
  # constructor because utils.R loads igraph after the package imports above.
  tidyr::crossing(
    scenario = cfg$scenarios,
    climate_model = cfg$climate_models,
    year = cfg$future_years
  ) |>
    mutate(
      dataset = "future",
      scenario_label = unname(cfg$scenario_labels[scenario]),
      .before = 1
    )
) |>
  mutate(
    job_id = row_number(),
    future_window = vapply(year, function(value) {
      match <- names(cfg$future_windows)[vapply(
        cfg$future_windows, function(years) value %in% years, logical(1)
      )]
      if (length(match) == 0L) NA_character_ else match[[1]]
    }, character(1)),
    output_file = mapply(
      network_bundle_file,
      dataset = dataset,
      year = year,
      scenario = scenario,
      climate_model = climate_model,
      USE.NAMES = FALSE
    ),
    summary_file = mapply(
      network_job_summary_file,
      dataset = dataset,
      year = year,
      scenario = scenario,
      climate_model = climate_model,
      USE.NAMES = FALSE
    )
  )

# Extract fitted layer parameters in the structure expected by the simulator.
selected_parameters <- function(resource) {
  lapply(fit$fits[[resource]], `[[`, "params")
}

scenario_identifier <- function(job) {
  if (job$dataset == "observed") return(sprintf("OBSERVED_%d", job$year))
  sprintf("%s_%s_%d", toupper(job$scenario), job$climate_model, job$year)
}

# Excluding the SSP from this seed gives matched stochastic draws for the same
# GCM, year, and resource under SSP2-4.5 and SSP3-7.0.
projection_seed <- function(job, resource) {
  model_index <- if (job$dataset == "observed") 0L else
    match(job$climate_model, cfg$climate_models)
  resource_offset <- if (resource == "watering") 5000000L else 0L
  as.integer(
    cfg$seeds$projection + model_index * 100000L + job$year * 10L +
      resource_offset
  )
}

# Load observed or SSP-GCM site-level NDVI for a single projection job.
load_ndvi <- function(job) {
  if (job$dataset == "observed") {
    return(observed$ndvi[c("grazing", "watering")])
  }
  readRDS(scenario_input_file(job$scenario, job$climate_model))$ndvi[
    c("grazing", "watering")
  ]
}

# Convert site NDVI to availability and classify the year's twelve regimes.
regime_information <- function(resource_ndvi, year) {
  availability <- combined_resource_availability(
    resource_ndvi,
    year,
    fit$availability_parameters
  )
  regimes <- assign_regime_signal_three(
    availability$combined,
    fit$thresholds
  )

  list(
    regimes = regimes,
    combined_availability = availability$combined,
    resource_availability = availability$by_resource
  )
}

# Compact one-mode summaries retained for Module 02 projection figures.
graph_summary <- function(graph) {
  distance <- edge_distance_vector(graph, distance_village)
  tibble(
    edges = igraph::ecount(graph),
    mean_degree = mean(igraph::degree(graph)),
    mean_edge_distance_km = if (length(distance) == 0L) NA_real_ else
      mean(distance),
    p95_edge_distance_km = if (length(distance) == 0L) NA_real_ else
      unname(quantile(distance, 0.95)),
    long_edges_20 = sum(distance >= 20)
  )
}

summarise_simulations <- function(simulations, resource, regime_info) {
  monthly <- bind_rows(lapply(seq_along(simulations), function(replicate) {
    bind_rows(lapply(month.abb, function(month) {
      graph_summary(simulations[[replicate]]$one_mode[[month]]) |>
        mutate(
          resource = resource,
          replicate = replicate,
          month = month,
          regime = unname(regime_info$regimes[[month]]),
          combined_availability = unname(
            regime_info$combined_availability[[month]]
          ),
          resource_availability = unname(
            regime_info$resource_availability[[resource]][[month]]
          ),
          .before = 1
        )
    }))
  }))

  annual <- bind_rows(lapply(seq_along(simulations), function(replicate) {
    graph_summary(simulations[[replicate]]$year_net) |>
      mutate(resource = resource, replicate = replicate, .before = 1)
  }))

  list(monthly = monthly, annual = annual)
}

# Skip a completed job only when its fit, replicate count, and scenario match.
outputs_are_current <- function(bundle_path, summary_path, job) {
  if (!file.exists(bundle_path) || !file.exists(summary_path)) return(FALSE)
  summary <- tryCatch(readRDS(summary_path), error = function(error) NULL)
  if (is.null(summary) || is.null(summary$metadata)) return(FALSE)

  metadata <- summary$metadata
  identical(metadata$fit_digest, fit_digest) &&
    identical(as.integer(metadata$n_replicates), cfg$projection_replicates) &&
    identical(as.integer(metadata$year), as.integer(job$year)) &&
    identical(metadata$dataset, job$dataset) &&
    identical(metadata$scenario %||% NA_character_, job$scenario %||% NA_character_) &&
    identical(metadata$climate_model %||% NA_character_,
              job$climate_model %||% NA_character_)
}

# Write atomically so an interrupted process cannot leave a valid-looking file.
save_object <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary, compress = TRUE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not move completed bundle to: ", path)
  }
}

# Simulate both resource networks for one year and save the Module 03 bundle
# alongside a compact, resumable summary used only for Module 02 aggregation.
project_job <- function(index) {
  job <- jobs[index, ]
  output_file <- job$output_file
  summary_file <- job$summary_file

  if (outputs_are_current(output_file, summary_file, job)) {
    return(tibble(
      job_id = job$job_id,
      scenario_id = scenario_identifier(job),
      output_file = output_file,
      summary_file = summary_file,
      status = "existing"
    ))
  }

  message(sprintf("[%d/%d] %s", index, nrow(jobs), scenario_identifier(job)))
  resource_ndvi <- load_ndvi(job)
  regime_info <- regime_information(resource_ndvi, job$year)

  bundle <- list(
    metadata = list(
      scenario_id = scenario_identifier(job),
      dataset = job$dataset,
      scenario = job$scenario,
      scenario_label = job$scenario_label,
      climate_model = job$climate_model,
      future_window = job$future_window,
      year = as.integer(job$year),
      regimes = regime_info$regimes,
      thresholds = fit$thresholds,
      combined_availability = regime_info$combined_availability,
      n_replicates = cfg$projection_replicates,
      fit_digest = fit_digest,
      common_random_numbers_across_ssps = TRUE
    ),
    vill_names = observed$village_names,
    grazing = list(
      cent_names = observed$grazing$resource_names,
      replicates = vector("list", cfg$projection_replicates)
    ),
    watering = list(
      cent_names = observed$watering$resource_names,
      replicates = vector("list", cfg$projection_replicates)
    )
  )

  monthly_summaries <- list()
  annual_summaries <- list()

  for (resource in c("grazing", "watering")) {
    availability_matrix <- resource_availability_matrix(
      resource_ndvi[[resource]],
      job$year,
      fit$availability_parameters[[resource]]
    )
    inputs <- build_resource_inputs(
      observed = observed,
      resource = resource,
      availability_matrix = availability_matrix,
      regimes = regime_info$regimes,
      config = cfg
    )
    simulations <- simulate_replicates(
      parameters = selected_parameters(resource),
      inputs = inputs,
      regimes = regime_info$regimes,
      n_reps = cfg$projection_replicates,
      seed = projection_seed(job, resource),
      n_cores = 1L
    )

    bundle[[resource]]$replicates <- lapply(simulations, `[[`, "bipartite")
    summary <- summarise_simulations(simulations, resource, regime_info)
    monthly_summaries[[resource]] <- summary$monthly
    annual_summaries[[resource]] <- summary$annual
  }

  compact_summary <- list(
    metadata = bundle$metadata,
    monthly = bind_rows(monthly_summaries),
    annual = bind_rows(annual_summaries)
  )
  save_object(bundle, output_file)
  save_object(compact_summary, summary_file)

  tibble(
    job_id = job$job_id,
    scenario_id = scenario_identifier(job),
    output_file = output_file,
    summary_file = summary_file,
    status = "written"
  )
}

message(
  "Projecting ", nrow(jobs), " year-level scenarios with ",
  cfg$projection_replicates, " replicates each on up to ", cfg$n_cores,
  " cores. Existing current bundles will be skipped."
)

job_results <- if (.Platform$OS.type == "unix" && cfg$n_cores > 1L) {
  parallel::mclapply(
    seq_len(nrow(jobs)), project_job,
    mc.cores = min(cfg$n_cores, nrow(jobs)), mc.preschedule = FALSE
  )
} else {
  lapply(seq_len(nrow(jobs)), project_job)
}

# A try-error marks a failed job; NULL marks a worker process that stopped
# without returning, which mclapply reports only as a warning.
failed <- vapply(job_results, function(result) {
  is.null(result) || inherits(result, "try-error")
}, logical(1))
if (any(failed)) {
  stop(sum(failed), " projection job(s) failed or their worker process stopped. ",
       "Re-run this script after checking the messages; completed bundles ",
       "will be skipped.")
}

bundle_index <- bind_rows(job_results) |>
  left_join(
    jobs |>
      select(job_id, dataset, scenario, scenario_label, climate_model,
             future_window, year),
    by = "job_id"
  ) |>
  mutate(
    fit_digest = fit_digest,
    n_replicates = cfg$projection_replicates
  ) |>
  arrange(dataset, scenario, climate_model, year)
write_csv(
  bundle_index,
  file.path(cfg$projection_summaries_dir, "network_bundle_index.csv")
)

message("Projection bundles written to: ", cfg$bundles_dir)
