#!/usr/bin/env Rscript
#' ####################################################################### #
#' 05 - Fit Final sdmTMB Model and Project Future NDVI
#'
#' Fit the selected NDVI model once to all observed data, then use that same
#' fit for every configured SSP-GCM projection. The fitted year smooth is held
#' at its 2018 value so projected change is driven by climate covariates rather
#' than extrapolation of an historical time trend.
#'
#' INPUTS:  derived/gimms_support/all_historical_observed.feather
#'          derived/gimms_support/future/{scenario}/all_future_{GCM}.feather
#'          outputs/selection/sdmtmb_best_model_observed.txt
#' OUTPUTS: outputs/models/final/model_fit.rds
#'          outputs/models/final/scaling.csv
#'          outputs/scenarios/{scenario}/projections/{GCM}/monthly_ndvi.feather
#'          outputs/scenarios/projection_index.csv
#' ####################################################################### #

pacman::p_load(
  "arrow",
  "tidyverse",
  "lubridate",
  "sf",
  "slider",
  "sdmTMB"
)

source(here::here("analysis", "01_ndvi_prediction", "_config.R"))

# --- Model covariates --------------------------------------------------------

add_spatial_coordinates <- function(df) {
  coords_sf <- df |>
    distinct(x, y) |>
    st_as_sf(coords = c("x", "y"), crs = 4326, remove = FALSE) |>
    st_transform(32736)
  coords_mat <- st_coordinates(coords_sf) / 1000

  coords_sf |>
    st_drop_geometry() |>
    mutate(
      x_utm_km = coords_mat[, 1],
      y_utm_km = coords_mat[, 2]
    )
}

add_historical_covariates <- function(df) {
  coords <- add_spatial_coordinates(df)

  df |>
    mutate(
      date = as.Date(date),
      year = year(date),
      month_num = month(date),
      month = factor(month_num, levels = 1:12),
      ele_km = ele / 1000
    ) |>
    left_join(coords, by = c("x", "y")) |>
    arrange(x, y, date) |>
    group_by(x, y) |>
    mutate(
      pr_lag1 = lag(pr),
      pr_roll3 = slide_dbl(pr, mean, .before = 2, .complete = TRUE),
      # Candidate models were compared on a common complete-case sample that
      # also supported a six-month rainfall term. Retain that same sample for
      # the final fit even though the selected formula does not use pr_roll6.
      pr_roll6 = slide_dbl(pr, mean, .before = 5, .complete = TRUE),
      tasmax_roll3 = slide_dbl(tasmax, mean, .before = 2, .complete = TRUE)
    ) |>
    ungroup() |>
    mutate(ndvi_beta = pmin(pmax(ndvi, 1e-4), 0.9999)) |>
    drop_na(
      ndvi_beta,
      pr_lag1,
      pr_roll3,
      pr_roll6,
      tasmax_roll3,
      ele_km,
      x_utm_km,
      y_utm_km
    )
}

add_future_covariates <- function(df) {
  coords <- add_spatial_coordinates(df)

  df |>
    mutate(
      date = as.Date(date),
      year = year(date),
      month_num = month(date),
      month = factor(month_num, levels = 1:12),
      ele_km = ele / 1000
    ) |>
    left_join(coords, by = c("x", "y")) |>
    arrange(x, y, date) |>
    group_by(x, y) |>
    mutate(
      pr_lag1 = lag(pr),
      pr_roll3 = slide_dbl(pr, mean, .before = 2, .complete = TRUE),
      pr_roll6 = slide_dbl(pr, mean, .before = 5, .complete = TRUE),
      tasmax_roll3 = slide_dbl(tasmax, mean, .before = 2, .complete = TRUE)
    ) |>
    ungroup() |>
    drop_na(
      pr_lag1,
      pr_roll3,
      pr_roll6,
      tasmax_roll3,
      ele_km,
      x_utm_km,
      y_utm_km
    )
}

apply_scaling <- function(df, scaling) {
  centres <- setNames(scaling$mean, scaling$variable)
  spreads <- setNames(scaling$sd, scaling$variable)

  df |>
    mutate(
      pr_lag1_sc = (pr_lag1 - centres[["pr_lag1"]]) /
        spreads[["pr_lag1"]],
      pr_roll3_sc = (pr_roll3 - centres[["pr_roll3"]]) /
        spreads[["pr_roll3"]],
      tasmax_roll3_sc = (tasmax_roll3 - centres[["tasmax_roll3"]]) /
        spreads[["tasmax_roll3"]],
      ele_km_sc = (ele_km - centres[["ele_km"]]) /
        spreads[["ele_km"]],
      year_sc = (year - centres[["year"]]) /
        spreads[["year"]]
    )
}

# --- Historical fit ---------------------------------------------------------

dir_data <- ndvi_config$gimms_support_dir
dir_outputs <- ndvi_config$outputs_dir
dir_selection <- file.path(dir_outputs, "selection")
dir_final_model <- file.path(dir_outputs, "models", "final")
dir_scenarios <- ndvi_config$scenario_outputs_root
dir.create(dir_final_model, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_scenarios, recursive = TRUE, showWarnings = FALSE)

historical_path <- file.path(dir_data, "all_historical_observed.feather")
selected_model_path <- file.path(
  dir_selection,
  "sdmtmb_best_model_observed.txt"
)

if (!file.exists(historical_path)) {
  stop("Missing historical feature table: ", historical_path)
}
if (!file.exists(selected_model_path)) {
  stop("Missing selected-model file: ", selected_model_path)
}

selected_model <- readLines(selected_model_path, warn = FALSE)[1]
if (!identical(selected_model, "sdmtmb_winner_t2")) {
  stop(
    "Script 05 implements the selected publication model ",
    "'sdmtmb_winner_t2', but script 03 selected '", selected_model, "'. ",
    "Review the model-selection result before projecting."
  )
}

cat("Preparing historical data for the final NDVI model...\n")
historical <- read_feather(historical_path) |>
  add_historical_covariates()

scaling_variables <- c(
  "pr_lag1",
  "pr_roll3",
  "tasmax_roll3",
  "ele_km",
  "year"
)
scaling_current <- tibble(
  variable = scaling_variables,
  mean = vapply(historical[scaling_variables], mean, numeric(1)),
  sd = vapply(historical[scaling_variables], sd, numeric(1))
)

fit_path <- file.path(dir_final_model, "model_fit.rds")
scaling_path <- file.path(dir_final_model, "scaling.csv")

if (file.exists(fit_path)) {
  if (!file.exists(scaling_path)) {
    stop(
      "A saved final fit exists without its scaling table: ", fit_path,
      "\nRestore the matching scaling.csv or refit the model."
    )
  }

  scaling_saved <- read_csv(scaling_path, show_col_types = FALSE) |>
    filter(variable %in% scaling_variables) |>
    arrange(match(variable, scaling_variables))
  scaling_check <- scaling_current |>
    arrange(match(variable, scaling_variables))

  scaling_matches <- identical(
    scaling_saved$variable,
    scaling_check$variable
  ) &&
    max(abs(scaling_saved$mean - scaling_check$mean)) < 1e-10 &&
    max(abs(scaling_saved$sd - scaling_check$sd)) < 1e-10

  if (!scaling_matches) {
    stop(
      "The historical scaling values differ from those used by the saved ",
      "final fit. Refit before projecting."
    )
  }

  scaling <- scaling_saved
  fit <- readRDS(fit_path)
  cat("Reusing the verified final model fit.\n")
} else {
  scaling <- scaling_current
  historical_scaled <- apply_scaling(historical, scaling)
  mesh <- make_mesh(
    historical_scaled,
    xy_cols = c("x_utm_km", "y_utm_km"),
    cutoff = 20
  )

  set.seed(20260507)
  fit <- sdmTMB(
    ndvi_beta ~
      s(pr_lag1_sc, k = 5) +
      s(ele_km_sc, k = 4) +
      t2(pr_roll3_sc, tasmax_roll3_sc, k = c(5, 5)) +
      s(year_sc, k = 5) +
      month,
    data = historical_scaled,
    mesh = mesh,
    family = Beta(link = "logit"),
    spatial = "on",
    spatiotemporal = "off",
    reml = TRUE,
    silent = TRUE
  )

  saveRDS(fit, fit_path)
  write_csv(scaling, scaling_path)
  cat("Saved one final model fit for all future projections.\n")
}

historical_scaled <- apply_scaling(historical, scaling)
fixed_year <- max(historical_scaled$year)
fixed_year_sc <- (fixed_year - scaling$mean[scaling$variable == "year"]) /
  scaling$sd[scaling$variable == "year"]
projection_start_year <- ndvi_config$projection_start_year

model_settings <- tibble(
  selected_model = selected_model,
  historical_start = as.character(min(historical_scaled$date)),
  historical_end = as.character(max(historical_scaled$date)),
  historical_rows = nrow(historical_scaled),
  historical_file_md5 = unname(tools::md5sum(historical_path)),
  mesh_cutoff_km = 20,
  year_effect = paste0("fixed at ", fixed_year),
  projection_start_year = projection_start_year,
  scenarios = paste(ndvi_config$future_scenarios, collapse = ";"),
  climate_models = paste(ndvi_config$chelsa_models, collapse = ";")
)
write_csv(model_settings, file.path(dir_final_model, "model_settings.csv"))

# --- Future projections -----------------------------------------------------

prepare_future_input <- function(model_name, scenario) {
  feature_path <- file.path(
    ndvi_future_features_dir(scenario),
    paste0("all_future_", model_name, ".feather")
  )
  if (!file.exists(feature_path)) {
    stop(
      "Missing future feature table: ", feature_path,
      "\nRun 01c_prepare_historical_feature_tables.R first."
    )
  }

  future <- read_feather(feature_path) |>
    add_future_covariates() |>
    apply_scaling(scaling) |>
    filter(year >= projection_start_year)
  future$year_sc <- fixed_year_sc
  future
}

prediction_blocks <- c(
  list(`2019-2040` = projection_start_year:2040),
  ndvi_config$analysis_windows
)
expected_projection_rows <- n_distinct(historical$x, historical$y) *
  12L * (max(unlist(prediction_blocks)) - projection_start_year + 1L)

read_valid_projection <- function(path) {
  if (!file.exists(path)) return(NULL)

  tryCatch({
    projection <- read_feather(path) |>
      mutate(date = as.Date(date))
    required <- c("x", "y", "date", "year", "month", "ndvi_pred")
    valid <- all(required %in% names(projection)) &&
      nrow(projection) == expected_projection_rows &&
      n_distinct(projection$x, projection$y, projection$date) ==
        expected_projection_rows &&
      min(projection$year) == projection_start_year &&
      max(projection$year) == max(unlist(prediction_blocks)) &&
      all(is.finite(projection$ndvi_pred)) &&
      all(projection$ndvi_pred > 0 & projection$ndvi_pred < 1)

    if (valid) projection else NULL
  }, error = function(e) NULL)
}

projection_index_row <- function(projection, scenario, model_name,
                                 output_path, status) {
  tibble(
    scenario = scenario,
    climate_model = model_name,
    status = status,
    rows = nrow(projection),
    cells = n_distinct(projection$x, projection$y),
    start = as.character(min(projection$date)),
    end = as.character(max(projection$date)),
    ndvi_min = min(projection$ndvi_pred),
    ndvi_max = max(projection$ndvi_pred),
    output = output_path
  )
}

project_gcm <- function(model_name) {
  map_dfr(ndvi_config$future_scenarios, function(scenario) {
    output_dir <- file.path(
      ndvi_scenario_outputs_dir(scenario, create = TRUE),
      "projections",
      model_name
    )
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_path <- file.path(output_dir, "monthly_ndvi.feather")

    existing <- read_valid_projection(output_path)
    if (!is.null(existing)) {
      cat("Reusing", scenario, "|", model_name, "\n")
      return(projection_index_row(
        existing,
        scenario,
        model_name,
        output_path,
        status = "reused"
      ))
    }

    cat("Projecting", scenario, "|", model_name, "\n")
    future <- prepare_future_input(model_name, scenario) |>
      mutate(.row_id = row_number())

    # The selected model has a static spatial field and no spatiotemporal
    # random effect, so predictions are row-independent. Predicting in the
    # four reporting blocks gives exactly the same estimates as one 2019-2100
    # call while avoiding a very large projection matrix in memory.
    prediction <- map_dfr(names(prediction_blocks), function(block_name) {
      block <- filter(future, year %in% prediction_blocks[[block_name]])
      cat("  ", block_name, "\n")
      tibble(
        .row_id = block$.row_id,
        ndvi_pred = predict(fit, newdata = block, type = "response")$est
      )
    }) |>
      arrange(.row_id)

    if (any(!is.finite(prediction$ndvi_pred)) ||
        any(prediction$ndvi_pred <= 0 | prediction$ndvi_pred >= 1)) {
      stop("Invalid NDVI prediction for ", scenario, " | ", model_name)
    }

    projection_output <- future |>
      transmute(
        x,
        y,
        date,
        year,
        month = month_num,
        ndvi_pred = prediction$ndvi_pred
      )
    write_feather(projection_output, output_path)

    projection_index_row(
      projection_output,
      scenario,
      model_name,
      output_path,
      status = "written"
    )
  })
}

gcm_models <- ndvi_config$chelsa_models
# sdmTMB's large spatial prediction matrix is memory-bound on this data set.
# Concurrent workers increased wall-clock time substantially, so the stable
# and reproducible publication default is one GCM at a time.
projection_parts <- lapply(gcm_models, project_gcm)

projection_index <- bind_rows(projection_parts) |>
  arrange(scenario, climate_model)
write_csv(
  projection_index,
  file.path(dir_scenarios, "projection_index.csv")
)

cat("\nCompleted", nrow(projection_index), "SSP-GCM projections.\n")
cat("Projection index:",
    file.path(dir_scenarios, "projection_index.csv"), "\n")
