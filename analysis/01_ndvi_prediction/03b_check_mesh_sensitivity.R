#!/usr/bin/env Rscript
#' ####################################################################### #
#' 03b - Check Spatial-Mesh Sensitivity
#'
#' Refit only the selected NDVI formula at 10, 20, and 30 km mesh cutoffs.
#' Compare temporal-holdout skill, spatial residual structure, estimated
#' spatial range, and fixed-response predictions. The 20 km holdout fit from
#' script 03 is reused; alternative fits are not retained unless a later
#' scientific decision makes one of them the main model.
#'
#' INPUTS:  outputs/selection/sdmtmb_best_model_observed.txt
#'          outputs/models/sdmtmb_winner_t2_observed_temporal_fit.rds
#'          derived/gimms_support/all_historical_observed.feather
#' OUTPUTS: outputs/sensitivity/mesh/mesh_comparison.csv
#'          outputs/sensitivity/mesh/mesh_response_comparison.csv
#' ####################################################################### #

pacman::p_load(
  "arrow",
  "tidyverse",
  "lubridate",
  "sf",
  "slider",
  "sdmTMB",
  "spdep"
)

source(here::here("analysis", "01_ndvi_prediction", "_config.R"))
source(here::here(
  "analysis", "01_ndvi_prediction", "R", "model_covariates.R"
))

selected_model <- readLines(
  file.path(
    ndvi_config$outputs_dir,
    "selection",
    "sdmtmb_best_model_observed.txt"
  ),
  warn = FALSE
)[1]
if (!identical(selected_model, "sdmtmb_winner_t2")) {
  stop("Mesh sensitivity is specified for sdmtmb_winner_t2; found ",
       selected_model)
}

historical <- read_feather(
  file.path(
    ndvi_config$gimms_support_dir,
    "all_historical_observed.feather"
  )
) |>
  ndvi_add_candidate_covariates()
training <- filter(
  historical,
  year < ndvi_config$temporal_holdout_start_year
)
holdout <- filter(
  historical,
  year >= ndvi_config$temporal_holdout_start_year
)
scaled <- ndvi_standardize_candidate_covariates(training, holdout)
training <- scaled$training
holdout <- scaled$holdout

output_dir <- file.path(
  ndvi_config$outputs_dir,
  "sensitivity",
  "mesh"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

reference_row <- training[1, , drop = FALSE]
reference_row$x <- median(training$x)
reference_row$y <- median(training$y)
reference_row$x_utm_km <- median(training$x_utm_km)
reference_row$y_utm_km <- median(training$y_utm_km)
reference_row$pr_lag1_sc <- median(training$pr_lag1_sc)
reference_row$pr_roll3_sc <- median(training$pr_roll3_sc)
reference_row$tasmax_roll3_sc <- median(training$tasmax_roll3_sc)
reference_row$ele_km_sc <- median(training$ele_km_sc)
reference_row$year_sc <- median(training$year_sc)

rain_scaled <- seq(
  quantile(training$pr_roll3_sc, 0.05),
  quantile(training$pr_roll3_sc, 0.95),
  length.out = 20
)
temperature_scaled <- seq(
  quantile(training$tasmax_roll3_sc, 0.05),
  quantile(training$tasmax_roll3_sc, 0.95),
  length.out = 20
)
response_index <- tidyr::crossing(
  rain_id = seq_along(rain_scaled),
  temperature_id = seq_along(temperature_scaled),
  month_num = 1:12
)
response_grid <- reference_row[rep(1, nrow(response_index)), , drop = FALSE]
response_grid$month_num <- response_index$month_num
response_grid$month <- factor(response_grid$month_num, levels = 1:12)
response_grid$pr_roll3_sc <- rain_scaled[response_index$rain_id]
response_grid$tasmax_roll3_sc <- temperature_scaled[
  response_index$temperature_id
]

rain_display <- rain_scaled * scaled$spreads[["pr_roll3"]] +
  scaled$centres[["pr_roll3"]]
temperature_display <- temperature_scaled *
  scaled$spreads[["tasmax_roll3"]] + scaled$centres[["tasmax_roll3"]]

mesh_cutoffs <- c(10, 20, 30)
comparison_parts <- list()
holdout_predictions <- list()
response_predictions <- list()

for (mesh_cutoff in mesh_cutoffs) {
  cat("\nEvaluating mesh cutoff:", mesh_cutoff, "km\n")
  mesh <- make_mesh(
    training,
    xy_cols = c("x_utm_km", "y_utm_km"),
    cutoff = mesh_cutoff
  )
  mesh_vertices <- mesh$mesh$n

  elapsed <- NA_real_
  if (mesh_cutoff == 20) {
    fit_path <- file.path(
      ndvi_config$outputs_dir,
      "models",
      paste0(selected_model, "_observed_temporal_fit.rds")
    )
    fit <- readRDS(fit_path)
    fit_source <- "reused selected fit"
  } else {
    set.seed(20260507 + mesh_cutoff)
    fit_time <- system.time({
      fit <- sdmTMB(
        ndvi_beta ~
          s(pr_lag1_sc, k = 5) +
          s(ele_km_sc, k = 4) +
          t2(pr_roll3_sc, tasmax_roll3_sc, k = c(5, 5)) +
          s(year_sc, k = 5) +
          month,
        data = training,
        mesh = mesh,
        family = Beta(link = "logit"),
        spatial = "on",
        spatiotemporal = "off",
        reml = TRUE,
        silent = TRUE
      )
    })
    elapsed <- unname(fit_time[["elapsed"]]) / 60
    fit_source <- "sensitivity refit"
  }

  holdout_prediction <- predict(
    fit,
    newdata = holdout,
    type = "response"
  )$est
  holdout_predictions[[as.character(mesh_cutoff)]] <- holdout_prediction

  monthly_prediction <- holdout |>
    mutate(predicted = holdout_prediction) |>
    group_by(date) |>
    summarise(
      observed = mean(ndvi),
      predicted = mean(predicted),
      .groups = "drop"
    )

  residual_cells <- holdout |>
    mutate(residual = ndvi - holdout_prediction) |>
    group_by(cell_id, x_utm_km, y_utm_km) |>
    summarise(mean_residual = mean(residual), .groups = "drop")
  neighbours <- residual_cells |>
    select(x_utm_km, y_utm_km) |>
    as.matrix() |>
    knearneigh(k = 4) |>
    knn2nb()
  weights <- nb2listw(neighbours, style = "W")
  moran_result <- moran.test(residual_cells$mean_residual, weights)

  response_prediction <- predict(
    fit,
    newdata = response_grid,
    type = "response",
    re_form = NA
  )$est
  response_predictions[[as.character(mesh_cutoff)]] <- response_prediction

  model_parameters <- fit$model$par
  kappa <- exp(unname(model_parameters[["ln_kappa"]]))
  tau <- exp(unname(model_parameters[["ln_tau_O"]]))

  comparison_parts[[as.character(mesh_cutoff)]] <- tibble(
    mesh_cutoff_km = mesh_cutoff,
    mesh_vertices = mesh_vertices,
    fit_source = fit_source,
    fit_minutes = elapsed,
    positive_definite_hessian = isTRUE(fit$pos_def_hessian),
    max_absolute_gradient = max(abs(fit$gradients)),
    spatial_range_km = sqrt(8) / kappa,
    spatial_sd = 1 / (sqrt(4 * pi) * tau * kappa),
    holdout_rmse = sqrt(mean((holdout$ndvi - holdout_prediction)^2)),
    holdout_mae = mean(abs(holdout$ndvi - holdout_prediction)),
    holdout_r_squared = cor(holdout$ndvi, holdout_prediction)^2,
    monthly_mean_rmse = sqrt(mean(
      (monthly_prediction$observed - monthly_prediction$predicted)^2
    )),
    monthly_mean_correlation = cor(
      monthly_prediction$observed,
      monthly_prediction$predicted
    ),
    holdout_residual_moran_i = unname(moran_result$estimate[[1]]),
    holdout_residual_moran_p = moran_result$p.value
  )

  write_csv(
    bind_rows(comparison_parts),
    file.path(output_dir, "mesh_comparison_partial.csv")
  )
  rm(fit)
  gc()
}

reference_holdout <- holdout_predictions[["20"]]
reference_response <- response_predictions[["20"]]
mesh_comparison <- bind_rows(comparison_parts) |>
  rowwise() |>
  mutate(
    holdout_rmse_vs_20km_predictions = sqrt(mean(
      (holdout_predictions[[as.character(mesh_cutoff_km)]] -
        reference_holdout)^2
    )),
    response_rmse_vs_20km = sqrt(mean(
      (response_predictions[[as.character(mesh_cutoff_km)]] -
        reference_response)^2
    )),
    response_max_abs_difference_vs_20km = max(abs(
      response_predictions[[as.character(mesh_cutoff_km)]] -
        reference_response
    ))
  ) |>
  ungroup()

response_comparison <- map_dfr(mesh_cutoffs, function(mesh_cutoff) {
  tibble(
    mesh_cutoff_km = mesh_cutoff,
    rain_id = response_index$rain_id,
    temperature_id = response_index$temperature_id,
    rainfall_mm = rain_display[response_index$rain_id],
    temperature_c = temperature_display[response_index$temperature_id],
    month = response_index$month_num,
    predicted_ndvi = response_predictions[[as.character(mesh_cutoff)]]
  )
}) |>
  group_by(
    mesh_cutoff_km,
    rain_id,
    temperature_id,
    rainfall_mm,
    temperature_c
  ) |>
  summarise(predicted_ndvi = mean(predicted_ndvi), .groups = "drop")

write_csv(
  mesh_comparison,
  file.path(output_dir, "mesh_comparison.csv")
)
write_csv(
  response_comparison,
  file.path(output_dir, "mesh_response_comparison.csv")
)
unlink(file.path(output_dir, "mesh_comparison_partial.csv"))

cat("\nMesh sensitivity complete.\n")
print(mesh_comparison)
