#' ####################################################################### #
#' 03 — Compare sdmTMB Candidate Models
#'
#' Fit a set of theory-driven sdmTMB models on GIMMS-support observed
#' history using temporal holdout validation. Treat RMSE differences below
#' 1e-4 as ties. Within that tied set, prefer the tensor-product interaction,
#' which lets rainfall and temperature use separate smoothness scales.
#'
#' Candidates vary in: rainfall term (rolling mean vs memory decomposition),
#' temperature metric (Tmax, Tmin, Tmean), year trend shape, temporal
#' dependence, and rain x heat interaction form. All share: Beta family,
#' spatial random field, month fixed effects.
#'
#' INPUTS:  derived/gimms_support/all_historical_observed.feather
#' OUTPUTS: outputs/selection/  - comparison table + selected model name
#'          outputs/models/     - selected temporal-holdout fit
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
source(here::here(
  "analysis", "01_ndvi_prediction", "R", "model_covariates.R"
))

# --- Helper functions --------------------------------------------------------

safe_cor_sq <- function(x, y) {
  if (length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(NA_real_)
  }
  cor(x, y)^2
}

compute_metrics <- function(obs, pred) {
  tibble(
    rmse = sqrt(mean((obs - pred)^2)),
    mae = mean(abs(obs - pred)),
    r2 = safe_cor_sq(obs, pred)
  )
}

compute_monthly_mean_skill <- function(df) {
  monthly_df <- df |>
    group_by(date) |>
    summarise(
      mean_obs = mean(ndvi),
      mean_pred = mean(pred),
      .groups = "drop"
    )

  tibble(
    monthly_mean_rmse = sqrt(mean((monthly_df$mean_obs - monthly_df$mean_pred)^2)),
    monthly_mean_mae = mean(abs(monthly_df$mean_obs - monthly_df$mean_pred)),
    monthly_mean_r = if (length(unique(monthly_df$mean_pred)) < 2) NA_real_ else cor(monthly_df$mean_obs, monthly_df$mean_pred)
  )
}

# --- Model fitting -----------------------------------------------------------

fit_candidate_model <- function(train_df, model_id, mesh_cutoff_km = 20, extra_time = NULL) {
  mesh <- make_mesh(train_df, xy_cols = c("x_utm_km", "y_utm_km"), cutoff = mesh_cutoff_km)

  if (model_id == "sdmtmb_baseline_season_year") {
    fit <- sdmTMB(
      ndvi_beta ~
        s(ele_km_sc, k = 4) +
        year_sc +
        month,
      data = train_df,
      mesh = mesh,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "off",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_linear_core_tmax") {
    fit <- sdmTMB(
      ndvi_beta ~ pr_roll3_sc + tasmax_roll3_sc + ele_km_sc + month,
      data = train_df,
      mesh = mesh,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "off",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_smooth_core_tmax") {
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        month,
      data = train_df,
      mesh = mesh,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "off",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_smooth_core_tmax_year") {
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        year_sc +
        month,
      data = train_df,
      mesh = mesh,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "off",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_smooth_core_tmax_syear_int_plus_prlag1") {
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(pr_lag1_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        s(pr_tasmax_inter, k = 5) +
        s(year_sc, k = 5) +
        month,
      data = train_df,
      mesh = mesh,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "off",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_smooth_core_tmax_year_ar1") {
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        year_sc +
        month,
      data = train_df,
      mesh = mesh,
      time = "time_id",
      extra_time = extra_time,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "ar1",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_smooth_core_tmax_year_int_plus_prlag1_ar1") {
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(pr_lag1_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        s(pr_tasmax_inter, k = 5) +
        year_sc +
        month,
      data = train_df,
      mesh = mesh,
      time = "time_id",
      extra_time = extra_time,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "ar1",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_smooth_core_tmax_year_int") {
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        s(pr_tasmax_inter, k = 5) +
        year_sc +
        month,
      data = train_df,
      mesh = mesh,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "off",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_smooth_core_tmax_year_surface_plus_prlag1") {
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(pr_lag1_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        s(pr_roll3_sc, tasmax_roll3_sc, k = 15) +
        year_sc +
        month,
      data = train_df,
      mesh = mesh,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "off",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_smooth_core_tmax_year_int_plus_prlag1") {
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(pr_lag1_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        s(pr_tasmax_inter, k = 5) +
        year_sc +
        month,
      data = train_df,
      mesh = mesh,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "off",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_smooth_core_tmax_year_int_plus_tmax_now") {
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(tasmax_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        s(pr_tasmax_inter, k = 5) +
        year_sc +
        month,
      data = train_df,
      mesh = mesh,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "off",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_smooth_core_tmin_year") {
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(tasmin_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        year_sc +
        month,
      data = train_df,
      mesh = mesh,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "off",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_smooth_core_tmean_year") {
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(tasmean_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        year_sc +
        month,
      data = train_df,
      mesh = mesh,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "off",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_smooth_memory_tmax_year_int") {
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll6_sc, k = 5) +
        s(pr_recent_anom_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        s(prrecent_tasmax_inter, k = 5) +
        year_sc +
        month,
      data = train_df,
      mesh = mesh,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "off",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_winner_iid") {
    # Winner's formula + IID spatiotemporal random field per time step.
    # Tests whether independent yearly spatial drift absorbs the residual
    # lag-1 ACF (~0.42 per cell) seen in the diagnostics figure.
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(pr_lag1_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        s(pr_tasmax_inter, k = 5) +
        s(year_sc, k = 5) +
        month,
      data = train_df,
      mesh = mesh,
      time = "time_id",
      extra_time = extra_time,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "iid",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_winner_rw1") {
    # Winner's formula + Random-Walk (RW1) spatiotemporal field. Smooths
    # year-to-year fields rather than reverting to the mean (AR1) — usually
    # safer for forecasting, per sdmTMB forecasting vignette.
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(pr_lag1_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        s(pr_tasmax_inter, k = 5) +
        s(year_sc, k = 5) +
        month,
      data = train_df,
      mesh = mesh,
      time = "time_id",
      extra_time = extra_time,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "rw",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_winner_ar1_fair") {
    # Winner's formula + AR1 (a "fair" AR1 — the original AR1 candidate had
    # a simpler formula; this one matches the winner so the comparison is
    # apples-to-apples).
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(pr_lag1_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        s(pr_tasmax_inter, k = 5) +
        s(year_sc, k = 5) +
        month,
      data = train_df,
      mesh = mesh,
      time = "time_id",
      extra_time = extra_time,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "ar1",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_winner_no_year_iid") {
    # Same as winner_iid but with year_sc dropped — tests the hypothesis
    # that the temporal random field overfits because it competes with
    # year for the long-term trend.
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(pr_lag1_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        s(pr_tasmax_inter, k = 5) +
        month,
      data = train_df,
      mesh = mesh,
      time = "time_id",
      extra_time = extra_time,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "iid",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_winner_no_year_rw1") {
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(pr_lag1_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        s(pr_tasmax_inter, k = 5) +
        month,
      data = train_df,
      mesh = mesh,
      time = "time_id",
      extra_time = extra_time,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "rw",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_winner_no_year_ar1") {
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll3_sc, k = 5) +
        s(pr_lag1_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        s(pr_tasmax_inter, k = 5) +
        month,
      data = train_df,
      mesh = mesh,
      time = "time_id",
      extra_time = extra_time,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "ar1",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_winner_t2") {
    # Winner's formula but the rain x heat interaction becomes a tensor
    # product surface t2(pr, tmax) instead of the smooth-of-product
    # s(pr * tmax). t2 is mgcv's anisotropic-stable alternative to s(x, y).
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_lag1_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        t2(pr_roll3_sc, tasmax_roll3_sc, k = c(5, 5)) +
        s(year_sc, k = 5) +
        month,
      data = train_df,
      mesh = mesh,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "off",
      reml = TRUE,
      silent = TRUE
    )
  } else if (model_id == "sdmtmb_smooth_memory_tmax_year_int_ar1") {
    fit <- sdmTMB(
      ndvi_beta ~
        s(pr_roll6_sc, k = 5) +
        s(pr_recent_anom_sc, k = 5) +
        s(tasmax_roll3_sc, k = 5) +
        s(ele_km_sc, k = 4) +
        s(prrecent_tasmax_inter, k = 5) +
        year_sc +
        month,
      data = train_df,
      mesh = mesh,
      time = "time_id",
      extra_time = extra_time,
      family = Beta(link = "logit"),
      spatial = "on",
      spatiotemporal = "ar1",
      reml = TRUE,
      silent = TRUE
    )
  } else {
    stop("Unknown model_id: ", model_id)
  }

  list(fit = fit, mesh = mesh)
}

# --- Run ---------------------------------------------------------------------

Dir.Data <- ndvi_config$gimms_support_dir
Dir.OutRoot <- ndvi_config$outputs_dir
Dir.Selection <- file.path(Dir.OutRoot, "selection")
Dir.Models <- file.path(Dir.OutRoot, "models")

historical_dataset_id <- "observed"
temporal_holdout_start_year <- ndvi_config$temporal_holdout_start_year
mesh_cutoff_km <- 20

walk(c(Dir.Selection, Dir.Models), dir.create, recursive = TRUE, showWarnings = FALSE)

candidate_models <- c(
  "sdmtmb_baseline_season_year",
  "sdmtmb_linear_core_tmax",
  "sdmtmb_smooth_core_tmax",
  "sdmtmb_smooth_core_tmax_year",
  "sdmtmb_smooth_core_tmax_syear_int_plus_prlag1",
  "sdmtmb_smooth_core_tmax_year_int",
  "sdmtmb_smooth_core_tmax_year_int_plus_prlag1",
  "sdmtmb_smooth_core_tmax_year_int_plus_tmax_now",
  "sdmtmb_smooth_core_tmin_year",
  "sdmtmb_smooth_core_tmean_year",
  "sdmtmb_smooth_memory_tmax_year_int",
  # Locked 2026-05 winner: tensor rain x heat surface plus frozen-year
  # projection in script 05. AR1/RW/IID and no-year variants are no longer
  # in the default comparison because they fit training data too closely and
  # degrade temporal-holdout prediction.
  "sdmtmb_winner_t2"
)
preferred_tied_model <- "sdmtmb_winner_t2"
if (!preferred_tied_model %in% candidate_models) {
  stop("Preferred tied-model structure is absent from candidate_models.")
}

cat("Loading GIMMS-support historical data for:", historical_dataset_id, "\n")
hist_df <- read_feather(file.path(Dir.Data, paste0("all_historical_", historical_dataset_id, ".feather"))) |>
  ndvi_add_candidate_covariates(include_time_id = TRUE)

if (!temporal_holdout_start_year %in% hist_df$year) {
  stop(
    "TEMPORAL_HOLDOUT_START_YEAR=", temporal_holdout_start_year,
    " is not present in the historical dataset."
  )
}

holdout_start_year <- temporal_holdout_start_year
train_df <- hist_df |>
  filter(year < holdout_start_year)
test_df <- hist_df |>
  filter(year >= holdout_start_year)

cat("Training years:", min(train_df$year), "-", max(train_df$year), "\n")
cat("Holdout years:", min(test_df$year), "-", max(test_df$year), "\n")
cat("Training rows:", nrow(train_df), "| Holdout rows:", nrow(test_df), "\n")

scaled <- ndvi_standardize_candidate_covariates(train_df, test_df)
train_df <- scaled$training
test_df <- scaled$holdout
extra_time_ids <- sort(unique(test_df$time_id[!test_df$time_id %in% train_df$time_id]))

comparison_df <- map_dfr(candidate_models, function(model_id) {
  cat("\n=============================================================\n")
  cat("Fitting sdmTMB candidate:", model_id, "\n")
  cat("=============================================================\n")

  result <- tryCatch({
    fitted_obj <- fit_candidate_model(
      train_df,
      model_id,
      mesh_cutoff_km = mesh_cutoff_km,
      extra_time = extra_time_ids
    )
    fit <- fitted_obj$fit

    train_pred <- predict(fit, newdata = train_df, type = "response")$est
    test_pred <- predict(fit, newdata = test_df, type = "response")$est

    train_metrics <- compute_metrics(train_df$ndvi, train_pred)
    test_metrics <- compute_metrics(test_df$ndvi, test_pred)
    monthly_skill <- compute_monthly_mean_skill(test_df |> mutate(pred = test_pred))

    saveRDS(
      fit,
      file.path(Dir.Models, paste0(model_id, "_", historical_dataset_id, "_temporal_fit.rds"))
    )

    tibble(
      historical_dataset_id = historical_dataset_id,
      holdout_start_year = holdout_start_year,
      holdout_end_year = max(test_df$year),
      model_id = model_id,
      train_rmse = train_metrics$rmse,
      train_mae = train_metrics$mae,
      train_r2 = train_metrics$r2,
      holdout_rmse = test_metrics$rmse,
      holdout_mae = test_metrics$mae,
      holdout_r2 = test_metrics$r2,
      holdout_monthly_mean_rmse = monthly_skill$monthly_mean_rmse,
      holdout_monthly_mean_mae = monthly_skill$monthly_mean_mae,
      holdout_monthly_mean_r = monthly_skill$monthly_mean_r
    )
  }, error = function(e) {
    cat("  FAILED:", conditionMessage(e), "\n")
    tibble(
      historical_dataset_id = historical_dataset_id,
      holdout_start_year = holdout_start_year,
      holdout_end_year = max(test_df$year),
      model_id = model_id,
      train_rmse = NA_real_, train_mae = NA_real_, train_r2 = NA_real_,
      holdout_rmse = NA_real_, holdout_mae = NA_real_, holdout_r2 = NA_real_,
      holdout_monthly_mean_rmse = NA_real_, holdout_monthly_mean_mae = NA_real_,
      holdout_monthly_mean_r = NA_real_
    )
  })

  result
})

rmse_tie_tolerance <- 1e-4
if (all(is.na(comparison_df$holdout_rmse))) {
  stop("All NDVI candidate models failed; no model can be selected.")
}
best_rmse <- min(comparison_df$holdout_rmse, na.rm = TRUE)

comparison_df <- comparison_df |>
  mutate(
    within_rmse_tie = !is.na(holdout_rmse) &
      holdout_rmse <= best_rmse + rmse_tie_tolerance,
    preferred_tied_structure = within_rmse_tie &
      model_id == preferred_tied_model
  ) |>
  arrange(
    desc(within_rmse_tie),
    desc(preferred_tied_structure),
    holdout_rmse,
    desc(holdout_r2),
    model_id
  ) |>
  mutate(selected = row_number() == 1L) |>
  select(-preferred_tied_structure)

n_failed <- sum(is.na(comparison_df$holdout_rmse))
if (n_failed > 0) cat("\n", n_failed, "candidate(s) failed to converge.\n")

best_model <- comparison_df |>
  filter(.data$selected) |>
  pull(.data$model_id)
if (length(best_model) != 1L) {
  stop("Candidate ranking did not identify exactly one selected model.")
}
write_csv(
  comparison_df,
  file.path(
    Dir.Selection,
    paste0("sdmtmb_candidate_comparison_", historical_dataset_id, ".csv")
  )
)
writeLines(
  best_model,
  file.path(
    Dir.Selection,
    paste0("sdmtmb_best_model_", historical_dataset_id, ".txt")
  )
)

# Candidate fits are expensive but only the selected holdout fit is needed for
# diagnostics. Remove the other large model objects after recording the full
# comparison table.
drop_models <- setdiff(candidate_models, best_model)
walk(drop_models, function(model_id) {
  unlink(file.path(
    Dir.Models,
    paste0(model_id, "_", historical_dataset_id, "_temporal_fit.rds")
  ))
})

cat("\nBest provisional model:", best_model, "\n")
cat("Saved comparison outputs to:", Dir.OutRoot, "\n")
