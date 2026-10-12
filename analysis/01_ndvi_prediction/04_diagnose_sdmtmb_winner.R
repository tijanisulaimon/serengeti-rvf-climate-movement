#!/usr/bin/env Rscript
#' ####################################################################### #
#' 04 - Validate the Selected sdmTMB Model
#'
#' Evaluate the selected model against the withheld 2015-2018 observations.
#' The diagnostic figure focuses on temporal prediction and residual structure;
#' a separate response figure shows each selected model component once.
#'
#' INPUTS:  outputs/selection/sdmtmb_best_model_observed.txt
#'          outputs/models/{selected_model}_observed_temporal_fit.rds
#'          derived/gimms_support/all_historical_observed.feather
#' OUTPUTS: outputs/diagnostics/selected_model_metrics.csv
#'          outputs/diagnostics/selected_model_effects.csv
#'          outputs/diagnostics/selected_model_sanity.txt
#'          outputs/plots/diagnostics/selected_model_validation.png
#'          outputs/plots/diagnostics/selected_model_responses.png
#' ####################################################################### #

pacman::p_load(
  "arrow",
  "tidyverse",
  "lubridate",
  "sf",
  "slider",
  "sdmTMB",
  "spdep",
  "patchwork"
)

source(here::here("analysis", "01_ndvi_prediction", "_config.R"))
source(here::here(
  "analysis", "01_ndvi_prediction", "R", "model_covariates.R"
))

# --- Recreate the model-selection data --------------------------------------

dir_data <- ndvi_config$gimms_support_dir
dir_outputs <- ndvi_config$outputs_dir
dir_models <- file.path(dir_outputs, "models")
dir_selection <- file.path(dir_outputs, "selection")
dir_diagnostics <- file.path(dir_outputs, "diagnostics")
dir_figures <- file.path(dir_outputs, "plots", "diagnostics")
dir.create(dir_diagnostics, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

selected_model <- readLines(
  file.path(dir_selection, "sdmtmb_best_model_observed.txt"),
  warn = FALSE
)[1]
fit_path <- file.path(
  dir_models,
  paste0(selected_model, "_observed_temporal_fit.rds")
)
if (!file.exists(fit_path)) stop("Missing temporal-holdout fit: ", fit_path)

historical <- read_feather(
  file.path(dir_data, "all_historical_observed.feather")
) |>
  ndvi_add_candidate_covariates()

holdout_start <- ndvi_config$temporal_holdout_start_year
training <- filter(historical, year < holdout_start)
holdout <- filter(historical, year >= holdout_start)
scaled <- ndvi_standardize_candidate_covariates(training, holdout)
training <- scaled$training
holdout <- scaled$holdout

fit <- readRDS(fit_path)
training$predicted <- predict(fit, newdata = training, type = "response")$est
holdout$predicted <- predict(fit, newdata = holdout, type = "response")$est
training$residual <- training$ndvi - training$predicted
holdout$residual <- holdout$ndvi - holdout$predicted

# --- Validation metrics and residual summaries ------------------------------

metric_row <- function(df, split) {
  tibble(
    split = split,
    n = nrow(df),
    rmse = sqrt(mean((df$ndvi - df$predicted)^2)),
    mae = mean(abs(df$ndvi - df$predicted)),
    r_squared = cor(df$ndvi, df$predicted)^2
  )
}

holdout_monthly <- holdout |>
  group_by(date) |>
  summarise(
    observed = mean(ndvi),
    predicted = mean(predicted),
    mean_residual = mean(residual),
    .groups = "drop"
  )

validation_metrics <- bind_rows(
  metric_row(training, "training"),
  metric_row(holdout, "temporal holdout")
)
monthly_metrics <- tibble(
  split = "temporal holdout, district monthly mean",
  n = nrow(holdout_monthly),
  rmse = sqrt(mean((holdout_monthly$observed -
    holdout_monthly$predicted)^2)),
  mae = mean(abs(holdout_monthly$observed -
    holdout_monthly$predicted)),
  r_squared = cor(
    holdout_monthly$observed,
    holdout_monthly$predicted
  )^2
)
validation_metrics <- bind_rows(validation_metrics, monthly_metrics)

holdout_cell_residual <- holdout |>
  group_by(cell_id, x, y, x_utm_km, y_utm_km) |>
  summarise(mean_residual = mean(residual), .groups = "drop")

neighbours <- holdout_cell_residual |>
  select(x_utm_km, y_utm_km) |>
  as.matrix() |>
  knearneigh(k = 4) |>
  knn2nb()
spatial_weights <- nb2listw(neighbours, style = "W")
moran_result <- moran.test(
  holdout_cell_residual$mean_residual,
  spatial_weights
)
moran_summary <- tibble(
  moran_i = unname(moran_result$estimate[[1]]),
  p_value = moran_result$p.value
)

acf_object <- acf(
  holdout_monthly$mean_residual,
  lag.max = 12,
  plot = FALSE
)
acf_summary <- tibble(
  lag_months = as.numeric(acf_object$lag)[-1],
  acf = as.numeric(acf_object$acf)[-1]
)
acf_bound <- 1.96 / sqrt(nrow(holdout_monthly))

pairwise_distance <- as.matrix(
  dist(holdout_cell_residual[, c("x_utm_km", "y_utm_km")])
)
pairwise_residual <- outer(
  holdout_cell_residual$mean_residual,
  holdout_cell_residual$mean_residual,
  function(a, b) (a - b)^2
)
upper_triangle <- upper.tri(pairwise_distance)
variogram <- tibble(
  distance_km = pairwise_distance[upper_triangle],
  semivariance = 0.5 * pairwise_residual[upper_triangle]
) |>
  mutate(
    distance_bin = cut(
      distance_km,
      breaks = seq(0, max(distance_km) + 5, by = 5),
      include.lowest = TRUE
    )
  ) |>
  group_by(distance_bin) |>
  summarise(
    distance_km = mean(distance_km),
    semivariance = mean(semivariance),
    pairs = n(),
    .groups = "drop"
  ) |>
  filter(pairs >= 30)

write_csv(
  validation_metrics,
  file.path(dir_diagnostics, "selected_model_metrics.csv")
)
write_csv(
  moran_summary,
  file.path(dir_diagnostics, "selected_model_moran.csv")
)
write_csv(
  acf_summary,
  file.path(dir_diagnostics, "selected_model_residual_acf.csv")
)

sanity_text <- capture.output(sanity(fit))
if (!length(sanity_text)) {
  sanity_text <- "sdmTMB::sanity() completed without capturable console text."
}
writeLines(
  sanity_text,
  file.path(dir_diagnostics, "selected_model_sanity.txt")
)

# --- Holdout-centred validation figure --------------------------------------

theme_diagnostic <- theme_minimal(base_size = 10) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

holdout_r2 <- validation_metrics$r_squared[
  validation_metrics$split == "temporal holdout"
]
holdout_rmse <- validation_metrics$rmse[
  validation_metrics$split == "temporal holdout"
]

observed_predicted_plot <- ggplot(
  holdout,
  aes(predicted, ndvi)
) +
  geom_point(alpha = 0.10, size = 0.45, color = "#0072B2") +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    color = "#B2182B"
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    hjust = -0.08,
    vjust = 1.3,
    label = sprintf("R2 = %.3f\nRMSE = %.3f", holdout_r2, holdout_rmse),
    size = 3.2
  ) +
  coord_equal() +
  theme_diagnostic +
  labs(
    title = "A  Point-level temporal holdout",
    x = "Predicted NDVI",
    y = "Observed NDVI"
  )

holdout_time_plot <- holdout_monthly |>
  select(date, observed, predicted) |>
  pivot_longer(
    c(observed, predicted),
    names_to = "series",
    values_to = "mean_ndvi"
  ) |>
  mutate(
    series = recode(
      series,
      observed = "Observed",
      predicted = "Predicted"
    )
  ) |>
  ggplot(aes(date, mean_ndvi, color = series)) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(
    values = c(Observed = "black", Predicted = "#0072B2"),
    name = NULL
  ) +
  theme_diagnostic +
  theme(legend.position = "top") +
  labs(
    title = "B  District monthly mean",
    x = NULL,
    y = "Mean NDVI"
  )

residual_limit <- max(abs(holdout_cell_residual$mean_residual))
residual_map_plot <- ggplot(
  holdout_cell_residual,
  aes(x, y, fill = mean_residual)
) +
  geom_tile() +
  coord_equal() +
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#1B7837",
    midpoint = 0,
    limits = c(-residual_limit, residual_limit),
    name = "Mean\nresidual"
  ) +
  theme_diagnostic +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "C  Spatial residual pattern", x = NULL, y = NULL)

variogram_plot <- ggplot(
  variogram,
  aes(distance_km, semivariance)
) +
  geom_point(
    aes(size = pairs),
    color = "#2166AC",
    alpha = 0.75,
    show.legend = FALSE
  ) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    span = 0.5,
    color = "#B2182B",
    linewidth = 0.8
  ) +
  theme_diagnostic +
  labs(
    title = "D  Holdout residual variogram",
    x = "Distance (km)",
    y = "Semivariance"
  )

acf_plot <- ggplot(acf_summary, aes(lag_months, acf)) +
  geom_hline(yintercept = 0, color = "grey45") +
  geom_hline(
    yintercept = c(-acf_bound, acf_bound),
    linetype = 2,
    color = "#B2182B",
    linewidth = 0.4
  ) +
  geom_segment(
    aes(xend = lag_months, yend = 0),
    color = "#2166AC",
    linewidth = 0.7
  ) +
  geom_point(color = "#2166AC", size = 1.6) +
  scale_x_continuous(breaks = 1:12) +
  theme_diagnostic +
  labs(
    title = "E  Monthly mean residual autocorrelation",
    x = "Lag (months)",
    y = "ACF"
  )

validation_figure <-
  (observed_predicted_plot | holdout_time_plot) /
  (residual_map_plot | variogram_plot) /
  acf_plot +
  plot_layout(heights = c(1, 1, 0.8))

ggsave(
  file.path(dir_figures, "selected_model_validation.png"),
  validation_figure,
  width = 11,
  height = 10.5,
  dpi = 300
)

# --- Conditional model responses -------------------------------------------

# Each one-dimensional curve averages population-level predictions over the
# 12 month effects while holding the other covariates at their training
# medians. Curves are descriptive model responses, not uncertainty intervals.
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

make_response_curve <- function(scaled_variable, raw_variable, label) {
  scaled_values <- seq(
    quantile(training[[scaled_variable]], 0.01),
    quantile(training[[scaled_variable]], 0.99),
    length.out = 80
  )
  grid_index <- tidyr::crossing(
    value_id = seq_along(scaled_values),
    month_num = 1:12
  )
  grid <- reference_row[rep(1, nrow(grid_index)), , drop = FALSE]
  grid$month_num <- grid_index$month_num
  grid$month <- factor(grid$month_num, levels = 1:12)
  grid[[scaled_variable]] <- scaled_values[grid_index$value_id]

  predictions <- predict(
    fit,
    newdata = grid,
    type = "response",
    re_form = NA
  )$est
  display_values <- scaled_values * scaled$spreads[[raw_variable]] +
    scaled$centres[[raw_variable]]

  tibble(
    response = label,
    value_id = grid_index$value_id,
    value = display_values[grid_index$value_id],
    predicted_ndvi = predictions
  ) |>
    group_by(response, value_id, value) |>
    summarise(predicted_ndvi = mean(predicted_ndvi), .groups = "drop")
}

response_curves <- bind_rows(
  make_response_curve("pr_lag1_sc", "pr_lag1", "Lagged rainfall (mm)"),
  make_response_curve("ele_km_sc", "ele_km", "Elevation (km)"),
  make_response_curve("year_sc", "year", "Year")
)

month_grid <- reference_row[rep(1, 12), , drop = FALSE]
month_grid$month_num <- 1:12
month_grid$month <- factor(1:12, levels = 1:12)
month_response <- tibble(
  month = 1:12,
  predicted_ndvi = predict(
    fit,
    newdata = month_grid,
    type = "response",
    re_form = NA
  )$est
)

rain_scaled <- seq(
  quantile(training$pr_roll3_sc, 0.05),
  quantile(training$pr_roll3_sc, 0.95),
  length.out = 40
)
temperature_scaled <- seq(
  quantile(training$tasmax_roll3_sc, 0.05),
  quantile(training$tasmax_roll3_sc, 0.95),
  length.out = 40
)
surface_index <- tidyr::crossing(
  rain_id = seq_along(rain_scaled),
  temperature_id = seq_along(temperature_scaled),
  month_num = 1:12
)
surface_grid <- reference_row[rep(1, nrow(surface_index)), , drop = FALSE]
surface_grid$month_num <- surface_index$month_num
surface_grid$month <- factor(surface_grid$month_num, levels = 1:12)
surface_grid$pr_roll3_sc <- rain_scaled[surface_index$rain_id]
surface_grid$tasmax_roll3_sc <- temperature_scaled[
  surface_index$temperature_id
]
surface_prediction <- predict(
  fit,
  newdata = surface_grid,
  type = "response",
  re_form = NA
)$est

rain_display <- rain_scaled * scaled$spreads[["pr_roll3"]] +
  scaled$centres[["pr_roll3"]]
temperature_display <- temperature_scaled *
  scaled$spreads[["tasmax_roll3"]] + scaled$centres[["tasmax_roll3"]]
response_surface <- tibble(
  rain_id = surface_index$rain_id,
  temperature_id = surface_index$temperature_id,
  rainfall_mm = rain_display[surface_index$rain_id],
  temperature_c = temperature_display[surface_index$temperature_id],
  predicted_ndvi = surface_prediction
) |>
  group_by(rain_id, temperature_id, rainfall_mm, temperature_c) |>
  summarise(predicted_ndvi = mean(predicted_ndvi), .groups = "drop")

effect_output <- bind_rows(
  response_curves |>
    transmute(
      component = response,
      x = value,
      y = NA_real_,
      predicted_ndvi
    ),
  month_response |>
    transmute(
      component = "Month",
      x = month,
      y = NA_real_,
      predicted_ndvi
    ),
  response_surface |>
    transmute(
      component = "Rainfall-temperature surface",
      x = rainfall_mm,
      y = temperature_c,
      predicted_ndvi
    )
)
write_csv(
  effect_output,
  file.path(dir_diagnostics, "selected_model_effects.csv")
)

response_support <- list(
  `Lagged rainfall (mm)` = training$pr_lag1,
  `Elevation (km)` = training$ele_km,
  Year = training$year
)

set.seed(103)
line_response_plots <- map(names(response_support), function(label) {
  response_data <- filter(response_curves, response == label)
  support_values <- sample(
    response_support[[label]],
    min(length(response_support[[label]]), 1500)
  )

  ggplot(response_data, aes(value, predicted_ndvi)) +
    geom_rug(
      data = tibble(value = support_values),
      aes(x = value),
      inherit.aes = FALSE,
      alpha = 0.06,
      color = "grey25"
    ) +
    geom_line(color = "#2166AC", linewidth = 0.9) +
    theme_diagnostic +
    labs(
      title = label,
      x = label,
      y = "Predicted NDVI (conditional)"
    )
})

month_plot <- ggplot(month_response, aes(month, predicted_ndvi)) +
  geom_line(color = "#2166AC", linewidth = 0.9) +
  geom_point(color = "#2166AC", size = 1.8) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme_diagnostic +
  labs(
    title = "Month effect",
    x = NULL,
    y = "Predicted NDVI (conditional)"
  )

set.seed(104)
surface_rainfall_range <- range(response_surface$rainfall_mm)
surface_temperature_range <- range(response_surface$temperature_c)
surface_support <- training |>
  transmute(
    rainfall_mm = pr_roll3,
    temperature_c = tasmax_roll3
  ) |>
  filter(
    between(
      rainfall_mm,
      surface_rainfall_range[[1]],
      surface_rainfall_range[[2]]
    ),
    between(
      temperature_c,
      surface_temperature_range[[1]],
      surface_temperature_range[[2]]
    )
  )
surface_support <- surface_support |>
  slice_sample(n = min(nrow(surface_support), 6000L))

surface_plot <- ggplot(
  response_surface,
  aes(rainfall_mm, temperature_c)
) +
  geom_raster(aes(fill = predicted_ndvi), interpolate = TRUE) +
  geom_contour(
    aes(z = predicted_ndvi),
    color = "white",
    alpha = 0.55,
    linewidth = 0.3
  ) +
  geom_density_2d(
    data = surface_support,
    aes(rainfall_mm, temperature_c),
    inherit.aes = FALSE,
    color = "black",
    alpha = 0.45,
    linewidth = 0.35,
    bins = 4
  ) +
  scale_fill_viridis_c(option = "C", name = "Predicted NDVI") +
  coord_cartesian(
    xlim = surface_rainfall_range,
    ylim = surface_temperature_range,
    expand = FALSE,
    clip = "on"
  ) +
  theme_diagnostic +
  theme(legend.position = "right") +
  labs(
    title = "Rainfall-temperature interaction",
    x = "Three-month mean rainfall (mm)",
    y = "Three-month mean maximum temperature (deg C)"
  )

response_figure <-
  wrap_plots(c(line_response_plots, list(month_plot)), ncol = 2) /
  surface_plot +
  plot_layout(heights = c(1, 1.35))

ggsave(
  file.path(dir_figures, "selected_model_responses.png"),
  response_figure,
  width = 11,
  height = 10,
  dpi = 300
)

cat("Selected-model validation complete.\n")
cat("Metrics:", file.path(dir_diagnostics, "selected_model_metrics.csv"), "\n")
