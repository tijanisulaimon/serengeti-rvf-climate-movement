# Shared covariate construction for NDVI candidate fitting and diagnostics.
# Keeping these transformations in one place prevents the selected-model
# diagnostics and mesh sensitivity from silently using a different sample.

ndvi_add_candidate_covariates <- function(df, include_time_id = FALSE) {
  coordinates_sf <- df |>
    dplyr::distinct(.data$x, .data$y) |>
    sf::st_as_sf(
      coords = c("x", "y"),
      crs = 4326,
      remove = FALSE
    ) |>
    sf::st_transform(32736)
  coordinates_km <- sf::st_coordinates(coordinates_sf) / 1000
  coordinates <- coordinates_sf |>
    sf::st_drop_geometry() |>
    dplyr::mutate(
      x_utm_km = coordinates_km[, 1],
      y_utm_km = coordinates_km[, 2]
    )

  covariates <- df |>
    dplyr::mutate(
      date = as.Date(.data$date),
      year = lubridate::year(.data$date),
      month_num = lubridate::month(.data$date),
      month = factor(.data$month_num, levels = 1:12),
      ele_km = .data$ele / 1000,
      dtr = .data$tasmax - .data$tasmin,
      tasmean = (.data$tasmax + .data$tasmin) / 2
    ) |>
    dplyr::left_join(coordinates, by = c("x", "y")) |>
    dplyr::arrange(.data$x, .data$y, .data$date) |>
    dplyr::group_by(.data$x, .data$y) |>
    dplyr::mutate(
      pr_lag1 = dplyr::lag(.data$pr, 1),
      pr_roll3 = slider::slide_dbl(
        .data$pr,
        mean,
        .before = 2,
        .complete = TRUE
      ),
      pr_roll6 = slider::slide_dbl(
        .data$pr,
        mean,
        .before = 5,
        .complete = TRUE
      ),
      tasmax_roll3 = slider::slide_dbl(
        .data$tasmax,
        mean,
        .before = 2,
        .complete = TRUE
      ),
      tasmin_roll3 = slider::slide_dbl(
        .data$tasmin,
        mean,
        .before = 2,
        .complete = TRUE
      ),
      tasmean_roll3 = slider::slide_dbl(
        .data$tasmean,
        mean,
        .before = 2,
        .complete = TRUE
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      pr_recent_anom = .data$pr_roll3 - .data$pr_roll6,
      cell_id = as.character(interaction(.data$x, .data$y, drop = TRUE)),
      ndvi_beta = pmin(pmax(.data$ndvi, 1e-4), 0.9999)
    )

  # Time indices are required only by candidate specifications that include a
  # temporal random effect. Construct them before complete-case filtering to
  # preserve the original monthly index.
  if (isTRUE(include_time_id)) {
    covariates <- covariates |>
      dplyr::mutate(
        time_id = match(.data$date, sort(unique(.data$date)))
      ) |>
      dplyr::relocate("time_id", .before = "ndvi_beta")
  }

  tidyr::drop_na(covariates)
}

ndvi_standardize_candidate_covariates <- function(training, holdout = NULL) {
  variables <- c(
    "pr_lag1",
    "pr_roll3",
    "pr_roll6",
    "pr_recent_anom",
    "tasmax",
    "tasmax_roll3",
    "tasmin_roll3",
    "tasmean_roll3",
    "ele_km",
    "year"
  )
  centres <- vapply(training[variables], mean, numeric(1))
  spreads <- vapply(training[variables], stats::sd, numeric(1))

  invalid_spreads <- names(spreads)[!is.finite(spreads) | spreads <= 0]
  if (length(invalid_spreads) > 0L) {
    stop(
      "Cannot standardise zero-variance or non-finite covariates: ",
      paste(invalid_spreads, collapse = ", ")
    )
  }

  scale_data <- function(df) {
    df |>
      dplyr::mutate(
        pr_lag1_sc = (.data$pr_lag1 - centres[["pr_lag1"]]) /
          spreads[["pr_lag1"]],
        pr_roll3_sc = (.data$pr_roll3 - centres[["pr_roll3"]]) /
          spreads[["pr_roll3"]],
        pr_roll6_sc = (.data$pr_roll6 - centres[["pr_roll6"]]) /
          spreads[["pr_roll6"]],
        pr_recent_anom_sc = (.data$pr_recent_anom -
          centres[["pr_recent_anom"]]) / spreads[["pr_recent_anom"]],
        tasmax_sc = (.data$tasmax - centres[["tasmax"]]) /
          spreads[["tasmax"]],
        tasmax_roll3_sc = (.data$tasmax_roll3 -
          centres[["tasmax_roll3"]]) / spreads[["tasmax_roll3"]],
        tasmin_roll3_sc = (.data$tasmin_roll3 -
          centres[["tasmin_roll3"]]) / spreads[["tasmin_roll3"]],
        tasmean_roll3_sc = (.data$tasmean_roll3 -
          centres[["tasmean_roll3"]]) / spreads[["tasmean_roll3"]],
        ele_km_sc = (.data$ele_km - centres[["ele_km"]]) /
          spreads[["ele_km"]],
        year_sc = (.data$year - centres[["year"]]) /
          spreads[["year"]],
        pr_tasmax_inter = .data$pr_roll3_sc * .data$tasmax_roll3_sc,
        prrecent_tasmax_inter = .data$pr_recent_anom_sc *
          .data$tasmax_roll3_sc
      )
  }

  list(
    training = scale_data(training),
    holdout = if (is.null(holdout)) NULL else scale_data(holdout),
    centres = centres,
    spreads = spreads,
    scaling = tibble::tibble(
      variable = variables,
      mean = unname(centres),
      sd = unname(spreads)
    )
  )
}
