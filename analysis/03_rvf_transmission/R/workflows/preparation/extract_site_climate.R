#!/usr/bin/env Rscript
#' Extract monthly rainfall and mean temperature at transmission sites.
#'
#' The transmission model uses rainfall for vector habitat capacity and mean
#' temperature for biting, mortality, and extrinsic incubation. Daily maximum
#' and minimum temperature are quality-controlled upstream but are not copied
#' into Module 03 because they do not enter the equations.

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(readr)
  library(sf)
  library(terra)
  library(tidyr)
})

source(here::here("analysis", "03_rvf_transmission", "R", "config.R"))
ensure_rvf_dirs()

cfg <- rvf_config
workspace_inputs <- load_workspace_inputs(cfg)

extract_monthly_values <- function(raster, points, years, value_name,
                                   method = "bilinear") {
  raster_dates <- as.Date(terra::time(raster))
  raster_years <- lubridate::year(raster_dates)
  missing_years <- setdiff(years, unique(raster_years))

  if (length(missing_years) > 0L) {
    stop(
      "Climate raster is missing requested year(s): ",
      paste(missing_years, collapse = ", ")
    )
  }

  keep <- which(raster_years %in% years)
  if (length(keep) != length(years) * 12L) {
    stop("Climate raster does not contain 12 monthly layers for every requested year.")
  }

  raster <- raster[[keep]]
  raster_dates <- as.Date(terra::time(raster))
  layer_ids <- paste0("D", format(raster_dates, "%Y_%m"))
  names(raster) <- layer_ids

  point_vector <- terra::vect(points)
  point_extent <- terra::ext(point_vector)
  crop_extent <- terra::ext(
    point_extent$xmin - 0.2, point_extent$xmax + 0.2,
    point_extent$ymin - 0.2, point_extent$ymax + 0.2
  )
  raster <- terra::crop(raster, crop_extent)
  extracted <- terra::extract(raster, point_vector, method = method, ID = FALSE)

  dates <- tibble(
    layer_id = layer_ids,
    date = raster_dates,
    year = lubridate::year(raster_dates),
    month_num = lubridate::month(raster_dates),
    month = month.abb[month_num]
  )

  bind_cols(
    static_site_metadata(points) |>
      select(all_of(c("location_name", "site_type"))),
    as.data.frame(extracted, check.names = FALSE)
  ) |>
    pivot_longer(starts_with("D"), names_to = "layer_id", values_to = value_name) |>
    left_join(dates, by = "layer_id") |>
    select(-all_of("layer_id"))
}

combine_climate <- function(points, years, rainfall_file, temperature_file,
                            dataset, scenario, scenario_label, climate_model) {
  rainfall <- terra::rast(rainfall_file)
  temperature <- terra::rast(temperature_file) - 273.15

  rainfall_data <- extract_monthly_values(
    rainfall, points, years, "pr", cfg$extraction_method
  )
  temperature_data <- extract_monthly_values(
    temperature, points, years, "tas", cfg$extraction_method
  )

  climate <- rainfall_data |>
    left_join(
      select(temperature_data, all_of(c("location_name", "site_type", "date", "tas"))),
      by = c("location_name", "site_type", "date")
    ) |>
    mutate(
      dataset = dataset,
      scenario = scenario,
      scenario_label = scenario_label,
      climate_model = climate_model
    ) |>
    arrange(.data$site_type, .data$year, .data$month_num, .data$location_name)

  if (any(!is.finite(climate$pr)) || any(!is.finite(climate$tas))) {
    stop("Non-finite climate values after extraction for ", climate_model)
  }

  climate
}

save_site_climate <- function(climate, site_type, output_file) {
  site_climate <- climate |>
    filter(.data$site_type == .env$site_type)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(site_climate, output_file)

  tibble(
    dataset = site_climate$dataset[[1]],
    scenario = site_climate$scenario[[1]],
    scenario_label = site_climate$scenario_label[[1]],
    climate_model = site_climate$climate_model[[1]],
    site_type = site_type,
    n_locations = n_distinct(site_climate$location_name),
    year_min = min(site_climate$year),
    year_max = max(site_climate$year),
    rainfall_min = min(site_climate$pr),
    rainfall_max = max(site_climate$pr),
    temperature_min = min(site_climate$tas),
    temperature_max = max(site_climate$tas),
    file = output_file
  )
}

site_sets <- list(
  vil = workspace_inputs$village_ndvi_observed,
  gra = workspace_inputs$grazing_ndvi_observed,
  wat = workspace_inputs$watering_ndvi_observed
)

# Extract all locations in one pass per climate dataset. This avoids opening
# each large NetCDF three times while retaining separate output files for the
# transmission site's three ecological types.
all_sites <- bind_rows(lapply(names(site_sets), function(site_type) {
  site_sets[[site_type]] |>
    mutate(site_type = site_type) |>
    select(all_of(c("location_name", "site_type", attr(site_sets[[site_type]], "sf_column"))))
}))

summary_rows <- list()

message("Observed climate")
observed_climate <- combine_climate(
  points = all_sites,
  years = cfg$current_years,
  rainfall_file = cfg$observed_pr_file,
  temperature_file = cfg$observed_tas_file,
  dataset = "observed",
  scenario = NA_character_,
  scenario_label = "Observed",
  climate_model = "OBSERVED_CURRENT"
)
for (site_type in names(site_sets)) {
  output_file <- file.path(cfg$climate_dir, "observed", paste0(site_type, ".rds"))
  summary_rows[[length(summary_rows) + 1L]] <- save_site_climate(
    observed_climate, site_type, output_file
  )
}

for (scenario in cfg$future_scenarios) {
  for (climate_model in cfg$climate_models) {
    message(scenario, " | ", climate_model)
    model_file <- function(variable) {
      file.path(
        cfg$chelsa_future_root, scenario,
        paste0(climate_model, "_", variable, ".nc")
      )
    }
    future_climate_data <- combine_climate(
      points = all_sites,
      years = cfg$future_years,
      rainfall_file = model_file("pr"),
      temperature_file = model_file("tas"),
      dataset = "future",
      scenario = scenario,
      scenario_label = unname(cfg$scenario_labels[[scenario]]),
      climate_model = climate_model
    )

    for (site_type in names(site_sets)) {
      output_file <- file.path(
        cfg$climate_dir, scenario,
        paste0(climate_model, "_", site_type, ".rds")
      )

      summary_rows[[length(summary_rows) + 1L]] <- save_site_climate(
        future_climate_data, site_type, output_file
      )
    }
  }
}

climate_index <- bind_rows(summary_rows) |>
  arrange(.data$dataset, .data$scenario, .data$climate_model, .data$site_type)
write_csv(climate_index, file.path(cfg$climate_dir, "climate_input_index.csv"))

message("Climate inputs written to: ", cfg$climate_dir)
