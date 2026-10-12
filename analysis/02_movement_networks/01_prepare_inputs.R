#!/usr/bin/env Rscript
# Prepare the observed calibration data and site-level NDVI for both SSPs.
#
# The participatory data contain monthly village-to-village contacts. They do
# not identify the shared resource sites. This script therefore preserves the
# observed graphs at their original one-mode scale and separately prepares the
# candidate village, grazing, and water locations used by the latent simulator.

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(dplyr)
  library(geodist)
  library(igraph)
  library(lubridate)
  library(readr)
  library(sf)
  library(stringr)
  library(terra)
})

source(here::here("analysis", "02_movement_networks", "_config.R"))
source(here::here("analysis", "02_movement_networks", "R", "utils.R"))

cfg <- mvt_config

month_column <- function(date) {
  paste0(month.abb[lubridate::month(date)], "_", lubridate::year(date))
}

check_complete_months <- function(data, years, label) {
  expected <- as.Date(sprintf(
    "%d-%02d-01", rep(years, each = 12), rep(seq_len(12), length(years))
  ))
  missing <- setdiff(expected, sort(unique(data$date)))
  if (length(missing) > 0L) {
    stop(label, " is missing: ", paste(format(missing, "%Y-%m"), collapse = ", "))
  }
}

normalise_name <- function(x) {
  x |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_replace_all("^_|_$", "")
}

optional_numeric <- function(data, candidates) {
  column <- candidates[candidates %in% names(data)][1]
  if (is.na(column)) return(rep(NA_real_, nrow(data)))
  suppressWarnings(as.numeric(data[[column]]))
}

# Combine the site-name patterns that mark restricted-access sites into one
# regular expression. The file holds one pattern per line and comes with the
# participatory data, so the patterns themselves never appear in this code.
read_restricted_site_patterns <- function(path) {
  if (!file.exists(path)) {
    stop("Missing ", path, ". This file is supplied with the participatory data.")
  }
  patterns <- trimws(readLines(path, warn = FALSE))
  patterns <- patterns[nzchar(patterns)]
  if (length(patterns) == 0L) {
    stop(path, " contains no site-name patterns.")
  }
  paste(patterns, collapse = "|")
}

load_site_templates <- function() {
  district <- sf::st_read(
    file.path(cfg$shape_dir, "Serengeti and Ngorongoro district shape files.shp"),
    quiet = TRUE
  ) |>
    dplyr::filter(.data$District_N == "Serengeti") |>
    sf::st_transform(4326)

  administration <- district |>
    dplyr::select(
      admin_district = District_N,
      admin_ward = Ward_Name,
      admin_village = Vil_Mtaa_N
    )

  villages <- readr::read_csv(
    file.path(cfg$shape_dir, "serengeti village centre xy.csv"),
    show_col_types = FALSE
  ) |>
    dplyr::rename(village_name = Village_name.x) |>
    sf::st_as_sf(coords = c("EASTING.x", "NORTHING.x"), crs = 32736) |>
    sf::st_transform(4326) |>
    dplyr::mutate(location_name = .data$village_name) |>
    dplyr::filter(!.data$location_name %in%
                    c("Bunda_Irimba", "Robanda", "Butiama_District")) |>
    sf::st_join(administration, join = sf::st_within, left = TRUE) |>
    dplyr::distinct(.data$location_name, .keep_all = TRUE) |>
    dplyr::select(
      village_name, location_name, admin_district,
      admin_ward, admin_village, geometry
    )

  restricted_regex <- read_restricted_site_patterns(
    cfg$restricted_site_patterns_file
  )

  grazing <- sf::st_read(
    file.path(cfg$shape_dir, "All grazing areas_centroid.shp"), quiet = TRUE
  ) |>
    sf::st_transform(4326) |>
    dplyr::filter(!is.na(.data$xcoord), !is.na(.data$ycoord)) |>
    dplyr::mutate(
      location_name = normalise_name(make.unique(.data$Location_C)),
      owner_village_raw = .data$Village_Na,
      owner_village_id = .data$Village_ID,
      owner_village_key = normalise_name(.data$owner_village_raw),
      is_restricted = stringr::str_detect(.data$location_name, restricted_regex)
    )
  grazing$area_size <- optional_numeric(
    grazing, c("Areas_Size", "Area_size", "area_size")
  )
  grazing <- grazing |>
    dplyr::select(
      location_name, owner_village_raw, owner_village_id,
      owner_village_key, area_size, type = Type,
      grazing_category = Grazing_Ca, is_restricted, geometry
    )

  # The canonical movement analysis uses the natural shared-water inventory.
  # The small set of constructed pumps and settlement pools is not used to
  # infer shared commutes because the observed one-mode data cannot identify
  # which of those sites generated a reported village contact.
  watering <- sf::st_read(
    file.path(cfg$shape_dir, "All water points.shp"), quiet = TRUE
  ) |>
    sf::st_transform(4326) |>
    dplyr::filter(!is.na(.data$xcoord), !is.na(.data$ycoord)) |>
    dplyr::mutate(
      location_name = normalise_name(make.unique(.data$Location_c)),
      owner_village_raw = .data$Village_Na,
      owner_village_id = .data$Village_ID,
      owner_village_key = normalise_name(.data$owner_village_raw),
      is_natural = !stringr::str_detect(
        .data$location_name,
        "water_pump|pool_spills?|perimeter_pool|mcu_town|nhc_village"
      ),
      is_restricted = stringr::str_detect(.data$location_name, restricted_regex)
    ) |>
    dplyr::filter(.data$is_natural)
  watering$area_size <- optional_numeric(
    watering, c("Area_size", "Areas_Size", "area_size")
  )
  watering <- watering |>
    dplyr::select(
      location_name, owner_village_raw, owner_village_id,
      owner_village_key, area_size, is_restricted, geometry
    )

  list(village = villages, grazing = grazing, watering = watering)
}

grid_to_raster <- function(data) {
  xyz <- data |>
    dplyr::select(x, y, ndvi_value) |>
    as.data.frame()
  names(xyz) <- c("x", "y", "value")
  terra::rast(xyz, type = "xyz", crs = "EPSG:4326")
}

extract_point_values <- function(raster, sites) {
  values <- terra::extract(raster, terra::vect(sites), method = "bilinear")[[2]]
  if (anyNA(values)) {
    nearest <- terra::extract(raster, terra::vect(sites), method = "simple")[[2]]
    values[is.na(values)] <- nearest[is.na(values)]
  }
  values
}

extract_site_ndvi <- function(monthly_grid, templates) {
  result <- templates
  for (date_value in sort(unique(monthly_grid$date))) {
    current_date <- as.Date(date_value, origin = "1970-01-01")
    raster <- grid_to_raster(
      dplyr::filter(monthly_grid, .data$date == .env$current_date)
    )
    column <- month_column(current_date)
    for (site_type in names(result)) {
      result[[site_type]][[column]] <- extract_point_values(
        raster, result[[site_type]]
      )
    }
  }
  result
}

read_observed_networks <- function(path, village_names) {
  excluded <- c("Bunda_Irimba", "Robanda", "Butiama_District")
  edges <- data.table::fread(path) |>
    as.data.frame() |>
    dplyr::transmute(
      from = .data$Village_name.x,
      to = .data$Village_name.y,
      month = month.abb[match(.data$Month, month.name)]
    ) |>
    dplyr::filter(!.data$from %in% excluded, !.data$to %in% excluded)

  unknown <- setdiff(unique(c(edges$from, edges$to)), village_names)
  if (length(unknown) > 0L) {
    stop("Movement data contain villages absent from the site file: ",
         paste(unknown, collapse = ", "))
  }

  networks <- lapply(month.abb, function(month_name) {
    month_edges <- dplyr::filter(
      edges, .data$month == .env$month_name
    )
    if (nrow(month_edges) == 0L) {
      graph <- igraph::make_empty_graph(n = length(village_names), directed = FALSE)
      igraph::V(graph)$name <- village_names
      return(graph)
    }
    make_graph(month_edges[c("from", "to")], village_names)
  })
  names(networks) <- month.abb
  networks
}

site_distance <- function(origins, destinations = origins) {
  origin_xy <- sf::st_coordinates(origins)
  destination_xy <- sf::st_coordinates(destinations)
  distance <- geodist::geodist(origin_xy, destination_xy, measure = "haversine") / 1000
  rownames(distance) <- origins$location_name
  colnames(distance) <- destinations$location_name
  distance
}

drop_ndvi_columns <- function(data) {
  pattern <- paste0("^(", paste(month.abb, collapse = "|"), ")_[0-9]{4}$")
  sf::st_drop_geometry(data) |>
    dplyr::select(-dplyr::matches(pattern))
}

templates <- load_site_templates()

historical <- arrow::read_feather(cfg$historical_ndvi_file) |>
  dplyr::transmute(
    x = .data$x, y = .data$y, date = as.Date(.data$date),
    year = lubridate::year(.data$date), ndvi_value = .data$ndvi
  ) |>
  dplyr::filter(.data$year %in% cfg$observed_years)
check_complete_months(historical, cfg$observed_years, "Historical NDVI")

message("Extracting observed site-level NDVI...")
observed_ndvi <- extract_site_ndvi(historical, templates)
village_names <- observed_ndvi$village$location_name

observed_inputs <- list(
  metadata = list(
    dataset = "observed",
    input_schema_version = cfg$input_schema_version,
    years = cfg$observed_years,
    ndvi_source = cfg$historical_ndvi_file,
    prepared_at = as.character(Sys.time())
  ),
  village_names = village_names,
  distance_village_village = site_distance(observed_ndvi$village),
  village_metadata = drop_ndvi_columns(observed_ndvi$village),
  grazing = list(
    distance_village_resource = site_distance(
      observed_ndvi$village, observed_ndvi$grazing
    ),
    is_restricted = observed_ndvi$grazing$is_restricted,
    resource_names = observed_ndvi$grazing$location_name,
    resource_metadata = drop_ndvi_columns(observed_ndvi$grazing)
  ),
  watering = list(
    distance_village_resource = site_distance(
      observed_ndvi$village, observed_ndvi$watering
    ),
    is_restricted = observed_ndvi$watering$is_restricted,
    resource_names = observed_ndvi$watering$location_name,
    resource_metadata = drop_ndvi_columns(observed_ndvi$watering)
  ),
  ndvi = observed_ndvi,
  observed_networks = list(
    grazing = read_observed_networks(cfg$movement_files[["grazing"]], village_names),
    watering = read_observed_networks(cfg$movement_files[["watering"]], village_names)
  )
)
saveRDS(observed_inputs, cfg$observed_inputs_file)

jobs <- tidyr::crossing(
  scenario = cfg$scenarios,
  climate_model = cfg$climate_models
)

prepare_scenario <- function(index) {
  job <- jobs[index, ]
  input_file <- file.path(
    cfg$future_ndvi_dir, job$scenario, "projections", job$climate_model,
    "monthly_ndvi.feather"
  )
  output_file <- scenario_input_file(job$scenario, job$climate_model)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(input_file)) stop("Missing Module 01 projection: ", input_file)
  if (file.exists(output_file) && file.mtime(output_file) >= file.mtime(input_file)) {
    existing <- tryCatch(readRDS(output_file), error = function(error) NULL)
    if (!is.null(existing) &&
        identical(existing$metadata$input_schema_version,
                  cfg$input_schema_version)) {
      message("Up to date: ", job$scenario, " / ", job$climate_model)
      return(existing$metadata)
    }
  }

  message("Preparing ", job$scenario, " / ", job$climate_model)
  future <- arrow::read_feather(input_file) |>
    dplyr::transmute(
      x = .data$x, y = .data$y, date = as.Date(.data$date),
      year = as.integer(.data$year), ndvi_value = .data$ndvi_pred
    ) |>
    dplyr::filter(.data$year %in% cfg$future_years)
  check_complete_months(
    future, cfg$future_years,
    paste(cfg$scenario_labels[[job$scenario]], job$climate_model)
  )

  site_ndvi <- extract_site_ndvi(future, templates)
  result <- list(
    metadata = list(
      dataset = "future",
      input_schema_version = cfg$input_schema_version,
      scenario = job$scenario,
      scenario_label = unname(cfg$scenario_labels[[job$scenario]]),
      climate_model = job$climate_model,
      years = cfg$future_years,
      source_file = input_file,
      prepared_at = as.character(Sys.time())
    ),
    ndvi = site_ndvi
  )
  saveRDS(result, output_file)
  result$metadata
}

n_cores <- min(cfg$n_cores, nrow(jobs))
future_metadata <- if (.Platform$OS.type == "unix" && n_cores > 1L) {
  parallel::mclapply(seq_len(nrow(jobs)), prepare_scenario,
                     mc.cores = n_cores, mc.preschedule = FALSE)
} else {
  lapply(seq_len(nrow(jobs)), prepare_scenario)
}

# A try-error marks a failed job; NULL marks a worker process that stopped
# without returning, which mclapply reports only as a warning.
failed <- vapply(future_metadata, function(result) {
  is.null(result) || inherits(result, "try-error")
}, logical(1))
if (any(failed)) {
  stop(sum(failed), " future input job(s) failed or their worker process ",
       "stopped. Re-run this script.")
}

input_index <- dplyr::bind_rows(lapply(future_metadata, function(metadata) {
  tibble::tibble(
    scenario = metadata$scenario,
    scenario_label = metadata$scenario_label,
    climate_model = metadata$climate_model,
    source_file = metadata$source_file,
    prepared_at = metadata$prepared_at
  )
}))
readr::write_csv(input_index, file.path(cfg$inputs_dir, "future_input_index.csv"))

message("Prepared observed inputs: ", cfg$observed_inputs_file)
message("Prepared ", nrow(input_index), " SSP-GCM input files.")
