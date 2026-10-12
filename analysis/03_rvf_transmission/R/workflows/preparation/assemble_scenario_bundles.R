#!/usr/bin/env Rscript
#' Assemble ODE-ready movement and climate scenarios.
#'
#' Each output contains one year of monthly commute matrices and site-level
#' rainfall and mean temperature for one SSP, GCM, and network replicate.
#' Module 02 remains the sole store of the original binary network draws.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

source(here::here("analysis", "03_rvf_transmission", "R", "config.R"))
ensure_rvf_dirs()

cfg <- rvf_config

row_normalise <- function(matrix, tolerance = 1e-12) {
  matrix <- as.matrix(matrix)
  matrix[!is.finite(matrix) | matrix < 0] <- 0
  totals <- rowSums(matrix)
  output <- matrix
  nonempty <- totals > tolerance
  output[nonempty, ] <- output[nonempty, , drop = FALSE] / totals[nonempty]
  output[!nonempty, ] <- 0
  output
}

validate_climate_table <- function(climate, path) {
  required_columns <- c("year", "month", "location_name", "pr", "tas")
  missing_columns <- setdiff(required_columns, names(climate))
  if (length(missing_columns) > 0L) {
    stop(
      "Site-climate file is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      "\nFile: ", path
    )
  }

  climate_key <- paste(climate$year, climate$month, climate$location_name, sep = "::")
  duplicated_key <- duplicated(climate_key)
  if (any(duplicated_key)) {
    stop(
      "Site-climate data contain duplicate year-month-location rows, for example: ",
      climate_key[which(duplicated_key)[1]],
      "\nFile: ", path
    )
  }

  invisible(climate)
}

# Select the requested climate year and order locations to match the movement
# matrices. Explicit `.env` lookup prevents dplyr from treating the requested
# year as another reference to the `year` column.
build_month_lists <- function(climate, scenario_year, target_names, site_type,
                              climate_model) {
  scenario_year <- as.integer(scenario_year)
  output <- lapply(month.abb, function(month_name) {
    monthly <- climate |>
      filter(
        .data$year == .env$scenario_year,
        .data$month == .env$month_name
      )

    if (nrow(monthly) == 0L) {
      stop(
        "No climate data for ", climate_model, ", ", scenario_year, ", ",
        month_name, ", site type ", site_type
      )
    }

    selected_years <- unique(as.integer(monthly$year))
    if (length(selected_years) != 1L || selected_years != scenario_year) {
      stop(
        "Climate-year selection failed for ", climate_model, ", requested ",
        scenario_year, ", site type ", site_type, ", month ", month_name,
        ". Selected year(s): ", paste(selected_years, collapse = ", ")
      )
    }

    duplicated_locations <- unique(monthly$location_name[duplicated(monthly$location_name)])
    if (length(duplicated_locations) > 0L) {
      stop(
        "Duplicate climate rows after year-month selection for ", site_type,
        ", for example: ", duplicated_locations[1]
      )
    }

    order_index <- match(target_names, monthly$location_name)
    if (anyNA(order_index)) {
      stop(
        "Climate data do not match the movement locations for ", site_type,
        ", for example: ", target_names[which(is.na(order_index))[1]]
      )
    }

    monthly <- monthly[order_index, , drop = FALSE]
    if (any(!is.finite(monthly$pr)) || any(!is.finite(monthly$tas))) {
      stop(
        "Non-finite climate values in ", climate_model, ", ",
        scenario_year, ", ", month_name
      )
    }
    monthly
  })
  names(output) <- month.abb
  output
}

# Summarise one scenario's climate independently of its movement-network
# replicate. The fingerprint uses monthly site-type means and is an integrity
# check, not a scientific summary or model input.
summarise_bundle_climate <- function(bundle_climate, expected_year) {
  site_types <- c("vil", "gra", "wat")
  fingerprint_values <- numeric()
  annual_means <- list()

  for (site_type in site_types) {
    monthly_tables <- bundle_climate[[site_type]]
    if (length(monthly_tables) != 12L) {
      stop("Expected 12 climate months for site type ", site_type, ".")
    }

    selected_years <- unique(unlist(lapply(monthly_tables, function(x) x$year)))
    if (length(selected_years) != 1L || as.integer(selected_years) != expected_year) {
      stop(
        "Bundle climate does not match metadata year ", expected_year,
        " for site type ", site_type, ". Selected year(s): ",
        paste(selected_years, collapse = ", ")
      )
    }

    monthly_means <- unlist(lapply(monthly_tables, function(x) {
      c(mean(x$tas), mean(x$pr))
    }))
    fingerprint_values <- c(fingerprint_values, monthly_means)

    annual_means[[paste0(site_type, "_mean_tas")]] <- mean(unlist(lapply(monthly_tables, `[[`, "tas")))
    annual_means[[paste0(site_type, "_mean_pr")]] <- mean(unlist(lapply(monthly_tables, `[[`, "pr")))
  }

  as_tibble(annual_means) |>
    mutate(
      climate_fingerprint = paste(
        formatC(fingerprint_values, digits = 10, format = "fg", flag = "#"),
        collapse = "|"
      )
    )
}

climate_file <- function(dataset, scenario, climate_model, site_type) {
  if (identical(dataset, "observed")) {
    file.path(cfg$climate_dir, "observed", paste0(site_type, ".rds"))
  } else {
    file.path(
      cfg$climate_dir, scenario,
      paste0(climate_model, "_", site_type, ".rds")
    )
  }
}

climate_cache <- new.env(parent = emptyenv())

get_climate <- function(dataset, scenario, climate_model, site_type) {
  key <- paste(dataset, scenario, climate_model, site_type, sep = "::")
  if (!exists(key, envir = climate_cache, inherits = FALSE)) {
    path <- climate_file(dataset, scenario, climate_model, site_type)
    if (!file.exists(path)) stop("Missing site-climate file: ", path)
    climate <- readRDS(path)
    validate_climate_table(climate, path)
    assign(key, climate, envir = climate_cache)
  }
  get(key, envir = climate_cache, inherits = FALSE)
}

network_index_file <- file.path(cfg$networks_dir, "network_scenario_index.csv")
if (!file.exists(network_index_file)) {
  stop("Run 01_prepare_inputs.R first: ", network_index_file)
}

network_index <- read_csv(network_index_file, show_col_types = FALSE)
if (nrow(network_index) == 0L) stop("The network scenario index is empty.")

# Dataset scope is declared once in R/config.R. The publication workflow
# rebuilds observed and future bundles together so their metadata stay aligned.
assemble_datasets <- cfg$execution$assemble_datasets
if (length(assemble_datasets) > 0L) {
  network_index <- network_index[network_index$dataset %in% assemble_datasets, , drop = FALSE]
  if (nrow(network_index) == 0L) {
    stop("No scenario-index rows match execution$assemble_datasets: ",
         paste(assemble_datasets, collapse = ", "))
  }
  message("Rebuilding scenario bundles for dataset(s): ", paste(assemble_datasets, collapse = ", "))
}

# Scenario files are derived inputs. Clear this cache only when all replacement
# inputs for the selected rebuild pass are present, immediately before
# rebuilding it. This lets observed and future bundles be rebuilt separately.
required_climate <- unique(unlist(lapply(seq_len(nrow(network_index)), function(i) {
  row <- network_index[i, ]
  scenario <- if (is.na(row$scenario)) NA_character_ else row$scenario
  c(
    climate_file(row$dataset, scenario, row$climate_model, "vil"),
    climate_file(row$dataset, scenario, row$climate_model, "gra"),
    climate_file(row$dataset, scenario, row$climate_model, "wat")
  )
}), use.names = FALSE))
if (any(!file.exists(required_climate))) {
  stop("Run 01_prepare_inputs.R before rebuilding scenario bundles.")
}

# Remove only the bundles this pass will regenerate, so a scoped rebuild (e.g.
# future) never deletes another dataset's bundles (e.g. the observed set the
# smoke test / calibration depends on).
regenerated <- file.path(cfg$scenarios_dir, paste0(network_index$scenario_id, ".rds"))
stale <- regenerated[file.exists(regenerated)]
if (length(stale) > 0L) unlink(stale)

scenario_rows <- vector("list", nrow(network_index))
climate_audit_cache <- new.env(parent = emptyenv())
loaded_source <- NULL
source_bundle <- NULL

for (i in seq_len(nrow(network_index))) {
  row <- network_index[i, ]
  if (!identical(loaded_source, row$source_bundle)) {
    source_bundle <- readRDS(row$source_bundle)
    loaded_source <- row$source_bundle
  }

  replicate_index <- as.integer(row$network_replicate)
  grazing <- lapply(source_bundle$grazing$replicates[[replicate_index]], row_normalise)
  watering <- lapply(source_bundle$watering$replicates[[replicate_index]], row_normalise)

  scenario <- if (is.na(row$scenario)) NA_character_ else row$scenario
  village_climate <- get_climate(row$dataset, scenario, row$climate_model, "vil")
  grazing_climate <- get_climate(row$dataset, scenario, row$climate_model, "gra")
  watering_climate <- get_climate(row$dataset, scenario, row$climate_model, "wat")

  metadata <- source_bundle$metadata
  metadata$scenario_id <- row$scenario_id
  metadata$scenario_base_id <- row$scenario_base_id
  metadata$dataset <- row$dataset
  metadata$scenario <- scenario
  metadata$scenario_label <- row$scenario_label
  metadata$climate_model <- row$climate_model
  metadata$scenario_year <- as.integer(row$year)
  metadata$future_window <- row$future_window
  metadata$network_replicate <- replicate_index
  metadata$source_network_bundle <- row$source_bundle

  bundle <- list(
    G_list = grazing,
    W_list = watering,
    climate = list(
      vil = build_month_lists(
        village_climate, row$year, source_bundle$vill_names,
        "village", row$climate_model
      ),
      gra = build_month_lists(
        grazing_climate, row$year, source_bundle$grazing$cent_names,
        "grazing", row$climate_model
      ),
      wat = build_month_lists(
        watering_climate, row$year, source_bundle$watering$cent_names,
        "watering", row$climate_model
      )
    ),
    n_vil = length(source_bundle$vill_names),
    n_gra = length(source_bundle$grazing$cent_names),
    n_wat = length(source_bundle$watering$cent_names),
    vil_names = source_bundle$vill_names,
    gra_names = source_bundle$grazing$cent_names,
    wat_names = source_bundle$watering$cent_names,
    metadata = metadata
  )

  # Climate is identical across network replicates for the same scenario year,
  # so record its integrity summary only once.
  audit_key <- paste(
    row$dataset, scenario, row$climate_model, as.integer(row$year),
    sep = "::"
  )
  if (!exists(audit_key, envir = climate_audit_cache, inherits = FALSE)) {
    climate_summary <- summarise_bundle_climate(
      bundle$climate,
      expected_year = as.integer(row$year)
    )
    audit_row <- bind_cols(
      tibble(
        dataset = row$dataset,
        scenario = scenario,
        scenario_label = row$scenario_label,
        climate_model = row$climate_model,
        future_window = row$future_window,
        year = as.integer(row$year)
      ),
      climate_summary
    )
    assign(audit_key, audit_row, envir = climate_audit_cache)
  }

  output_file <- file.path(cfg$scenarios_dir, paste0(row$scenario_id, ".rds"))
  saveRDS(bundle, output_file, compress = "gzip")

  scenario_rows[[i]] <- tibble(
    scenario_id = row$scenario_id,
    scenario_base_id = row$scenario_base_id,
    dataset = row$dataset,
    scenario = scenario,
    scenario_label = row$scenario_label,
    climate_model = row$climate_model,
    future_window = row$future_window,
    year = as.integer(row$year),
    network_replicate = replicate_index,
    n_vil = bundle$n_vil,
    n_gra = bundle$n_gra,
    n_wat = bundle$n_wat,
    file = output_file
  )

  if (i == 1L || i %% 500L == 0L || i == nrow(network_index)) {
    message("Assembled ", i, " / ", nrow(network_index), " scenarios")
  }
}

scenario_index <- bind_rows(scenario_rows)
# Merge with any existing index so a scoped rebuild preserves the other
# dataset's rows -- observed and future must both stay listed for downstream.
index_path <- file.path(cfg$scenarios_dir, "scenario_bundle_index.csv")
if (file.exists(index_path)) {
  existing_index <- read_csv(index_path, show_col_types = FALSE)
  existing_index <- existing_index[!existing_index$dataset %in% unique(scenario_index$dataset), , drop = FALSE]
  scenario_index <- bind_rows(existing_index, scenario_index)
}
scenario_index <- dplyr::arrange(
  scenario_index, .data$dataset, .data$scenario, .data$climate_model,
  .data$year, .data$network_replicate
)

# Preserve the other dataset during a scoped observed/future rebuild, then
# require a one-to-one match between the scenario index and climate audit.
climate_audit <- bind_rows(as.list(climate_audit_cache))
audit_path <- file.path(cfg$scenarios_dir, "scenario_climate_integrity.csv")
if (file.exists(audit_path)) {
  existing_audit <- read_csv(audit_path, show_col_types = FALSE)
  existing_audit <- existing_audit[
    !existing_audit$dataset %in% unique(climate_audit$dataset),
    , drop = FALSE
  ]
  climate_audit <- bind_rows(existing_audit, climate_audit)
}

expected_climate_keys <- scenario_index |>
  distinct(
    .data$dataset, .data$scenario, .data$scenario_label,
    .data$climate_model, .data$future_window, .data$year
  )
climate_audit <- climate_audit |>
  semi_join(
    expected_climate_keys,
    by = c("dataset", "scenario", "scenario_label", "climate_model", "future_window", "year")
  ) |>
  arrange(.data$dataset, .data$scenario, .data$climate_model, .data$year)

missing_audit_rows <- anti_join(
  expected_climate_keys,
  climate_audit,
  by = c("dataset", "scenario", "scenario_label", "climate_model", "future_window", "year")
)
if (nrow(missing_audit_rows) > 0L) {
  stop("Climate-integrity audit is missing ", nrow(missing_audit_rows), " scenario-year row(s).")
}

# An exact repeat of all monthly site-type means in different years is an
# assembly error in this pipeline, not plausible interannual climate variation.
repeated_fingerprints <- climate_audit |>
  group_by(
    .data$dataset, .data$scenario, .data$climate_model,
    .data$climate_fingerprint
  ) |>
  summarise(
    n_years = n_distinct(.data$year),
    years = paste(sort(unique(.data$year)), collapse = ", "),
    .groups = "drop"
  ) |>
  filter(.data$n_years > 1L)
if (nrow(repeated_fingerprints) > 0L) {
  first_repeat <- repeated_fingerprints[1, ]
  stop(
    "Climate-integrity audit found an identical annual climate sequence in ",
    first_repeat$n_years, " years for dataset=", first_repeat$dataset,
    ", scenario=", first_repeat$scenario,
    ", climate_model=", first_repeat$climate_model,
    ". Years: ", first_repeat$years
  )
}

write_csv(scenario_index, index_path)
write_csv(climate_audit, audit_path)
message("Scenario bundles written to: ", cfg$scenarios_dir)
message("Climate-integrity audit written to: ", audit_path)
