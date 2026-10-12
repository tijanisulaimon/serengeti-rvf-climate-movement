#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
})

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

# Reload `R/config.R` on demand so helper scripts pick up the current
# defaults without needing to manage a shared global object.
load_rvf_config <- function() {
  source(here::here("analysis", "03_rvf_transmission", "R", "config.R"))
  rvf_config
}

rvf_resolve_site_profile_id <- function(site_profile,
                                        cfg = load_rvf_config()) {
  if (is.null(site_profile)) return(NULL)
  site_profile <- as.character(site_profile)[1]
  if (!nzchar(site_profile)) return(site_profile)

  aliases <- cfg$site_profile_aliases %||% list()
  canonical <- aliases[[site_profile]] %||% site_profile
  if (!canonical %in% names(cfg$site_profiles)) {
    allowed <- sort(unique(c(names(cfg$site_profiles), names(aliases))))
    stop(
      "Unknown site-capacity profile: ", site_profile,
      "\nKnown profiles/aliases: ", paste(allowed, collapse = ", ")
    )
  }
  canonical
}

rvf_resolve_site_profile_ids <- function(site_profiles,
                                         cfg = load_rvf_config()) {
  if (is.null(site_profiles)) return(NULL)
  unique(vapply(site_profiles, rvf_resolve_site_profile_id, character(1), cfg = cfg))
}

# Row-normalise a matrix but leave empty rows as zeros.
normalize_rows_scenario <- function(M, tol = 1e-12) {
  rs <- rowSums(M, na.rm = TRUE)
  out <- M
  ok <- rs > tol
  out[ok, ] <- out[ok, , drop = FALSE] / rs[ok]
  out[!ok, ] <- 0
  out
}

## Read the scenario index written by 01_prepare_inputs.R.
list_rvf_scenarios <- function(cfg = load_rvf_config()) {
  index_file <- file.path(cfg$scenarios_dir, "scenario_bundle_index.csv")
  read.csv(index_file, stringsAsFactors = FALSE)
}

# Format replicate tags consistently on disk and in output tables.
rvf_replicate_suffix <- function(network_replicate) {
  sprintf("REP%03d", as.integer(network_replicate))
}

# Add a replicate tag to one or more base scenario ids.
rvf_with_replicate <- function(base_scenario_id, network_replicate) {
  # Scenario filenames on disk are now replicate-tagged so one base climate
  # year can be re-run under multiple network realisations.
  paste0(base_scenario_id, "_", rvf_replicate_suffix(network_replicate))
}

## Strip the replicate tag when collapsing back to a base scenario.
rvf_strip_replicate <- function(scenario_id) {
  sub("_REP[0-9]+$", "", scenario_id)
}

# Shared index filter so all runner scripts subset scenarios the same way.
rvf_select_scenarios <- function(index = list_rvf_scenarios(),
                                 dataset = NULL,
                                 climate_model = NULL,
                                 year = NULL,
                                 network_replicate = NULL,
                                 scenario_base_id = NULL,
                                 scenario = NULL) {
  # Lightweight scenario-index filter used by the runner scripts so they
  # all pick observed / future / replicate subsets in the same way.
  out <- index

  if (!is.null(dataset)) {
    out <- out[out$dataset %in% dataset, , drop = FALSE]
  }
  # SSP filter. The index carries a `scenario` column (ssp245 / ssp370) for
  # future rows so a single window/GCM can be picked per SSP; observed rows
  # have scenario = NA and are unaffected when `scenario` is left NULL.
  if (!is.null(scenario) && "scenario" %in% names(out)) {
    out <- out[out$scenario %in% scenario, , drop = FALSE]
  }
  if (!is.null(climate_model) && "climate_model" %in% names(out)) {
    out <- out[out$climate_model %in% climate_model, , drop = FALSE]
  }
  if (!is.null(year)) {
    out <- out[out$year %in% year, , drop = FALSE]
  }
  if (!is.null(network_replicate) && "network_replicate" %in% names(out)) {
    out <- out[out$network_replicate %in% network_replicate, , drop = FALSE]
  }
  if (!is.null(scenario_base_id) && "scenario_base_id" %in% names(out)) {
    out <- out[out$scenario_base_id %in% scenario_base_id, , drop = FALSE]
  }

  order_cols <- intersect(c("dataset", "climate_model", "year", "network_replicate", "scenario_id"), names(out))
  if (length(order_cols) > 0 && nrow(out) > 0) {
    out <- out[do.call(order, unname(out[order_cols])), , drop = FALSE]
  }
  out
}

## Default waning means effectively durable post-infection immunity.
rvf_is_default_waning <- function(waning_years) {
  waning_years <- suppressWarnings(as.numeric(waning_years))
  !is.finite(waning_years) || waning_years <= 0
}

## Human-readable label for output filenames and tables.
rvf_waning_label <- function(waning_years) {
  if (rvf_is_default_waning(waning_years)) "lifelong" else sprintf("%sy", as.character(waning_years))
}

## Validate movement matrices, climate tables, and metadata before ODE assembly.
validate_rvf_scenario <- function(scenario) {
  if (!is.list(scenario)) stop("Scenario must be a list.")
  required_components <- c(
    "G_list", "W_list", "climate", "n_vil", "n_gra", "n_wat",
    "vil_names", "gra_names", "wat_names", "metadata"
  )
  missing_components <- setdiff(required_components, names(scenario))
  if (length(missing_components) > 0L) {
    stop(
      "Scenario is missing required component(s): ",
      paste(missing_components, collapse = ", ")
    )
  }

  dimensions <- as.integer(c(
    scenario$n_vil, scenario$n_gra, scenario$n_wat
  ))
  if (length(dimensions) != 3L || anyNA(dimensions) || any(dimensions < 1L)) {
    stop("Scenario site counts must be positive integers.")
  }
  names(dimensions) <- c("vil", "gra", "wat")
  if (length(scenario$vil_names) != dimensions[["vil"]] ||
      length(scenario$gra_names) != dimensions[["gra"]] ||
      length(scenario$wat_names) != dimensions[["wat"]]) {
    stop("Scenario site-name vectors do not match the declared site counts.")
  }

  n_steps <- length(scenario$G_list)
  if (n_steps < 1L || n_steps %% 12L != 0L) {
    stop("Scenario must contain one or more complete 12-month forcing cycles.")
  }
  forcing_lists <- list(
    G = scenario$G_list,
    W = scenario$W_list,
    vil = scenario$climate$vil,
    gra = scenario$climate$gra,
    wat = scenario$climate$wat
  )
  invalid_lengths <- names(forcing_lists)[vapply(
    forcing_lists, length, integer(1)
  ) != n_steps]
  if (length(invalid_lengths) > 0L) {
    stop(
      "Scenario forcing lists have inconsistent lengths: ",
      paste(invalid_lengths, collapse = ", ")
    )
  }

  forcing_names <- lapply(forcing_lists, names)
  named_lists <- !vapply(forcing_names, is.null, logical(1))
  if (any(named_lists) &&
      !all(vapply(forcing_names[named_lists], identical, logical(1),
                  forcing_names[[which(named_lists)[1]]]))) {
    stop("Scenario movement and climate forcing labels are not aligned.")
  }

  validate_matrix <- function(matrix, expected_rows, expected_cols, label) {
    matrix <- as.matrix(matrix)
    if (!identical(dim(matrix), c(expected_rows, expected_cols))) {
      stop(
        label, " has dimensions ", paste(dim(matrix), collapse = " x "),
        "; expected ", expected_rows, " x ", expected_cols, "."
      )
    }
    if (any(!is.finite(matrix)) || any(matrix < 0)) {
      stop(label, " contains non-finite or negative movement weights.")
    }
    row_totals <- rowSums(matrix)
    normalised <- abs(row_totals - 1) <= 1e-8
    empty <- abs(row_totals) <= 1e-8
    if (any(!(normalised | empty))) {
      stop(label, " contains rows that sum to neither zero nor one.")
    }
  }

  validate_climate <- function(table, expected_rows, label) {
    required_columns <- c("location_name", "year", "month", "pr", "tas")
    missing_columns <- setdiff(required_columns, names(table))
    if (length(missing_columns) > 0L) {
      stop(
        label, " is missing climate column(s): ",
        paste(missing_columns, collapse = ", ")
      )
    }
    if (nrow(table) != expected_rows ||
        length(unique(table$location_name)) != expected_rows) {
      stop(label, " does not contain exactly one row per expected location.")
    }
    if (any(!is.finite(table$pr)) || any(!is.finite(table$tas))) {
      stop(label, " contains non-finite rainfall or temperature values.")
    }
  }

  for (step in seq_len(n_steps)) {
    step_label <- names(scenario$G_list)[step] %||% as.character(step)
    validate_matrix(
      scenario$G_list[[step]],
      dimensions[["vil"]],
      dimensions[["gra"]],
      paste0("G_list[[", step_label, "]]")
    )
    validate_matrix(
      scenario$W_list[[step]],
      dimensions[["vil"]],
      dimensions[["wat"]],
      paste0("W_list[[", step_label, "]]")
    )
    validate_climate(
      scenario$climate$vil[[step]],
      dimensions[["vil"]],
      paste0("climate$vil[[", step_label, "]]")
    )
    validate_climate(
      scenario$climate$gra[[step]],
      dimensions[["gra"]],
      paste0("climate$gra[[", step_label, "]]")
    )
    validate_climate(
      scenario$climate$wat[[step]],
      dimensions[["wat"]],
      paste0("climate$wat[[", step_label, "]]")
    )
  }

  # One-year bundles must contain climate from the year recorded in metadata.
  # Multi-year sequence objects do not carry `scenario_year` and are validated
  # when their component bundles are loaded.
  expected_year <- suppressWarnings(as.integer(scenario$metadata$scenario_year %||% NA_integer_))
  if (length(expected_year) == 1L && is.finite(expected_year)) {
    for (site_type in c("vil", "gra", "wat")) {
      for (month_index in seq_len(n_steps)) {
        monthly_climate <- scenario$climate[[site_type]][[month_index]]
        if (!"year" %in% names(monthly_climate)) {
          stop("Scenario climate is missing its year column for site type ", site_type, ".")
        }
        selected_years <- unique(as.integer(monthly_climate$year))
        if (length(selected_years) != 1L || selected_years != expected_year) {
          stop(
            "Scenario climate year does not match metadata: expected ", expected_year,
            ", found ", paste(selected_years, collapse = ", "),
            " for site type ", site_type, "."
          )
        }
      }
    }
  }
  invisible(TRUE)
}

## Read one processed scenario bundle from disk.
load_rvf_scenario <- function(scenario_id, cfg = load_rvf_config(), validate = TRUE) {
  scenario_file <- file.path(cfg$scenarios_dir, paste0(scenario_id, ".rds"))

  if (!file.exists(scenario_file)) {
    stop("Scenario not found: ", scenario_file)
  }

  scenario <- readRDS(scenario_file)

  if (validate) {
    validate_rvf_scenario(scenario)
  }

  scenario
}

## Concatenate multiple one-year scenarios into a repeated monthly sequence.
build_rvf_sequence_scenario <- function(scenario_ids, cfg = load_rvf_config(), validate = TRUE,
                                        scenario_type = "climate_sequence",
                                        sequence_label = NULL) {
  bundles <- lapply(scenario_ids, function(id) load_rvf_scenario(id, cfg = cfg, validate = validate))

  ref <- bundles[[1]]
  for (b in bundles[-1]) {
    stopifnot(b$n_vil == ref$n_vil, b$n_gra == ref$n_gra, b$n_wat == ref$n_wat)
    stopifnot(identical(b$vil_names, ref$vil_names))
    stopifnot(identical(b$gra_names, ref$gra_names))
    stopifnot(identical(b$wat_names, ref$wat_names))
  }

  G_list <- list()
  W_list <- list()
  climate_vil <- list()
  climate_gra <- list()
  climate_wat <- list()

  for (bundle in bundles) {
    year_value <- bundle$metadata$scenario_year
    for (m in month.abb) {
      key <- paste0(m, "_", year_value)
      G_list[[key]] <- bundle$G_list[[m]]
      W_list[[key]] <- bundle$W_list[[m]]
      climate_vil[[key]] <- bundle$climate$vil[[m]]
      climate_gra[[key]] <- bundle$climate$gra[[m]]
      climate_wat[[key]] <- bundle$climate$wat[[m]]
    }
  }

  scenario <- list(
    G_list = G_list,
    W_list = W_list,
    climate = list(vil = climate_vil, gra = climate_gra, wat = climate_wat),
    n_vil = ref$n_vil,
    n_gra = ref$n_gra,
    n_wat = ref$n_wat,
    vil_names = ref$vil_names,
    gra_names = ref$gra_names,
    wat_names = ref$wat_names,
    metadata = list(
      scenario_type = scenario_type,
      sequence_label = sequence_label %||% paste(scenario_ids, collapse = ";"),
      component_scenarios = scenario_ids,
      sequence_years = vapply(bundles, function(x) x$metadata$scenario_year, numeric(1)),
      sequence_length_months = length(G_list)
    )
  )

  if (validate) {
    validate_rvf_scenario(scenario)
  }

  scenario
}

rvf_complete_replicates_for_years <- function(index = list_rvf_scenarios(),
                                              dataset,
                                              years,
                                              climate_model = NULL,
                                              network_replicate = NULL,
                                              scenario = NULL) {
  out <- rvf_select_scenarios(
    index = index,
    dataset = dataset,
    climate_model = climate_model,
    year = years,
    network_replicate = network_replicate,
    scenario = scenario
  )
  if (nrow(out) == 0 || !"network_replicate" %in% names(out)) {
    return(integer(0))
  }
  out |>
    dplyr::group_by(.data$network_replicate) |>
    dplyr::summarise(n_years = dplyr::n_distinct(.data$year), .groups = "drop") |>
    dplyr::filter(.data$n_years == length(unique(years))) |>
    dplyr::pull(.data$network_replicate) |>
    sort()
}

## Replace each year's monthly networks with that year's mean commute matrix.
annualize_rvf_scenario_networks <- function(scenario, validate = TRUE) {
  step_names <- names(scenario$G_list)
  if (is.null(step_names)) {
    step_names <- paste0("step_", seq_along(scenario$G_list))
    names(scenario$G_list) <- step_names
    names(scenario$W_list) <- step_names
    names(scenario$climate$vil) <- step_names
    names(scenario$climate$gra) <- step_names
    names(scenario$climate$wat) <- step_names
  }

  year_tag <- sub("^[A-Za-z]{3}_?", "", step_names)
  year_group <- ifelse(grepl("^[0-9]{4}$", year_tag), year_tag, "all")

  new_G <- scenario$G_list
  new_W <- scenario$W_list

  for (grp in unique(year_group)) {
    idx <- which(year_group == grp)
    G_mean <- Reduce("+", scenario$G_list[idx]) / length(idx)
    W_mean <- Reduce("+", scenario$W_list[idx]) / length(idx)
    G_mean <- normalize_rows_scenario(G_mean)
    W_mean <- normalize_rows_scenario(W_mean)

    for (j in idx) {
      new_G[[j]] <- G_mean
      new_W[[j]] <- W_mean
    }
  }

  scenario$G_list <- new_G
  scenario$W_list <- new_W
  scenario$metadata$network_mode <- "annualized"

  if (validate) {
    validate_rvf_scenario(scenario)
  }

  scenario
}

## Apply either the raw monthly networks or the annualized alternative.
rvf_apply_network_mode <- function(scenario, network_mode = c("monthly", "annualized"), validate = TRUE) {
  network_mode <- match.arg(network_mode)

  if (network_mode == "monthly") {
    scenario$metadata$network_mode <- "monthly"
    if (validate) validate_rvf_scenario(scenario)
    return(scenario)
  }

  annualize_rvf_scenario_networks(scenario, validate = validate)
}

# --------------------------------------------------------------------------
# Calibration recommendation helpers.
#
# Script 02 writes a per-profile `calibrate_recommended.csv`. These helpers
# load the matching profile and rainfall-response calibration for downstream
# analyses, avoiding hand-copied capacity scales.
# --------------------------------------------------------------------------

rvf_calibration_output_dir <- function(cfg = load_rvf_config(), suffix = NULL) {
  suffix <- trimws(suffix %||% cfg$sim$calibration_source_suffix %||% "")
  subdir <- "calibrate_seroprevalence"
  if (nzchar(suffix)) {
    suffix <- gsub("[^A-Za-z0-9_.-]+", "_", suffix)
    subdir <- paste0(subdir, "_", suffix)
  }
  file.path(cfg$results_dir, subdir)
}

rvf_calibration_recommended_path <- function(cfg = load_rvf_config(), suffix = NULL) {
  file.path(rvf_calibration_output_dir(cfg, suffix), "calibrate_recommended.csv")
}

rvf_read_calibration_recommendations <- function(cfg = load_rvf_config(), suffix = NULL) {
  path <- rvf_calibration_recommended_path(cfg, suffix)
  if (!file.exists(path)) {
    stop(
      "No calibration recommendation table found at: ", path,
      "\nRun 02_calibrate_vector_capacity.R first."
    )
  }

  scenario_index_path <- file.path(cfg$scenarios_dir, "scenario_bundle_index.csv")
  if (file.exists(scenario_index_path) && file.info(path)$mtime < file.info(scenario_index_path)$mtime) {
    stop(
      "Calibration recommendations predate the current scenario bundles: ", path,
      "\nRerun 02_calibrate_vector_capacity.R before any downstream transmission analysis."
    )
  }
  read.csv(path, stringsAsFactors = FALSE)
}

# Keep the primary rainfall formulation in canonical output directories and
# tag structural sensitivities so they cannot overwrite manuscript inputs.
rvf_rainfall_output_subdir <- function(base_subdir,
                                       rainfall_model,
                                       cfg = load_rvf_config()) {
  primary_model <- cfg$primary_rainfall_model %||% cfg$sim$rainfall_model
  if (identical(rainfall_model, primary_model)) {
    return(base_subdir)
  }
  paste0(base_subdir, "_", gsub("[^A-Za-z0-9_.-]+", "_", rainfall_model))
}

rvf_profile_calibration_recommendation <- function(site_profile,
                                                   rainfall_model,
                                                   p_mobile = NULL,
                                                   cfg = load_rvf_config(),
                                                   suffix = NULL,
                                                   require_target_band = FALSE) {
  site_profile <- rvf_resolve_site_profile_id(site_profile, cfg = cfg)
  recs <- rvf_read_calibration_recommendations(cfg, suffix)
  if ("site_profile" %in% names(recs)) {
    recs$site_profile <- vapply(
      recs$site_profile,
      function(x) rvf_resolve_site_profile_id(x, cfg = cfg),
      character(1)
    )
  }
  recs <- recs[
    recs$site_profile == site_profile & recs$rainfall_model == rainfall_model,
    ,
    drop = FALSE
  ]
  if (nrow(recs) == 0L) {
    stop(
      "No calibration recommendation for site_profile='", site_profile,
      "' and rainfall_model='", rainfall_model, "'.",
      "\nAdd the profile and rainfall response to execution$calibration in R/config.R, then rerun script 02."
    )
  }

  if (!is.null(p_mobile)) {
    p_mobile <- suppressWarnings(as.numeric(p_mobile)[1])
    if (!is.finite(p_mobile)) {
      stop("p_mobile must be finite when selecting a profile calibration.")
    }
    if (!"p_mobile" %in% names(recs)) {
      stop(
        "Calibration recommendation table does not record p_mobile, so it predates the p_mobile=1 baseline update.",
        "\nRerun 02_calibrate_vector_capacity.R before downstream analyses."
      )
    }
    recs$p_mobile <- suppressWarnings(as.numeric(recs$p_mobile))
    recs <- recs[is.finite(recs$p_mobile) & abs(recs$p_mobile - p_mobile) < 1e-12, , drop = FALSE]
    if (nrow(recs) == 0L) {
      stop(
        "No calibration recommendation for site_profile='", site_profile,
        "', rainfall_model='", rainfall_model, "', and p_mobile=", signif(p_mobile, 4), ".",
        "\nRerun script 02 with the requested p_mobile baseline before applying profile calibration."
      )
    }
  }

  recs$target_distance <- suppressWarnings(as.numeric(recs$target_distance))
  recs$target_mid_distance <- suppressWarnings(as.numeric(recs$target_mid_distance))
  recs <- recs[order(recs$target_distance, recs$target_mid_distance), , drop = FALSE]
  rec <- recs[1, , drop = FALSE]
  rec_sero <- if ("mean_seroprevalence_calibration" %in% names(rec)) {
    as.numeric(rec$mean_seroprevalence_calibration[[1]])
  } else {
    as.numeric(rec$mean_seroprevalence_last5y[[1]])
  }

  in_band <- isTRUE(as.logical(rec$target_band[[1]]))
  if (require_target_band && !in_band) {
    stop(
      "Recommended calibration for site_profile='", site_profile,
      "' is outside the serology target band: calibration_seroprevalence=",
      signif(rec_sero, 4),
      "%. Use this as a rejected/stress-test profile, or rerun 06 with a wider biologically justified grid."
    )
  }

  rec
}

rvf_apply_profile_calibration <- function(sim,
                                          cfg = load_rvf_config(),
                                          suffix = NULL,
                                          require_target_band = FALSE,
                                          verbose = TRUE) {
  rec <- rvf_profile_calibration_recommendation(
    site_profile = sim$site_profile %||% "grazing_floodplain_dominant",
    rainfall_model = sim$rainfall_model %||% "quadratic_rainfall",
    p_mobile = sim$p_mobile %||% cfg$sim$p_mobile %||% 1,
    cfg = cfg,
    suffix = suffix,
    require_target_band = require_target_band
  )

  sim$site_profile <- rec$site_profile[[1]]
  sim$rainfall_model <- rec$rainfall_model[[1]]
  sim$capacity_scale <- as.numeric(rec$capacity_scale[[1]])
  sim$profile_calibration_source <- basename(rvf_calibration_recommended_path(cfg, suffix))
  sim$profile_calibration_target_band <- as.logical(rec$target_band[[1]])
  sim$profile_calibration_endemic_regime <- rec$endemic_regime[[1]]

  if (verbose) {
    rec_sero <- if ("mean_seroprevalence_calibration" %in% names(rec)) {
      as.numeric(rec$mean_seroprevalence_calibration[[1]])
    } else {
      as.numeric(rec$mean_seroprevalence_last5y[[1]])
    }
    message(sprintf(
      "Using profile calibration: profile=%s rainfall=%s capacity_scale=%.3g mean_sero=%.2f%% target_band=%s",
      sim$site_profile,
      sim$rainfall_model,
      sim$capacity_scale,
      rec_sero,
      as.character(sim$profile_calibration_target_band)
    ))
  }
  sim
}

# --------------------------------------------------------------------------
# Baseline-file resolution helpers.
#
# Script 03 writes observed-baseline outputs under a run tag that
# encodes whether transmission_scale / capacity_scale / rainfall_model /
# waning_years differ from the configured unscaled default
# (trans=1, cap=1, primary rainfall response, lifelong immunity). The helpers below
# reconstructs the same tag from cfg$sim so downstream scripts can find the
# right files without hard-coding a particular calibration.
# --------------------------------------------------------------------------

rvf_baseline_run_tag <- function(sim, cfg = load_rvf_config()) {
  if (!is.null(sim[["sim"]]) && is.list(sim[["sim"]])) {
    sim <- sim[["sim"]]
  }

  transmission_scale <- sim[["transmission_scale"]] %||% 1
  capacity_scale <- sim[["capacity_scale"]] %||% 1
  rainfall_model <- sim[["rainfall_model"]] %||% cfg$primary_rainfall_model
  site_profile <- sim[["site_profile"]] %||% "balanced"
  site_profile <- rvf_resolve_site_profile_id(site_profile)
  p_mobile <- sim[["p_mobile"]] %||% 1
  waning_years <- sim[["waning_years"]] %||% Inf
  baseline_years <- sim[["baseline_years"]] %||% sim[["baseline_year"]] %||% NA_integer_
  external_introduction_per_year <- sim[["external_introduction_per_year"]] %||% 0

  is_default_transmission <- abs(transmission_scale - 1) < 1e-12
  is_default_capacity <- abs(capacity_scale - 1) < 1e-12
  is_default_scale <- is_default_transmission && is_default_capacity
  is_default_rainfall <- identical(rainfall_model, cfg$primary_rainfall_model)
  is_default_waning <- rvf_is_default_waning(waning_years)

  parts <- "observed_cycle_endemic"
  if (!is_default_transmission) {
    parts <- c(parts, sprintf("t%.2f", transmission_scale))
  }
  if (!is_default_capacity) {
    parts <- c(parts, sprintf("k%.2f", capacity_scale))
  }
  if (!is_default_rainfall) {
    parts <- c(parts, paste0("rf_", rainfall_model))
  }
  if (!identical(site_profile, "balanced")) {
    parts <- c(parts, paste0("sp_", gsub("[^A-Za-z0-9]+", "_", site_profile)))
  }
  p_mobile <- suppressWarnings(as.numeric(p_mobile)[1])
  if (is.finite(p_mobile)) {
    parts <- c(parts, sprintf("pm%.2f", p_mobile))
  }
  if (length(baseline_years) > 1) {
    parts <- c(parts, sprintf("y%d-%d", min(baseline_years, na.rm = TRUE), max(baseline_years, na.rm = TRUE)))
  } else if (length(baseline_years) == 1 && is.finite(baseline_years)) {
    parts <- c(parts, sprintf("y%d", as.integer(baseline_years)))
  }
  if (!is_default_waning) {
    parts <- c(parts, paste0("wane", rvf_waning_label(waning_years)))
  }
  if (is.finite(external_introduction_per_year) && external_introduction_per_year > 0) {
    parts <- c(parts, paste0("intro", gsub("\\.", "p", as.character(external_introduction_per_year)), "py"))
  }
  paste(parts, collapse = "_")
}

rvf_baseline_subdir <- function(sim) {
  if (!is.null(sim[["sim"]]) && is.list(sim[["sim"]])) {
    sim <- sim[["sim"]]
  }

  if (identical(sim[["network_mode"]], "annualized")) {
    "observed_cycle_endemic_annualized"
  } else {
    "observed_cycle_endemic"
  }
}

# Compact baseline tables used by downstream summaries.
rvf_baseline_dir <- function(sim, cfg = load_rvf_config()) {
  file.path(cfg$results_dir, rvf_baseline_subdir(sim))
}

# Nested or expensive baseline objects that are not manuscript tables.
rvf_baseline_intermediate_dir <- function(sim, cfg = load_rvf_config()) {
  file.path(cfg$intermediate_dir, rvf_baseline_subdir(sim))
}

rvf_baseline_state_path <- function(sim, cfg = load_rvf_config()) {
  file.path(
    rvf_baseline_intermediate_dir(sim, cfg),
    paste0(rvf_baseline_run_tag(sim), "_baseline_state.rds")
  )
}

rvf_baseline_summary_path <- function(sim, cfg = load_rvf_config()) {
  file.path(
    rvf_baseline_dir(sim, cfg),
    paste0(rvf_baseline_run_tag(sim), "_summary.csv")
  )
}

rvf_baseline_summary_aggregate_path <- function(sim, cfg = load_rvf_config()) {
  file.path(
    rvf_baseline_dir(sim, cfg),
    paste0(rvf_baseline_run_tag(sim), "_summary_aggregate.csv")
  )
}

rvf_baseline_annual_system_path <- function(sim, cfg = load_rvf_config()) {
  file.path(
    rvf_baseline_dir(sim, cfg),
    paste0(rvf_baseline_run_tag(sim), "_annual_system.csv")
  )
}

rvf_load_baseline_summary <- function(sim, cfg = load_rvf_config()) {
  summary_path <- rvf_baseline_summary_path(sim, cfg)
  if (file.exists(summary_path)) {
    return(readr::read_csv(summary_path, show_col_types = FALSE))
  }

  # Backward-compatible fallback for archived runs saved before the lean-output
  # baseline summary became mandatory.
  results <- rvf_load_baseline_run(sim, cfg)
  dplyr::bind_rows(lapply(results, `[[`, "summary"))
}

rvf_baseline_total_hosts <- function(sim, cfg = load_rvf_config()) {
  if (!file.exists(cfg$mvt_inputs_file)) {
    stop("Module 02 observed inputs are missing: ", cfg$mvt_inputs_file)
  }
  inputs <- readRDS(cfg$mvt_inputs_file)
  n_villages <- length(inputs$village_names %||% character())
  hosts_per_village <- as.numeric(sim$host_population_per_village)
  if (n_villages < 1L || !is.finite(hosts_per_village) || hosts_per_village <= 0) {
    stop("Could not derive the total cattle population from current inputs.")
  }
  n_villages * hosts_per_village
}

# Pull a per-replicate vector of "cumulative incidence in the last 5 endemic
# years" from 05's compact summary. Used by future-scenario scripts to
# express cum-inc results as a percentage change against the current baseline.
# Returns one value per replicate so callers can compute their own quantiles.
rvf_baseline_cum_incidence_last5y <- function(sim, cfg = load_rvf_config()) {
  summary_df <- rvf_load_baseline_summary(sim, cfg)
  if ("cum_incidence_last5y" %in% names(summary_df)) {
    out <- summary_df$cum_incidence_last5y
    if ("network_replicate" %in% names(summary_df)) {
      names(out) <- as.character(summary_df$network_replicate)
    }
    return(out)
  }

  # Archived-output fallback. Current runs should not need the heavy run RDS.
  results <- rvf_load_baseline_run(sim, cfg)

  vapply(results, function(rep_res) {
    ts <- rep_res$time_series
    if (is.null(ts) || nrow(ts) == 0) return(NA_real_)

    endemic_years <- as.numeric(rep_res$summary$years_endemic[1])
    if (!is.finite(endemic_years)) endemic_years <- max(ts$year_since_seed, na.rm = TRUE)

    last5_window <- ts$year_since_seed >= (endemic_years - 5)
    if (!any(last5_window)) return(NA_real_)

    cum_end <- max(ts$cumulative_incidence, na.rm = TRUE)
    cum_start_of_last5 <- min(ts$cumulative_incidence[last5_window], na.rm = TRUE)
    cum_end - cum_start_of_last5
  }, numeric(1))
}

# Shared loader for the observed-baseline run object. Both helpers below use it; pulling
# it out keeps the file-existence error in one place and avoids two readRDS
# calls when an analysis script needs both cum-inc and sero baselines.
rvf_load_baseline_run <- function(sim, cfg = load_rvf_config()) {
  run_path <- file.path(
    rvf_baseline_intermediate_dir(sim, cfg),
    paste0(rvf_baseline_run_tag(sim), "_run.rds")
  )
  if (!file.exists(run_path)) {
    stop(
      "No observed-baseline run.rds at: ", run_path,
      "\nRun 03_run_observed_baseline.R first to produce the current baseline.",
      "\nFor lean-output runs, rerun script 03 with save_baseline_run_rds = TRUE ",
      "only when this optional diagnostic object is required."
    )
  }
  readRDS(run_path)
}

# Per-replicate baseline seroprevalence trajectory (medians/min/max across
# replicates at each output time). Returns a tibble with one row per
# year_since_seed step, ready to overlay as a reference curve on 08's
# future-scenario sero plots. Sero values are in percent.
rvf_baseline_seroprevalence_timeseries <- function(sim, cfg = load_rvf_config()) {
  results <- rvf_load_baseline_run(sim, cfg)

  per_rep <- lapply(seq_along(results), function(i) {
    ts <- results[[i]]$time_series
    if (is.null(ts) || nrow(ts) == 0) return(NULL)
    tibble::tibble(
      rep_index = i,
      year_since_seed = ts$year_since_seed,
      seroprevalence = ts$seroprevalence
    )
  })
  per_rep <- per_rep[!vapply(per_rep, is.null, logical(1))]
  if (length(per_rep) == 0) {
    stop("Baseline run.rds had no usable replicate timeseries.")
  }

  combined <- dplyr::bind_rows(per_rep)
  combined |>
    dplyr::group_by(year_since_seed) |>
    dplyr::summarise(
      sero_mid = stats::median(seroprevalence, na.rm = TRUE),
      sero_lo = min(seroprevalence, na.rm = TRUE),
      sero_hi = max(seroprevalence, na.rm = TRUE),
      n_reps = dplyr::n(),
      .groups = "drop"
    )
}

# Single scalar summary for adding a horizontal reference line: the baseline
# median seroprevalence across the last 5 endemic years, averaged across
# replicates first then taking the rep-level median.
rvf_baseline_seroprevalence_last5y <- function(sim, cfg = load_rvf_config()) {
  summary_df <- rvf_load_baseline_summary(sim, cfg)
  if ("mean_seroprevalence_last5y" %in% names(summary_df)) {
    per_rep <- summary_df$mean_seroprevalence_last5y
    return(list(
      per_rep = per_rep,
      median = stats::median(per_rep, na.rm = TRUE),
      lo = min(per_rep, na.rm = TRUE),
      hi = max(per_rep, na.rm = TRUE)
    ))
  }

  results <- rvf_load_baseline_run(sim, cfg)

  per_rep <- vapply(results, function(rep_res) {
    ts <- rep_res$time_series
    if (is.null(ts) || nrow(ts) == 0) return(NA_real_)
    endemic_years <- as.numeric(rep_res$summary$years_endemic[1])
    if (!is.finite(endemic_years)) endemic_years <- max(ts$year_since_seed, na.rm = TRUE)
    late <- ts$year_since_seed >= (endemic_years - 5)
    if (!any(late)) return(NA_real_)
    mean(ts$seroprevalence[late], na.rm = TRUE)
  }, numeric(1))

  list(
    per_rep = per_rep,
    median = stats::median(per_rep, na.rm = TRUE),
    lo = min(per_rep, na.rm = TRUE),
    hi = max(per_rep, na.rm = TRUE)
  )
}
