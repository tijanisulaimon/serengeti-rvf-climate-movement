#!/usr/bin/env Rscript

#' 09_summarise_spatial_seasonal_patterns.R
#' -------------------------------------------------------------------------
#' Test whether village-level endemic burden becomes less spatially
#' concentrated under future climate and movement windows. The comparison
#' uses the p_mobile = 1 observed scenario from script 04 and the matched
#' future-window results from script 06.
#'
#' Two calibrated exposure profiles are retained: the primary
#' grazing/floodplain profile and the contrasting peridomestic profile. Gini
#' coefficients quantify concentration without imposing an arbitrary burden
#' threshold. A degree--burden diagnostic tests whether villages that share
#' resources with more other villages also carry greater endemic burden. Maps
#' show the ensemble-median annual burden for the primary profile only; the
#' concentration and degree-association figures provide cross-profile checks.
#'
#' No transmission simulations are run here.
#' -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(sf)
  library(tidyr)
})

source(here::here("analysis", "03_rvf_transmission", "R", "transmission_model.R"))
source(here::here("analysis", "03_rvf_transmission", "R", "plot_helpers.R"))

cfg <- load_rvf_config()
sim <- cfg$sim

if (!identical(sim$network_mode, "monthly")) {
  stop("The publication spatial summary requires monthly movement networks.")
}

primary_profile <- rvf_resolve_site_profile_id(
  sim$future_site_profile %||% "grazing_floodplain_dominant"
)
profiles <- unique(c(primary_profile, "peridomestic_dominant"))
profile_plot_labels <- c(
  grazing_floodplain_dominant = "Shared grazing/floodplain exposure",
  peridomestic_dominant = "Village-associated exposure"
)

out_dir <- file.path(cfg$results_dir, "spatial_patterns")
figures_dir <- cfg$figures_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Function: read the current and future village summaries for one profile.
# The current subset is restricted to p_mobile = 1 so its residence-time
# structure matches the future-window experiment.
read_profile_spatial_outputs <- function(profile) {
  mobility_dir <- file.path(
    cfg$intermediate_dir,
    paste0(
      "mobility_capacity_cold_start_",
      profile,
      "_pmobile_capacity"
    )
  )

  future_name <- if (identical(profile, primary_profile)) {
    "future_cold_endemic"
  } else {
    paste0("future_cold_endemic_profile_", profile)
  }
  future_name <- rvf_rainfall_output_subdir(
    future_name,
    sim$rainfall_model,
    cfg
  )
  future_dir <- file.path(cfg$intermediate_dir, future_name)

  mobility_path <- file.path(mobility_dir, "mobility_capacity_runs.rds")
  future_path <- file.path(
    future_dir,
    "future_cold_endemic_village_terminal.rds"
  )
  missing_paths <- c(mobility_path, future_path)[
    !file.exists(c(mobility_path, future_path))
  ]
  if (length(missing_paths) > 0L) {
    stop(
      "Missing spatial input(s) for ", profile, ":\n  ",
      paste(missing_paths, collapse = "\n  "),
      "\nRun scripts 04 and 06 for this profile before script 09."
    )
  }

  mobility <- readRDS(mobility_path)
  if (is.null(mobility$village_terminal)) {
    stop("Script-04 output lacks village_terminal: ", mobility_path)
  }

  current <- mobility$village_terminal |>
    filter(abs(.data$p_mobile - 1) < 1e-12) |>
    transmute(
      site_profile = .data$site_profile,
      profile_label = unname(profile_plot_labels[.data$site_profile]),
      scenario = "current",
      scenario_label = "Current",
      future_window = "Current",
      future_year = NA_integer_,
      window_start = NA_integer_,
      window_end = NA_integer_,
      climate_model = "Observed",
      .data$network_replicate,
      .data$village_id,
      .data$annual_new_infections
    )

  future <- readRDS(future_path) |>
    filter(
      .data$scenario %in% cfg$future_scenarios,
      .data$climate_model %in% cfg$climate_models
    ) |>
    transmute(
      site_profile = .data$site_profile,
      profile_label = unname(profile_plot_labels[.data$site_profile]),
      .data$scenario,
      .data$scenario_label,
      .data$future_window,
      .data$future_year,
      .data$window_start,
      .data$window_end,
      .data$climate_model,
      .data$network_replicate,
      .data$village_id,
      .data$annual_new_infections
    )

  if (!identical(unique(current$site_profile), profile) ||
      !identical(unique(future$site_profile), profile)) {
    stop("Profile metadata do not match the requested profile: ", profile)
  }

  bind_rows(current, future)
}

village_burden <- bind_rows(lapply(profiles, read_profile_spatial_outputs))

# Every run should contain the same 94 home villages. This gate catches
# incomplete joins or interrupted script-06 outputs before spatial summaries.
run_counts <- village_burden |>
  count(
    .data$site_profile,
    .data$scenario,
    .data$future_window,
    .data$climate_model,
    .data$network_replicate,
    name = "n_villages"
  )
if (any(run_counts$n_villages != 94L)) {
  stop("At least one spatial run does not contain all 94 villages.")
}
if (any(!is.finite(village_burden$annual_new_infections)) ||
    any(village_burden$annual_new_infections < 0)) {
  stop("Village infection burdens must be finite and non-negative.")
}

# Require the complete publication grid. Counting villages alone would allow
# an older 10-replicate result to pass after the configuration moved to 20.
observed_spatial_runs <- village_burden |>
  filter(.data$scenario == "current") |>
  distinct(.data$site_profile, .data$network_replicate)
expected_observed_spatial_runs <- crossing(
  site_profile = profiles,
  network_replicate = sort(unique(cfg$analysis_network_replicates))
)
future_spatial_runs <- village_burden |>
  filter(.data$scenario %in% cfg$future_scenarios) |>
  distinct(
    .data$site_profile,
    .data$scenario,
    .data$future_window,
    .data$climate_model,
    .data$network_replicate
  )
expected_future_spatial_runs <- crossing(
  site_profile = profiles,
  scenario = cfg$future_scenarios,
  future_window = names(cfg$future_windows),
  climate_model = cfg$climate_models,
  network_replicate = sort(unique(cfg$analysis_network_replicates))
)
if (!setequal(
  observed_spatial_runs,
  expected_observed_spatial_runs
) || !setequal(
  future_spatial_runs,
  expected_future_spatial_runs
)) {
  stop(
    "Spatial summaries contain an incomplete, duplicated, or stale ",
    "profile-by-scenario-by-replicate grid."
  )
}

# -------- Village network degree ----------------------------------------
# A village-level R0 would be the natural quantity to relate to network degree,
# but none is available for this periodically forced model, in which villages
# share vector sites. We instead test the directly observable relationship
# between one-mode degree and annual infections. This is a descriptive
# association, not an estimate of movement's causal effect.
scenario_index_path <- file.path(
  cfg$scenarios_dir,
  "scenario_bundle_index.csv"
)
if (!file.exists(scenario_index_path)) {
  stop("Scenario bundle index not found: ", scenario_index_path)
}

window_range_to_id <- setNames(
  names(cfg$future_windows),
  vapply(
    cfg$future_windows,
    function(years) sprintf("%d-%d", min(years), max(years)),
    character(1)
  )
)

scenario_index <- read_csv(scenario_index_path, show_col_types = FALSE) |>
  filter(
    .data$network_replicate %in% cfg$analysis_network_replicates,
    (.data$dataset == "observed" & .data$year %in% sim$baseline_years) |
      (
        .data$dataset == "future" &
          .data$scenario %in% cfg$future_scenarios &
          .data$climate_model %in% cfg$climate_models &
          .data$year %in% unlist(cfg$future_windows, use.names = FALSE)
      )
  ) |>
  arrange(
    .data$dataset,
    .data$scenario,
    .data$climate_model,
    .data$year,
    .data$network_replicate
  ) |>
  mutate(
    # Align Module 02's calendar-range labels with Module 03's window IDs.
    future_window = if_else(
      .data$dataset == "observed",
      "Current",
      unname(window_range_to_id[.data$future_window])
    )
  )

if (anyNA(scenario_index$future_window)) {
  stop("At least one network scenario could not be assigned to a forcing window.")
}

missing_scenarios <- scenario_index$file[!file.exists(scenario_index$file)]
if (length(missing_scenarios) > 0L) {
  stop(
    "Processed scenario bundle(s) are missing, for example: ",
    missing_scenarios[[1]]
  )
}

# Function: project one monthly bipartite matrix to village adjacency.
# A binary one-mode edge joins two villages when both use at least one of the
# same resource sites. Matrix values are row-normalised use weights, so only
# their positive/non-positive status is used for this structural summary.
one_mode_adjacency <- function(bipartite_matrix) {
  adjacency <- tcrossprod(bipartite_matrix > 0) > 0
  diag(adjacency) <- FALSE
  adjacency
}

# Function: calculate village network and climate metrics for one year/draw.
# Combined degree is the union of grazing- and watering-mediated connections,
# not the sum, so a village pair sharing both resource types is counted once.
residence_time <- rvf_time_budget_profile(no_movement = FALSE)$routine
if (!isTRUE(all.equal(sum(residence_time), 1))) {
  stop("The residence-time fractions must sum to one.")
}

read_scenario_spatial_metrics <- function(index_row) {
  scenario <- readRDS(index_row$file[[1]])
  if (length(scenario$G_list) != 12L || length(scenario$W_list) != 12L) {
    stop("Scenario does not contain 12 monthly networks: ", index_row$file[[1]])
  }

  monthly_degree <- bind_rows(lapply(seq_along(scenario$G_list), function(month) {
    grazing_matrix <- scenario$G_list[[month]]
    watering_matrix <- scenario$W_list[[month]]
    grazing_adjacency <- one_mode_adjacency(grazing_matrix)
    watering_adjacency <- one_mode_adjacency(watering_matrix)

    home_temperature <- scenario$climate$vil[[month]]$tas
    grazing_temperature <- scenario$climate$gra[[month]]$tas
    watering_temperature <- scenario$climate$wat[[month]]$tas
    if (any(!is.finite(c(
      home_temperature,
      grazing_temperature,
      watering_temperature
    )))) {
      stop("Non-finite temperature in scenario: ", index_row$file[[1]])
    }

    # Match the ODE treatment of missing boundary rainfall as a dry month.
    home_rainfall <- scenario$climate$vil[[month]]$pr
    grazing_rainfall <- scenario$climate$gra[[month]]$pr
    watering_rainfall <- scenario$climate$wat[[month]]$pr
    home_rainfall[!is.finite(home_rainfall)] <- 0
    grazing_rainfall[!is.finite(grazing_rainfall)] <- 0
    watering_rainfall[!is.finite(watering_rainfall)] <- 0

    # Row-normalised matrices give each village's mean climate at the resource
    # sites it uses. The 15/8/1 weights then provide a simple daily
    # residence-time climate summary; they do not represent a fitted climate
    # effect or replace the model's non-linear monthly trait functions.
    used_grazing_temperature <- as.vector(
      grazing_matrix %*% grazing_temperature
    )
    used_watering_temperature <- as.vector(
      watering_matrix %*% watering_temperature
    )
    used_grazing_rainfall <- as.vector(grazing_matrix %*% grazing_rainfall)
    used_watering_rainfall <- as.vector(watering_matrix %*% watering_rainfall)

    tibble(
      village_id = seq_len(scenario$n_vil),
      grazing_degree = rowSums(grazing_adjacency),
      watering_degree = rowSums(watering_adjacency),
      combined_degree = rowSums(grazing_adjacency | watering_adjacency),
      home_temperature = home_temperature,
      home_rainfall = home_rainfall,
      residence_temperature =
        residence_time[["t_v"]] * home_temperature +
        residence_time[["t_g"]] * used_grazing_temperature +
        residence_time[["t_w"]] * used_watering_temperature,
      residence_rainfall =
        residence_time[["t_v"]] * home_rainfall +
        residence_time[["t_g"]] * used_grazing_rainfall +
        residence_time[["t_w"]] * used_watering_rainfall
    )
  }))

  monthly_degree |>
    group_by(.data$village_id) |>
    summarise(
      # Mean number of other villages sharing a resource in a month.
      grazing_degree = mean(.data$grazing_degree),
      watering_degree = mean(.data$watering_degree),
      combined_degree = mean(.data$combined_degree),
      home_temperature = mean(.data$home_temperature),
      home_rainfall = mean(.data$home_rainfall),
      residence_temperature = mean(.data$residence_temperature),
      residence_rainfall = mean(.data$residence_rainfall),
      .groups = "drop"
    ) |>
    mutate(
      dataset = index_row$dataset[[1]],
      scenario = if_else(
        index_row$dataset[[1]] == "observed",
        "current",
        index_row$scenario[[1]]
      ),
      future_window = index_row$future_window[[1]],
      climate_model = if_else(
        index_row$dataset[[1]] == "observed",
        "Observed",
        index_row$climate_model[[1]]
      ),
      year = index_row$year[[1]],
      network_replicate = index_row$network_replicate[[1]],
      .before = 1
    )
}

# Each row below is already averaged over the 12 months of one forcing year.
# Average again over the complete 10-year observed cycle or matched 20-year
# future window because the transmission burden uses those repeated cycles.
spatial_metrics_by_year <- bind_rows(lapply(
  seq_len(nrow(scenario_index)),
  function(i) read_scenario_spatial_metrics(scenario_index[i, ])
))

village_spatial_metrics <- spatial_metrics_by_year |>
  group_by(
    .data$scenario,
    .data$future_window,
    .data$climate_model,
    .data$network_replicate,
    .data$village_id
  ) |>
  summarise(
    forcing_years = n_distinct(.data$year),
    grazing_degree = mean(.data$grazing_degree),
    watering_degree = mean(.data$watering_degree),
    combined_degree = mean(.data$combined_degree),
    home_temperature = mean(.data$home_temperature),
    home_rainfall = mean(.data$home_rainfall),
    residence_temperature = mean(.data$residence_temperature),
    residence_rainfall = mean(.data$residence_rainfall),
    .groups = "drop"
  )

expected_years <- if_else(
  village_spatial_metrics$scenario == "current",
  length(sim$baseline_years),
  unique(vapply(cfg$future_windows, length, integer(1)))
)
if (any(village_spatial_metrics$forcing_years != expected_years)) {
  stop("At least one village-degree summary has an incomplete forcing cycle.")
}

degree_join_keys <- c(
  "scenario",
  "future_window",
  "climate_model",
  "network_replicate",
  "village_id"
)
village_burden_metrics <- village_burden |>
  left_join(village_spatial_metrics, by = degree_join_keys)
if (anyNA(village_burden_metrics$combined_degree)) {
  stop("Village burden could not be matched to all network-degree summaries.")
}

# Spearman correlation is scale-free and allows a monotonic, non-linear
# relationship. Separate grazing and watering estimates are kept for
# interpretation, and the figure uses their one-mode union as a single
# combined degree.
degree_association_by_run <- village_burden_metrics |>
  group_by(
    .data$site_profile,
    .data$profile_label,
    .data$scenario,
    .data$scenario_label,
    .data$future_window,
    .data$future_year,
    .data$window_start,
    .data$window_end,
    .data$climate_model,
    .data$network_replicate
  ) |>
  summarise(
    n_villages = n(),
    grazing_degree_rho = cor(
      .data$grazing_degree,
      .data$annual_new_infections,
      method = "spearman"
    ),
    watering_degree_rho = cor(
      .data$watering_degree,
      .data$annual_new_infections,
      method = "spearman"
    ),
    combined_degree_rho = cor(
      .data$combined_degree,
      .data$annual_new_infections,
      method = "spearman"
    ),
    home_temperature_rho = cor(
      .data$home_temperature,
      .data$annual_new_infections,
      method = "spearman"
    ),
    home_rainfall_rho = cor(
      .data$home_rainfall,
      .data$annual_new_infections,
      method = "spearman"
    ),
    residence_temperature_rho = cor(
      .data$residence_temperature,
      .data$annual_new_infections,
      method = "spearman"
    ),
    residence_rainfall_rho = cor(
      .data$residence_rainfall,
      .data$annual_new_infections,
      method = "spearman"
    ),
    .groups = "drop"
  )

association_columns <- c(
  "combined_degree_rho",
  "home_temperature_rho",
  "home_rainfall_rho",
  "residence_temperature_rho",
  "residence_rainfall_rho"
)
if (any(!is.finite(as.matrix(degree_association_by_run[association_columns])))) {
  stop("A village-level burden association is undefined for at least one run.")
}

degree_association_summary <- degree_association_by_run |>
  group_by(
    .data$site_profile,
    .data$profile_label,
    .data$scenario,
    .data$scenario_label,
    .data$future_window,
    .data$future_year,
    .data$window_start,
    .data$window_end
  ) |>
  summarise(
    n_runs = n(),
    combined_rho_mid = median(.data$combined_degree_rho),
    combined_rho_lo = as.numeric(
      quantile(.data$combined_degree_rho, 0.10)
    ),
    combined_rho_hi = as.numeric(
      quantile(.data$combined_degree_rho, 0.90)
    ),
    residence_temperature_rho_mid = median(.data$residence_temperature_rho),
    residence_temperature_rho_lo = as.numeric(
      quantile(.data$residence_temperature_rho, 0.10)
    ),
    residence_temperature_rho_hi = as.numeric(
      quantile(.data$residence_temperature_rho, 0.90)
    ),
    residence_rainfall_rho_mid = median(.data$residence_rainfall_rho),
    residence_rainfall_rho_lo = as.numeric(
      quantile(.data$residence_rainfall_rho, 0.10)
    ),
    residence_rainfall_rho_hi = as.numeric(
      quantile(.data$residence_rainfall_rho, 0.90)
    ),
    .groups = "drop"
  )

# Function: Gini coefficient of non-negative village burdens.
# Zero denotes equal burden across villages; values approaching one indicate
# concentration in a small number of villages.
gini_coefficient <- function(x) {
  x <- sort(x[is.finite(x) & x >= 0])
  n <- length(x)
  if (n < 2L || sum(x) <= 0) {
    return(NA_real_)
  }
  (2 * sum(seq_len(n) * x) - (n + 1) * sum(x)) / (n * sum(x))
}

# Calculate concentration within each GCM x movement-network realisation.
# The effective-village count is the inverse Simpson concentration of burden:
# 1 means one village carries all burden and 94 means perfectly even burden.
spatial_by_run <- village_burden |>
  group_by(
    .data$site_profile,
    .data$profile_label,
    .data$scenario,
    .data$scenario_label,
    .data$future_window,
    .data$future_year,
    .data$window_start,
    .data$window_end,
    .data$climate_model,
    .data$network_replicate
  ) |>
  summarise(
    n_villages = n(),
    total_annual_infections = sum(.data$annual_new_infections),
    gini = gini_coefficient(.data$annual_new_infections),
    top_decile_burden_share = {
      ordered_burden <- sort(
        .data$annual_new_infections,
        decreasing = TRUE
      )
      sum(head(ordered_burden, ceiling(n() * 0.10))) / sum(ordered_burden)
    },
    effective_villages = {
      burden <- .data$annual_new_infections
      sum(burden)^2 / sum(burden^2)
    },
    .groups = "drop"
  )

# P10-P90 intervals represent variation across network replicates for the
# current baseline and across GCM x network-replicate combinations in future.
spatial_summary <- spatial_by_run |>
  group_by(
    .data$site_profile,
    .data$profile_label,
    .data$scenario,
    .data$scenario_label,
    .data$future_window,
    .data$future_year,
    .data$window_start,
    .data$window_end
  ) |>
  summarise(
    n_runs = n(),
    gini_mid = median(.data$gini),
    gini_lo = as.numeric(quantile(.data$gini, 0.10)),
    gini_hi = as.numeric(quantile(.data$gini, 0.90)),
    top_decile_share_mid = median(.data$top_decile_burden_share),
    top_decile_share_lo = as.numeric(
      quantile(.data$top_decile_burden_share, 0.10)
    ),
    top_decile_share_hi = as.numeric(
      quantile(.data$top_decile_burden_share, 0.90)
    ),
    effective_villages_mid = median(.data$effective_villages),
    effective_villages_lo = as.numeric(
      quantile(.data$effective_villages, 0.10)
    ),
    effective_villages_hi = as.numeric(
      quantile(.data$effective_villages, 0.90)
    ),
    .groups = "drop"
  ) |>
  arrange(.data$site_profile, .data$scenario, .data$future_year)

# Express each village's burden relative to the total within the same model
# run. These shares isolate spatial concentration from the much larger future
# totals and therefore provide the appropriate quantity for comparative maps.
village_burden <- village_burden |>
  group_by(
    .data$site_profile,
    .data$scenario,
    .data$future_window,
    .data$climate_model,
    .data$network_replicate
  ) |>
  mutate(
    total_annual_infections = sum(.data$annual_new_infections),
    annual_infection_share =
      .data$annual_new_infections / .data$total_annual_infections
  ) |>
  ungroup()

if (any(!is.finite(village_burden$annual_infection_share)) ||
    any(village_burden$annual_infection_share < 0)) {
  stop("Village infection shares must be finite and non-negative.")
}

# Median village burdens and shares support the maps and provide a compact
# numerical record without duplicating the full simulation output.
village_summary <- village_burden |>
  group_by(
    .data$site_profile,
    .data$profile_label,
    .data$scenario,
    .data$scenario_label,
    .data$future_window,
    .data$future_year,
    .data$window_start,
    .data$window_end,
    .data$village_id
  ) |>
  summarise(
    n_runs = n(),
    annual_infections_mid = median(.data$annual_new_infections),
    annual_infections_lo = as.numeric(
      quantile(.data$annual_new_infections, 0.10)
    ),
    annual_infections_hi = as.numeric(
      quantile(.data$annual_new_infections, 0.90)
    ),
    annual_share_mid = median(.data$annual_infection_share),
    annual_share_lo = as.numeric(
      quantile(.data$annual_infection_share, 0.10)
    ),
    annual_share_hi = as.numeric(
      quantile(.data$annual_infection_share, 0.90)
    ),
    .groups = "drop"
  )

write_csv(
  spatial_by_run,
  file.path(out_dir, "spatial_concentration_by_run.csv")
)
write_csv(
  spatial_summary,
  file.path(out_dir, "spatial_concentration_summary.csv")
)
write_csv(
  village_summary,
  file.path(out_dir, "village_burden_summary.csv")
)
write_csv(
  degree_association_by_run,
  file.path(out_dir, "village_burden_correlations_by_run.csv")
)
write_csv(
  degree_association_summary,
  file.path(out_dir, "village_burden_correlation_summary.csv")
)

# -------- Spatial-concentration figure -----------------------------------
# The inverse-Simpson effective number is easier to interpret than Gini in the
# main narrative: it is the equivalent number of villages contributing evenly
# to the total annual burden, bounded above by the 94 modelled villages. Future
# intervals pool GCM and network-replicate uncertainty (10th--90th percentiles).
future_extent <- spatial_summary |>
  filter(.data$scenario != "current") |>
  mutate(time_label = sprintf("%d-%d", .data$window_start, .data$window_end))

current_extent <- spatial_summary |>
  filter(.data$scenario == "current") |>
  select(
    site_profile,
    profile_label,
    effective_villages_mid,
    effective_villages_lo,
    effective_villages_hi,
    gini_mid,
    gini_lo,
    gini_hi
  )

current_plot <- crossing(
  site_profile = profiles,
  scenario = cfg$future_scenarios
) |>
  left_join(current_extent, by = "site_profile") |>
  mutate(
    scenario_label = unname(cfg$scenario_labels[.data$scenario]),
    time_label = "Current"
  )

time_levels <- c(
  "Current",
  vapply(
    cfg$future_windows,
    function(years) sprintf("%d-%d", min(years), max(years)),
    character(1)
  )
)
concentration_plot_data <- bind_rows(current_plot, future_extent) |>
  mutate(
    profile_label = factor(
      .data$profile_label,
      levels = unname(profile_plot_labels[profiles])
    ),
    scenario_label = factor(
      .data$scenario_label,
      levels = unname(cfg$scenario_labels[cfg$future_scenarios])
    ),
    time_label = factor(.data$time_label, levels = time_levels)
  )

p_concentration <- ggplot(
  concentration_plot_data,
  aes(
    x = .data$time_label,
    y = .data$effective_villages_mid,
    group = 1
  )
) +
  geom_line(colour = "#0072B2", linewidth = 0.75) +
  geom_linerange(
    aes(
      ymin = .data$effective_villages_lo,
      ymax = .data$effective_villages_hi
    ),
    colour = "#0072B2",
    linewidth = 0.65
  ) +
  geom_point(colour = "#0072B2", size = 2.5) +
  facet_grid(profile_label ~ scenario_label) +
  scale_y_continuous(
    limits = c(0, length(unique(village_burden$village_id))),
    breaks = seq(0, 90, 15),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = NULL,
    y = "Effective number of villages\ncontributing to annual infections"
  ) +
  theme_rvf(base_size = 12.5) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black", size = 10.5),
    axis.text.x = element_text(angle = 20, hjust = 1),
    axis.title = element_text(colour = "black", size = 12.5),
    strip.text = element_text(colour = "black", size = 11, face = "bold"),
    panel.grid.minor = element_blank()
  )

rvf_ggsave(
  file.path(figures_dir, "06_spatial_concentration.png"),
  p_concentration,
  width = 11.2,
  height = 7.2
)

# -------- Network degree--burden association figure ---------------------
# Plot the run-level Spearman association rather than pooling villages across
# GCMs or network draws. This preserves the matched simulation unit and shows
# whether the association is consistent across network uncertainty.
future_degree_gcm <- degree_association_by_run |>
  filter(.data$scenario != "current") |>
  group_by(
    .data$site_profile,
    .data$profile_label,
    .data$scenario,
    .data$scenario_label,
    .data$future_window,
    .data$future_year,
    .data$window_start,
    .data$window_end,
    .data$climate_model
  ) |>
  summarise(
    rho_mid = median(.data$combined_degree_rho),
    rho_lo = as.numeric(quantile(.data$combined_degree_rho, 0.10)),
    rho_hi = as.numeric(quantile(.data$combined_degree_rho, 0.90)),
    .groups = "drop"
  ) |>
  mutate(
    time_label = sprintf("%d-%d", .data$window_start, .data$window_end)
  )

current_degree <- degree_association_by_run |>
  filter(.data$scenario == "current") |>
  group_by(.data$site_profile, .data$profile_label) |>
  summarise(
    rho_mid = median(.data$combined_degree_rho),
    rho_lo = as.numeric(quantile(.data$combined_degree_rho, 0.10)),
    rho_hi = as.numeric(quantile(.data$combined_degree_rho, 0.90)),
    .groups = "drop"
  )

current_degree_plot <- crossing(
  site_profile = profiles,
  scenario = cfg$future_scenarios,
  climate_model = cfg$climate_models
) |>
  left_join(current_degree, by = "site_profile") |>
  mutate(
    scenario_label = unname(cfg$scenario_labels[.data$scenario]),
    time_label = "Current"
  )

degree_plot_data <- bind_rows(current_degree_plot, future_degree_gcm) |>
  mutate(
    profile_label = factor(
      .data$profile_label,
      levels = unname(profile_plot_labels[profiles])
    ),
    scenario_label = factor(
      .data$scenario_label,
      levels = unname(cfg$scenario_labels[cfg$future_scenarios])
    ),
    time_label = factor(.data$time_label, levels = time_levels)
  )

p_degree_association <- ggplot(
  degree_plot_data,
  aes(
    x = .data$time_label,
    y = .data$rho_mid,
    colour = .data$climate_model,
    group = .data$climate_model
  )
) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.45) +
  geom_line(linewidth = 0.65) +
  geom_linerange(
    aes(ymin = .data$rho_lo, ymax = .data$rho_hi),
    linewidth = 0.55
  ) +
  geom_point(size = 2.2) +
  facet_grid(profile_label ~ scenario_label) +
  scale_colour_gcm() +
  scale_y_continuous(
    limits = c(-0.75, 0.40),
    breaks = seq(-0.75, 0.25, 0.25)
  ) +
  labs(
    x = NULL,
    y = "Spearman correlation between\nmean village degree and annual infections"
  ) +
  theme_rvf(base_size = 12.5) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black", size = 10.5),
    axis.text.x = element_text(angle = 20, hjust = 1),
    axis.title = element_text(colour = "black", size = 12.5),
    strip.text = element_text(colour = "black", size = 11, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 10.5),
    panel.grid.minor = element_blank()
  )

rvf_ggsave(
  file.path(figures_dir, "06_network_degree_burden_association.png"),
  p_degree_association,
  width = 11.2,
  height = 7.2
)

# -------- Spatial-driver sanity diagnostic ------------------------------
# Raw climate means are descriptive checks only. The transmission model uses
# monthly non-linear rainfall and temperature responses, so these correlations
# should not be interpreted as fitted climate effects or causal attribution.
driver_labels <- c(
  combined_degree_rho = "Mean network degree",
  residence_temperature_rho = "Mean temperature",
  residence_rainfall_rho = "Mean rainfall"
)
driver_colours <- c(
  "Mean network degree" = "#6A3D9A",
  "Mean temperature" = "#D55E00",
  "Mean rainfall" = "#0072B2"
)

driver_by_run <- degree_association_by_run |>
  select(
    "site_profile",
    "profile_label",
    "scenario",
    "scenario_label",
    "future_window",
    "window_start",
    "window_end",
    all_of(names(driver_labels))
  ) |>
  pivot_longer(
    cols = all_of(names(driver_labels)),
    names_to = "driver",
    values_to = "rho"
  ) |>
  mutate(driver = unname(driver_labels[.data$driver]))

future_driver <- driver_by_run |>
  filter(.data$scenario != "current") |>
  group_by(
    .data$site_profile,
    .data$profile_label,
    .data$scenario,
    .data$scenario_label,
    .data$future_window,
    .data$window_start,
    .data$window_end,
    .data$driver
  ) |>
  summarise(
    rho_mid = median(.data$rho),
    rho_lo = as.numeric(quantile(.data$rho, 0.10)),
    rho_hi = as.numeric(quantile(.data$rho, 0.90)),
    .groups = "drop"
  ) |>
  mutate(time_label = sprintf("%d-%d", .data$window_start, .data$window_end))

current_driver <- driver_by_run |>
  filter(.data$scenario == "current") |>
  group_by(.data$site_profile, .data$profile_label, .data$driver) |>
  summarise(
    rho_mid = median(.data$rho),
    rho_lo = as.numeric(quantile(.data$rho, 0.10)),
    rho_hi = as.numeric(quantile(.data$rho, 0.90)),
    .groups = "drop"
  )

current_driver_plot <- crossing(
  current_driver,
  scenario = cfg$future_scenarios
) |>
  mutate(
    scenario_label = unname(cfg$scenario_labels[.data$scenario]),
    time_label = "Current"
  )

driver_plot_data <- bind_rows(current_driver_plot, future_driver) |>
  mutate(
    profile_label = factor(
      .data$profile_label,
      levels = unname(profile_plot_labels[profiles])
    ),
    scenario_label = factor(
      .data$scenario_label,
      levels = unname(cfg$scenario_labels[cfg$future_scenarios])
    ),
    time_label = factor(.data$time_label, levels = time_levels),
    driver = factor(.data$driver, levels = unname(driver_labels))
  )

p_driver_diagnostic <- ggplot(
  driver_plot_data,
  aes(
    x = .data$time_label,
    y = .data$rho_mid,
    colour = .data$driver,
    group = .data$driver
  )
) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.45) +
  geom_line(linewidth = 0.65) +
  geom_linerange(
    aes(ymin = .data$rho_lo, ymax = .data$rho_hi),
    linewidth = 0.55
  ) +
  geom_point(size = 2.2) +
  facet_grid(profile_label ~ scenario_label) +
  scale_colour_manual(values = driver_colours) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  labs(
    x = NULL,
    y = "Spearman correlation with annual infections",
    colour = NULL
  ) +
  theme_rvf(base_size = 12.5) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black", size = 10.5),
    axis.text.x = element_text(angle = 20, hjust = 1),
    axis.title = element_text(colour = "black", size = 12.5),
    strip.text = element_text(colour = "black", size = 11, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 10.5),
    panel.grid.minor = element_blank()
  )

if (rvf_write_diagnostic_plots()) {
  diagnostic_dir <- file.path(figures_dir, "diagnostics")
  rvf_ggsave(
    file.path(diagnostic_dir, "spatial_driver_correlation_diagnostic.png"),
    p_driver_diagnostic,
    width = 11.2,
    height = 7.2
  )
}

# -------- Primary-profile spatial map ------------------------------------
# Village IDs follow the order stored in the scenario bundles. Coordinates
# are projected UTM zone 36S values from the movement-data source.
baseline_scenario <- load_rvf_scenario(
  rvf_with_replicate(
    sprintf("OBSERVED_CURRENT_%d", sim$baseline_years[[1]]),
    sort(unique(cfg$analysis_network_replicates))[[1]]
  ),
  cfg = cfg,
  validate = TRUE
)
village_coordinates <- read_csv(
  here::here(
    "data",
    "shapefiles",
    "resource_areas",
    "serengeti village centre xy.csv"
  ),
  show_col_types = FALSE
) |>
  setNames(c("village_name", "easting", "northing")) |>
  right_join(
    tibble(
      village_id = seq_along(baseline_scenario$vil_names),
      village_name = baseline_scenario$vil_names
    ),
    by = "village_name"
  )
if (nrow(village_coordinates) != 94L ||
    anyNA(village_coordinates[c("easting", "northing")])) {
  stop("Village coordinates do not match the 94 modelled villages.")
}

district_boundary <- st_read(
  here::here(
    "data",
    "shapefiles",
    "resource_areas",
    "Serengeti and Ngorongoro district shape files.shp"
  ),
  quiet = TRUE
) |>
  filter(.data$District_N == "Serengeti") |>
  summarise() |>
  st_transform(32736)

# Protected land covers about seven tenths of the district and holds almost no
# village herds, which is why the modelled villages sit in the northwest. Same
# layer as the manuscript burden maps, built by
# analysis/prepare_protected_areas.R from the World Database on
# Protected Areas.
protected_areas <- st_read(
  here::here("data", "shapefiles", "serengeti_protected_areas.gpkg"),
  quiet = TRUE
) |>
  st_transform(32736)

primary_villages <- village_summary |>
  filter(.data$site_profile == primary_profile)
current_map <- primary_villages |>
  filter(.data$scenario == "current") |>
  select(
    "village_id",
    "annual_infections_mid"
  )
map_scenarios <- tibble(
  scenario = cfg$future_scenarios,
  scenario_label = unname(cfg$scenario_labels[cfg$future_scenarios])
)
current_map <- crossing(
  current_map,
  map_scenarios
) |>
  mutate(time_label = "Current")

future_map <- primary_villages |>
  filter(.data$scenario != "current") |>
  transmute(
    .data$village_id,
    .data$annual_infections_mid,
    .data$scenario,
    .data$scenario_label,
    time_label = sprintf("%d-%d", .data$window_start, .data$window_end)
  )

map_data <- bind_rows(current_map, future_map) |>
  left_join(village_coordinates, by = "village_id") |>
  mutate(
    scenario_label = factor(
      .data$scenario_label,
      levels = unname(cfg$scenario_labels[cfg$future_scenarios])
    ),
    time_label = factor(.data$time_label, levels = time_levels)
  )

p_map <- ggplot() +
  geom_sf(
    data = district_boundary,
    fill = "white",
    colour = "grey25",
    linewidth = 0.55
  ) +
  geom_sf(
    data = protected_areas,
    fill = "grey88",
    colour = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = map_data,
    aes(
      x = .data$easting,
      y = .data$northing,
      colour = .data$annual_infections_mid
    ),
    size = 2.25
  ) +
  facet_grid(scenario_label ~ time_label) +
  scale_colour_viridis_c(
    option = "H",
    trans = "sqrt",
    breaks = c(0, 50, 100, 200, 300),
    labels = scales::label_number(accuracy = 1),
    name = "Annual cattle\ninfections"
  ) +
  coord_sf(crs = 32736, datum = NA) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 12) +
  theme(
    text = element_text(colour = "black"),
    strip.text = element_text(colour = "black", size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey94", colour = "grey65"),
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 10.5),
    legend.title = element_text(colour = "black", size = 10.5),
    legend.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  ) +
  guides(
    colour = guide_colourbar(
      title.position = "top",
      barwidth = grid::unit(8, "cm")
    )
  )

rvf_ggsave(
  file.path(figures_dir, "06_spatial_burden_map.png"),
  p_map,
  width = 11.6,
  height = 6.8
)

message("Spatial-pattern summaries written to: ", out_dir)
