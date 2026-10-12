#!/usr/bin/env Rscript
#' Summarise the seasonal adult-vector abundance implied by Module 03.
#'
#' This diagnostic isolates the total adult-vector equation from the infection
#' compartments. That reduction is exact because infection only redistributes
#' adult mosquitoes among S, E, and I; it does not change total recruitment or
#' mortality. The resulting vector-to-host ratios are model-implied quantities,
#' not field estimates of mosquito density.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

source(here::here("analysis", "03_rvf_transmission", "R", "model_parameters.R"))

cfg <- load_rvf_config()
out_dir <- file.path(cfg$results_dir, "diagnostics", "vector_abundance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The diagnostic uses the two profiles carried through the future analysis.
profiles <- c(
  "grazing_floodplain_dominant",
  "peridomestic_dominant"
)
rainfall_model <- cfg$primary_rainfall_model

calibration <- read_csv(
  rvf_calibration_recommended_path(cfg),
  show_col_types = FALSE
) |>
  filter(
    .data$rainfall_model == .env$rainfall_model,
    .data$site_profile %in% .env$profiles
  ) |>
  select(
    "site_profile",
    "site_profile_label",
    "rainfall_model",
    "capacity_scale"
  )

if (nrow(calibration) != length(profiles)) {
  stop("The primary-rainfall calibration table does not contain both profiles.")
}

# Function: advance total adult abundance through one interval with constant
# carrying capacity and mortality. Below K this is the exact logistic solution;
# above K, recruitment is zero until natural mortality brings abundance below K.
advance_adults <- function(M, K, mu, recruitment_rate, interval_days) {
  M <- pmax(as.numeric(M), 0)
  K <- pmax(as.numeric(K), 1)
  mu <- pmax(as.numeric(mu), 0)
  updated <- numeric(length(M))

  above_capacity <- M > K
  if (any(above_capacity)) {
    above_index <- which(above_capacity)
    time_to_capacity <- log(M[above_index] / K[above_index]) /
      pmax(mu[above_index], .Machine$double.eps)
    remains_above <- time_to_capacity >= interval_days

    if (any(remains_above)) {
      index <- above_index[remains_above]
      updated[index] <- M[index] * exp(-mu[index] * interval_days)
    }

    if (any(!remains_above)) {
      index <- above_index[!remains_above]
      remaining_days <- interval_days - time_to_capacity[!remains_above]
      net_growth <- recruitment_rate - mu[index]
      density_coefficient <- recruitment_rate / K[index]
      starting_abundance <- K[index]

      updated[index] <- ifelse(
        abs(net_growth) < 1e-12,
        starting_abundance /
          (1 + density_coefficient * starting_abundance * remaining_days),
        net_growth * starting_abundance /
          (
            density_coefficient * starting_abundance +
              (net_growth - density_coefficient * starting_abundance) *
                exp(-net_growth * remaining_days)
          )
      )
    }
  }

  if (any(!above_capacity)) {
    index <- which(!above_capacity)
    net_growth <- recruitment_rate - mu[index]
    density_coefficient <- recruitment_rate / K[index]
    starting_abundance <- M[index]

    updated[index] <- ifelse(
      abs(net_growth) < 1e-12,
      starting_abundance /
        (1 + density_coefficient * starting_abundance * interval_days),
      net_growth * starting_abundance /
        (
          density_coefficient * starting_abundance +
            (net_growth - density_coefficient * starting_abundance) *
              exp(-net_growth * interval_days)
        )
    )
  }

  pmax(updated, 0)
}

# Function: return cattle present at each site. Resource-site totals use the
# same row-normalised commute matrices and p_mobile = 1 reference as the model.
hosts_at_sites <- function(params, month_index, site_type) {
  if (site_type == "Village") {
    return(rep(params$N_h_const, params$n_vil))
  }

  mobile_hosts <- rep(params$N_h_const, params$n_vil)
  if (site_type == "Grazing") {
    return(as.vector(t(params$G_list[[month_index]]) %*% mobile_hosts))
  }
  as.vector(t(params$W_list[[month_index]]) %*% mobile_hosts)
}

# Function: calculate the recurrent adult-vector cycle for one forcing sequence.
# Daily updates give a time-weighted monthly mean, and four complete forcing
# cycles are more than sufficient for the short-lived adult population to settle.
summarise_vector_cycle <- function(params, cell, profile_label) {
  site_inputs <- list(
    Village = list(K = params$K_v, mu = params$mu_m_v),
    Grazing = list(K = params$K_g, mu = params$mu_m_g),
    Watering = list(K = params$K_w, mu = params$mu_m_w)
  )
  n_months <- length(params$time_labels)
  cycle_days <- as.integer(round(n_months * params$days_per_step))
  n_cycles <- 4L

  monthly_rows <- list()
  convergence_rows <- list()

  for (site_type in names(site_inputs)) {
    site <- site_inputs[[site_type]]
    # Host loads change only at monthly forcing boundaries. Precomputing them
    # avoids repeating large village-to-resource matrix products every day.
    monthly_host_loads <- lapply(
      seq_len(n_months),
      hosts_at_sites,
      params = params,
      site_type = site_type
    )
    first_mu <- site$mu[[1]]
    equilibrium_fraction <- pmax(0, 1 - first_mu / params$r_m)
    adults <- site$K[[1]] * equilibrium_fraction
    cycle_end_states <- vector("list", n_cycles)

    # Only the final forcing cycle is retained for seasonal summaries.
    final_daily_ratio <- numeric(cycle_days)
    final_daily_capacity_ratio <- numeric(cycle_days)
    final_daily_abundance <- numeric(cycle_days)
    final_month <- integer(cycle_days)

    for (cycle_index in seq_len(n_cycles)) {
      for (day_index in seq_len(cycle_days)) {
        month_index <- floor((day_index - 1) / params$days_per_step) + 1L
        month_index <- min(month_index, n_months)

        adults <- advance_adults(
          M = adults,
          K = site$K[[month_index]],
          mu = site$mu[[month_index]],
          recruitment_rate = params$r_m,
          interval_days = 1
        )

        if (cycle_index == n_cycles) {
          hosts <- monthly_host_loads[[month_index]]
          occupied <- is.finite(hosts) & hosts >= 1

          # This is the cattle-weighted mean of site-specific M/N ratios.
          # Sites with no visiting cattle do not contribute to exposure.
          final_daily_ratio[day_index] <-
            sum(adults[occupied]) / sum(hosts[occupied])
          final_daily_capacity_ratio[day_index] <-
            sum(site$K[[month_index]][occupied]) / sum(hosts[occupied])
          final_daily_abundance[day_index] <- mean(adults[occupied])
          final_month[day_index] <- ((month_index - 1L) %% 12L) + 1L
        }
      }
      cycle_end_states[[cycle_index]] <- adults
    }

    previous <- cycle_end_states[[n_cycles - 1L]]
    current <- cycle_end_states[[n_cycles]]
    convergence_rows[[site_type]] <- tibble(
      scenario = cell$scenario,
      scenario_label = cell$scenario_label,
      climate_model = cell$climate_model,
      future_window = cell$future_window,
      site_profile = params$site_profile,
      site_profile_label = profile_label,
      site_type = site_type,
      max_relative_cycle_difference = max(
        abs(current - previous) / pmax(abs(current), 1),
        na.rm = TRUE
      )
    )

    monthly_rows[[site_type]] <- tibble(
      month_num = final_month,
      vector_host_ratio = final_daily_ratio,
      capacity_host_ratio = final_daily_capacity_ratio,
      adult_vectors_per_occupied_site = final_daily_abundance
    ) |>
      group_by(.data$month_num) |>
      summarise(
        vector_host_ratio = mean(.data$vector_host_ratio, na.rm = TRUE),
        capacity_host_ratio = mean(.data$capacity_host_ratio, na.rm = TRUE),
        adult_vectors_per_occupied_site = mean(
          .data$adult_vectors_per_occupied_site,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) |>
      mutate(
        scenario = cell$scenario,
        scenario_label = cell$scenario_label,
        climate_model = cell$climate_model,
        future_window = cell$future_window,
        site_profile = params$site_profile,
        site_profile_label = profile_label,
        site_type = site_type
      )
  }

  list(
    monthly = bind_rows(monthly_rows),
    convergence = bind_rows(convergence_rows)
  )
}

# Observed forcing is one 10-year cycle. Future forcing comprises each complete
# SSP x GCM x 20-year window. Network replicate 1 is sufficient here because
# the aim is to expose model-implied vector seasonality, not quantify network
# uncertainty already propagated in the main analyses.
diagnostic_cells <- tibble(
  scenario = "observed",
  scenario_label = "Observed",
  climate_model = "Observed",
  future_window = "2009-2018",
  years = list(cfg$current_years)
)

for (scenario_id in cfg$future_scenarios) {
  for (climate_model in cfg$climate_models) {
    for (future_window in names(cfg$future_windows)) {
      diagnostic_cells <- bind_rows(
        diagnostic_cells,
        tibble(
          scenario = scenario_id,
          scenario_label = unname(cfg$scenario_labels[[scenario_id]]),
          climate_model = climate_model,
          future_window = future_window,
          years = list(cfg$future_windows[[future_window]])
        )
      )
    }
  }
}

monthly_results <- list()
convergence_results <- list()
result_index <- 0L

for (cell_index in seq_len(nrow(diagnostic_cells))) {
  cell <- diagnostic_cells[cell_index, ]
  years <- unlist(cell$years[[1]], use.names = FALSE)
  scenario_base_ids <- if (cell$scenario == "observed") {
    sprintf("OBSERVED_CURRENT_%d", years)
  } else {
    sprintf(
      "%s_%s_%d",
      toupper(cell$scenario),
      cell$climate_model,
      years
    )
  }
  scenario_ids <- rvf_with_replicate(scenario_base_ids, 1L)
  scenario <- build_rvf_sequence_scenario(
    scenario_ids,
    cfg = cfg,
    validate = TRUE,
    scenario_type = "vector_abundance_diagnostic",
    sequence_label = paste(
      cell$scenario_label,
      cell$climate_model,
      cell$future_window,
      sep = "_"
    )
  )

  for (profile_index in seq_len(nrow(calibration))) {
    profile <- calibration[profile_index, ]
    calibrated <- rvf_adjust_calibrated_params(
      capacity_scale = profile$capacity_scale,
      site_profile = profile$site_profile,
      rainfall_model = profile$rainfall_model,
      transmission_scale = 1
    )
    params <- rvf_assemble_params(
      scenario = scenario,
      calibrated = calibrated,
      use_temp_dependence = TRUE,
      p_mobile = 1,
      capacity_reference_p_mobile = 1,
      verbose = FALSE
    )

    result_index <- result_index + 1L
    result <- summarise_vector_cycle(
      params = params,
      cell = cell,
      profile_label = profile$site_profile_label
    )
    monthly_results[[result_index]] <- result$monthly
    convergence_results[[result_index]] <- result$convergence
  }
}

monthly <- bind_rows(monthly_results) |>
  mutate(
    month = factor(month.abb[.data$month_num], levels = month.abb),
    site_type = factor(
      .data$site_type,
      levels = c("Village", "Grazing", "Watering")
    ),
    site_profile_label = factor(
      .data$site_profile_label,
      levels = calibration$site_profile_label
    )
  )
convergence <- bind_rows(convergence_results)

write_csv(monthly, file.path(out_dir, "model_implied_vector_seasonality.csv"))
write_csv(convergence, file.path(out_dir, "vector_cycle_convergence.csv"))

# The publication figure compares observed and late-century seasonal cycles.
# Future bands span the four GCM-specific 20-year means; the observed line is
# the mean seasonal cycle over 2009-2018.
observed_plot <- monthly |>
  filter(.data$scenario == "observed") |>
  transmute(
    .data$site_profile_label,
    .data$site_type,
    .data$month,
    forcing = "Observed 2009-2018",
    mid = .data$vector_host_ratio,
    low = .data$vector_host_ratio,
    high = .data$vector_host_ratio
  )

future_plot <- monthly |>
  filter(.data$scenario != "observed", .data$future_window == "2100") |>
  group_by(
    .data$site_profile_label,
    .data$site_type,
    .data$month,
    .data$scenario_label
  ) |>
  summarise(
    mid = median(.data$vector_host_ratio, na.rm = TRUE),
    low = min(.data$vector_host_ratio, na.rm = TRUE),
    high = max(.data$vector_host_ratio, na.rm = TRUE),
    .groups = "drop"
  ) |>
  rename(forcing = "scenario_label")

plot_data <- bind_rows(observed_plot, future_plot) |>
  mutate(
    profile_display = case_when(
      grepl("Grazing/floodplain", as.character(.data$site_profile_label)) ~
        "Shared grazing/floodplain",
      grepl("Peridomestic|Village-associated", as.character(.data$site_profile_label)) ~
        "Village-associated",
      TRUE ~ as.character(.data$site_profile_label)
    ),
    profile_display = factor(
      .data$profile_display,
      levels = c(
        "Shared grazing/floodplain",
        "Village-associated"
      )
    ),
    forcing = factor(
      .data$forcing,
      levels = c("Observed 2009-2018", "SSP2-4.5", "SSP3-7.0")
    )
  )

forcing_colours <- c(
  "Observed 2009-2018" = "black",
  "SSP2-4.5" = "#0072B2",
  "SSP3-7.0" = "#D55E00"
)

vector_plot <- ggplot(
  plot_data,
  aes(
    x = .data$month,
    y = .data$mid,
    colour = .data$forcing,
    fill = .data$forcing,
    group = .data$forcing
  )
) +
  geom_ribbon(
    data = filter(plot_data, .data$forcing != "Observed 2009-2018"),
    aes(ymin = .data$low, ymax = .data$high),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(linewidth = 0.75) +
  facet_grid(
    rows = vars(.data$profile_display),
    cols = vars(.data$site_type),
    scales = "free_y"
  ) +
  scale_colour_manual(values = forcing_colours, drop = FALSE) +
  scale_fill_manual(values = forcing_colours, drop = FALSE, guide = "none") +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = NULL,
    y = "Cattle-weighted adult vector-to-host ratio",
    colour = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(colour = "black", size = 9.2),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(colour = "black", size = 10.8),
    legend.position = "top",
    legend.text = element_text(colour = "black", size = 9.5),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", face = "bold", size = 9.5),
    panel.spacing = grid::unit(0.9, "lines"),
    plot.margin = margin(5.5, 16, 5.5, 5.5)
  )

ggsave(
  file.path(cfg$figures_dir, "model_implied_vector_seasonality.pdf"),
  vector_plot,
  width = 7.6,
  height = 5.2,
  device = grDevices::cairo_pdf
)

max_cycle_difference <- max(
  convergence$max_relative_cycle_difference,
  na.rm = TRUE
)
if (!is.finite(max_cycle_difference) || max_cycle_difference > 1e-6) {
  warning(
    "Adult-vector cycles did not meet the 1e-6 convergence tolerance; max = ",
    signif(max_cycle_difference, 3)
  )
}

message("Vector-abundance diagnostics written to: ", out_dir)
message("Maximum final-cycle relative difference: ", signif(max_cycle_difference, 3))
