#!/usr/bin/env Rscript

# Build the cross-module environmental figure and copy the small set of
# canonical analysis figures used by the active Nature Communications draft.
# Run this after the Module 01, 02, and 03 summary scripts have completed.

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tibble)
})

source(here::here("analysis", "03_rvf_transmission", "R", "config.R"))

manuscript_dir <- here::here("manuscript", "serengeti_rvf")
main_figure_dir <- file.path(manuscript_dir, "figures", "main")
supplementary_figure_dir <- file.path(
  manuscript_dir,
  "figures",
  "supplementary"
)
supplementary_table_dir <- file.path(manuscript_dir, "tables", "supplementary")

dir.create(main_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplementary_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplementary_table_dir, recursive = TRUE, showWarnings = FALSE)

annual_ndvi_file <- here::here(
  "analysis", "01_ndvi_prediction", "outputs", "comparison", "summaries",
  "annual_ndvi_ensemble.csv"
)
historical_ndvi_file <- here::here(
  "analysis", "01_ndvi_prediction", "derived", "gimms_support",
  "all_historical_observed.feather"
)
network_summary_file <- here::here(
  "analysis", "02_movement_networks", "outputs", "projections", "summaries",
  "future_network_ensemble_summary.csv"
)

required_summary_files <- c(
  annual_ndvi_file,
  historical_ndvi_file,
  network_summary_file
)
missing_summary_files <- required_summary_files[!file.exists(required_summary_files)]
if (length(missing_summary_files) > 0L) {
  stop(
    "Missing summary input(s):\n - ",
    paste(missing_summary_files, collapse = "\n - ")
  )
}

# Panel A places the observed annual series before the four-GCM projections.
# Future lines are ensemble medians and ribbons span the GCM range.
annual_ndvi <- read_csv(annual_ndvi_file, show_col_types = FALSE) |>
  mutate(
    scenario_label = factor(
      scenario_label,
      levels = c("SSP2-4.5", "SSP3-7.0")
    )
  )

historical_annual_ndvi <- read_feather(
  historical_ndvi_file,
  col_select = c("date", "ndvi")
) |>
  mutate(year = as.integer(format(as.Date(date), "%Y"))) |>
  group_by(year) |>
  summarise(mean_ndvi = mean(ndvi, na.rm = TRUE), .groups = "drop")

scenario_colours <- c("SSP2-4.5" = "#0072B2", "SSP3-7.0" = "#D55E00")

ndvi_panel <- ggplot() +
  geom_line(
    data = historical_annual_ndvi,
    aes(x = year, y = mean_ndvi, colour = "Observed 1982-2018"),
    linewidth = 0.75
  ) +
  geom_vline(
    xintercept = 2018.5,
    colour = "grey45",
    linewidth = 0.5,
    linetype = 2
  ) +
  geom_ribbon(
    data = annual_ndvi,
    aes(
      x = year,
      ymin = min_ndvi,
      ymax = max_ndvi,
      fill = scenario_label
    ),
    inherit.aes = FALSE,
    alpha = 0.22,
    colour = NA
  ) +
  geom_line(
    data = annual_ndvi,
    aes(x = year, y = median_ndvi, colour = scenario_label),
    linewidth = 0.8
  ) +
  scale_colour_manual(
    values = c("Observed 1982-2018" = "black", scenario_colours),
    breaks = c("Observed 1982-2018", "SSP2-4.5", "SSP3-7.0"),
    name = NULL
  ) +
  scale_fill_manual(values = scenario_colours, guide = "none") +
  scale_x_continuous(
    breaks = c(1982, 2000, 2020, 2040, 2060, 2080, 2100)
  ) +
  labs(x = NULL, y = "Mean annual NDVI") +
  theme_classic(base_size = 13) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    legend.position = "top"
  )

# Panel B reports the ensemble median change in monthly one-mode mean degree.
# Bars retain the range across the four GCM-specific medians. We avoid joining
# the discrete future windows because the connecting segments would imply a
# continuous trajectory that was not simulated.
network_degree <- read_csv(network_summary_file, show_col_types = FALSE) |>
  filter(metric == "mean_degree") |>
  mutate(
    scenario_label = factor(
      scenario_label,
      levels = c("SSP2-4.5", "SSP3-7.0")
    ),
    future_window = factor(
      future_window,
      levels = c("2041-2060", "2061-2080", "2081-2100")
    ),
    resource = factor(
      resource,
      levels = c("grazing", "watering"),
      labels = c("Grazing", "Watering")
    )
  )

network_panel <- ggplot(
  network_degree,
  aes(
    x = future_window,
    y = ensemble_median_change_percent,
    colour = scenario_label,
    shape = resource,
    group = interaction(scenario_label, resource)
  )
) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.55) +
  geom_linerange(
    aes(
      ymin = gcm_minimum_change_percent,
      ymax = gcm_maximum_change_percent
    ),
    position = position_dodge(width = 0.62),
    linewidth = 0.65
  ) +
  geom_point(position = position_dodge(width = 0.62), size = 2.9) +
  scale_colour_manual(values = scenario_colours, name = NULL) +
  scale_shape_manual(
    values = c("Grazing" = 16, "Watering" = 17),
    name = NULL
  ) +
  guides(colour = "none") +
  labs(
    x = NULL,
    y = "Change in mean village\ncontact degree (%)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "top"
  )

environmental_figure <- ndvi_panel + network_panel +
  plot_layout(widths = c(1.05, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    plot.tag = element_text(face = "bold", size = 15)
  )

ggsave(
  file.path(main_figure_dir, "env_pathway_combined.png"),
  environmental_figure,
  width = 12,
  height = 5.8,
  dpi = 300,
  bg = "white"
)

# Copy only figures referenced by the active main paper or Supplement. The
# analysis outputs remain canonical; these are stable publication-facing copies.
main_asset_map <- tribble(
  ~source, ~destination,
  here::here("analysis", "03_rvf_transmission", "outputs", "figures", "03_pmobile_all_profiles_relative_incidence.png"),
  file.path(main_figure_dir, "rvf_mobility_all_profiles.png"),
  here::here("analysis", "03_rvf_transmission", "outputs", "figures", "04_future_window_annual_infections_profiles_fold.png"),
  file.path(main_figure_dir, "rvf_future_profiles.png"),
  here::here("analysis", "03_rvf_transmission", "outputs", "figures", "06_spatial_broadening.png"),
  file.path(main_figure_dir, "rvf_spatial_broadening.png"),
  here::here("analysis", "03_rvf_transmission", "outputs", "figures", "08_attribution_scenario_change.png"),
  file.path(main_figure_dir, "rvf_attribution.png")
)

supplementary_asset_map <- tribble(
  ~source, ~destination,
  here::here("analysis", "01_ndvi_prediction", "outputs", "plots", "diagnostics", "selected_model_validation.png"),
  file.path(supplementary_figure_dir, "ndvi_model_validation.png"),
  here::here("analysis", "01_ndvi_prediction", "outputs", "plots", "diagnostics", "selected_model_responses.png"),
  file.path(supplementary_figure_dir, "ndvi_model_responses.png"),
  here::here("analysis", "01_ndvi_prediction", "outputs", "comparison", "figures", "climate_support_by_ssp.png"),
  file.path(supplementary_figure_dir, "ndvi_climate_support_by_ssp.png"),
  here::here("analysis", "01_ndvi_prediction", "outputs", "comparison", "figures", "seasonal_ndvi_by_ssp.png"),
  file.path(supplementary_figure_dir, "ndvi_seasonal_by_ssp.png"),
  here::here("analysis", "02_movement_networks", "outputs", "validation", "movement_regime_definition.pdf"),
  file.path(supplementary_figure_dir, "movement_regime_definition.pdf"),
  here::here("analysis", "02_movement_networks", "outputs", "validation", "movement_model_validation.pdf"),
  file.path(supplementary_figure_dir, "movement_model_validation.pdf"),
  here::here("analysis", "02_movement_networks", "outputs", "projections", "figures", "future_network_change.pdf"),
  file.path(supplementary_figure_dir, "movement_future_change.pdf"),
  here::here("analysis", "02_movement_networks", "outputs", "projections", "figures", "future_regime_change.pdf"),
  file.path(supplementary_figure_dir, "movement_future_regimes.pdf"),
  here::here("analysis", "03_rvf_transmission", "outputs", "figures", "vector_temperature_traits.pdf"),
  file.path(supplementary_figure_dir, "rvf_vector_temperature_traits.pdf"),
  here::here("analysis", "03_rvf_transmission", "outputs", "figures", "rainfall_capacity_responses.pdf"),
  file.path(supplementary_figure_dir, "rvf_rainfall_capacity_responses.pdf"),
  here::here("analysis", "03_rvf_transmission", "outputs", "figures", "model_implied_vector_seasonality.pdf"),
  file.path(supplementary_figure_dir, "rvf_vector_seasonality.pdf"),
  here::here("analysis", "03_rvf_transmission", "outputs", "figures", "07_seasonal_infectious_cycle.png"),
  file.path(supplementary_figure_dir, "rvf_seasonal_timing.png"),
  here::here("analysis", "03_rvf_transmission", "outputs", "figures", "04_rainfall_response_sensitivity.png"),
  file.path(supplementary_figure_dir, "rvf_rainfall_sensitivity.png")
)

asset_map <- bind_rows(main_asset_map, supplementary_asset_map)

missing_assets <- asset_map$source[!file.exists(asset_map$source)]
if (length(missing_assets) > 0L) {
  stop("Missing manuscript asset(s):\n - ", paste(missing_assets, collapse = "\n - "))
}

# Do not copy figures from an earlier 5- or 10-replicate run into the active
# manuscript. These compact metadata checks cover every Module 03 analysis
# represented by the publication assets.
expected_network_replicates <- length(
  unique(rvf_config$analysis_network_replicates)
)
expected_future_runs_per_summary <-
  length(rvf_config$climate_models) * expected_network_replicates

calibration_selected_file <- file.path(
  rvf_config$results_dir,
  "calibrate_seroprevalence",
  "calibrate_recommended.csv"
)
mobility_summary_file <- file.path(
  rvf_config$results_dir,
  "mobility_capacity_profile_comparison",
  "mobility_profile_summary.csv"
)
future_summary_files <- file.path(
  rvf_config$results_dir,
  c(
    "future_cold_endemic",
    "future_cold_endemic_profile_peridomestic_dominant",
    "future_cold_endemic_quadratic_rainfall"
  ),
  "future_cold_endemic_window_burden_aggregate.csv"
)
attribution_summary_files <- file.path(
  rvf_config$results_dir,
  c(
    "climate_movement_attribution_cold_start",
    "climate_movement_attribution_cold_start_profile_peridomestic_dominant"
  ),
  "climate_movement_attribution_decomposition_aggregate.csv"
)
spatial_summary_file <- file.path(
  rvf_config$results_dir,
  "spatial_patterns",
  "spatial_concentration_summary.csv"
)
seasonal_summary_file <- file.path(
  rvf_config$results_dir,
  "seasonal_cycles",
  "seasonal_characteristics_summary.csv"
)
module03_summary_files <- c(
  calibration_selected_file,
  mobility_summary_file,
  future_summary_files,
  attribution_summary_files,
  spatial_summary_file,
  seasonal_summary_file
)
missing_module03_summaries <- module03_summary_files[
  !file.exists(module03_summary_files)
]
if (length(missing_module03_summaries) > 0L) {
  stop(
    "Missing current Module 03 summary file(s):\n - ",
    paste(missing_module03_summaries, collapse = "\n - ")
  )
}

calibration_check <- read_csv(
  calibration_selected_file,
  show_col_types = FALSE
)
mobility_check <- read_csv(mobility_summary_file, show_col_types = FALSE)
future_checks <- lapply(
  future_summary_files,
  read_csv,
  show_col_types = FALSE
)
attribution_checks <- lapply(
  attribution_summary_files,
  read_csv,
  show_col_types = FALSE
)

protocol_mismatch <- function(data, window_field) {
  if (!"n_network_replicates" %in% names(data) ||
      any(data$n_network_replicates != expected_network_replicates)) {
    return(TRUE)
  }
  !window_field %in% names(data) ||
    any(data[[window_field]] != rvf_config$execution$reporting_years)
}

if (protocol_mismatch(calibration_check, "calibration_window_years") ||
    protocol_mismatch(mobility_check, "terminal_summary_years") ||
    any(vapply(
      future_checks,
      protocol_mismatch,
      window_field = "window_length_years",
      logical(1)
    )) ||
    any(vapply(
      attribution_checks,
      protocol_mismatch,
      window_field = "terminal_summary_years",
      logical(1)
    ))) {
  stop(
    "Module 03 manuscript assets do not match the configured ",
    expected_network_replicates, "-replicate, ",
    rvf_config$execution$reporting_years, "-year reporting protocol."
  )
}

spatial_check <- read_csv(
  spatial_summary_file,
  show_col_types = FALSE
)
seasonal_check <- read_csv(
  seasonal_summary_file,
  show_col_types = FALSE
)
if (any(spatial_check$n_runs[spatial_check$scenario == "current"] !=
        expected_network_replicates) ||
    any(spatial_check$n_runs[spatial_check$scenario != "current"] !=
        expected_future_runs_per_summary) ||
    any(seasonal_check$n_runs[seasonal_check$scenario == "current"] !=
        expected_network_replicates) ||
    any(seasonal_check$n_runs[seasonal_check$scenario != "current"] !=
        expected_future_runs_per_summary)) {
  stop("Spatial or seasonal manuscript assets use a stale replicate grid.")
}

copied <- file.copy(asset_map$source, asset_map$destination, overwrite = TRUE)
if (!all(copied)) stop("At least one manuscript asset could not be copied.")

# Regenerate the climate-input table, then verify every external SI table. The
# remaining short tables are maintained as auditable LaTeX fragments beside it.
source(
  here::here("analysis", "build_supp_climate_input_table.R"),
  local = TRUE
)

# Generate the vegetation-model comparison directly from the temporal-holdout
# output. Formula labels below are a readable rendering of the exact candidate
# definitions in 03_fit_sdmtmb_candidates.R; elevation and month are common to
# every candidate and are stated once in the table caption.
ndvi_comparison_file <- here::here(
  "analysis", "01_ndvi_prediction", "outputs", "selection",
  "sdmtmb_candidate_comparison_observed.csv"
)
if (!file.exists(ndvi_comparison_file)) {
  stop("Missing vegetation candidate comparison: ", ndvi_comparison_file)
}

ndvi_formula_labels <- c(
  sdmtmb_winner_t2 = "$t_2(P_3,T_{\\max,3})+s(P_{t-1})+s(Y)$",
  sdmtmb_smooth_core_tmax_syear_int_plus_prlag1 = "$s(P_3)+s(P_{t-1})+s(T_{\\max,3})+s(P_3T_{\\max,3})+s(Y)$",
  sdmtmb_smooth_core_tmax_year_int_plus_prlag1 = "$s(P_3)+s(P_{t-1})+s(T_{\\max,3})+s(P_3T_{\\max,3})+Y$",
  sdmtmb_smooth_core_tmax_year_int_plus_tmax_now = "$s(P_3)+s(T_{\\max,3})+s(T_{\\max,t})+s(P_3T_{\\max,3})+Y$",
  sdmtmb_smooth_core_tmax_year = "$s(P_3)+s(T_{\\max,3})+Y$",
  sdmtmb_smooth_core_tmax_year_int = "$s(P_3)+s(T_{\\max,3})+s(P_3T_{\\max,3})+Y$",
  sdmtmb_smooth_memory_tmax_year_int = "$s(P_6)+s(\\Delta P)+s(T_{\\max,3})+s(\\Delta P\\,T_{\\max,3})+Y$",
  sdmtmb_smooth_core_tmean_year = "$s(P_3)+s(T_{\\mathrm{mean},3})+Y$",
  sdmtmb_smooth_core_tmax = "$s(P_3)+s(T_{\\max,3})$",
  sdmtmb_smooth_core_tmin_year = "$s(P_3)+s(T_{\\min,3})+Y$",
  sdmtmb_linear_core_tmax = "$P_3+T_{\\max,3}$",
  sdmtmb_baseline_season_year = "$Y$"
)

ndvi_candidates <- read_csv(ndvi_comparison_file, show_col_types = FALSE) |>
  mutate(specification = unname(ndvi_formula_labels[.data$model_id]))
if (!"selected" %in% names(ndvi_candidates)) {
  ndvi_candidates <- ndvi_candidates |>
    mutate(selected = row_number() == 1L)
}
if (any(is.na(ndvi_candidates$specification))) {
  stop(
    "Missing manuscript label for vegetation candidate(s): ",
    paste(ndvi_candidates$model_id[is.na(ndvi_candidates$specification)],
          collapse = ", ")
  )
}

ndvi_table_rows <- character(nrow(ndvi_candidates))
for (row_index in seq_len(nrow(ndvi_candidates))) {
  candidate <- ndvi_candidates[row_index, ]
  values <- c(
    candidate$specification,
    sprintf("%.4f", candidate$holdout_rmse),
    sprintf("%.3f", candidate$holdout_r2),
    sprintf("%.3f", candidate$holdout_monthly_mean_r),
    sprintf("%.4f", candidate$holdout_monthly_mean_rmse)
  )
  if (isTRUE(candidate$selected[[1]])) {
    values <- paste0("\\textbf{", values, "}")
  }
  ndvi_table_rows[[row_index]] <- paste0(
    "    ", paste(values, collapse = " & "), " \\\\"
  )
}

ndvi_table_lines <- c(
  "\\begin{table}[htbp]",
  "  \\centering",
  "  \\caption{\\textbf{Candidate vegetation-model comparison.} Twelve",
  "  \\texttt{sdmTMB} specifications were trained on monthly observations before",
  "  2015 and evaluated on the 2015--2018 temporal holdout. The two leading",
  "  specifications differed in RMSE by less than $10^{-4}$ and were treated as",
  "  tied. The selected tensor-product model is shown first; it permits rainfall",
  "  and temperature to have separate smoothness scales. All candidates used a",
  "  Beta response, a Mat\\'{e}rn spatial random",
  "  field, month fixed effects, and elevation. $P_3$ and $P_6$ are",
  "  three- and six-month rainfall means, $\\Delta P=P_3-P_6$, $P_{t-1}$ is",
  "  one-month lagged rainfall, temperature subscripts 3 and $t$ denote a",
  "  three-month mean and the current month, and $Y$ is year.}",
  "  \\label{tab:supp_ndvi_candidate_comparison}",
  "  \\small",
  "  \\setlength{\\tabcolsep}{2.5pt}",
  "  \\begin{tabularx}{0.98\\textwidth}{@{}>{\\raggedright\\arraybackslash}Xcccc@{}}",
  "    \\toprule",
  "    Specification & RMSE & $r^2$ & Monthly $r$ & Monthly RMSE \\\\",
  "    \\midrule",
  ndvi_table_rows,
  "    \\bottomrule",
  "  \\end{tabularx}",
  "\\end{table}"
)
writeLines(
  ndvi_table_lines,
  file.path(supplementary_table_dir, "ndvi_candidate_comparison.tex")
)

# Regenerate the movement-model fit table from the final observed-cycle checks.
# The checks use new random seeds but the same 2017 data used for calibration;
# they therefore measure model adequacy rather than out-of-sample prediction.
movement_validation_file <- here::here(
  "analysis", "02_movement_networks", "outputs", "validation",
  "regime_network_validation.csv"
)
if (!file.exists(movement_validation_file)) {
  stop("Missing movement-model validation summary: ", movement_validation_file)
}
movement_validation <- read_csv(
  movement_validation_file,
  show_col_types = FALSE
)

movement_resources <- c(grazing = "Grazing", watering = "Watering")
movement_regimes <- c(routine = "Routine", transition = "Transition", drought = "Drought")

movement_fit_rows <- character()
for (resource_id in names(movement_resources)) {
  for (regime_id in names(movement_regimes)) {
    edge_row <- movement_validation |>
      filter(
        .data$resource == resource_id,
        .data$regime == regime_id,
        .data$metric == "edges"
      )
    distance_row <- movement_validation |>
      filter(
        .data$resource == resource_id,
        .data$regime == regime_id,
        .data$metric == "p95_dist"
      )
    if (nrow(edge_row) != 1L || nrow(distance_row) != 1L) {
      stop(
        "Expected one edge and one p95-distance row for ",
        resource_id, " / ", regime_id, "."
      )
    }

    resource_label <- if (regime_id == "routine") {
      movement_resources[[resource_id]]
    } else {
      ""
    }
    movement_fit_rows <- c(
      movement_fit_rows,
      sprintf(
        "    %s & %s & %.0f & %.0f [%.0f, %.0f] & %.1f & %.1f [%.1f, %.1f] \\\\",
        resource_label,
        movement_regimes[[regime_id]],
        edge_row$observed,
        edge_row$median,
        edge_row$q10,
        edge_row$q90,
        distance_row$observed,
        distance_row$median,
        distance_row$q10,
        distance_row$q90
      )
    )
  }
}

writeLines(
  c(
    "\\begin{table}[htbp]",
    "  \\centering",
    "  \\caption{\\textbf{Observed and simulated movement-network summaries.} Simulated values are medians across $100$ new stochastic realisations from the fitted point-estimate model, with $10$th--$90$th percentile intervals in brackets. Observed values are regime-mean monthly edge counts and 95th-percentile edge distances from the 2017 participatory-mapping network. These checks use the same data as calibration and therefore assess model adequacy rather than out-of-sample prediction.}",
    "  \\label{tab:movement_ppc}",
    "  \\small",
    "  \\setlength{\\tabcolsep}{4pt}",
    "  \\begin{tabular}{@{}llrrrr@{}}",
    "    \\toprule",
    "    & & \\multicolumn{2}{c}{Monthly edges} & \\multicolumn{2}{c}{\\begin{tabular}[c]{@{}c@{}}95th-percentile edge\\\\distance (km)\\end{tabular}} \\\\",
    "    \\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
    "    Resource & Regime & Observed & Simulated & Observed & Simulated \\\\",
    "    \\midrule",
    movement_fit_rows,
    "    \\bottomrule",
    "  \\end{tabular}",
    "\\end{table}"
  ),
  file.path(supplementary_table_dir, "movement_fit_quantile_summary.tex")
)

# Regenerate fitted movement parameters from the canonical calibration output.
# Transition gamma is derived as gamma_yr * m_tr and is not listed separately.
movement_estimates_file <- here::here(
  "analysis", "02_movement_networks", "outputs", "calibration",
  "point_estimates.csv"
)
if (!file.exists(movement_estimates_file)) {
  stop("Missing movement point estimates: ", movement_estimates_file)
}
movement_estimates <- read_csv(
  movement_estimates_file,
  show_col_types = FALSE
)

movement_parameter_rows <- tribble(
  ~layer_id, ~layer_label, ~parameter_id, ~parameter_label,
  "year_round", "Year-round", "lambda", "$\\lambda_{\\mathrm{yr}}$",
  "year_round", "", "gamma", "$\\gamma_{\\mathrm{yr}}$",
  "year_round", "", "beta_z", "$\\beta_{z,\\mathrm{yr}}$",
  "transition", "Transition", "lambda", "$\\lambda_{\\mathrm{tr}}$",
  "transition", "", "m_sh", "$m_{\\mathrm{tr}}$",
  "transition", "", "beta_z", "$\\beta_{z,\\mathrm{tr}}$",
  "drought", "Drought", "lambda_local", "$\\lambda_{\\mathrm{local}}$",
  "drought", "", "gamma", "$\\gamma_{\\mathrm{dr}}$",
  "drought", "", "beta_z", "$\\beta_{z,\\mathrm{dr}}$",
  "drought", "", "p_keep", "$p_{\\mathrm{keep}}$",
  "drought", "", "lambda_tail", "$\\lambda_{\\mathrm{tail}}$",
  "drought", "", "tau", "$\\tau$"
)

movement_estimate_rows <- vapply(
  seq_len(nrow(movement_parameter_rows)),
  function(i) {
    parameter_row <- movement_parameter_rows[i, ]
    estimates <- movement_estimates |>
      filter(
        .data$layer == parameter_row$layer_id,
        .data$parameter == parameter_row$parameter_id
      )
    grazing_estimate <- filter(estimates, .data$resource == "grazing")
    watering_estimate <- filter(estimates, .data$resource == "watering")
    if (nrow(grazing_estimate) != 1L || nrow(watering_estimate) != 1L) {
      stop(
        "Expected one grazing and one watering estimate for ",
        parameter_row$layer_id, " / ", parameter_row$parameter_id, "."
      )
    }
    sprintf(
      "    %s & %s & %.3f & %.3f \\\\",
      parameter_row$layer_label,
      parameter_row$parameter_label,
      grazing_estimate$estimate,
      watering_estimate$estimate
    )
  },
  character(1)
)

# Add visual spacing between the three fitted layers.
movement_estimate_rows <- append(movement_estimate_rows, "    \\addlinespace", after = 3L)
movement_estimate_rows <- append(movement_estimate_rows, "    \\addlinespace", after = 7L)

writeLines(
  c(
    "\\begin{table}[htbp]",
    "  \\centering",
    "  \\caption{\\textbf{Fitted movement-model point estimates.} Estimates were obtained by sequential multi-start Nelder--Mead optimisation against the 2017 observed village-network summaries. The transition distance-decay coefficient is derived as $\\gamma_{\\mathrm{tr}}=\\gamma_{\\mathrm{yr}}m_{\\mathrm{tr}}$.}",
    "  \\label{tab:supp_movement_point_estimates}",
    "  \\small",
    "  \\begin{tabular}{llrr}",
    "    \\toprule",
    "    Layer & Parameter & Grazing & Watering \\\\",
    "    \\midrule",
    movement_estimate_rows,
    "    \\bottomrule",
    "  \\end{tabular}",
    "\\end{table}"
  ),
  file.path(supplementary_table_dir, "movement_point_estimates.tex")
)

# Build a compact table of the selected calibration results. The full response
# curves remain in the analysis outputs for diagnostics; the paper needs the
# selected capacity scales and the seroprevalence they produced.
calibration_selected <- read_csv(
  calibration_selected_file,
  show_col_types = FALSE
)
profile_display <- c(
  grazing_floodplain_dominant = "Shared grazing/floodplain",
  peridomestic_dominant = "Village-associated",
  balanced = "Balanced",
  resource_edge_dominant = "Resource-edge",
  waterpoint_extreme = "Waterpoint-extreme"
)
profile_order <- c(
  "Shared grazing/floodplain",
  "Village-associated",
  "Balanced",
  "Resource-edge",
  "Waterpoint-extreme"
)
calibration_table <- calibration_selected |>
  transmute(
    profile = unname(profile_display[.data$site_profile]),
    rainfall_model = .data$rainfall_model,
    capacity_scale = .data$capacity_scale,
    seroprevalence = sprintf(
      "%.2f [%.2f, %.2f]",
      .data$mean_seroprevalence_calibration,
      .data$mean_seroprevalence_calibration_p10,
      .data$mean_seroprevalence_calibration_p90
    )
  ) |>
  mutate(profile = factor(.data$profile, levels = profile_order)) |>
  arrange(.data$profile)

calibration_rows <- vapply(profile_order, function(profile_name) {
  profile_data <- filter(calibration_table, .data$profile == profile_name)
  hill <- filter(profile_data, .data$rainfall_model == "saturating_rainfall")
  quadratic <- filter(profile_data, .data$rainfall_model == "quadratic_rainfall")
  sprintf(
    "%s & %.2f & %s & %.2f & %s \\\\",
    profile_name,
    hill$capacity_scale,
    hill$seroprevalence,
    quadratic$capacity_scale,
    quadratic$seroprevalence
  )
}, character(1))

writeLines(
  c(
    "\\begin{table}[htbp]",
    "\\centering",
    sprintf(
      "\\caption{\\textbf{Selected calibration results for the five exposure profiles.} Mean cattle seroprevalence is followed by the 10th--90th percentile across %d movement-network replicates. All selected mean values fall within the 4--10\\%% reference range.}",
      expected_network_replicates
    ),
    "\\label{tab:supp_calibration_profiles}",
    "\\small",
    "\\begin{tabular}{lcccc}",
    "\\toprule",
    "& \\multicolumn{2}{c}{Hill} & \\multicolumn{2}{c}{Quadratic} \\\\",
    "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
    "Profile & $c_K$ & Seroprevalence [P10--P90] (\\%) & $c_K$ & Seroprevalence [P10--P90] (\\%) \\\\",
    "\\midrule",
    calibration_rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ),
  file.path(supplementary_table_dir, "calibration_profile_results.tex")
)

# Confirm every table this script writes was written. The two hand-maintained
# Supplement tables, ndvi_mesh_sensitivity.tex and
# movement_calibration_weights.tex, are not checked here because they live only
# in the manuscript source; a missing one stops the LaTeX build instead. This
# lets the script run in a copy of the code that has no manuscript folder.
required_tables <- file.path(
  supplementary_table_dir,
  c(
    "supp_climate_input_summary.tex",
    "ndvi_candidate_comparison.tex",
    "movement_fit_quantile_summary.tex",
    "movement_point_estimates.tex",
    "calibration_profile_results.tex"
  )
)
missing_tables <- required_tables[!file.exists(required_tables)]
if (length(missing_tables) > 0L) {
  stop("Missing manuscript table(s):\n - ", paste(missing_tables, collapse = "\n - "))
}

message(
  "Built the environmental pathway figure and refreshed ",
  nrow(asset_map),
  " figure assets and ",
  length(required_tables),
  " verified supplementary table assets under: ",
  manuscript_dir
)
