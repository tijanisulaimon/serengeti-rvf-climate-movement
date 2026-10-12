#!/usr/bin/env Rscript

#' 03_run_observed_baseline.R
#' -------------------------------------------------------------------------
#' Establish the endemic baseline used by the downstream scenario summaries.
#' Observed baseline: cycles the configured observed CHELSA years (default
#' 2009--2018) through burn-in and endemic phases. This keeps real
#' interannual variability instead of calibrating to a single current year.
#' Saves compact summaries, a final-state bundle, and the baseline convergence
#' plot. Full per-replicate ODE output and extra diagnostics are opt-in because
#' they can be several GB or duplicate downstream figures.
#'
#' Parameters are read from `cfg$sim` (see `R/config.R`). To do a non-default
#' run, edit the `sim$...` assignments in the CONFIG block below.
#' -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
})

source(here::here("analysis", "03_rvf_transmission", "R", "transmission_model.R"))
source(here::here("analysis", "03_rvf_transmission", "R", "cycle_diagnostics.R"))

cfg <- load_rvf_config()

# -------- Analysis specification -----------------------------------------
sim <- cfg$sim
if (exists("analysis_profile", inherits = FALSE)) {
  sim$site_profile <- analysis_profile
}
if (exists("analysis_rainfall_model", inherits = FALSE)) {
  sim$rainfall_model <- analysis_rainfall_model
}
if (isTRUE(cfg$execution$use_calibrated_profiles)) {
  sim <- rvf_apply_profile_calibration(
    sim,
    cfg = cfg,
    require_target_band = isTRUE(cfg$execution$require_calibration_target)
  )
}

# Calibration, baseline, mobility, and future comparisons use the same
# establishment assay: one infectious cow is introduced in every village.
if (!identical(sim$seed_village, "all") || sim$n_seed != 1L) {
  stop("The publication baseline requires one infectious cow in every village.")
}

# Keep compact tables separate from nested simulation and cycle objects.
out_dir <- rvf_baseline_dir(sim, cfg)
intermediate_dir <- rvf_baseline_intermediate_dir(sim, cfg)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

# Observed multi-year baseline. `baseline_year` is retained only as the
# legacy single-year reference; the current default uses `baseline_years`.
current_years <- sim$baseline_years %||% sim$baseline_year
scenario_base_ids <- sprintf("OBSERVED_CURRENT_%d", current_years)
# Pull observed scenario bundles for the chosen baseline years and retained
# analysis replicates, then run one coherent climate/network sequence per
# replicate.
requested_network_replicates <- if (exists("analysis_network_replicates", inherits = FALSE)) {
  analysis_network_replicates
} else {
  cfg$analysis_network_replicates
}
scenario_index <- rvf_select_scenarios(
  dataset = "observed",
  year = current_years,
  network_replicate = requested_network_replicates
)

if (nrow(scenario_index) == 0) {
  stop("No observed scenario bundles found for year(s): ",
       paste(current_years, collapse = ", "),
       "\nRun 01_prepare_inputs.R first.")
}

network_replicates <- rvf_complete_replicates_for_years(
  index = scenario_index,
  dataset = "observed",
  years = current_years,
  network_replicate = requested_network_replicates
)
if (length(network_replicates) == 0 || any(is.na(network_replicates))) {
  stop("No complete observed network replicate set found for years: ",
       paste(current_years, collapse = ", "))
}
if (!identical(
  sort(unique(as.integer(network_replicates))),
  sort(unique(as.integer(requested_network_replicates)))
)) {
  stop(
    "The observed baseline input contains an incomplete or stale ",
    "network-replicate set."
  )
}

expected_ids <- unlist(lapply(network_replicates, function(rep_idx) {
  rvf_with_replicate(scenario_base_ids, rep_idx)
}), use.names = FALSE)
missing_ids <- setdiff(expected_ids, scenario_index$scenario_id)
if (length(missing_ids) > 0) {
  stop("Missing observed scenario bundle(s) for replicate-aware baseline run: ",
       paste(missing_ids, collapse = ", "))
}

run_tag <- rvf_baseline_run_tag(sim)

# -------- Vector time series ---------------------------------------------
#' Aggregate vector S/E/I across sites to mean vector count and mean
#' infectious fraction per site type (village / grazing / watering).
compute_vector_summary <- function(output_matrix, params, t_seed = 0) {
  # +1L because deSolve places time in column 1.
  idx <- params$idx
  S_v_cols <- idx$S_m_v + 1L; E_v_cols <- idx$E_m_v + 1L; I_v_cols <- idx$I_m_v + 1L
  S_g_cols <- idx$S_m_g + 1L; E_g_cols <- idx$E_m_g + 1L; I_g_cols <- idx$I_m_g + 1L
  S_w_cols <- idx$S_m_w + 1L; E_w_cols <- idx$E_m_w + 1L; I_w_cols <- idx$I_m_w + 1L

  M_v <- output_matrix[, S_v_cols] + output_matrix[, E_v_cols] + output_matrix[, I_v_cols]
  M_g <- output_matrix[, S_g_cols] + output_matrix[, E_g_cols] + output_matrix[, I_g_cols]
  M_w <- output_matrix[, S_w_cols] + output_matrix[, E_w_cols] + output_matrix[, I_w_cols]
  I_v <- output_matrix[, I_v_cols]
  I_g <- output_matrix[, I_g_cols]
  I_w <- output_matrix[, I_w_cols]

  # Floor avoids divide-by-near-zero at dry-season sites.
  inf_frac_v <- I_v / pmax(M_v, 1e-6)
  inf_frac_g <- I_g / pmax(M_g, 1e-6)
  inf_frac_w <- I_w / pmax(M_w, 1e-6)

  tibble(
    time = output_matrix[, 1],
    step_label = vapply(
      output_matrix[, 1],
      time_label_from_time,
      character(1),
      time_labels = params$time_labels,
      days_per_step = params$days_per_step
    ),
    year_since_seed = (output_matrix[, 1] - t_seed) / 365,
    village_vector_mean    = rowMeans(M_v),
    village_vector_per_host = rowMeans(M_v) / params$N_h_const,
    grazing_vector_mean    = rowMeans(M_g),
    watering_vector_mean   = rowMeans(M_w),
    village_inf_frac_mean  = rowMeans(inf_frac_v) * 100,
    grazing_inf_frac_mean  = rowMeans(inf_frac_g) * 100,
    watering_inf_frac_mean = rowMeans(inf_frac_w) * 100
  )
}

run_one_replicate <- function(rep_idx) {
  # Keep the network replicate fixed across all component years in the
  # sequence so the burn-in / endemic run represents one coherent network
  # realisation rather than a mixture of replicate draws.
  scenario_ids <- rvf_with_replicate(scenario_base_ids, rep_idx)
  scenario_label <- sprintf("OBSERVED_%d_%d_CYCLE", min(current_years), max(current_years))
  scenario <- build_rvf_sequence_scenario(
    scenario_ids,
    cfg = cfg,
    validate = TRUE,
    scenario_type = "observed_current_cycle",
    sequence_label = scenario_label
  )
  scenario <- rvf_apply_network_mode(scenario, network_mode = sim$network_mode)

  calibrated <- rvf_adjust_calibrated_params(
    transmission_scale = sim$transmission_scale,
    capacity_scale = sim$capacity_scale,
    site_profile = sim$site_profile,
    rainfall_model = sim$rainfall_model
  )

  params <- rvf_assemble_params(
    scenario = scenario,
    calibrated = calibrated,
    use_temp_dependence = TRUE,
    p_mobile = sim$p_mobile,
    waning_years = sim$waning_years,
    external_introduction_per_year = sim$external_introduction_per_year,
    verbose = FALSE
  )

  res <- rvf_run_simulation(
    params = params,
    years_burn = sim$years_burn,
    years_endemic = sim$years_endemic,
    output_step = sim$output_step,
    n_seed = sim$n_seed,
    seed_village = sim$seed_village,
    init_state = NULL,
    verbose = FALSE
  )
  annual <- rvf_annual_metrics(res$output_endemic, params, t_seed = res$t_seed)

  ts <- res$time_series |>
    mutate(
      scenario = scenario_label,
      scenario_id = if (length(scenario_ids) == 1) scenario_ids else paste(scenario_ids, collapse = ";"),
      scenario_base_id = if (length(scenario_base_ids) == 1) scenario_base_ids else paste(scenario_base_ids, collapse = ";"),
      network_replicate = rep_idx,
      network_mode = sim$network_mode
    )
  annual_system <- annual$system_year |>
    mutate(
      scenario = scenario_label,
      scenario_id = if (length(scenario_ids) == 1) scenario_ids else paste(scenario_ids, collapse = ";"),
      scenario_base_id = if (length(scenario_base_ids) == 1) scenario_base_ids else paste(scenario_base_ids, collapse = ";"),
      network_replicate = rep_idx,
      network_mode = sim$network_mode,
      p_mobile = sim$p_mobile,
      transmission_scale = sim$transmission_scale,
      capacity_scale = sim$capacity_scale,
      site_profile = calibrated$site_profile,
      site_profile_label = calibrated$site_profile_label,
      rainfall_model = sim$rainfall_model,
      rainfall_label = rvf_rainfall_model_label(sim$rainfall_model, cfg = cfg)
    )

  late <- ts |> filter(year_since_seed >= (sim$years_endemic - 5))
  # Saved for optional baseline-state diagnostics; the headline mobility and
  # future analyses cold-start their own scenarios.
  baseline_state <- unname(res$output_endemic[nrow(res$output_endemic), -1])
  names(baseline_state) <- colnames(res$output_endemic)[-1]

  vec_ts <- compute_vector_summary(res$output_endemic, params, t_seed = res$t_seed) |>
    mutate(
      scenario_id = if (length(scenario_ids) == 1) scenario_ids else paste(scenario_ids, collapse = ";"),
      scenario_base_id = if (length(scenario_base_ids) == 1) scenario_base_ids else paste(scenario_base_ids, collapse = ";"),
      network_replicate = rep_idx,
      network_mode = sim$network_mode
    )
  late_vec <- vec_ts |> filter(year_since_seed >= (sim$years_endemic - 5))
  final_analysis_year <- max(annual_system$analysis_year, na.rm = TRUE)
  annual_last5 <- annual_system |>
    filter(.data$analysis_year >= final_analysis_year - 4L)

  summary <- tibble(
    scenario = scenario_label,
    scenario_id = if (length(scenario_ids) == 1) scenario_ids else paste(scenario_ids, collapse = ";"),
    scenario_base_id = if (length(scenario_base_ids) == 1) scenario_base_ids else paste(scenario_base_ids, collapse = ";"),
    network_replicate = rep_idx,
    network_mode = sim$network_mode,
    years_burn = sim$years_burn,
    years_endemic = sim$years_endemic,
    p_mobile = sim$p_mobile,
    transmission_scale = sim$transmission_scale,
    capacity_scale = sim$capacity_scale,
    site_profile = calibrated$site_profile,
    site_profile_label = calibrated$site_profile_label,
    m_max_v = calibrated$m_max_v,
    m_max_g = calibrated$m_max_g,
    m_max_w = calibrated$m_max_w,
    rainfall_model = sim$rainfall_model,
    rainfall_label = rvf_rainfall_model_label(sim$rainfall_model, cfg = cfg),
    waning_years = as.numeric(sim$waning_years),
    seed_village = sim$seed_village,
    n_seed = sim$n_seed,
    external_introduction_per_year = sim$external_introduction_per_year,
    final_seroprevalence = tail(ts$seroprevalence, 1),
    mean_seroprevalence_last5y = mean(late$seroprevalence, na.rm = TRUE),
    sd_seroprevalence_last5y = sd(late$seroprevalence, na.rm = TRUE),
    cv_seroprevalence_last5y = sd(late$seroprevalence, na.rm = TRUE) /
      mean(late$seroprevalence, na.rm = TRUE),
    final_prevalence = tail(ts$prevalence, 1),
    mean_prevalence_last5y = mean(late$prevalence, na.rm = TRUE),
    annual_new_infections_last5y = mean(annual_last5$annual_new_infections, na.rm = TRUE),
    cum_incidence_last5y = sum(annual_last5$annual_new_infections, na.rm = TRUE),
    mean_village_vector_last5y = mean(late_vec$village_vector_mean, na.rm = TRUE),
    mean_village_inf_frac_last5y = mean(late_vec$village_inf_frac_mean, na.rm = TRUE)
  )

  list(
    summary = summary,
    annual_system = annual_system,
    time_series = ts,
    vector_time_series = vec_ts,
    run = res,
    baseline_bundle = list(
      network_replicate = rep_idx,
      scenario_ids = scenario_ids,
      baseline_state = baseline_state,
      params = params
    )
  )
}

replicate_jobs <- seq_along(network_replicates)
n_cores <- rvf_max_cores(
  default = cfg$execution$max_cores,
  n_tasks = length(replicate_jobs)
)
save_full_run_rds <- isTRUE(cfg$execution$output$save_baseline_run_rds)
message(sprintf("Running endemic baseline across %d replicate(s) on %d core(s).",
                length(replicate_jobs), n_cores))

# Reuse cached endemic runs only when a full run RDS was deliberately saved.
# Normal publication runs save the compact baseline-state bundle instead.
endemic_runs_cache <- file.path(intermediate_dir, paste0(run_tag, "_run.rds"))
if (isTRUE(cfg$execution$reuse$observed_baseline) &&
    file.exists(endemic_runs_cache)) {
  message(sprintf(
    "Reusing cached endemic runs from %s (no re-simulation).",
    endemic_runs_cache
  ))
  results <- readRDS(endemic_runs_cache)
} else {
  results <- rvf_parallel_lapply(replicate_jobs, function(i) {
    rep_idx <- network_replicates[i]
    message(sprintf("Running endemic baseline %d/%d: replicate %03d",
                    i, length(network_replicates), rep_idx))
    run_one_replicate(rep_idx)
  }, n_cores = n_cores)
}

summary_df <- bind_rows(lapply(results, `[[`, "summary")) |>
  arrange(network_replicate)
annual_system_all <- bind_rows(lapply(results, `[[`, "annual_system")) |>
  arrange(network_replicate, analysis_year)
ts_all <- bind_rows(lapply(results, `[[`, "time_series"))
vec_ts_all <- bind_rows(lapply(results, `[[`, "vector_time_series"))

# A periodically forced system need not converge to a constant state. Use four
# 10-year cycles for the final 40-year reporting block and compare it with the
# preceding four cycles at matching phases. This gives the observed baseline the
# same reporting duration as the two-cycle future summaries.
cycle_years <- length(current_years)
reporting_years <- as.integer(cfg$execution$reporting_years)
if (reporting_years %% cycle_years != 0L) {
  stop("The reporting duration must contain complete observed forcing cycles.")
}
terminal_cycles <- reporting_years %/% cycle_years
cycle_diagnostics <- rvf_build_cycle_diagnostics(
  annual_system = annual_system_all,
  time_series = ts_all,
  group_cols = c("scenario", "network_replicate", "network_mode"),
  forcing_cycle_years = cycle_years,
  years_endemic = sim$years_endemic,
  terminal_cycles = terminal_cycles
)
complete_observed_cycles <- cycle_diagnostics$protocol$complete_cycles
recurrent_state_summary <- cycle_diagnostics$convergence_summary
periodicity_summary <- cycle_diagnostics$periodicity_by_run |>
  count(.data$recurrent_pattern, name = "n_network_replicates")

# Keep a compact aggregate table for quick inspection, but preserve the
# replicate-level CSV because that is the uncertainty-carrying output.
summary_aggregate <- summary_df |>
  summarise(
    n_network_replicates = n(),
    mean_seroprevalence_last5y_median = median(mean_seroprevalence_last5y, na.rm = TRUE),
    mean_seroprevalence_last5y_min = min(mean_seroprevalence_last5y, na.rm = TRUE),
    mean_seroprevalence_last5y_max = max(mean_seroprevalence_last5y, na.rm = TRUE),
    mean_prevalence_last5y_median = median(mean_prevalence_last5y, na.rm = TRUE),
    annualized_note = "raw replicate-level summaries are in the companion *_summary.csv file"
  ) |>
  mutate(
    scenario = sprintf("OBSERVED_%d_%d_CYCLE", min(current_years), max(current_years)),
    baseline_years = paste(current_years, collapse = ","),
    network_mode = sim$network_mode,
    years_burn = sim$years_burn,
    years_endemic = sim$years_endemic,
    p_mobile = sim$p_mobile,
    transmission_scale = sim$transmission_scale,
    capacity_scale = sim$capacity_scale,
    site_profile = sim$site_profile,
    site_profile_label = rvf_site_profile_label(sim$site_profile, cfg = cfg),
    rainfall_model = sim$rainfall_model,
    rainfall_label = rvf_rainfall_model_label(sim$rainfall_model, cfg = cfg),
    external_introduction_per_year = sim$external_introduction_per_year,
    waning_label = rvf_waning_label(sim$waning_years),
    complete_observed_cycles = complete_observed_cycles,
    recurrent_state_supported = recurrent_state_summary$recurrent_state_supported[[1]],
    recurrent_state_supported_fraction =
      recurrent_state_summary$recurrent_state_supported_fraction[[1]]
  )

print(summary_df)
cat("\nComplete-cycle recurrent-state check:\n")
print(recurrent_state_summary)
cat("\nPhase-matched recurrent pattern by network replicate:\n")
print(periodicity_summary)

# -------- Supplementary equilibrium figure (3 panels) --------------------
# Publication-facing baseline figure. Off-protocol summaries (the phase-sensitive
# last-5-year reference line) and the village-only vector panels are dropped;
# bands are the median and 10th-90th percentile ACROSS network replicates.
report_window <- cycle_diagnostics$protocol$terminal_window_years
final_start   <- sim$years_endemic - report_window
prev_start    <- final_start - report_window
cycle_breaks  <- seq(cycle_years, sim$years_endemic, by = cycle_years)
pctl <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE))

# Panel 1 -- full seroprevalence trajectory: initial transient + recurrent
# state, with the final 40-year reporting window shaded.
sero_traj <- ts_all |>
  group_by(year_since_seed) |>
  summarise(mid = median(seroprevalence, na.rm = TRUE),
            lo  = pctl(seroprevalence, 0.10),
            hi  = pctl(seroprevalence, 0.90), .groups = "drop")
p_traj <- ggplot(sero_traj, aes(year_since_seed, mid)) +
  annotate("rect", xmin = final_start, xmax = sim$years_endemic,
           ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.55) +
  geom_vline(xintercept = cycle_breaks, colour = "grey88", linewidth = 0.2) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue3", alpha = 0.22) +
  geom_line(linewidth = 0.6, colour = "steelblue4") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Years since seeding", y = "Seroprevalence (%)")

# Panel 2 -- recurrent-state comparison: overlay the preceding and final
# reporting blocks at matching phases. Each 40-year block contains four full
# observed-climate cycles and two repetitions of a possible period-two pattern.
block_labels <- c("Preceding 40-year block", "Final 40-year block")
block_overlay <- ts_all |>
  filter(year_since_seed >= prev_start,
         year_since_seed < sim$years_endemic) |>
  mutate(is_final = year_since_seed >= final_start,
         block = factor(if_else(is_final, block_labels[2], block_labels[1]), levels = block_labels),
         year_in_block = year_since_seed - if_else(is_final, final_start, prev_start)) |>
  group_by(block, year_in_block) |>
  summarise(mid = median(seroprevalence, na.rm = TRUE),
            lo = pctl(seroprevalence, 0.10), hi = pctl(seroprevalence, 0.90), .groups = "drop")
p_overlay <- ggplot(block_overlay, aes(year_in_block, mid, colour = block, fill = block)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.14, colour = NA) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("grey55", "steelblue4")) +
  scale_fill_manual(values = c("grey55", "steelblue4")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.title = element_blank()) +
  labs(x = "Year within 40-year block", y = "Seroprevalence (%)")

# Panel 3 -- terminal seasonal pattern. The shared diagnostic helper first
# computes time-weighted monthly means within each network replicate; the band
# below then represents variation among network replicates.
seasonal_host <- cycle_diagnostics$seasonal_by_run |>
  mutate(month = factor(.data$month, levels = month.abb)) |>
  group_by(month) |>
  summarise(mid = median(.data$mean_prevalence, na.rm = TRUE),
            lo = pctl(.data$mean_prevalence, 0.10),
            hi = pctl(.data$mean_prevalence, 0.90),
            .groups = "drop")
p_season <- ggplot(seasonal_host, aes(month, mid, group = 1)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "firebrick2", alpha = 0.15) +
  geom_line(colour = "firebrick3", linewidth = 0.8) +
  geom_point(colour = "firebrick3", size = 1.8) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "Host prevalence (%)")

if (rvf_write_diagnostic_plots()) {
  equilibrium_fig <- (p_traj / p_overlay / p_season) +
    plot_layout(heights = c(1.15, 1, 1))
  diagnostic_dir <- file.path(cfg$figures_dir, "diagnostics")
  dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    file.path(diagnostic_dir, paste0(run_tag, "_equilibrium.png")),
    equilibrium_fig,
    width = 8,
    height = 10,
    dpi = 300
  )
}

# Optional village-vector seasonal diagnostic (not a manuscript panel; villages
# only). Built and written only on request to avoid redundant outputs.
if (rvf_write_diagnostic_plots()) {
  baseline_params <- results[[1]]$baseline_bundle$params
  forcing_months <- baseline_params$time_labels
  vec_forcing <- tibble(
    step_label = factor(forcing_months, levels = forcing_months),
    rainfall_village_mm = vapply(forcing_months,
      function(m) mean(baseline_params$rainfall_vil[[m]], na.rm = TRUE), numeric(1)))
  seasonal_vector <- vec_ts_all |>
    filter(year_since_seed >= final_start) |>
    mutate(step_label = factor(step_label, levels = forcing_months)) |>
    group_by(step_label) |>
    summarise(vph_mid = median(village_vector_per_host, na.rm = TRUE),
              vph_lo = pctl(village_vector_per_host, 0.10),
              vph_hi = pctl(village_vector_per_host, 0.90), .groups = "drop") |>
    left_join(vec_forcing, by = "step_label")
  p_dv <- ggplot(seasonal_vector, aes(step_label, vph_mid, group = 1)) +
    geom_ribbon(aes(ymin = vph_lo, ymax = vph_hi), fill = "darkgreen", alpha = 0.14) +
    geom_line(colour = "darkgreen", linewidth = 0.8) + geom_point(colour = "darkgreen", size = 1.8) +
    theme_bw(base_size = 11) + labs(x = NULL, y = "Village vectors per host")
  p_dr <- ggplot(seasonal_vector, aes(step_label, rainfall_village_mm)) +
    geom_col(fill = "steelblue3", alpha = 0.75) +
    theme_bw(base_size = 11) + labs(x = NULL, y = "Village rainfall (mm)")
  ggsave(file.path(out_dir, paste0(run_tag, "_seasonal_vector_diagnostic.png")),
         p_dv / p_dr, width = 9, height = 7, dpi = 150)
} else {
  message("Skipping the optional village-vector seasonal diagnostic (see execution$output in R/config.R).")
}

# -------- Persist outputs ------------------------------------------------
write_derived_csvs <- rvf_write_derived_csvs()
if (write_derived_csvs) {
  write_csv(ts_all,          file.path(out_dir, paste0(run_tag, "_timeseries.csv")))
  write_csv(vec_ts_all,      file.path(out_dir, paste0(run_tag, "_vector_timeseries.csv")))
} else {
  message("Skipping bulky derived baseline CSVs (see execution$output in R/config.R).")
}
write_csv(summary_df,          file.path(out_dir, paste0(run_tag, "_summary.csv")))
write_csv(annual_system_all,   file.path(out_dir, paste0(run_tag, "_annual_system.csv")))
write_csv(summary_aggregate,   file.path(out_dir, paste0(run_tag, "_summary_aggregate.csv")))
write_csv(
  recurrent_state_summary,
  file.path(out_dir, paste0(run_tag, "_recurrent_state_summary.csv"))
)
saveRDS(
  cycle_diagnostics,
  file.path(intermediate_dir, paste0(run_tag, "_cycle_diagnostics.rds"))
)
if (save_full_run_rds) {
  saveRDS(results, file.path(intermediate_dir, paste0(run_tag, "_run.rds")))
} else {
  message("Skipping the full baseline run RDS (see execution$output in R/config.R).")
}

# Baseline state bundle: final state + sim config for downstream rehydration.
# Scalar fields are duplicated out of `sim` for backward compatibility with
# downstream scripts that expect them at the top level.
if (isTRUE(cfg$execution$output$save_baseline_state)) {
  saveRDS(
    list(
      baseline_states = setNames(lapply(results, `[[`, "baseline_bundle"),
                                 sprintf("rep_%03d", network_replicates)),
      sim = sim,
      baseline_years = current_years,
      network_replicates = network_replicates,
      network_mode = sim$network_mode,
      years_burn = sim$years_burn,
      years_endemic = sim$years_endemic,
      complete_observed_cycles = complete_observed_cycles,
      recurrent_state_summary = recurrent_state_summary,
      p_mobile = sim$p_mobile,
      transmission_scale = sim$transmission_scale,
      capacity_scale = sim$capacity_scale,
      site_profile = sim$site_profile,
      site_profile_label = rvf_site_profile_label(sim$site_profile, cfg = cfg),
      rainfall_model = sim$rainfall_model,
      rainfall_label = rvf_rainfall_model_label(sim$rainfall_model, cfg = cfg),
      external_introduction_per_year = sim$external_introduction_per_year
    ),
    file.path(intermediate_dir, paste0(run_tag, "_baseline_state.rds"))
  )
} else {
  message("Skipping optional baseline-state bundles (see execution$output in R/config.R).")
}

cat("Observed-cycle endemic run complete.\n")
