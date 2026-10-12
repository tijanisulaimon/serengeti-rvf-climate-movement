#!/usr/bin/env Rscript
#' Quick check that the transmission model runs, in about 5 seconds.
#'
#' Runs one short observed scenario to confirm that scenario loading,
#' parameter assembly, ODE integration, and output writing still work.

suppressPackageStartupMessages({
  library(readr)
})

source(here::here("analysis", "03_rvf_transmission", "R", "transmission_model.R"))

cfg <- load_rvf_config()
sim <- cfg$sim

# Check the shared force-of-infection helper against the equations that were
# previously implemented independently in the solver and diagnostics. The
# boundary p_mobile = 0 and the reference p_mobile = 1 are included because
# they use different resource-site branches.
N_h_test <- c(100, 120, 80)
I_h_test <- c(4, 9, 2)
I_m_v_test <- c(20, 15, 8)
I_m_g_test <- c(30, 10)
I_m_w_test <- c(12, 18)
G_test <- rbind(c(1, 0), c(0.25, 0.75), c(0, 1))
W_test <- rbind(c(0.5, 0.5), c(1, 0), c(0.2, 0.8))
a_v_test <- c(0.20, 0.21, 0.22)
a_g_test <- c(0.19, 0.23)
a_w_test <- c(0.18, 0.24)
time_test <- list(t_v = 15 / 24, t_g = 8 / 24, t_w = 1 / 24)
alpha_vh_test <- 0.4
alpha_hv_test <- 0.6
external_test <- c(0, 1e-5, 0)

for (p_test in c(0, 0.4, 1)) {
  observed <- rvf_compute_foi_components(
    N_h = N_h_test,
    I_h = I_h_test,
    I_m_v = I_m_v_test,
    I_m_g = I_m_g_test,
    I_m_w = I_m_w_test,
    G = G_test,
    W = W_test,
    a_v = a_v_test,
    a_g = a_g_test,
    a_w = a_w_test,
    time_budget = time_test,
    p_mobile = p_test,
    alpha_vh = alpha_vh_test,
    alpha_hv = alpha_hv_test,
    external_hazard = external_test
  )

  N_g <- as.vector(t(G_test) %*% (p_test * N_h_test))
  I_g <- as.vector(t(G_test) %*% (p_test * I_h_test))
  N_w <- as.vector(t(W_test) %*% (p_test * N_h_test))
  I_w <- as.vector(t(W_test) %*% (p_test * I_h_test))
  village_foi <- a_v_test * alpha_vh_test * I_m_v_test / N_h_test

  if (p_test == 0) {
    grazing_foi <- rep(0, length(N_h_test))
    watering_foi <- rep(0, length(N_h_test))
    grazing_vector_foi <- rep(0, length(N_g))
    watering_vector_foi <- rep(0, length(N_w))
  } else {
    grazing_site_foi <- a_g_test * time_test$t_g * alpha_vh_test *
      I_m_g_test / N_g
    watering_site_foi <- a_w_test * time_test$t_w * alpha_vh_test *
      I_m_w_test / N_w
    grazing_foi <- p_test * as.vector(G_test %*% grazing_site_foi)
    watering_foi <- p_test * as.vector(W_test %*% watering_site_foi)
    grazing_vector_foi <- a_g_test * p_test * time_test$t_g *
      alpha_hv_test * I_g / N_g
    watering_vector_foi <- a_w_test * p_test * time_test$t_w *
      alpha_hv_test * I_w / N_w
  }

  expected_village <- ((1 - p_test) + p_test * time_test$t_v) * village_foi
  expected_total <- expected_village + grazing_foi + watering_foi + external_test
  expected_village_vector <- a_v_test *
    ((1 - p_test) + p_test * time_test$t_v) *
    alpha_hv_test * I_h_test / N_h_test

  stopifnot(
    isTRUE(all.equal(observed$lambda_h, expected_total, tolerance = 1e-12)),
    isTRUE(all.equal(observed$lambda_h_village, expected_village, tolerance = 1e-12)),
    isTRUE(all.equal(observed$lambda_h_grazing, grazing_foi, tolerance = 1e-12)),
    isTRUE(all.equal(observed$lambda_h_watering, watering_foi, tolerance = 1e-12)),
    isTRUE(all.equal(
      observed$lambda_hv_village,
      expected_village_vector,
      tolerance = 1e-12
    )),
    isTRUE(all.equal(
      observed$lambda_hv_grazing,
      grazing_vector_foi,
      tolerance = 1e-12
    )),
    isTRUE(all.equal(
      observed$lambda_hv_watering,
      watering_vector_foi,
      tolerance = 1e-12
    ))
  )
}

stopifnot(identical(formals(rvf_assemble_params)$capacity_reference_p_mobile, 1))
cat("FOI and capacity-reference checks passed.\n")

# Smoke-test the first year of the configured CHELSA historical baseline.
smoke_year <- min(sim$baseline_years %||% sim$baseline_year)
smoke_index <- rvf_select_scenarios(
  dataset = "observed",
  year = smoke_year,
  network_replicate = cfg$representative_replicate
)
if (nrow(smoke_index) == 0) {
  stop("No scenario found for the quick check: observed ", smoke_year, " replicate ", cfg$representative_replicate,
       ". Rebuild scenario bundles first.")
}
scenario_id <- smoke_index$scenario_id[1]

out_dir <- file.path(cfg$cache_dir, "smoke_tests")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Fresh cold-start run with the current default parameter block.
scenario <- load_rvf_scenario(scenario_id)
params <- rvf_assemble_params(
  scenario = scenario,
  use_temp_dependence = TRUE,
  p_mobile = sim$p_mobile,
  waning_years = sim$waning_years,
  verbose = TRUE
)

res <- rvf_run_simulation(
  params = params,
  years_burn = 1,
  years_endemic = 2,
  output_step = 7,
  n_seed = 25,
  seed_village = 10,
  init_state = NULL,
  verbose = TRUE
)

# Numerical and accounting invariants for the complete ODE output. These checks
# guard the quantities used by every downstream burden and convergence summary.
out <- as.matrix(res$output_endemic)
state <- out[, -1, drop = FALSE]
if (any(!is.finite(state))) {
  stop("Smoke-test ODE output contains non-finite state values.")
}
if (min(state) < -1e-5) {
  stop("Smoke-test ODE output contains materially negative state values.")
}

state_columns <- function(indices) indices + 1L
host_total <- out[, state_columns(params$idx$S_h), drop = FALSE] +
  out[, state_columns(params$idx$E_h), drop = FALSE] +
  out[, state_columns(params$idx$I_h), drop = FALSE] +
  out[, state_columns(params$idx$R_h), drop = FALSE]
host_target <- matrix(
  params$N_h_const,
  nrow = nrow(host_total),
  ncol = params$n_vil
)
if (max(abs(host_total - host_target)) > 1e-3) {
  stop("Host population is not conserved within the expected ODE tolerance.")
}

for (diagnostic_indices in list(params$idx$C_h, params$idx$H_h)) {
  diagnostic_state <- out[, state_columns(diagnostic_indices), drop = FALSE]
  if (min(apply(diagnostic_state, 2, diff)) < -1e-6) {
    stop("A cumulative incidence or hazard state decreases during integration.")
  }
}

capacity_values <- unlist(c(params$K_v, params$K_g, params$K_w))
trait_values <- unlist(c(
  params$epsilon_m_v, params$epsilon_m_g, params$epsilon_m_w,
  params$mu_m_v, params$mu_m_g, params$mu_m_w,
  params$a_m_v, params$a_m_g, params$a_m_w
))
if (any(!is.finite(capacity_values)) || any(capacity_values <= 0)) {
  stop("Climate-dependent vector carrying capacities are not finite and positive.")
}
if (any(!is.finite(trait_values)) || any(trait_values <= 0)) {
  stop("Temperature-dependent vector rates are not finite and positive.")
}

recruitment_check <- rvf_adult_recruitment(
  M = c(0, 50, 100, 150),
  K = rep(100, 4),
  r_m = params$r_m
)
if (any(recruitment_check < 0) ||
    recruitment_check[[3]] != 0 ||
    recruitment_check[[4]] != 0) {
  stop("Bounded adult recruitment does not stop at or above carrying capacity.")
}

annual <- rvf_annual_metrics(out, params, t_seed = res$t_seed)
cumulative_incidence <- out[, state_columns(params$idx$C_h), drop = FALSE]
integrated_total <- sum(
  cumulative_incidence[nrow(cumulative_incidence), ] -
    cumulative_incidence[1, ]
)
annual_total <- sum(annual$system_year$annual_new_infections)
if (!isTRUE(all.equal(integrated_total, annual_total, tolerance = 1e-8))) {
  stop("Annual infection summaries do not reconcile with cumulative incidence.")
}
cat("ODE state and annual-accounting checks passed.\n")

# Save both the full object and the compact weekly summary.
saveRDS(res, file.path(out_dir, paste0(scenario_id, "_smoke_test.rds")))
if (rvf_write_derived_csvs()) {
  write_csv(res$time_series, file.path(out_dir, paste0(scenario_id, "_smoke_test_timeseries.csv")))
} else {
  message("Skipping the quick-check time-series CSV. The full result is saved in the RDS.")
}

cat("Smoke test complete for ", scenario_id, "\n", sep = "")
