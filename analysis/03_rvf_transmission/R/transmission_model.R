#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(deSolve)
  library(dplyr)
})

source(here::here("analysis", "03_rvf_transmission", "R", "model_parameters.R"))

# Minimum number of hosts at a resource site below which that site is
# treated as depopulated for infection-pressure calculations (avoids
# divide-by-very-small near-zero mobile-host populations during burn-in).
MIN_HOSTS <- 1.0

rvf_nonnegative_state <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  pmax(x, 0)
}

rvf_safe_ratio <- function(num, den, min_den = 0) {
  den <- as.numeric(den)
  out <- rep(0, length(den))
  ok <- is.finite(den) & den >= min_den
  out[ok] <- rvf_nonnegative_state(num)[ok] / den[ok]
  out[!is.finite(out)] <- 0
  out
}

rvf_adult_recruitment <- function(M, K, r_m) {
  M <- rvf_nonnegative_state(M)
  K <- pmax(rvf_nonnegative_state(K), 1)
  r_m <- suppressWarnings(as.numeric(r_m)[1])
  if (!is.finite(r_m) || is.na(r_m) || r_m < 0) r_m <- 0

  # Carrying capacity limits emergence into the adult susceptible class.
  # When adult abundance exceeds the current habitat capacity, recruitment
  # stops rather than becoming negative; existing adults then decline through
  # temperature-dependent mortality.
  r_m * M * pmax(0, 1 - M / K)
}

rvf_ode_controls <- function(cfg = load_rvf_config()) {
  controls <- cfg$execution$ode
  values <- c(controls$rtol, controls$atol, controls$maxsteps)
  if (length(values) != 3L || any(!is.finite(values)) ||
      controls$rtol <= 0 || controls$atol <= 0 || controls$maxsteps < 1) {
    stop("ODE controls in R/config.R must be finite and strictly positive.")
  }

  # These tolerances are below the scale of reported serology and annual
  # burden summaries, but avoid excessive tiny-step LSODA work in long runs.
  list(
    rtol = as.numeric(controls$rtol),
    atol = as.numeric(controls$atol),
    maxsteps = as.integer(controls$maxsteps)
  )
}

# Function: compute the host and vector forces of infection for one model state.
# This helper is shared by the ODE and the reporting diagnostics so the
# decomposed FOI cannot drift away from the rates integrated by the solver.
# N_h / I_h: total and infectious cattle by home village.
# I_m_*: infectious adult vectors at villages, grazing sites, and water points.
# G / W: row-normalised village-to-resource matrices for the current month.
rvf_compute_foi_components <- function(N_h,
                                       I_h,
                                       I_m_v,
                                       I_m_g,
                                       I_m_w,
                                       G,
                                       W,
                                       a_v,
                                       a_g,
                                       a_w,
                                       time_budget,
                                       p_mobile,
                                       alpha_vh,
                                       alpha_hv,
                                       external_hazard = NULL) {
  n_vil <- length(N_h)
  n_gra <- length(I_m_g)
  n_wat <- length(I_m_w)

  p_mob <- suppressWarnings(as.numeric(p_mobile)[1])
  if (!is.finite(p_mob)) p_mob <- 0
  p_mob <- max(0, min(1, p_mob))

  # Aggregate the mobile sub-herd at each selected resource site. Multiplying
  # both total and infectious cattle by p_mobile preserves prevalence while
  # retaining the correct site-level vector-to-host denominator.
  if (p_mob < 1e-10) {
    N_h_at_gra <- rep(0, n_gra)
    I_h_at_gra <- rep(0, n_gra)
    N_h_at_wat <- rep(0, n_wat)
    I_h_at_wat <- rep(0, n_wat)
  } else {
    N_h_mobile <- p_mob * N_h
    I_h_mobile <- p_mob * I_h
    N_h_at_gra <- as.vector(t(G) %*% N_h_mobile)
    I_h_at_gra <- as.vector(t(G) %*% I_h_mobile)
    N_h_at_wat <- as.vector(t(W) %*% N_h_mobile)
    I_h_at_wat <- as.vector(t(W) %*% I_h_mobile)
  }

  prev_vil <- rvf_safe_ratio(I_h, N_h, min_den = 1)
  prev_gra <- rvf_safe_ratio(I_h_at_gra, N_h_at_gra, min_den = MIN_HOSTS)
  prev_wat <- rvf_safe_ratio(I_h_at_wat, N_h_at_wat, min_den = MIN_HOSTS)

  # Vector-to-host FOI is first calculated for a host spending the full day at
  # a site, then weighted by the mobile fraction and daily residence times.
  lambda_vh_vil <- a_v * alpha_vh *
    rvf_safe_ratio(I_m_v, N_h, min_den = 1)
  if (p_mob > 1e-10) {
    lambda_vh_gra_site <- a_g * time_budget$t_g * alpha_vh *
      rvf_safe_ratio(I_m_g, N_h_at_gra, min_den = MIN_HOSTS)
    lambda_vh_wat_site <- a_w * time_budget$t_w * alpha_vh *
      rvf_safe_ratio(I_m_w, N_h_at_wat, min_den = MIN_HOSTS)
    lambda_vh_gra_mobile <- as.vector(G %*% lambda_vh_gra_site)
    lambda_vh_wat_mobile <- as.vector(W %*% lambda_vh_wat_site)
  } else {
    lambda_vh_gra_mobile <- rep(0, n_vil)
    lambda_vh_wat_mobile <- rep(0, n_vil)
  }

  lambda_h_village <-
    (1 - p_mob) * lambda_vh_vil +
    p_mob * time_budget$t_v * lambda_vh_vil
  lambda_h_grazing <- p_mob * lambda_vh_gra_mobile
  lambda_h_watering <- p_mob * lambda_vh_wat_mobile
  if (is.null(external_hazard)) external_hazard <- 0
  external_hazard <- as.numeric(external_hazard)
  if (length(external_hazard) == 1L) {
    external_hazard <- rep(external_hazard, n_vil)
  }
  if (length(external_hazard) != n_vil || any(!is.finite(external_hazard))) {
    stop("external_hazard must be finite and have length 1 or n_vil.")
  }
  lambda_h_external <- pmax(0, external_hazard)
  lambda_h_total <- pmax(
    0,
    lambda_h_village + lambda_h_grazing +
      lambda_h_watering + lambda_h_external
  )

  # Host-to-vector FOI uses the herd-average fraction of each day during which
  # cattle are present at a site and the infectious-cattle prevalence there.
  q_v <- (1 - p_mob) + p_mob * time_budget$t_v
  q_g <- p_mob * time_budget$t_g
  q_w <- p_mob * time_budget$t_w

  list(
    p_mobile = p_mob,
    lambda_h = lambda_h_total,
    lambda_h_village = lambda_h_village,
    lambda_h_grazing = lambda_h_grazing,
    lambda_h_watering = lambda_h_watering,
    lambda_h_external = lambda_h_external,
    lambda_hv_village = a_v * q_v * alpha_hv * prev_vil,
    lambda_hv_grazing = a_g * q_g * alpha_hv * prev_gra,
    lambda_hv_watering = a_w * q_w * alpha_hv * prev_wat
  )
}

#' RHS of the endemic RVF ODE system.
#'
#' Integrated by deSolve::ode(). `t` is time in days since the simulation
#' start; `y` is the flattened state vector (9 states per village + 3
#' states per grazing site + 3 states per watering site -- see
#' `rvf_get_state_indices()` in `R/model_parameters.R` for the layout);
#' `parms` is the full precomputed parameter bundle produced by
#' `rvf_assemble_params()`.
#'
#' The function evaluates month-indexed climate forcing (K, mu_m, eps_m,
#' and biting rate), reconstructs the four compartmental populations (host SEIR +
#' vector SEI at each of three site types), computes the two forces of
#' infection (host->vector, vector->host), then returns the time
#' derivative vector dy.
rvf_model_endemic <- function(t, y, parms) {
  with(parms, {
    # -------- 1. Pull month-indexed forcing at current time ----------
    # `time_label_from_time()` converts absolute simulation time (days)
    # into the month label (Jan..Dec) under which the monthly-resolution
    # climate and network data are stored. This is what makes the
    # climate forcing seasonal without needing explicit Fourier terms.
    month <- time_label_from_time(t, time_labels, days_per_step)

    # Carrying capacities at villages / grazing / watering for this
    # month (vectors of length n_vil / n_gra / n_wat). Driven by
    # the configured monthly-rainfall response family in the parameter
    # bundle (`hill` baseline or named rainfall-shape sensitivities).
    K_v_t <- K_v[[month]]
    K_g_t <- K_g[[month]]
    K_w_t <- K_w[[month]]

    # Row-normalised village->grazing and village->watering commute
    # matrices for this month (n_vil x n_gra, n_vil x n_wat).
    G <- G_list[[month]]
    W <- W_list[[month]]

    # EIP rate (per day) at each site type for this month (Cecilia
    # Table 1, Barker 2013 linear-in-T form). Higher T -> faster
    # extrinsic incubation.
    eps_m_v <- epsilon_m_v[[month]]
    eps_m_g <- epsilon_m_g[[month]]
    eps_m_w <- epsilon_m_w[[month]]

    # Adult vector mortality (per day), Cecilia Table 1 row i=1:
    # quadratic in T, U-shaped with minimum around 22 degC.
    mu_v <- mu_m_v[[month]]
    mu_g <- mu_m_g[[month]]
    mu_w <- mu_m_w[[month]]

    # Biting rate (bites/vector/day). Temperature-dependent per Cecilia
    # Table 1: a(T) = (1 + 0.17) / g(T) where g(T) is the
    # gonotrophic-cycle duration. Vectors at each site type use their
    # own local temperature, so a_v / a_g / a_w are independent vectors.
    a_v <- a_m_v[[month]]
    a_g <- a_m_g[[month]]
    a_w <- a_m_w[[month]]
    tb_month <- rvf_time_budget_at(
      list(t_v = t_v, t_g = t_g, t_w = t_w),
      month
    )

    # -------- 2. Unpack state vector into compartments --------------
    # Host compartments by village. Tiny negative values can appear from the
    # adaptive solver near zero; clamp them before any rate calculation so they
    # cannot create negative prevalence or FOI.
    S_h <- rvf_nonnegative_state(y[idx$S_h])
    E_h <- rvf_nonnegative_state(y[idx$E_h])
    I_h <- rvf_nonnegative_state(y[idx$I_h])
    R_h <- rvf_nonnegative_state(y[idx$R_h])
    # Total host population per village. `pmax(..., 1)` avoids a
    # pathological divide-by-zero if a village somehow empties (should
    # not happen given births, but defensive).
    N_h <- pmax(S_h + E_h + I_h + R_h, 1)

    # Vector compartments at each village (SEI, no R for mosquitoes).
    S_m_v <- rvf_nonnegative_state(y[idx$S_m_v])
    E_m_v <- rvf_nonnegative_state(y[idx$E_m_v])
    I_m_v <- rvf_nonnegative_state(y[idx$I_m_v])
    M_v <- S_m_v + E_m_v + I_m_v

    # Vector compartments at each grazing site.
    S_m_g <- rvf_nonnegative_state(y[idx$S_m_g])
    E_m_g <- rvf_nonnegative_state(y[idx$E_m_g])
    I_m_g <- rvf_nonnegative_state(y[idx$I_m_g])
    M_g <- S_m_g + E_m_g + I_m_g

    # Vector compartments at each watering site.
    S_m_w <- rvf_nonnegative_state(y[idx$S_m_w])
    E_m_w <- rvf_nonnegative_state(y[idx$E_m_w])
    I_m_w <- rvf_nonnegative_state(y[idx$I_m_w])
    M_w <- S_m_w + E_m_w + I_m_w

    # -------- 3. Forces of infection under daily commuting ----------
    # The same helper is used after integration to decompose reported FOI.
    # Keeping one implementation makes the solver and diagnostics auditable.
    external_hazard <- if (exists("external_intro_hazard", inherits = FALSE)) {
      external_intro_hazard
    } else {
      rep(0, n_vil)
    }
    foi <- rvf_compute_foi_components(
      N_h = N_h,
      I_h = I_h,
      I_m_v = I_m_v,
      I_m_g = I_m_g,
      I_m_w = I_m_w,
      G = G,
      W = W,
      a_v = a_v,
      a_g = a_g,
      a_w = a_w,
      time_budget = tb_month,
      p_mobile = p_mobile,
      alpha_vh = alpha_vh,
      alpha_hv = alpha_hv,
      external_hazard = external_hazard
    )
    p_mob <- foi$p_mobile
    lambda_h <- foi$lambda_h
    lambda_hv_vil <- foi$lambda_hv_village
    lambda_hv_gra <- foi$lambda_hv_grazing
    lambda_hv_wat <- foi$lambda_hv_watering

    # -------- 4. Host SEIR ODE -------------------------------------
    # Constant per-village population driven by `N_h_const` birth rate
    # (so S + E + I + R fluctuates minimally around N_h_const).
    # `omega = 0` is the default (life-long immunity); if set to 1/
    # (waning_years * 365) the R compartment feeds back to S.
    births <- mu_h * N_h_const
    dS_h <- births + omega * R_h - lambda_h * S_h - mu_h * S_h
    dE_h <- lambda_h * S_h - epsilon_h * E_h - mu_h * E_h
    dI_h <- epsilon_h * E_h - gamma_h * I_h - mu_h * I_h
    dR_h <- gamma_h * I_h - omega * R_h - mu_h * R_h
    # Cumulative incidence and hazard diagnostics. `C_h` counts new
    # infections entering E (lambda_h * S_h), while `H_h` integrates the
    # per-susceptible host force of infection. Neither feeds back into the
    # dynamics, but both make annual FOI summaries independent of the
    # chosen ODE output interval.
    dC_h <- lambda_h * S_h
    dH_h <- lambda_h

    # -------- 5. Vector SEI ODE at villages -------------------------
    # Adult-only compression of immature-stage limitation: climate-driven K
    # caps new adult emergence into S_m, while existing adults die through
    # temperature-dependent mortality. This avoids imposing larval habitat
    # collapse as an extra mortality term on exposed/infectious adults.
    births_v <- rvf_adult_recruitment(M_v, K_v_t, r_m)
    dS_m_v <- births_v - lambda_hv_vil * S_m_v - mu_v * S_m_v
    dE_m_v <- lambda_hv_vil * S_m_v - eps_m_v * E_m_v - mu_v * E_m_v
    dI_m_v <- eps_m_v * E_m_v - mu_v * I_m_v

    if (p_mob < 1e-10) {
      # In sedentary mode, resource sites have no cattle exposure and cannot
      # feed back to hosts. Freeze their vector states at zero so arbitrary
      # network-replicate resource capacities cannot affect adaptive solver
      # step choices.
      dS_m_g <- rep(0, n_gra)
      dE_m_g <- rep(0, n_gra)
      dI_m_g <- rep(0, n_gra)
      dS_m_w <- rep(0, n_wat)
      dE_m_w <- rep(0, n_wat)
      dI_m_w <- rep(0, n_wat)
    } else {
      # -------- 6. Vector SEI ODE at grazing sites ------------------
      births_g <- rvf_adult_recruitment(M_g, K_g_t, r_m)
      dS_m_g <- births_g - lambda_hv_gra * S_m_g - mu_g * S_m_g
      dE_m_g <- lambda_hv_gra * S_m_g - eps_m_g * E_m_g - mu_g * E_m_g
      dI_m_g <- eps_m_g * E_m_g - mu_g * I_m_g

      # -------- 7. Vector SEI ODE at watering sites -----------------
      births_w <- rvf_adult_recruitment(M_w, K_w_t, r_m)
      dS_m_w <- births_w - lambda_hv_wat * S_m_w - mu_w * S_m_w
      dE_m_w <- lambda_hv_wat * S_m_w - eps_m_w * E_m_w - mu_w * E_m_w
      dI_m_w <- eps_m_w * E_m_w - mu_w * I_m_w
    }

    dy <- rep(0, length(y))
    dy[idx$S_h] <- dS_h; dy[idx$E_h] <- dE_h; dy[idx$I_h] <- dI_h; dy[idx$R_h] <- dR_h; dy[idx$C_h] <- dC_h; dy[idx$H_h] <- dH_h
    dy[idx$S_m_v] <- dS_m_v; dy[idx$E_m_v] <- dE_m_v; dy[idx$I_m_v] <- dI_m_v
    dy[idx$S_m_g] <- dS_m_g; dy[idx$E_m_g] <- dE_m_g; dy[idx$I_m_g] <- dI_m_g
    dy[idx$S_m_w] <- dS_m_w; dy[idx$E_m_w] <- dE_m_w; dy[idx$I_m_w] <- dI_m_w

    list(dy)
  })
}

rvf_compute_host_foi_state <- function(y, params, t) {
  idx <- params$idx
  month <- time_label_from_time(t, params$time_labels, params$days_per_step)

  G <- params$G_list[[month]]
  W <- params$W_list[[month]]
  # Temperature-dependent biting rates at this month for each site type,
  # matching the lookup inside the ODE right-hand side so that the FoI
  # decomposition reconstructs the same infection pressure the solver saw.
  a_v <- params$a_m_v[[month]]
  a_g <- params$a_m_g[[month]]
  a_w <- params$a_m_w[[month]]
  tb_month <- rvf_time_budget_at(
    list(t_v = params$t_v, t_g = params$t_g, t_w = params$t_w),
    month
  )

  # Mirror the ODE's non-negative rate calculations so diagnostics decompose
  # exactly the same force-of-infection logic used during integration.
  S_h <- rvf_nonnegative_state(y[idx$S_h])
  E_h <- rvf_nonnegative_state(y[idx$E_h])
  I_h <- rvf_nonnegative_state(y[idx$I_h])
  R_h <- rvf_nonnegative_state(y[idx$R_h])
  N_h <- pmax(S_h + E_h + I_h + R_h, 1)

  I_m_v <- rvf_nonnegative_state(y[idx$I_m_v])
  I_m_g <- rvf_nonnegative_state(y[idx$I_m_g])
  I_m_w <- rvf_nonnegative_state(y[idx$I_m_w])

  foi <- rvf_compute_foi_components(
    N_h = N_h,
    I_h = I_h,
    I_m_v = I_m_v,
    I_m_g = I_m_g,
    I_m_w = I_m_w,
    G = G,
    W = W,
    a_v = a_v,
    a_g = a_g,
    a_w = a_w,
    time_budget = tb_month,
    p_mobile = params$p_mobile,
    alpha_vh = params$alpha_vh,
    alpha_hv = params$alpha_hv,
    external_hazard = params$external_intro_hazard %||% rep(0, params$n_vil)
  )

  foi[c(
    "lambda_h", "lambda_h_village", "lambda_h_grazing",
    "lambda_h_watering", "lambda_h_external"
  )]
}

rvf_annual_metrics <- function(out, params, t_seed = 0) {
  state_cols <- function(idx) idx + 1L

  if (nrow(out) < 2) {
    stop("Need at least two output rows to compute annual metrics.")
  }

  C_h <- out[, state_cols(params$idx$C_h), drop = FALSE]
  H_h <- out[, state_cols(params$idx$H_h), drop = FALSE]
  I_h <- out[, state_cols(params$idx$I_h), drop = FALSE]
  host_pop <- as.numeric(params[["N_h_const"]])
  if (length(host_pop) == 1L) host_pop <- rep(host_pop, params$n_vil)
  host_pop <- pmax(host_pop, .Machine$double.eps)
  interval_rows <- vector("list", nrow(out) - 1L)
  system_rows <- vector("list", nrow(out) - 1L)

  # Total incidence / hazard comes from cumulative ODE states (exact over the
  # solver output interval). Source-specific FOI is reconstructed at interval
  # starts for decomposition, so use it as attribution rather than a second
  # independent fitted target.
  for (i in seq_len(nrow(out) - 1L)) {
    t0 <- out[i, 1]
    t1 <- out[i + 1L, 1]
    dt <- t1 - t0
    if (!is.finite(dt) || dt <= 0) {
      stop("Output times must be strictly increasing for annual metrics.")
    }
    analysis_year <- floor((t0 - t_seed) / 365) + 1L
    foi <- rvf_compute_host_foi_state(out[i, -1], params, t0)
    new_infections <- pmax(0, C_h[i + 1L, ] - C_h[i, ])
    cumulative_hazard <- pmax(0, H_h[i + 1L, ] - H_h[i, ])
    mean_daily_foi <- cumulative_hazard / dt
    infectious_host_days <- ((pmax(0, I_h[i, ]) + pmax(0, I_h[i + 1L, ])) / 2) * dt
    mean_I_h <- infectious_host_days / dt

    interval_rows[[i]] <- tibble(
      analysis_year = analysis_year,
      village_id = seq_len(params$n_vil),
      dt = dt,
      new_infections = as.numeric(new_infections),
      infectious_host_days = as.numeric(infectious_host_days),
      mean_infectious_hosts = as.numeric(mean_I_h),
      mean_prevalence = as.numeric(100 * mean_I_h / host_pop),
      cumulative_hazard = as.numeric(cumulative_hazard),
      lambda_h = as.numeric(mean_daily_foi),
      lambda_h_village = foi$lambda_h_village,
      lambda_h_grazing = foi$lambda_h_grazing,
      lambda_h_watering = foi$lambda_h_watering,
      lambda_h_external = foi$lambda_h_external
    )

    system_rows[[i]] <- tibble(
      analysis_year = analysis_year,
      dt = dt,
      annual_new_infections_interval = sum(new_infections),
      infectious_host_days_interval = sum(infectious_host_days),
      mean_infectious_hosts = sum(infectious_host_days) / dt,
      mean_prevalence = 100 * sum(infectious_host_days) / (sum(host_pop) * dt),
      mean_cumulative_hazard = mean(cumulative_hazard),
      mean_lambda_h = mean(mean_daily_foi),
      mean_lambda_h_village = mean(foi$lambda_h_village),
      mean_lambda_h_grazing = mean(foi$lambda_h_grazing),
      mean_lambda_h_watering = mean(foi$lambda_h_watering),
      mean_lambda_h_external = mean(foi$lambda_h_external)
    )
  }

  village_interval <- bind_rows(interval_rows)
  system_interval <- bind_rows(system_rows)

  village_year <- village_interval |>
    group_by(analysis_year, village_id) |>
    summarise(
      days = sum(dt),
      annual_new_infections = sum(new_infections),
      annual_infectious_host_days = sum(infectious_host_days),
      annual_mean_infectious_hosts = annual_infectious_host_days / days,
      annual_mean_prevalence = sum(mean_prevalence * dt) / days,
      annual_mean_daily_foi = sum(lambda_h * dt) / sum(dt),
      annual_mean_daily_foi_village = sum(lambda_h_village * dt) / sum(dt),
      annual_mean_daily_foi_grazing = sum(lambda_h_grazing * dt) / sum(dt),
      annual_mean_daily_foi_watering = sum(lambda_h_watering * dt) / sum(dt),
      annual_mean_daily_foi_external = sum(lambda_h_external * dt) / sum(dt),
      annual_cum_hazard = sum(cumulative_hazard),
      annual_foi_probability = 1 - exp(-annual_cum_hazard),
      annual_cum_hazard_village = sum(lambda_h_village * dt),
      annual_cum_hazard_grazing = sum(lambda_h_grazing * dt),
      annual_cum_hazard_watering = sum(lambda_h_watering * dt),
      annual_cum_hazard_external = sum(lambda_h_external * dt),
      .groups = "drop"
    )

  system_year <- system_interval |>
    group_by(analysis_year) |>
    summarise(
      days = sum(dt),
      annual_new_infections = sum(annual_new_infections_interval),
      annual_infectious_host_days = sum(infectious_host_days_interval),
      annual_mean_infectious_hosts = annual_infectious_host_days / days,
      annual_mean_prevalence = sum(mean_prevalence * dt) / days,
      mean_annual_cum_hazard = sum(mean_cumulative_hazard),
      annual_mean_daily_foi = mean_annual_cum_hazard / sum(dt),
      annual_mean_daily_foi_village = sum(mean_lambda_h_village * dt) / sum(dt),
      annual_mean_daily_foi_grazing = sum(mean_lambda_h_grazing * dt) / sum(dt),
      annual_mean_daily_foi_watering = sum(mean_lambda_h_watering * dt) / sum(dt),
      annual_mean_daily_foi_external = sum(mean_lambda_h_external * dt) / sum(dt),
      .groups = "drop"
    )

  system_foi_year <- village_year |>
    group_by(analysis_year) |>
    summarise(
      mean_annual_foi_probability = mean(annual_foi_probability),
      median_annual_foi_probability = median(annual_foi_probability),
      p90_annual_foi_probability = as.numeric(stats::quantile(annual_foi_probability, 0.9)),
      mean_annual_infectious_host_days = mean(annual_infectious_host_days),
      median_annual_infectious_host_days = median(annual_infectious_host_days),
      p90_annual_infectious_host_days = as.numeric(stats::quantile(annual_infectious_host_days, 0.9)),
      mean_annual_mean_prevalence = mean(annual_mean_prevalence),
      median_annual_mean_prevalence = median(annual_mean_prevalence),
      p90_annual_mean_prevalence = as.numeric(stats::quantile(annual_mean_prevalence, 0.9)),
      .groups = "drop"
    )

  system_year <- system_year |>
    left_join(system_foi_year, by = "analysis_year")

  distribution_year <- village_year |>
    group_by(analysis_year) |>
    summarise(
      mean_annual_mean_daily_foi = mean(annual_mean_daily_foi),
      median_annual_mean_daily_foi = median(annual_mean_daily_foi),
      p90_annual_mean_daily_foi = as.numeric(stats::quantile(annual_mean_daily_foi, 0.9)),
      mean_annual_cum_hazard = mean(annual_cum_hazard),
      median_annual_cum_hazard = median(annual_cum_hazard),
      p90_annual_cum_hazard = as.numeric(stats::quantile(annual_cum_hazard, 0.9)),
      mean_annual_foi_probability = mean(annual_foi_probability),
      median_annual_foi_probability = median(annual_foi_probability),
      p90_annual_foi_probability = as.numeric(stats::quantile(annual_foi_probability, 0.9)),
      .groups = "drop"
    )

  list(
    village_year = village_year,
    system_year = system_year,
    distribution_year = distribution_year
  )
}

rvf_get_time_series <- function(out, params, t_seed = 0) {
  state_cols <- function(idx) idx + 1L
  state_matrix <- function(idx) {
    mat <- as.matrix(out[, state_cols(idx), drop = FALSE])
    mat[!is.finite(mat)] <- 0
    pmax(mat, 0)
  }

  S_h <- state_matrix(params$idx$S_h)
  E_h <- state_matrix(params$idx$E_h)
  I_h <- state_matrix(params$idx$I_h)
  R_h <- state_matrix(params$idx$R_h)
  C_h <- state_matrix(params$idx$C_h)
  H_h <- state_matrix(params$idx$H_h)
  host_total <- pmax(rowSums(S_h + E_h + I_h + R_h), 1)

  data.frame(
    time = out[, 1],
    year_since_seed = (out[, 1] - t_seed) / 365,
    step_label = vapply(out[, 1], time_label_from_time, character(1), time_labels = params$time_labels, days_per_step = params$days_per_step),
    seroprevalence = rowSums(R_h) / host_total * 100,
    prevalence = rowSums(I_h) / host_total * 100,
    infectious_hosts = rowSums(I_h),
    cumulative_incidence = rowSums(C_h),
    mean_cumulative_hazard = rowMeans(H_h),
    mean_foi_probability_since_start = rowMeans(1 - exp(-H_h))
  )
}

rvf_output_times <- function(t_start, t_end, output_step, boundary_origin = t_start,
                             forcing_step = 365 / 12) {
  if (!is.finite(output_step) || output_step <= 0) {
    stop("output_step must be a positive number of days.")
  }
  if (!is.finite(t_start) || !is.finite(t_end) || t_end < t_start) {
    stop("Invalid output time range.")
  }

  base_times <- seq(t_start, t_end, by = output_step)
  # Include exact annual boundaries relative to the analysis origin so annual
  # hazard/incidence summaries are made of 365-day bins, not weekly fragments.
  first_boundary <- boundary_origin + ceiling((t_start - boundary_origin) / 365) * 365
  last_boundary <- boundary_origin + floor((t_end - boundary_origin) / 365) * 365
  year_boundaries <- if (first_boundary <= last_boundary) {
    seq(first_boundary, last_boundary, by = 365)
  } else {
    numeric(0)
  }

  # The RHS forcing is piecewise monthly. Include exact forcing boundaries
  # so LSODA lands on discontinuities instead of interpolating over them on
  # the requested output grid.
  first_forcing_boundary <- ceiling(t_start / forcing_step) * forcing_step
  last_forcing_boundary <- floor(t_end / forcing_step) * forcing_step
  forcing_boundaries <- if (first_forcing_boundary <= last_forcing_boundary) {
    seq(first_forcing_boundary, last_forcing_boundary, by = forcing_step)
  } else {
    numeric(0)
  }

  sort(unique(c(base_times, year_boundaries, forcing_boundaries, t_end)))
}

rvf_run_simulation <- function(params,
                               years_burn = 5,
                               years_endemic = 20,
                               output_step = 7,
                               n_seed = 10,
                               seed_village = 1,
                               init_state = NULL,
                               verbose = TRUE) {
  ode_controls <- rvf_ode_controls()
  if (is.null(init_state)) {
    y0 <- rvf_create_initial_state(params)
    times_burn <- rvf_output_times(0, max(years_burn, 1/365) * 365, output_step, boundary_origin = 0)

    if (years_burn > 0) {
      out_burn <- deSolve::ode(
        y = y0, times = times_burn, func = rvf_model_endemic,
        parms = params, method = "lsoda",
        rtol = ode_controls$rtol, atol = ode_controls$atol,
        maxsteps = ode_controls$maxsteps, verbose = FALSE
      )
      y_end <- out_burn[nrow(out_burn), -1]
      t_seed <- max(times_burn)
      out_burn_history <- out_burn
    } else {
      y_end <- y0
      t_seed <- 0
      out_burn_history <- NULL
    }

    y_ready <- if (n_seed > 0) rvf_seed_infection(y_end, params, n_seed, seed_village) else y_end
  } else {
    y_ready <- rvf_upgrade_state_vector(init_state, params)
    t_seed <- 0
    out_burn_history <- NULL
  }

  # C_h and H_h are diagnostics, not biological state. Reset them at the
  # start of every analysed endemic/future run so continuation states from
  # previous simulations do not pollute cumulative-incidence or hazard outputs.
  y_ready[params$idx$C_h] <- 0
  y_ready[params$idx$H_h] <- 0

  if (verbose) {
    q <- rvf_calculate_q_factors(params$p_mobile, list(t_v = params$t_v, t_g = params$t_g, t_w = params$t_w))
    cat(sprintf("Running fresh RVF simulation | p_mobile=%.2f | %s\n",
                params$p_mobile, rvf_format_q_factors(q)))
  }

  times_endemic <- rvf_output_times(t_seed, t_seed + years_endemic * 365, output_step, boundary_origin = t_seed)
  out_endemic <- deSolve::ode(
    y = y_ready, times = times_endemic, func = rvf_model_endemic,
    parms = params, method = "lsoda",
    rtol = ode_controls$rtol, atol = ode_controls$atol,
    maxsteps = ode_controls$maxsteps, verbose = FALSE
  )

  ts <- rvf_get_time_series(out_endemic, params, t_seed = t_seed)

  list(
    output_burn = out_burn_history,
    output_endemic = out_endemic,
    t_seed = t_seed,
    time_series = ts,
    params = params
  )
}
