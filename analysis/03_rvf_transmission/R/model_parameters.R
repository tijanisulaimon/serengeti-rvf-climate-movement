#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

source(here::here("analysis", "03_rvf_transmission", "R", "scenario_io.R"))

rvf_get_rainfall_model_spec <- function(rainfall_model = NULL,
                                        cfg = load_rvf_config()) {
  if (is.null(rainfall_model)) {
    rainfall_model <- cfg$sim$rainfall_model %||% cfg$primary_rainfall_model
  }

  spec <- cfg$rainfall_models[[rainfall_model]]
  if (is.null(spec)) {
    stop("Unknown rainfall model: ", rainfall_model)
  }

  c(list(rainfall_model = rainfall_model), spec)
}

rvf_rainfall_model_label <- function(rainfall_model = NULL,
                                     cfg = load_rvf_config()) {
  spec <- rvf_get_rainfall_model_spec(rainfall_model = rainfall_model, cfg = cfg)
  spec$label
}

rvf_get_site_profile_spec <- function(site_profile = NULL,
                                      cfg = load_rvf_config()) {
  if (is.null(site_profile)) {
    site_profile <- cfg$sim$site_profile %||% "balanced"
  }
  site_profile <- rvf_resolve_site_profile_id(site_profile, cfg = cfg)

  spec <- cfg$site_profiles[[site_profile]]
  if (is.null(spec)) {
    stop("Unknown site-capacity profile: ", site_profile)
  }

  ratios <- spec$ratios
  required <- c("village", "grazing", "watering")
  if (is.null(ratios) || any(!required %in% names(ratios))) {
    stop("Site profile must define ratios named village, grazing, watering: ", site_profile)
  }
  ratios <- as.numeric(ratios[required])
  if (any(!is.finite(ratios)) || any(ratios <= 0)) {
    stop("Site profile ratios must be finite and positive: ", site_profile)
  }

  base_mean <- cfg$site_profile_base_mean_m_max %||% mean(c(2, 6, 30))
  m_max <- ratios / mean(ratios) * base_mean
  names(m_max) <- required

  list(
    site_profile = site_profile,
    site_profile_label = spec$label,
    site_profile_role = spec$role %||% NA_character_,
    site_profile_ratio_v = ratios[[1]],
    site_profile_ratio_g = ratios[[2]],
    site_profile_ratio_w = ratios[[3]],
    site_profile_base_mean_m_max = base_mean,
    profile_m_max_v = m_max[["village"]],
    profile_m_max_g = m_max[["grazing"]],
    profile_m_max_w = m_max[["watering"]]
  )
}

rvf_site_profile_label <- function(site_profile = NULL,
                                   cfg = load_rvf_config()) {
  spec <- rvf_get_site_profile_spec(site_profile = site_profile, cfg = cfg)
  spec$site_profile_label
}

rvf_rainfall_response_table <- function(rainfall_models = NULL,
                                        rainfall_mm = seq(0, 300, by = 5),
                                        cfg = load_rvf_config()) {
  if (is.null(rainfall_models)) {
    rainfall_models <- cfg$sim$rainfall_models_to_compare %||% cfg$sim$rainfall_model
  }

  bind_rows(lapply(rainfall_models, function(model_name) {
    spec <- rvf_get_rainfall_model_spec(rainfall_model = model_name, cfg = cfg)
    tibble(
      rainfall_model = model_name,
      rainfall_label = spec$label,
      rainfall_mm = rainfall_mm,
      rainfall_fraction = rainfall_response_fraction(
        R = rainfall_mm,
        rainfall_model = model_name,
        R_half = spec$R_half %||% 60,
        n = spec$n_hill %||% 1.25,
        min_K = spec$min_K,
        rainfall_lower = spec$rainfall_lower %||% 20,
        rainfall_upper = spec$rainfall_upper %||% 500
      )
    )
  }))
}

rvf_get_calibrated_params <- function(rainfall_model = NULL,
                                      site_profile = NULL,
                                      cfg = load_rvf_config()) {
  N_h_const <- as.numeric(cfg$sim$host_population_per_village)
  if (!is.finite(N_h_const) || N_h_const <= 0) {
    stop("host_population_per_village in R/config.R must be positive.")
  }
  eq_frac <- 0.50

  # Pooled-vector maximum vector:host ratios at saturated carrying capacity.
  # These are competence-weighted effective contact-capacity assumptions, not
  # measured mosquito densities. `site_profile` controls their relative
  # allocation across villages, grazing sites, and watering sites.
  site_spec <- rvf_get_site_profile_spec(site_profile = site_profile, cfg = cfg)
  m_max_v <- site_spec$profile_m_max_v
  m_max_g <- site_spec$profile_m_max_g
  m_max_w <- site_spec$profile_m_max_w
  rainfall_spec <- rvf_get_rainfall_model_spec(rainfall_model = rainfall_model, cfg = cfg)

  list(
    site_profile = site_spec$site_profile,
    site_profile_label = site_spec$site_profile_label,
    site_profile_role = site_spec$site_profile_role,
    site_profile_ratio_v = site_spec$site_profile_ratio_v,
    site_profile_ratio_g = site_spec$site_profile_ratio_g,
    site_profile_ratio_w = site_spec$site_profile_ratio_w,
    profile_m_max_v = site_spec$profile_m_max_v,
    profile_m_max_g = site_spec$profile_m_max_g,
    profile_m_max_w = site_spec$profile_m_max_w,
    rainfall_model = rainfall_spec$rainfall_model,
    rainfall_label = rainfall_spec$label,
    R_half = rainfall_spec$R_half %||% NA_real_,
    n_hill = rainfall_spec$n_hill %||% NA_real_,
    min_K = rainfall_spec$min_K,
    rainfall_lower = rainfall_spec$rainfall_lower %||% 20,
    rainfall_upper = rainfall_spec$rainfall_upper %||% 500,
    rainfall_reference = rainfall_spec$reference %||% NA_character_,
    K0_v = m_max_v * N_h_const / eq_frac,
    K0_g = m_max_g * N_h_const / eq_frac,
    K0_w = m_max_w * N_h_const / eq_frac,
    m_max_v = m_max_v,
    m_max_g = m_max_g,
    m_max_w = m_max_w,
    r_m = 0.05,
    # Transmission probabilities aligned with Cecilia et al. (2020) Table 1:
    # host-to-vector alpha_hv = 0.6, vector-to-host alpha_vh = 0.4.
    # (Previously 0.65 / 0.60, a Cavalerie 2015 / Gachohi 2016 blend;
    # updated 2026-04-19 to keep the vector block internally consistent
    # with the Cecilia-anchored pooled-vector parameterisation.)
    alpha_hv = 0.6,
    alpha_vh = 0.4,
    N_h_const = N_h_const
  )
}

rvf_adjust_calibrated_params <- function(params = rvf_get_calibrated_params(),
                                         transmission_scale = 1,
                                         capacity_scale = 1,
                                         village_capacity_scale = 1,
                                         grazing_capacity_scale = 1,
                                         watering_capacity_scale = 1,
                                         site_profile = NULL,
                                         rainfall_model = NULL,
                                         alpha_hv = NULL,
                                         alpha_vh = NULL,
                                         m_max_v = NULL,
                                         m_max_g = NULL,
                                         m_max_w = NULL,
                                         r_m = NULL,
                                         R_half = NULL,
                                         n_hill = NULL,
                                         min_K = NULL,
                                         rainfall_lower = NULL,
                                         rainfall_upper = NULL,
                                         cfg = load_rvf_config()) {
  eq_frac <- params$m_max_v * params$N_h_const / params$K0_v
  scales <- c(transmission_scale, capacity_scale, village_capacity_scale, grazing_capacity_scale, watering_capacity_scale)
  if (any(!is.finite(scales)) || any(scales < 0)) {
    stop("Transmission and capacity scale parameters must be finite and non-negative.")
  }

  if (!is.null(site_profile)) {
    site_spec <- rvf_get_site_profile_spec(site_profile = site_profile, cfg = cfg)
    params$site_profile <- site_spec$site_profile
    params$site_profile_label <- site_spec$site_profile_label
    params$site_profile_role <- site_spec$site_profile_role
    params$site_profile_ratio_v <- site_spec$site_profile_ratio_v
    params$site_profile_ratio_g <- site_spec$site_profile_ratio_g
    params$site_profile_ratio_w <- site_spec$site_profile_ratio_w
    params$profile_m_max_v <- site_spec$profile_m_max_v
    params$profile_m_max_g <- site_spec$profile_m_max_g
    params$profile_m_max_w <- site_spec$profile_m_max_w
    params$m_max_v <- site_spec$profile_m_max_v
    params$m_max_g <- site_spec$profile_m_max_g
    params$m_max_w <- site_spec$profile_m_max_w
  }

  if (!is.null(rainfall_model)) {
    rainfall_spec <- rvf_get_rainfall_model_spec(rainfall_model = rainfall_model, cfg = cfg)
    params$rainfall_model <- rainfall_spec$rainfall_model
    params$rainfall_label <- rainfall_spec$label
    params$R_half <- rainfall_spec$R_half %||% NA_real_
    params$n_hill <- rainfall_spec$n_hill %||% NA_real_
    params$min_K <- rainfall_spec$min_K
    params$rainfall_lower <- rainfall_spec$rainfall_lower %||% 20
    params$rainfall_upper <- rainfall_spec$rainfall_upper %||% 500
    params$rainfall_reference <- rainfall_spec$reference %||% NA_character_
  }

  params$alpha_hv <- if (is.null(alpha_hv)) params$alpha_hv * transmission_scale else alpha_hv
  params$alpha_vh <- if (is.null(alpha_vh)) params$alpha_vh * transmission_scale else alpha_vh

  params$m_max_v <- if (is.null(m_max_v)) params$m_max_v * capacity_scale * village_capacity_scale else m_max_v
  params$m_max_g <- if (is.null(m_max_g)) params$m_max_g * capacity_scale * grazing_capacity_scale else m_max_g
  params$m_max_w <- if (is.null(m_max_w)) params$m_max_w * capacity_scale * watering_capacity_scale else m_max_w

  params$K0_v <- params$m_max_v * params$N_h_const / eq_frac
  params$K0_g <- params$m_max_g * params$N_h_const / eq_frac
  params$K0_w <- params$m_max_w * params$N_h_const / eq_frac

  if (!is.null(r_m)) params$r_m <- r_m
  if (!is.null(R_half)) params$R_half <- R_half
  if (!is.null(n_hill)) params$n_hill <- n_hill
  if (!is.null(min_K)) params$min_K <- min_K
  if (!is.null(rainfall_lower)) params$rainfall_lower <- rainfall_lower
  if (!is.null(rainfall_upper)) params$rainfall_upper <- rainfall_upper
  if (any(!is.finite(c(params$alpha_hv, params$alpha_vh, params$m_max_v, params$m_max_g, params$m_max_w, params$r_m))) ||
      any(c(params$alpha_hv, params$alpha_vh, params$m_max_v, params$m_max_g, params$m_max_w, params$r_m) < 0)) {
    stop("Adjusted transmission, vector-scale, and vector-growth parameters must be finite and non-negative.")
  }

  params
}

rvf_get_fixed_params <- function(waning_years = Inf) {
  # `waning_years` = mean post-infection immunity (years). Baseline = Inf
  # (no waning / effectively durable immunity after infection). Finite
  # waning is retained as a sensitivity toggle if a good endemic fit cannot
  # be obtained under the life-long-immunity baseline.
  waning_years <- suppressWarnings(as.numeric(waning_years))
  omega <- if (is.finite(waning_years) && waning_years > 0) {
    1 / (waning_years * 365)
  } else 0

  list(
    gamma_h = 1 / 6,         # infectious period: 6 days (Cecilia et al. 2020)
    epsilon_h = 1 / 2,       # intrinsic incubation: 2 days (Gaff et al. 2007)
    omega = omega,           # waning rate per day; 0 = lifelong immunity
    mu_h = 1 / (8 * 365)     # cattle demographic turnover: 8-year lifespan
  )
}

rvf_normalize_time_regime <- function(x) {
  x <- as.character(x)
  x[!nzchar(x) | is.na(x)] <- "routine"
  x[x %in% c("wet", "non_drought", "nondrought", "year_round", "normal")] <- "routine"
  x[x == "shoulder"] <- "transition"  # compatibility with early Module 02 outputs
  x[!x %in% c("routine", "transition", "drought")] <- "routine"
  x
}

rvf_time_budget_profile <- function(no_movement = FALSE) {
  if (no_movement) {
    return(list(
      routine = c(t_g = 0, t_w = 0, t_v = 1),
      transition = c(t_g = 0, t_w = 0, t_v = 1),
      drought = c(t_g = 0, t_w = 0, t_v = 1)
    ))
  }

  canonical <- c(t_g = 8.0 / 24, t_w = 1.0 / 24, t_v = 15.0 / 24)
  list(
    # Daily exposure budget used for the site-profile definitions:
    # 15 h at/near the village, 8 h grazing, 1 h watering. Regime labels are
    # retained for compatibility with scenario metadata; the current analysis
    # does not impose additional seasonal variation in residence time.
    routine = canonical,
    transition = canonical,
    drought = canonical
  )
}

rvf_get_time_budget <- function(no_movement = FALSE,
                                regimes = NULL,
                                time_labels = NULL) {
  profile <- rvf_time_budget_profile(no_movement = no_movement)

  if (is.null(time_labels)) {
    if (is.null(regimes)) {
      return(as.list(profile$routine))
    }
    time_labels <- names(regimes) %||% names(profile)
  }

  if (is.null(regimes)) {
    regimes <- rep("routine", length(time_labels))
  }
  if (length(regimes) == 1L && length(time_labels) > 1L) {
    regimes <- rep(regimes, length(time_labels))
  }
  if (length(regimes) != length(time_labels)) {
    stop("Time-budget regimes must have length 1 or match time_labels.")
  }

  regimes <- rvf_normalize_time_regime(regimes)
  names(regimes) <- time_labels

  t_mat <- do.call(rbind, lapply(regimes, function(regime) profile[[regime]]))
  rownames(t_mat) <- time_labels

  list(
    t_g = setNames(as.numeric(t_mat[, "t_g"]), time_labels),
    t_w = setNames(as.numeric(t_mat[, "t_w"]), time_labels),
    t_v = setNames(as.numeric(t_mat[, "t_v"]), time_labels),
    regime = regimes,
    profile = profile
  )
}

rvf_time_budget_at <- function(tb, month) {
  pick <- function(x) {
    if (length(x) == 1L) return(as.numeric(x)[1])
    if (!is.null(names(x)) && month %in% names(x)) return(as.numeric(x[[month]]))
    as.numeric(x[[1]])
  }

  list(
    t_v = pick(tb$t_v),
    t_g = pick(tb$t_g),
    t_w = pick(tb$t_w)
  )
}

rvf_calculate_q_factors <- function(p_mobile, tb = rvf_get_time_budget()) {
  p_mobile <- suppressWarnings(as.numeric(p_mobile)[1])
  if (!is.finite(p_mobile)) p_mobile <- 0
  p_mobile <- max(0, min(1, p_mobile))
  list(
    q_v = (1 - p_mobile) * 1.0 + p_mobile * tb$t_v,
    q_g = p_mobile * tb$t_g,
    q_w = p_mobile * tb$t_w
  )
}

rvf_effective_profile_shares <- function(site_profile = NULL,
                                         p_mobile = NULL,
                                         cfg = load_rvf_config(),
                                         regimes = c("routine", "transition", "drought")) {
  if (is.null(p_mobile)) {
    p_mobile <- cfg$sim$p_mobile %||% 1
  }
  spec <- rvf_get_site_profile_spec(site_profile = site_profile, cfg = cfg)
  ratios <- c(
    village = spec$site_profile_ratio_v,
    grazing = spec$site_profile_ratio_g,
    watering = spec$site_profile_ratio_w
  )
  regimes <- rvf_normalize_time_regime(regimes)
  profile <- rvf_time_budget_profile(no_movement = FALSE)

  bind_rows(lapply(regimes, function(regime) {
    tb <- as.list(profile[[regime]])
    q <- c(
      village = (1 - p_mobile) + p_mobile * tb$t_v,
      grazing = p_mobile * tb$t_g,
      watering = p_mobile * tb$t_w
    )
    effective_pressure <- ratios * q
    shares <- effective_pressure / sum(effective_pressure)

    tibble(
      site_profile = spec$site_profile,
      site_profile_label = spec$site_profile_label,
      regime = regime,
      p_mobile = p_mobile,
      q_village = unname(q[["village"]]),
      q_grazing = unname(q[["grazing"]]),
      q_watering = unname(q[["watering"]]),
      effective_share_village = unname(shares[["village"]]),
      effective_share_grazing = unname(shares[["grazing"]]),
      effective_share_watering = unname(shares[["watering"]])
    )
  }))
}

rvf_effective_profile_share_wide <- function(site_profile = NULL,
                                             p_mobile = NULL,
                                             cfg = load_rvf_config()) {
  shares <- rvf_effective_profile_shares(
    site_profile = site_profile,
    p_mobile = p_mobile,
    cfg = cfg
  )
  out <- list()
  for (i in seq_len(nrow(shares))) {
    regime <- shares$regime[[i]]
    out[[paste0("effective_share_village_", regime)]] <- shares$effective_share_village[[i]]
    out[[paste0("effective_share_grazing_", regime)]] <- shares$effective_share_grazing[[i]]
    out[[paste0("effective_share_watering_", regime)]] <- shares$effective_share_watering[[i]]
  }
  out
}

rvf_format_q_factors <- function(q) {
  fmt <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) return("NA")
    if (length(unique(round(x, 6))) == 1L) return(sprintf("%.3f", x[1]))
    sprintf("%.3f-%.3f", min(x), max(x))
  }

  sprintf("q_v=%s q_g=%s q_w=%s", fmt(q$q_v), fmt(q$q_g), fmt(q$q_w))
}

# Temperature-dependent vector traits follow the Cecilia et al. table-1
# forms, but the pooled-vector model uses the T-only mortality row. The
# humidity-adjusted alternative became numerically unstable in Serengeti
# conditions and was effectively being clamped back to the same shape.

# Extrinsic incubation rate (per day). Linear in mean daily T (degC);
# Barker et al. (2013), used for BOTH species rows in Cecilia Table 1.
epsilon_m_temp <- function(T) {
  pmax(0.0071 * T - 0.1038, 1e-8)
}

eip_days <- function(T) {
  1 / epsilon_m_temp(T)
}

# Adult vector lifespan (days). This is the temperature-only mortality
# formulation in Cecilia et al. (2020), Table 1, vector group i = 1
# (Ae. vexans arabiensis/Cx. tritaeniorhynchus parameter set). It is used as
# a pooled-vector proxy, not as a Cx. pipiens-specific curve.
lifespan_days <- function(T) {
  quad <- pmax(0.000148 * T^2 - 0.00667 * T + 0.1, 1e-8)
  1 / quad
}

mu_m_temp <- function(T) {
  1 / lifespan_days(T)
}

# Biting rate (bites per vector per day). Cecilia Table 1 uses
# a(T) = (1 + c) / g(T), where c = 0.17 is the proportion of double
# blood meals (Ba et al. 2006) and g(T) is the gonotrophic-cycle length
# in days, 1 / (0.0173 * (T - 9.6)) (Madder et al. 1983). Biting stops
# below 9.6 degC; clamp at zero keeps the rate physical.
biting_rate_temp <- function(T, c = 0.17) {
  (1 + c) * 0.0173 * pmax(T - 9.6, 0)
}

# Monthly-rainfall responses for the pooled-vector carrying-capacity block.
# Both are habitat-capacity proxies rather than explicit breeding-site
# hydrology. Temperature is deliberately excluded here because it already
# controls biting, mortality, and extrinsic incubation in the transmission
# equations.
quadratic_rainfall_response <- function(R, lower = 20, upper = 500,
                                        min_K = 0.05) {
  R <- pmax(0, as.numeric(R))
  R[!is.finite(R)] <- 0

  # Rescale the bounded quadratic to 0-1 before applying the dry-season
  # carrying-capacity fraction. Rainfall outside [lower, upper] contributes 0.
  response <- rep(0, length(R))
  within_range <- R > lower & R < upper
  response[within_range] <-
    (R[within_range] - lower) * (upper - R[within_range]) /
    (((upper - lower)^2) / 4)
  min_K + (1 - min_K) * response
}

saturating_rainfall_response <- function(R, R_half = 60, n = 1.25,
                                         min_K = 0.05) {
  R <- pmax(0, as.numeric(R))
  R[!is.finite(R)] <- 0
  response <- (R^n) / (R^n + R_half^n)
  min_K + (1 - min_K) * response
}

rainfall_response_fraction <- function(R,
                                       rainfall_model = "quadratic_rainfall",
                                       R_half = 60,
                                       n = 1.25,
                                       min_K = 0.05,
                                       rainfall_lower = 20,
                                       rainfall_upper = 500) {
  if (identical(rainfall_model, "quadratic_rainfall")) {
    return(quadratic_rainfall_response(
      R = R,
      lower = rainfall_lower,
      upper = rainfall_upper,
      min_K = min_K
    ))
  }

  if (identical(rainfall_model, "saturating_rainfall")) {
    return(saturating_rainfall_response(
      R = R,
      R_half = R_half,
      n = n,
      min_K = min_K
    ))
  }

  stop("Unknown rainfall model: ", rainfall_model)
}

carrying_capacity <- function(R,
                              K0,
                              rainfall_model = "quadratic_rainfall",
                              R_half = 60,
                              n = 1.25,
                              min_K = 0.05,
                              rainfall_lower = 20,
                              rainfall_upper = 500) {
  K0 * rainfall_response_fraction(
    R = R,
    rainfall_model = rainfall_model,
    R_half = R_half,
    n = n,
    min_K = min_K,
    rainfall_lower = rainfall_lower,
    rainfall_upper = rainfall_upper
  )
}

normalize_rows <- function(M, tol = 1e-12) {
  M <- as.matrix(M)
  # Network rows are probabilities. Treat missing / non-finite / negative
  # entries as unavailable routes before normalising, so a single bad cell
  # cannot propagate NA through the ODE.
  M[!is.finite(M)] <- 0
  M[M < 0] <- 0
  rs <- rowSums(M, na.rm = TRUE)
  out <- M
  ok <- rs > tol
  out[ok, ] <- out[ok, , drop = FALSE] / rs[ok]
  out[!ok, ] <- 0
  out
}

time_label_from_time <- function(t, time_labels, days_per_step = 365 / 12) {
  idx <- (floor(t / days_per_step) %% length(time_labels)) + 1
  time_labels[idx]
}

rvf_get_state_indices <- function(n_vil, n_gra, n_wat) {
  states_per_vil <- 9

  S_h <- seq(1, by = states_per_vil, length.out = n_vil)
  E_h <- seq(2, by = states_per_vil, length.out = n_vil)
  I_h <- seq(3, by = states_per_vil, length.out = n_vil)
  R_h <- seq(4, by = states_per_vil, length.out = n_vil)
  C_h <- seq(5, by = states_per_vil, length.out = n_vil)
  H_h <- seq(6, by = states_per_vil, length.out = n_vil)

  S_m_v <- seq(7, by = states_per_vil, length.out = n_vil)
  E_m_v <- seq(8, by = states_per_vil, length.out = n_vil)
  I_m_v <- seq(9, by = states_per_vil, length.out = n_vil)

  offset_gra <- n_vil * states_per_vil
  S_m_g <- offset_gra + seq(1, by = 3, length.out = n_gra)
  E_m_g <- offset_gra + seq(2, by = 3, length.out = n_gra)
  I_m_g <- offset_gra + seq(3, by = 3, length.out = n_gra)

  offset_wat <- offset_gra + n_gra * 3
  S_m_w <- offset_wat + seq(1, by = 3, length.out = n_wat)
  E_m_w <- offset_wat + seq(2, by = 3, length.out = n_wat)
  I_m_w <- offset_wat + seq(3, by = 3, length.out = n_wat)

  list(
    S_h = S_h, E_h = E_h, I_h = I_h, R_h = R_h, C_h = C_h, H_h = H_h,
    S_m_v = S_m_v, E_m_v = E_m_v, I_m_v = I_m_v,
    S_m_g = S_m_g, E_m_g = E_m_g, I_m_g = I_m_g,
    S_m_w = S_m_w, E_m_w = E_m_w, I_m_w = I_m_w,
    total_states = offset_wat + n_wat * 3
  )
}

rvf_get_legacy_state_indices <- function(n_vil, n_gra, n_wat) {
  states_per_vil <- 8

  S_h <- seq(1, by = states_per_vil, length.out = n_vil)
  E_h <- seq(2, by = states_per_vil, length.out = n_vil)
  I_h <- seq(3, by = states_per_vil, length.out = n_vil)
  R_h <- seq(4, by = states_per_vil, length.out = n_vil)
  C_h <- seq(5, by = states_per_vil, length.out = n_vil)

  S_m_v <- seq(6, by = states_per_vil, length.out = n_vil)
  E_m_v <- seq(7, by = states_per_vil, length.out = n_vil)
  I_m_v <- seq(8, by = states_per_vil, length.out = n_vil)

  offset_gra <- n_vil * states_per_vil
  S_m_g <- offset_gra + seq(1, by = 3, length.out = n_gra)
  E_m_g <- offset_gra + seq(2, by = 3, length.out = n_gra)
  I_m_g <- offset_gra + seq(3, by = 3, length.out = n_gra)

  offset_wat <- offset_gra + n_gra * 3
  S_m_w <- offset_wat + seq(1, by = 3, length.out = n_wat)
  E_m_w <- offset_wat + seq(2, by = 3, length.out = n_wat)
  I_m_w <- offset_wat + seq(3, by = 3, length.out = n_wat)

  list(
    S_h = S_h, E_h = E_h, I_h = I_h, R_h = R_h, C_h = C_h,
    S_m_v = S_m_v, E_m_v = E_m_v, I_m_v = I_m_v,
    S_m_g = S_m_g, E_m_g = E_m_g, I_m_g = I_m_g,
    S_m_w = S_m_w, E_m_w = E_m_w, I_m_w = I_m_w,
    total_states = offset_wat + n_wat * 3
  )
}

rvf_upgrade_state_vector <- function(y, params) {
  if (length(y) == params$idx$total_states) {
    return(y)
  }

  legacy_idx <- rvf_get_legacy_state_indices(params$n_vil, params$n_gra, params$n_wat)
  if (length(y) != legacy_idx$total_states) {
    stop("State vector length does not match current or legacy layout: ", length(y))
  }

  out <- rep(0, params$idx$total_states)
  fields <- c(
    "S_h", "E_h", "I_h", "R_h", "C_h",
    "S_m_v", "E_m_v", "I_m_v",
    "S_m_g", "E_m_g", "I_m_g",
    "S_m_w", "E_m_w", "I_m_w"
  )

  for (field in fields) {
    out[params$idx[[field]]] <- y[legacy_idx[[field]]]
  }

  out
}

rvf_validate_scenario_inputs <- function(scenario) {
  validate_rvf_scenario(scenario)

  time_labels <- names(scenario$G_list) %||% seq_along(scenario$G_list)

  for (m in time_labels) {
    stopifnot(nrow(scenario$G_list[[m]]) == scenario$n_vil)
    stopifnot(ncol(scenario$G_list[[m]]) == scenario$n_gra)
    stopifnot(nrow(scenario$W_list[[m]]) == scenario$n_vil)
    stopifnot(ncol(scenario$W_list[[m]]) == scenario$n_wat)

    stopifnot(nrow(scenario$climate$vil[[m]]) == scenario$n_vil)
    stopifnot(nrow(scenario$climate$gra[[m]]) == scenario$n_gra)
    stopifnot(nrow(scenario$climate$wat[[m]]) == scenario$n_wat)
  }

  invisible(TRUE)
}

rvf_precompute_temp_params <- function(climate, n_vil, n_gra, n_wat, use_temp_dependence = TRUE) {
  # Build per-month, per-site-type vectors for every vector trait that the
  # ODE looks up by month: EIP rate (eps), mortality (mu), and biting
  # rate (a). All three depend only on monthly mean temperature `tas`.
  epsilon_m_v <- epsilon_m_g <- epsilon_m_w <- list()
  mu_m_v <- mu_m_g <- mu_m_w <- list()
  a_m_v <- a_m_g <- a_m_w <- list()
  temp_summary <- list()
  time_labels <- names(climate$vil) %||% seq_along(climate$vil)

  # Fill NA in a vector with the vector's own mean so downstream
  # vectorised ops stay finite. If the whole vector is NA we return
  # zeros, which lets the simulation continue with a cold / non-biting
  # site rather than silently producing NaNs that spread via network ops.
  fill_na_with_mean <- function(x) {
    x[!is.finite(x)] <- NA_real_
    if (all(is.na(x))) return(rep(0, length(x)))
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  }

  for (m in time_labels) {
    T_vil <- fill_na_with_mean(climate$vil[[m]]$tas[seq_len(n_vil)])
    T_gra <- fill_na_with_mean(climate$gra[[m]]$tas[seq_len(n_gra)])
    T_wat <- fill_na_with_mean(climate$wat[[m]]$tas[seq_len(n_wat)])

    if (!use_temp_dependence) {
      T_vil[] <- 25; T_gra[] <- 25; T_wat[] <- 25
    }

    temp_summary[[m]] <- c(vil = mean(T_vil), gra = mean(T_gra), wat = mean(T_wat))

    epsilon_m_v[[m]] <- epsilon_m_temp(T_vil)
    epsilon_m_g[[m]] <- epsilon_m_temp(T_gra)
    epsilon_m_w[[m]] <- epsilon_m_temp(T_wat)

    mu_m_v[[m]] <- mu_m_temp(T_vil)
    mu_m_g[[m]] <- mu_m_temp(T_gra)
    mu_m_w[[m]] <- mu_m_temp(T_wat)

    a_m_v[[m]] <- biting_rate_temp(T_vil)
    a_m_g[[m]] <- biting_rate_temp(T_gra)
    a_m_w[[m]] <- biting_rate_temp(T_wat)
  }

  list(
    epsilon_m_v = epsilon_m_v,
    epsilon_m_g = epsilon_m_g,
    epsilon_m_w = epsilon_m_w,
    mu_m_v = mu_m_v,
    mu_m_g = mu_m_g,
    mu_m_w = mu_m_w,
    a_m_v = a_m_v,
    a_m_g = a_m_g,
    a_m_w = a_m_w,
    temp_summary = temp_summary
  )
}

# Function: assemble the complete parameter bundle used by the ODE solver.
# scenario: monthly climate and village-resource matrices for one run.
# calibrated: profile-specific vector-capacity and transmission settings.
# p_mobile: fraction of each herd following the daily resource commute.
# capacity_reference_p_mobile: fixed mobile fraction used only to scale
#   resource-site capacity; it should not vary across mobility counterfactuals.
rvf_assemble_params <- function(scenario,
                                calibrated = NULL,
                                no_movement = FALSE,
                                use_temp_dependence = TRUE,
                                p_mobile = 1.0,
                                capacity_reference_p_mobile = 1.0,
                                waning_years = Inf,     # baseline = no waning
                                external_introduction_per_year = 0,
                                verbose = FALSE) {
  rvf_validate_scenario_inputs(scenario)

  p_mobile <- suppressWarnings(as.numeric(p_mobile)[1])
  if (!is.finite(p_mobile)) {
    stop("p_mobile must be finite.")
  }
  p_mobile <- max(0, min(1, p_mobile))
  capacity_reference_p_mobile <- suppressWarnings(
    as.numeric(capacity_reference_p_mobile)[1]
  )
  if (!is.finite(capacity_reference_p_mobile)) {
    stop("capacity_reference_p_mobile must be finite.")
  }
  capacity_reference_p_mobile <- max(0, min(1, capacity_reference_p_mobile))

  if (is.null(calibrated)) calibrated <- rvf_get_calibrated_params()
  fixed <- rvf_get_fixed_params(waning_years = waning_years)
  tb <- rvf_get_time_budget(
    no_movement = no_movement,
    regimes = scenario$metadata$regimes %||% NULL,
    time_labels = names(scenario$G_list) %||% paste0("step_", seq_along(scenario$G_list))
  )
  external_introduction_per_year <- suppressWarnings(as.numeric(external_introduction_per_year)[1])
  if (!is.finite(external_introduction_per_year) || external_introduction_per_year < 0) {
    stop("external_introduction_per_year must be finite and non-negative.")
  }

  n_vil <- scenario$n_vil
  n_gra <- scenario$n_gra
  n_wat <- scenario$n_wat
  time_labels <- names(scenario$G_list) %||% paste0("step_", seq_along(scenario$G_list))

  G_list <- lapply(scenario$G_list, normalize_rows)
  W_list <- lapply(scenario$W_list, normalize_rows)

  # Resource-site K0 values are scaled to a reference mobile-host load. During
  # p_mobile counterfactuals this can be frozen at the calibrated baseline, so
  # reducing mobility changes host residence/aggregation without also changing
  # the vector habitat baseline.
  N_h_mobile <- capacity_reference_p_mobile * calibrated$N_h_const
  mean_N_h_gra <- rep(0, n_gra)
  mean_N_h_wat <- rep(0, n_wat)

  # Time fractions enter later in the FOI terms; keeping them out here avoids
  # applying exposure duration twice.
  for (m in time_labels) {
    mean_N_h_gra <- mean_N_h_gra + as.vector(t(G_list[[m]]) %*% rep(N_h_mobile, n_vil))
    mean_N_h_wat <- mean_N_h_wat + as.vector(t(W_list[[m]]) %*% rep(N_h_mobile, n_vil))
  }
  mean_N_h_gra <- mean_N_h_gra / length(time_labels)
  mean_N_h_wat <- mean_N_h_wat / length(time_labels)

  eq_frac_approx <- 0.50
  K0_v_vec <- rep(calibrated$K0_v, n_vil)
  K0_g_vec <- pmax(calibrated$m_max_g * mean_N_h_gra / eq_frac_approx, 100)
  K0_w_vec <- pmax(calibrated$m_max_w * mean_N_h_wat / eq_frac_approx, 100)

  K_v <- K_g <- K_w <- list()
  rainfall_vil <- rainfall_gra <- rainfall_wat <- list()

  for (m in time_labels) {
    R_vil <- scenario$climate$vil[[m]]$pr[seq_len(n_vil)]
    R_gra <- scenario$climate$gra[[m]]$pr[seq_len(n_gra)]
    R_wat <- scenario$climate$wat[[m]]$pr[seq_len(n_wat)]

    # A few boundary points have missing rainfall in some months.
    # Treat missing values as dry months so the carrying-capacity floor applies
    # rather than propagating NA through the ODE state.
    R_vil[!is.finite(R_vil)] <- 0
    R_gra[!is.finite(R_gra)] <- 0
    R_wat[!is.finite(R_wat)] <- 0

    rainfall_vil[[m]] <- R_vil
    rainfall_gra[[m]] <- R_gra
    rainfall_wat[[m]] <- R_wat

    K_v[[m]] <- carrying_capacity(
      R = R_vil,
      K0 = K0_v_vec,
      rainfall_model = calibrated$rainfall_model,
      R_half = calibrated$R_half,
      n = calibrated$n_hill,
      min_K = calibrated$min_K,
      rainfall_lower = calibrated$rainfall_lower,
      rainfall_upper = calibrated$rainfall_upper
    )
    K_g[[m]] <- carrying_capacity(
      R = R_gra,
      K0 = K0_g_vec,
      rainfall_model = calibrated$rainfall_model,
      R_half = calibrated$R_half,
      n = calibrated$n_hill,
      min_K = calibrated$min_K,
      rainfall_lower = calibrated$rainfall_lower,
      rainfall_upper = calibrated$rainfall_upper
    )
    K_w[[m]] <- carrying_capacity(
      R = R_wat,
      K0 = K0_w_vec,
      rainfall_model = calibrated$rainfall_model,
      R_half = calibrated$R_half,
      n = calibrated$n_hill,
      min_K = calibrated$min_K,
      rainfall_lower = calibrated$rainfall_lower,
      rainfall_upper = calibrated$rainfall_upper
    )
  }

  temp_params <- rvf_precompute_temp_params(
    climate = scenario$climate,
    n_vil = n_vil,
    n_gra = n_gra,
    n_wat = n_wat,
    use_temp_dependence = use_temp_dependence
  )

  idx <- rvf_get_state_indices(n_vil, n_gra, n_wat)

  out <- list(
    n_vil = n_vil,
    n_gra = n_gra,
    n_wat = n_wat,
    idx = idx,
    time_labels = time_labels,
    days_per_step = 365 / 12,
    G_list = G_list,
    W_list = W_list,
    K_v = K_v,
    K_g = K_g,
    K_w = K_w,
    rainfall_vil = rainfall_vil,
    rainfall_gra = rainfall_gra,
    rainfall_wat = rainfall_wat,
    epsilon_m_v = temp_params$epsilon_m_v,
    epsilon_m_g = temp_params$epsilon_m_g,
    epsilon_m_w = temp_params$epsilon_m_w,
    mu_m_v = temp_params$mu_m_v,
    mu_m_g = temp_params$mu_m_g,
    mu_m_w = temp_params$mu_m_w,
    # Per-site, per-month biting rates use the Cecilia et al. (2020) form.
    a_m_v = temp_params$a_m_v,
    a_m_g = temp_params$a_m_g,
    a_m_w = temp_params$a_m_w,
    temp_summary = temp_params$temp_summary,
    use_temp_dependence = use_temp_dependence,
    r_m = calibrated$r_m,
    alpha_hv = calibrated$alpha_hv,
    alpha_vh = calibrated$alpha_vh,
    rainfall_model = calibrated$rainfall_model,
    rainfall_label = calibrated$rainfall_label,
    gamma_h = fixed$gamma_h,
    epsilon_h = fixed$epsilon_h,
    omega = fixed$omega,
    mu_h = fixed$mu_h,
    N_h_const = calibrated$N_h_const,
    # Sensitivity-only external pressure: expected successful exposed-host
    # introductions per year across the whole village system. It is converted
    # to a uniform daily per-susceptible hazard so it preserves host mass and
    # can represent unmodelled trade/gifting/neighbouring-system importation.
    external_introduction_per_year = external_introduction_per_year,
    external_intro_hazard = rep(
      external_introduction_per_year / (365 * n_vil * calibrated$N_h_const),
      n_vil
    ),
    t_v = tb$t_v,
    t_g = tb$t_g,
    t_w = tb$t_w,
    time_budget_regime = tb$regime,
    time_budget_profile = tb$profile,
    p_mobile = p_mobile,
    capacity_reference_p_mobile = capacity_reference_p_mobile,
    metadata = scenario$metadata,
    # Site-profile identity carried forward so downstream summaries (07/08/09)
    # can group by profile. Without these, params$site_profile is NULL and
    # tibble() silently drops the column, breaking the aggregation group_by.
    site_profile = calibrated$site_profile %||% NA_character_,
    site_profile_label = calibrated$site_profile_label %||% NA_character_,
    site_profile_role = calibrated$site_profile_role %||% NA_character_,
    site_profile_ratio_v = calibrated$site_profile_ratio_v %||% NA_real_,
    site_profile_ratio_g = calibrated$site_profile_ratio_g %||% NA_real_,
    site_profile_ratio_w = calibrated$site_profile_ratio_w %||% NA_real_,
    # Profile-level m_max values (post-capacity-scale); the per-site K0_*_vec
    # already absorbs these, but exposing them on params allows summaries to
    # report the effective vector ceilings without re-running the helper.
    m_max_v = calibrated$m_max_v,
    m_max_g = calibrated$m_max_g,
    m_max_w = calibrated$m_max_w
  )

  if (verbose) {
    q <- rvf_calculate_q_factors(p_mobile, tb)
    cat(sprintf("Assembled fresh RVF params | n_vil=%d n_gra=%d n_wat=%d\n", n_vil, n_gra, n_wat))
    cat(sprintf("p_mobile=%.2f | %s\n", p_mobile, rvf_format_q_factors(q)))
    cat(sprintf("Temp mean range= %.1f to %.1f C\n",
                min(unlist(temp_params$temp_summary), na.rm = TRUE),
                max(unlist(temp_params$temp_summary), na.rm = TRUE)))
  }

  out
}

rvf_create_initial_state <- function(params) {
  y0 <- rep(0, params$idx$total_states)
  y0[params$idx$S_h] <- params$N_h_const

  mu_m_mean <- mean(unlist(params$mu_m_v), na.rm = TRUE)
  eq_frac <- if (params$r_m > mu_m_mean) (1 - mu_m_mean / params$r_m) else 0.5
  m0 <- params$time_labels[[1]]

  y0[params$idx$S_m_v] <- params$K_v[[m0]] * eq_frac
  p_mob <- suppressWarnings(as.numeric(params$p_mobile)[1])
  if (!is.finite(p_mob)) p_mob <- 0
  if (p_mob > 1e-10) {
    y0[params$idx$S_m_g] <- params$K_g[[m0]] * eq_frac
    y0[params$idx$S_m_w] <- params$K_w[[m0]] * eq_frac
  }

  y0
}

rvf_seed_infection <- function(y, params, n_seed = 10, village = 1) {
  n_seed <- suppressWarnings(as.numeric(n_seed)[1])
  if (!is.finite(n_seed) || n_seed < 0) {
    stop("n_seed must be a non-negative finite number.")
  }

  villages <- if (identical(village, "all")) seq_len(params$n_vil) else village
  if (any(!is.finite(villages)) || any(villages < 1) || any(villages > params$n_vil)) {
    stop("Seed village index out of range.")
  }

  for (v in villages) {
    to_seed <- min(n_seed, y[params$idx$S_h[v]])
    y[params$idx$S_h[v]] <- y[params$idx$S_h[v]] - to_seed
    y[params$idx$I_h[v]] <- y[params$idx$I_h[v]] + to_seed
  }
  y
}
