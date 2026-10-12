#!/usr/bin/env Rscript
# Fit the three movement layers sequentially for grazing and watering networks.
# Each optimisation uses several dispersed starts and averages its loss over
# replicate simulations. The observed targets are one-mode village networks;
# the latent bipartite networks are projected before their fit is evaluated.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

source(here::here("analysis", "02_movement_networks", "_config.R"))
source(here::here("analysis", "02_movement_networks", "R", "utils.R"))
source(here::here("analysis", "02_movement_networks", "R", "simulator.R"))

cfg <- mvt_config
options(movement_n_cores = cfg$n_cores)

observed <- readRDS(cfg$observed_inputs_file)
targets_file <- file.path(cfg$calibration_dir, "targets.rds")
targets <- readRDS(targets_file)
regimes <- targets$regimes
dist_vv <- observed$distance_village_village

checkpoint_file <- file.path(cfg$calibration_dir, "fit_checkpoint.rds")
checkpoint_signature <- list(
  estimator_version = 4L,
  targets = unname(tools::md5sum(targets_file)),
  replicates = cfg$fit_replicates,
  starts = cfg$fit_starts,
  max_iterations = cfg$fit_max_iterations,
  use_restriction = cfg$use_restriction,
  monthly_route_turnover_fraction = cfg$monthly_route_turnover_fraction,
  drought_tail_turnover_fraction = cfg$drought_tail_turnover_fraction
)
checkpoint <- if (file.exists(checkpoint_file)) {
  candidate <- tryCatch(readRDS(checkpoint_file), error = function(error) NULL)
  if (!is.null(candidate) &&
      identical(candidate$signature, checkpoint_signature)) {
    candidate
  } else {
    list(signature = checkpoint_signature, fits = list())
  }
} else {
  list(signature = checkpoint_signature, fits = list())
}

# Save each fitted layer immediately so a long multi-resource fit can resume.
save_layer_checkpoint <- function(resource, layer, fitted_layer) {
  updated <- checkpoint
  resource_fits <- updated$fits[[resource]] %||% list()
  resource_fits[[layer]] <- fitted_layer
  updated$fits[[resource]] <- resource_fits
  checkpoint <<- updated
  saveRDS(checkpoint, checkpoint_file)
}

# Return a completed layer only when it matches the current targets and settings.
checkpointed_layer <- function(resource, layer) {
  checkpoint$fits[[resource]][[layer]] %||% NULL
}

# Select the lowest finite loss and retain every start for convergence auditing.
select_best_start <- function(results) {
  valid <- vapply(results, function(result) {
    is.finite(result$loss) && result$loss < non_finite_penalty
  }, logical(1))
  if (!any(valid)) {
    stop("All optimisation starts returned invalid losses.")
  }
  best_id <- which.min(vapply(results, `[[`, numeric(1), "loss"))
  best <- results[[best_id]]
  best$starts <- results
  best$best_start <- best_id
  best
}

# Require the configured number of starts explicitly. A mismatch should stop
# the fit rather than silently dropping a newly added exploratory start.
configured_starts <- function(starts, resource, layer) {
  expected_starts <- unname(cfg$fit_starts[[layer]])
  if (is.null(expected_starts) || length(expected_starts) != 1L) {
    stop("No configured start count for layer: ", layer)
  }
  if (length(starts) != expected_starts) {
    stop(
      resource, " ", layer, " defines ", length(starts),
      " starts, but the configured count is ", expected_starts, "."
    )
  }
  starts
}

non_finite_penalty <- 1e12
finite_objective <- function(loss_function) {
  function(parameters) {
    value <- tryCatch(
      loss_function(parameters),
      error = function(error) NA_real_
    )
    if (length(value) != 1L || !is.finite(value)) non_finite_penalty else value
  }
}

# Map unconstrained optimisation values to finite positive model parameters.
bounded_positive <- function(log_value, upper) {
  value <- exp(pmin(log_value, log(upper)))
  pmin(upper, pmax(1e-8, value))
}

# --- Loss and simulation summaries -----------------------------------------

scale_vec <- function(v) pmax(abs(v), 1, na.rm = TRUE)

# Normalised weighted squared error between two named numeric vectors.
wsse <- function(sim, obs, weights) {
  common <- Reduce(intersect, list(names(sim), names(obs), names(weights)))
  if (length(common) == 0L) {
    stop("The fitting loss found no shared named network metrics.")
  }
  sc <- scale_vec(obs[common])
  sum(weights[common] * ((sim[common] - obs[common]) / sc)^2, na.rm = TRUE)
}

# Assemble site distances, access restrictions, and monthly availability for one
# resource network using the same conventions as projection.
resource_model_inputs <- function(resource) {
  build_resource_inputs(
    observed = observed,
    resource = resource,
    availability_matrix = targets$observed_availability_matrices[[resource]],
    regimes = regimes,
    config = cfg
  )
}

# Summarise density and spatial reach at the observed one-mode network scale.
graph_fit_summary <- function(g) {
  d <- edge_distance_vector(g, dist_vv)
  c(
    edges = as.numeric(igraph::ecount(g)),
    mean_dist = if (length(d) == 0) NA_real_ else mean(d, na.rm = TRUE),
    p95_dist = if (length(d) == 0) NA_real_ else unname(stats::quantile(d, 0.95, na.rm = TRUE)),
    long_edges_15 = sum(d >= 15, na.rm = TRUE),
    long_edges_20 = sum(d >= 20, na.rm = TRUE),
    long_edges_25 = sum(d >= 25, na.rm = TRUE),
    long_edges_30 = sum(d >= 30, na.rm = TRUE)
  )
}

observed_year_summary <- function(resource) {
  targets[[resource]]$annual
}

# Average annual-union summaries over stochastic fitting replicates.
average_year_summary <- function(sims) {
  mat <- do.call(rbind, lapply(sims, function(sim) graph_fit_summary(sim$year_net)))
  colMeans(mat, na.rm = TRUE)
}

# --- Layer 1: year-round ---------------------------------------------------

pack1 <- function(th) c(log(th$lambda), log(th$gamma), th$beta_z)
unpack1 <- function(v) list(
  lambda = bounded_positive(v[1], 20),
  gamma = bounded_positive(v[2], 5),
  beta_z = v[3]
)

# Fit routine resource use with transition and drought additions switched off.
fit_year_round <- function(resource) {
  obs <- targets[[resource]]$summaries$routine   # 3-vec: edges, mean_dist, p95_dist
  # Routine long edges are sparse, so density and distance shape receive most
  # of the weight. Long-distance reach is fitted more strongly in drought.
  weights <- c(edges = 4.0, mean_dist = 1.0, p95_dist = 1.0,
               long_edges_15 = 0.75, long_edges_20 = 0.25)

  # Treat all months as routine so only the year-round layer contributes.
  regimes_all_routine <- setNames(rep("routine", 12), month.abb)
  theta_shim <- list(
    transition = list(lambda = 0, gamma = 1, beta_z = 0, p_keep = 1),
    drought  = list(lambda_local = 0, gamma = 1, beta_z = 0,
                    p_keep = 1, lambda_tail = 0, tau = 1)
  )

  loss_fn <- function(v) {
    th <- unpack1(v)
    inputs <- resource_model_inputs(resource)
    theta_layers <- list(
      year_round = th,
      transition   = theta_shim$transition,
      drought    = theta_shim$drought
    )
    sims <- simulate_replicates(theta_layers, inputs, regimes_all_routine,
                                n_reps = cfg$fit_replicates,
                                seed = cfg$seeds$fit_routine)
    summ <- average_layer_summaries(sims, dist_vv)$year_round
    wsse(summ, obs, weights)
  }

  starts <- list(
    list(lambda = 2.0, gamma = 0.7, beta_z = -1.5),
    list(lambda = 1.5, gamma = 1.0,  beta_z = -0.5),
    list(lambda = 2.5, gamma = 0.5,  beta_z = -2.0),
    list(lambda = 2.0, gamma = 0.4,  beta_z = -1.0),
    list(lambda = 1.0, gamma = 1.4,  beta_z = 0.0),
    list(lambda = 3.0, gamma = 0.25, beta_z = -2.5)
  )
  starts <- configured_starts(starts, resource, "year_round")
  objective <- finite_objective(loss_fn)

  results <- lapply(seq_along(starts), function(i) {
    v0 <- pack1(starts[[i]])
    l0 <- objective(v0)
    if (l0 >= non_finite_penalty) stop("Invalid year-round starting value ", i)
    res <- stats::optim(v0, objective, method = "Nelder-Mead",
                        control = list(maxit = cfg$fit_max_iterations,
                                       reltol = 1e-4))
    list(
      start_id = i,
      params = unpack1(res$par),
      loss = res$value,
      initial_loss = l0,
      convergence = res$convergence,
      evaluations = unname(res$counts[["function"]]),
      message = res$message %||% NA_character_
    )
  })
  select_best_start(results)
}

# --- Layer 2: transition (year-round frozen) --------------------------------

# Constrain the transition distance decay to be no stronger than the routine
# decay. This permits equal or longer typical journeys as resources decline.
pack2 <- function(th) {
  multiplier <- pmin(1 - 1e-6, pmax(1e-6, th$m_sh))
  c(log(th$lambda), qlogis(multiplier), th$beta_z)
}
unpack2 <- function(v) list(lambda = bounded_positive(v[1], 20),
                            m_sh = stats::plogis(v[2]),
                            beta_z = v[3])

# Fit only the additional transition routes, holding routine parameters fixed.
# `m_sh` is the internal multiplier that weakens distance decay in transition.
fit_transition <- function(resource, yr_params) {
  inputs <- resource_model_inputs(resource)
  obs <- targets[[resource]]$summaries$transition
  # Transition fit emphasises contact density and distance distribution.
  weights <- c(edges = 4.0, mean_dist = 1.0, p95_dist = 1.0,
               long_edges_15 = 0.75, long_edges_20 = 0.25)

  # Retain the observed regime labels so the transition objective is averaged
  # over the seven transition months, matching the observed target exactly.
  regimes_fit <- regimes
  theta_drought_shim <- list(lambda_local = 0, gamma = 1, beta_z = 0,
                             p_keep = 1, lambda_tail = 0, tau = 1)

  loss_fn <- function(v) {
    th <- unpack2(v)
    th$gamma  <- yr_params$gamma * th$m_sh
    theta_layers <- list(
      year_round = yr_params,
      transition   = th,
      drought    = theta_drought_shim
    )
    sims <- simulate_replicates(theta_layers, inputs, regimes_fit,
                                n_reps = cfg$fit_replicates,
                                seed = cfg$seeds$fit_transition)
    summ <- average_layer_summaries(sims, dist_vv)$transition
    wsse(summ, obs, weights)
  }

  starts <- switch(
    resource,
    grazing = list(
      list(lambda = 0.13, m_sh = 0.55, beta_z = -1.0),
      list(lambda = 0.6022, m_sh = pmin(0.95, 0.4909 / yr_params$gamma),
           beta_z = -1.749),
      list(lambda = 1.0, m_sh = 0.75, beta_z = -1.5),
      list(lambda = 0.7, m_sh = 0.55, beta_z = -2.0),
      list(lambda = 0.3, m_sh = 0.35, beta_z = -1.0),
      list(lambda = 1.4, m_sh = 0.90, beta_z = -0.5),
      list(lambda = 0.5, m_sh = 0.25, beta_z = -2.5)
    ),
    watering = list(
      list(lambda = 0.45, m_sh = 0.36, beta_z = -2.65),
      list(lambda = 0.9581, m_sh = pmin(0.95, 0.3564 / yr_params$gamma),
           beta_z = -2.131),
      list(lambda = 1.5, m_sh = 0.90, beta_z = -0.5),
      list(lambda = 0.7, m_sh = 0.55, beta_z = -2.0),
      list(lambda = 0.35, m_sh = 0.40, beta_z = -1.0),
      list(lambda = 1.1, m_sh = 0.70, beta_z = -1.5),
      list(lambda = 0.5, m_sh = 0.25, beta_z = -2.5)
    )
  )
  starts <- configured_starts(starts, resource, "transition")
  objective <- finite_objective(loss_fn)

  results <- lapply(seq_along(starts), function(i) {
    v0 <- pack2(starts[[i]])
    l0 <- objective(v0)
    if (l0 >= non_finite_penalty) stop("Invalid transition starting value ", i)
    res <- stats::optim(v0, objective, method = "Nelder-Mead",
                        control = list(maxit = cfg$fit_max_iterations,
                                       reltol = 1e-4))
    p <- unpack2(res$par)
    p$gamma <- yr_params$gamma * p$m_sh
    list(
      start_id = i,
      params = p,
      loss = res$value,
      initial_loss = l0,
      convergence = res$convergence,
      evaluations = unname(res$counts[["function"]]),
      message = res$message %||% NA_character_
    )
  })
  select_best_start(results)
}

# --- Layer 3: drought (year-round + transition frozen) ----------------------

# Local and heavy-tailed drought choices have separate Poisson intensities.
# This lets density and rare long-distance contacts be fitted independently.
pack3 <- function(th) {
  scaled_keep <- (th$p_keep - 0.3) / 0.7
  scaled_keep <- pmin(1 - 1e-6, pmax(1e-6, scaled_keep))
  c(log(th$lambda_local), log(th$gamma), th$beta_z,
    qlogis(scaled_keep), log(th$lambda_tail), log(th$tau))
}
unpack3 <- function(v) list(
                            lambda_local = bounded_positive(v[1], 20),
                            gamma = bounded_positive(v[2], 5),
                            beta_z = v[3],
                            p_keep = 0.3 + 0.7 * stats::plogis(v[4]),
                            lambda_tail = bounded_positive(v[5], 20),
                            tau = pmin(5, pmax(0.1, exp(v[6]))))

# Fit local and heavy-tailed drought additions while preserving the two fitted
# lower layers. Monthly targets govern typical drought structure; the lower-
# weighted annual union preserves information about rare long-distance edges.
fit_drought <- function(resource, yr_params, sh_params) {
  inputs <- resource_model_inputs(resource)
  obs <- targets[[resource]]$summaries$drought
  # Drought fit gives greater weight to spatial reach. Very rare 25 and 30 km
  # contacts enter through the annual union rather than unstable monthly means.
  weights <- c(edges = 0.5, mean_dist = 0.5, p95_dist = 4.0,
               long_edges_15 = 5.0, long_edges_20 = 5.0)
  obs_year <- observed_year_summary(resource)
  year_weights <- c(edges = cfg$fit_year_union_weight,
                    p95_dist = cfg$fit_year_union_weight,
                    long_edges_15 = cfg$fit_year_union_weight,
                    long_edges_20 = cfg$fit_year_union_weight,
                    long_edges_25 = cfg$fit_year_union_weight,
                    long_edges_30 = cfg$fit_year_union_weight)
  loss_fn <- function(v) {
    th <- unpack3(v)
    theta_layers <- list(
      year_round = yr_params,
      transition   = sh_params,
      drought    = th
    )
    sims <- simulate_replicates(theta_layers, inputs, regimes,
                                n_reps = cfg$fit_replicates,
                                seed = cfg$seeds$fit_drought)
    summ <- average_layer_summaries(sims, dist_vv)$drought
    wsse(summ, obs, weights) +
      wsse(average_year_summary(sims), obs_year, year_weights)
  }

  # Starts span local-heavy and tail-heavy drought behaviour.
  starts <- switch(
    resource,
    grazing = list(
      list(lambda_local = 0.19, gamma = 0.08, beta_z = -0.05,
           p_keep = 0.98, lambda_tail = 0.56, tau = 3.2),
      list(lambda_local = 0.4, gamma = 0.5, beta_z = -1.0,
           p_keep = 0.9, lambda_tail = 0.20, tau = 1.5),
      list(lambda_local = 0.3, gamma = 0.4, beta_z = -0.5,
           p_keep = 0.85, lambda_tail = 0.30, tau = 1.0),
      list(lambda_local = 0.5, gamma = 0.6, beta_z = -1.5,
           p_keep = 0.92, lambda_tail = 0.30, tau = 0.5),
      list(lambda_local = 0.2, gamma = 0.35, beta_z = -0.8,
           p_keep = 0.95, lambda_tail = 0.45, tau = 0.75),
      list(lambda_local = 0.6, gamma = 0.8, beta_z = -2.0,
           p_keep = 0.75, lambda_tail = 0.55, tau = 0.35),
      list(lambda_local = 0.15, gamma = 0.25, beta_z = 0,
           p_keep = 0.98, lambda_tail = 0.70, tau = 0.60)
    ),
    watering = list(
      list(lambda_local = 0.47, gamma = 0.16, beta_z = -1.0,
           p_keep = 0.98, lambda_tail = 0.47, tau = 2.2),
      list(lambda_local = 0.7, gamma = 0.4, beta_z = -0.7,
           p_keep = 0.9, lambda_tail = 0.25, tau = 1.5),
      list(lambda_local = 0.5, gamma = 0.3, beta_z = -0.5,
           p_keep = 0.85, lambda_tail = 0.40, tau = 1.0),
      list(lambda_local = 0.9, gamma = 0.5, beta_z = -1.0,
           p_keep = 0.92, lambda_tail = 0.30, tau = 0.5),
      list(lambda_local = 0.4, gamma = 0.25, beta_z = -0.8,
           p_keep = 0.95, lambda_tail = 0.50, tau = 0.75),
      list(lambda_local = 1.1, gamma = 0.7, beta_z = -1.8,
           p_keep = 0.75, lambda_tail = 0.60, tau = 0.35),
      list(lambda_local = 0.25, gamma = 0.18, beta_z = 0,
           p_keep = 0.98, lambda_tail = 0.75, tau = 0.60)
    )
  )
  starts <- configured_starts(starts, resource, "drought")
  objective <- finite_objective(loss_fn)

  results <- lapply(seq_along(starts), function(i) {
    v0 <- pack3(starts[[i]])
    l0 <- objective(v0)
    if (l0 >= non_finite_penalty) stop("Invalid drought starting value ", i)
    res <- stats::optim(v0, objective, method = "Nelder-Mead",
                        control = list(maxit = cfg$fit_max_iterations,
                                       reltol = 1e-4))
    list(
      start_id = i,
      params = unpack3(res$par),
      loss = res$value,
      initial_loss = l0,
      convergence = res$convergence,
      evaluations = unname(res$counts[["function"]]),
      message = res$message %||% NA_character_
    )
  })
  select_best_start(results)
}

# --- Run per resource ------------------------------------------------------

# Fit the three ecological layers in order for grazing or watering contacts.
calibrate_resource <- function(resource) {
  cat(sprintf("\n=========== %s ===========\n", resource))

  cat("--- Layer 1 (year-round) ---\n")
  fit1 <- checkpointed_layer(resource, "year_round")
  if (is.null(fit1)) {
    fit1 <- fit_year_round(resource)
    save_layer_checkpoint(resource, "year_round", fit1)
  } else {
    cat("  using completed checkpoint\n")
  }
  cat(sprintf("  loss = %.4f | lambda=%.3f gamma=%.3f beta_z=%.3f\n",
              fit1$loss, fit1$params$lambda, fit1$params$gamma,
              fit1$params$beta_z))

  cat("--- Layer 2 (transition residue, year-round frozen) ---\n")
  fit2 <- checkpointed_layer(resource, "transition")
  if (is.null(fit2)) {
    fit2 <- fit_transition(resource, fit1$params)
    save_layer_checkpoint(resource, "transition", fit2)
  } else {
    cat("  using completed checkpoint\n")
  }
  cat(sprintf(
    "  loss = %.4f | lambda=%.3f m_sh=%.3f gamma(eff)=%.3f beta_z=%.3f\n",
    fit2$loss, fit2$params$lambda, fit2$params$m_sh, fit2$params$gamma,
    fit2$params$beta_z))

  cat("--- Layer 3 (drought residue, year-round + transition frozen) ---\n")
  fit3 <- checkpointed_layer(resource, "drought")
  if (is.null(fit3)) {
    fit3 <- fit_drought(resource, fit1$params, fit2$params)
    save_layer_checkpoint(resource, "drought", fit3)
  } else {
    cat("  using completed checkpoint\n")
  }
  cat(sprintf(
    "  loss = %.4f | lambda_local=%.3f gamma=%.3f beta_z=%.3f lambda_tail=%.3f tau=%.3f p_keep=%.3f\n",
    fit3$loss, fit3$params$lambda_local, fit3$params$gamma,
    fit3$params$beta_z, fit3$params$lambda_tail, fit3$params$tau,
    fit3$params$p_keep))

  list(year_round = fit1, transition = fit2, drought = fit3)
}

fits <- list(
  grazing  = calibrate_resource("grazing"),
  watering = calibrate_resource("watering")
)

fit_file <- file.path(cfg$calibration_dir, "fit.rds")
saveRDS(
  list(
    fits = fits,
    targets_digest = unname(tools::md5sum(targets_file)),
    regimes = targets$regimes,
    thresholds = targets$thresholds,
    availability_parameters = targets$availability_parameters,
    observed_availability_matrices = targets$observed_availability_matrices,
    use_restriction = cfg$use_restriction,
    monthly_route_turnover_fraction = cfg$monthly_route_turnover_fraction,
    drought_tail_turnover_fraction = cfg$drought_tail_turnover_fraction,
    fit_settings = list(
      replicates = cfg$fit_replicates,
      starts = cfg$fit_starts,
      max_iterations = cfg$fit_max_iterations,
      seeds = cfg$seeds[c("fit_routine", "fit_transition", "fit_drought")]
    ),
    observation_scale = "one-mode village networks"
  ),
  fit_file
)

# A long table is easier to inspect and does not require empty columns for
# parameters that only occur in one layer.
estimate_rows <- bind_rows(lapply(names(fits), function(resource) {
  bind_rows(lapply(names(fits[[resource]]), function(layer) {
    selected_fit <- fits[[resource]][[layer]]
    enframe(unlist(selected_fit$params), name = "parameter", value = "estimate") |>
      mutate(
        resource = resource,
        layer = layer,
        selected_loss = selected_fit$loss,
        .before = 1
      )
  }))
}))

estimate_file <- file.path(cfg$calibration_dir, "point_estimates.csv")
write_csv(
  estimate_rows |> mutate(across(where(is.numeric), ~ signif(.x, 6))),
  estimate_file
)

start_rows <- dplyr::bind_rows(lapply(names(fits), function(resource) {
  dplyr::bind_rows(lapply(names(fits[[resource]]), function(layer) {
    fit_layer <- fits[[resource]][[layer]]
    dplyr::bind_rows(lapply(fit_layer$starts, function(start) {
      tibble::tibble(
        resource = resource,
        layer = layer,
        start_id = start$start_id,
        selected = start$start_id == fit_layer$best_start,
        initial_loss = start$initial_loss,
        final_loss = start$loss,
        convergence = start$convergence,
        evaluations = start$evaluations,
        message = start$message
      )
    }))
  }))
}))
diagnostics_file <- file.path(cfg$calibration_dir, "multistart_diagnostics.csv")
write_csv(start_rows, diagnostics_file)

if (file.exists(checkpoint_file)) unlink(checkpoint_file)

cat("\nFit complete. Outputs:\n")
cat("  ", fit_file, "\n")
cat("  ", estimate_file, "\n")
cat("  ", diagnostics_file, "\n")
