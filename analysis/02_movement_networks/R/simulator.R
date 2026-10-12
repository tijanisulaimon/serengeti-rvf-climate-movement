#!/usr/bin/env Rscript
# Three-layer village-resource movement-network simulator.
#
# A mostly persistent year-round route pool is supplemented by transition
# routes and, in drought months, local and heavy-tailed search routes. Monthly
# networks use subsets of these yearly pools according to resource severity.
# A small common fraction of routine and transition routes, and a larger fixed
# fraction of the drought tail, are redrawn between months.

suppressPackageStartupMessages(library(igraph))
source(here::here("analysis", "02_movement_networks", "R", "utils.R"))

# --- Destination-choice kernels --------------------------------------------

clean_availability <- function(availability, length_out) {
  if (is.null(availability)) return(rep(1, length_out))
  availability <- as.numeric(availability)
  availability[!is.finite(availability)] <- 0
  pmin(1, pmax(0, availability))
}

normalise_log_weights <- function(log_weights) {
  finite <- is.finite(log_weights)
  if (!any(finite)) return(rep(1 / length(log_weights), length(log_weights)))
  weights <- rep(0, length(log_weights))
  weights[finite] <- exp(log_weights[finite] - max(log_weights[finite]))
  if (sum(weights) <= 0) return(rep(1 / length(weights), length(weights)))
  weights / sum(weights)
}

# Routine and local choices decline exponentially with distance and favour
# resource sites in proportion to their current NDVI-derived availability.
exponential_choice_probabilities <- function(distance, is_restricted,
                                             gamma, beta_z,
                                             availability) {
  availability <- clean_availability(availability, length(distance))
  log_availability <- ifelse(availability > 0, log(availability), -Inf)
  normalise_log_weights(
    -gamma * distance + beta_z * is_restricted + log_availability
  )
}

# The drought-search kernel has a heavier distance tail and uses sqrt(A), so
# moderately suitable distant sites are not excluded as strongly as local sites.
tail_choice_probabilities <- function(distance, is_restricted,
                                      beta_z, tau,
                                      availability) {
  availability <- clean_availability(availability, length(distance))
  log_availability <- ifelse(availability > 0, 0.5 * log(availability), -Inf)
  normalise_log_weights(
    -tau * log1p(pmax(distance, 0)) +
      beta_z * is_restricted + log_availability
  )
}

draw_choices <- function(probability, number) {
  if (number <= 0L) return(integer(0))
  if (number >= length(probability)) return(seq_along(probability))
  # Extreme distance-decay proposals can underflow negligible weights to zero.
  # A tiny floor preserves their ordering while keeping weighted sampling valid.
  probability <- pmax(probability, sqrt(.Machine$double.xmin))
  probability <- probability / sum(probability)
  sample.int(length(probability), number, replace = FALSE, prob = probability)
}

empty_bipartite <- function(village_names, resource_names) {
  matrix(
    0L,
    nrow = length(village_names),
    ncol = length(resource_names),
    dimnames = list(village_names, resource_names)
  )
}

# Draw the persistent routine route pool, including at least one site per village.
draw_year_round <- function(parameters, distances, is_restricted,
                            village_names, availability) {
  number_villages <- nrow(distances)
  number_resources <- ncol(distances)
  choices_per_village <- pmin(
    pmax(1L, stats::rpois(number_villages, parameters$lambda)),
    number_resources
  )
  network <- empty_bipartite(village_names, colnames(distances))

  for (village in seq_len(number_villages)) {
    probability <- exponential_choice_probabilities(
      distances[village, ], is_restricted,
      parameters$gamma, parameters$beta_z, availability
    )
    selected <- draw_choices(probability, choices_per_village[[village]])
    network[village, selected] <- 1L
  }
  network
}

# Add resource-seeking routes without replacing the supplied routine network.
add_transition_choices <- function(parameters, network, distances,
                                   is_restricted, availability) {
  additions <- stats::rpois(nrow(distances), parameters$lambda)

  for (village in seq_len(nrow(distances))) {
    available_cells <- which(network[village, ] == 0L)
    number <- min(additions[[village]], length(available_cells))
    if (number == 0L) next

    probability <- exponential_choice_probabilities(
      distances[village, available_cells],
      is_restricted[available_cells],
      parameters$gamma,
      parameters$beta_z,
      availability[available_cells]
    )
    selected <- available_cells[draw_choices(probability, number)]
    network[village, selected] <- 1L
  }
  network
}

# Thin transition routes according to fitted persistence and site availability.
retain_transition_choices <- function(network, probability, availability) {
  probability <- pmin(1, pmax(0, probability))
  retention <- matrix(
    probability * availability,
    nrow(network), ncol(network), byrow = TRUE
  )
  keep <- matrix(stats::runif(length(network)), nrow(network), ncol(network)) <
    retention
  network[network == 1L & !keep] <- 0L
  network
}

# Replace a small fraction of selected routes while preserving each village's
# number of resource choices. Replacement sites follow the same fitted kernel
# and current-month availability as the corresponding ecological layer.
redraw_monthly_routes <- function(network, fraction, parameters, distances,
                                  is_restricted, availability) {
  fraction <- pmin(1, pmax(0, fraction))
  if (fraction <= 0) return(network)

  updated <- network
  for (village in seq_len(nrow(network))) {
    selected <- which(network[village, ] == 1L)
    if (length(selected) == 0L) next

    dropped <- selected[stats::runif(length(selected)) < fraction]
    if (length(dropped) == 0L) next

    # Exclude dropped sites from immediate reselection so that every accepted
    # redraw represents an actual change of resource route.
    alternatives <- setdiff(which(network[village, ] == 0L), dropped)
    number <- min(length(dropped), length(alternatives))
    if (number == 0L) next

    probability <- exponential_choice_probabilities(
      distances[village, alternatives],
      is_restricted[alternatives],
      parameters$gamma,
      parameters$beta_z,
      availability[alternatives]
    )
    replacement <- alternatives[draw_choices(probability, number)]
    changed <- dropped[seq_len(number)]
    updated[village, changed] <- 0L
    updated[village, replacement] <- 1L
  }
  updated
}

# Retain part of the transition pool, then add distance-limited drought routes.
add_drought_local_choices <- function(parameters, network, distances,
                                      is_restricted, availability) {
  network <- retain_transition_choices(
    network, parameters$p_keep, availability
  )
  additions <- stats::rpois(nrow(distances), parameters$lambda)

  for (village in seq_len(nrow(distances))) {
    available_cells <- which(network[village, ] == 0L)
    number <- min(additions[[village]], length(available_cells))
    if (number == 0L) next

    probability <- exponential_choice_probabilities(
      distances[village, available_cells],
      is_restricted[available_cells],
      parameters$gamma,
      parameters$beta_z,
      availability[available_cells]
    )
    selected <- available_cells[draw_choices(probability, number)]
    network[village, selected] <- 1L
  }
  network
}

# Add rare long-distance drought-search routes using the heavy-tailed kernel.
add_drought_tail_choices <- function(parameters, network, distances,
                                     is_restricted, availability) {
  additions <- stats::rpois(nrow(distances), parameters$lambda)

  for (village in seq_len(nrow(distances))) {
    available_cells <- which(network[village, ] == 0L)
    number <- min(additions[[village]], length(available_cells))
    if (number == 0L) next

    probability <- tail_choice_probabilities(
      distances[village, available_cells],
      is_restricted[available_cells],
      parameters$beta_z,
      parameters$tau,
      availability[available_cells]
    )
    selected <- available_cells[draw_choices(probability, number)]
    network[village, selected] <- 1L
  }
  network
}

union_bipartite <- function(network, persistent_base) {
  result <- (network | persistent_base) + 0L
  storage.mode(result) <- "integer"
  dimnames(result) <- dimnames(persistent_base)
  result
}

# --- Climate and resource inputs -------------------------------------------

# Convert prepared data into the common input object used for fitting,
# validation, and future projection.
build_resource_inputs <- function(observed, resource, availability_matrix,
                                  regimes, config) {
  distances <- observed[[resource]]$distance_village_resource
  availability <- resource_availability_inputs_from_matrix(
    availability_matrix,
    regimes = regimes,
    resource_names = colnames(distances)
  )

  list(
    village_names = observed$village_names,
    distances = distances,
    is_restricted = if (isTRUE(config$use_restriction)) {
      as.numeric(observed[[resource]]$is_restricted)
    } else {
      rep(0, ncol(distances))
    },
    availability_by_regime = availability$by_regime,
    availability_by_month = availability$by_month,
    severity_strength = config$availability_severity_strength,
    severity_cap = config$availability_severity_cap,
    severity_memory = config$availability_severity_memory,
    route_turnover = config$monthly_route_turnover_fraction,
    tail_turnover = config$drought_tail_turnover_fraction
  )
}

availability_for_regime <- function(inputs, regime) {
  availability <- inputs$availability_by_regime[[regime]]
  if (is.null(availability)) {
    availability <- rowMeans(
      do.call(cbind, inputs$availability_by_month), na.rm = TRUE
    )
  }
  clean_availability(availability, ncol(inputs$distances))
}

availability_for_month <- function(inputs, month, regime) {
  availability <- inputs$availability_by_month[[month]] %||%
    inputs$availability_by_regime[[regime]]
  clean_availability(availability, ncol(inputs$distances))
}

previous_month <- function(month) {
  index <- match(month, month.abb)
  month.abb[if (index == 1L) 12L else index - 1L]
}

# Scale extra route intensity when a month is drier than its regime average.
# This varies route use within a regime without fitting twelve monthly models.
severity_multiplier <- function(inputs, month, regime) {
  if (!regime %in% c("transition", "drought")) return(1)
  if (inputs$severity_strength <= 0) return(1)

  regime_availability <- availability_for_regime(inputs, regime)
  current <- availability_for_month(inputs, month, regime)
  effective <- if (identical(inputs$severity_memory, "lag_min")) {
    pmin(current, availability_for_month(inputs, previous_month(month), regime))
  } else {
    current
  }

  relative_shortfall <- pmax(
    0,
    mean(regime_availability) - mean(effective)
  ) / pmax(mean(regime_availability), 0.05)
  pmin(inputs$severity_cap, 1 + inputs$severity_strength * relative_shortfall)
}

maximum_severity <- function(inputs, months, regime) {
  if (length(months) == 0L) return(1)
  max(vapply(months, severity_multiplier, numeric(1),
             inputs = inputs, regime = regime))
}

# Activate a severity-dependent subset of non-routine routes while always
# preserving the year-round base.
activate_route_pool <- function(pool, persistent_base, activity) {
  activity <- pmin(1, pmax(0, activity))
  network <- persistent_base
  additional <- pool == 1L & persistent_base == 0L
  if (activity >= 0.999) {
    network[additional] <- 1L
  } else if (any(additional)) {
    active <- matrix(stats::runif(length(pool)), nrow(pool), ncol(pool)) < activity
    network[additional & active] <- 1L
  }
  storage.mode(network) <- "integer"
  dimnames(network) <- dimnames(persistent_base)
  network
}

# Form an annual village contact network from monthly co-use edges.
union_one_mode_graphs <- function(monthly_graphs, village_names) {
  keys <- unique(unlist(lapply(monthly_graphs, edge_keys), use.names = FALSE))
  if (length(keys) == 0L) return(make_graph(data.frame(), village_names))

  endpoints <- strsplit(keys, "--", fixed = TRUE)
  make_graph(
    data.frame(
      from = vapply(endpoints, `[`, character(1), 1),
      to = vapply(endpoints, `[`, character(1), 2)
    ),
    village_names
  )
}

# --- Annual simulation ------------------------------------------------------

# Draw persistent route pools once per simulated year, then activate monthly
# subsets according to regime and within-regime resource severity.
simulate_year <- function(parameters, inputs, regimes, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  village_names <- inputs$village_names
  distances <- inputs$distances
  is_restricted <- inputs$is_restricted

  year_round <- draw_year_round(
    parameters$year_round,
    distances,
    is_restricted,
    village_names,
    availability_for_regime(inputs, "routine")
  )

  transition_months <- month.abb[
    regimes[month.abb] %in% c("transition", "drought")
  ]
  drought_months <- month.abb[regimes[month.abb] == "drought"]
  maximum_transition <- maximum_severity(
    inputs, transition_months, "transition"
  )
  maximum_drought <- maximum_severity(inputs, drought_months, "drought")

  transition_parameters <- parameters$transition
  transition_parameters$lambda <-
    transition_parameters$lambda * maximum_transition
  transition_pool <- if (length(transition_months) > 0L) {
    union_bipartite(
      add_transition_choices(
        transition_parameters,
        year_round,
        distances,
        is_restricted,
        availability_for_regime(inputs, "transition")
      ),
      year_round
    )
  } else {
    year_round
  }

  drought_parameters <- parameters$drought
  local_parameters <- list(
    lambda = drought_parameters$lambda_local * maximum_drought,
    gamma = drought_parameters$gamma,
    beta_z = drought_parameters$beta_z,
    p_keep = pmin(1, pmax(0.3, drought_parameters$p_keep))
  )
  drought_availability <- availability_for_regime(inputs, "drought")
  drought_local_pool <- if (length(drought_months) > 0L) {
    union_bipartite(
      add_drought_local_choices(
        local_parameters,
        transition_pool,
        distances,
        is_restricted,
        drought_availability
      ),
      year_round
    )
  } else {
    transition_pool
  }

  tail_turnover <- pmin(1, pmax(0, inputs$tail_turnover))
  stable_tail_parameters <- list(
    lambda = drought_parameters$lambda_tail * maximum_drought *
      (1 - tail_turnover),
    beta_z = drought_parameters$beta_z,
    tau = pmin(5, pmax(0.1, drought_parameters$tau))
  )
  drought_stable_pool <- if (length(drought_months) > 0L) {
    union_bipartite(
      add_drought_tail_choices(
        stable_tail_parameters,
        drought_local_pool,
        distances,
        is_restricted,
        drought_availability
      ),
      year_round
    )
  } else {
    drought_local_pool
  }

  bipartite <- list()
  one_mode <- list()
  for (month in month.abb) {
    regime <- unname(regimes[[month]])
    network <- if (regime == "routine") {
      year_round
    } else if (regime == "transition") {
      activity <- severity_multiplier(inputs, month, "transition") /
        maximum_transition
      activate_route_pool(transition_pool, year_round, activity)
    } else {
      month_severity <- severity_multiplier(inputs, month, "drought")
      stable_network <- activate_route_pool(
        drought_stable_pool, year_round, month_severity / maximum_drought
      )
      changing_tail <- list(
        lambda = drought_parameters$lambda_tail * month_severity * tail_turnover,
        beta_z = drought_parameters$beta_z,
        tau = pmin(5, pmax(0.1, drought_parameters$tau))
      )
      union_bipartite(
        add_drought_tail_choices(
          changing_tail,
          stable_network,
          distances,
          is_restricted,
          availability_for_month(inputs, month, "drought")
        ),
        year_round
      )
    }

    if (regime %in% c("routine", "transition") &&
        inputs$route_turnover > 0) {
      kernel_parameters <- if (regime == "routine") {
        parameters$year_round
      } else {
        parameters$transition
      }
      network <- redraw_monthly_routes(
        network,
        fraction = inputs$route_turnover,
        parameters = kernel_parameters,
        distances = distances,
        is_restricted = is_restricted,
        availability = availability_for_month(inputs, month, regime)
      )
    }

    bipartite[[month]] <- network
    one_mode[[month]] <- project_one_mode(network, village_names)
  }

  list(
    bipartite = bipartite,
    one_mode = one_mode,
    year_net = union_one_mode_graphs(one_mode, village_names),
    regimes = regimes
  )
}

# Generate independent stochastic network years, optionally in parallel.
simulate_replicates <- function(parameters, inputs, regimes,
                                n_reps, seed = NULL,
                                n_cores = getOption("movement_n_cores", 1L)) {
  seeds <- if (is.null(seed)) {
    rep(list(NULL), n_reps)
  } else {
    as.list(as.integer(seed + seq_len(n_reps)))
  }
  run_one <- function(current_seed) {
    simulate_year(parameters, inputs, regimes, seed = current_seed)
  }

  results <- if (.Platform$OS.type == "unix" && n_cores > 1L) {
    parallel::mclapply(
      seeds, run_one,
      mc.cores = min(n_cores, n_reps), mc.preschedule = FALSE
    )
  } else {
    lapply(seeds, run_one)
  }

  # A try-error marks a failed replicate; NULL marks a worker process that
  # stopped without returning, which mclapply reports only as a warning.
  failed <- vapply(results, function(result) {
    is.null(result) || inherits(result, "try-error")
  }, logical(1))
  if (any(failed)) {
    first <- results[[which(failed)[1]]]
    stop("One or more replicate simulations failed: ",
         if (is.null(first)) "a worker process stopped without returning a result"
         else as.character(first))
  }
  results
}

# The fitting loss is evaluated at the observed one-mode scale.
simulation_layer_summaries <- function(simulation, village_distances) {
  summarise_graph <- function(graph) {
    distance <- edge_distance_vector(graph, village_distances)
    c(
      edges = igraph::ecount(graph),
      mean_dist = if (length(distance) == 0L) NA_real_ else mean(distance),
      p95_dist = if (length(distance) == 0L) NA_real_ else
        unname(stats::quantile(distance, 0.95)),
      long_edges_15 = sum(distance >= 15),
      long_edges_20 = sum(distance >= 20),
      long_edges_25 = sum(distance >= 25),
      long_edges_30 = sum(distance >= 30)
    )
  }

  summary_template <- stats::setNames(
    numeric(7),
    c(
      "edges", "mean_dist", "p95_dist", "long_edges_15",
      "long_edges_20", "long_edges_25", "long_edges_30"
    )
  )

  setNames(lapply(c("routine", "transition", "drought"), function(regime) {
    months <- names(simulation$regimes)[simulation$regimes == regime]
    if (length(months) == 0L) return(summary_template + NA_real_)
    summaries <- vapply(
      months,
      function(month) summarise_graph(simulation$one_mode[[month]]),
      summary_template
    )
    rowMeans(summaries, na.rm = TRUE)
  }), c("year_round", "transition", "drought"))
}

# Average observed-scale layer summaries over stochastic replicates for fitting.
average_layer_summaries <- function(simulations, village_distances) {
  summaries <- lapply(
    simulations, simulation_layer_summaries,
    village_distances = village_distances
  )
  setNames(lapply(c("year_round", "transition", "drought"), function(layer) {
    values <- do.call(rbind, lapply(summaries, `[[`, layer))
    colMeans(values, na.rm = TRUE)
  }), c("year_round", "transition", "drought"))
}
