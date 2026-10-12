#!/usr/bin/env Rscript
# Shared data, availability, and graph helpers for Module 02.

suppressPackageStartupMessages({
  library(igraph)
  library(sf)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

month_column_name <- function(month, year) paste0(month, "_", year)

extract_year_ndvi_matrix <- function(data, year) {
  columns <- month_column_name(month.abb, year)
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0L) {
    stop("Missing NDVI columns for ", year, ": ", paste(missing, collapse = ", "))
  }

  matrix <- as.matrix(sf::st_drop_geometry(data)[, columns, drop = FALSE])
  rownames(matrix) <- data$location_name
  colnames(matrix) <- month.abb
  matrix
}

# --- NDVI-derived resource availability ------------------------------------

# Map NDVI to relative resource availability on [minimum, 1]. The midpoint
# and slope are fixed before fitting movement parameters.
availability_curve <- function(ndvi, midpoint, slope, minimum = 0.02) {
  minimum <- pmin(0.95, pmax(0, minimum))
  slope <- pmax(slope, .Machine$double.eps)
  minimum + (1 - minimum) * stats::plogis((ndvi - midpoint) / slope)
}

calibrate_availability_curve <- function(resource_ndvi, years,
                                         reference_months,
                                         target_unavailable,
                                         minimum = 0.02,
                                         slope = NA_real_) {
  reference_months <- intersect(reference_months, month.abb)
  if (length(reference_months) == 0L) {
    stop("At least one valid reference month is required.")
  }

  columns <- as.vector(outer(reference_months, years, month_column_name))
  columns <- intersect(columns, names(resource_ndvi))
  dry_ndvi <- unlist(
    sf::st_drop_geometry(resource_ndvi)[, columns, drop = FALSE],
    use.names = FALSE
  )
  dry_ndvi <- dry_ndvi[is.finite(dry_ndvi)]
  if (length(dry_ndvi) == 0L) stop("No finite NDVI values in reference months.")

  if (!is.finite(slope) || slope <= 0) {
    quantiles <- stats::quantile(
      dry_ndvi, c(0.10, 0.90), names = FALSE, na.rm = TRUE
    )
    slope <- max(diff(quantiles) / 6, 0.015)
  }

  target_available <- 1 - pmin(0.95, pmax(0.01, target_unavailable))
  difference_from_target <- function(midpoint) {
    mean(availability_curve(dry_ndvi, midpoint, slope, minimum)) -
      target_available
  }
  midpoint <- stats::uniroot(
    difference_from_target,
    interval = range(dry_ndvi) + c(-1, 1)
  )$root
  realised <- availability_curve(dry_ndvi, midpoint, slope, minimum)

  list(
    midpoint = midpoint,
    slope = slope,
    minimum = minimum,
    target_unavailable = target_unavailable,
    realised_reference_unavailable = mean(1 - realised),
    reference_months = reference_months,
    reference_years = years
  )
}

# Calibrate separate availability curves for grazing and watering sites using
# the working dry-season availability deficits in the reference months.
calibrate_resource_availability_params <- function(resource_ndvi, years,
                                                   reference_months,
                                                   target_unavailable,
                                                   minimum = 0.02,
                                                   slope = NA_real_) {
  parameters <- lapply(names(resource_ndvi), function(resource) {
    target <- if (is.list(target_unavailable)) {
      target_unavailable[[resource]]
    } else {
      target_unavailable
    }
    if (is.null(target)) {
      stop("Missing unavailable-resource target for ", resource, ".")
    }

    calibrate_availability_curve(
      resource_ndvi[[resource]],
      years = years,
      reference_months = reference_months,
      target_unavailable = target,
      minimum = minimum,
      slope = slope
    )
  })
  names(parameters) <- names(resource_ndvi)
  parameters
}

resource_availability_matrix <- function(resource_ndvi, year, parameters) {
  ndvi <- extract_year_ndvi_matrix(resource_ndvi, year)
  availability <- availability_curve(
    ndvi,
    midpoint = parameters$midpoint,
    slope = parameters$slope,
    minimum = parameters$minimum
  )
  dimnames(availability) <- dimnames(ndvi)
  availability
}

resource_availability_cycle_mean_matrix <- function(resource_ndvi, years,
                                                    parameters) {
  matrices <- lapply(years, function(year) {
    resource_availability_matrix(resource_ndvi, year, parameters)
  })
  Reduce(`+`, matrices) / length(matrices)
}

resource_availability_inputs_from_matrix <- function(availability_matrix,
                                                     regimes,
                                                     resource_names) {
  availability_matrix <- availability_matrix[
    resource_names, month.abb, drop = FALSE
  ]
  annual_mean <- rowMeans(availability_matrix, na.rm = TRUE)

  by_month <- setNames(lapply(month.abb, function(month) {
    availability_matrix[, month]
  }), month.abb)

  by_regime <- setNames(lapply(c("routine", "transition", "drought"), function(regime) {
    months <- names(regimes)[regimes == regime]
    if (length(months) == 0L) return(annual_mean)
    rowMeans(availability_matrix[, months, drop = FALSE], na.rm = TRUE)
  }), c("routine", "transition", "drought"))

  list(
    matrix = availability_matrix,
    annual_mean = annual_mean,
    by_month = by_month,
    by_regime = by_regime
  )
}

combined_resource_availability <- function(resource_ndvi, year,
                                           availability_parameters) {
  by_resource <- lapply(names(resource_ndvi), function(resource) {
    matrix <- resource_availability_matrix(
      resource_ndvi[[resource]], year, availability_parameters[[resource]]
    )
    setNames(colMeans(matrix, na.rm = TRUE), month.abb)
  })
  names(by_resource) <- names(resource_ndvi)

  list(
    combined = Reduce(`+`, by_resource) / length(by_resource),
    by_resource = by_resource
  )
}

# Average site-level availability across the observed climate cycle before
# assigning the calendar-month labels used for fitting observed networks.
combined_resource_availability_cycle_mean <- function(resource_ndvi, years,
                                                      availability_parameters) {
  by_resource <- lapply(names(resource_ndvi), function(resource) {
    matrix <- resource_availability_cycle_mean_matrix(
      resource_ndvi[[resource]], years, availability_parameters[[resource]]
    )
    setNames(colMeans(matrix, na.rm = TRUE), month.abb)
  })
  names(by_resource) <- names(resource_ndvi)

  list(
    combined = Reduce(`+`, by_resource) / length(by_resource),
    by_resource = by_resource
  )
}

calibrate_availability_regime_thresholds_cycle_current <- function(
    resource_ndvi, availability_parameters, years,
    drought_percentile, transition_percentile) {
  signals <- unlist(lapply(years, function(year) {
    combined_resource_availability(
      resource_ndvi, year, availability_parameters
    )$combined
  }), use.names = FALSE)

  stats::setNames(
    as.numeric(stats::quantile(
      signals,
      probs = c(drought_percentile, transition_percentile),
      names = FALSE,
      na.rm = TRUE
    )),
    c("drought", "transition")
  )
}

# Apply the fixed drought and transition thresholds to a named monthly series.
assign_regime_signal_three <- function(signal, thresholds) {
  regimes <- vapply(signal, function(value) {
    if (!is.finite(value)) return("routine")
    if (value < thresholds[["drought"]]) return("drought")
    if (value < thresholds[["transition"]]) return("transition")
    "routine"
  }, character(1))
  names(regimes) <- names(signal)
  regimes
}

# --- One-mode network helpers ----------------------------------------------

edge_key <- function(from, to) paste(sort(c(from, to)), collapse = "--")

edge_keys <- function(graph) {
  if (igraph::ecount(graph) == 0L) return(character(0))
  edges <- igraph::as_edgelist(graph, names = TRUE)
  apply(edges, 1, function(edge) edge_key(edge[1], edge[2]))
}

make_graph <- function(edges, vertices) {
  if (nrow(edges) == 0L) {
    graph <- igraph::make_empty_graph(n = length(vertices), directed = FALSE)
    igraph::V(graph)$name <- vertices
    return(graph)
  }
  igraph::graph_from_data_frame(edges, directed = FALSE, vertices = vertices) |>
    igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE)
}

project_one_mode <- function(bipartite_matrix, village_names) {
  # Villages are connected when they use at least one resource in the same
  # month; edge multiplicity is deliberately discarded at the observed scale.
  adjacency <- bipartite_matrix %*% t(bipartite_matrix)
  diag(adjacency) <- 0
  adjacency[adjacency > 0] <- 1
  graph <- igraph::graph_from_adjacency_matrix(
    adjacency, mode = "undirected", diag = FALSE
  )
  igraph::V(graph)$name <- village_names
  graph
}

edge_distance_vector <- function(graph, village_distances) {
  if (igraph::ecount(graph) == 0L) return(numeric(0))
  edges <- igraph::as_edgelist(graph, names = TRUE)
  from <- match(edges[, 1], rownames(village_distances))
  to <- match(edges[, 2], colnames(village_distances))
  as.numeric(village_distances[cbind(from, to)])
}

jaccard <- function(first, second) {
  if (length(first) == 0L && length(second) == 0L) return(1)
  if (length(first) == 0L || length(second) == 0L) return(0)
  length(intersect(first, second)) / length(union(first, second))
}
