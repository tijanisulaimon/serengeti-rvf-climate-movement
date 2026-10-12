#!/usr/bin/env Rscript
#' Module 03 configuration: paths, defaults, and small helpers shared by
#' the RVF transmission scripts. Upstream inputs now come directly from
#' module 02's projected bipartite bundles plus `workspace_inputs.RData`.

suppressPackageStartupMessages({
  library(readr)
})

source(here::here("_project_paths.R"))
source(here::here("analysis", "_future_climate_config.R"))

# Null coalescing helper used throughout the config and loaders.
`%||%` <- function(x, y) if (is.null(x)) y else x

# Create a directory if needed, then return the path for inline config use.
make_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  path
}

# -------------------------------------------------------------------------
# Publication-run settings
# -------------------------------------------------------------------------
# These are the only routine execution choices a user should need to edit.
# The numbered scripts run with these values when called with plain `Rscript`.
# Keep exploratory changes on a separate branch or restore these values before
# regenerating manuscript outputs.
publication_network_replicates <- seq_len(20L)
publication_reporting_years <- 40L

rvf_execution <- list(
  max_cores = 8L,
  # All publication comparisons use the same final 40-year duration. This is
  # four complete observed 10-year cycles or two complete future 20-year cycles.
  reporting_years = publication_reporting_years,
  # Module 02 materialised 100 network realisations for every observed and
  # future climate year. Module 03 uses the first 20 throughout the publication
  # rerun so uncertainty summaries are directly comparable across questions.
  network_replicates = publication_network_replicates,
  analysis_replicates = publication_network_replicates,
  attribution_replicates = publication_network_replicates,
  assemble_datasets = c("observed", "future"),
  use_calibrated_profiles = TRUE,
  require_calibration_target = TRUE,
  strict_preflight = FALSE,
  calibration = list(
    profiles = c(
      "peridomestic_dominant",
      "balanced",
      "grazing_floodplain_dominant",
      "resource_edge_dominant",
      "waterpoint_extreme"
    ),
    rainfall_models = c("saturating_rainfall", "quadratic_rainfall"),
    replicates = publication_network_replicates,
    years_burn = NULL,              # NULL inherits the model setting below
    years_endemic = NULL,
    terminal_window_years = publication_reporting_years,
    seed_mode = "all_villages",
    n_seed = 1L,
    seed_village = 10L,             # used only for an explicit invasion assay
    output_suffix = "",
    resume_existing = FALSE,
    use_checkpoints = TRUE,
    plan_only = FALSE,
    write_headline_figures = TRUE
  ),
  baseline_runs = list(
    list(
      site_profile = "grazing_floodplain_dominant",
      rainfall_model = "saturating_rainfall",
      network_replicates = publication_network_replicates
    ),
    list(
      site_profile = "peridomestic_dominant",
      rainfall_model = "saturating_rainfall",
      network_replicates = publication_network_replicates
    ),
    list(
      site_profile = "grazing_floodplain_dominant",
      rainfall_model = "quadratic_rainfall",
      network_replicates = publication_network_replicates
    )
  ),
  mobility_profiles = c(
    "peridomestic_dominant",
    "balanced",
    "grazing_floodplain_dominant",
    "resource_edge_dominant",
    "waterpoint_extreme"
  ),
  future_runs = list(
    list(
      site_profile = "grazing_floodplain_dominant",
      rainfall_model = "saturating_rainfall",
      network_replicates = publication_network_replicates
    ),
    list(
      site_profile = "peridomestic_dominant",
      rainfall_model = "saturating_rainfall",
      network_replicates = publication_network_replicates
    ),
    list(
      site_profile = "grazing_floodplain_dominant",
      rainfall_model = "quadratic_rainfall",
      network_replicates = publication_network_replicates
    )
  ),
  attribution_profiles = c(
    "grazing_floodplain_dominant",
    "peridomestic_dominant"
  ),
  smoke = list(
    mobility = FALSE,
    future = FALSE,
    attribution = FALSE
  ),
  reuse = list(
    calibration = FALSE,
    observed_baseline = FALSE,
    attribution = FALSE
  ),
  output = list(
    write_derived_csvs = FALSE,
    write_diagnostic_plots = FALSE,
    save_run_rds = TRUE,
    save_baseline_run_rds = FALSE,
    save_baseline_state = FALSE,
    save_heavy_diagnostics = FALSE
  ),
  ode = list(
    rtol = 1e-6,
    atol = 1e-8,
    maxsteps = 50000L
  )
)

# Cap local parallel work so long replicate sweeps stay responsive on a laptop.
rvf_max_cores <- function(default = rvf_execution$max_cores, n_tasks = Inf) {
  requested <- suppressWarnings(as.integer(default))
  if (!is.finite(requested) || is.na(requested) || requested < 1L) {
    requested <- rvf_execution$max_cores
  }

  n_tasks <- suppressWarnings(as.integer(n_tasks))
  if (!is.finite(n_tasks) || is.na(n_tasks) || n_tasks < 1L) n_tasks <- requested

  max(1L, min(as.integer(requested), rvf_execution$max_cores, as.integer(n_tasks)))
}

rvf_parallel_lapply <- function(X, FUN, ..., n_cores = rvf_max_cores(n_tasks = length(X)),
                                mc_preschedule = FALSE) {
  if (length(X) == 0L) return(list())
  n_cores <- rvf_max_cores(default = n_cores, n_tasks = length(X))

  if (.Platform$OS.type == "windows" || n_cores <= 1L) {
    return(lapply(X, FUN, ...))
  }

  # mclapply reports a worker process that stops without a result, or a task
  # that fails, only as a warning, and returns NULL or a try-error in its place.
  # Stop instead, so a lost replicate can never drop silently out of a summary.
  results <- withCallingHandlers(
    parallel::mclapply(X, FUN, ..., mc.cores = n_cores, mc.preschedule = mc_preschedule),
    warning = function(w) {
      if (grepl("did not deliver|resulted in an error|encountered error",
                conditionMessage(w))) {
        stop("A parallel task failed or its worker process stopped (",
             conditionMessage(w), "). Rerun the script.", call. = FALSE)
      }
    }
  )
  failed <- vapply(results, inherits, logical(1), "try-error")
  if (any(failed)) {
    stop(sum(failed), " of ", length(X), " parallel task(s) failed. First error: ",
         conditionMessage(attr(results[[which(failed)[1]]], "condition")), call. = FALSE)
  }
  results
}

rvf_write_derived_csvs <- function() {
  isTRUE(rvf_execution$output$write_derived_csvs)
}

rvf_write_diagnostic_plots <- function() {
  isTRUE(rvf_execution$output$write_diagnostic_plots)
}

rvf_save_heavy_diagnostics <- function() {
  isTRUE(rvf_execution$output$save_heavy_diagnostics)
}

rvf_config <- local({
  # Root of module 03's data and outputs — wrapped in make_dir() so a
  # fresh checkout has all required directories on first source().
  data_dir <- make_dir(proj$rvf_data)
  processed_dir <- make_dir(file.path(data_dir, "processed"))
  outputs_dir <- make_dir(file.path(proj$rvf_dir, "outputs"))
  intermediate_dir <- make_dir(file.path(outputs_dir, "intermediate"))
  results_dir <- make_dir(file.path(outputs_dir, "results"))
  figures_dir <- make_dir(file.path(outputs_dir, "figures"))
  logs_dir <- make_dir(file.path(outputs_dir, "logs"))
  cache_dir <- make_dir(file.path(outputs_dir, "cache"))

  future_windows <- list(
    `2050` = 2041:2060,
    `2080` = 2061:2080,
    `2100` = 2081:2100
  )
  future_years <- sort(unique(unlist(future_windows, use.names = FALSE)))
  primary_rainfall_model <- "saturating_rainfall"
  # Recent decade of CHELSA-monthly historical coverage (1981-2018) used as
  # the "current baseline" for cycling. Was 2016:2022 (CHIRPS-era); switched
  # in 2026-05 when historical climate moved to CHELSA-monthly.
  baseline_years <- 2009:2018
  list(
    execution = rvf_execution,
    # --- module directories -----------------------------------------------
    data_dir      = data_dir,
    processed_dir = processed_dir,
    networks_dir  = make_dir(file.path(processed_dir, "networks")),
    climate_dir   = make_dir(file.path(processed_dir, "climate")),
    scenarios_dir = make_dir(file.path(processed_dir, "scenarios")),
    # Generated-output contract. Expensive or nested R objects are kept apart
    # from compact tables so users can inspect or publish `results/` without
    # copying simulation caches. `outputs_dir` is retained as the common root;
    # active workflows use the five category-specific paths below.
    outputs_dir      = outputs_dir,
    intermediate_dir = intermediate_dir,
    results_dir      = results_dir,
    figures_dir      = figures_dir,
    logs_dir         = logs_dir,
    cache_dir        = cache_dir,
    notes_dir     = file.path(proj$rvf_dir, "notes"),

    # --- inputs from module 02 --------------------------------------------
    # Module 03 reads Module 02's saved outputs directly. The bundle index
    # lists each scenario bundle with its SSP metadata.
    mvt_outputs_dir = proj$mvt_outputs,
    mvt_inputs_file = file.path(proj$mvt_outputs, "inputs", "observed_inputs.rds"),
    mvt_bundle_index = file.path(
      proj$mvt_outputs, "projections", "summaries", "network_bundle_index.csv"
    ),
    mvt_bundles_dir = file.path(proj$mvt_outputs, "projections", "bundles"),

    # --- observed climate data (single source, not per-GCM) ---------------
    # Historical climate is now CHELSA-monthly (Karger et al. 2025) across
    # all four variables (tas/tasmax/tasmin/pr). One dataset across all
    # scenarios. Pairs cleanly with the CHELSA-CMIP6 future projections
    # (delta-change downscaling from the same CHELSA V2.1 backbone).
    # Was CHIRPS + CHIRTS-ERA5 before 2026-05.
    observed_pr_file   = file.path(proj$root, "data", "processed",
                                   "chelsa_historical",
                                   "CHELSA_monthly_pr_1981_2018.nc"),
    observed_tas_file  = file.path(proj$root, "data", "processed",
                                   "chelsa_historical",
                                   "CHELSA_monthly_tas_1981_2018.nc"),
    observed_tmax_file = file.path(proj$root, "data", "processed",
                                   "chelsa_historical",
                                   "CHELSA_monthly_tasmax_1981_2018.nc"),
    observed_tmin_file = file.path(proj$root, "data", "processed",
                                   "chelsa_historical",
                                   "CHELSA_monthly_tasmin_1981_2018.nc"),
    future_scenarios = future_climate$scenarios,
    primary_future_scenario = future_climate$primary_scenario,
    scenario_labels = future_climate$scenario_labels,
    chelsa_future_root = file.path(
      proj$root, "data", "processed", "chelsa_cmip6_future"
    ),

    # --- scenario / parameter defaults ------------------------------------
    current_years = baseline_years,                   # observed years used for current baseline cycling
    future_windows = future_windows,                  # exemplar label -> climate years
    future_years  = future_years,                     # all future years to extract / assemble
    # Use the same four-GCM ensemble for both SSPs. MPI-ESM1-2-LR is retained
    # upstream for quality assessment but excluded from reporting because its
    # SSP2-4.5 mean-temperature channel failed the internal consistency check.
    climate_models = future_climate$reporting_gcms,

    # Reference year used by 01_prepare_inputs.R for any
    # threshold/standardisation against observed climate; separate from
    # module 02's reference year for regime classification (2017), which
    # is baked into 02's saved bundles.
    threshold_reference_year = 2016L,

    # Module 02 saves 100 networks per scenario; Module 03 materialises and
    # analyses the first 20 for the publication rerun.
    network_replicates = rvf_execution$network_replicates,
    analysis_network_replicates = rvf_execution$analysis_replicates,
    # Backward-compatible fallback for scripts / tests that still want one
    # deterministic replicate explicitly.
    representative_replicate = 1L,

    # Raster-to-point extraction method for climate values at
    # village/grazing/watering coordinates (01_prepare_inputs.R).
    extraction_method = "bilinear",

    # Rainfall-only carrying-capacity responses. Temperature affects biting,
    # mortality, and extrinsic incubation elsewhere in the model, but does not
    # also multiply habitat capacity. This avoids counting an uncertain
    # cross-species thermal-density response twice.
    rainfall_models = list(
      quadratic_rainfall = list(
        label = "Quadratic rainfall response",
        min_K = 0.05,
        rainfall_lower = 20,
        rainfall_upper = 500,
        reference = "Quadratic response adapted from Caldwell et al. 2021 to monthly Serengeti rainfall"
      ),
      saturating_rainfall = list(
        label = "Saturating rainfall response",
        min_K = 0.05,
        R_half = 60,
        n_hill = 1.25,
        reference = "Provisional pooled-vector response without an assumed high-rainfall decline"
      )
    ),
    primary_rainfall_model = primary_rainfall_model,

    # --- site-capacity profiles ------------------------------------------
    # Maximum vector:host ratios at saturated carrying capacity are not
    # identifiable from the available literature. Treat them as named
    # ecological profiles, each specifying relative vector pressure at
    # villages, grazing sites, and watering sites. The helper in
    # `R/model_parameters.R` rescales the ratios to the same mean m_max
    # as the prior 2:6:30 reference (mean = 12.67), so profile controls
    # spatial allocation while `capacity_scale` controls overall abundance.
    #
    # *Profiles are tuned to span the exposure-weighted FOI share space*, not
    # the raw m_max-ratio space. The default baseline now assumes every herd
    # follows the daily village-grazing-watering itinerary (p_mobile = 1) with
    # a 15 h / 8 h / 1 h exposure budget. Per-hour effective vector pressure
    # must therefore be much higher at short-duration water interfaces before
    # water can dominate total exposure. The ratios below invert this exposure
    # budget to hit the indicated mobile-herd exposure shares.
    site_profile_base_mean_m_max = mean(c(2, 6, 30)),
    site_profiles = list(
      peridomestic_dominant = list(
        label = "Village-associated exposure (75/20/5)",
        ratios = c(village = 1, grazing = 0.50, watering = 1.00),
        role = "peridomestic upper bound; reduced movement should tend to amplify residence exposure"
      ),
      balanced = list(
        label = "Balanced exposure (33/33/33)",
        ratios = c(village = 1, grazing = 1.875, watering = 15.0),
        role = "neutral exposure midpoint after applying the 15/8/1 h time budget"
      ),
      grazing_floodplain_dominant = list(
        label = "Grazing/floodplain-dominant exposure (20/65/15)",
        ratios = c(village = 1, grazing = 6.09375, watering = 11.25),
        role = "agropastoral/floodplain Aedes-and-pasture exposure"
      ),
      resource_edge_dominant = list(
        label = "Resource-edge-dominant exposure (10/55/35)",
        ratios = c(village = 1, grazing = 10.3125, watering = 52.5),
        role = "water-associated vector habitat spilling into grazing/floodplain exposure"
      ),
      waterpoint_extreme = list(
        label = "Waterpoint-extreme exposure (10/20/70)",
        ratios = c(village = 1, grazing = 3.75, watering = 105),
        role = "stress-test upper bound where one-hour water exposure dominates"
      )
    ),
    site_profile_aliases = list(
      village_extreme = "peridomestic_dominant",
      grazing_dominant = "grazing_floodplain_dominant",
      water_dominant = "resource_edge_dominant",
      water_extreme = "waterpoint_extreme"
    ),

    # --- simulation defaults (used by scripts 02-10) ----------------------
    # One place to edit the ODE / calibration knobs. Scripts read
    # `cfg$sim$field` rather than using env vars. Override in a script
    # by redefining just the field(s) you need.
    sim = list(
      baseline_year  = 2017L,     # stationary baseline = module 02's network ref year
      baseline_years = baseline_years,
      years_burn     = 3L,        # disease-free vector burn-in only. The host starts at its
                                  #   all-susceptible disease-free equilibrium and the vector
                                  #   pool re-equilibrates within one season, so 2-3 y suffices;
                                  #   the long transient that matters (immunity build-up) is
                                  #   entirely post-seeding and is unaffected by this length.
      years_endemic  = 150L,      # post-seed horizon for recurrent endemic dynamics.
                                  #   Each analysis checks matched forcing cycles rather than
                                  #   assuming convergence to a constant equilibrium.
      summary_burnin_years = 110L, # retained as a fallback for older calibration code;
                                   #   the publication reporting duration is now set
                                   #   explicitly under execution$reporting_years.
      output_step    = 7,         # weekly output (days between ODE samples)
      host_population_per_village = 5000L,
      p_mobile       = 1.0,       # baseline assumes all cattle follow the daily resource itinerary
      n_seed         = 1L,        # default: one infectious cow per seeded village after vector burn-in
      seed_village   = "all",     # default capacity assay seeds all villages; invasion tests can override
      external_introduction_per_year = 0, # system-wide exogenous infection pressure; sensitivity only
      external_introduction_sensitivity = c(0, 0.1, 0.5, 1),
      transmission_scale = 1.0,   # fixed beta_hv / beta_vh literature baseline; not tuned in the new calibration
      # FOI-targeted profile set (2026-05-12 redesign): five profiles now span
      # the exposure-weighted FOI share space rather than the raw m_max-ratio space, so
      # peridomestic_dominant / balanced / grazing_floodplain_dominant / resource_edge_dominant /
      # waterpoint_extreme test genuinely distinct transmission ecologies after
      # cattle residence time is applied. The grazing/floodplain/resource-
      # associated archetype remains the primary working profile. With
      # p_mobile = 1, reductions in p_mobile are now
      # interpreted as increasing the non-commuting village-resident fraction.
      #
      # Scripts 03, 04, 06, and 08 read the canonical calibration table by
      # default. The fallback below is used only before a calibration exists or
      # for an explicitly uncalibrated diagnostic.
      capacity_scale     = 0.49, # fallback only; script 02's calibration table is authoritative
      site_profile       = "grazing_floodplain_dominant",
      # Future climate comparisons use two prespecified contrasting exposure
      # profiles. The full five-profile set is reserved for the mobility
      # experiment rather than multiplying every future scenario by all five
      # latent vector landscapes.
      future_site_profile = "grazing_floodplain_dominant",
      # The canonical calibration directory contains both rainfall responses;
      # downstream helpers select the row matching the configured response.
      calibration_source_suffix = "",
      # The saturating Hill response is the primary pooled-vector formulation;
      # the quadratic response is retained as a structural sensitivity.
      rainfall_model = primary_rainfall_model,
      rainfall_models_to_compare = c(
        "saturating_rainfall",
        "quadratic_rainfall"
      ),
      sero_target_min = 4,        # lower inter-epidemic cattle anchor: de Glanville et al. 2022, PLoS Negl Trop Dis 16:e0010871
      sero_target_max = 10,       # pragmatic upper bound spanning Salekwa et al. 2019, Am J Trop Med Hyg 101:459-466 (7.7% overall)
      sero_target_mid = 7,        # ranking midpoint inside the literature-informed 4-10% target band
      waning_years   = Inf,       # baseline = no waning; finite waning is sensitivity-only if needed for fit
      network_mode   = "monthly"  # "monthly" or "annualized"
    )
  )
})

# Create the standard module-03 output tree for a clean checkout.
ensure_rvf_dirs <- function(cfg = rvf_config) {
  dirs <- c(cfg$data_dir, cfg$processed_dir, cfg$networks_dir,
            cfg$climate_dir, cfg$scenarios_dir, cfg$outputs_dir,
            cfg$intermediate_dir, cfg$results_dir, cfg$figures_dir,
            cfg$logs_dir, cfg$cache_dir)
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

# Load module 02's site objects into a private environment so callers get
# the inputs they need without polluting the workspace.
load_workspace_inputs <- function(cfg = rvf_config) {
  inputs <- readRDS(cfg$mvt_inputs_file)
  list(
    village_ndvi_observed = inputs$ndvi$village,
    grazing_ndvi_observed = inputs$ndvi$grazing,
    watering_ndvi_observed = inputs$ndvi$watering
  )
}

# Drop time-varying raster columns and keep static site metadata only.
static_site_metadata <- function(sf_obj) {
  df <- sf::st_drop_geometry(sf_obj)
  keep <- !grepl("^[A-Z][a-z]{2}_[0-9]{4}$", names(df))
  df[, keep, drop = FALSE]
}
