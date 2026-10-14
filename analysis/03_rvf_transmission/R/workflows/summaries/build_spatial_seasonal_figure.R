#!/usr/bin/env Rscript

# Build the main spatial-broadening figure from the summaries produced earlier
# in script 09. One absolute current map provides the reference distribution;
# future maps show the change from that reference so spatial broadening does
# not have to be inferred by comparing colours across seven absolute maps. A
# separate concentration panel quantifies redistribution across villages.
# No transmission simulations are run here.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(sf)
  library(tidyr)
})

source(here::here("analysis", "03_rvf_transmission", "R", "transmission_model.R"))
source(here::here("analysis", "03_rvf_transmission", "R", "plot_helpers.R"))

cfg <- load_rvf_config()
sim <- cfg$sim
primary_profile <- rvf_resolve_site_profile_id(
  sim$future_site_profile %||% "grazing_floodplain_dominant"
)

spatial_dir <- file.path(cfg$results_dir, "spatial_patterns")
figures_dir <- cfg$figures_dir
village_path <- file.path(spatial_dir, "village_burden_summary.csv")
concentration_path <- file.path(
  spatial_dir,
  "spatial_concentration_summary.csv"
)
required_paths <- c(village_path, concentration_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop(
    "Missing spatial summary input(s):\n  ",
    paste(missing_paths, collapse = "\n  ")
  )
}

village_summary <- read_csv(village_path, show_col_types = FALSE)
concentration_summary <- read_csv(
  concentration_path,
  show_col_types = FALSE
)
if (!all(c("annual_infections_mid", "annual_share_mid") %in% names(village_summary))) {
  stop(
    "The village summary lacks annual infection shares. ",
    "Rerun script 09 with the current spatial-summary code."
  )
}

# Village IDs follow the order stored in the processed scenario bundles.
baseline_scenario <- load_rvf_scenario(
  rvf_with_replicate(
    sprintf("OBSERVED_CURRENT_%d", sim$baseline_years[[1]]),
    sort(unique(cfg$analysis_network_replicates))[[1]]
  ),
  cfg = cfg,
  validate = TRUE
)

village_coordinates <- read_csv(
  here::here(
    "data",
    "shapefiles",
    "resource_areas",
    "serengeti village centre xy.csv"
  ),
  show_col_types = FALSE
) |>
  setNames(c("village_name", "easting", "northing")) |>
  right_join(
    tibble(
      village_id = seq_along(baseline_scenario$vil_names),
      village_name = baseline_scenario$vil_names
    ),
    by = "village_name"
  )

if (nrow(village_coordinates) != 94L ||
    anyNA(village_coordinates[c("easting", "northing")])) {
  stop("Village coordinates do not match the 94 modelled villages.")
}

district_boundary <- st_read(
  here::here(
    "data",
    "shapefiles",
    "resource_areas",
    "Serengeti and Ngorongoro district shape files.shp"
  ),
  quiet = TRUE
) |>
  filter(.data$District_N == "Serengeti") |>
  summarise() |>
  st_transform(32736)

# Serengeti National Park and the Ikorongo and Grumeti game reserves cover about
# seven tenths of the district and hold no village herds, which is why the
# modelled villages sit in the northwest. Drawing them stops a reader taking the
# empty part of the map for missing data. Built by
# analysis/prepare_protected_areas.R from the World Database on Protected
# Areas, which requires attribution in the caption.
protected_areas <- st_read(
  here::here("data", "shapefiles", "serengeti_protected_areas.gpkg"),
  quiet = TRUE
) |>
  st_transform(32736)

# -------- A-B. Observed reference and all six future maps ----------------
profile_map_data <- village_summary |>
  filter(.data$site_profile == primary_profile) |>
  mutate(
    annual_infections_per_1000 =
      1000 * .data$annual_infections_mid /
        as.numeric(sim$host_population_per_village)
  ) |>
  left_join(village_coordinates, by = "village_id")

current_map_data <- profile_map_data |>
  filter(.data$scenario == "current")

current_reference <- current_map_data |>
  select(
    "village_id",
    current_infections_per_1000 = "annual_infections_per_1000"
  )

future_map_data <- profile_map_data |>
  filter(.data$scenario %in% cfg$future_scenarios) |>
  left_join(current_reference, by = "village_id") |>
  mutate(
    change_per_1000 =
      .data$annual_infections_per_1000 - .data$current_infections_per_1000,
    scenario_label = factor(
      rvf_ssp_label(.data$scenario),
      levels = unname(cfg$scenario_labels[cfg$future_scenarios])
    ),
    time_label = factor(
      sprintf("%d-%d", .data$window_start, .data$window_end),
      levels = vapply(
        cfg$future_windows,
        function(years) sprintf("%d-%d", min(years), max(years)),
        character(1)
      )
    )
  )

if (nrow(current_map_data) != 94L ||
    nrow(future_map_data) != 94L * length(cfg$future_scenarios) *
      length(cfg$future_windows)) {
  stop("The spatial maps do not contain the expected village-scenario cells.")
}
if (any(!is.finite(future_map_data$change_per_1000))) {
  stop("Future-minus-current village changes contain non-finite values.")
}

# Observed burden and future changes have separate scales because they answer
# different questions. The same change scale is retained across all six
# future maps.
current_limit <- ceiling(max(current_map_data$annual_infections_per_1000) / 5) * 5
current_breaks <- seq(0, current_limit, by = 5)
change_limit <- ceiling(max(future_map_data$change_per_1000) / 10) * 10
change_limit <- max(change_limit, 10)
change_breaks <- seq(0, change_limit, by = 10)

p_current <- ggplot() +
  geom_sf(
    data = district_boundary,
    fill = "white",
    colour = "grey25",
    linewidth = 0.5
  ) +
  geom_sf(
    data = protected_areas,
    aes(fill = "Protected area"),
    colour = "grey70",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = c("Protected area" = "grey88"), name = NULL) +
  ggnewscale::new_scale_fill() +
  geom_point(
    data = current_map_data,
    # Filled circles with a thin outline: pale low values stay visible against
    # the white district fill instead of disappearing into it.
    aes(
      x = .data$easting,
      y = .data$northing,
      fill = .data$annual_infections_per_1000
    ),
    shape = 21,
    colour = "grey30",
    stroke = 0.12,
    size = 2.1
  ) +
  scale_fill_viridis_c(
    # Absolute burden (this panel) and future change (next panel) are distinct
    # quantities on different ranges, so they use different hue families: blue
    # here, red there. Both ramps run pale -> dark, so darker always means more
    # and only the hue signals which quantity is being shown.
    option = "mako",
    direction = -1,
    limits = c(0, current_limit),
    breaks = current_breaks,
    labels = scales::label_number(accuracy = 1),
    name = "Simulated annual infections per 1,000 cattle"
  ) +
  coord_sf(crs = 32736, datum = NA) +
  labs(title = "Observed climate", x = NULL, y = NULL) +
  theme_void(base_size = 11.5) +
  theme(
    text = element_text(colour = "black"),
    plot.title = element_text(
      colour = "black",
      size = 11,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 3)
    ),
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 8.5),
    legend.title = element_text(colour = "black", size = 9)
  ) +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      barwidth = grid::unit(3.2, "cm"),
      barheight = grid::unit(0.25, "cm")
    )
  )

p_future <- ggplot() +
  geom_sf(
    data = district_boundary,
    fill = "white",
    colour = "grey25",
    linewidth = 0.45
  ) +
  geom_sf(
    data = protected_areas,
    fill = "grey88",
    colour = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = future_map_data,
    aes(
      x = .data$easting,
      y = .data$northing,
      fill = .data$change_per_1000
    ),
    shape = 21,
    colour = "grey30",
    stroke = 0.1,
    size = 1.75
  ) +
  facet_grid(scenario_label ~ time_label) +
  scale_fill_gradientn(
    # All village-level future changes are non-negative. A sequential scale
    # therefore uses the full colour range and preserves the common 0-based
    # reference across all six future maps. The yellow-orange-red progression
    # gives stronger perceptual separation than a single-hue red ramp, while
    # keeping red reserved for the largest increases.
    colours = c(
      "#FFFFCC", "#FED976", "#FEB24C", "#FD8D3C",
      "#FC4E2A", "#E31A1C", "#BD0026"
    ),
    limits = c(0, change_limit),
    breaks = change_breaks,
    labels = scales::label_number(accuracy = 1),
    name = "Increase from observed-climate baseline (per 1,000 cattle)"
  ) +
  coord_sf(crs = 32736, datum = NA) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 11.5) +
  theme(
    text = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10.2, face = "bold"),
    strip.text.y.right = element_text(angle = 0),
    panel.spacing = grid::unit(0.14, "cm"),
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 9),
    legend.title = element_text(colour = "black", size = 9.5),
    legend.margin = margin(t = -2)
  ) +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      barwidth = grid::unit(6.2, "cm"),
      barheight = grid::unit(0.28, "cm")
    )
  )

# -------- C. Intuitive concentration summary across both profiles -------
# The top ten villages are the ceiling of 10% of 94 villages. If burden were
# perfectly even, they would contribute 10/94 = 10.6% of annual infections.
profile_labels <- c(
  grazing_floodplain_dominant = "Shared grazing/floodplain exposure",
  peridomestic_dominant = "Village-associated exposure"
)
time_levels <- c(
  "Observed climate",
  vapply(
    cfg$future_windows,
    function(years) sprintf("%d-%d", min(years), max(years)),
    character(1)
  )
)

current_concentration <- concentration_summary |>
  filter(.data$scenario == "current") |>
  mutate(time_label = "Observed climate")

future_concentration <- concentration_summary |>
  filter(.data$scenario %in% cfg$future_scenarios) |>
  mutate(time_label = sprintf("%d-%d", .data$window_start, .data$window_end))

concentration_data <- bind_rows(current_concentration, future_concentration) |>
  mutate(
    profile_label = factor(
      unname(profile_labels[.data$site_profile]),
      levels = unname(profile_labels[c(
        "grazing_floodplain_dominant",
        "peridomestic_dominant"
      )])
    ),
    scenario_label = factor(
      rvf_ssp_label(.data$scenario),
      levels = unname(cfg$scenario_labels[cfg$future_scenarios])
    ),
    time_label = factor(.data$time_label, levels = time_levels)
  )

current_points <- concentration_data |>
  filter(.data$scenario == "current")
future_points <- concentration_data |>
  filter(.data$scenario %in% cfg$future_scenarios)
dodge_ssp <- position_dodge(width = 0.42)

p_concentration <- ggplot() +
  geom_hline(
    yintercept = 10 / 94,
    colour = "grey55",
    linewidth = 0.45
  ) +
  geom_linerange(
    data = current_points,
    aes(
      x = .data$time_label,
      ymin = .data$top_decile_share_lo,
      ymax = .data$top_decile_share_hi
    ),
    colour = "black",
    linewidth = 0.65
  ) +
  geom_point(
    data = current_points,
    aes(x = .data$time_label, y = .data$top_decile_share_mid),
    colour = "black",
    shape = 18,
    size = 3
  ) +
  geom_linerange(
    data = future_points,
    aes(
      x = .data$time_label,
      ymin = .data$top_decile_share_lo,
      ymax = .data$top_decile_share_hi,
      colour = .data$scenario_label,
      group = .data$scenario_label
    ),
    position = dodge_ssp,
    linewidth = 0.65
  ) +
  geom_point(
    data = future_points,
    aes(
      x = .data$time_label,
      y = .data$top_decile_share_mid,
      colour = .data$scenario_label,
      shape = .data$scenario_label,
      group = .data$scenario_label
    ),
    position = dodge_ssp,
    size = 2.65
  ) +
  facet_wrap(vars(profile_label), nrow = 1) +
  scale_colour_ssp() +
  scale_shape_manual(values = c("SSP2-4.5" = 16, "SSP3-7.0" = 17), name = NULL) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    breaks = seq(0.10, 0.40, by = 0.05)
  ) +
  coord_cartesian(ylim = c(0.08, 0.38)) +
  labs(
    x = NULL,
    y = "Share of district infections in\n10 highest-burden villages",
    colour = NULL
  ) +
  theme_classic(base_size = 11.5) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black", size = 9.5),
    axis.text.x = element_text(angle = 18, hjust = 1),
    axis.title = element_text(colour = "black", size = 11),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10.5, face = "bold"),
    legend.position = "top",
    legend.text = element_text(colour = "black", size = 9.5),
    legend.margin = margin(b = -3),
    plot.margin = margin(t = 2)
  )

map_row <- p_current + p_future +
  plot_layout(widths = c(1, 3.05))

spatial_broadening_figure <- map_row / p_concentration +
  plot_layout(heights = c(1.65, 1)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(face = "bold", size = 14))
  )

rvf_ggsave(
  file.path(figures_dir, "06_spatial_broadening.png"),
  spatial_broadening_figure,
  width = 12,
  height = 9.2
)

message("Spatial-broadening figure written to: ", figures_dir)
