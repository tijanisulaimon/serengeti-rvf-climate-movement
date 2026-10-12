#!/usr/bin/env Rscript
#' Plot the climate functions used by the transmission model.
#'
#' The figures are suitable for supplementary material. They contain no
#' manuscript-specific titles; interpretation belongs in the figure captions.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tidyr)
})

source(here::here("analysis", "03_rvf_transmission", "R", "model_parameters.R"))

cfg <- load_rvf_config()
output_dir <- file.path(cfg$results_dir, "diagnostics", "climate_vector_functions")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

climate_index <- read_csv(
  file.path(cfg$climate_dir, "climate_input_index.csv"),
  show_col_types = FALSE
)
# Find each table inside this copy's climate folder rather than at the stored
# path, so the index still works after the project is moved.
climate_files <- file.path(
  cfg$climate_dir, basename(dirname(climate_index$file)), basename(climate_index$file)
)
climate <- bind_rows(lapply(climate_files, readRDS)) |>
  mutate(
    period = case_when(
      .data$dataset == "observed" ~ "Observed 2009-2018",
      .data$year <= 2060 ~ "2041-2060",
      .data$year <= 2080 ~ "2061-2080",
      TRUE ~ "2081-2100"
    ),
    pathway = if_else(
      .data$dataset == "observed", "Observed", .data$scenario_label
    )
  )

# These quantiles document the climate range over which the functions are used.
climate_support <- climate |>
  group_by(.data$pathway, .data$period) |>
  summarise(
    rainfall_q01 = quantile(.data$pr, 0.01),
    rainfall_q10 = quantile(.data$pr, 0.10),
    rainfall_median = median(.data$pr),
    rainfall_q90 = quantile(.data$pr, 0.90),
    rainfall_q99 = quantile(.data$pr, 0.99),
    rainfall_max = max(.data$pr),
    temperature_q01 = quantile(.data$tas, 0.01),
    temperature_median = median(.data$tas),
    temperature_q99 = quantile(.data$tas, 0.99),
    temperature_max = max(.data$tas),
    .groups = "drop"
  )
write_csv(climate_support, file.path(output_dir, "climate_function_support.csv"))

# Translate the climate support into the trait values the model implies.
# These are deterministic transformations of monthly mean temperature, not
# fitted or locally observed mosquito traits.
trait_support <- climate_support |>
  transmute(
    .data$pathway,
    .data$period,
    temperature_q01 = .data$temperature_q01,
    temperature_median = .data$temperature_median,
    temperature_q99 = .data$temperature_q99,
    biting_q01 = biting_rate_temp(.data$temperature_q01),
    biting_median = biting_rate_temp(.data$temperature_median),
    biting_q99 = biting_rate_temp(.data$temperature_q99),
    eip_days_q01_temperature = eip_days(.data$temperature_q01),
    eip_days_median_temperature = eip_days(.data$temperature_median),
    eip_days_q99_temperature = eip_days(.data$temperature_q99),
    lifespan_days_q01_temperature = lifespan_days(.data$temperature_q01),
    lifespan_days_median_temperature = lifespan_days(.data$temperature_median),
    lifespan_days_q99_temperature = lifespan_days(.data$temperature_q99)
  )
write_csv(trait_support, file.path(output_dir, "temperature_trait_support.csv"))

trait_levels <- c(
  "Biting rate\n(bites per vector-day)",
  "Extrinsic incubation\n(days)",
  "Adult mortality\n(per day)"
)

temperature_grid <- tibble(temperature = seq(18, 32, by = 0.05)) |>
  mutate(
    `Biting rate\n(bites per vector-day)` = biting_rate_temp(.data$temperature),
    `Extrinsic incubation\n(days)` = eip_days(.data$temperature),
    `Adult mortality\n(per day)` = mu_m_temp(.data$temperature)
  ) |>
  pivot_longer(
    -all_of("temperature"),
    names_to = "trait",
    values_to = "value"
  ) |>
  mutate(trait = factor(.data$trait, levels = trait_levels))

temperature_support <- climate_support |>
  filter(
    (.data$pathway == "Observed") |
      (.data$period == "2081-2100" & .data$pathway %in% c("SSP2-4.5", "SSP3-7.0"))
  ) |>
  mutate(
    support = factor(
      if_else(.data$pathway == "Observed", .data$period,
              paste(.data$pathway, .data$period)),
      levels = c(
        "Observed 2009-2018",
        "SSP2-4.5 2081-2100",
        "SSP3-7.0 2081-2100"
      )
    )
  )

support_colours <- c(
  "Observed 2009-2018" = "#7A7A7A",
  "SSP2-4.5 2081-2100" = "#0072B2",
  "SSP3-7.0 2081-2100" = "#D55E00"
)

temperature_plot <- ggplot() +
  geom_rect(
    data = temperature_support,
    aes(
      xmin = .data$temperature_q01,
      xmax = .data$temperature_q99,
      ymin = -Inf,
      ymax = Inf,
      fill = .data$support
    ),
    # The three ranges overlap across most of the temperature axis, so the
    # fill has to be light enough to stack without hiding the trait curve.
    alpha = 0.18,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = temperature_grid,
    aes(x = .data$temperature, y = .data$value),
    colour = "black",
    linewidth = 0.75
  ) +
  facet_wrap(vars(.data$trait), scales = "free_y", nrow = 1) +
  scale_fill_manual(values = support_colours, drop = FALSE) +
  scale_x_continuous(breaks = seq(18, 32, by = 2)) +
  labs(x = "Monthly mean temperature (degree C)", y = NULL, fill = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(colour = "black", size = 9),
    axis.title = element_text(colour = "black", size = 10),
    legend.position = "top",
    legend.text = element_text(colour = "black", size = 9),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", face = "bold", size = 10),
    panel.spacing.x = grid::unit(1, "lines")
  )

ggsave(
  file.path(cfg$figures_dir, "vector_temperature_traits.pdf"),
  temperature_plot,
  width = 7.2,
  height = 2.8,
  device = grDevices::cairo_pdf
)

rainfall_grid <- rvf_rainfall_response_table(
  rainfall_mm = seq(0, 600, by = 1),
  cfg = cfg
) |>
  mutate(
    rainfall_label = factor(
      .data$rainfall_label,
      levels = vapply(
        cfg$sim$rainfall_models_to_compare,
        rvf_rainfall_model_label,
        character(1),
        cfg = cfg
      )
    )
  )

rainfall_plot <- ggplot(
  rainfall_grid,
  aes(
    x = .data$rainfall_mm,
    y = .data$rainfall_fraction,
    colour = .data$rainfall_label
  )
) +
  geom_line(linewidth = 0.85) +
  scale_colour_manual(values = c("black", "#0072B2")) +
  scale_x_continuous(breaks = seq(0, 600, by = 100), expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(0, 1.02),
    breaks = seq(0, 1, by = 0.2),
    expand = c(0, 0)
  ) +
  labs(
    x = "Monthly rainfall (mm)",
    y = "Fraction of maximum carrying capacity",
    colour = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(colour = "black", size = 9),
    axis.title = element_text(colour = "black", size = 10),
    legend.position = "top",
    legend.text = element_text(colour = "black", size = 9)
  )

rainfall_support <- climate_support |>
  filter(
    .data$pathway == "Observed" |
      (.data$period == "2081-2100" &
        .data$pathway %in% c("SSP2-4.5", "SSP3-7.0"))
  ) |>
  mutate(
    support = factor(
      if_else(
        .data$pathway == "Observed",
        "Observed 2009-2018",
        paste(.data$pathway, .data$period)
      ),
      levels = rev(c(
        "Observed 2009-2018",
        "SSP2-4.5 2081-2100",
        "SSP3-7.0 2081-2100"
      ))
    )
  )

# Thin lines show the 1st-99th percentiles, thick lines the 10th-90th
# percentiles, and points the median monthly rainfall. This displays the model's
# climate support without adding a dense or bandwidth-dependent density plot.
rainfall_support_plot <- ggplot(
  rainfall_support,
  aes(y = .data$support, colour = .data$support)
) +
  geom_segment(
    aes(
      x = .data$rainfall_q01,
      xend = .data$rainfall_q99,
      yend = .data$support
    ),
    linewidth = 0.6
  ) +
  geom_segment(
    aes(
      x = .data$rainfall_q10,
      xend = .data$rainfall_q90,
      yend = .data$support
    ),
    linewidth = 2.3
  ) +
  geom_point(aes(x = .data$rainfall_median), size = 2.2) +
  scale_colour_manual(values = support_colours, guide = "none") +
  scale_x_continuous(
    limits = c(0, 600),
    breaks = seq(0, 600, by = 100),
    expand = c(0, 0)
  ) +
  labs(x = "Monthly rainfall (mm)", y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(colour = "black", size = 9),
    axis.title = element_text(colour = "black", size = 10),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

rainfall_mechanism_plot <- rainfall_plot / rainfall_support_plot +
  plot_layout(heights = c(2.1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 11))

ggsave(
  file.path(cfg$figures_dir, "rainfall_capacity_responses.pdf"),
  rainfall_mechanism_plot,
  width = 6.5,
  height = 5.2,
  device = grDevices::cairo_pdf
)

message("Climate-function diagnostics written to: ", output_dir)
