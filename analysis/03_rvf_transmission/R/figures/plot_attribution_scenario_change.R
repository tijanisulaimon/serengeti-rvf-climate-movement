#!/usr/bin/env Rscript
# Replot the climate-movement attribution as scenario changes relative to the
# current-climate, current-movement cell.
#
# The published version split the decomposition into direct-climate, movement
# and interaction components across two panels on different scales, which made
# the reader reassemble the total mentally. This version plots what each
# scenario actually produces, all relative to cell A:
#
#   climate_contribution  = B - A   future climate, current movement
#   movement_contribution = C - A   current climate, future movement
#   total_change          = D - A   both future
#
# Those are the quantities already stored by run_attribution_profile.R, so no
# new computation is involved. The interaction term is the residual
# (D - B - C + A) and is reported in the caption as a numerical check rather
# than drawn, because it has no scenario of its own.
#
# Plotting climate-only against both-future on one axis is the point: the two
# series nearly coincide, which is the finding.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(stringr)
})

source(here::here("analysis", "03_rvf_transmission", "R", "config.R"))
cfg <- rvf_config

results_dir <- cfg$results_dir
# The window summary collapses across GCMs, so its min and max equal the median.
# The aggregate file keeps one row per GCM and stores the raw 2x2 cells, which
# is what the GCM range has to be computed from.
aggregate_files <- list.files(
  results_dir,
  pattern = "climate_movement_attribution_decomposition_aggregate\\.csv$",
  recursive = TRUE, full.names = TRUE
)
stopifnot(length(aggregate_files) > 0)

per_gcm <- bind_rows(lapply(aggregate_files, read_csv, show_col_types = FALSE))

profiles_kept <- c("grazing_floodplain_dominant", "peridomestic_dominant")
profile_labels <- c(
  grazing_floodplain_dominant = "Shared grazing/floodplain",
  peridomestic_dominant       = "Village-associated"
)
window_labels <- c(
  `2050` = "2041-2060", `2080` = "2061-2080", `2100` = "2081-2100"
)

# Each scenario contrast measured against cell A, then summarised across GCMs.
scenario_change <- per_gcm |>
  filter(.data$site_profile %in% profiles_kept) |>
  # Burden is reported per 1,000 cattle throughout the manuscript, so the
  # district-wide counts are divided by the 470 thousand modelled cattle.
  transmute(
    scenario_label, future_window, site_profile, climate_model,
    `Both future`          = (.data$D_mid - .data$A_mid) / 470,
    `Future climate only`  = (.data$B_mid - .data$A_mid) / 470,
    `Future movement only` = (.data$C_mid - .data$A_mid) / 470
  ) |>
  pivot_longer(
    c(`Both future`, `Future climate only`, `Future movement only`),
    names_to = "series", values_to = "change"
  ) |>
  group_by(scenario_label, future_window, site_profile, series) |>
  summarise(
    change_mid = median(.data$change),
    change_lo  = min(.data$change),
    change_hi  = max(.data$change),
    .groups = "drop"
  ) |>
  mutate(
    profile_panel = factor(
      unname(profile_labels[.data$site_profile]),
      levels = unname(profile_labels)
    ),
    window_label = factor(
      unname(window_labels[as.character(.data$future_window)]),
      levels = unname(window_labels)
    ),
    series = factor(
      series,
      levels = c("Both future", "Future climate only", "Future movement only")
    )
  )

# Okabe-Ito, matching the framework figure. Both-future and climate-only sit
# almost on top of each other, so they need clearly different shapes as well.
series_colours <- c(
  "Both future"          = "#000000",
  "Future climate only"  = "#D55E00",
  "Future movement only" = "#0072B2"
)
series_shapes <- c(
  "Both future" = 16, "Future climate only" = 17, "Future movement only" = 15
)

p <- ggplot(
  scenario_change,
  aes(window_label, change_mid, colour = series, shape = series, group = series)
) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_linerange(
    aes(ymin = change_lo, ymax = change_hi),
    position = position_dodge(width = 0.5), linewidth = 0.7
  ) +
  geom_point(position = position_dodge(width = 0.5), size = 2.6) +
  facet_grid(profile_panel ~ scenario_label) +
  scale_colour_manual(values = series_colours, name = NULL) +
  scale_shape_manual(values = series_shapes, name = NULL) +
  labs(
    x = "Future climate window",
    y = "Change in annual infections per 1,000 cattle\nrelative to current climate and movement"
  ) +
  theme_classic(base_size = 12.5) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black", size = 10.5),
    axis.title = element_text(colour = "black", size = 12),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 11, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 11)
  )

out <- file.path(cfg$figures_dir, "08_attribution_scenario_change.png")
ggsave(out, p, width = 9.6, height = 6.4, dpi = 400)
message("Written: ", out)

# Numerical check reported in the caption rather than drawn.
check <- per_gcm |>
  filter(.data$site_profile %in% profiles_kept) |>
  summarise(
    climate_share = sprintf("%.1f-%.1f%%", min(climate_share_of_total_mid) * 100,
                            max(climate_share_of_total_mid) * 100),
    movement = sprintf("%.0f to %.0f", min(movement_contribution_mid),
                       max(movement_contribution_mid)),
    interaction = sprintf("%.0f to %.0f", min(interaction_contribution_mid),
                          max(interaction_contribution_mid))
  )
print(as.data.frame(check))
