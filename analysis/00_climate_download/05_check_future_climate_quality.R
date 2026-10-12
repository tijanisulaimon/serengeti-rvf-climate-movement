#!/usr/bin/env Rscript

# Check whether the processed future climate fields are complete, physically
# coherent, and suitable for the two-SSP ensemble. The main diagnostic compares
# monthly mean temperature (`tas`) with the midpoint of `tasmin` and `tasmax`.
# An abrupt change in that difference indicates a channel-specific processing
# artefact rather than a plausible change in climate.

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

source(here::here("analysis", "_future_climate_config.R"))

processed_root <- here::here(
  "data", "processed", "chelsa_cmip6_future"
)
output_dir <- here::here(
  "analysis", "00_climate_download", "outputs", "climate_quality"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expected_dates <- seq(
  as.Date("2015-01-15"),
  as.Date("2100-12-15"),
  by = "month"
)

# Extract spatial monthly means and cell-level validity checks for one GCM.
summarise_gcm <- function(scenario, gcm) {
  paths <- file.path(
    processed_root,
    scenario,
    paste0(gcm, c("_tas.nc", "_tasmax.nc", "_tasmin.nc", "_pr.nc"))
  )
  names(paths) <- c("tas", "tasmax", "tasmin", "pr")

  missing_files <- paths[!file.exists(paths)]
  if (length(missing_files)) {
    stop(
      "Missing processed climate file(s) for ", scenario, " | ", gcm,
      ": ", paste(basename(missing_files), collapse = ", ")
    )
  }

  climate <- lapply(paths, terra::rast)
  layer_counts <- vapply(climate, terra::nlyr, numeric(1))
  if (any(layer_counts != length(expected_dates))) {
    stop(
      "Unexpected layer count for ", scenario, " | ", gcm, ": ",
      paste(names(layer_counts), layer_counts, sep = "=", collapse = ", ")
    )
  }

  geometry_matches <- terra::compareGeom(
    climate$tas,
    climate$tasmax,
    climate$tasmin,
    climate$pr,
    stopOnError = FALSE
  )
  if (!geometry_matches) {
    stop("Climate-variable geometry mismatch for ", scenario, " | ", gcm)
  }

  dates <- as.Date(terra::time(climate$tas))
  if (!identical(dates, expected_dates)) {
    stop("Unexpected monthly time axis for ", scenario, " | ", gcm)
  }

  # Check all months at 256 evenly spaced cells. This deterministic sample is
  # sufficient to detect the known channel discontinuity; downstream modules
  # calculate exact summaries at villages and resource locations.
  sample_cells <- unique(as.integer(round(seq(
    1,
    terra::ncell(climate$tas),
    length.out = 256
  ))))
  tas_sample <- as.matrix(terra::extract(climate$tas, sample_cells))
  tasmax_sample <- as.matrix(terra::extract(climate$tasmax, sample_cells))
  tasmin_sample <- as.matrix(terra::extract(climate$tasmin, sample_cells))
  pr_sample <- as.matrix(terra::extract(climate$pr, sample_cells))

  tas_c <- colMeans(tas_sample, na.rm = TRUE) - 273.15
  tasmax_c <- colMeans(tasmax_sample, na.rm = TRUE) - 273.15
  tasmin_c <- colMeans(tasmin_sample, na.rm = TRUE) - 273.15
  pr_mm <- colMeans(pr_sample, na.rm = TRUE)

  sampled_temperature_order_violations <- sum(
    tasmin_sample > tas_sample | tas_sample > tasmax_sample,
    na.rm = TRUE
  )
  sampled_negative_rainfall_values <- sum(pr_sample < 0, na.rm = TRUE)

  monthly <- tibble(
    scenario = scenario,
    climate_model = gcm,
    date = dates,
    year = as.integer(format(dates, "%Y")),
    month = as.integer(format(dates, "%m")),
    tas_c = tas_c,
    tasmax_c = tasmax_c,
    tasmin_c = tasmin_c,
    temperature_midpoint_c = (tasmax_c + tasmin_c) / 2,
    pr_mm_month = pr_mm
  ) |>
    mutate(
      tas_minus_midpoint_c = tas_c - temperature_midpoint_c
    )

  checks <- tibble(
    scenario = scenario,
    climate_model = gcm,
    n_months = nrow(monthly),
    start_date = min(monthly$date),
    end_date = max(monthly$date),
    sampled_cells = length(sample_cells),
    sampled_temperature_order_violations =
      sampled_temperature_order_violations,
    sampled_negative_rainfall_values = sampled_negative_rainfall_values,
    all_monthly_means_finite = all(
      is.finite(
        unlist(
          monthly[c("tas_c", "tasmax_c", "tasmin_c", "pr_mm_month")]
        )
      )
    )
  )

  list(monthly = monthly, checks = checks)
}

results <- list()
result_index <- 0L
for (scenario in future_climate$scenarios) {
  for (gcm in future_climate$download_gcms) {
    message("Checking ", scenario, " | ", gcm)
    result_index <- result_index + 1L
    results[[result_index]] <- summarise_gcm(scenario, gcm)
  }
}

monthly_summary <- bind_rows(lapply(results, `[[`, "monthly"))
base_checks <- bind_rows(lapply(results, `[[`, "checks"))

annual_summary <- monthly_summary |>
  group_by(scenario, climate_model, year) |>
  summarise(
    tas_c = mean(tas_c),
    temperature_midpoint_c = mean(temperature_midpoint_c),
    tas_minus_midpoint_c = mean(tas_minus_midpoint_c),
    pr_mm_year = sum(pr_mm_month),
    .groups = "drop"
  ) |>
  group_by(scenario, climate_model) |>
  arrange(year, .by_group = TRUE) |>
  mutate(
    annual_tas_change_c = tas_c - lag(tas_c),
    annual_midpoint_change_c =
      temperature_midpoint_c - lag(temperature_midpoint_c),
    annual_mismatch_change_c =
      tas_minus_midpoint_c - lag(tas_minus_midpoint_c)
  ) |>
  ungroup()

same_month_changes <- monthly_summary |>
  group_by(scenario, climate_model, month) |>
  arrange(year, .by_group = TRUE) |>
  mutate(
    tas_same_month_change_c = tas_c - lag(tas_c),
    mismatch_same_month_change_c =
      tas_minus_midpoint_c - lag(tas_minus_midpoint_c)
  ) |>
  ungroup() |>
  group_by(scenario, climate_model) |>
  summarise(
    max_abs_same_month_tas_change_c = max(
      abs(tas_same_month_change_c), na.rm = TRUE
    ),
    max_abs_same_month_mismatch_change_c = max(
      abs(mismatch_same_month_change_c), na.rm = TRUE
    ),
    .groups = "drop"
  )

continuity_checks <- annual_summary |>
  group_by(scenario, climate_model) |>
  summarise(
    min_annual_mean_tas_c = min(tas_c),
    max_annual_mean_tas_c = max(tas_c),
    min_annual_mean_pr_mm_month = min(pr_mm_year / 12),
    max_annual_mean_pr_mm_month = max(pr_mm_year / 12),
    max_abs_annual_tas_change_c = max(abs(annual_tas_change_c), na.rm = TRUE),
    max_abs_annual_midpoint_change_c = max(
      abs(annual_midpoint_change_c), na.rm = TRUE
    ),
    max_abs_annual_mismatch_change_c = max(
      abs(annual_mismatch_change_c), na.rm = TRUE
    ),
    tas_step_2080_2081_c =
      tas_c[year == 2081][1] - tas_c[year == 2080][1],
    midpoint_step_2080_2081_c =
      temperature_midpoint_c[year == 2081][1] -
      temperature_midpoint_c[year == 2080][1],
    mismatch_step_2080_2081_c =
      tas_minus_midpoint_c[year == 2081][1] -
      tas_minus_midpoint_c[year == 2080][1],
    .groups = "drop"
  ) |>
  left_join(
    same_month_changes,
    by = c("scenario", "climate_model")
  )

# Compare headline four-GCM summaries with and without MPI. This quantifies
# influence but is not itself an exclusion rule: a coherent GCM may reasonably
# lie outside the other models' envelope.
ensemble_influence <- annual_summary |>
  group_by(scenario, year) |>
  summarise(
    tas_mean_5 = mean(tas_c),
    pr_mean_5 = mean(pr_mm_year),
    tas_mean_4 = mean(tas_c[climate_model != "MPI_ESM1_2_LR"]),
    pr_mean_4 = mean(pr_mm_year[climate_model != "MPI_ESM1_2_LR"]),
    .groups = "drop"
  ) |>
  mutate(
    tas_5_minus_4_c = tas_mean_5 - tas_mean_4,
    pr_5_minus_4_mm = pr_mean_5 - pr_mean_4
  )

influence_checks <- ensemble_influence |>
  group_by(scenario) |>
  summarise(
    max_abs_tas_ensemble_shift_c = max(abs(tas_5_minus_4_c)),
    max_abs_pr_ensemble_shift_mm = max(abs(pr_5_minus_4_mm)),
    late_century_tas_ensemble_shift_c = mean(
      tas_5_minus_4_c[year >= 2081]
    ),
    late_century_pr_ensemble_shift_mm = mean(
      pr_5_minus_4_mm[year >= 2081]
    ),
    .groups = "drop"
  )

quality_checks <- base_checks |>
  left_join(
    continuity_checks,
    by = c("scenario", "climate_model")
  ) |>
  left_join(influence_checks, by = "scenario") |>
  mutate(
    channel_discontinuity_flag =
      max_abs_same_month_mismatch_change_c > 1,
    physical_validity_flag =
      sampled_temperature_order_violations > 0 |
      sampled_negative_rainfall_values > 0 |
      !all_monthly_means_finite,
    review_status = case_when(
      physical_validity_flag ~ "fail_physical_validity",
      channel_discontinuity_flag ~ "inspect_temperature_discontinuity",
      TRUE ~ "pass"
    )
  )

write_csv(
  monthly_summary,
  file.path(output_dir, "future_climate_monthly_summary.csv")
)
write_csv(
  annual_summary,
  file.path(output_dir, "future_climate_annual_summary.csv")
)
write_csv(
  quality_checks,
  file.path(output_dir, "future_climate_quality_checks.csv")
)
write_csv(
  ensemble_influence,
  file.path(output_dir, "future_climate_ensemble_influence.csv")
)

scenario_labels <- unname(
  future_climate$scenario_labels[future_climate$scenarios]
)
names(scenario_labels) <- future_climate$scenarios

plot_data <- annual_summary |>
  select(
    scenario,
    climate_model,
    year,
    tas_c,
    tas_minus_midpoint_c,
    pr_mm_year
  ) |>
  pivot_longer(
    cols = c(tas_c, tas_minus_midpoint_c, pr_mm_year),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    scenario = factor(
      scenario,
      levels = future_climate$scenarios,
      labels = scenario_labels
    ),
    metric = factor(
      metric,
      levels = c("tas_c", "tas_minus_midpoint_c", "pr_mm_year"),
      labels = c(
        "Annual mean temperature (°C)",
        "Mean temperature minus min/max midpoint (°C)",
        "Annual precipitation (mm)"
      )
    )
  )

quality_plot <- ggplot(
  plot_data,
  aes(x = year, y = value, colour = climate_model)
) +
  geom_line(linewidth = 0.45, alpha = 0.9) +
  facet_grid(metric ~ scenario, scales = "free_y") +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Year", y = NULL, colour = "GCM") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(colour = "black"),
    legend.text = element_text(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    strip.text = element_text(colour = "black", face = "bold")
  )

ggsave(
  file.path(output_dir, "future_climate_gcm_quality.png"),
  quality_plot,
  width = 10,
  height = 8,
  dpi = 320
)

report_lines <- c(
  "# Future Climate Quality Assessment",
  "",
  "## Dataset",
  "",
  paste0(
    "Processed CHELSA-CMIP6 monthly fields for ",
    length(future_climate$scenarios), " SSPs, ",
    length(future_climate$download_gcms), " GCMs, four variables, and ",
    length(expected_dates), " months (2015-2100). Monthly climate summaries ",
    "use 256 evenly spaced cells across the download domain."
  ),
  "",
  "## Decision Rule",
  "",
  paste0(
    "A GCM is flagged for a channel discontinuity when the year-on-year ",
    "same-calendar-month change in `tas - (tasmin + tasmax)/2` exceeds 1°C. ",
    "This tests internal channel coherence, not whether a GCM is warmer or ",
    "wetter than the ensemble."
  ),
  "",
  "## Results",
  "",
  paste0(
    "- Physical-validity failures: ",
    sum(quality_checks$physical_validity_flag), "."
  ),
  paste0(
    "- Temperature-channel discontinuity flags: ",
    sum(quality_checks$channel_discontinuity_flag), "."
  ),
  "",
  "Flagged combinations:",
  ""
)

flagged <- quality_checks |>
  filter(review_status != "pass")
if (nrow(flagged) == 0L) {
  report_lines <- c(report_lines, "- None.")
} else {
  report_lines <- c(
    report_lines,
    paste0(
      "- ", flagged$scenario, " | ", flagged$climate_model,
      ": ", flagged$review_status,
      " (maximum same-month mismatch change = ",
      round(flagged$max_abs_same_month_mismatch_change_c, 2), "°C)."
    )
  )
}

report_lines <- c(
  report_lines,
  "",
  "## Interpretation",
  "",
  paste0(
    "A model should not be excluded solely because it lies outside the other ",
    "four GCMs. Exclusion requires evidence of missing data, physical ",
    "invalidity, or a channel-specific discontinuity. Ensemble influence is ",
    "reported separately so the scientific effect of retaining MPI can be ",
    "judged transparently."
  ),
  "",
  "## MPI-ESM1-2-LR Decision",
  "",
  paste0(
    "MPI-ESM1-2-LR fails the SSP2-4.5 internal-coherence screen but passes ",
    "under SSP3-7.0. It is therefore excluded from the analysis ensemble for ",
    "both pathways. This keeps the GCM composition identical across SSPs and ",
    "avoids relying on a technical MPI-LR substitute whose SSP2-4.5 `tas` ",
    "channel is not physically valid. The downloaded and processed MPI files ",
    "are retained as auditable quality-assessment inputs."
  )
)

writeLines(
  report_lines,
  file.path(output_dir, "future_climate_quality_report.md")
)

message("Climate quality assessment written to: ", output_dir)
