#!/usr/bin/env Rscript
# Build the supplementary climate table from the site-level files used by
# Module 03. The table is regenerated from analysis inputs, not maintained by
# hand, so pathway and GCM summaries stay consistent with the simulations.

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
})

repo_root <- here::here()
table_dir <- file.path(
  repo_root, "manuscript", "serengeti_rvf", "tables", "supplementary"
)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
climate_dir <- file.path(
  repo_root, "analysis", "03_rvf_transmission", "data", "processed", "climate"
)
index_path <- file.path(climate_dir, "climate_input_index.csv")

# Find each table inside this copy's climate folder rather than at the stored
# path, so the index still works after the project is moved.
climate_index <- read_csv(index_path, show_col_types = FALSE) |>
  mutate(file = file.path(
    climate_dir, basename(dirname(.data$file)), basename(.data$file)
  ))
reported_models <- c(
  "GFDL_ESM4", "IPSL_CM6A_LR", "MRI_ESM2_0", "UKESM1_0_LL"
)
model_labels <- c(
  OBSERVED_CURRENT = "Observed",
  GFDL_ESM4 = "GFDL-ESM4",
  IPSL_CM6A_LR = "IPSL-CM6A-LR",
  MRI_ESM2_0 = "MRI-ESM2-0",
  UKESM1_0_LL = "UKESM1-0-LL"
)
model_table_labels <- c(
  OBSERVED_CURRENT = "--",
  GFDL_ESM4 = "GFDL",
  IPSL_CM6A_LR = "IPSL",
  MRI_ESM2_0 = "MRI",
  UKESM1_0_LL = "UKESM"
)

# Read only the files used in the reported four-GCM ensemble. Site-months are
# pooled across village, grazing, and watering locations, matching the climate
# support passed to the transmission model.
climate <- climate_index |>
  filter(
    .data$climate_model == "OBSERVED_CURRENT" |
      .data$climate_model %in% reported_models
  ) |>
  select("file", "site_type") |>
  pmap_dfr(function(file, site_type) {
    if (!file.exists(file)) stop("Missing climate input: ", file)
    readRDS(file) |> mutate(site_type = site_type)
  }) |>
  mutate(
    pathway = if_else(.data$dataset == "observed", "Observed", .data$scenario_label),
    period = case_when(
      .data$dataset == "observed" ~ "2009--2018",
      .data$year <= 2060 ~ "2041--2060",
      .data$year <= 2080 ~ "2061--2080",
      TRUE ~ "2081--2100"
    )
  )

summary_table <- climate |>
  group_by(.data$pathway, .data$period, .data$climate_model) |>
  summarise(
    n_site_months = n(),
    median_rainfall_mm = median(.data$pr, na.rm = TRUE),
    p99_rainfall_mm = as.numeric(quantile(.data$pr, 0.99, na.rm = TRUE)),
    months_above_200_percent = 100 * mean(.data$pr > 200, na.rm = TRUE),
    months_above_400_percent = 100 * mean(.data$pr > 400, na.rm = TRUE),
    mean_temperature_c = mean(.data$tas, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    pathway = factor(.data$pathway, levels = c("Observed", "SSP2-4.5", "SSP3-7.0")),
    climate_model = factor(
      .data$climate_model,
      levels = c("OBSERVED_CURRENT", reported_models)
    ),
    gcm = unname(model_labels[as.character(.data$climate_model)]),
    gcm_table = unname(model_table_labels[as.character(.data$climate_model)])
  ) |>
  arrange(.data$pathway, .data$period, .data$climate_model)

# Suppress repeated pathway and period labels to keep the 25-row table easy to
# scan. Numeric precision follows the resolution useful for interpretation.
latex_rows <- summary_table |>
  mutate(
    pathway_text = if_else(duplicated(.data$pathway), "", as.character(.data$pathway)),
    period_group = interaction(.data$pathway, .data$period, drop = TRUE),
    period_text = if_else(duplicated(.data$period_group), "", .data$period),
    row = sprintf(
      "%s & %s & %s & %.1f & %.1f & %.1f & %.2f & %.1f \\\\",
      .data$pathway_text,
      .data$period_text,
      .data$gcm_table,
      .data$median_rainfall_mm,
      .data$p99_rainfall_mm,
      .data$months_above_200_percent,
      .data$months_above_400_percent,
      .data$mean_temperature_c
    )
  ) |>
  pull(.data$row)

latex_table <- c(
  "\\begin{table}[htbp]",
  "  \\centering",
  paste0(
    "  \\caption{\\textbf{Observed and future climate at modelled village, grazing, and watering locations.} ",
    "Values summarise all site-months used in the transmission analyses. Future rows are shown separately for each SSP, GCM, and 20-year window. ",
    "GFDL, IPSL, MRI, and UKESM denote GFDL-ESM4, IPSL-CM6A-LR, MRI-ESM2-0, and UKESM1-0-LL, respectively.}"
  ),
  "  \\label{tab:supp_climate_input_summary}",
  "  \\small",
  "  \\setlength{\\tabcolsep}{2.5pt}",
  "  \\renewcommand{\\arraystretch}{0.94}",
  "  \\begin{tabularx}{\\textwidth}{@{}ll>{\\raggedright\\arraybackslash}Xrrrrr@{}}",
  "    \\toprule",
  paste0(
    "    Pathway & Period & GCM & \\begin{tabular}[c]{@{}c@{}}Median rain\\\\(mm/month)\\end{tabular} & ",
    "\\begin{tabular}[c]{@{}c@{}}p99 rain\\\\(mm/month)\\end{tabular} & ",
    "\\begin{tabular}[c]{@{}c@{}}Months $>200$ mm\\\\(\\%)\\end{tabular} & ",
    "\\begin{tabular}[c]{@{}c@{}}Months $>400$ mm\\\\(\\%)\\end{tabular} & ",
    "\\begin{tabular}[c]{@{}c@{}}Mean temp.\\\\($^\\circ$C)\\end{tabular} \\\\"
  ),
  "    \\midrule",
  paste0("    ", latex_rows),
  "    \\bottomrule",
  "  \\end{tabularx}",
  "\\end{table}"
)

writeLines(
  latex_table,
  file.path(table_dir, "supp_climate_input_summary.tex")
)

message("Updated supplementary climate table in: ", table_dir)
