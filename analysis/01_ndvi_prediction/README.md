# Vegetation (NDVI)

Fits a spatial model of monthly NDVI against climate on the GIMMS grid, then
projects NDVI from 2019 to 2100 under SSP2-4.5 and SSP3-7.0 for four climate
models.

| Script, in run order | What it does |
|---|---|
| `01a_download_gimms_ndvi.R` | Downloads GIMMS NDVI (needs a NASA Earthdata account) |
| `01c_prepare_historical_feature_tables.R` | Builds the model inputs from NDVI, climate and elevation, observed and future |
| `03_fit_sdmtmb_candidates.R` | Fits and compares the candidate models |
| `03b_check_mesh_sensitivity.R` | Refits the chosen model with 10, 20 and 30 km mesh spacing |
| `04_diagnose_sdmtmb_winner.R` | Checks the chosen model on held-out data |
| `05_fit_and_project_sdmtmb.R` | Fits the final model and projects NDVI for each pathway and climate model |
| `06_plot_multimodel_projections.R` | Summarises and plots the projections |

The final model is fitted once to observed data and used for every projection,
with its year effect held at 2018. Rerunning script 05 reuses the saved fit and
any completed projections.

Settings are in `_config.R`, and pathways and climate models in
`analysis/_future_climate_config.R`. Model inputs go to `derived/` and results
to `outputs/`. The movement model reads the projections in
`outputs/scenarios/<pathway>/projections/<climate model>/monthly_ndvi.feather`.
