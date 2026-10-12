# Rift Valley fever transmission

Runs a cattle-mosquito transmission model across villages, grazing areas and
water points, using the simulated movement networks and site climate.

| Script, in run order | What it does |
|---|---|
| `01_prepare_inputs.R` | Joins each simulated network to the site climate, one input file per scenario year and network |
| `02_calibrate_vector_capacity.R` | Calibrates mosquito capacity for each exposure profile to 4–10% cattle seroprevalence |
| `03_run_observed_baseline.R` | Runs the model under observed climate and movement |
| `04_run_mobility_scenarios.R` | Varies the fraction of cattle that use shared grazing and water sites |
| `05_summarise_mobility_scenarios.R` | Summarises script 04 |
| `06_run_future_scenarios.R` | Runs the model under both pathways and four climate models |
| `07_summarise_future_scenarios.R` | Compares future with current burden over matched 40-year periods |
| `08_run_climate_movement_attribution.R` | Separates the effects of climate change and network change |
| `09_summarise_spatial_seasonal_patterns.R` | Summarises where and when transmission occurs |
| `10_build_manuscript_outputs.R` | Builds the remaining figures and the headline tables |

`tests/quick_model_check.R` runs one short simulation, in about 5 seconds, to
check the model works after scripts 01 and 02.

Run settings are in the `rvf_execution` block of `R/config.R`: cores, network
replicates, calibration profiles, and a short test mode (`smoke`) for scripts
04, 06 and 08. Model parameters are in `R/model_parameters.R`. Results go to
`outputs/`.
