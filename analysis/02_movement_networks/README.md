# Livestock movement networks

Fits a model of where herds graze and drink to the village contact networks
recorded by participatory mapping, then simulates monthly networks under
observed and projected vegetation for the transmission model.

The recorded data are village-to-village contacts. The model therefore
simulates which grazing areas and water points each village uses, and is
judged on the contacts those choices imply.

| Script, in run order | What it does |
|---|---|
| `01_prepare_inputs.R` | Prepares site distances and NDVI at each village and site (needs the field data) |
| `02_define_targets.R` | Sets resource availability, the monthly conditions (routine, transition, drought) and the observed network summaries to fit (needs the field data) |
| `03_fit_model.R` | Fits the model, separately for grazing and watering |
| `04_validate_model.R` | Compares new simulations with the observed networks (needs the field data) |
| `05_project_networks.R` | Simulates 100 networks for each observed and projected year. Rerun it to resume |
| `06_summarise_projections.R` | Summarises the projected networks |

Settings are in `_config.R`, and shared functions in `R/simulator.R` and
`R/utils.R`. Outputs go to `outputs/`. Keep `outputs/projections/job_summaries/`,
because script 06 reads it and script 05 needs it to resume.
