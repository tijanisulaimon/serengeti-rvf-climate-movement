# Climate data

Downloads observed CHELSA monthly climate for 1981–2018 and generates
CHELSA-CMIP6 projections under SSP2-4.5 and SSP3-7.0.

| Script, in run order | What it does |
|---|---|
| `01_download_chelsa_historical.R` | Downloads observed monthly rainfall and temperature |
| `04_download_future_scenarios.sh` | Generates the future fields by running `02_generate_chelsa_cmip6_future.py` for each pathway and climate model |
| `03_process_chelsa_cmip6_future.R` | Joins the yearly future files into one file per variable and climate model |
| `05_check_future_climate_quality.R` | Checks the future fields for gaps, temperature consistency and valid rainfall |

The future downloads need the Python packages in `requirements-chelsa.txt`.
Set them up once:

```bash
python3 -m venv analysis/00_climate_download/.venv
analysis/00_climate_download/.venv/bin/python -m pip install -r analysis/00_climate_download/requirements-chelsa.txt
```

Rerunning a script skips files that are already complete.

Pathways, climate models and years are set in `analysis/_future_climate_config.R`.
Downloads go to `data/raw/`, processed files to `data/processed/`, and the
quality checks to `outputs/climate_quality/`.

Five climate models are downloaded and four are used. MPI-ESM1-2-LR is only
checked, because its SSP2-4.5 mean temperature does not agree with its minimum
and maximum temperature.
