#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${CHELSA_PYTHON:-${SCRIPT_DIR}/.venv/bin/python}"
DOWNLOADER="${SCRIPT_DIR}/02_generate_chelsa_cmip6_future.py"
MAX_JOBS="${CHELSA_MAX_JOBS:-2}"

read -r -a SCENARIOS <<< "${CHELSA_SCENARIOS:-ssp245 ssp370}"
read -r -a GCMS <<< "${CHELSA_GCMS:-GFDL-ESM4 IPSL-CM6A-LR MPI-ESM1-2-LR MRI-ESM2-0 UKESM1-0-LL}"

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "CHELSA Python environment not found: ${PYTHON_BIN}" >&2
  echo "Create it using the setup commands in analysis/00_climate_download/README.md." >&2
  exit 2
fi

if ! [[ "${MAX_JOBS}" =~ ^[1-6]$ ]]; then
  echo "CHELSA_MAX_JOBS must be an integer from 1 to 6." >&2
  exit 2
fi

printf "Scenarios: %s\n" "${SCENARIOS[*]}"
printf "GCMs: %s\n" "${GCMS[*]}"
printf "Concurrent downloads: %s\n" "${MAX_JOBS}"

for scenario in "${SCENARIOS[@]}"; do
  for gcm in "${GCMS[@]}"; do
    printf "%s\t%s\n" "${gcm}" "${scenario}"
  done
done |
  xargs -n 2 -P "${MAX_JOBS}" "${PYTHON_BIN}" "${DOWNLOADER}"
