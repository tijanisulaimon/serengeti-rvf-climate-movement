#!/usr/bin/env python3
"""
02 - Generate CHELSA-CMIP6 future climate as monthly time series for the
Serengeti study domain. One CHELSA-derived NetCDF per month per variable per
GCM is produced for every year in 2015-2100 (matching the existing CHELSA-CMIP6
archive structure of monthly time series).

Follows Example 2 of the chelsa_cmip6 README: loop year-by-year, calling
chelsa_cmip6() with fefps = YYYY-01-01 and fefpe = YYYY-12-31. Handles the
UKESM1-0-LL 360-day-calendar case (last-of-year = -12-30 instead of -12-31).

Usage (one GCM and scenario per process):
  python 02_generate_chelsa_cmip6_future.py GFDL-ESM4 ssp245
  python 02_generate_chelsa_cmip6_future.py UKESM1-0-LL ssp370
  ...

Idempotent: a year is skipped if its expected NetCDFs already exist.
Transient download failures are retried three times with increasing delays.
Raw files are isolated under:
  data/raw/chelsa_cmip6_future/{scenario}/{GCM}/{year}/

Citation: Karger, D.N. et al. (2023). chelsa-cmip6 1.0: a python package
to create high resolution bioclimatic variables based on CHELSA ver. 2.1
and CMIP6 data. Ecography. doi:10.1111/ecog.06535
"""

import os
import sys
import time
import logging
from pathlib import Path

from chelsa_cmip6.GetClim import chelsa_cmip6

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RAW_ROOT = PROJECT_ROOT / "data" / "raw" / "chelsa_cmip6_future"
LOG_ROOT = PROJECT_ROOT / "analysis" / "00_climate_download" / "logs"
RAW_ROOT.mkdir(parents=True, exist_ok=True)
LOG_ROOT.mkdir(parents=True, exist_ok=True)

# Bbox: existing Serengeti crop (34.22-35.27, -2.50, -1.42) + 0.2 deg buffer.
XMIN, XMAX = 34.02, 35.47
YMIN, YMAX = -2.70, -1.22

# CHELSA V2.1 standard reference period.
REFPS = "1981-01-15"
REFPE = "2010-12-15"

# Years to generate. Existing archive covered 2015-2100 (1032 months); we
# match that. If interrupted, restart picks up from cached year-files.
START_YEAR = int(os.getenv("CHELSA_START_YEAR", "2015"))
END_YEAR = int(os.getenv("CHELSA_END_YEAR", "2100"))

# Four ISIMIP3b core GCMs plus MPI-ESM1-2-LR as a technical substitute for
# the core MPI-ESM1-2-HR model.
# Tuple: (institution_id, calendar_last_day, member_id, use_esgf)
# Pangeo-hosting verified by direct query of cmip6-zarr-consolidated-stores.csv.
# Constraint: chelsa_cmip6 reuses one institution_id for BOTH the historical
# reference query (CMIP) and the future scenario query (ScenarioMIP).
# Models therefore need a single institution that hosts both branches under
# r1i1p1f1 Amon.
#   - MPI-ESM1-2-HR is split-hosted (CMIP under MPI-M, ScenarioMIP under DKRZ),
#     so chelsa_cmip6 can NOT use it. Falling back to MPI-ESM1-2-LR, which is
#     single-institution (MPI-M) and is the chelsa_cmip6 README's example.
#   - UKESM1-0-LL uses member_id=r1i1p1f2 and a 360-day calendar.
GCM_REGISTRY = {
    "GFDL-ESM4":      ("NOAA-GFDL", 31, "r1i1p1f1", False),
    "IPSL-CM6A-LR":   ("IPSL",      31, "r1i1p1f1", False),
    "MPI-ESM1-2-LR":  ("MPI-M",     31, "r1i1p1f1", False),  # LR fallback (HR is split-hosted)
    "MRI-ESM2-0":     ("MRI",       31, "r1i1p1f1", False),
    "UKESM1-0-LL":    ("MOHC",      30, "r1i1p1f2", False),  # 360-day cal + f2 member
}

SUPPORTED_EXPERIMENTS = {"ssp126", "ssp245", "ssp370", "ssp585"}
TABLE_ID = "Amon"
ACTIVITY_ID = "ScenarioMIP"
MAX_ATTEMPTS = 3
RETRY_DELAYS_SECONDS = (30, 90)

# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------

def setup_logger(source_id: str, experiment_id: str) -> logging.Logger:
    log_dir = LOG_ROOT / experiment_id
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"chelsa_cmip6_{source_id}.log"
    logger = logging.getLogger(f"chelsa.{experiment_id}.{source_id}")
    logger.setLevel(logging.INFO)
    if logger.handlers:
        return logger
    fmt = logging.Formatter(
        f"%(asctime)s [{source_id} %(levelname)s] %(message)s"
    )
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    fh = logging.FileHandler(log_path, mode="a")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    return logger


def year_output_dir(source_id: str, experiment_id: str, year: int) -> Path:
    d = RAW_ROOT / experiment_id / source_id / f"{year}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def expected_year_complete(out_dir: Path, experiment_id: str, year: int) -> bool:
    """Check that all four core future-variable files exist for this year.

    Each output directory also contains historical reference-period files, so
    variable and scenario names alone are not sufficient evidence that the
    requested future year completed.
    """
    need = ("tas", "tasmax", "tasmin", "pr")
    year_span = f"_{year}-01-01_{year}-12-"
    nc_files = [p.name for p in out_dir.glob("*.nc") if p.stat().st_size > 0]
    for var in need:
        if not any(
            f"_{var}_" in f
            and f"_{experiment_id}_" in f
            and year_span in f
            for f in nc_files
        ):
            return False
    return True


def run_one_year(source_id: str, institution_id: str, member_id: str,
                 experiment_id: str, year: int, last_day: int, use_esgf: bool,
                 logger: logging.Logger) -> bool:
    out_dir = year_output_dir(source_id, experiment_id, year)
    if expected_year_complete(out_dir, experiment_id, year):
        logger.info("skip %s %d (already complete: %d NetCDFs)",
                    source_id, year, len(list(out_dir.glob("*.nc"))))
        return True

    fefps = f"{year}-01-01"
    fefpe = f"{year}-12-{last_day:02d}"
    logger.info("run  %s %d  (%s..%s, member=%s, esgf=%s) -> %s",
                source_id, year, fefps, fefpe, member_id, use_esgf, out_dir)
    t0 = time.time()
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            chelsa_cmip6(
                activity_id=ACTIVITY_ID, table_id=TABLE_ID,
                experiment_id=experiment_id,
                institution_id=institution_id,
                source_id=source_id,
                member_id=member_id,
                refps=REFPS, refpe=REFPE,
                fefps=fefps, fefpe=fefpe,
                xmin=XMIN, xmax=XMAX,
                ymin=YMIN, ymax=YMAX,
                output=str(out_dir) + os.sep,
                use_esgf=use_esgf,
            )
            break
        except Exception as exc:
            if attempt == MAX_ATTEMPTS:
                logger.error(
                    "fail %s %d after %d attempts: %s",
                    source_id, year, MAX_ATTEMPTS, exc
                )
                return False
            delay = RETRY_DELAYS_SECONDS[attempt - 1]
            logger.warning(
                "retry %s %d after attempt %d/%d failed: %s; waiting %ds",
                source_id, year, attempt, MAX_ATTEMPTS, exc, delay
            )
            time.sleep(delay)

    # CHELSA can return without an exception while one variable is still
    # missing. Treat that case as a failed year so a later resume retries it.
    if not expected_year_complete(out_dir, experiment_id, year):
        logger.error(
            "incomplete %s %d after download: missing one or more core "
            "variables (tas, tasmax, tasmin, pr)",
            source_id,
            year,
        )
        return False

    logger.info("done %s %d in %.1f min",
                source_id, year, (time.time() - t0) / 60.0)
    return True


def main(argv: list[str]) -> int:
    if len(argv) not in (2, 3) or argv[1] not in GCM_REGISTRY:
        print(__doc__)
        print("Available GCMs:", ", ".join(sorted(GCM_REGISTRY.keys())))
        print("Available scenarios:", ", ".join(sorted(SUPPORTED_EXPERIMENTS)))
        return 2
    source_id = argv[1]
    experiment_id = (
        argv[2] if len(argv) == 3
        else os.getenv("CHELSA_SCENARIO", "ssp245")
    ).lower()
    if experiment_id not in SUPPORTED_EXPERIMENTS:
        print(
            f"Unsupported scenario {experiment_id!r}; choose one of "
            f"{', '.join(sorted(SUPPORTED_EXPERIMENTS))}",
            file=sys.stderr,
        )
        return 2
    if START_YEAR > END_YEAR:
        print("CHELSA_START_YEAR must not exceed CHELSA_END_YEAR", file=sys.stderr)
        return 2

    institution_id, last_day, member_id, use_esgf = GCM_REGISTRY[source_id]
    log = setup_logger(source_id, experiment_id)
    log.info("Starting monthly-timeseries generation: %s, %s (%s), %d-%d",
             experiment_id, source_id, member_id, START_YEAR, END_YEAR)
    log.info("Bbox: %.2f..%.2f lon, %.2f..%.2f lat",
             XMIN, XMAX, YMIN, YMAX)

    overall_t0 = time.time()
    n_done = n_failed = 0
    for year in range(START_YEAR, END_YEAR + 1):
        ok = run_one_year(source_id, institution_id, member_id, experiment_id, year,
                          last_day, use_esgf, log)
        if ok:
            n_done += 1
        else:
            n_failed += 1
            # Persist past failures rather than abort: a single year failing
            # is usually a transient Pangeo timeout that the next run-through
            # of the script will retry.
    elapsed_h = (time.time() - overall_t0) / 3600.0
    log.info("Finished %s. years_ok=%d, years_failed=%d, elapsed=%.2f h",
             source_id, n_done, n_failed, elapsed_h)
    return 0 if n_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
