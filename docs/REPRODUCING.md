# Reproducing the analysis

This repository reproduces the figures and tables of *Interacting effects of
climate on vectors and livestock movements drive Rift Valley fever infection
risk*, for Serengeti District, Tanzania. The analysis runs in four stages, each
in its own folder under `analysis/`. A final step gathers the figures under
their publication names and regenerates the Supplement tables.

Run every command below from the project root.

## Requirements

### System

- **R 4.5.3.** The code was tested with R 4.5.3 on macOS 27 (Apple M4, 10 cores,
  32 GB RAM). Linux and Windows have not been tested.
- **About 60 GB of disk** for downloaded and prepared inputs.
- **A multi-core machine** for the movement and transmission simulations. Set
  the core counts to no more than your machine has: `n_cores` in
  `analysis/02_movement_networks/_config.R` (6 as supplied) and `max_cores` in
  the `rvf_execution` block of `analysis/03_rvf_transmission/R/config.R` (8 as
  supplied). The same block sets the number of network replicates.
- **Python 3**, for the future-climate downloads only. Its packages are listed in
  `analysis/00_climate_download/requirements-chelsa.txt`.
- **A NASA Earthdata account** for the GIMMS NDVI download. Supply it through
  `EARTHDATA_BEARER_TOKEN`, `EARTHDATA_NETRC`, or `EARTHDATA_USERNAME` with
  `EARTHDATA_PASSWORD`. Credentials are read from the environment so they never
  enter the repository.

### Compilers and system libraries

`renv.lock` records the package versions that produced the published results.
CRAN distributes ready-built packages only for current versions, so restoring
the recorded versions compiles many packages from source (89 of the 147 in
October 2026). Install these tools before restoring.

**macOS**

- Xcode Command Line Tools: `xcode-select --install`.
- CRAN's GNU Fortran compiler, which installs to `/opt/gfortran`, from
  <https://mac.r-project.org/tools/>. `deSolve`, `igraph`, `mvtnorm` and `TMB`
  need it. R does not find a Homebrew `gfortran`.
- The OpenMP runtime for Apple clang, from <https://mac.r-project.org/openmp/>.
  `data.table` needs it.
- System libraries from Homebrew:

  ```bash
  brew install pkgconf udunits gdal geos proj harfbuzz fribidi freetype libpng jpeg-turbo libtiff webp
  ```

**Debian or Ubuntu** (not tested)

```bash
sudo apt-get install build-essential gfortran cmake pkg-config libudunits2-dev libgdal-dev libgeos-dev libproj-dev libsqlite3-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libjpeg-dev libtiff-dev libwebp-dev libxml2-dev libcurl4-openssl-dev libssl-dev libfontconfig1-dev
```

**Windows** (not tested): Rtools45 provides the compilers and most of the
libraries.

### R packages

Start R in the project folder and run:

```r
install.packages("renv")
renv::activate()
```

Restart R in the same folder, then run:

```r
renv::restore(exclude = "arrow")
install.packages("arrow")
```

`renv` installs the recorded versions into a library inside the project folder,
so your own R library is left unchanged. `arrow` is the one exception. Its
recorded version (21.0.0.1) looks for a prebuilt Arrow C++ library that was
never published, then falls back to a long source build. The current CRAN
release installs without compiling. The code uses `arrow` only to read and
write feather tables, and `arrow` 25.0.1 reads all 17 of the project's feather
files identically to 21.0.0.1 (checked October 2026).

## Inputs that no script downloads

Some inputs must be put in place by hand before the stage that reads them.
[`DATA_ACCESS.md`](DATA_ACCESS.md) lists each file, where it goes and how to
obtain it. Folder locations are set in `_project_paths.R`.

The files from the participatory field study cannot be shared publicly. Without
them, the climate and vegetation stages run in full. The protected-area step
needs the district boundaries, and the movement and transmission stages need the
movement data, from which their inputs are built.

## Run order

### 1. Climate

Historical climate, then the future scenarios. The shell script launches the
Python download worker for each pathway and climate model.

```bash
Rscript analysis/00_climate_download/01_download_chelsa_historical.R
bash analysis/00_climate_download/04_download_future_scenarios.sh
Rscript analysis/00_climate_download/03_process_chelsa_cmip6_future.R
Rscript analysis/00_climate_download/05_check_future_climate_quality.R
```

### 2. Vegetation

Fits the NDVI candidates, checks the selected model, then projects it under both
pathways and all four climate models.

```bash
Rscript analysis/01_ndvi_prediction/01a_download_gimms_ndvi.R
Rscript analysis/01_ndvi_prediction/01c_prepare_historical_feature_tables.R
Rscript analysis/01_ndvi_prediction/03_fit_sdmtmb_candidates.R
Rscript analysis/01_ndvi_prediction/03b_check_mesh_sensitivity.R
Rscript analysis/01_ndvi_prediction/04_diagnose_sdmtmb_winner.R
Rscript analysis/01_ndvi_prediction/05_fit_and_project_sdmtmb.R
Rscript analysis/01_ndvi_prediction/06_plot_multimodel_projections.R
```

### 3. Movement networks

Fits the resource-selection model to the observed networks, validates it, then
projects monthly networks under observed and future vegetation.

```bash
Rscript analysis/02_movement_networks/01_prepare_inputs.R
Rscript analysis/02_movement_networks/02_define_targets.R
Rscript analysis/02_movement_networks/03_fit_model.R
Rscript analysis/02_movement_networks/04_validate_model.R
Rscript analysis/02_movement_networks/05_project_networks.R
Rscript analysis/02_movement_networks/06_summarise_projections.R
```

### 4. Protected areas

Downloads the World Database on Protected Areas for Tanzania and clips it to the
district. The spatial figure in the next stage reads the result. The release
month is pinned in the script, because the download address changes with each
release.

```bash
Rscript analysis/prepare_protected_areas.R
```

### 5. Transmission

Prepare the inputs and calibrate first. Scripts 02, 03, 04, 06 and 08 run
simulations. Script 10 builds the remaining figures and the headline tables.

```bash
Rscript analysis/03_rvf_transmission/01_prepare_inputs.R
Rscript analysis/03_rvf_transmission/02_calibrate_vector_capacity.R
Rscript analysis/03_rvf_transmission/03_run_observed_baseline.R
Rscript analysis/03_rvf_transmission/04_run_mobility_scenarios.R
Rscript analysis/03_rvf_transmission/05_summarise_mobility_scenarios.R
Rscript analysis/03_rvf_transmission/06_run_future_scenarios.R
Rscript analysis/03_rvf_transmission/07_summarise_future_scenarios.R
Rscript analysis/03_rvf_transmission/08_run_climate_movement_attribution.R
Rscript analysis/03_rvf_transmission/09_summarise_spatial_seasonal_patterns.R
Rscript analysis/03_rvf_transmission/10_build_manuscript_outputs.R
```

To check an installation before the long runs, run one short simulation after
scripts 01 and 02. It takes about 5 seconds on the test machine:

```bash
Rscript analysis/03_rvf_transmission/tests/quick_model_check.R
```

### 6. Check the stages agree, then gather the figures and tables

The check reads the completed outputs of all four stages and confirms they share
the same pathways, climate models and reporting windows. The second script builds
the combined vegetation and network figure, copies every other figure under its
publication name, and regenerates five of the Supplement tables. It writes them
into `manuscript/serengeti_rvf/`, creating the folder.

```bash
Rscript analysis/check_pipeline.R
Rscript analysis/build_manuscript_assets.R
```

Two Supplement tables, the vegetation mesh sensitivity and the movement
calibration weights, are maintained by hand in the manuscript rather than written
by this script.

## Where the outputs go

Each stage writes to an `outputs/` folder inside its own directory. Step 6 then
gathers the figures and tables used by the paper into `manuscript/serengeti_rvf/`.
The analysis outputs stay the source, and that folder holds copies.
