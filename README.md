# Where livestock encounter vectors determines climate-driven Rift Valley fever risk

Code for the paper of the same title.

**Correspondence**: [tijani@aims.ac.za](mailto:tijani@aims.ac.za)

## Contents

The code links projected climate, vegetation, livestock movement and Rift
Valley fever virus transmission in Serengeti District, Tanzania, in four stages
that run in order.

```text
climate -> vegetation (NDVI) -> livestock movement networks -> RVF transmission
```

- `analysis/`: all code. Each stage folder has its own README with a local
  run order.
  - `00_climate_download/`: downloads and processes CHELSA historical climate and
    CHELSA-CMIP6 projections under SSP2-4.5 and SSP3-7.0.
  - `01_ndvi_prediction/`: fits monthly NDVI to climate and projects it under both
    pathways and four climate models.
  - `02_movement_networks/`: fits a resource-driven livestock movement model to
    observed village contact networks, then projects monthly networks from the
    vegetation projections.
  - `03_rvf_transmission/`: runs the vector-host transmission model on the
    projected networks and site climate.
  - **prepare_protected_areas.R**: builds the protected-area map layer used by
    the spatial figure.
  - **check_pipeline.R**: checks that the four stages share the same pathways,
    climate models and reporting windows.
  - **build_manuscript_assets.R**: builds the combined vegetation and network
    figure, gathers the other figures under their publication names, and
    regenerates the Supplement tables.
- `docs/`
  - **REPRODUCING.md**: software requirements and the full run order.
  - **DATA_ACCESS.md**: every input and how to obtain it.
- `renv.lock`: the R package versions used.

## Getting Started

1. Install the compilers and system libraries listed under Requirements in
   [`docs/REPRODUCING.md`](docs/REPRODUCING.md). Restoring the recorded package
   versions compiles many packages from source. Then restore the R packages
   from `renv.lock`, using R 4.5.3:

   ```r
   install.packages("renv")
   renv::activate()
   # restart R in the project folder, then:
   renv::restore(exclude = "arrow")
   install.packages("arrow")
   ```

2. Obtain the inputs listed in [`docs/DATA_ACCESS.md`](docs/DATA_ACCESS.md). None
   of the input data are stored in this repository.
3. Run the stages in the order given in
   [`docs/REPRODUCING.md`](docs/REPRODUCING.md).

## Licence

MIT. See [`LICENSE`](LICENSE).
