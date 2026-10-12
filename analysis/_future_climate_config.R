# Shared settings for future climate analyses.
#
# Edit this file when the scenario set, GCM ensemble, or projection period
# changes. Downstream R scripts should not require scenario arguments at the
# command line.

future_climate <- list(
  scenarios = c("ssp245", "ssp370"),
  primary_scenario = "ssp245",
  scenario_labels = c(
    ssp245 = "SSP2-4.5",
    ssp370 = "SSP3-7.0",
    ssp585 = "SSP5-8.5"
  ),
  download_gcms = c(
    "GFDL_ESM4",
    "IPSL_CM6A_LR",
    "MPI_ESM1_2_LR",
    "MRI_ESM2_0",
    "UKESM1_0_LL"
  ),
  # MPI-ESM1-2-LR is retained only for quality assessment. Its SSP2-4.5 `tas`
  # channel is internally inconsistent with `tasmin` and `tasmax`, although
  # SSP3-7.0 passes the same check. A common four-GCM set is used for both SSPs
  # so pathway comparisons do not depend on different model compositions.
  reporting_gcms = c(
    "GFDL_ESM4",
    "IPSL_CM6A_LR",
    "MRI_ESM2_0",
    "UKESM1_0_LL"
  ),
  climate_years = 2015:2100,
  projection_start_year = 2019L,
  # Common 20-year reporting windows used by vegetation, movement, and
  # transmission analyses. Define them once so cross-module comparisons use
  # exactly the same years.
  analysis_windows = list(
    `2041-2060` = 2041:2060,
    `2061-2080` = 2061:2080,
    `2081-2100` = 2081:2100
  )
)

if (!future_climate$primary_scenario %in% future_climate$scenarios) {
  stop("The primary future scenario must be included in future_climate$scenarios.")
}

if (!all(future_climate$reporting_gcms %in% future_climate$download_gcms)) {
  stop("Every reporting GCM must also be listed in future_climate$download_gcms.")
}
