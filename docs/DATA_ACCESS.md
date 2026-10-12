# Data access

The analysis reads three kinds of input. Some are downloaded by the scripts,
some are public but must be downloaded by hand, and some come from the
participatory field study. None are stored in this repository.

## Downloaded by the scripts

| Source | Script | Notes |
|---|---|---|
| CHELSA monthly climate, 1981--2018 | `analysis/00_climate_download/01_download_chelsa_historical.R` | |
| CHELSA-CMIP6 projections, SSP2-4.5 and SSP3-7.0 | `analysis/00_climate_download/04_download_future_scenarios.sh` | Needs the Python environment in `requirements-chelsa.txt` |
| GIMMS NDVI | `analysis/01_ndvi_prediction/01a_download_gimms_ndvi.R` | Needs a NASA Earthdata account |
| World Database on Protected Areas, Tanzania | `analysis/prepare_protected_areas.R` | The release month is pinned in the script |

## Public, downloaded by hand

| File | Place it in | Source |
|---|---|---|
| `gadm41_TZA_2_pk.rds` | `data/shapefiles/` | GADM 4.1 boundaries for Tanzania at level 2, downloaded with `geodata::gadm("tanzania", level = 2)`. The package saves it inside a `gadm/` subfolder, so move it up into `data/shapefiles/` |
| `serengeti_ele.tif` | `data/elevation/` | SRTM GL3 elevation at 3 arc-seconds, about 90 m, from OpenTopography. Downloaded with `elevatr::get_elev_raster()` using `src = "gl3"`, `z = 12` and `prj = "EPSG:4326"`, for the Serengeti District polygon from the GADM file above. OpenTopography requires a free API key |

The `geodata` and `elevatr` packages are needed only for these two one-off
downloads, so they are not pinned in `renv.lock`.

## From the participatory field study

These come from participatory mapping with livestock keepers in Serengeti
District (Ekwem et al. 2023). The movement files record monthly village-to-village
contacts for grazing and for watering separately. They do not identify which
grazing area or water point produced each contact.

| File | Place it in | Content |
|---|---|---|
| `Local_movt_Gza_tj.csv` | `data/observed_movements/` | Monthly grazing contacts, 2017 |
| `Local_movt_wp_tj.csv` | `data/observed_movements/` | Monthly watering contacts, 2017 |
| `serengeti village centre xy.csv` | `data/shapefiles/resource_areas/` | Village centre coordinates |
| `All grazing areas_centroid.shp` | `data/shapefiles/resource_areas/` | Grazing-area centroids |
| `All water points.shp` | `data/shapefiles/resource_areas/` | Water-point locations |
| `Serengeti and Ngorongoro district shape files.shp` | `data/shapefiles/resource_areas/` | Village, ward and district boundaries from the 2012 census, supplied with the field data. The National Bureau of Statistics publishes district and ward boundaries but not village boundaries. Outside movement step 01, only the Serengeti District outline is used |
| `restricted_site_patterns.txt` | `data/shapefiles/resource_areas/` | Site-name patterns marking sites whose land designation restricts livestock access, one per line |

These files contain sensitive information on the locations of communal
livestock resources and cannot be shared publicly due to government
regulations. Summary statistics of the networks have been published in
Ekwem et al. (2023).
