#!/usr/bin/env Rscript
#' Build the Serengeti District protected-area layer used in the burden maps.
#'
#' Village burden maps show cattle confined to the unprotected part of the
#' district, and a reader cannot tell why the villages cluster to the northwest
#' without seeing the protected land. This script fetches the Tanzanian
#' protected areas, keeps the three that create that pattern, clips them to the
#' district, and writes a small layer the figure scripts can read.
#'
#' Source is the World Database on Protected Areas, which requires attribution
#' wherever the layer is shown.
#'
#'   UNEP-WCMC and IUCN (2026). Protected Planet: The World Database on
#'   Protected Areas (WDPA). Cambridge, UK. www.protectedplanet.net
#'
#' The national download is about 13 MB and is not kept. Only the clipped layer
#' is written, which is a few hundred kilobytes.
#'
#' Rerun this only when refreshing the WDPA release. The monthly stamp in the
#' download URL changes, so `wdpa_release` below has to move with it.

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
})

wdpa_release <- "Aug2026"
out_file <- here::here("data", "shapefiles", "serengeti_protected_areas.gpkg")

# Ikona WMA and Loliondo GCA are left out on purpose. Wildlife management areas
# are community land where settlement is allowed, and eight modelled villages
# sit inside Ikona, so shading it would misrepresent them as uninhabited. The
# Loliondo game controlled area reaches only 13 square kilometres into the
# district and would not be visible at print size.
keep_areas <- c("Serengeti National Park", "Ikorongo G.R.", "Grumeti G.R.")

scratch <- tempfile("wdpa_")
dir.create(scratch)
zip_path <- file.path(scratch, "wdpa_tza.zip")

message("Downloading WDPA Tanzania, release ", wdpa_release, " ...")
download.file(
  sprintf(
    "https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_%s_Public_TZA_shp.zip",
    wdpa_release
  ),
  zip_path, mode = "wb", quiet = TRUE
)
unzip(zip_path, exdir = scratch)

# The country download nests three shapefile bundles, split by geometry source.
inner <- list.files(scratch, pattern = "_shp_[0-9]+\\.zip$", full.names = TRUE)
for (z in inner) unzip(z, exdir = file.path(scratch, tools::file_path_sans_ext(basename(z))))

polygon_files <- list.files(
  scratch, pattern = "polygons\\.shp$", recursive = TRUE, full.names = TRUE
)
protected <- bind_rows(lapply(polygon_files, st_read, quiet = TRUE)) |>
  st_make_valid()

# The district outline is the union of its village polygons, which is what the
# burden maps already draw.
district <- st_read(
  here::here(
    "data", "shapefiles", "resource_areas",
    "Serengeti and Ngorongoro district shape files.shp"
  ),
  quiet = TRUE
) |>
  filter(.data$District_N == "Serengeti") |>
  summarise() |>
  st_transform(32736) |>
  st_make_valid()

serengeti_protected <- protected |>
  filter(.data$NAME %in% keep_areas) |>
  st_transform(32736) |>
  st_intersection(district) |>
  group_by(.data$NAME, .data$DESIG_ENG) |>
  summarise(.groups = "drop")

serengeti_protected$area_km2 <-
  as.numeric(st_area(serengeti_protected)) / 1e6

if (nrow(serengeti_protected) != length(keep_areas)) {
  stop("Expected ", length(keep_areas), " protected areas, found ",
       nrow(serengeti_protected), ". Check the WDPA release and area names.")
}

st_write(serengeti_protected, out_file, delete_dsn = TRUE, quiet = TRUE)

district_km2 <- as.numeric(st_area(district)) / 1e6
message("Written: ", out_file)
print(as.data.frame(st_drop_geometry(serengeti_protected)), digits = 4, row.names = FALSE)
message(sprintf(
  "Protected share of the district: %.0f%% of %.0f km2",
  100 * sum(serengeti_protected$area_km2) / district_km2, district_km2
))
