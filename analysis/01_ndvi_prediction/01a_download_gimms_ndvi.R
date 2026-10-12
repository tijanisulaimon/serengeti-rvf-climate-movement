#!/usr/bin/env Rscript
#' ####################################################################### #
#' 01a — Download GIMMS-3G+ NDVI
#'
#' Query NASA CMR for the GIMMS-3G+ collection (ORNL DAAC, 1982-2022)
#' and download all NetCDF granules into data/raw/GIMMS_3Gplus/.
#' Run once; set GIMMS_OVERWRITE=true to re-download.
#'
#' Auth: ~/.netrc | EARTHDATA_NETRC | EARTHDATA_BEARER_TOKEN |
#'       EARTHDATA_USERNAME + EARTHDATA_PASSWORD
#'
#' INPUTS:  NASA CMR API
#' OUTPUTS: data/raw/GIMMS_3Gplus/*.nc
#' ####################################################################### #

suppressPackageStartupMessages({
  library(jsonlite)
  library(curl)
  library(readr)
  library(dplyr)
  library(purrr)
})

source(here::here("analysis", "01_ndvi_prediction", "_config.R"))

cfg <- ndvi_config
overwrite <- tolower(Sys.getenv("GIMMS_OVERWRITE", unset = "false")) %in% c("1", "true", "yes", "y")
dry_run <- tolower(Sys.getenv("GIMMS_DRY_RUN", unset = "false")) %in% c("1", "true", "yes", "y")
max_files <- suppressWarnings(as.integer(Sys.getenv("GIMMS_MAX_FILES", unset = "")))

cmr_url <- paste0(
  cfg$gimms_cmr_granules_url,
  "?collection_concept_id=",
  cfg$gimms_collection_concept_id,
  "&page_size=2000"
)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

extract_granule_rows <- function(entry) {
  links <- entry$links %||% list()
  rows <- purrr::map_dfr(links, function(x) {
    href <- x$href %||% ""
    if (!grepl("\\.(nc|nc4)(\\?|$)", href, ignore.case = TRUE)) {
      return(NULL)
    }

    tibble::tibble(
      url = href,
      rel = x$rel %||% NA_character_,
      link_title = x$title %||% NA_character_
    )
  })

  if (!nrow(rows)) {
    return(NULL)
  }

  rows |>
    mutate(
      title = entry$title %||% NA_character_,
      updated = entry$updated %||% NA_character_,
      link_type = dplyr::case_when(
        grepl("^https://data\\.ornldaac\\.earthdata\\.nasa\\.gov/", url) ~ "ornl_download",
        grepl("^https://opendap\\.earthdata\\.nasa\\.gov/", url) ~ "opendap",
        grepl("^s3://", url) ~ "s3",
        TRUE ~ "other"
      ),
      link_rank = dplyr::case_when(
        link_type == "ornl_download" ~ 1L,
        link_type == "opendap" ~ 2L,
        link_type == "s3" ~ 3L,
        TRUE ~ 4L
      )
    ) |>
    select(title, updated, url, rel, link_title, link_type, link_rank)
}

payload <- jsonlite::fromJSON(cmr_url, simplifyVector = FALSE)
entries <- payload$feed$entry %||% list()

granules <- purrr::map(entries, extract_granule_rows) |>
  purrr::compact() |>
  bind_rows() |>
  arrange(link_rank, url) |>
  mutate(
    filename = basename(sub("\\?.*$", "", url)),
    dest = file.path(cfg$gimms_raw_dir, filename)
  ) |>
  filter(link_type == "ornl_download") |>
  distinct(filename, .keep_all = TRUE)

if (!nrow(granules)) {
  stop("No NetCDF download URLs found for GIMMS-3G+ via CMR.")
}

if (!is.na(max_files) && max_files > 0 && max_files < nrow(granules)) {
  granules <- granules[seq_len(max_files), , drop = FALSE]
}

dir.create(cfg$gimms_raw_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(granules, cfg$gimms_granule_list)

message("GIMMS granule list written to: ", cfg$gimms_granule_list)
message("Granules listed: ", nrow(granules))

if (dry_run) {
  message("Dry run requested; no files downloaded.")
  quit(save = "no", status = 0)
}

token <- Sys.getenv("EARTHDATA_BEARER_TOKEN", unset = "")
username <- Sys.getenv("EARTHDATA_USERNAME", unset = "")
password <- Sys.getenv("EARTHDATA_PASSWORD", unset = "")
netrc_env <- Sys.getenv("EARTHDATA_NETRC", unset = "")
netrc_path <- if (nzchar(netrc_env)) path.expand(netrc_env) else path.expand("~/.netrc")
cookie_path <- file.path(tempdir(), "earthdata_cookies.txt")

auth_mode <- if (nzchar(token)) {
  "bearer_token"
} else if (nzchar(username) && nzchar(password)) {
  "username_password"
} else if (file.exists(netrc_path)) {
  "netrc"
} else {
  "missing"
}

if (identical(auth_mode, "missing")) {
  stop(
    "No Earthdata authentication found.\n",
    "Provide one of:\n",
    "  1. ~/.netrc with: machine urs.earthdata.nasa.gov login <user> password <pass>\n",
    "  2. EARTHDATA_NETRC=/path/to/.netrc\n",
    "  3. EARTHDATA_BEARER_TOKEN\n",
    "  4. EARTHDATA_USERNAME and EARTHDATA_PASSWORD\n"
  )
}

message("Earthdata auth mode: ", auth_mode)
if (identical(auth_mode, "netrc")) {
  message("Using .netrc file: ", netrc_path)
}

build_handle <- function() {
  h <- curl::new_handle(
    followlocation = TRUE,
    cookiejar = cookie_path,
    cookiefile = cookie_path,
    useragent = "resourceDrivenMoveModel/ndvi_prediction GIMMS downloader"
  )

  if (identical(auth_mode, "bearer_token")) {
    curl::handle_setopt(
      h,
      httpheader = c(Authorization = paste("Bearer", token)),
      unrestricted_auth = 1L
    )
  } else if (identical(auth_mode, "username_password")) {
    curl::handle_setopt(
      h,
      username = username,
      password = password,
      unrestricted_auth = 1L
    )
  } else if (identical(auth_mode, "netrc")) {
    curl::handle_setopt(
      h,
      netrc = 1L,
      netrc_file = netrc_path
    )
  }

  h
}

for (i in seq_len(nrow(granules))) {
  dest <- granules$dest[[i]]
  url <- granules$url[[i]]

  if (file.exists(dest) && !overwrite) {
    message("Skipping existing file: ", basename(dest))
    next
  }

  message("Downloading GIMMS granule ", i, "/", nrow(granules), ": ", basename(dest))
  tryCatch(
    curl::curl_download(url, destfile = dest, handle = build_handle(), quiet = FALSE),
    error = function(e) {
      stop(
        "Failed while downloading ", basename(dest), ".\n",
        "Current auth mode: ", auth_mode, "\n",
        "Original error: ", conditionMessage(e), "\n",
        "If you see HTTP 401 from urs.earthdata.nasa.gov, try ~/.netrc first."
      )
    }
  )
}

message("GIMMS raw download complete.")
