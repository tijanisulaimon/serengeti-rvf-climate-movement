#' Shared plotting conventions for Module 03 figures.
#'
#' One place for the publication style so every manuscript figure looks the
#' same and figure captions live in LaTeX rather than baked into the plot.
#' Figure scripts should source this and build on `theme_rvf()` + the SSP/GCM
#' scales below, rather than re-styling ad hoc.
#'
#'   source(here::here("analysis", "03_rvf_transmission", "R", "plot_helpers.R"))

suppressPackageStartupMessages(library(ggplot2))

# --- Two-SSP palette + labels ----------------------------------------------
# The canonical scheme used across every Module 03 figure. SSP2-4.5 is the
# primary (central) pathway; SSP3-7.0 is the high pathway. Colour-blind safe.
rvf_ssp_levels  <- c("SSP2-4.5", "SSP3-7.0")
rvf_ssp_colours <- c("SSP2-4.5" = "#2c7fb8", "SSP3-7.0" = "#d95f0e")

# Map internal scenario codes (ssp245 / ssp370) to display labels, returning a
# factor with a fixed order so panels and legends are consistent everywhere.
rvf_ssp_label <- function(x) {
  x <- as.character(x)
  x[x %in% c("ssp245", "SSP245")] <- "SSP2-4.5"
  x[x %in% c("ssp370", "SSP370")] <- "SSP3-7.0"
  factor(x, levels = rvf_ssp_levels)
}

# --- Publication theme -----------------------------------------------------
# No titles/subtitles: interpretation belongs in the LaTeX caption. A classic
# white background keeps axes and uncertainty marks legible at journal size.
theme_rvf <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      axis.line         = element_line(colour = "black", linewidth = 0.45),
      axis.ticks        = element_line(colour = "black", linewidth = 0.45),
      axis.text         = element_text(colour = "black"),
      axis.title        = element_text(colour = "black"),
      strip.background  = element_blank(),
      strip.text        = element_text(colour = "black", face = "bold"),
      legend.position   = "bottom",
      legend.title      = element_blank(),
      legend.text       = element_text(colour = "black")
    )
}

# --- Consistent scales -----------------------------------------------------
scale_colour_ssp <- function(...) scale_colour_manual(values = rvf_ssp_colours, drop = FALSE, ...)
scale_fill_ssp   <- function(...) scale_fill_manual(values = rvf_ssp_colours, drop = FALSE, ...)

# Four-GCM reporting ensemble (MPI excluded upstream). Qualitative, cb-friendly.
rvf_gcm_colours <- c(
  GFDL_ESM4    = "#e41a1c", IPSL_CM6A_LR = "#377eb8",
  MRI_ESM2_0   = "#4daf4a", UKESM1_0_LL  = "#984ea3"
)
rvf_gcm_labels <- c(
  GFDL_ESM4 = "GFDL-ESM4",
  IPSL_CM6A_LR = "IPSL-CM6A-LR",
  MRI_ESM2_0 = "MRI-ESM2-0",
  UKESM1_0_LL = "UKESM1-0-LL"
)
scale_colour_gcm <- function(...) {
  scale_colour_manual(
    values = rvf_gcm_colours,
    breaks = names(rvf_gcm_labels),
    labels = unname(rvf_gcm_labels),
    ...
  )
}

# --- Default two-SSP layout ------------------------------------------------
# SSP on facet columns (optionally a second variable on rows). Expects an
# `ssp_label` column produced with rvf_ssp_label().
facet_ssp <- function(rows = NULL, scales = "free_y", ...) {
  if (is.null(rows)) {
    facet_wrap(ggplot2::vars(.data$ssp_label), nrow = 1, scales = scales, ...)
  } else {
    facet_grid(rows = rows, cols = ggplot2::vars(.data$ssp_label), scales = scales, ...)
  }
}

# --- Consistent save -------------------------------------------------------
# Single width/dpi so figures scale uniformly; height set per figure. Manuscript
# figures go to outputs/figures/; diagnostics belong in their own script dirs.
rvf_ggsave <- function(filename, plot, width = 8, height = 5, dpi = 300) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename, plot, width = width, height = height, dpi = dpi)
}
