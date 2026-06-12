# Nature Portfolio figure theme and utilities
# Unified visual style for all main and Extended Data figures.
# v3 2026-04-27: final-size typography, Nature display dimensions,
#   restrained semantic palettes, axis-only framing, vector-first exports.

library(ggplot2)
library(cowplot)
library(patchwork)
library(scales)
library(RColorBrewer)

# ---------------------------------------------------------------------------
#  1. Colour semantics: one meaning, one colour family
# ---------------------------------------------------------------------------

# Shared figure palettes. Set PERTUSSIS_FIGURE_PALETTE=earth to preview the
# alternate scheme; publication is the default for manuscript rendering.
figure_palette_options <- list(
  publication = list(
    discrete = c(
      "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
      "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF"
    ),
    roles = c(
      red = "#E64B35FF",
      teal = "#4DBBD5FF",
      green = "#00A087FF",
      blue = "#3C5488FF",
      orange = "#F39B7FFF",
      peach = "#FFB84DFF",
      purple = "#8491B4FF",
      yellow = "#FFB84DFF",
      brown = "#7E6148FF",
      tan = "#B09C85FF",
      grey = "#8491B4FF",
      light_grey = "#EDEAE6FF",
      mid_grey = "#B09C85FF",
      dark_grey = "#7E6148FF",
      black = "#373634FF"
    )
  ),
  earth = list(
    discrete = c(
      "#785D37FF", "#62BBA5FF", "#FFB84DFF", "#AAA488FF", "#B2432FFF",
      "#3A6589FF", "#9B5672FF", "#908150FF", "#373634FF"
    ),
    roles = c(
      red = "#B2432FFF",
      teal = "#62BBA5FF",
      green = "#AAA488FF",
      blue = "#3A6589FF",
      orange = "#FFB84DFF",
      peach = "#908150FF",
      purple = "#9B5672FF",
      yellow = "#FFB84DFF",
      brown = "#785D37FF",
      tan = "#908150FF",
      grey = "#AAA488FF",
      light_grey = "#E8E4DAFF",
      mid_grey = "#908150FF",
      dark_grey = "#785D37FF",
      black = "#373634FF"
    )
  )
)

FIGURE_PALETTE_NAME <- Sys.getenv("PERTUSSIS_FIGURE_PALETTE", unset = "publication")
if (!FIGURE_PALETTE_NAME %in% names(figure_palette_options)) {
  warning(sprintf("Unknown PERTUSSIS_FIGURE_PALETTE='%s'; using publication.", FIGURE_PALETTE_NAME))
  FIGURE_PALETTE_NAME <- "publication"
}
figure_palette <- figure_palette_options[[FIGURE_PALETTE_NAME]]
figure_discrete <- figure_palette$discrete
figure_discrete_at <- function(i) {
  unname(figure_discrete[((i - 1) %% length(figure_discrete)) + 1])
}

# Backwards-compatible semantic palette used across the scripts.
npg_colors <- figure_palette$roles

figure_warm_seq <- c(
  "#FFE1BDFF", "#FFD6A3FF", "#FFC47DFF", "#FFB04FFF",
  "#FF9B21FF", "#F28500FF", "#DE7A00FF"
)
figure_teal_seq <- c(
  "#CCE5E5FF", "#99CCCCFF", "#66B2B2FF", "#329999FF",
  "#008080FF", "#006666FF", "#004C4CFF"
)

FIGURE_BASE_SIZE <- 6.5
FIGURE_TAG_SIZE <- 8
FIGURE_LEGEND_TEXT_SIZE <- 5
FIGURE_LEGEND_TITLE_SIZE <- 5.5
FIGURE_ANNOT_SIZE <- 1.85
FIGURE_TEXT_COLOUR <- unname(npg_colors["black"])
FIGURE_INK <- unname(npg_colors["black"])
FIGURE_MUTED_TEXT <- unname(npg_colors["dark_grey"])
FIGURE_DARK_GREY <- unname(npg_colors["dark_grey"])
FIGURE_MID_GREY <- unname(npg_colors["mid_grey"])
FIGURE_GREY <- unname(npg_colors["grey"])
FIGURE_LIGHT_GREY <- unname(npg_colors["light_grey"])
FIGURE_RULE_COLOUR <- "#D8D2C7FF"
FIGURE_PANEL_FILL <- "#F5F3F0FF"
FIGURE_MAP_FILL <- "#F0F7F6FF"

# PRN status — used across the active main-figure contract
prn_colors <- c(
  "Intact"       = unname(npg_colors["blue"]),
  "Disrupted"    = unname(npg_colors["red"]),
  "Insufficient" = unname(npg_colors["grey"]),
  "Insufficient data" = unname(npg_colors["grey"]),
  "Reference" = unname(npg_colors["black"])
)

prn_binary_colors <- c(
  "Intact" = unname(npg_colors["blue"]),
  "Disrupted" = unname(npg_colors["red"]),
  "Non-interpretable / uncertain" = unname(npg_colors["grey"])
)

# Disruption mechanism — used in Figures 2, 3, and 5
mechanism_colors <- c(
  "IS481"          = unname(npg_colors["red"]),
  "IS481 insertion" = unname(npg_colors["red"]),
  "IS481 1043 bp" = unname(npg_colors["red"]),
  "Inversion"      = unname(npg_colors["blue"]),
  "Inversion / rearrangement" = unname(npg_colors["blue"]),
  "Rearrangement" = unname(npg_colors["blue"]),
  "Other"          = unname(npg_colors["green"]),
  "Other disruptions" = unname(npg_colors["green"]),
  "Mixed/Unknown"  = unname(npg_colors["mid_grey"]),
  "Insufficient data" = unname(npg_colors["grey"])
)

architecture_colors <- c(
  "IS481 gap1043" = unname(npg_colors["red"]),
  "IS481 1,043-bp" = unname(npg_colors["red"]),
  "Rearr. cov58" = unname(npg_colors["blue"]),
  "Rearr. cov91" = unname(npg_colors["teal"]),
  "Rearrangement family" = unname(npg_colors["blue"]),
  "gap1042" = unname(npg_colors["peach"]),
  "gap1045" = unname(npg_colors["purple"]),
  "gap54" = unname(npg_colors["dark_grey"]),
  "gap204" = unname(npg_colors["green"]),
  "IS481 gap1041" = unname(npg_colors["orange"]),
  "Other gap1041" = unname(npg_colors["mid_grey"]),
  "Other gap1044" = unname(npg_colors["grey"]),
  "Other gap1040" = unname(npg_colors["tan"]),
  "Rearr. cov94" = unname(npg_colors["teal"]),
  "Other insertion-like" = unname(npg_colors["green"]),
  "Insufficient" = unname(npg_colors["light_grey"]),
  "Other" = unname(npg_colors["grey"])
)

country_role_colors <- c(
  "USA" = figure_discrete_at(8),
  "NZL" = figure_discrete_at(1),
  "JPN" = figure_discrete_at(4),
  "AUS" = figure_discrete_at(3),
  "GBR" = figure_discrete_at(6),
  "CHN" = figure_discrete_at(5),
  "FRA" = figure_discrete_at(2),
  "BRA" = figure_discrete_at(9),
  "CZE" = figure_discrete_at(10),
  "FIN" = figure_discrete_at(7)
)

# Evidence / support tier — used in Figures 2, 3
support_colors <- c(
  "strong"         = unname(npg_colors["green"]),
  "moderate"       = unname(npg_colors["peach"]),
  "weak"           = unname(npg_colors["grey"]),
  "limited"        = unname(npg_colors["teal"]),
  "not_applicable" = unname(npg_colors["light_grey"]),
  "none"           = unname(npg_colors["light_grey"])
)

evidence_colors <- c(
  "Read-backed"       = unname(npg_colors["green"]),
  "Long-read / hybrid" = unname(npg_colors["blue"]),
  "Published overlap"  = unname(npg_colors["peach"]),
  "Assembly only"      = unname(npg_colors["light_grey"])
)

# Vaccine-programme exposure — supplementary ecology diagnostics
exposure_colors <- c(
  "aP (PRN+)" = unname(npg_colors["orange"]),
  "aP (PRN-)" = unname(npg_colors["blue"]),
  "wP / pre-aP" = unname(npg_colors["dark_grey"]),
  "Mixed / transition" = unname(npg_colors["purple"]),
  "Uncertain" = unname(npg_colors["grey"])
)

bmj_program_colors <- c(
  "wp_only_or_pre_ap"                = figure_discrete_at(9),
  "routine_ap_prn_negative"          = figure_discrete_at(4),
  "routine_ap_mixed"                 = figure_discrete_at(3),
  "routine_ap_prn_positive"          = figure_discrete_at(5),
  "transition_mixed_within_period"   = figure_discrete_at(6),
  "routine_ap_unknown"               = unname(npg_colors["grey"])
)

programme_epoch_colors <- c(
  "wP only" = figure_discrete_at(9),
  "aP + PRN" = figure_discrete_at(5),
  "aP - PRN" = figure_discrete_at(4),
  "Mixed PRN" = figure_discrete_at(3),
  "Transition" = figure_discrete_at(6),
  "Uncertain" = unname(npg_colors["grey"])
)

direction_colors <- c(
  "Upward" = unname(npg_colors["red"]),
  "Downward" = unname(npg_colors["blue"]),
  "No change" = unname(npg_colors["grey"]),
  "Uncertain" = unname(npg_colors["grey"])
)

identifiability_colors <- c(
  "Identified" = unname(npg_colors["green"]),
  "Bounded" = unname(npg_colors["peach"]),
  "Needs new data" = unname(npg_colors["grey"]),
  "Not identified" = unname(npg_colors["grey"])
)

bmj_status_colors <- c(
  "eligible_wp_only_or_pre_ap"                = unname(bmj_program_colors["wp_only_or_pre_ap"]),
  "eligible_routine_ap_prn_negative"          = unname(bmj_program_colors["routine_ap_prn_negative"]),
  "eligible_routine_ap_mixed"                 = unname(bmj_program_colors["routine_ap_mixed"]),
  "eligible_routine_ap_prn_positive"          = unname(bmj_program_colors["routine_ap_prn_positive"]),
  "eligible_transition_mixed_within_period"   = unname(bmj_program_colors["transition_mixed_within_period"]),
  "missing_ipw_response"                      = FIGURE_PANEL_FILL,
  "below_interpretable_threshold_5"           = FIGURE_GREY,
  "missing_reported_cases"                    = FIGURE_MID_GREY
)

# Sequential ramps
teal_seq <- colorRampPalette(figure_teal_seq)(100)
orange_seq <- colorRampPalette(figure_warm_seq)(100)
blue_seq <- teal_seq
red_seq <- orange_seq
grey_seq <- colorRampPalette(c("#FFFFFF", FIGURE_PANEL_FILL, FIGURE_GREY, FIGURE_DARK_GREY))(100)

# ---------------------------------------------------------------------------
#  2. Nature display dimensions (mm -> inches internally)
# ---------------------------------------------------------------------------

NC_SINGLE_COL  <- 89  / 25.4   # 3.504 in
NC_1_5_COL     <- 120 / 25.4   # 4.724 in
NC_WIDE_COL    <- 136 / 25.4   # 5.354 in
NC_DOUBLE_COL  <- 183 / 25.4   # 7.205 in
NC_MAX_HEIGHT  <- 150 / 25.4   # 6.693 in; leaves caption/legend space
NC_PAGE_HEIGHT <- 247 / 25.4   # 9.724 in; absolute page depth reference

# ---------------------------------------------------------------------------
#  3. theme_nature — main panel theme (axis-only framing, 7 pt base)
# ---------------------------------------------------------------------------

theme_nature <- function(base_size = FIGURE_BASE_SIZE, base_family = "") {
  half_line <- base_size / 2
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Typography
      text             = element_text(colour = FIGURE_TEXT_COLOUR, size = base_size),
      plot.title       = element_text(face = "bold", size = rel(1.05), hjust = 0,
                                      margin = margin(b = half_line)),
      plot.subtitle    = element_text(size = rel(0.9), colour = FIGURE_MUTED_TEXT,
                                      margin = margin(b = half_line)),
      plot.caption     = element_text(size = rel(0.8), colour = FIGURE_MUTED_TEXT, hjust = 0),
      plot.tag         = element_text(face = "bold", size = FIGURE_TAG_SIZE),
      plot.tag.position = c(0, 1),
      plot.margin      = margin(2.5, 2.5, 2.5, 2.5),

      # Axes — clean, thin lines at edges only
      axis.line        = element_line(colour = FIGURE_INK, linewidth = 0.28),
      axis.text        = element_text(size = rel(0.86), colour = FIGURE_TEXT_COLOUR),
      axis.title       = element_text(size = rel(0.92), colour = FIGURE_TEXT_COLOUR),
      axis.ticks       = element_line(colour = FIGURE_INK, linewidth = 0.28),
      axis.ticks.length = unit(1.5, "pt"),

      # Panel — no border, no major grid
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.spacing    = unit(5, "pt"),

      # Strips (facets)
      strip.background = element_rect(fill = FIGURE_PANEL_FILL, colour = NA),
      strip.text       = element_text(face = "bold", size = rel(0.84),
                                      margin = margin(2, 0, 2, 0)),

      # Legend — compact, bottom
      legend.background = element_rect(fill = NA, colour = NA),
      legend.key        = element_rect(fill = NA, colour = NA),
      legend.key.size   = unit(6, "pt"),
      legend.title      = element_text(face = "bold", size = rel(0.82)),
      legend.text       = element_text(size = rel(0.76)),
      legend.position   = "bottom",
      legend.box        = "horizontal",
      legend.spacing.x  = unit(2.5, "pt"),
      legend.margin     = margin(0, 0, 0, 0)
    )
}

# ---------------------------------------------------------------------------
#  3b. theme_nature_tree — for phylogenetic panels (no axes, no grid)
# ---------------------------------------------------------------------------

theme_nature_tree <- function(base_size = FIGURE_BASE_SIZE, base_family = "") {
  theme_nature(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line        = element_blank(),
      axis.text        = element_blank(),
      axis.title       = element_blank(),
      axis.ticks       = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border     = element_blank(),
      legend.position  = "none"
    )
}

# ---------------------------------------------------------------------------
#  3c. theme_nature_map — for cartographic panels
# ---------------------------------------------------------------------------

theme_nature_map <- function(base_size = FIGURE_BASE_SIZE, base_family = "") {
  theme_nature(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line  = element_blank(),
      axis.text  = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(colour = FIGURE_RULE_COLOUR, fill = NA, linewidth = 0.28)
    )
}

# ---------------------------------------------------------------------------
#  3d. Multipanel finishing themes
# ---------------------------------------------------------------------------

theme_nature_multipanel <- function(base_size = FIGURE_BASE_SIZE) {
  theme(
    plot.margin = margin(2, 2, 2, 2),
    plot.tag = element_text(face = "bold", size = FIGURE_TAG_SIZE, colour = FIGURE_TEXT_COLOUR),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(3, "pt"),
    legend.spacing.y = unit(1.5, "pt"),
    legend.key.size = unit(6, "pt"),
    legend.title = element_text(face = "bold", size = rel(0.76)),
    legend.text = element_text(size = rel(0.74)),
    strip.text = element_text(face = "bold", size = rel(0.80), margin = margin(2, 0, 2, 0))
  )
}

theme_nature_matrix <- function(base_size = FIGURE_BASE_SIZE, base_family = "") {
  theme_nature(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = rel(0.82), colour = FIGURE_TEXT_COLOUR),
      axis.title = element_blank(),
      legend.position = "bottom"
    )
}

theme_nature_extended_data <- function(base_size = FIGURE_BASE_SIZE) {
  theme_nature_multipanel(base_size = base_size) +
    theme(
      plot.subtitle = element_blank(),
      legend.box = "horizontal",
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.x = unit(3, "pt"),
      legend.spacing.y = unit(1.5, "pt"),
      legend.key.size = unit(6, "pt")
    )
}

figure_annotation_theme <- function(base_size = FIGURE_BASE_SIZE) {
  theme(
    plot.tag = element_text(face = "bold", size = FIGURE_TAG_SIZE, colour = FIGURE_TEXT_COLOUR),
    plot.margin = margin(2, 2, 2, 2)
  )
}

nature_tag_levels <- function(n = 8) {
  list(letters[seq_len(n)])
}

# ---------------------------------------------------------------------------
#  4. Utility functions
# ---------------------------------------------------------------------------

format_pval <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

get_support_script_dir <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      if (!is.null(frame$ofile)) {
        return(normalizePath(frame$ofile, winslash = "/", mustWork = FALSE))
      }
      ""
    },
    character(1)
  )
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0) {
    return(dirname(tail(frame_files, 1)))
  }
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)))
  }
  normalizePath(".", winslash = "/", mustWork = FALSE)
}

locate_repo_root <- function(start_dir = getwd()) {
  current <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  repeat {
    if (dir.exists(file.path(current, ".git")) &&
        dir.exists(file.path(current, "manuscript")) &&
        dir.exists(file.path(current, "workflow"))) {
      return(current)
    }
    parent <- normalizePath(file.path(current, ".."), winslash = "/", mustWork = FALSE)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }
  stop(sprintf("Could not locate repository root from: %s", start_dir))
}

get_figure_root <- function() {
  env_root <- Sys.getenv("PERTUSSIS_FIGURE_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = FALSE))
  }
  repo_env <- Sys.getenv("PERTUSSIS_REPO_ROOT", unset = "")
  repo_root <- if (nzchar(repo_env)) {
    normalizePath(repo_env, winslash = "/", mustWork = FALSE)
  } else {
    locate_repo_root(get_support_script_dir())
  }
  normalizePath(file.path(repo_root, "manuscript", "figures"), winslash = "/", mustWork = FALSE)
}

resolve_output_dir <- function(subdir = "main") {
  out_dir <- normalizePath(file.path(get_figure_root(), "outputs", subdir), winslash = "/", mustWork = FALSE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_dir
}

# ---------------------------------------------------------------------------
#  5. NC-compliant save helpers
# ---------------------------------------------------------------------------

save_figure <- function(plot, filename,
                        width  = NC_DOUBLE_COL,
                        height = NC_MAX_HEIGHT,
                        units  = "in",
                        dpi    = 300,
                        subdir = "main") {
  if (identical(units, "in")) {
    if (width > NC_DOUBLE_COL) {
      message(sprintf(
        "Width %.2f in exceeds Nature double-column width %.2f in; saving at requested size for legibility. Redesign the panel count before final submission.",
        width, NC_DOUBLE_COL
      ))
    }
    if (height > NC_MAX_HEIGHT) {
      message(sprintf(
        "Height %.2f in exceeds the 170 mm target %.2f in; saving at requested size for legibility. Repackage oversized audit panels as supplementary figures or supplementary data before final submission.",
        height, NC_MAX_HEIGHT
      ))
    }
  }
  out_dir <- resolve_output_dir(subdir)
  path <- file.path(out_dir, filename)
  ggsave(path, plot, width = width, height = height, units = units, dpi = dpi, bg = "white")
  message(sprintf("Saved: %s (%.2f x %.2f in, %d dpi)", path, width, height, dpi))
  invisible(path)
}

save_nc_pdf <- function(plot, filename,
                        width  = NC_DOUBLE_COL,
                        height = NC_MAX_HEIGHT,
                        subdir = "main") {
  save_figure(plot, filename, width = width, height = height, dpi = 600, subdir = subdir)
}

save_nc_png <- function(plot, filename,
                        width  = NC_DOUBLE_COL,
                        height = NC_MAX_HEIGHT,
                        subdir = "main") {
  save_figure(plot, filename, width = width, height = height, dpi = 300, subdir = subdir)
}

save_ed_pdf <- function(plot, filename,
                        width  = NC_DOUBLE_COL,
                        height = NC_MAX_HEIGHT) {
  save_figure(plot, filename, width = width, height = height, dpi = 600, subdir = "extended_data")
}

save_ed_png <- function(plot, filename,
                        width  = NC_DOUBLE_COL,
                        height = NC_MAX_HEIGHT) {
  save_figure(plot, filename, width = width, height = height, dpi = 300, subdir = "extended_data")
}

save_bmj_figure <- function(plot, filename, width = 7.2, height = 9, units = "in", res = 300) {
  save_figure(plot, filename, width = width, height = height, dpi = res, subdir = "bmj")
}

# ---------------------------------------------------------------------------
#  6. Data helpers (kept for backward compat)
# ---------------------------------------------------------------------------

validate_columns <- function(df, required_cols, label = "data") {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf("[%s] Missing required columns: %s", label, paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}

data_quality <- function(df, label = "data") {
  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("Rows: %d | Columns: %d\n", nrow(df), ncol(df)))
  cat(sprintf("Missing values: %d (%.1f%%)\n",
              sum(is.na(df)),
              100 * sum(is.na(df)) / (nrow(df) * ncol(df))))
  invisible(df)
}
