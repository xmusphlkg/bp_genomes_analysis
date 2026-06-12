#!/usr/bin/env Rscript
# Extended Data Fig. 11: archive identifiability and sampling DAG

get_script_dir <- function() {
  frame_files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) normalizePath(frame$ofile, winslash = "/", mustWork = FALSE) else ""
  }, character(1))
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0) return(dirname(tail(frame_files, 1)))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)))
  normalizePath(".", winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir()

source(file.path(script_dir, "..", "lib", "theme_nature.R"))

library(ggplot2)
library(patchwork)

make_nodes <- function(panel) {
  if (panel == "target") {
    data.frame(
      node = c(
        "Vaccine programme\nand product history",
        "Population immunity,\nboosters and age mix",
        "Lineage background\nand importation",
        "Antimicrobial selection",
        "Calendar time and\nrespiratory-epidemic context",
        "Circulating PRN phenotype\nprevalence",
        "Observed pertussis cases"
      ),
      x = c(1.2, 1.2, 1.2, 1.2, 1.2, 4.5, 7.6),
      y = c(6.0, 4.8, 3.6, 2.4, 1.2, 3.6, 3.6),
      fill = c("context", "context", "context", "context", "context", "target", "measured")
    )
  } else {
    data.frame(
      node = c(
        "Vaccine programme\nand product history",
        "Lineage background\nand antimicrobial selection",
        "Observed pertussis cases",
        "Surveillance intensity,\noutbreak investigation,\nand sequencing choice",
        "Public upload and\nstudy-block clustering",
        "Raw-read availability,\nassembly quality,\nand contig fragmentation",
        "PRN locus\ninterpretability",
        "Observed public-genome\nprn disruption"
      ),
      x = c(1.6, 1.8, 3.9, 4.1, 6.3, 6.3, 8.9, 11.0),
      y = c(6.0, 3.9, 5.1, 2.7, 5.1, 2.7, 1.3, 4.1),
      fill = c("context", "context", "measured", "process", "process", "process", "measured", "estimand")
    )
  }
}

make_edges <- function(panel) {
  if (panel == "target") {
    tibble::tribble(
      ~x, ~y, ~xend, ~yend, ~curvature,
      2.0, 6.0, 3.6, 3.9, -0.18,
      2.0, 4.8, 3.6, 3.75, -0.08,
      2.0, 3.6, 3.6, 3.6, 0,
      2.0, 2.4, 3.6, 3.45, 0.08,
      5.4, 3.6, 6.8, 3.6, 0,
      8.1, 3.6, 8.6, 3.6, 0
    )
  } else {
    tibble::tribble(
      ~x, ~y, ~xend, ~yend, ~curvature,
      2.4, 6.0, 3.3, 5.3, -0.12,
      2.6, 3.9, 3.5, 2.95, 0.10,
      4.8, 5.1, 5.5, 5.1, 0,
      4.9, 2.7, 5.5, 2.7, 0,
      6.9, 5.1, 10.0, 4.4, 0.08,
      6.9, 2.7, 8.1, 1.8, 0.05,
      9.4, 1.3, 10.3, 3.3, 0.12,
      4.3, 5.1, 10.1, 4.5, 0.10
    )
  }
}

palette <- c(
  context = alpha(unname(npg_colors["blue"]), 0.30),
  process = alpha(unname(npg_colors["peach"]), 0.34),
  measured = alpha(unname(npg_colors["green"]), 0.28),
  target = alpha(unname(npg_colors["red"]), 0.24),
  estimand = alpha(unname(npg_colors["teal"]), 0.46)
)

build_panel <- function(panel) {
  nodes <- make_nodes(panel)
  edges <- make_edges(panel)

  plot <- ggplot()
  for (idx in seq_len(nrow(edges))) {
    plot <- plot + geom_curve(
      data = edges[idx, , drop = FALSE],
      aes(x = x, y = y, xend = xend, yend = yend),
      curvature = edges$curvature[[idx]],
      arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
      linewidth = 0.35,
      colour = FIGURE_MUTED_TEXT
    )
  }

  plot +
    geom_label(
      data = nodes,
      aes(x = x, y = y, label = node, fill = fill),
      linewidth = 0.22,
      size = 2.3,
      label.padding = unit(0.15, "lines"),
      colour = FIGURE_INK,
      lineheight = 0.95,
      fontface = "plain"
    ) +
    scale_fill_manual(values = palette, guide = "none") +
    coord_cartesian(xlim = c(0.2, 12.1), ylim = c(0.2, 6.8), clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(6, 8, 6, 8)
    )
}

p_target <- build_panel("target")

p_archive <- build_panel("archive")

fig <- p_target / p_archive +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = FIGURE_TAG_SIZE, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.margin = margin(2, 2, 2, 2)
  )

save_ed_pdf(fig, "Extended_Data_Fig_11_Identifiability_DAG.pdf", height = 5.8)
save_ed_png(fig, "Extended_Data_Fig_11_Identifiability_DAG.png", height = 5.8)
