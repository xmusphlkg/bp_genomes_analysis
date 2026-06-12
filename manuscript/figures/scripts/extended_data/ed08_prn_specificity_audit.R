#!/usr/bin/env Rscript
# Extended Data Fig. 8: PRN-locus structural specificity audit

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
source(file.path(script_dir, "..", "lib", "data_utils.R"))
source(file.path(script_dir, "..", "lib", "manuscript_rebuild_helpers.R"))

library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(scales)

controls <- safe_load(file.path(FIGURE_DATA_DIR, "selected_country", "prn_specificity_negative_control.tsv"), "PRN specificity controls") %>%
  mutate(
    signal_positive_fraction_among_interpretable = tidyr::replace_na(as.numeric(signal_positive_fraction_among_interpretable), 0),
    signal_rate_per_1000_interpretable_genomes = tidyr::replace_na(as.numeric(signal_rate_per_1000_interpretable_genomes), 0),
    n_signal_positive = as.numeric(n_signal_positive),
    locus_label = if_else(locus == "prn", "PRN", locus_label),
    locus_group = case_when(
      locus == "prn" ~ "PRN",
      locus_category %in% c("structure_matched_autotransporter", "pertactin_homologous_autotransporter") ~ "Structure-matched",
      TRUE ~ "Other comparator"
    )
  )

specificity_colors <- c(
  "PRN" = unname(npg_colors["red"]),
  "Structure-matched" = unname(npg_colors["blue"]),
  "Other comparator" = FIGURE_GREY
)

global <- controls %>%
  filter(table_scope == "global_overlap_frame") %>%
  mutate(locus_label = factor(locus_label, levels = locus_label[order(signal_positive_fraction_among_interpretable)]))

pA <- ggplot(global, aes(signal_positive_fraction_among_interpretable, locus_label, fill = locus_group)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.18) +
  scale_x_continuous(labels = percent) +
  scale_fill_manual(values = specificity_colors, name = NULL) +
  labs(x = "Signal-positive fraction", y = NULL) +
  theme_nature() +
  theme(legend.position = "bottom")

pB <- ggplot(global, aes(signal_rate_per_1000_interpretable_genomes, locus_label, fill = locus_group)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.18) +
  scale_x_continuous(trans = "pseudo_log", breaks = c(0, 10, 100, 500), labels = comma) +
  scale_fill_manual(values = specificity_colors, guide = "none") +
  labs(x = "Signal-positive calls per 1,000 interpretable genomes", y = NULL) +
  theme_nature()

selected_loci <- c("PRN", "FHA/FhaB", "TcfA", "BrkA", "SphB1", "Phg", "Vag8", "BapC", "Fim2", "Fim3")
selected_countries <- c("USA", "NZL", "JPN", "AUS")

selected <- controls %>%
  filter(table_scope == "selected_country_epoch_overlap_frame", !is.na(country_iso3), country_iso3 %in% selected_countries) %>%
  group_by(country_iso3, locus_label) %>%
  summarise(
    signal_positive_fraction_among_interpretable = max(tidyr::replace_na(signal_positive_fraction_among_interpretable, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    country_iso3 = factor(country_iso3, levels = selected_countries),
    locus_label = factor(locus_label, levels = rev(selected_loci))
  ) %>%
  filter(!is.na(locus_label)) %>%
  tidyr::complete(
    country_iso3,
    locus_label,
    fill = list(signal_positive_fraction_among_interpretable = 0)
  )

pC <- ggplot(selected, aes(country_iso3, locus_label, fill = signal_positive_fraction_among_interpretable)) +
  geom_tile(colour = "white", linewidth = 0.28) +
  scale_fill_gradientn(
    colours = red_seq,
    labels = percent_format(accuracy = 1),
    name = "Signal fraction",
    guide = guide_colourbar(
      direction = "vertical",
      barheight = grid::unit(46, "pt"),
      barwidth = grid::unit(4, "pt")
    )
  ) +
  theme_nature_matrix() +
  theme(
    axis.text.x = element_text(face = "bold"),
    legend.position = "right",
    legend.key.height = grid::unit(46, "pt"),
    legend.title = element_text(size = 5.8, face = "bold"),
    legend.text = element_text(size = 5.4)
  )

ed8 <- (pA | pB) / free(pC) +
  plot_layout(heights = c(0.95, 1.05)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = FIGURE_TAG_SIZE, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot",
    plot.margin = margin(3, 3, 3, 3)
  )

save_ed_pdf(ed8, "Extended_Data_Fig_08_PRN_Specificity_Audit.pdf", height = NC_MAX_HEIGHT)
save_ed_png(ed8, "Extended_Data_Fig_08_PRN_Specificity_Audit.png", height = NC_MAX_HEIGHT)
