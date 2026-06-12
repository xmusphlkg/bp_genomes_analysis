#!/usr/bin/env Rscript
# Extended Data Fig. 1: sample-frame reconciliation and selection audit

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
library(tidyr)

frame <- safe_load(file.path(BASE_DIR, "manuscript", "supplementary", "Supplementary_Table_33_Sample_Frame_Reconciliation.tsv"), "sample-frame reconciliation") %>%
  mutate(
    n_total = as.numeric(n_total),
    n_prn_interpretable = as.numeric(n_prn_interpretable),
    n_prn_disrupted = as.numeric(n_prn_disrupted)
  )

score <- safe_load(data_paths$rr_selection_scorecard, "selection scorecard") %>%
  mutate(
    stage1_primary_default = as_logical_flag(stage1_primary_default),
    stage2_triangulated_default = as_logical_flag(stage2_triangulated_default),
    n_prn_interpretable = as.numeric(n_prn_interpretable),
    interpretability_fraction = as.numeric(interpretability_fraction)
  )

country_levels <- ordered_countries(unique(score$country_iso3))

pA <- frame %>%
  filter(sample_frame %in% c("retained_country_total", "broad_descriptive_prevalence_panel", "selected_country_epoch_analysis")) %>%
  group_by(country_iso3, sample_frame) %>%
  summarise(n_total = sum(n_total, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    country_iso3 = factor(country_iso3, levels = rev(country_levels)),
    sample_frame = recode(sample_frame,
      retained_country_total = "Retained country total",
      broad_descriptive_prevalence_panel = "Descriptive panel",
      selected_country_epoch_analysis = "Selected-country epochs"
    )
  ) %>%
  ggplot(aes(n_total, country_iso3, fill = sample_frame)) +
  geom_col(position = position_dodge2(width = 0.7, preserve = "single"), width = 0.62, colour = "white", linewidth = 0.15) +
  scale_x_continuous(labels = comma) +
  scale_fill_manual(values = c("Retained country total" = FIGURE_DARK_GREY, "Descriptive panel" = npg_colors["blue"], "Selected-country epochs" = npg_colors["green"]), name = NULL) +
  labs(x = "Genomes", y = NULL) +
  theme_nature() +
  theme(legend.position = "bottom")

pB <- score %>%
  filter(country_iso3 %in% country_levels) %>%
  transmute(
    country_iso3 = factor(country_iso3, levels = rev(country_levels)),
    `Stage 1 primary` = if_else(stage1_primary_default, "yes", "no"),
    `Stage 2 triangulated` = if_else(stage2_triangulated_default, "yes", "no"),
    `Primary only` = if_else(selection_state == "primary_only", "yes", "no")
  ) %>%
  pivot_longer(-country_iso3, names_to = "status", values_to = "value") %>%
  ggplot(aes(status, country_iso3)) +
  geom_tile(aes(fill = value), colour = "white", linewidth = 0.32) +
  scale_fill_manual(values = c("yes" = npg_colors["green"], "no" = FIGURE_LIGHT_GREY), guide = "none") +
  theme_nature_matrix() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), axis.text.y = element_text(face = "bold"))

pC <- score %>%
  filter(country_iso3 %in% country_levels) %>%
  mutate(country_iso3 = factor(country_iso3, levels = rev(country_levels))) %>%
  ggplot(aes(n_prn_interpretable, interpretability_fraction, label = country_iso3)) +
  geom_point(aes(fill = selection_state), shape = 21, size = 2.7, colour = FIGURE_INK, stroke = 0.22) +
  ggrepel::geom_text_repel(size = 2.0, min.segment.length = 0, segment.size = 0.15, max.overlaps = Inf) +
  scale_x_continuous(trans = "log10", breaks = c(1, 8, 25, 100, 500, 1000), labels = comma) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  scale_fill_manual(values = c(primary_and_triangulated = npg_colors["green"], primary_only = npg_colors["peach"], context_only = FIGURE_GREY), name = "Selection") +
  labs(x = "PRN-interpretable genomes", y = "Interpretability fraction") +
  theme_nature() +
  theme(legend.position = "bottom")

ed1 <- (pA | pB) / pC +
  plot_layout(heights = c(1, 1), widths = c(1.1, 0.9)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = FIGURE_TAG_SIZE, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot",
    plot.margin = margin(3, 3, 3, 3)
  )

save_ed_pdf(ed1, "Extended_Data_Fig_01_Sample_Frame_Reconciliation.pdf", height = NC_MAX_HEIGHT)
save_ed_png(ed1, "Extended_Data_Fig_01_Sample_Frame_Reconciliation.png", height = NC_MAX_HEIGHT)
