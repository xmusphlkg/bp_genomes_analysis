#!/usr/bin/env Rscript
# Supplementary Figure 14: structural recurrence evidence compendium.

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
library(stringr)
library(tidyr)

parse_top3_share <- function(x) {
  out <- str_match(x, "share=([0-9.]+)")[, 2]
  as.numeric(out)
}

hierarchy <- safe_load(file.path(FIGURE_DATA_DIR, "event_definition_hierarchy_sensitivity.tsv"), "event-definition hierarchy") %>%
  mutate(
    dominant_definition_share = as.numeric(dominant_definition_share),
    top3_definition_share = as.numeric(top3_definition_share),
    event_definition_label = case_when(
      event_definition_tier == "exact_breakpoint_orientation_tsd" ~ "Exact event grammar",
      event_definition_tier == "breakpoint_window_5bp" ~ "Breakpoint window",
      event_definition_tier == "mechanism_approximate_architecture" ~ "Approx. architecture",
      event_definition_tier == "broad_mechanism_class" ~ "Broad mechanism",
      TRUE ~ str_wrap(event_definition_label, 18)
    )
  ) %>%
  mutate(event_definition_label = factor(event_definition_label, levels = event_definition_label))

grammar <- safe_load(file.path(FIGURE_DATA_DIR, "structural_grammar_evidence.tsv"), "structural grammar evidence") %>%
  mutate(
    dominant_event_share = as.numeric(dominant_event_share),
    top3_share = parse_top3_share(top3_share_or_summary),
    evidence_label = case_when(
      evidence_layer == "raw_structurally_resolved_event_burden" ~ "Raw resolved genomes",
      collapse_or_weighting_rule == "country_x_st" ~ "Country + ST collapse",
      collapse_or_weighting_rule == "country_x_st_x_ptxp_fim_signature" ~ "Country + lineage collapse",
      collapse_or_weighting_rule == "st_only" ~ "ST-only collapse",
      evidence_layer == "origin_package_collapse" ~ "Tree-package collapse",
      evidence_layer == "within_origin_concentration" ~ "Within-package concentration",
      collapse_or_weighting_rule == "study_block_equalized" ~ "Equal study block",
      collapse_or_weighting_rule == "drop_largest_block_naive" ~ "Drop largest block",
      TRUE ~ str_replace_all(collapse_or_weighting_rule, "_", " ")
    )
  ) %>%
  filter(!is.na(dominant_event_share)) %>%
  mutate(evidence_label = factor(evidence_label, levels = rev(evidence_label)))

events <- safe_load(file.path(FIGURE_DATA_DIR, "event_specific_acquisition_summary.tsv"), "event-specific acquisition summary") %>%
  mutate(
    rank_by_genome_burden = as.numeric(rank_by_genome_burden),
    sample_count = as.numeric(sample_count),
    sample_share_among_structurally_resolved = as.numeric(sample_share_among_structurally_resolved),
    acquisition_package_count = as.numeric(acquisition_package_count),
    n_country_year_cells = as.numeric(n_country_year_cells),
    event_display = paste0(rank_by_genome_burden, ". ", event_label(prn_event_id))
  ) %>%
  arrange(desc(sample_count)) %>%
  slice_head(n = 10) %>%
  mutate(event_display = factor(event_display, levels = rev(event_display)))

junction <- safe_load(file.path(FIGURE_DATA_DIR, "prn_junction_confidence_matrix.tsv"), "junction confidence matrix") %>%
  mutate(
    sample_count = as.numeric(sample_count),
    event_display = paste0(row_number(), ". ", event_label(prn_event_id)),
    confidence_display = case_when(
      confidence_tier == "tier_1_read_backed_tsd_recovered" ~ "Read-backed TSD",
      confidence_tier == "tier_2_public_longread_or_hybrid" ~ "Long-read / hybrid",
      confidence_tier == "tier_3_assembly_or_rule_supported" ~ "Assembly/rule",
      confidence_tier == "tier_4_validation_unresolved" ~ "Unresolved",
      TRUE ~ str_replace_all(confidence_tier, "_", " ")
    )
  ) %>%
  arrange(desc(sample_count)) %>%
  slice_head(n = 10) %>%
  mutate(event_display = factor(event_display, levels = rev(event_display)))

pA <- hierarchy %>%
  pivot_longer(
    cols = c(dominant_definition_share, top3_definition_share),
    names_to = "statistic",
    values_to = "share"
  ) %>%
  mutate(statistic = recode(
    statistic,
    dominant_definition_share = "Dominant event",
    top3_definition_share = "Top three events"
  )) %>%
  ggplot(aes(event_definition_label, share, group = statistic, colour = statistic)) +
  geom_line(linewidth = 0.35) +
  geom_point(size = 1.8) +
  scale_y_continuous(labels = percent, limits = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
  scale_colour_manual(values = c(
    "Dominant event" = unname(npg_colors["red"]),
    "Top three events" = unname(npg_colors["blue"])
  ), name = NULL) +
  labs(x = NULL, y = "Share of resolved disruptions") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 28, hjust = 1), legend.position = "bottom")

pB <- grammar %>%
  ggplot(aes(dominant_event_share, evidence_label)) +
  geom_col(width = 0.62, fill = npg_colors["red"], alpha = 0.86, colour = "white", linewidth = 0.14) +
  geom_point(aes(x = top3_share), size = 2.0, fill = npg_colors["blue"], shape = 21, colour = FIGURE_INK, stroke = 0.18, na.rm = TRUE) +
  scale_x_continuous(labels = percent, limits = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
  labs(x = "Dominant share; points mark top-three share", y = NULL) +
  theme_nature(base_size = 6.0) +
  theme(axis.text.y = element_text(size = 4.8, lineheight = 0.88))

pC <- events %>%
  ggplot(aes(sample_share_among_structurally_resolved, event_display)) +
  geom_point(aes(size = acquisition_package_count, fill = n_country_year_cells),
             shape = 21, colour = FIGURE_INK, stroke = 0.18, alpha = 0.92) +
  scale_x_continuous(labels = percent, limits = c(0, 0.70), expand = expansion(mult = c(0, 0.07))) +
  scale_size_continuous(range = c(1.4, 5.0), breaks = c(0, 1, 3, 9), name = "Minimum\npackages") +
  scale_fill_gradient(low = FIGURE_LIGHT_GREY, high = npg_colors["green"], name = "Country-year\ncells") +
  labs(x = "Share of resolved disruptions", y = NULL) +
  theme_nature(base_size = 6.0) +
  guides(size = "none", fill = "none") +
  theme(axis.text.y = element_text(size = 4.8), legend.position = "none")

pD <- junction %>%
  ggplot(aes(sample_count, event_display, fill = confidence_display)) +
  geom_col(width = 0.66, colour = "white", linewidth = 0.14) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.06))) +
  scale_fill_manual(
    values = c(
      "Read-backed TSD" = unname(npg_colors["green"]),
      "Long-read / hybrid" = unname(npg_colors["blue"]),
      "Assembly/rule" = unname(npg_colors["peach"]),
      "Unresolved" = FIGURE_GREY
    ),
    name = "Confidence tier"
  ) +
  labs(x = "Genomes in event class", y = NULL) +
  theme_nature(base_size = 6.0) +
  guides(fill = "none") +
  theme(axis.text.y = element_text(size = 4.8), legend.position = "none")

ed14 <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(0.94, 1.08), widths = c(1.02, 1.08), guides = "keep") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = FIGURE_TAG_SIZE, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot",
    plot.margin = margin(3, 3, 3, 3),
    legend.position = "bottom"
  )

save_ed_pdf(ed14, "Extended_Data_Fig_14_Structural_Recurrence_Compendium.pdf", height = NC_MAX_HEIGHT)
save_ed_png(ed14, "Extended_Data_Fig_14_Structural_Recurrence_Compendium.png", height = NC_MAX_HEIGHT)
