#!/usr/bin/env Rscript
# Supplementary Figure 15: validation and caller-sensitivity compendium.

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

concordance <- safe_load(file.path(FIGURE_DATA_DIR, "published_overlap_concordance.tsv"), "published-overlap concordance") %>%
  mutate(
    n_compared_rows = as.numeric(n_compared_rows),
    concordance_fraction = as.numeric(concordance_fraction)
  ) %>%
  filter(
    metric_name == "prn_status_concordance",
    summary_level == "country",
    !is.na(country_iso3),
    country_iso3 != "",
    n_compared_rows >= 10
  ) %>%
  arrange(concordance_fraction, desc(n_compared_rows)) %>%
  mutate(country_iso3 = factor(country_iso3, levels = country_iso3))

validation <- safe_load(file.path(FIGURE_DATA_DIR, "caller_validation_sensitivity_summary.tsv"), "caller validation sensitivity") %>%
  mutate(
    fraction = suppressWarnings(as.numeric(fraction)),
    n_compared_rows = as.numeric(n_compared_rows),
    evidence_display = case_when(
      evidence_layer == "PRN status" ~ "PRN status",
      evidence_layer == "broad mechanism" ~ "Broad mechanism",
      evidence_layer == "manual/read junction evidence tiers" ~ "Junction tier",
      evidence_layer == "event-definition hierarchy" ~ "Event hierarchy",
      evidence_layer == "full caller threshold grid" ~ "Threshold grid",
      TRUE ~ str_wrap(evidence_layer, 18)
    ),
    audit_display = case_when(
      audit_section == "caller_validation" ~ "Caller validation",
      audit_section == "read_validation" ~ "Manual/read validation",
      audit_section == "threshold_sensitivity" ~ "Threshold sensitivity",
      audit_section == "false_negative_risk" ~ "False-negative risk",
      TRUE ~ str_replace_all(audit_section, "_", " ")
    )
  ) %>%
  filter(!is.na(fraction)) %>%
  mutate(evidence_display = factor(evidence_display, levels = rev(evidence_display)))

threshold <- safe_load(file.path(FIGURE_DATA_DIR, "prn_threshold_grid_full.tsv"), "full threshold grid") %>%
  mutate(
    hsp_min_pident = as.numeric(hsp_min_pident),
    locus_qcov_threshold = as.numeric(locus_qcov_threshold),
    dominant_event_share = as.numeric(dominant_event_share),
    top3_event_share = as.numeric(top3_event_share),
    status_changed_vs_manuscript_fraction = as.numeric(status_changed_vs_manuscript_fraction),
    is_support_profile = str_to_title(is_support_profile)
  )

external <- safe_load(file.path(FIGURE_DATA_DIR, "biology_bridge_external_context.tsv"), "external phenotype context") %>%
  mutate(
    n_total_or_frame = as.numeric(n_total_or_frame),
    n_prn_deficient = as.numeric(n_prn_deficient),
    deficient_fraction = n_prn_deficient / n_total_or_frame,
    context_label = case_when(
      str_detect(context_id, "japan") ~ "Japan 1990-2010",
      str_detect(context_id, "usa") ~ "USA 1935-2012",
      str_detect(context_id, "australia") ~ "Australia 2008-2012",
      str_detect(context_id, "europe") ~ "Europe 1998-2015",
      str_detect(context_id, "france") ~ "France 2023-2024",
      str_detect(context_id, "belgium_2022") ~ "Belgium 2022-2023",
      str_detect(context_id, "belgium_2000") ~ "Belgium 2000-2023",
      TRUE ~ str_wrap(country_or_region, 16)
    )
  ) %>%
  filter(!is.na(deficient_fraction), !is.na(n_total_or_frame), n_total_or_frame > 0) %>%
  arrange(deficient_fraction) %>%
  mutate(context_label = factor(context_label, levels = context_label))

pA <- concordance %>%
  ggplot(aes(concordance_fraction, country_iso3)) +
  geom_col(width = 0.68, fill = npg_colors["green"], alpha = 0.82, colour = "white", linewidth = 0.14) +
  geom_vline(xintercept = 0.93, linetype = "dashed", linewidth = 0.25, colour = FIGURE_MID_GREY) +
  scale_x_continuous(labels = percent, limits = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
  labs(x = "Caller/published PRN-status concordance", y = NULL) +
  theme_nature(base_size = 6.0) +
  guides(fill = "none") +
  theme(axis.text.y = element_text(size = 5.1), legend.position = "none")

pB <- validation %>%
  ggplot(aes(fraction, evidence_display)) +
  geom_segment(aes(x = 0, xend = fraction, yend = evidence_display), linewidth = 0.45, colour = FIGURE_GREY) +
  geom_point(aes(fill = audit_display, size = n_compared_rows), shape = 21, colour = FIGURE_INK, stroke = 0.18, alpha = 0.92) +
  scale_x_continuous(labels = percent, limits = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
  scale_fill_manual(values = c(
    "Caller validation" = unname(npg_colors["blue"]),
    "Manual/read validation" = unname(npg_colors["green"]),
    "Threshold sensitivity" = unname(npg_colors["peach"]),
    "False-negative risk" = unname(npg_colors["grey"])
  ), name = NULL) +
  scale_size_continuous(range = c(1.8, 4.2), name = "Compared\nrecords") +
  labs(x = "Audit-support fraction", y = NULL) +
  guides(fill = "none", size = "none") +
  theme_nature(base_size = 6.0) +
  theme(legend.position = "none")

pC <- threshold %>%
  select(grid_id, hsp_min_pident, locus_qcov_threshold, is_support_profile,
         dominant_event_share, top3_event_share) %>%
  pivot_longer(c(dominant_event_share, top3_event_share), names_to = "statistic", values_to = "share") %>%
  mutate(statistic = recode(
    statistic,
    dominant_event_share = "Dominant event",
    top3_event_share = "Top three events"
  )) %>%
  ggplot(aes(locus_qcov_threshold, share, colour = is_support_profile, linetype = statistic)) +
  geom_line(linewidth = 0.34) +
  geom_point(size = 1.2) +
  facet_wrap(~hsp_min_pident, nrow = 1, labeller = labeller(hsp_min_pident = function(x) paste0("HSP ", x, "%"))) +
  scale_y_continuous(labels = percent, limits = c(0.55, 1.0), expand = expansion(mult = c(0.02, 0.04))) +
  scale_x_continuous(breaks = c(90, 92.5, 95, 97.5, 99)) +
  scale_colour_manual(values = c(
    Relaxed = unname(npg_colors["teal"]),
    Default = unname(npg_colors["blue"]),
    Strict = unname(npg_colors["red"])
  ), name = "IS profile") +
  scale_linetype_manual(values = c("Dominant event" = "solid", "Top three events" = "dashed"), name = NULL) +
  labs(x = "Locus coverage threshold (%)", y = "Share of structural calls") +
  theme_nature(base_size = 5.8) +
  guides(colour = "none", linetype = "none") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), strip.text = element_text(size = 5.3), legend.position = "none")

pD <- external %>%
  ggplot(aes(deficient_fraction, context_label)) +
  geom_col(fill = npg_colors["red"], alpha = 0.82, width = 0.66, colour = "white", linewidth = 0.14) +
  geom_text(aes(label = comma(n_total_or_frame)), hjust = -0.08, size = 1.55, colour = FIGURE_MUTED_TEXT) +
  scale_x_continuous(labels = percent, limits = c(0, 1.05), expand = expansion(mult = c(0, 0.04))) +
  labs(x = "Published PRN-deficient fraction; labels show tested isolates", y = NULL) +
  theme_nature(base_size = 6.0) +
  theme(axis.text.y = element_text(size = 4.9, lineheight = 0.86))

ed15 <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1.08), widths = c(1.02, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = FIGURE_TAG_SIZE, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot",
    plot.margin = margin(3, 3, 3, 3),
    legend.position = "bottom"
  )

save_ed_pdf(ed15, "Extended_Data_Fig_15_Validation_Sensitivity_Compendium.pdf", height = NC_MAX_HEIGHT)
save_ed_png(ed15, "Extended_Data_Fig_15_Validation_Sensitivity_Compendium.png", height = NC_MAX_HEIGHT)
