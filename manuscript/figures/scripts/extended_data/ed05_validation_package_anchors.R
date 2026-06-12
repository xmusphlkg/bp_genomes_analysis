#!/usr/bin/env Rscript
# Extended Data Fig. 5: validation anchors and package-level rerun audit

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

validation <- safe_load(data_paths$rr_validation_matrix, "selected-country validation matrix") %>%
  mutate(
    country_sample_count = as.numeric(country_sample_count),
    event = event_label(prn_event_id),
    support_scope = if_else(is.na(support_scope) | support_scope == "", "none", support_scope)
  )

packages <- safe_load(data_paths$rr_origin_package_summary, "origin package summary") %>%
  mutate(
    n_disrupted_descendants = as.numeric(n_disrupted_descendants),
    n_total_descendants = as.numeric(n_total_descendants),
    origin_year = as.numeric(origin_year),
    origin_package_hard_anchor = tidyr::replace_na(as_logical_flag(origin_package_hard_anchor), FALSE),
    dominant_event_hard_anchor = tidyr::replace_na(as_logical_flag(dominant_event_hard_anchor), FALSE)
  )

local <- safe_load(file.path(FIGURE_DATA_DIR, "local_rooted_package_tree_summary.tsv"), "local rooted package tree summary") %>%
  mutate(across(c(local_tip_count, local_disrupted_tip_count, local_fitch_origin_events, target_package_tip_retention_fraction), as.numeric))

pA <- validation %>%
  filter(!str_detect(event, "Insufficient")) %>%
  mutate(
    country_iso3 = factor(country_iso3, levels = rev(ordered_countries(country_iso3))),
    event = factor(event, levels = rev(unique(event)))
  ) %>%
  ggplot(aes(event, country_iso3)) +
  geom_point(aes(size = country_sample_count, fill = support_scope), shape = 21, colour = FIGURE_INK, stroke = 0.2) +
  scale_size_area(max_size = 5.4, name = "Country event count") +
  scale_fill_manual(
    values = c(same_country = unname(npg_colors["green"]), cross_country = unname(npg_colors["blue"]), none = FIGURE_LIGHT_GREY),
    labels = c(same_country = "Same country", cross_country = "Cross country", none = "None"),
    name = "Anchor scope"
  ) +
  labs(x = NULL, y = NULL) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 5.6, face = "plain"), legend.position = "bottom")

pB <- packages %>%
  filter(!is.na(origin_id)) %>%
  mutate(origin_id = factor(origin_id, levels = origin_id[order(n_disrupted_descendants)])) %>%
  ggplot(aes(n_disrupted_descendants, origin_id)) +
  geom_segment(aes(x = 0, xend = n_total_descendants, yend = origin_id), colour = FIGURE_RULE_COLOUR, linewidth = 0.8) +
  geom_point(aes(fill = origin_package_hard_anchor), size = 2.8, shape = 21, colour = FIGURE_INK, stroke = 0.22) +
  scale_fill_manual(values = c(`TRUE` = unname(npg_colors["green"]), `FALSE` = FIGURE_GREY), labels = c(`TRUE` = "Yes", `FALSE` = "No"), name = "Hard anchor") +
  labs(x = "Disrupted descendants; grey segment = total descendants", y = NULL) +
  theme_nature()

pC <- local %>%
  mutate(origin_id = factor(origin_id, levels = origin_id[order(local_disrupted_tip_count)])) %>%
  ggplot(aes(local_disrupted_tip_count, origin_id)) +
  geom_segment(aes(x = 0, xend = local_tip_count, yend = origin_id), colour = FIGURE_RULE_COLOUR, linewidth = 0.8) +
  geom_point(aes(fill = local_origin_consistency_status), shape = 21, size = 2.8, colour = FIGURE_INK, stroke = 0.22) +
  scale_fill_manual(
    values = c(
      single_origin_consistent = unname(npg_colors["green"]),
      split_across_multiple_local_origins = unname(npg_colors["peach"]),
      no_local_origin_call = FIGURE_GREY
    ),
    labels = c(
      single_origin_consistent = "Single origin",
      split_across_multiple_local_origins = "Split origins",
      no_local_origin_call = "No local call"
    ),
    name = "Local rerun"
  ) +
  labs(x = "Local disrupted tips; grey segment = local tips", y = NULL) +
  theme_nature()

ed5 <- pA / (pB | pC) +
  plot_layout(heights = c(1.05, 0.95)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = FIGURE_TAG_SIZE, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot",
    plot.margin = margin(3, 3, 3, 3)
  )

save_ed_pdf(ed5, "Extended_Data_Fig_05_Validation_Package_Anchors.pdf", height = NC_MAX_HEIGHT)
save_ed_png(ed5, "Extended_Data_Fig_05_Validation_Package_Anchors.png", height = NC_MAX_HEIGHT)
