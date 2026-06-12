#!/usr/bin/env Rscript
# Extended Data Fig. 3: tree representativeness

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

rep <- safe_load(file.path(FIGURE_DATA_DIR, "asr_representativeness_adjustment_summary.tsv"), "ASR representativeness adjustment") %>%
  mutate(across(c(full_interpretable_fraction_largest_gap, tree_subset_fraction_largest_gap,
                  max_absolute_fraction_gap, mean_absolute_fraction_gap, full_interpretable_count,
                  tree_subset_count, tree_subset_fraction, median_fitch_origin_events,
                  min_fitch_origin_events, max_fitch_origin_events), as.numeric))

gap <- rep %>% filter(row_type == "representativeness_gap_summary")
detail <- rep %>% filter(!is.na(comparison_dimension), !is.na(category))
detail <- detail %>%
  group_by(comparison_dimension) %>%
  mutate(
    full_interpretable_fraction = full_interpretable_count / sum(full_interpretable_count, na.rm = TRUE),
    tree_subset_fraction = tree_subset_count / sum(tree_subset_count, na.rm = TRUE)
  ) %>%
  ungroup()
balanced <- rep %>% filter(row_type == "composition_balanced_resampling_summary")

pA <- gap %>%
  mutate(scenario_or_dimension = factor(scenario_or_dimension, levels = scenario_or_dimension[order(max_absolute_fraction_gap)])) %>%
  ggplot(aes(max_absolute_fraction_gap, scenario_or_dimension)) +
  geom_col(fill = npg_colors["blue"], width = 0.65, colour = "white", linewidth = 0.15) +
  geom_text(aes(label = largest_gap_category), hjust = -0.05, size = 2.1) +
  scale_x_continuous(labels = percent, limits = c(0, max(gap$max_absolute_fraction_gap, na.rm = TRUE) * 1.18)) +
  labs(x = "Largest full-versus-tree fraction gap", y = NULL) +
  theme_nature()

pB <- detail %>%
  filter(comparison_dimension %in% c("country", "year_band", "mechanism_group")) %>%
  group_by(comparison_dimension) %>%
  slice_max(abs(tree_subset_fraction - full_interpretable_fraction), n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(category = factor(category, levels = rev(unique(category)))) %>%
  ggplot(aes(full_interpretable_fraction, tree_subset_fraction, label = category)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, linetype = "dashed", colour = FIGURE_MID_GREY) +
  geom_point(size = 2.1, fill = npg_colors["green"], shape = 21, colour = FIGURE_INK, stroke = 0.18) +
  ggrepel::geom_text_repel(size = 1.8, min.segment.length = 0, segment.size = 0.12, max.overlaps = Inf) +
  facet_wrap(~comparison_dimension, nrow = 1, scales = "free") +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(labels = percent) +
  labs(x = "Full interpretable fraction", y = "Tree subset fraction") +
  theme_nature(base_size = FIGURE_BASE_SIZE)

pC <- balanced %>%
  mutate(scenario_or_dimension = factor(scenario_or_dimension, levels = scenario_or_dimension[order(median_fitch_origin_events)])) %>%
  ggplot(aes(median_fitch_origin_events, scenario_or_dimension)) +
  geom_segment(aes(x = min_fitch_origin_events, xend = max_fitch_origin_events, yend = scenario_or_dimension), colour = FIGURE_GREY, linewidth = 0.8) +
  geom_point(size = 2.8, fill = npg_colors["red"], shape = 21, colour = FIGURE_INK, stroke = 0.22) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = FIGURE_MID_GREY, linewidth = 0.25) +
  labs(x = "Fitch origin count after balancing", y = NULL) +
  theme_nature()

ed3 <- pA + free(pC) + pB +
  plot_layout(design = "AB\nCC", heights = c(0.86, 1.1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = FIGURE_TAG_SIZE, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot",
    plot.margin = margin(3, 3, 3, 3)
  )

save_ed_pdf(ed3, "Extended_Data_Fig_03_Tree_Representativeness.pdf", height = NC_MAX_HEIGHT)
save_ed_png(ed3, "Extended_Data_Fig_03_Tree_Representativeness.png", height = NC_MAX_HEIGHT)
