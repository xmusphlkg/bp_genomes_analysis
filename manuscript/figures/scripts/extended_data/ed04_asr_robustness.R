#!/usr/bin/env Rscript
# Extended Data Fig. 4: ASR robustness and stochastic mapping

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

registry <- safe_load(file.path(FIGURE_DATA_DIR, "asr_scenario_registry.tsv"), "ASR scenario registry") %>%
  mutate(
    fitch_origin_events = as.numeric(fitch_origin_events),
    largest_disrupted_clade_share = as.numeric(largest_disrupted_clade_share),
    scenario_class = stringr::str_replace_all(scenario_class, "_", " ")
  )

clone <- safe_load(file.path(FIGURE_DATA_DIR, "asr_one_global_clone_summary.tsv"), "one-global-clone summary") %>%
  mutate(across(c(n_scenarios, n_reject_one_global_clone_fitch, min_fitch_origin_events,
                  median_fitch_origin_events, max_fitch_origin_events), as.numeric))

simmap <- safe_load(file.path(FIGURE_DATA_DIR, "asr_stochastic_mapping_summary.tsv"), "stochastic mapping summary") %>%
  mutate(across(c(stochastic_origin_count_median, stochastic_origin_count_lower_95,
                  stochastic_origin_count_upper_95, proportion_origin_count_gt_1), as.numeric),
         scenario_class = stringr::str_replace_all(scenario_class, "_", " "),
         scenario_label = scenario_id %>%
           stringr::str_replace_all("_", " ") %>%
           stringr::str_replace_all("composition filtered", "composition-filtered") %>%
           stringr::str_replace_all("reference rooted", "reference-rooted") %>%
           stringr::str_replace_all("midpoint rooted", "midpoint-rooted") %>%
           stringr::str_replace_all("support ge ", "support >= "))

study <- safe_load(file.path(FIGURE_DATA_DIR, "study_dependence", "asr_study_block_resampling.tsv"), "ASR study-block resampling") %>%
  filter(row_type == "replicate") %>%
  mutate(fitch_origin_events = as.numeric(fitch_origin_events))

pA <- registry %>%
  filter(!is.na(fitch_origin_events)) %>%
  ggplot(aes(fitch_origin_events, scenario_class)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.25, colour = FIGURE_MID_GREY) +
  geom_point(aes(fill = scenario_class), shape = 21, size = 2.0, colour = FIGURE_INK, stroke = 0.18, alpha = 0.85) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = "Fitch origin count", y = NULL) +
  theme_nature()

pB <- clone %>%
  filter(scenario_class != "overall") %>%
  mutate(
    reject_fraction = n_reject_one_global_clone_fitch / n_scenarios,
    scenario_class = stringr::str_replace_all(scenario_class, "_", " ")
  ) %>%
  ggplot(aes(reject_fraction, reorder(scenario_class, reject_fraction))) +
  geom_col(fill = npg_colors["green"], width = 0.68, colour = "white", linewidth = 0.15) +
  scale_x_continuous(labels = percent, limits = c(0, 1)) +
  labs(x = "Fitch scenarios rejecting one-clone", y = NULL) +
  theme_nature()

pC <- simmap %>%
  ggplot(aes(stochastic_origin_count_median, reorder(scenario_label, stochastic_origin_count_median))) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.25, colour = FIGURE_MID_GREY) +
  geom_segment(aes(x = stochastic_origin_count_lower_95, xend = stochastic_origin_count_upper_95,
                   yend = reorder(scenario_label, stochastic_origin_count_median)),
               colour = FIGURE_GREY, linewidth = 0.65) +
  geom_point(size = 2.1, fill = npg_colors["purple"], shape = 21, colour = FIGURE_INK, stroke = 0.18) +
  labs(x = "Stochastic-map origin count", y = NULL) +
  theme_nature(base_size = FIGURE_BASE_SIZE) +
  theme(axis.text.y = element_text(size = 4.5, lineheight = 0.85))

pD <- study %>%
  count(fitch_origin_events, name = "replicates") %>%
  ggplot(aes(fitch_origin_events, replicates)) +
  geom_col(width = 0.82, fill = npg_colors["blue"], colour = "white", linewidth = 0.18) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.25, colour = FIGURE_MID_GREY) +
  scale_x_continuous(breaks = pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = pretty_breaks(n = 4), expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Study-block balanced Fitch origins", y = "Replicates") +
  theme_nature()

ed4 <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1), widths = c(1.1, 0.9)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = FIGURE_TAG_SIZE, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot",
    plot.margin = margin(3, 3, 3, 3)
  )

save_ed_pdf(ed4, "Extended_Data_Fig_04_ASR_Robustness.pdf", height = NC_MAX_HEIGHT)
save_ed_png(ed4, "Extended_Data_Fig_04_ASR_Robustness.png", height = NC_MAX_HEIGHT)
