#!/usr/bin/env Rscript
# Extended Data Fig. 6: country-by-epoch architecture turnover

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

turnover <- safe_load(file.path(FIGURE_DATA_DIR, "selected_country", "country_epoch_architecture_turnover.tsv"), "architecture turnover") %>%
  mutate(
    epoch_event_count = as.numeric(epoch_event_count),
    epoch_total_disrupted = as.numeric(epoch_total_disrupted),
    epoch_event_share = as.numeric(epoch_event_share),
    epoch_class = short_epoch_type(epoch_type),
    epoch_display = stringr::str_wrap(epoch_label, width = 18),
    event = event_label(prn_event_id)
  )

summary <- safe_load(file.path(FIGURE_DATA_DIR, "selected_country", "country_epoch_architecture_turnover_summary.tsv"), "architecture turnover summary") %>%
  mutate(
    architecture_total_variation_distance = as.numeric(architecture_total_variation_distance),
    previous_dominant_share = as.numeric(previous_dominant_share),
    next_dominant_share = as.numeric(next_dominant_share)
  )

country_levels <- c("USA", "NZL", "JPN", "AUS")
event_cols <- c(
  "IS481 gap1043" = unname(npg_colors["red"]),
  "Rearrangement cov58" = unname(npg_colors["blue"]),
  "Rearrangement cov91" = unname(npg_colors["teal"]),
  "gap1042" = unname(npg_colors["peach"]),
  "Insufficient" = FIGURE_GREY,
  "None observed" = FIGURE_PANEL_FILL
)

pA <- turnover %>%
  filter(country_iso3 %in% country_levels) %>%
  mutate(
    country_iso3 = factor(country_iso3, levels = country_levels),
    epoch_display = factor(epoch_display, levels = unique(epoch_display[order(country_iso3, epoch_label)]))
  ) %>%
  ggplot(aes(epoch_display, epoch_event_count, fill = event)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.18) +
  facet_wrap(~country_iso3, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = event_cols, na.value = FIGURE_LIGHT_GREY, name = "Architecture") +
  labs(x = NULL, y = "Disrupted genomes") +
  theme_nature(base_size = 6.4) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, face = "plain"), legend.position = "bottom")

pB <- summary %>%
  filter(comparison_type == "primary_selected_country_pair", country_iso3 %in% country_levels) %>%
  mutate(country_iso3 = factor(country_iso3, levels = rev(country_levels))) %>%
  ggplot(aes(architecture_total_variation_distance, country_iso3)) +
  geom_col(fill = npg_colors["purple"], width = 0.65, colour = "white", linewidth = 0.18) +
  scale_x_continuous(labels = percent, limits = c(0, 1)) +
  labs(x = "Architecture total variation distance", y = NULL) +
  theme_nature()

pC <- summary %>%
  filter(comparison_type == "primary_selected_country_pair", country_iso3 %in% country_levels) %>%
  transmute(
    country_iso3 = factor(country_iso3, levels = rev(country_levels)),
    Previous = previous_dominant_share,
    Next = next_dominant_share
  ) %>%
  pivot_longer(-country_iso3, names_to = "epoch", values_to = "share") %>%
  ggplot(aes(share, country_iso3, colour = epoch)) +
  geom_point(size = 2.4) +
  scale_x_continuous(labels = percent, limits = c(0, 1), expand = expansion(mult = c(0.04, 0.10))) +
  scale_colour_manual(values = c(Previous = FIGURE_MID_GREY, Next = unname(npg_colors["red"])), name = "Dominant share") +
  labs(x = "Dominant event share", y = NULL) +
  theme_nature() +
  theme(legend.position = "bottom", plot.margin = margin(3, 8, 3, 3))

ed6 <- pA / (pB | pC) +
  plot_layout(heights = c(1.12, 0.9)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = FIGURE_TAG_SIZE, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot",
    plot.margin = margin(3, 3, 3, 3)
  )

save_ed_pdf(ed6, "Extended_Data_Fig_06_Architecture_Turnover.pdf", height = NC_MAX_HEIGHT)
save_ed_png(ed6, "Extended_Data_Fig_06_Architecture_Turnover.png", height = NC_MAX_HEIGHT)
