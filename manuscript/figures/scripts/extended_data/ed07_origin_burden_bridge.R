#!/usr/bin/env Rscript
# Extended Data Fig. 7: amplification heterogeneity after origin or detection

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

burden <- safe_load(file.path(FIGURE_DATA_DIR, "selected_country", "origin_burden_prevalence_shift.tsv"), "origin burden prevalence shift") %>%
  mutate(
    n_local_origin_packages = as.numeric(n_local_origin_packages),
    total_package_descendants = as.numeric(total_package_descendants),
    post_minus_pre_ipw_prevalence = as.numeric(post_minus_pre_ipw_prevalence),
    has_local_origin_package = as_logical_flag(has_local_origin_package),
    direction = direction_label(prevalence_direction)
  )

relative <- safe_load(data_paths$rr_relative_year_plot_data, "relative-year plot data") %>%
  mutate(
    relative_year = as.numeric(relative_year),
    ipw_prevalence = as.numeric(ipw_prevalence),
    n_genomes_prn_interpretable = as.numeric(n_genomes_prn_interpretable),
    event_type = recode(event_type, first_local_origin = "First local origin", first_prn_detection = "First detection")
  )

country_levels <- c("USA", "NZL", "JPN", "AUS", "CHN", "FRA", "BRA")

pA <- burden %>%
  filter(country_iso3 %in% country_levels) %>%
  mutate(country_iso3 = factor(country_iso3, levels = rev(country_levels))) %>%
  ggplot(aes(total_package_descendants, country_iso3, fill = has_local_origin_package)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.18) +
  scale_fill_manual(values = c(`TRUE` = unname(npg_colors["green"]), `FALSE` = FIGURE_GREY), guide = "none") +
  labs(x = "Total package descendants", y = NULL) +
  theme_nature()

pB <- burden %>%
  filter(country_iso3 %in% country_levels, !is.na(post_minus_pre_ipw_prevalence), !is.na(total_package_descendants)) %>%
  ggplot(aes(post_minus_pre_ipw_prevalence, total_package_descendants, label = country_iso3)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = FIGURE_MID_GREY) +
  geom_point(aes(fill = direction, size = pmax(n_local_origin_packages, 1)), shape = 21, colour = FIGURE_INK, stroke = 0.22, alpha = 0.92) +
  ggrepel::geom_text_repel(size = 2.0, min.segment.length = 0, segment.size = 0.13, max.overlaps = Inf) +
  scale_fill_manual(
    values = direction_colors,
    breaks = c("Uncertain", "Upward", "Downward", "No change"),
    labels = c("Uncertain", "Up", "Down", "Flat"),
    name = "Direction"
  ) +
  scale_size_area(max_size = 4.6, breaks = 1:5, name = "Local packages") +
  scale_x_continuous(labels = percent, expand = expansion(mult = c(0.04, 0.12))) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.10))) +
  guides(
    size = guide_legend(nrow = 1, title.position = "left", order = 1),
    fill = guide_legend(nrow = 1, title.position = "left", order = 2)
  ) +
  labs(x = "Post-minus-pre IPW disrupted fraction", y = "Package descendants") +
  theme_nature() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.justification = "left",
    legend.key.size = unit(5, "pt"),
    legend.spacing.x = unit(1.5, "pt"),
    legend.text = element_text(size = 5.1),
    legend.title = element_text(face = "bold", size = 5.2),
    plot.margin = margin(3, 9, 3, 3)
  )

pC <- relative %>%
  filter(country_iso3 %in% c("USA", "NZL", "JPN", "AUS", "CHN"), !is.na(ipw_prevalence)) %>%
  mutate(country_iso3 = factor(country_iso3, levels = c("USA", "NZL", "JPN", "AUS", "CHN"))) %>%
  ggplot(aes(relative_year, ipw_prevalence, colour = event_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = FIGURE_MID_GREY, linewidth = 0.25) +
  geom_line(linewidth = 0.42, na.rm = TRUE) +
  geom_point(aes(size = pmax(n_genomes_prn_interpretable, 1)), na.rm = TRUE) +
  facet_wrap(~country_iso3, nrow = 1) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  scale_x_continuous(breaks = -3:3) +
  scale_colour_manual(values = c("First local origin" = unname(npg_colors["red"]), "First detection" = unname(npg_colors["blue"])), name = NULL) +
  scale_size_area(max_size = 4, guide = "none") +
  labs(x = "Relative year", y = "IPW disrupted fraction") +
  theme_nature(base_size = FIGURE_BASE_SIZE) +
  theme(legend.position = "bottom")

ed7 <- (pA | pB) / pC +
  plot_layout(heights = c(0.95, 1.05), widths = c(0.85, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = FIGURE_TAG_SIZE, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot",
    plot.margin = margin(3, 9, 3, 3)
  )

save_ed_pdf(ed7, "Extended_Data_Fig_07_Origin_Burden_Bridge.pdf", height = NC_MAX_HEIGHT)
save_ed_png(ed7, "Extended_Data_Fig_07_Origin_Burden_Bridge.png", height = NC_MAX_HEIGHT)
