#!/usr/bin/env Rscript
# Extended Data Fig. 9: support-only product-aware programme diagnostics

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

models <- safe_load(data_paths$fig5_models, "programme exposure models") %>%
  mutate(across(c(effect_estimate, ci_lower, ci_upper), as.numeric))

loo <- safe_load(data_paths$fig5_leave_one_out, "leave-one-country-out exposure models") %>%
  mutate(effect_estimate = as.numeric(effect_estimate), same_direction_as_primary = as_logical_flag(same_direction_as_primary))

coverage <- safe_load(data_paths$fig5_formulation, "formulation coverage") %>%
  mutate(
    known_prn_fraction = as.numeric(known_prn_fraction),
    role_product_metadata_fraction = as.numeric(role_product_metadata_fraction),
    mean_primary_prn_positive_share = as.numeric(mean_primary_prn_positive_share)
  )

class_summary <- safe_load(data_paths$fig5_programme_class_summary, "programme class summary") %>%
  mutate(
    pooled_ipw_prevalence = as.numeric(pooled_ipw_prevalence),
    mean_period_ipw_prevalence = as.numeric(mean_period_ipw_prevalence)
  )

ladder <- tibble::tribble(
  ~exposure, ~layer, ~x,
  "DTP3", "coverage", 1,
  "V1", "timing", 2,
  "V2", "PRN", 3,
  "V3", "products", 4
)

pA <- ggplot(ladder, aes(x, 1)) +
  annotate("segment", x = 1, xend = 4, y = 1, yend = 1, colour = FIGURE_RULE_COLOUR, linewidth = 1.0) +
  geom_point(aes(fill = exposure), shape = 21, size = 5, colour = FIGURE_INK, stroke = 0.25) +
  geom_text(aes(label = exposure), nudge_y = 0.16, size = 2.1, fontface = "bold") +
  geom_text(aes(label = layer), nudge_y = -0.16, size = 1.9, colour = FIGURE_MUTED_TEXT) +
  scale_fill_manual(values = c("DTP3" = FIGURE_GREY, "V1" = unname(npg_colors["blue"]), "V2" = unname(npg_colors["peach"]), "V3" = unname(npg_colors["green"])), guide = "none") +
  coord_cartesian(xlim = c(0.6, 4.4), ylim = c(0.55, 1.45), clip = "off") +
  theme_void(base_size = FIGURE_BASE_SIZE)

pB <- models %>%
  filter(panel_id %in% c("primary_exposure_comparison", "cluster_robust_v3", "cluster_robust_v2")) %>%
  mutate(model_label = str_replace_all(sensitivity_label, "_", " "),
         model_label = factor(model_label, levels = rev(model_label))) %>%
  ggplot(aes(effect_estimate, model_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = FIGURE_MID_GREY) +
  geom_segment(aes(x = ci_lower, xend = ci_upper, yend = model_label), colour = FIGURE_GREY, linewidth = 0.8) +
  geom_point(aes(fill = focal_exposure_family), shape = 21, size = 2.6, colour = FIGURE_INK, stroke = 0.22) +
  scale_fill_manual(values = c(v1 = unname(npg_colors["blue"]), v2 = unname(npg_colors["peach"]), v3 = unname(npg_colors["green"]), dtp3 = FIGURE_GREY), name = "Exposure") +
  labs(x = "Log-odds effect", y = NULL) +
  theme_nature() +
  theme(legend.position = "bottom")

pC <- loo %>%
  filter(focal_exposure_family == "v3") %>%
  mutate(excluded_country_iso3 = factor(excluded_country_iso3, levels = excluded_country_iso3[order(effect_estimate)])) %>%
  ggplot(aes(effect_estimate, excluded_country_iso3, fill = same_direction_as_primary)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = FIGURE_MID_GREY) +
  scale_fill_manual(values = c(`TRUE` = unname(npg_colors["green"]), `FALSE` = unname(npg_colors["red"])), name = "Same direction") +
  labs(x = "Leave-one-country-out V3 effect", y = "Excluded country") +
  theme_nature() +
  theme(legend.position = "bottom")

pD <- coverage %>%
  filter(country_iso3 %in% c("USA", "NZL", "JPN", "AUS", "CHN", "GBR", "FRA", "BRA")) %>%
  mutate(country_iso3 = factor(country_iso3, levels = rev(c("USA", "NZL", "JPN", "AUS", "CHN", "GBR", "FRA", "BRA")))) %>%
  ggplot(aes(known_prn_fraction, country_iso3)) +
  geom_segment(aes(x = 0, xend = role_product_metadata_fraction, yend = country_iso3), colour = FIGURE_RULE_COLOUR, linewidth = 0.8) +
  geom_point(aes(fill = dominant_prn_in_vaccine_curated), shape = 21, size = 2.8, colour = FIGURE_INK, stroke = 0.22) +
  scale_x_continuous(labels = percent, limits = c(0, 1)) +
  scale_fill_manual(values = c(yes = unname(npg_colors["orange"]), no = unname(npg_colors["blue"]), mixed = unname(npg_colors["peach"]), unknown = FIGURE_GREY), name = "Dominant PRN") +
  labs(x = "Known PRN formulation", y = NULL) +
  theme_nature() +
  theme(legend.position = "bottom")

pE <- if ("program_formulation_class" %in% names(class_summary) && "pooled_ipw_prevalence" %in% names(class_summary)) {
  class_summary %>%
    mutate(program_class_label = recode(program_formulation_class,
      routine_ap_prn_positive = "aP + PRN",
      routine_ap_mixed = "Mixed aP",
      routine_ap_prn_negative = "aP - PRN",
      wp_only_or_pre_ap = "wP / pre-aP"
    )) %>%
    ggplot(aes(pooled_ipw_prevalence, reorder(program_class_label, pooled_ipw_prevalence))) +
    geom_col(fill = unname(npg_colors["blue"]), width = 0.65, colour = "white", linewidth = 0.15) +
    scale_x_continuous(labels = percent) +
    labs(x = "Mean IPW prn-disruption fraction", y = NULL) +
    theme_nature()
  } else {
  tibble::tibble(program_class_label = "Programme summary unavailable", pooled_ipw_prevalence = 0) %>%
    ggplot(aes(pooled_ipw_prevalence, program_class_label)) +
    geom_col(fill = FIGURE_LIGHT_GREY, width = 0.65, colour = "white", linewidth = 0.15) +
    scale_x_continuous(labels = percent, limits = c(0, 1)) +
    labs(x = "Mean IPW prn-disruption fraction", y = NULL) +
    theme_nature()
}

ed9 <- (pA | pB) / (pC | pD | pE) +
  plot_layout(heights = c(0.9, 1.1), widths = c(1, 1, 0.8)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = FIGURE_TAG_SIZE, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot",
    plot.margin = margin(3, 3, 3, 3)
  )

save_ed_pdf(ed9, "Extended_Data_Fig_09_Ecology_Sidecar.pdf", height = NC_MAX_HEIGHT)
save_ed_png(ed9, "Extended_Data_Fig_09_Ecology_Sidecar.png", height = NC_MAX_HEIGHT)
