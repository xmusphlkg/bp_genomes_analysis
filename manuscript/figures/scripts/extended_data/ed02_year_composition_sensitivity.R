#!/usr/bin/env Rscript
# Extended Data Fig. 2: missingness, year composition, and study-dependence audit

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

year_sens <- safe_load(file.path(FIGURE_DATA_DIR, "selected_country", "selected_country_year_sensitivity_summary.tsv"), "year sensitivity") %>%
  mutate(across(where(is.character), ~na_if(.x, ""))) %>%
  mutate(
    baseline_naive_delta = as.numeric(baseline_naive_delta),
    max_leave_one_year_naive_delta_change = as.numeric(max_leave_one_year_naive_delta_change),
    pre_epoch_max_single_year_share = as.numeric(pre_epoch_max_single_year_share),
    post_epoch_max_single_year_share = as.numeric(post_epoch_max_single_year_share)
  )

model <- safe_load(file.path(FIGURE_DATA_DIR, "selected_country", "prn_interpretability_model.tsv"), "interpretability model") %>%
  mutate(across(c(model_auc, cv_model_auc, standardized_coefficient), as.numeric))

dr <- safe_load(file.path(FIGURE_DATA_DIR, "selected_country", "selected_country_dr_missingness_summary.tsv"), "DR missingness") %>%
  mutate(across(starts_with("delta_"), as.numeric))

tipping <- safe_load(file.path(FIGURE_DATA_DIR, "selected_country", "selected_country_missingness_tipping_summary.tsv"), "tipping summary") %>%
  mutate(
    baseline_delta_interpretable_only = as.numeric(baseline_delta_interpretable_only),
    sign_change_observed_within_grid = as_logical_flag(sign_change_observed_within_grid)
  )

bootstrap <- safe_load(file.path(FIGURE_DATA_DIR, "study_dependence", "selected_country_block_bootstrap.tsv"), "block bootstrap") %>%
  mutate(
    observed_delta_naive_prevalence = as.numeric(observed_delta_naive_prevalence),
    bootstrap_delta_naive_prevalence_lower_95 = as.numeric(bootstrap_delta_naive_prevalence_lower_95),
    bootstrap_delta_naive_prevalence_upper_95 = as.numeric(bootstrap_delta_naive_prevalence_upper_95),
    bootstrap_delta_naive_prevalence_median = as.numeric(bootstrap_delta_naive_prevalence_median)
  )

transport <- safe_load(file.path(FIGURE_DATA_DIR, "selected_country", "selected_country_read_linked_transportability_ledger.tsv"), "read-linked transportability") %>%
  filter(comparison_dimension == "country") %>%
  mutate(
    full_interpretable_fraction = as.numeric(full_interpretable_fraction),
    read_linked_interpretable_fraction = as.numeric(read_linked_interpretable_fraction),
    category = factor(category, levels = rev(ordered_countries(category)))
  )

country_order <- rev(c("USA", "NZL", "JPN", "AUS"))

pA <- year_sens %>%
  mutate(country_iso3 = factor(country_iso3, levels = country_order)) %>%
  ggplot(aes(baseline_naive_delta, country_iso3)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = FIGURE_MID_GREY) +
  geom_point(aes(fill = baseline_naive_sign), shape = 21, size = 2.7, colour = FIGURE_INK, stroke = 0.22) +
  geom_segment(aes(x = baseline_naive_delta - max_leave_one_year_naive_delta_change,
                   xend = baseline_naive_delta + max_leave_one_year_naive_delta_change,
                   yend = country_iso3), colour = FIGURE_GREY, linewidth = 0.65) +
  scale_x_continuous(labels = percent) +
  scale_fill_manual(values = c(increase = unname(npg_colors["red"]), decrease = unname(npg_colors["blue"])), guide = "none") +
  labs(x = "Contrast delta", y = NULL) +
  theme_nature()

pB <- year_sens %>%
  transmute(
    country_iso3 = factor(country_iso3, levels = country_order),
    `Pre epoch` = pre_epoch_max_single_year_share,
    `Post epoch` = post_epoch_max_single_year_share
  ) %>%
  pivot_longer(-country_iso3, names_to = "epoch", values_to = "share") %>%
  ggplot(aes(share, country_iso3, colour = epoch)) +
  geom_point(size = 2.4) +
  scale_x_continuous(labels = percent, limits = c(0, 0.6)) +
  scale_colour_manual(values = c("Pre epoch" = unname(npg_colors["grey"]), "Post epoch" = unname(npg_colors["red"])), name = NULL) +
  labs(x = "Largest year share", y = NULL) +
  theme_nature() +
  theme(legend.position = "bottom")

pC <- model %>%
  filter(row_type == "feature_coefficient", model_variant == "full_model") %>%
  mutate(feature = factor(feature, levels = feature[order(standardized_coefficient)])) %>%
  ggplot(aes(standardized_coefficient, feature, fill = standardized_coefficient > 0)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = FIGURE_MID_GREY) +
  scale_fill_manual(values = c(`TRUE` = unname(npg_colors["green"]), `FALSE` = unname(npg_colors["blue"])), guide = "none") +
  labs(x = "Coefficient", y = NULL) +
  theme_nature()

pD <- model %>%
  filter(row_type == "model_diagnostic", !is.na(cv_model_auc)) %>%
  mutate(model_variant = factor(model_variant, levels = model_variant[order(cv_model_auc)])) %>%
  ggplot(aes(cv_model_auc, model_variant)) +
  geom_point(size = 2.5, fill = npg_colors["green"], shape = 21, colour = FIGURE_INK, stroke = 0.2) +
  scale_x_continuous(limits = c(0.7, 1.0)) +
  labs(x = "CV AUC", y = NULL) +
  theme_nature()

dr_long <- dr %>%
  select(country_iso3, starts_with("delta_")) %>%
  pivot_longer(-country_iso3, names_to = "estimator", values_to = "delta") %>%
  filter(!is.na(delta)) %>%
  mutate(country_iso3 = factor(country_iso3, levels = country_order))

pE <- ggplot(dr_long, aes(delta, country_iso3)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = FIGURE_MID_GREY) +
  geom_point(alpha = 0.62, size = 1.8, colour = npg_colors["purple"]) +
  labs(x = "Estimator contrast", y = NULL) +
  theme_nature()

pF <- bootstrap %>%
  mutate(country_iso3 = factor(country_iso3, levels = country_order)) %>%
  ggplot(aes(bootstrap_delta_naive_prevalence_median, country_iso3)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = FIGURE_MID_GREY) +
  geom_segment(aes(x = bootstrap_delta_naive_prevalence_lower_95, xend = bootstrap_delta_naive_prevalence_upper_95, yend = country_iso3),
               colour = FIGURE_GREY, linewidth = 0.8) +
  geom_point(size = 2.5, fill = npg_colors["peach"], shape = 21, colour = FIGURE_INK, stroke = 0.22) +
  scale_x_continuous(labels = percent) +
  labs(x = "Bootstrap contrast", y = NULL) +
  theme_nature()

pG <- tipping %>%
  mutate(country_iso3 = factor(country_iso3, levels = country_order)) %>%
  ggplot(aes(baseline_delta_interpretable_only, country_iso3, fill = sign_change_observed_within_grid)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = FIGURE_MID_GREY) +
  geom_point(shape = 21, size = 2.7, colour = FIGURE_INK, stroke = 0.22) +
  scale_fill_manual(values = c(`TRUE` = unname(npg_colors["peach"]), `FALSE` = unname(npg_colors["green"])), name = "Grid sign change") +
  scale_x_continuous(labels = percent) +
  labs(x = "Missingness-bound contrast", y = NULL) +
  theme_nature() +
  theme(legend.position = "bottom")

pH <- transport %>%
  filter(category %in% c("USA", "NZL", "JPN", "AUS", "CHN", "GBR", "FRA")) %>%
  pivot_longer(c(full_interpretable_fraction, read_linked_interpretable_fraction), names_to = "frame", values_to = "fraction") %>%
  mutate(frame = recode(frame, full_interpretable_fraction = "Full interpretable", read_linked_interpretable_fraction = "Read linked")) %>%
  ggplot(aes(fraction, category, colour = frame)) +
  geom_point(size = 2.2) +
  scale_x_continuous(labels = percent) +
  scale_colour_manual(values = c("Full interpretable" = unname(npg_colors["grey"]), "Read linked" = unname(npg_colors["green"])), name = NULL) +
  labs(x = "Country composition", y = NULL) +
  theme_nature() +
  theme(legend.position = "bottom")

ed2 <- (pA | pB | pC | pD) / (pE | pF | pG | pH) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = FIGURE_TAG_SIZE, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot",
    plot.margin = margin(3, 3, 3, 3)
  )

save_ed_pdf(ed2, "Extended_Data_Fig_02_Year_Composition_Sensitivity.pdf", height = NC_MAX_HEIGHT)
save_ed_png(ed2, "Extended_Data_Fig_02_Year_Composition_Sensitivity.png", height = NC_MAX_HEIGHT)
