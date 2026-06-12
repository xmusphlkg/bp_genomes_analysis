#!/usr/bin/env Rscript
# Extended Data Fig. 10: readiness-conditional USA focal-country sidecar

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

readiness <- safe_load(file.path(FIGURE_DATA_DIR, "dynamic_identifiability_report.tsv"), "dynamic identifiability") %>%
  mutate(across(c(full_mechanistic_readiness, recovery_readiness, fit_readiness, event_study_readiness), as_logical_flag))

summary <- safe_load(file.path(FIGURE_DATA_DIR, "dynamic_transmission_advantage_summary.tsv"), "dynamic summary") %>%
  mutate(across(c(effect_estimate, ci_lower, ci_upper, effect_ratio, effect_ratio_ci_lower, effect_ratio_ci_upper), as.numeric))

fit <- safe_load(file.path(FIGURE_DATA_DIR, "dynamic_fit_summary.tsv"), "dynamic fit summary") %>%
  mutate(n_obs = as.numeric(n_obs))

model_input <- safe_load(file.path(FIGURE_DATA_DIR, "dynamic_model_input.tsv"), "dynamic model input") %>%
  mutate(
    date = as.Date(date),
    cases = as.numeric(cases),
    annual_ipw_prevalence = as.numeric(annual_ipw_prevalence)
  )

recovery <- safe_load(file.path(FIGURE_DATA_DIR, "focal_country_recovery_summary.tsv"), "focal-country recovery summary") %>%
  mutate(
    year = as.numeric(year),
    current_interpretable = as.numeric(current_interpretable),
    reconciled_interpretable = as.numeric(reconciled_interpretable),
    max_attainable_interpretable = as.numeric(max_attainable_interpretable)
  )

country_levels <- c("USA", "CHN", "JPN")

ready_long <- readiness %>%
  select(country_iso3, full_mechanistic_readiness, recovery_readiness, fit_readiness, event_study_readiness) %>%
  pivot_longer(-country_iso3, names_to = "criterion", values_to = "ready") %>%
  mutate(
    country_iso3 = factor(country_iso3, levels = rev(c("USA", "CHN", "JPN"))),
    criterion = recode(criterion,
      full_mechanistic_readiness = "Full mechanistic",
      recovery_readiness = "Recovery",
      fit_readiness = "Model fit",
      event_study_readiness = "Event-study"
    )
  )

pA <- ggplot(ready_long, aes(criterion, country_iso3)) +
  geom_tile(aes(fill = ready), colour = "white", linewidth = 0.32) +
  scale_fill_manual(values = c(`TRUE` = unname(npg_colors["green"]), `FALSE` = FIGURE_LIGHT_GREY), guide = "none") +
  theme_nature_matrix() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, face = "plain"), axis.text.y = element_text(face = "bold"))

fit_status_raw <- fit %>%
  filter(country_iso3 %in% country_levels) %>%
  transmute(
    country_iso3,
    branch = case_when(
      analysis_branch == "half_mechanistic_main" ~ "Overlap NB",
      analysis_branch == "event_study" ~ "Event study",
      grepl("full_mechanistic", analysis_branch) ~ "Full mechanistic",
      TRUE ~ "Other"
    ),
    status = case_when(
      status == "ok" ~ "OK",
      status == "failed" ~ "Failed",
      status == "not_run" ~ "Not run",
      status == "skipped" ~ "Skipped",
      TRUE ~ "Other"
    )
  ) %>%
  filter(branch != "Other") %>%
  group_by(country_iso3, branch) %>%
  summarise(status = first(status), .groups = "drop")

pB <- expand_grid(country_iso3 = country_levels, branch = c("Overlap NB", "Event study", "Full mechanistic")) %>%
  left_join(fit_status_raw, by = c("country_iso3", "branch")) %>%
  mutate(
    country_iso3 = factor(country_iso3, levels = rev(country_levels)),
    branch = factor(branch, levels = c("Overlap NB", "Event study", "Full mechanistic")),
    status = factor(coalesce(status, "Not attempted"), levels = c("OK", "Failed", "Skipped", "Not run", "Not attempted"))
  ) %>%
  ggplot(aes(branch, country_iso3)) +
  geom_tile(aes(fill = status), colour = "white", linewidth = 0.32) +
  scale_fill_manual(
    values = c(
      "OK" = unname(npg_colors["green"]),
      "Failed" = unname(npg_colors["red"]),
      "Skipped" = unname(npg_colors["peach"]),
      "Not run" = FIGURE_MID_GREY,
      "Not attempted" = FIGURE_LIGHT_GREY
    ),
    name = "Status"
  ) +
  labs(x = NULL, y = NULL) +
  theme_nature_matrix() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, face = "plain"), axis.text.y = element_text(face = "bold"), legend.position = "bottom")

pC <- model_input %>%
  filter(country_iso3 %in% country_levels, !is.na(cases)) %>%
  mutate(country_iso3 = factor(country_iso3, levels = country_levels)) %>%
  ggplot(aes(date, cases, colour = country_iso3)) +
  geom_line(linewidth = 0.36, alpha = 0.9) +
  facet_wrap(~country_iso3, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = country_role_colors, guide = "none") +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = "Monthly reported cases") +
  theme_nature(base_size = FIGURE_BASE_SIZE) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), strip.text = element_text(face = "bold"))

pD <- recovery %>%
  filter(country_iso3 %in% country_levels) %>%
  select(country_iso3, year, current_interpretable, reconciled_interpretable, max_attainable_interpretable) %>%
  pivot_longer(
    c(current_interpretable, reconciled_interpretable, max_attainable_interpretable),
    names_to = "series",
    values_to = "interpretable"
  ) %>%
  mutate(
    country_iso3 = factor(country_iso3, levels = country_levels),
    series = recode(series,
      current_interpretable = "Current",
      reconciled_interpretable = "Reconciled",
      max_attainable_interpretable = "Attainable"
    )
  ) %>%
  ggplot(aes(year, interpretable, colour = series)) +
  geom_line(linewidth = 0.34, na.rm = TRUE) +
  geom_point(size = 0.85, na.rm = TRUE) +
  facet_wrap(~country_iso3, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c(Current = FIGURE_MID_GREY, Reconciled = unname(npg_colors["green"]), Attainable = unname(npg_colors["peach"])), name = NULL) +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = "Interpretable genomes") +
  theme_nature(base_size = FIGURE_BASE_SIZE) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom", strip.text = element_text(face = "bold"))

pE <- summary %>%
  filter(row_type == "event_study_control") %>%
  mutate(country_iso3 = factor(country_iso3, levels = rev(c("CHN", "JPN")))) %>%
  ggplot(aes(effect_estimate, country_iso3, fill = status)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.25, colour = FIGURE_MID_GREY) +
  geom_segment(aes(x = ci_lower, xend = ci_upper, yend = country_iso3), colour = FIGURE_GREY, linewidth = 0.8, na.rm = TRUE) +
  geom_point(shape = 21, size = 2.8, colour = FIGURE_INK, stroke = 0.22, na.rm = TRUE) +
  scale_fill_manual(values = c(ok = unname(npg_colors["peach"]), skipped = FIGURE_LIGHT_GREY), name = "Status") +
  scale_x_continuous(labels = number_format(accuracy = 1), expand = expansion(mult = c(0.04, 0.08))) +
  labs(x = "Post/pre descriptive ratio", y = NULL) +
  theme_nature() +
  theme(legend.position = "bottom")

ed10 <- (pA | pB) / pC / (pD | pE) +
  plot_layout(heights = c(0.78, 1.02, 0.95), widths = c(1.35, 0.9)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = FIGURE_TAG_SIZE, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot",
    plot.margin = margin(3, 3, 3, 3)
  )

save_ed_pdf(ed10, "Extended_Data_Fig_10_USA_Focal_Country.pdf", height = NC_MAX_HEIGHT)
save_ed_png(ed10, "Extended_Data_Fig_10_USA_Focal_Country.png", height = NC_MAX_HEIGHT)
