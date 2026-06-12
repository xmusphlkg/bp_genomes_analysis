# Data Loading, Validation, and Transformation Utilities
# Redesigned 2026-04-03: strict validation, no synthetic data

library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(stringr)

# Active figure scripts should read only frozen figure/submission assets or
# small tracked state, never the live module output roots.
get_support_script_dir <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      if (!is.null(frame$ofile)) {
        return(normalizePath(frame$ofile, winslash = "/", mustWork = FALSE))
      }
      ""
    },
    character(1)
  )
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0) {
    return(dirname(tail(frame_files, 1)))
  }
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)))
  }
  normalizePath(".", winslash = "/", mustWork = FALSE)
}

locate_repo_root <- function(start_dir = getwd()) {
  current <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  repeat {
    if (dir.exists(file.path(current, ".git")) &&
        dir.exists(file.path(current, "manuscript")) &&
        dir.exists(file.path(current, "workflow"))) {
      return(current)
    }
    parent <- normalizePath(file.path(current, ".."), winslash = "/", mustWork = FALSE)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }
  stop(sprintf("Could not locate repository root from: %s", start_dir))
}

get_repo_root <- function() {
  env_root <- Sys.getenv("PERTUSSIS_REPO_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = FALSE))
  }
  locate_repo_root(get_support_script_dir())
}

get_figure_root <- function() {
  env_root <- Sys.getenv("PERTUSSIS_FIGURE_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = FALSE))
  }
  normalizePath(file.path(get_repo_root(), "manuscript", "figures"), winslash = "/", mustWork = FALSE)
}

BASE_DIR <- get_repo_root()
FIGURE_ROOT <- get_figure_root()
FIGURE_DATA_DIR <- file.path(BASE_DIR, "manuscript", "figure_data")
SUBMISSION_DATA_DIR <- file.path(BASE_DIR, "manuscript", "submission_data")
STATE_MANIFEST_DIR <- file.path(BASE_DIR, "state", "manifest")

data_paths <- list(
  selected_country_figure_data_dir = file.path(FIGURE_DATA_DIR, "selected_country"),
  selected_country_docs_dir = file.path(BASE_DIR, "docs", "manuscript_rebuild", "selected_country"),
  selected_country_curation_dir = file.path(BASE_DIR, "manuscript", "curation", "selected_country"),
  prn_country_year = file.path(FIGURE_DATA_DIR, "fig01_prn_country_year_summary.tsv"),
  prn_mechanisms = file.path(FIGURE_DATA_DIR, "fig02_prn_mechanism_summary.tsv"),
  fig1_landscape = file.path(FIGURE_DATA_DIR, "figure1_data_landscape.tsv"),
  fig3_tree_segments = file.path(FIGURE_DATA_DIR, "figure3_workflow_tree_segments.tsv"),
  fig3_tree_nodes = file.path(FIGURE_DATA_DIR, "figure3_workflow_tree_nodes.tsv"),
  fig3_asr_summary = file.path(FIGURE_DATA_DIR, "figure3_workflow_asr_sensitivity.tsv"),
  fig3_asr_resampling = file.path(FIGURE_DATA_DIR, "figure3_workflow_asr_resampling.tsv"),
  fig4_lineage_collapse = file.path(FIGURE_DATA_DIR, "lineage_collapsed_event_table.tsv"),
  fig4_origin_collapse = file.path(FIGURE_DATA_DIR, "origin_collapsed_event_table.tsv"),
  fig4_validation_matrix = file.path(FIGURE_DATA_DIR, "architecture_origin_validation_matrix.tsv"),
  fig4_local_package_trees = file.path(FIGURE_DATA_DIR, "local_rooted_package_tree_summary.tsv"),
  fig4_dynamic_counterfactual = file.path(FIGURE_DATA_DIR, "dynamic_counterfactual_summary.tsv"),
  fig4_dynamic_ident = file.path(FIGURE_DATA_DIR, "dynamic_identifiability_report.tsv"),
  fig5_models = file.path(FIGURE_DATA_DIR, "figure5_association_model_panels.tsv"),
  fig5_leave_one_out = file.path(FIGURE_DATA_DIR, "figure5_leave_one_country_out_summary.tsv"),
  fig5_formulation = file.path(FIGURE_DATA_DIR, "figure5_formulation_coverage.tsv"),
  fig5_programme_class_summary = file.path(FIGURE_DATA_DIR, "supplementary_programme_class_summary.tsv"),
  fig5_detectability_summary = file.path(FIGURE_DATA_DIR, "prn_event_class_detectability.tsv"),
  fig5_detectability_detail = file.path(FIGURE_DATA_DIR, "prn_event_class_detectability_detail.tsv"),
  fig7_transmission_summary = file.path(FIGURE_DATA_DIR, "dynamic_transmission_advantage_summary.tsv"),
  fig7_transmission_predictions = file.path(FIGURE_DATA_DIR, "dynamic_transmission_advantage_predictions.tsv"),
  fig7_structural_concentration = file.path(FIGURE_DATA_DIR, "structural_event_concentration.tsv"),
  supplementary_asr_sensitivity = file.path(FIGURE_DATA_DIR, "fig04_asr_sensitivity.tsv"),
  supplementary_representativeness_audit = file.path(FIGURE_DATA_DIR, "supplementary_programme_representativeness_audit.tsv"),
  supplementary_country_influence = file.path(FIGURE_DATA_DIR, "supplementary_programme_country_influence.tsv"),
  supplementary_panel_eligibility = file.path(FIGURE_DATA_DIR, "supplementary_programme_panel_eligibility.tsv"),
  supplementary_precedence_conflicts = file.path(FIGURE_DATA_DIR, "supplementary_programme_precedence_conflicts.tsv"),
  rr_program_history_manifest = file.path(FIGURE_DATA_DIR, "selected_country", "country_program_history_manifest.tsv"),
  rr_selection_scorecard = file.path(FIGURE_DATA_DIR, "selected_country", "country_selection_scorecard.tsv"),
  rr_epoch_eligibility = file.path(FIGURE_DATA_DIR, "selected_country", "country_epoch_eligibility.tsv"),
  rr_epoch_prevalence = file.path(FIGURE_DATA_DIR, "selected_country", "country_epoch_prn_prevalence.tsv"),
  rr_epoch_bounds = file.path(FIGURE_DATA_DIR, "selected_country", "country_epoch_bounds.tsv"),
  rr_epoch_contrast = file.path(FIGURE_DATA_DIR, "selected_country", "country_epoch_contrast_summary.tsv"),
  rr_origin_amplification = file.path(FIGURE_DATA_DIR, "selected_country", "selected_country_origin_amplification.tsv"),
  rr_detection_amplification = file.path(FIGURE_DATA_DIR, "selected_country", "selected_country_detection_amplification.tsv"),
  rr_origin_package_summary = file.path(FIGURE_DATA_DIR, "selected_country", "selected_country_origin_package_summary.tsv"),
  rr_structure_reuse = file.path(FIGURE_DATA_DIR, "selected_country", "selected_country_structure_reuse.tsv"),
  rr_validation_matrix = file.path(FIGURE_DATA_DIR, "selected_country", "selected_country_validation_matrix.tsv"),
  rr_selection_synthesis = file.path(FIGURE_DATA_DIR, "selected_country", "cross_country_selection_synthesis.tsv"),
  rr_evidence_grid = file.path(FIGURE_DATA_DIR, "selected_country", "cross_country_evidence_grid.tsv"),
  rr_mechanism_bridge = file.path(FIGURE_DATA_DIR, "selected_country", "cross_country_mechanism_bridge.tsv"),
  rr_relative_year_plot_data = file.path(FIGURE_DATA_DIR, "selected_country", "selected_country_relative_year_plot_data.tsv"),
  manifest = file.path(STATE_MANIFEST_DIR, "manifest.tsv")
)

safe_load <- function(path, label = "file") {
  if (!file.exists(path)) stop(sprintf("[MISSING] %s not found: %s", label, path))
  df <- read_tsv(path, show_col_types = FALSE, na = c("", "NA"))
  message(sprintf("[OK] %s: %d rows x %d cols", label, nrow(df), ncol(df)))
  df
}

validate_prn_summary <- function(df) {
  df$year <- suppressWarnings(as.integer(df$year))
  df$n_genomes_total <- suppressWarnings(as.numeric(df$n_genomes_total))
  df$n_genomes_prn_interpretable <- suppressWarnings(as.numeric(df$n_genomes_prn_interpretable))
  df$n_prn_intact <- suppressWarnings(as.numeric(df$n_prn_intact))
  df$n_prn_disrupted <- suppressWarnings(as.numeric(df$n_prn_disrupted))
  df$frac_prn_disrupted <- suppressWarnings(as.numeric(df$frac_prn_disrupted))
  for (col in c("frac_prn_disrupted_ci_lower", "frac_prn_disrupted_ci_upper")) {
    if (!col %in% names(df)) {
      df[[col]] <- NA_real_
    }
  }
  required <- c("country_iso3", "year", "n_genomes_total", "n_genomes_prn_interpretable",
                "n_prn_intact", "n_prn_disrupted", "frac_prn_disrupted",
                "frac_prn_disrupted_ci_lower", "frac_prn_disrupted_ci_upper")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) stop(sprintf("PRN summary missing columns: %s", paste(missing, collapse = ", ")))
  year_values <- suppressWarnings(as.numeric(df$year))
  message(sprintf("[VALIDATED] PRN summary: %d countries, years %s-%s, %d interpretable genomes",
                  n_distinct(df$country_iso3),
                  format(min(year_values, na.rm = TRUE), trim = TRUE, scientific = FALSE),
                  format(max(year_values, na.rm = TRUE), trim = TRUE, scientific = FALSE),
                  as.integer(sum(df$n_genomes_prn_interpretable, na.rm=TRUE))))
  df
}

validate_re_trajectories <- function(df) {
  required <- c("country", "year", "time_index", "re_estimate", "re_ci_lower", "re_ci_upper",
                "raw_incidence", "smoothed_incidence", "quality_flag")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) stop(sprintf("RE trajectories missing columns: %s", paste(missing, collapse=", ")))
  flag_counts <- df %>% count(quality_flag) %>% arrange(desc(n))
  message("[RE QUALITY FLAGS]")
  for (i in seq_len(nrow(flag_counts))) {
    message(sprintf("  %s: %d (%.1f%%)", flag_counts$quality_flag[i], flag_counts$n[i],
                    100*flag_counts$n[i]/nrow(df)))
  }
  df_ok <- df %>% filter(quality_flag == "OK")
  if (nrow(df_ok) == 0) stop("NO RE trajectory data passes quality check.")
  extreme_re <- df_ok %>% filter(re_estimate > 10 | re_estimate < 0.01)
  if (nrow(extreme_re) > 0) {
    warning(sprintf("%d OK rows have extreme Re values", nrow(extreme_re)))
    df_ok <- df_ok %>% filter(re_estimate >= 0.01, re_estimate <= 10)
  }
  message(sprintf("[VALIDATED] RE trajectories: %d OK rows, %d countries",
                  nrow(df_ok), n_distinct(df_ok$country)))
  df_ok
}

validate_re_summary <- function(df) {
  required <- c("country", "year", "mean_re", "median_re", "num_weeks_valid")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) stop(sprintf("RE summary missing columns: %s", paste(missing, collapse=", ")))
  df_valid <- df %>% filter(!is.na(mean_re) & mean_re > 0 & num_weeks_valid > 0)
  empty_count <- nrow(df) - nrow(df_valid)
  if (empty_count > 0) message(sprintf("[FILTERED] Removed %d empty/invalid rows (%.1f%%)",
    empty_count, 100*empty_count/nrow(df)))
  message(sprintf("[VALIDATED] RE summary: %d country-year estimates, %d countries",
                  nrow(df_valid), n_distinct(df_valid$country)))
  df_valid
}

validate_association_models <- function(df) {
  required <- c("panel_id", "model_id", "estimate_term", "effect_estimate",
                "ci_lower", "ci_upper", "p_value", "n_obs", "n_countries")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) stop(sprintf("Association models missing columns: %s", paste(missing, collapse=", ")))
  if ("standard_glm_warning_types" %in% names(df)) {
    warned <- df %>% filter(!is.na(standard_glm_warning_types) & standard_glm_warning_types != "")
    if (nrow(warned) > 0) message(sprintf("[WARNINGS] %d model terms have GLM warnings", nrow(warned)))
  }
  low_n <- df %>% filter(n_obs < 10)
  if (nrow(low_n) > 0) warning(sprintf("%d model terms have n_obs < 10", nrow(low_n)))
  message(sprintf("[VALIDATED] Association models: %d terms, %d countries, n_obs range [%d, %d]",
                  nrow(df), max(df$n_countries, na.rm=TRUE),
                  min(df$n_obs, na.rm=TRUE), max(df$n_obs, na.rm=TRUE)))
  df
}

aggregate_prn_by_year <- function(df) {
  df %>% filter(n_genomes_prn_interpretable > 0) %>%
    group_by(year) %>%
    summarise(n_genomes = sum(n_genomes_total, na.rm=TRUE),
              n_interpretable = sum(n_genomes_prn_interpretable, na.rm=TRUE),
              n_intact = sum(n_prn_intact, na.rm=TRUE),
              n_disrupted = sum(n_prn_disrupted, na.rm=TRUE),
              frac_disrupted = n_disrupted / n_interpretable, .groups="drop") %>%
    filter(n_interpretable >= 5)
}

aggregate_prn_by_country <- function(df, top_n = 15) {
  df %>% filter(n_genomes_prn_interpretable > 0) %>%
    group_by(country_iso3) %>%
    summarise(n_genomes = sum(n_genomes_total, na.rm=TRUE),
              n_interpretable = sum(n_genomes_prn_interpretable, na.rm=TRUE),
              n_intact = sum(n_prn_intact, na.rm=TRUE),
              n_disrupted = sum(n_prn_disrupted, na.rm=TRUE),
              frac_disrupted = n_disrupted / n_interpretable,
              year_min = min(year, na.rm=TRUE), year_max = max(year, na.rm=TRUE),
              .groups="drop") %>%
    filter(n_interpretable >= 5) %>% arrange(desc(n_interpretable)) %>% slice_head(n=top_n)
}
