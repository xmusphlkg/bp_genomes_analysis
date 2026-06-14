#!/usr/bin/env Rscript
# Run all Extended Data figures from a stable script-directory working root

get_script_dir <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      if (!is.null(frame$ofile)) {
        return(normalizePath(frame$ofile))
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
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))))
  }
  normalizePath(".")
}

script_dir <- get_script_dir()
figure_root <- normalizePath(file.path(script_dir, ".."))
repo_root <- normalizePath(file.path(figure_root, "..", ".."))
Sys.setenv(
  PERTUSSIS_FIGURE_ROOT = figure_root,
  PERTUSSIS_REPO_ROOT = repo_root
)
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
scripts_root <- normalizePath(file.path(figure_root, "scripts"))
setwd(scripts_root)

source("lib/theme_nature.R")
source("lib/data_utils.R")

cat("\n========================================\n")
cat("  PERTUSSIS EXTENDED DATA PIPELINE\n")
cat("========================================\n")

required_pkgs <- c(
  "ggplot2", "cowplot", "scales", "RColorBrewer",
  "dplyr", "readr", "tidyr", "lubridate", "stringr",
  "maps", "ggrepel", "patchwork"
)
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' not installed. Run: install.packages('%s')", pkg, pkg))
  }
}

figures <- list(
  list(label = "Extended Data Fig. 1", script = "extended_data/ed01_sample_frame_reconciliation.R", name = "Sample-frame reconciliation"),
  list(label = "Extended Data Fig. 2", script = "extended_data/ed02_year_composition_sensitivity.R", name = "Sensitivity and missingness diagnostics"),
  list(label = "Extended Data Fig. 3", script = "extended_data/ed03_tree_representativeness.R", name = "Tree representativeness"),
  list(label = "Extended Data Fig. 4", script = "extended_data/ed04_asr_robustness.R", name = "ASR robustness"),
  list(label = "Extended Data Fig. 5", script = "extended_data/ed05_validation_package_anchors.R", name = "Validation matrix and package anchors"),
  list(label = "Extended Data Fig. 6", script = "extended_data/ed06_architecture_turnover.R", name = "Architecture turnover"),
  list(label = "Extended Data Fig. 7", script = "extended_data/ed07_origin_burden_bridge.R", name = "Origin burden bridge"),
  list(label = "Extended Data Fig. 8", script = "extended_data/ed08_prn_specificity_audit.R", name = "PRN-specificity audit"),
  list(label = "Extended Data Fig. 9", script = "extended_data/ed09_ecology_sidecar.R", name = "Ecology sidecar"),
  list(label = "Extended Data Fig. 10", script = "extended_data/ed10_usa_focal_country.R", name = "USA focal-country sidecar"),
  list(label = "Extended Data Fig. 11", script = "extended_data/ed11_identifiability_dag.R", name = "Archive identifiability DAG"),
  list(label = "Extended Data Fig. 12", script = "extended_data/ed12_country_year_missingness_audit.R", name = "Country-year missingness and block audit"),
  list(label = "Extended Data Fig. 13", script = "extended_data/ed13_fig01_context_panels.R", name = "Figure 1 context panels"),
  list(label = "Extended Data Fig. 14", script = "extended_data/ed14_structural_recurrence_compendium.R", name = "Structural recurrence compendium"),
  list(label = "Extended Data Fig. 15", script = "extended_data/ed15_validation_sensitivity_compendium.R", name = "Validation and sensitivity compendium")
)

results <- vector("list", length(figures))
for (i in seq_along(figures)) {
  fig <- figures[[i]]
  cat(sprintf("\n========== %s: %s ==========\n", fig$label, fig$name))
  tryCatch(
    {
      source(fig$script, local = new.env(parent = globalenv()))
      results[[i]] <- list(status = "SUCCESS", name = fig$name)
    },
    error = function(e) {
      cat(sprintf("ERROR in %s (%s): %s\n", fig$label, fig$name, e$message))
      results[[i]] <<- list(status = "FAILED", name = fig$name, error = e$message)
    }
  )
}

cat("\n\n========================================\n")
cat("  EXECUTION SUMMARY\n")
cat("========================================\n")
for (result in results) {
  symbol <- ifelse(result$status == "SUCCESS", "[OK]", "[FAIL]")
  cat(sprintf("  %s %s\n", symbol, result$name))
  if (!is.null(result$error)) {
    cat(sprintf("    Error: %s\n", result$error))
  }
}
success_count <- sum(vapply(results, function(result) identical(result$status, "SUCCESS"), logical(1)))
cat(sprintf("\n%d/%d Extended Data figures completed.\n", success_count, length(figures)))
