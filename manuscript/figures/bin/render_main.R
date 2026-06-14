#!/usr/bin/env Rscript
# Run all active manuscript figures from a stable script-directory working root

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

cat("\n========================================\n")
cat("  PERTUSSIS MANUSCRIPT FIGURE PIPELINE\n")
cat("========================================\n")

required_pkgs <- c(
  "ggplot2", "ggforce", "cowplot", "scales", "RColorBrewer",
  "dplyr", "readr", "tidyr", "lubridate", "stringr",
  "maps", "ggrepel", "patchwork", "ape", "ggtree"
)
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' not installed. Run: install.packages('%s')", pkg, pkg))
  }
}
cat(sprintf("All %d required packages are available.\n", length(required_pkgs)))

cat("\nLoading shared libraries...\n")
source("lib/theme_nature.R")
source("lib/data_utils.R")
cat("Shared libraries loaded.\n")

figures <- list(
  list(label = "Figure 1", script = "main/fig01_public_genome_atlas.R", name = "Public-Genome Atlas"),
  list(label = "Figure 2", script = "main/fig02_prn_structural_solution_space.R", name = "PRN Structural Solution Space"),
  list(label = "Figure 3", script = "main/fig03_repeated_origin_phylogeny.R", name = "Repeated-Acquisition Phylogeny"),
  list(label = "Figure 4", script = "main/fig04_country_programme_amplification.R", name = "Archive Context and Sampling Heterogeneity"),
  list(label = "Figure 5", script = "main/fig05_validation_synthesis.R", name = "Validation, Specificity and Phenotype Synthesis")
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
cat(sprintf("\n%d/%d figures completed.\n", success_count, length(figures)))
