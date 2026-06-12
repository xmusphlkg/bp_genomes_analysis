#!/usr/bin/env Rscript
# Run binary stochastic mapping for registered ASR tree frames.

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(readr)
  library(dplyr)
  library(purrr)
  library(parallel)
})

`%||%` <- function(x, y) {
  if (length(x) == 0 || is.null(x) || is.na(x)) y else x
}

script_dir <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      if (!is.null(frame$ofile)) return(normalizePath(frame$ofile))
      ""
    },
    character(1)
  )
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0) return(dirname(tail(frame_files, 1)))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))))
  normalizePath(".")
}

repo_root <- normalizePath(file.path(script_dir(), "..", "..", ".."))
figure_data_dir <- file.path(repo_root, "manuscript", "figure_data")
supp_dir <- file.path(repo_root, "manuscript", "supplementary")

manifest_path <- file.path(repo_root, "outputs", "workflow", "manifest", "manifest.tsv")
scenario_registry_path <- file.path(figure_data_dir, "asr_scenario_registry.tsv")

as_bool <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "y", "t")
}

is_truthy <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "y", "t")
}

count_intact_to_disrupted <- function(sim_tree) {
  sum(vapply(
    sim_tree$maps,
    function(edge_map) {
      states <- names(edge_map)
      if (length(states) < 2) return(0L)
      sum(head(states, -1) == "intact" & tail(states, -1) == "disrupted")
    },
    integer(1)
  ))
}

state_for_tree <- function(tree, manifest) {
  state_lookup <- manifest %>%
    transmute(
      assembly_accession,
      state = case_when(
        prn_interpretable_bool & prn_disrupted_bool ~ "disrupted",
        prn_interpretable_bool & !prn_disrupted_bool ~ "intact",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(state), assembly_accession != "") %>%
    distinct(assembly_accession, .keep_all = TRUE)
  states <- setNames(state_lookup$state, state_lookup$assembly_accession)
  x <- states[tree$tip.label]
  names(x) <- tree$tip.label
  if ("Reference" %in% tree$tip.label) x["Reference"] <- "intact"
  x
}

make_fixed_er_q <- function(rate = 1.0) {
  q <- matrix(
    c(-rate, rate, rate, -rate),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("intact", "disrupted"), c("intact", "disrupted"))
  )
  q
}

split_nsim <- function(nsim, cores) {
  n_chunks <- max(1L, min(as.integer(nsim), as.integer(cores)))
  base <- floor(nsim / n_chunks)
  remainder <- nsim %% n_chunks
  chunks <- rep(base, n_chunks)
  if (remainder > 0) chunks[seq_len(remainder)] <- chunks[seq_len(remainder)] + 1L
  chunks[chunks > 0]
}

run_simmap_chunks <- function(tree, x, nsim, cores, scenario_seed, simmap_model, simmap_q, root_prior) {
  chunks <- split_nsim(nsim, cores)
  worker <- function(chunk_index) {
    chunk_nsim <- chunks[[chunk_index]]
    set.seed(scenario_seed + chunk_index)
    if (is.null(simmap_q)) {
      sims <- make.simmap(
        tree,
        x,
        model = simmap_model,
        nsim = chunk_nsim,
        pi = root_prior,
        message = FALSE
      )
    } else {
      sims <- make.simmap(
        tree,
        x,
        model = simmap_model,
        Q = simmap_q,
        nsim = chunk_nsim,
        pi = root_prior,
        message = FALSE
      )
    }
    sim_list <- if (inherits(sims, "multiSimmap")) sims else list(sims)
    vapply(sim_list, count_intact_to_disrupted, integer(1))
  }
  counts_by_chunk <- if (.Platform$OS.type == "unix" && cores > 1L && length(chunks) > 1L) {
    mclapply(seq_along(chunks), worker, mc.cores = cores, mc.preschedule = TRUE)
  } else {
    lapply(seq_along(chunks), worker)
  }
  unlist(counts_by_chunk, use.names = FALSE)
}

run_scenario <- function(scenario_id, scenario_class, tree_path, manifest, nsim = 50L, cores = 1L, simmap_model, simmap_q, root_prior, simmap_model_label) {
  cat(sprintf("  - %s: nsim=%d, cores=%d\n", scenario_id, nsim, cores))
  flush.console()
  if (!file.exists(tree_path)) {
    return(list(
      summary = tibble(
        scenario_id = scenario_id,
        scenario_class = scenario_class,
        tree_path = tree_path,
        nsim_requested = nsim,
        parallel_cores = cores,
        simmap_model = simmap_model_label,
        nsim_completed = 0L,
        tip_count = NA_integer_,
        disrupted_tip_count = NA_integer_,
        stochastic_origin_count_mean = NA_real_,
        stochastic_origin_count_median = NA_real_,
        stochastic_origin_count_lower_95 = NA_real_,
        stochastic_origin_count_upper_95 = NA_real_,
        stochastic_origin_count_min = NA_real_,
        stochastic_origin_count_max = NA_real_,
        proportion_origin_count_gt_1 = NA_real_,
        notes = "tree_missing"
      ),
      replicates = tibble()
    ))
  }
  tree <- read.tree(tree_path)
  x <- state_for_tree(tree, manifest)
  missing_tips <- names(x)[is.na(x)]
  if (length(missing_tips) > 0) {
    tree <- drop.tip(tree, missing_tips)
    x <- x[tree$tip.label]
  }
  if (!is.null(tree$edge.length)) {
    tree$edge.length[is.na(tree$edge.length) | tree$edge.length <= 0] <- 1e-8
  }
  x <- factor(x, levels = c("intact", "disrupted"))
  if (length(unique(x)) < 2) {
    return(list(
      summary = tibble(
        scenario_id = scenario_id,
        scenario_class = scenario_class,
        tree_path = tree_path,
        nsim_requested = nsim,
        parallel_cores = cores,
        simmap_model = simmap_model_label,
        nsim_completed = 0L,
        tip_count = length(x),
        disrupted_tip_count = sum(x == "disrupted"),
        stochastic_origin_count_mean = NA_real_,
        stochastic_origin_count_median = NA_real_,
        stochastic_origin_count_lower_95 = NA_real_,
        stochastic_origin_count_upper_95 = NA_real_,
        stochastic_origin_count_min = NA_real_,
        stochastic_origin_count_max = NA_real_,
        proportion_origin_count_gt_1 = NA_real_,
        notes = "single_state_after_missing_tip_drop"
      ),
      replicates = tibble()
    ))
  }
  scenario_seed <- 20260418 + sum(utf8ToInt(scenario_id))
  counts <- tryCatch(
    run_simmap_chunks(
      tree,
      x,
      nsim = nsim,
      cores = cores,
      scenario_seed = scenario_seed,
      simmap_model = simmap_model,
      simmap_q = simmap_q,
      root_prior = root_prior
    ),
    error = function(e) e
  )
  if (inherits(counts, "error")) {
    return(list(
      summary = tibble(
        scenario_id = scenario_id,
        scenario_class = scenario_class,
        tree_path = tree_path,
        nsim_requested = nsim,
        parallel_cores = cores,
        simmap_model = simmap_model_label,
        nsim_completed = 0L,
        tip_count = length(x),
        disrupted_tip_count = sum(x == "disrupted"),
        stochastic_origin_count_mean = NA_real_,
        stochastic_origin_count_median = NA_real_,
        stochastic_origin_count_lower_95 = NA_real_,
        stochastic_origin_count_upper_95 = NA_real_,
        stochastic_origin_count_min = NA_real_,
        stochastic_origin_count_max = NA_real_,
        proportion_origin_count_gt_1 = NA_real_,
        notes = paste("make.simmap_failed:", counts$message)
      ),
      replicates = tibble()
    ))
  }
  replicates <- tibble(
    scenario_id = scenario_id,
    scenario_class = scenario_class,
    simmap_model = simmap_model_label,
    stochastic_map_replicate = seq_along(counts),
    stochastic_origin_count = counts
  )
  summary <- tibble(
    scenario_id = scenario_id,
    scenario_class = scenario_class,
    tree_path = tree_path,
    nsim_requested = nsim,
    parallel_cores = cores,
    simmap_model = simmap_model_label,
    nsim_completed = length(counts),
    tip_count = length(x),
    disrupted_tip_count = sum(x == "disrupted"),
    stochastic_origin_count_mean = mean(counts),
    stochastic_origin_count_median = median(counts),
    stochastic_origin_count_lower_95 = unname(quantile(counts, 0.025)),
    stochastic_origin_count_upper_95 = unname(quantile(counts, 0.975)),
    stochastic_origin_count_min = min(counts),
    stochastic_origin_count_max = max(counts),
    proportion_origin_count_gt_1 = mean(counts > 1),
    notes = paste(
      "binary make.simmap",
      simmap_model_label,
      "model; missing/non-interpretable tips dropped before mapping; simulations split across parallel chunks"
    )
  )
  list(summary = summary, replicates = replicates)
}

build_resampling_scenarios <- function(resampling_root, nsim, require_registry = FALSE) {
  scheme_dirs <- sort(Sys.glob(file.path(resampling_root, "*_balanced")))
  rows <- list()
  missing_requirements <- character()
  registry_ids <- character()
  if (file.exists(scenario_registry_path)) {
    registry <- read_tsv(scenario_registry_path, show_col_types = FALSE)
    if ("scenario_id" %in% names(registry)) {
      registry_ids <- unique(as.character(registry$scenario_id))
      registry_ids <- registry_ids[nzchar(registry_ids)]
    }
  } else if (require_registry) {
    stop(
      paste(
        "Registered ASR scenario registry is missing:",
        scenario_registry_path,
        "Run ms_16_build_analysis_upgrade_sidecars.py before manuscript stochastic mapping."
      )
    )
  }
  for (scheme_dir in scheme_dirs) {
    replicate_dirs <- sort(Sys.glob(file.path(scheme_dir, "replicate_*")))
    for (replicate_dir in replicate_dirs) {
      tree_path <- file.path(replicate_dir, "rooted_ml_tree.reference_rooted.nwk")
      origin_path <- file.path(replicate_dir, "origin_events.tsv")
      tip_path <- file.path(replicate_dir, "tip_states.tsv")
      scenario_id <- paste0(basename(scheme_dir), "_", basename(replicate_dir))
      if (!file.exists(tree_path) || !file.exists(origin_path) || !file.exists(tip_path)) {
        missing_requirements <- c(missing_requirements, scenario_id)
        next
      }
      if (length(registry_ids) > 0 && !(scenario_id %in% registry_ids)) {
        missing_requirements <- c(missing_requirements, paste0(scenario_id, "::unregistered"))
        next
      }
      rows[[length(rows) + 1L]] <- tibble(
        scenario_id = scenario_id,
        scenario_class = paste0("resampling_", basename(scheme_dir)),
        tree_path = tree_path,
        nsim = nsim
      )
    }
  }
  if (length(missing_requirements) > 0 && require_registry) {
    stop(
      paste(
        "Resampling stochastic-mapping scenarios failed validation:",
        paste(missing_requirements, collapse = ", ")
      )
    )
  }
  if (length(rows) == 0) {
    return(tibble(
      scenario_id = character(),
      scenario_class = character(),
      tree_path = character(),
      nsim = integer()
    ))
  }
  bind_rows(rows)
}

resolve_output_paths <- function(full_run) {
  if (full_run) {
    return(list(
      summary = file.path(figure_data_dir, "asr_stochastic_mapping_summary.tsv"),
      replicates = file.path(figure_data_dir, "asr_stochastic_mapping_replicates.tsv"),
      supplementary = file.path(supp_dir, "Supplementary_Table_47_ASR_Stochastic_Mapping_Summary.tsv"),
      metadata = file.path(figure_data_dir, "asr_stochastic_mapping_metadata.tsv")
    ))
  }
  dev_dir <- file.path(figure_data_dir, "dev")
  list(
    summary = file.path(dev_dir, "asr_stochastic_mapping_summary.smoke.tsv"),
    replicates = file.path(dev_dir, "asr_stochastic_mapping_replicates.smoke.tsv"),
    supplementary = file.path(dev_dir, "Supplementary_Table_47_ASR_Stochastic_Mapping_Summary.smoke.tsv"),
    metadata = file.path(dev_dir, "asr_stochastic_mapping_metadata.smoke.tsv")
  )
}

detected_cores <- parallel::detectCores(logical = TRUE)
default_cores <- max(1L, detected_cores - 2L)
requested_cores <- suppressWarnings(as.integer(Sys.getenv("ASR_SIMMAP_CORES", default_cores)))
if (is.na(requested_cores) || requested_cores < 1L) requested_cores <- default_cores
simmap_cores <- min(requested_cores, detected_cores)
full_run <- is_truthy(Sys.getenv("ASR_SIMMAP_FULL", "1"))
smoke_run <- is_truthy(Sys.getenv("ASR_SIMMAP_SMOKE", "0"))
if (smoke_run) full_run <- FALSE
if (!full_run && !smoke_run) {
  stop(
    paste(
      "Reduced stochastic-mapping runs now require ASR_SIMMAP_SMOKE=1.",
      "Default manuscript-facing behavior is the full registered grid."
    )
  )
}
model_mode <- tolower(Sys.getenv("ASR_SIMMAP_MODEL", "empirical_ard"))
allow_fixed_er <- is_truthy(Sys.getenv("ASR_SIMMAP_ALLOW_FIXED_ER", "0"))
fixed_er_rate <- suppressWarnings(as.numeric(Sys.getenv("ASR_SIMMAP_FIXED_ER_RATE", "20000")))
if (is.na(fixed_er_rate) || fixed_er_rate <= 0) fixed_er_rate <- 1.0
if (model_mode %in% c("ard", "ml_ard", "empirical_ard")) {
  simmap_model <- "ARD"
  simmap_q <- NULL
  root_prior <- "equal"
  simmap_model_label <- "ARD_empirical_Q"
} else if (model_mode %in% c("er", "fixed_er")) {
  if (!allow_fixed_er) {
    stop(
      paste(
        "Fixed-ER stochastic mapping is now development-only.",
        "Set ASR_SIMMAP_ALLOW_FIXED_ER=1 to run it explicitly."
      )
    )
  }
  simmap_model <- "ER"
  simmap_q <- make_fixed_er_q(fixed_er_rate)
  root_prior <- "equal"
  simmap_model_label <- paste0("fixed_ER_Q_rate_", fixed_er_rate, "_root_equal")
} else {
  stop(paste("Unsupported ASR_SIMMAP_MODEL:", model_mode))
}

default_primary_nsim <- if (full_run) 1000L else 100L
default_secondary_nsim <- if (full_run) 500L else 50L
default_resampling_nsim <- if (full_run) 200L else 25L
primary_nsim <- suppressWarnings(as.integer(Sys.getenv("ASR_SIMMAP_PRIMARY_NSIM", default_primary_nsim)))
secondary_nsim <- suppressWarnings(as.integer(Sys.getenv("ASR_SIMMAP_SECONDARY_NSIM", default_secondary_nsim)))
resampling_nsim <- suppressWarnings(as.integer(Sys.getenv("ASR_SIMMAP_RESAMPLING_NSIM", default_resampling_nsim)))
if (is.na(primary_nsim) || primary_nsim < 1L) primary_nsim <- default_primary_nsim
if (is.na(secondary_nsim) || secondary_nsim < 1L) secondary_nsim <- default_secondary_nsim
if (is.na(resampling_nsim) || resampling_nsim < 1L) resampling_nsim <- default_resampling_nsim
output_paths <- resolve_output_paths(full_run)

cat(sprintf(
  "Running ASR stochastic mapping sidecar: detected_cores=%d, using_cores=%d, full_run=%s, smoke_run=%s\n",
  detected_cores,
  simmap_cores,
  full_run,
  smoke_run
))
cat(sprintf("Stochastic mapping model: %s\n", simmap_model_label))
flush.console()

manifest <- read_tsv(manifest_path, show_col_types = FALSE) %>%
  mutate(
    prn_interpretable_bool = as_bool(prn_interpretable),
    prn_disrupted_bool = as_bool(prn_disrupted)
  )

scenario_rows <- tibble(
  scenario_id = c(
    "primary_reference_rooted",
    "composition_filtered_reference_rooted",
    "support_ge_70_reference_rooted",
    "support_ge_90_reference_rooted",
    "composition_filtered_midpoint",
    "unpruned_midpoint"
  ),
  scenario_class = c(
    "primary",
    "composition_filtered",
    "support_threshold",
    "support_threshold",
    "rooting_sensitivity",
    "rooting_sensitivity"
  ),
  tree_path = file.path(
    repo_root,
    c(
      "outputs/workflow/asr/rooted_ml_tree.reference_rooted.nwk",
      "outputs/workflow/asr_sensitivity/composition_filtered/rooted_ml_tree.reference_rooted.nwk",
      "outputs/workflow/asr_sensitivity/support_70/rooted_ml_tree.reference_rooted.nwk",
      "outputs/workflow/asr_sensitivity/support_90/rooted_ml_tree.reference_rooted.nwk",
      "outputs/workflow/asr_rooting_sensitivity/composition_filtered_midpoint/rooted_ml_tree.midpoint_rooted.nwk",
      "outputs/workflow/asr_rooting_sensitivity/unpruned_midpoint/rooted_ml_tree.midpoint_rooted.nwk"
    )
  ),
  nsim = c(primary_nsim, rep(secondary_nsim, 5))
)

resampling_rows <- build_resampling_scenarios(
  file.path(repo_root, "outputs/workflow/asr_resampling"),
  nsim = resampling_nsim,
  require_registry = full_run
)
if (nrow(resampling_rows) > 0) {
  scenario_rows <- bind_rows(scenario_rows, resampling_rows)
}

results <- pmap(
  list(scenario_rows$scenario_id, scenario_rows$scenario_class, scenario_rows$tree_path, scenario_rows$nsim),
  ~ run_scenario(..1, ..2, ..3, manifest, ..4, simmap_cores, simmap_model, simmap_q, root_prior, simmap_model_label)
)

summary <- bind_rows(lapply(results, `[[`, "summary"))
replicates <- bind_rows(lapply(results, `[[`, "replicates"))
metadata <- tibble(
  run_timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  full_run = full_run,
  smoke_run = smoke_run,
  manuscript_supported = full_run && !smoke_run && model_mode %in% c("ard", "ml_ard", "empirical_ard"),
  simmap_model = simmap_model_label,
  simmap_model_mode = model_mode,
  fixed_er_rate = if (model_mode %in% c("er", "fixed_er")) fixed_er_rate else NA_real_,
  primary_nsim = primary_nsim,
  secondary_nsim = secondary_nsim,
  resampling_nsim = resampling_nsim,
  scenario_count = nrow(summary),
  output_summary_path = output_paths$summary,
  output_replicate_path = output_paths$replicates,
  notes = if (full_run) {
    "canonical_manuscript_facing_full_registered_run"
  } else {
    "development_smoke_run_written_to_dev_output_paths_only"
  }
)

dir.create(dirname(output_paths$summary), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_paths$supplementary), recursive = TRUE, showWarnings = FALSE)
write_tsv(summary, output_paths$summary)
write_tsv(summary, output_paths$supplementary)
write_tsv(replicates, output_paths$replicates)
write_tsv(metadata, output_paths$metadata)

cat(sprintf("Wrote %d stochastic-mapping scenario summaries and %d replicate rows\n", nrow(summary), nrow(replicates)))
