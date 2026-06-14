#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

get_script_dir <- function() {
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

locate_repo_root <- function(start_dir) {
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

script_dir <- get_script_dir()
repo_root <- locate_repo_root(script_dir)
supp_dir <- file.path(repo_root, "manuscript", "supplementary")
figure_data_dir <- file.path(repo_root, "manuscript", "figure_data")

read_tsv_quiet <- function(path) {
  read_tsv(path, show_col_types = FALSE, na = c("", "NA"))
}

first_nonempty <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

top_counts_string <- function(x, n = 3) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    return("")
  }
  tab <- sort(table(x), decreasing = TRUE)
  head_tab <- head(tab, n)
  paste(sprintf("%s (%d)", names(head_tab), as.integer(head_tab)), collapse = "; ")
}

event_defs <- read_tsv_quiet(file.path(supp_dir, "Supplementary_Table_9_prn_Event_Definitions.tsv")) %>%
  mutate(
    sample_count = suppressWarnings(as.integer(sample_count)),
    country_count = suppressWarnings(as.integer(country_count)),
    year_min = suppressWarnings(as.integer(year_min)),
    year_max = suppressWarnings(as.integer(year_max))
  )

mechanism_calls <- read_tsv_quiet(file.path(figure_data_dir, "fig02_prn_mechanism_calls.tsv")) %>%
  mutate(
    year = suppressWarnings(as.integer(year)),
    prn_interpretable = prn_mechanism_call %in% c(
      "intact",
      "coding_disrupted_is481",
      "coding_disrupted_inversion_or_rearrangement",
      "coding_disrupted_other"
    ),
    prn_disrupted = prn_mechanism_call %in% c(
      "coding_disrupted_is481",
      "coding_disrupted_inversion_or_rearrangement",
      "coding_disrupted_other"
    )
  )

published_overlap <- read_tsv_quiet(file.path(figure_data_dir, "published_overlap_annotation.tsv")) %>%
  mutate(
    year = suppressWarnings(as.integer(year)),
    prn_interpretable = if_else(is.na(prn_interpretable), FALSE, prn_interpretable),
    prn_disrupted = if_else(is.na(prn_disrupted), FALSE, prn_disrupted)
  )

phylo_context <- read_tsv_quiet(file.path(figure_data_dir, "figure3_global_phylogeny_context.tsv")) %>%
  filter(panel_id == "tip_metadata") %>%
  transmute(
    sample_id_canonical,
    phylo_mlst_st = na_if(as.character(mlst_st), "")
  )

study_blocks <- read_tsv_quiet(file.path(supp_dir, "Supplementary_Table_52_Study_Block_Assignment_and_Dominance.tsv")) %>%
  filter(row_type == "sample_block_assignment") %>%
  mutate(
    year = suppressWarnings(as.integer(year)),
    prn_interpretable = if_else(is.na(prn_interpretable), FALSE, prn_interpretable),
    prn_disrupted = if_else(is.na(prn_disrupted), FALSE, prn_disrupted)
  )

country_lookup <- published_overlap %>%
  transmute(
    country_iso3,
    country_name = na_if(as.character(published_country), "")
  ) %>%
  filter(!is.na(country_iso3), country_iso3 != "", !is.na(country_name)) %>%
  group_by(country_iso3) %>%
  summarise(country_name = first_nonempty(country_name), .groups = "drop")

disrupted_event_defs <- event_defs %>%
  filter(mechanism_call %in% c(
    "coding_disrupted_is481",
    "coding_disrupted_inversion_or_rearrangement",
    "coding_disrupted_other"
  )) %>%
  transmute(
    prn_event_id,
    mechanism_call,
    event_subcategory,
    sample_count,
    country_count,
    year_min,
    year_max,
    validation_level,
    evidence_type,
    tsd_or_flank_sequence_status,
    phenotype_evidence_tier = case_when(
      event_subcategory %in% c("IS481 insertion", "Inversion / rearrangement") &
        validation_level == "read_backed_supported" &
        tsd_or_flank_sequence_status == "target_site_duplication_recovered" ~ "Tier 1a junction-class phenotype bridge",
      event_subcategory %in% c("IS481 insertion", "Inversion / rearrangement") ~ "Tier 1b lesion-class bridge only",
      tsd_or_flank_sequence_status == "target_site_duplication_recovered" ~ "Tier 2 genome-disruption plausible",
      validation_level %in% c("read_backed_supported", "public_longread_or_hybrid_assembly") ~ "Tier 2 genome-disruption plausible",
      TRUE ~ "Tier 3 genome-only disruption"
    ),
    phenotype_inference = case_when(
      phenotype_evidence_tier == "Tier 1a junction-class phenotype bridge" ~ "Read-backed junction or event-class evidence places this architecture within lesion classes repeatedly linked to PRN non-expression, but the present archive still lacks genome-by-genome expression assays.",
      phenotype_evidence_tier == "Tier 1b lesion-class bridge only" ~ "The lesion class is externally linked to PRN non-expression, but this exact event class is not protein-tested here and has weaker junction-level support in the present archive summary.",
      phenotype_evidence_tier == "Tier 2 genome-disruption plausible" ~ "Coding disruption is structurally compatible with PRN loss and has strong genome evidence, but no exact event-class expression bridge was identified for this archive summary.",
      TRUE ~ "Genome architecture supports disruption, but phenotype inference remains indirect and event-specific protein evidence is currently absent."
    ),
    supporting_external_context = case_when(
      event_subcategory == "IS481 insertion" ~ "Japan 1990-2010, United States 1935-2012, Europe 1998-2015, Belgium 2000-2023",
      event_subcategory == "Inversion / rearrangement" ~ "Europe 1998-2015, Belgium 2000-2023",
      TRUE ~ "External phenotype bridge supports the broader disrupted-locus interpretation rather than this exact event label"
    ),
    caution_note = case_when(
      phenotype_evidence_tier == "Tier 1a junction-class phenotype bridge" ~ "Treat as junction-class bridged, not genome-by-genome protein proven.",
      phenotype_evidence_tier == "Tier 1b lesion-class bridge only" ~ "Treat as lesion-class bridged only; prioritize read-backed or expression validation.",
      TRUE ~ "Retain as genome-defined disruption until paired PRN expression assays are available."
    )
  )

intact_summary <- mechanism_calls %>%
  filter(prn_mechanism_call == "intact") %>%
  summarise(
    sample_count = n(),
    country_count = n_distinct(country_iso3[country_iso3 != ""]),
    year_min = suppressWarnings(min(year, na.rm = TRUE)),
    year_max = suppressWarnings(max(year, na.rm = TRUE))
  )

phenotype_tiers <- bind_rows(
  disrupted_event_defs,
  tibble(
    prn_event_id = "prn_evt_intact",
    mechanism_call = "intact",
    event_subcategory = "Intact locus",
    sample_count = intact_summary$sample_count,
    country_count = intact_summary$country_count,
    year_min = intact_summary$year_min,
    year_max = intact_summary$year_max,
    validation_level = "not_applicable",
    evidence_type = "genome_intact_boundary",
    tsd_or_flank_sequence_status = "not_applicable",
    phenotype_evidence_tier = "Boundary only, not expression-proven",
    phenotype_inference = "A genome-intact prn locus does not prove PRN production at the protein level in this archive.",
    supporting_external_context = "Australia 2008-2012 and Belgium 2000-2023 report PRN-negative isolates without an obvious genome-only prn lesion.",
    caution_note = "Do not treat intact-locus calls as direct PRN-positive phenotype assignments."
  )
) %>%
  arrange(
    factor(
      phenotype_evidence_tier,
      levels = c(
        "Tier 1a junction-class phenotype bridge",
        "Tier 1b lesion-class bridge only",
        "Tier 2 genome-disruption plausible",
        "Tier 3 genome-only disruption",
        "Boundary only, not expression-proven"
      )
    ),
    desc(sample_count),
    prn_event_id
  )

write_tsv(
  phenotype_tiers,
  file.path(supp_dir, "Supplementary_Table_62_Event_Class_Phenotype_Evidence_Tiers.tsv")
)

country_year_block_audit <- study_blocks %>%
  filter(!is.na(year), !is.na(country_iso3), country_iso3 != "") %>%
  group_by(country_iso3, year) %>%
  group_modify(~{
    base_tab <- sort(table(.x$base_block_id), decreasing = TRUE)
    sub_tab <- sort(table(.x$subblock_id), decreasing = TRUE)
    dominant_base <- names(base_tab)[1]
    dominant_sub <- names(sub_tab)[1]
    dominant_base_rows <- .x %>% filter(base_block_id == dominant_base)
    tibble(
      country_name = first_nonempty(.x$country_name),
      n_total_records = nrow(.x),
      n_prn_interpretable = sum(.x$prn_interpretable, na.rm = TRUE),
      n_prn_noninterpretable = sum(!.x$prn_interpretable, na.rm = TRUE),
      interpretability_fraction = sum(.x$prn_interpretable, na.rm = TRUE) / nrow(.x),
      n_prn_disrupted = sum(.x$prn_disrupted & .x$prn_interpretable, na.rm = TRUE),
      disrupted_fraction_among_interpretable = if_else(
        sum(.x$prn_interpretable, na.rm = TRUE) > 0,
        sum(.x$prn_disrupted & .x$prn_interpretable, na.rm = TRUE) / sum(.x$prn_interpretable, na.rm = TRUE),
        NA_real_
      ),
      n_base_blocks = n_distinct(.x$base_block_id),
      n_subblocks = n_distinct(.x$subblock_id),
      dominant_base_block_id = dominant_base,
      dominant_base_block_n_total = as.integer(base_tab[1]),
      dominant_base_block_share_total = as.integer(base_tab[1]) / nrow(.x),
      dominant_base_block_share_interpretable = if_else(
        sum(.x$prn_interpretable, na.rm = TRUE) > 0,
        sum(dominant_base_rows$prn_interpretable, na.rm = TRUE) / sum(.x$prn_interpretable, na.rm = TRUE),
        NA_real_
      ),
      dominant_subblock_id = dominant_sub,
      dominant_subblock_n_total = as.integer(sub_tab[1]),
      dominant_subblock_share_total = as.integer(sub_tab[1]) / nrow(.x),
      audit_flag = case_when(
        nrow(.x) < 3 ~ "sparse_cell",
        as.integer(base_tab[1]) / nrow(.x) >= 0.8 ~ "base_block_dominated",
        sum(!.x$prn_interpretable, na.rm = TRUE) / nrow(.x) >= 0.5 ~ "low_interpretability",
        TRUE ~ "distributed_or_moderate"
      )
    )
  }) %>%
  ungroup() %>%
  left_join(country_lookup, by = "country_iso3", suffix = c("", "_lookup")) %>%
  mutate(country_name = coalesce(country_name, country_name_lookup)) %>%
  select(-country_name_lookup) %>%
  arrange(desc(n_total_records), country_iso3, year)

write_tsv(
  country_year_block_audit,
  file.path(supp_dir, "Supplementary_Table_63_Country_Year_Interpretability_Study_Block_Audit.tsv")
)

write_tsv(
  country_year_block_audit,
  file.path(figure_data_dir, "epidemiology_revision_country_year_audit.tsv")
)

event_anchor_input <- published_overlap %>%
  filter(
    prn_mechanism_call %in% c(
      "coding_disrupted_is481",
      "coding_disrupted_inversion_or_rearrangement",
      "coding_disrupted_other"
    ),
    prn_disrupted,
    !is.na(prn_event_id),
    prn_event_id != "",
    prn_event_id != "prn_evt_intact"
  ) %>%
  left_join(
    event_defs %>% select(prn_event_id, event_subcategory, mechanism_call, sample_count, country_count, year_min, year_max),
    by = "prn_event_id"
  ) %>%
  left_join(phylo_context, by = "sample_id_canonical") %>%
  mutate(
    mlst_st_clean = na_if(as.character(mlst_st), ""),
    ptxP_label = na_if(as.character(ptxP_label), ""),
    fim3_label = na_if(as.character(fim3_label), ""),
    background_display_label = na_if(as.character(background_display_label), ""),
    published_sublineage_label = na_if(as.character(published_sublineage_label), ""),
    country_year_cell = if_else(!is.na(country_iso3) & !is.na(year), sprintf("%s:%d", country_iso3, year), NA_character_)
  )

event_anchor_summary <- event_anchor_input %>%
  group_by(prn_event_id, event_subcategory, mechanism_call, sample_count, country_count, year_min, year_max) %>%
  summarise(
    n_genomes = n(),
    n_countries_observed = n_distinct(country_iso3[!is.na(country_iso3) & country_iso3 != ""]),
    n_country_year_cells = n_distinct(country_year_cell[!is.na(country_year_cell)]),
    n_mlst_st = n_distinct(mlst_st_clean[!is.na(mlst_st_clean)]),
    top_mlst_st = top_counts_string(mlst_st_clean),
    n_ptxP_labels = n_distinct(ptxP_label[!is.na(ptxP_label)]),
    top_ptxP_labels = top_counts_string(ptxP_label),
    n_fim3_labels = n_distinct(fim3_label[!is.na(fim3_label)]),
    top_fim3_labels = top_counts_string(fim3_label),
    n_background_profiles = n_distinct(background_display_label[!is.na(background_display_label)]),
    top_background_profiles = top_counts_string(background_display_label),
    n_published_sublineages = n_distinct(published_sublineage_label[!is.na(published_sublineage_label)]),
    top_published_sublineages = top_counts_string(published_sublineage_label),
    top_country_year_cells = top_counts_string(country_year_cell),
    anchor_interpretation = case_when(
      n_country_year_cells >= 3 | n_background_profiles >= 2 ~ "Multiple lineage or sampling anchors; not consistent with one sampled country-year clone.",
      n_countries_observed >= 2 ~ "Observed in more than one country, but additional phylogenetic anchoring remains sparse.",
      TRUE ~ "Sparse archive anchor; interpret with stronger dependence on local validation and block audits."
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(n_genomes), prn_event_id)

write_tsv(
  event_anchor_summary,
  file.path(supp_dir, "Supplementary_Table_64_Recurrent_Event_Lineage_Country_Year_Anchors.tsv")
)

write_tsv(
  event_anchor_summary,
  file.path(figure_data_dir, "epidemiology_revision_event_anchors.tsv")
)

message("Wrote Supplementary Tables 62-64 and figure-data sidecars for the epidemiology revision audit.")
