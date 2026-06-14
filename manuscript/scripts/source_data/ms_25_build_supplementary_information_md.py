#!/usr/bin/env python3
"""Build the Communications Biology Supplementary Information markdown.

Evidence ledgers that are needed for interpretation are rendered as readable
markdown tables in the Supplementary Information. Wide ledgers are reduced to
reader-facing columns for the SI text. Figure-source or diagnostic TSVs are
tracked in a disposition manifest rather than repeated in the SI text.
"""

from __future__ import annotations

import csv
import re
import shutil
from pathlib import Path
from typing import Dict, Iterable, List

ROOT = Path(__file__).resolve().parents[3]
MANUSCRIPT_DIR = ROOT / "manuscript"
SUPPLEMENTARY_DIR = MANUSCRIPT_DIR / "supplementary"
AUDIT_LEDGER_DIR = MANUSCRIPT_DIR / "submission_data" / "audit_ledgers" / "supplementary_table_sources"
TEXT_DIR = MANUSCRIPT_DIR / "text"
TEMPLATE_DIR = TEXT_DIR / "templates"
TABLE_TEMPLATE = TEMPLATE_DIR / "supplementary_table_temp.md"
FIGURE_TEMPLATE = TEMPLATE_DIR / "supplementary_fig_temp.md"
FINAL_OUT = TEXT_DIR / "commsbio_supplementary_information.md"
DISPOSITION_OUT = TEXT_DIR / "commsbio_supplementary_table_disposition.tsv"
PDF_TABLE_COLUMN_LIMIT = 14

OPTIMIZED_COLUMN_SUBSETS: Dict[str, List[str]] = {
    "Supplementary_Table_3_independent_origins.tsv": [
        "origin_id",
        "n_tips_total",
        "n_tips_disrupted",
        "n_countries",
        "first_year",
        "last_year",
        "major_mlst_st",
        "dominant_prn_mechanism",
        "branch_support",
        "origin_support_score",
    ],
    "Supplementary_Table_5_prn_Event_Definitions.tsv": [
        "prn_event_id",
        "mechanism_call",
        "event_subcategory",
        "sample_count",
        "country_count",
        "year_min",
        "year_max",
        "insertion_subject_gap_bp",
        "is_element_name",
        "call_confidence",
        "validation_level",
        "validation_status",
    ],
    "Supplementary_Table_6_ASR_Scenario_Registry.tsv": [
        "scenario_id",
        "scenario_class",
        "analysis_frame",
        "rooting_mode",
        "disrupted_tip_count",
        "fitch_origin_events",
        "pastml_origin_events",
        "mk_origin_count_mean",
        "mk_origin_count_lower_95",
        "mk_origin_count_upper_95",
        "largest_disrupted_clade_share",
        "max_disrupted_tips_per_origin",
        "rejects_one_global_clone_fitch",
        "rejects_one_global_clone_mk95",
    ],
    "Supplementary_Table_10_Event_Class_Phenotype_Evidence_Tiers.tsv": [
        "prn_event_id",
        "mechanism_call",
        "event_subcategory",
        "sample_count",
        "country_count",
        "year_min",
        "year_max",
        "validation_level",
        "phenotype_evidence_tier",
        "phenotype_inference",
        "caution_note",
    ],
    "Supplementary_Table_8_Junction_Confidence_Matrix.tsv": [
        "prn_event_id",
        "event_label",
        "mechanism_call",
        "event_subcategory",
        "sample_count",
        "country_count",
        "year_min",
        "year_max",
        "breakpoint_coordinate_basis",
        "breakpoint_left",
        "breakpoint_right",
        "representative_tsd_direct_repeats",
        "confidence_tier",
        "junction_interpretation",
    ],
    "Supplementary_Table_9_IS481_Target_Site_Accessibility.tsv": [
        "locus",
        "locus_label",
        "analysis_role",
        "structural_match_class",
        "locus_length_bp",
        "exact_ACTAGG_or_reverse_complement_count",
        "exact_ACTAGG_or_reverse_complement_per_kb",
        "hamming_distance_le1_to_ACTAGG_or_reverse_complement_count",
        "hamming_distance_le1_to_ACTAGG_or_reverse_complement_per_kb",
        "nearest_reference_IS481_distance_bp",
        "n_reference_IS481_features_within_10kb",
        "n_reference_IS481_features_within_50kb",
        "distance_from_observed_breakpoint_to_nearest_exact_target_bp",
        "interpretation_note",
    ],
    "Supplementary_Table_11_Recurrent_Event_Lineage_Country_Year_Anchors.tsv": [
        "prn_event_id",
        "event_subcategory",
        "mechanism_call",
        "sample_count",
        "country_count",
        "year_min",
        "year_max",
        "n_mlst_st",
        "top_mlst_st",
        "n_ptxP_labels",
        "top_ptxP_labels",
        "top_country_year_cells",
        "anchor_interpretation",
    ],
    "Supplementary_Table_12_Event_Specific_Acquisition_Packages.tsv": [
        "prn_event_id",
        "mechanism_call",
        "sample_count",
        "country_count",
        "n_country_year_cells",
        "n_mlst_st",
        "top_mlst_st",
        "top_country_year_cells",
        "acquisition_package_count",
        "non_singleton_package_count",
        "largest_package_disrupted_tips",
        "validation_level",
        "representative_tsd_direct_repeats",
        "event_specific_interpretation",
    ],
}

PUBLISHED_TABLES: List[Dict[str, str]] = [
    {
        "pub_no": "1",
        "source_file": "Supplementary_Table_1_Dataset_Composition.tsv",
        "title": "Genome collection and analysis-set composition",
        "rationale": "compact cohort and analysis-frame definition needed for reader orientation",
    },
    {
        "pub_no": "2",
        "source_file": "Supplementary_Table_2_prn_mechanism_classification.tsv",
        "title": "Mechanistic classes of prn disruption",
        "rationale": "definition table required to interpret event-class terminology",
    },
    {
        "pub_no": "3",
        "source_file": "Supplementary_Table_3_independent_origins.tsv",
        "title": "Phylogenetic evidence for minimum tree-level prn-disruption acquisition packages",
        "rationale": "compact evidence ledger for the central tree-level recurrence claim",
    },
    {
        "pub_no": "4",
        "source_file": "Supplementary_Table_4_Cohort_Flow_and_Tree_Selection.tsv",
        "title": "Cohort-flow filters and tree-selection criteria",
        "rationale": "compact audit trail linking raw genomes to the analysis tree",
        "legacy_source_file": "Supplementary_Table_7_Cohort_Flow_and_Tree_Selection.tsv",
    },
    {
        "pub_no": "5",
        "source_file": "Supplementary_Table_5_prn_Event_Definitions.tsv",
        "title": "Operational definitions for recurrent prn events",
        "rationale": "definition table required to reproduce event classification rules",
        "legacy_source_file": "Supplementary_Table_9_prn_Event_Definitions.tsv",
    },
    {
        "pub_no": "6",
        "source_file": "Supplementary_Table_6_ASR_Scenario_Registry.tsv",
        "title": "Ancestral-state reconstruction scenarios",
        "rationale": "compact sensitivity registry bounding ancestral-state reconstruction choices",
        "legacy_source_file": "Supplementary_Table_44_ASR_Scenario_Registry.tsv",
    },
    {
        "pub_no": "7",
        "source_file": "Supplementary_Table_7_Published_Overlap_Concordance.tsv",
        "title": "Concordance with published prn-status annotations",
        "rationale": "external concordance evidence that is not fully conveyed in figures",
        "legacy_source_file": "Supplementary_Table_55_Published_Overlap_Concordance.tsv",
    },
    {
        "pub_no": "8",
        "source_file": "Supplementary_Table_8_Junction_Confidence_Matrix.tsv",
        "title": "Junction-level support for recurrent prn-disruption events",
        "rationale": "technical confidence table separating read-supported and inferred events",
        "legacy_source_file": "Supplementary_Table_58_Junction_Confidence_Matrix.tsv",
    },
    {
        "pub_no": "9",
        "source_file": "Supplementary_Table_9_IS481_Target_Site_Accessibility.tsv",
        "title": "Comparator-locus context and IS481 target-site accessibility",
        "rationale": "compact sequence-context evidence for comparator matching and IS481-associated event opportunities",
        "legacy_source_file": "Supplementary_Table_59_IS481_Target_Site_Accessibility.tsv",
    },
    {
        "pub_no": "10",
        "source_file": "Supplementary_Table_10_Event_Class_Phenotype_Evidence_Tiers.tsv",
        "title": "Phenotype-evidence tiers for recurrent prn event classes",
        "rationale": "compact biology bridge while keeping causal language bounded",
        "legacy_source_file": "Supplementary_Table_62_Event_Class_Phenotype_Evidence_Tiers.tsv",
    },
    {
        "pub_no": "11",
        "source_file": "Supplementary_Table_11_Recurrent_Event_Lineage_Country_Year_Anchors.tsv",
        "title": "Lineage, country and year anchors for recurrent prn events",
        "rationale": "compact anchor table for recurrent event interpretation",
        "legacy_source_file": "Supplementary_Table_64_Recurrent_Event_Lineage_Country_Year_Anchors.tsv",
    },
    {
        "pub_no": "12",
        "source_file": "Supplementary_Table_12_Event_Specific_Acquisition_Packages.tsv",
        "title": "Event-specific minimum tree-level acquisition package support for recurrent prn events",
        "rationale": "compact event-specific recurrence summary supporting the revised Fig. 3d",
        "legacy_source_file": "Supplementary_Table_65_Event_Specific_Acquisition_Packages.tsv",
    },
]

COLUMN_RENAMES: Dict[str, str] = {
    "acquisition_package_count": "No. acquisition packages",
    "acquisition_package_ids": "Acquisition package IDs",
    "analysis_frame": "Analysis frame",
    "analysis_role": "Analysis role",
    "anchor_interpretation": "Anchor interpretation",
    "ancestral_state": "Ancestral state",
    "audit_section": "Audit section",
    "bp_category": "Breakpoint category",
    "branch_support": "Branch support",
    "breakpoint_coordinate_basis": "Breakpoint coordinate basis",
    "breakpoint_left": "Left breakpoint",
    "breakpoint_right": "Right breakpoint",
    "call_confidence": "Call confidence",
    "caution_note": "Caution note",
    "clade_id": "Clade ID",
    "collapse_or_weighting_rule": "Collapse or weighting rule",
    "confidence_tier": "Confidence tier",
    "concordance_fraction": "Concordance fraction",
    "comparison_or_risk": "Comparison or risk",
    "country_count": "No. countries",
    "country_iso3": "Country (ISO3)",
    "definition_limitations": "Definition limitations",
    "delta_from_previous": "Delta from previous stage",
    "descendant_state": "Descendant state",
    "distance_from_observed_breakpoint_to_nearest_exact_target_bp": "Distance to nearest exact target site (bp)",
    "dominant_event_collapsed_burden": "Dominant event collapsed burden",
    "dominant_event_genome_count": "Dominant event genome count",
    "dominant_event_hard_anchor": "Dominant event hard anchor",
    "dominant_event_id": "Dominant event ID",
    "dominant_event_label": "Dominant event label",
    "dominant_event_longread_exemplar": "Dominant event long-read exemplar",
    "dominant_event_rank_or_tie": "Dominant event rank or tie",
    "dominant_event_share": "Dominant event share",
    "dominant_event_supporting_read_or_public_longread": "Dominant event read or public long-read support",
    "dominant_event_supporting_validation_rows": "Dominant event validation rows",
    "dominant_event_validation_level": "Dominant event validation level",
    "dominant_prn_event_id": "Dominant prn event ID",
    "dominant_prn_mechanism": "Dominant prn mechanism",
    "earliest_year": "Earliest year",
    "event_count": "No. events",
    "event_definition_rule": "Event definition rule",
    "event_label": "Event label",
    "event_specific_interpretation": "Event-specific interpretation",
    "event_subcategory": "Event subcategory",
    "evidence_alignment": "Evidence alignment",
    "evidence_flags": "Evidence flags",
    "evidence_layer": "Evidence layer",
    "evidence_type": "Evidence type",
    "exact_ACTAGG_or_reverse_complement_count": "Exact ACTAGG or reverse-complement sites",
    "exact_ACTAGG_or_reverse_complement_per_kb": "Exact ACTAGG or reverse-complement sites per kb",
    "example_assembly_accession": "Example assembly accession",
    "example_contig_id": "Example contig ID",
    "example_gap_end": "Example gap end",
    "example_gap_start": "Example gap start",
    "example_sample_id_canonical": "Example sample ID",
    "example_sequencing_tech": "Example sequencing technology",
    "exemplar_replacement_applied": "Exemplar replacement applied",
    "exemplar_selection_rule": "Exemplar selection rule",
    "fitch_origin_events": "Fitch acquisition packages",
    "first_year": "First year",
    "followup_class": "Follow-up class",
    "hamming_distance_le1_to_ACTAGG_or_reverse_complement_count": "Sites within one mismatch of ACTAGG or reverse complement",
    "hamming_distance_le1_to_ACTAGG_or_reverse_complement_per_kb": "Sites within one mismatch per kb",
    "hsp_min_pident": "Minimum HSP identity (%)",
    "hit_orientation": "Hit orientation",
    "hit_support_tier": "Hit support tier",
    "inference_method": "Inference method",
    "insertion_subject_gap_bp": "Insertion gap size (bp)",
    "interpretation": "Interpretation",
    "interpretation_note": "Interpretation note",
    "is_definitive_disrupted": "Definitive disrupted call",
    "is_element_name": "Insertion sequence element",
    "is_interpretable": "Interpretable call",
    "is_support_profile": "IS support profile",
    "is_uncertain_fragmented": "Uncertain or fragmented call",
    "junction_interpretation": "Junction interpretation",
    "largest_disrupted_clade_share": "Largest disrupted-clade share",
    "last_year": "Last year",
    "latest_year": "Latest year",
    "largest_package_disrupted_tips": "Largest package disrupted tips",
    "locus": "Locus",
    "locus_label": "Locus label",
    "locus_length_bp": "Locus length (bp)",
    "locus_qcov_threshold": "Locus coverage threshold (%)",
    "longread_exemplar": "Long-read exemplar",
    "major_23s_status": "Major 23S status",
    "major_background_label": "Major background label",
    "major_fhaB2400_5550_label": "Major fhaB2400-5550 label",
    "major_fim3_label": "Major fim3 label",
    "major_lineage": "Major lineage",
    "major_lineage_source": "Major lineage source",
    "main_result": "Main result",
    "major_mlst_st": "Major MLST sequence type",
    "major_ptxP_label": "Major ptxP label",
    "max_disrupted_tips_per_origin": "Maximum disrupted tips per package",
    "max_total_clipped_reads": "Maximum clipped reads",
    "median_package_disrupted_tips": "Median package disrupted tips",
    "mechanism_call": "Mechanism call",
    "metric_name": "Metric",
    "mk_origin_count_lower_95": "Mk package lower 95% CI",
    "mk_origin_count_mean": "Mk package mean",
    "mk_origin_count_upper_95": "Mk package upper 95% CI",
    "n_background_profiles": "No. background profiles",
    "n_compared_rows": "No. compared records",
    "n_concordant": "No. concordant records",
    "n_countries": "No. countries",
    "n_countries_observed": "No. observed countries",
    "n_country_year_cells": "No. country-year cells",
    "n_fim3_labels": "No. fim3 labels",
    "n_genomes": "No. genomes",
    "n_insufficient": "No. insufficient calls",
    "n_is481": "No. IS481-disrupted calls",
    "n_mlst_st": "No. MLST sequence types",
    "n_missing_year": "No. records missing year",
    "n_overlap_rows": "No. overlapping records",
    "n_other_disruption": "No. other-disruption calls",
    "n_prn_disrupted": "No. prn-disrupted genomes",
    "n_prn_intact": "No. prn-intact genomes",
    "n_ptxP_labels": "No. ptxP labels",
    "n_published_sublineages": "No. published sublineages",
    "n_rearrangement": "No. rearrangement calls",
    "n_reference_IS481_features_within_10kb": "No. reference IS481 features within 10 kb",
    "n_reference_IS481_features_within_50kb": "No. reference IS481 features within 50 kb",
    "n_rows": "No. records",
    "n_structural_disrupted": "No. structural disrupted calls",
    "n_supporting_rows": "No. supporting records",
    "n_tips_disrupted": "No. disrupted tips",
    "n_tips_total": "No. tips",
    "nearest_reference_IS481_distance_bp": "Nearest reference IS481 distance (bp)",
    "non_singleton_package_count": "No. non-singleton packages",
    "notes": "Notes",
    "observed_breakpoint_flanking_sequence_25bp": "Observed breakpoint flanking sequence (25 bp)",
    "observed_gap1043_breakpoint_left": "Observed gap1043 left breakpoint",
    "observed_gap1043_breakpoint_right": "Observed gap1043 right breakpoint",
    "orientation": "Orientation",
    "origin_confidence_tier": "Origin confidence tier",
    "origin_confidence_tier_rank": "Origin confidence tier rank",
    "origin_evidence_class": "Origin evidence class",
    "origin_id": "Origin ID",
    "origin_n_disrupted_tips": "No. disrupted tips in origin",
    "origin_package_hard_anchor": "Origin package hard anchor",
    "origin_support_score": "Origin support score",
    "pastml_origin_events": "PastML acquisition packages",
    "phenotype_evidence_tier": "Phenotype evidence tier",
    "phenotype_inference": "Phenotype inference",
    "phylo_tree_id": "Phylogenetic tree ID",
    "primary_use_in_manuscript": "Primary manuscript use",
    "priority_origin_ids": "Priority origin IDs",
    "prn_disrupted_pct": "prn-disrupted genomes (%)",
    "prn_event_id": "prn event ID",
    "prn_mechanism_call": "prn mechanism call",
    "public_data_recovery_status": "Public-data recovery status",
    "rank_by_genome_burden": "Rank by genome burden",
    "read_locus_end": "Read-locus end",
    "read_locus_start": "Read-locus start",
    "read_reference_record": "Read reference record",
    "read_support_class": "Read-support class",
    "read_support_classes": "Read-support classes",
    "read_validation_status": "Read-validation status",
    "read_validation_statuses": "Read-validation statuses",
    "recovery_plan_status": "Recovery-plan status",
    "reference_end_1based": "Reference end (1-based)",
    "reference_start_1based": "Reference start (1-based)",
    "rejects_one_global_clone_fitch": "Rejects one-global-clone model (Fitch)",
    "rejects_one_global_clone_mk95": "Rejects one-global-clone model (Mk 95% CI)",
    "repo_prn_negative_and_published_prn_negative": "Repository prn-negative, published prn-negative",
    "repo_prn_negative_and_published_prn_positive": "Repository prn-negative, published prn-positive",
    "repo_prn_positive_and_published_prn_negative": "Repository prn-positive, published prn-negative",
    "repo_prn_positive_and_published_prn_positive": "Repository prn-positive, published prn-positive",
    "representative_assembly_accession": "Representative assembly accession",
    "representative_country_iso3": "Representative country (ISO3)",
    "representative_sample_id_canonical": "Representative sample ID",
    "representative_supporting_read_or_public_longread": "Representative read or public long-read support",
    "representative_tsd_direct_repeats": "Representative TSD direct repeats",
    "representative_validation_level": "Representative validation level",
    "representative_year": "Representative year",
    "rooting_frame": "Rooting frame",
    "rooting_mode": "Rooting mode",
    "rule_or_definition": "Rule or definition",
    "sample_count": "No. samples",
    "sample_fraction_all": "Fraction of all samples",
    "sample_fraction_within_mechanism": "Fraction within mechanism class",
    "sample_share_among_structurally_resolved": "Share among structurally resolved disruptions",
    "scenario_class": "Scenario class",
    "scenario_id": "Scenario ID",
    "scenario_source": "Scenario source",
    "sequencing_tech": "Sequencing technology",
    "sister_clade_id": "Sister clade ID",
    "singleton_package_count": "No. singleton packages",
    "source_file": "Source file",
    "stage_id": "Stage ID",
    "stage_name": "Stage name",
    "status_changed_vs_manuscript_fraction": "Status changed versus manuscript fraction",
    "strand": "Strand",
    "structural_match_class": "Structural match class",
    "summary_level": "Summary level",
    "supporting_external_context": "Supporting external context",
    "supporting_read_count": "No. supporting reads",
    "supporting_read_or_public_longread": "Read or public long-read support",
    "supporting_table": "Supporting table",
    "supporting_validation_rows": "No. supporting validation rows",
    "target_tsd_motif": "Target TSD motif",
    "tip_count": "No. tips",
    "top3_share_or_summary": "Top-three share or summary",
    "top_background_profiles": "Top background profiles",
    "top_country_year_cells": "Top country-year cells",
    "top_fim3_labels": "Top fim3 labels",
    "top_mlst_st": "Top MLST sequence types",
    "top_ptxP_labels": "Top ptxP labels",
    "top_published_sublineages": "Top published sublineages",
    "total_genomes": "No. genomes",
    "tree_representative_assembly_accession": "Tree representative assembly accession",
    "tree_representative_country_iso3": "Tree representative country (ISO3)",
    "tree_representative_sample_id_canonical": "Tree representative sample ID",
    "tree_representative_supporting_read_or_public_longread": "Tree representative read or public long-read support",
    "tree_representative_validation_level": "Tree representative validation level",
    "tree_representative_year": "Tree representative year",
    "tsd_direct_repeats": "TSD direct repeats",
    "tsd_or_flank_sequence_status": "TSD or flank sequence status",
    "tsd_supported_validation_rows": "TSD-supported validation rows",
    "validation_level": "Validation level",
    "validation_priority": "Validation priority",
    "validation_status": "Validation status",
    "year_max": "Latest year",
    "year_min": "Earliest year",
}

VALUE_RENAMES: Dict[str, str] = {
    "all_disrupted_descendants": "all disrupted descendants",
    "assembly_coordinate_only": "assembly-coordinate evidence only",
    "assembly_high_or_moderate": "high- or moderate-quality assembly evidence",
    "assembly_only": "assembly-level evidence only",
    "best_hit_is481": "best hit to IS481",
    "bp_insertion_like": "insertion-like breakpoint pattern",
    "bp_opposite_strand_or_complex": "opposite-strand or complex breakpoint pattern",
    "bp_within_contig": "within-contig breakpoint pattern",
    "branch_support_not_available": "branch support not available",
    "coding_disrupted_inversion_or_rearrangement": "coding-disrupting inversion or rearrangement",
    "coding_disrupted_is481": "coding-disrupting IS481 insertion",
    "coding_disrupted_other": "other coding-disrupting event",
    "combined_manuscript_manifest": "combined manuscript manifest",
    "composition_filtered": "composition-filtered",
    "composition_filtered_midpoint_rooted": "composition-filtered, midpoint-rooted",
    "composition_filtered_reference_rooted_primary": "composition-filtered, reference-rooted primary scenario",
    "composition_pruned": "composition-pruned",
    "composition_pruned_primary_asr_tree_nonreference_state_breakdown": "composition-pruned primary ASR tree non-reference state breakdown",
    "composition_pruned_primary_asr_tree_total_tips": "composition-pruned primary ASR tree total tips",
    "composition_pruned_quality_frame": "composition-pruned quality frame",
    "core_full_alignment_total": "core full-alignment total",
    "core_repeated_origin_package_support": "core repeated-origin package support",
    "country_balanced": "country-balanced",
    "country_x_st": "country and sequence-type collapse",
    "country_x_st_x_ptxp_fim_signature": "country, sequence-type, ptxP and fim signature collapse",
    "descendant_disrupted_tip_count": "descendant disrupted-tip count",
    "descendant_tip_count": "descendant tip count",
    "descendant_tips": "descendant tips",
    "disrupted_multi_hsp": "multi-HSP disrupted alignment",
    "dominant_1043bp_architecture": "dominant 1,043-bp architecture",
    "dominant_event_selection_scope": "dominant-event selection scope",
    "dominant_prn_event_sample_count": "dominant prn event sample count",
    "drop_largest_block_naive": "drop-largest-block naive sensitivity",
    "event_catalog_not_rooting_specific": "event catalogue not rooting-specific",
    "example_assembly": "example assembly",
    "excluded_pre_gubbins_missingness": "excluded before Gubbins because of missingness",
    "exploratory_ledger_not_core_package_support": "exploratory ledger, not core package support",
    "gap_sequence_extracted": "gap sequence extracted",
    "genome_intact_boundary": "genome-intact boundary",
    "include_in_snippy_ctg": "included in Snippy contig mode",
    "insertion_gap_ge_50bp": "insertion gap >=50 bp",
    "insertion_like": "insertion-like",
    "insufficient_data": "insufficient data",
    "is1002_ismapper_only": "IS1002 detected by ISMapper only",
    "is481_ismapper_only": "IS481 detected by ISMapper only",
    "is481_ismapper_panisa": "IS481 supported by ISMapper and PanISa",
    "is481_panisa_only": "IS481 detected by PanISa only",
    "large_surface_or_secreted_acellular_vaccine_adhesin": "large surface or secreted acellular-vaccine adhesin",
    "lineage_proxy_collapse": "lineage-proxy collapse",
    "mad_cli_not_available": "MAD command-line tool not available",
    "mad_rooting_feasibility": "MAD rooting feasibility check",
    "mapped_from_legacy_mirror_or_alt_accession": "mapped from a legacy mirror or alternate accession",
    "missing_fraction": "missing fraction",
    "moderate_is_hit": "moderate IS hit",
    "no_download_plan_match": "no matching download plan",
    "no_local_origin_call": "no local-origin call",
    "no_prn_is_signal_detected": "no prn IS signal detected",
    "no_prn_local_is_signal": "no local prn IS signal",
    "no_supported_is_hit": "no supported IS hit",
    "not_applicable": "not applicable",
    "not_recovered_current_public_data": "not recovered from available public data",
    "not_run": "not run",
    "origin_exemplar": "origin exemplar",
    "origin_exemplar_is_longread_anchored_while_same_event_is_read_backed_elsewhere": "origin exemplar is long-read anchored, while the same event is read-backed elsewhere",
    "origin_exemplar_matches_dominant_event_evidence_tier": "origin exemplar matches the dominant-event evidence tier",
    "origin_followup_exemplar": "origin follow-up exemplar",
    "origin_linked_event": "origin-linked event",
    "origin_package_collapse": "origin-package collapse",
    "origin_package_hard_anchor": "origin-package hard anchor",
    "origin_packages": "origin packages",
    "other_disruption": "other disruption",
    "other_or_unspecified": "other or unspecified",
    "pertactin_homologous_autotransporter": "pertactin-homologous autotransporter",
    "pertactin_target_locus": "pertactin target locus",
    "pre_gubbins_missingness": "pre-Gubbins missingness",
    "primary_asr": "primary ASR",
    "primary_asr_origin_packages": "primary ASR origin packages",
    "primary_structure_matched_pseudo_control": "primary structure-matched pseudo-control",
    "primary_target_locus": "primary target locus",
    "prn_insufficient_or_uncertain_subset": "prn-insufficient or uncertain subset",
    "prn_interpretable_subset": "prn-interpretable subset",
    "prn_mechanism_broad_concordance": "broad prn-mechanism concordance",
    "prn_status_concordance": "prn-status concordance",
    "public_assembly_universe": "public assembly universe",
    "public_longread_or_hybrid_anchor_present": "public long-read or hybrid-assembly anchor present",
    "public_longread_or_hybrid_assembly": "public long-read or hybrid assembly",
    "public_longread_or_hybrid_exemplar_present": "public long-read or hybrid-assembly exemplar present",
    "public_raw_read_assembly": "public raw-read assembly",
    "rank_1": "rank 1",
    "rank_1_after_dropping_largest_block": "rank 1 after dropping the largest block",
    "raw_read_linked_subset": "raw-read-linked subset",
    "raw_reads_available": "raw reads available",
    "raw_structurally_resolved_event_burden": "raw structurally resolved event burden",
    "read_backed_candidate": "read-backed candidate",
    "read_backed_or_candidate_available": "read-backed or candidate evidence available",
    "read_backed_supported": "read-backed support",
    "read_backed_targeted_validation": "read-backed targeted validation",
    "read_interval_without_tsd": "read interval without recovered target-site duplication",
    "read_or_longread_anchored_non_singleton": "read- or long-read-anchored non-singleton origin",
    "read_reference": "read-reference coordinate",
    "read_validation_table": "read-validation table",
    "read_validation_unresolved": "read validation unresolved",
    "rearrangement_breakpoint": "rearrangement breakpoint",
    "recoverable_paired_illumina": "recoverable paired-end Illumina reads",
    "reference_support_threshold": "reference support threshold",
    "representative_validation_scope": "representative validation scope",
    "resampling_country_balanced": "country-balanced resampling",
    "resampling_study_block_balanced": "study-block-balanced resampling",
    "resampling_time_balanced": "time-balanced resampling",
    "rooting_sensitivity": "rooting sensitivity",
    "rule_high": "high-confidence rule-based call",
    "rule_medium": "medium-confidence rule-based call",
    "sample_level": "sample-level",
    "secondary_structure_matched_pseudo_control": "secondary structure-matched pseudo-control",
    "secondary_vaccine_antigen_pseudo_control": "secondary vaccine-antigen pseudo-control",
    "single_origin_consistent": "single-origin consistent",
    "singleton_lower_confidence": "singleton, lower-confidence origin",
    "snippy_ctg": "Snippy contig mode",
    "snippy_ctg_plan": "Snippy contig-mode plan",
    "snippy_ctg_qc_pass": "Snippy contig-mode QC pass",
    "split_across_multiple_local_origins": "split across multiple local origins",
    "st_only": "sequence-type-only collapse",
    "strong_is_hit": "strong IS hit",
    "study_block_balanced": "study-block-balanced",
    "study_block_equalized": "study-block equalized",
    "study_block_stress_test": "study-block stress test",
    "study_weighted_top_event": "study-weighted top event",
    "support_threshold": "support-threshold",
    "support_threshold_sensitivity": "support-threshold sensitivity",
    "supported_candidate": "supported candidate",
    "supported_concordant": "supported concordant",
    "target_site_duplication_recovered": "target-site duplication recovered",
    "targeted_read_validation_completed_read_backed": "targeted read validation completed with read-backed support",
    "tier_1_read_backed_tsd_recovered": "Tier 1, read-backed with recovered target-site duplication",
    "tier_2_public_longread_or_hybrid": "Tier 2, public long-read or hybrid assembly",
    "tier_3_assembly_or_rule_supported": "Tier 3, assembly- or rule-supported",
    "tier_4_validation_unresolved": "Tier 4, validation unresolved",
    "time_balanced": "time-balanced",
    "tip_package": "tip package",
    "tip_states": "tip states",
    "top_collapsed_burden": "top collapsed burden",
    "top_collapsed_burden_events": "top collapsed-burden events",
    "top_collapsed_burden_tied_with": "top collapsed burden tied with",
    "top_three_events": "top three events",
    "transition_scan_on_fitch_parsimony_ml_tree": "transition scan on the Fitch-parsimony maximum-likelihood tree",
    "type_V_autotransporter_colonization_factor": "type V autotransporter colonization factor",
    "type_V_autotransporter_serine_protease": "type V autotransporter serine protease",
    "type_V_autotransporter_surface_virulence_factor": "type V autotransporter surface virulence factor",
    "type_V_autotransporter_with_reference_pseudogene_caveat": "type V autotransporter with reference-pseudogene caveat",
    "uncertain_fragmented_assembly": "uncertain fragmented assembly",
    "unpruned_asr_tree_nonreference_state_breakdown": "unpruned ASR tree non-reference state breakdown",
    "unpruned_asr_tree_total_tips": "unpruned ASR tree total tips",
    "unpruned_comparability_frame": "unpruned comparability frame",
    "unpruned_midpoint_rooted": "unpruned, midpoint-rooted",
    "unpruned_reference_rooted_comparability": "unpruned, reference-rooted comparability frame",
    "unpruned_support_ge_70": "unpruned, support >=70",
    "unpruned_support_ge_90": "unpruned, support >=90",
    "within_contig": "within-contig",
    "within_origin_concentration": "within-origin concentration",
    "workflow_ml_tree_composition_filtered": "composition-filtered maximum-likelihood tree",
}

SOURCE_FILE_LABELS: Dict[str, str] = {
    "bp_public_genome_manifest.tsv": "public genome manifest",
    "bp_public_genome_qc_manifest.tsv": "public genome QC manifest",
    "bp_prn_mechanism_calls.tsv": "prn mechanism-call table",
    "bp_prn_read_validation.tsv": "prn read-validation table",
    "pre_gubbins_missingness.tsv": "pre-Gubbins missingness table",
    "snippy_ctg_plan.tsv": "Snippy contig-mode plan",
    "tip_states.tsv": "ASR tip-state table",
}

FIGURE_SOURCE_OR_DUPLICATIVE = {
    4, 5, 6, 8, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23,
    25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
    41, 42, 43, 46, 47, 48, 49, 50, 51, 53, 54, 56, 60, 61, 63,
}


def extract_table_number(path: Path) -> int:
    match = re.search(r"Supplementary_Table_(\d+)_", path.name)
    if not match:
        raise ValueError(f"Could not parse supplementary table number from {path.name}")
    return int(match.group(1))


def read_tsv(path: Path) -> List[List[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.reader(handle, delimiter="\t"))


def count_data_rows(path: Path) -> int:
    rows = read_tsv(path)
    return max(len(rows) - 1, 0)


def count_columns(path: Path) -> int:
    rows = read_tsv(path)
    return len(rows[0]) if rows else 0


def title_from_table_file(path: Path) -> str:
    number = extract_table_number(path)
    prefix = f"Supplementary_Table_{number}_"
    stem = path.stem
    if stem.startswith(prefix):
        stem = stem[len(prefix):]
    title = humanize_snake_phrase(stem)
    return title[:1].upper() + title[1:]


def humanize_header(value: str) -> str:
    if value in COLUMN_RENAMES:
        return COLUMN_RENAMES[value]

    tokens = value.replace("_", " ").replace("  ", " ").strip().split(" ")
    fixed_tokens = []
    for token in tokens:
        if token.lower() in {"id", "ids"}:
            fixed_tokens.append(token.upper())
        elif token.lower() in {"iso3", "mlst", "asr", "tsd", "is481", "ci"}:
            fixed_tokens.append(token.upper())
        elif token.lower() == "prn":
            fixed_tokens.append("prn")
        else:
            fixed_tokens.append(token)
    if not fixed_tokens:
        return value
    return " ".join(fixed_tokens[:1]).capitalize() + (" " + " ".join(fixed_tokens[1:]) if len(fixed_tokens) > 1 else "")


def event_display_label(value: str) -> str | None:
    lower = value.lower()
    if lower == "prn_evt_intact":
        return "intact prn"
    if "fragmented" in lower or "insufficient" in lower:
        return "fragmented or insufficient prn call"
    if "cov58" in lower:
        return "rearrangement cov58"
    if "cov91" in lower:
        return "rearrangement cov91"
    if "cov94" in lower:
        return "rearrangement cov94"
    for gap in ("1045", "1044", "1043", "1042", "1041", "1040", "204", "54"):
        if f"gap{gap}" in lower:
            prefix = "IS481" if "is481" in lower else "other insertion-like"
            return f"{prefix} gap{gap}"
    return None


def humanize_snake_phrase(value: str) -> str:
    tokens = value.replace("_", " ").split()
    replacements = {
        "asr": "ASR",
        "bp": "breakpoint",
        "ci": "CI",
        "ctg": "contig",
        "fim": "fim",
        "fitch": "Fitch",
        "gubbins": "Gubbins",
        "hsp": "HSP",
        "illumina": "Illumina",
        "ipw": "IPW",
        "is481": "IS481",
        "mad": "MAD",
        "mk": "Mk",
        "ml": "maximum-likelihood",
        "mlst": "MLST",
        "panisa": "PanISa",
        "pastml": "PastML",
        "prn": "prn",
        "ptxp": "ptxP",
        "qc": "QC",
        "snippy": "Snippy",
        "st": "sequence type",
        "tsd": "target-site duplication",
    }
    return " ".join(replacements.get(token.lower(), token) for token in tokens)


def normalize_known_terms(value: str) -> str:
    value = re.sub(r"\bptxP_(\d+)\b", r"ptxP\1", value)
    value = re.sub(r"\bfim3_(\d+)\b", r"fim3-\1", value)
    value = re.sub(r"\borigin_(\d+)\b", r"Origin \1", value)
    value = re.sub(r"\bnode_(\d+)\b", r"Node \1", value)
    value = re.sub(r"\bstudy_block_balanced_replicate_(\d+)\b", r"study-block-balanced replicate \1", value)
    value = re.sub(r"\bcountry_balanced_replicate_(\d+)\b", r"country-balanced replicate \1", value)
    value = re.sub(r"\btime_balanced_replicate_(\d+)\b", r"time-balanced replicate \1", value)
    return value


def is_official_identifier(value: str) -> bool:
    return bool(re.match(r"^(GCA|GCF|SAMN|SRR|ERR|DRR|CP|NC)_?[A-Za-z0-9.:-]+$", value))


def humanize_source_path(value: str) -> str:
    basename = Path(value).name
    if basename in SOURCE_FILE_LABELS:
        return SOURCE_FILE_LABELS[basename]
    if basename.endswith(".tsv"):
        basename = basename[:-4]
    return humanize_snake_phrase(basename)


def humanize_token(token: str) -> str:
    token = token.strip()
    if not token:
        return token
    if token in VALUE_RENAMES:
        return VALUE_RENAMES[token]
    if is_official_identifier(token):
        return token
    event_label = event_display_label(token)
    if event_label:
        return event_label
    token = normalize_known_terms(token)
    if "_" in token and not is_official_identifier(token):
        return humanize_snake_phrase(token)
    return token


def humanize_value(value: str) -> str:
    value = value.strip()
    if value == "":
        return ""
    if value in {"True", "TRUE", "true"}:
        return "Yes"
    if value in {"False", "FALSE", "false"}:
        return "No"
    if "/" in value and " " not in value and not value.startswith(("http://", "https://")):
        return humanize_source_path(value)
    if is_official_identifier(value):
        return value
    if value in VALUE_RENAMES:
        return VALUE_RENAMES[value]
    if re.fullmatch(r"[A-Za-z0-9_.-]+", value):
        event_label = event_display_label(value)
        if event_label:
            return event_label
    if ">" in value and " " not in value:
        return " > ".join(humanize_value(part) for part in value.split(">"))
    if ";" in value:
        return "; ".join(humanize_value(part) for part in value.split(";"))
    if "|" in value:
        return " | ".join(humanize_value(part) for part in value.split("|"))
    if "==" in value and not value.startswith(("http://", "https://")):
        key, rest = value.split("==", 1)
        return f"{humanize_token(key)} = {humanize_value(rest)}"
    single_eq_match = re.match(r"^([A-Za-z0-9_ -]+)\s*=\s*([^=;,]+)$", value)
    if single_eq_match and not value.startswith(("http://", "https://")):
        key = single_eq_match.group(1)
        rest = single_eq_match.group(2)
        return f"{humanize_token(key)} = {humanize_value(rest)}"

    value = normalize_known_terms(value)
    value = re.sub(
        r"\b[A-Za-z][A-Za-z0-9]+(?:_[A-Za-z0-9]+)+\b",
        lambda match: humanize_token(match.group(0)),
        value,
    )
    return value


def escape_markdown_cell(value: str) -> str:
    value = value.replace("\r\n", "\n").replace("\r", "\n")
    value = value.replace("\n", "<br>")
    value = value.replace("|", r"\|")
    return value.strip()


def render_markdown_table(source_file: str) -> str:
    path = SUPPLEMENTARY_DIR / source_file
    rows = read_tsv(path)
    if not rows:
        return "_No rows in source TSV._"

    raw_header = rows[0]
    optimized_columns = OPTIMIZED_COLUMN_SUBSETS.get(source_file)
    if optimized_columns:
        missing = [column for column in optimized_columns if column not in raw_header]
        if missing:
            raise ValueError(f"Optimized column subset for {source_file} has missing columns: {missing}")
        keep_indices = [raw_header.index(column) for column in optimized_columns]
        rows = [
            [row[index] if index < len(row) else "" for index in keep_indices]
            for row in rows
        ]
        raw_header = rows[0]

    header = [escape_markdown_cell(humanize_header(cell)) for cell in raw_header]
    if len(header) > PDF_TABLE_COLUMN_LIMIT:
        raise ValueError(
            f"{source_file} has {len(header)} columns after optimization; "
            f"limit is {PDF_TABLE_COLUMN_LIMIT}"
        )

    body_rows = [
        [
            escape_markdown_cell(humanize_value(cell))
            for cell in row
        ]
        for row in rows[1:]
    ]
    width = len(header)

    normalized_rows: List[List[str]] = []
    for row in body_rows:
        if len(row) < width:
            row = row + [""] * (width - len(row))
        elif len(row) > width:
            row = row[:width]
        normalized_rows.append(row)

    lines = [
        "| " + " | ".join(header) + " |",
        "| " + " | ".join(["---"] * width) + " |",
    ]
    for row in normalized_rows:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def official_table_filenames() -> set[str]:
    return {table["source_file"] for table in PUBLISHED_TABLES}


def source_candidates_for_table(table: Dict[str, str]) -> list[Path]:
    source_file = table["source_file"]
    legacy_file = table.get("legacy_source_file", source_file)
    candidates = [
        SUPPLEMENTARY_DIR / source_file,
        SUPPLEMENTARY_DIR / legacy_file,
        AUDIT_LEDGER_DIR / source_file,
        AUDIT_LEDGER_DIR / legacy_file,
    ]
    return candidates


def copy_submission_table(source: Path, target: Path) -> None:
    """Copy a table into the official submission namespace.

    Historical audit ledgers may have trailing blank TSV fields. The official
    machine-readable tables use explicit NA cells so Git whitespace checks stay
    clean without changing the archived source ledgers.
    """

    rows = read_tsv(source)
    if not rows:
        target.write_text("", encoding="utf-8")
        return

    width = len(rows[0])
    normalized_rows: list[list[str]] = []
    for row in rows:
        row = row[:width] + [""] * max(width - len(row), 0)
        normalized_rows.append([cell if cell != "" else "NA" for cell in row])

    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(normalized_rows)


def prepare_curated_table_namespace() -> None:
    """Keep only official submission tables at the supplementary root.

    Upstream sidecar scripts still emit many machine-readable audit ledgers
    using historical Supplementary_Table_* names. For submission packaging,
    those wide ledgers are retained under submission_data/audit_ledgers and
    the supplementary root is repopulated with the compact official tables.
    """

    AUDIT_LEDGER_DIR.mkdir(parents=True, exist_ok=True)
    official_files = official_table_filenames()

    for path in sorted(SUPPLEMENTARY_DIR.glob("Supplementary_Table_*.tsv"), key=extract_table_number):
        if path.name in official_files:
            continue
        target = AUDIT_LEDGER_DIR / path.name
        if target.exists():
            target.unlink()
        shutil.move(str(path), str(target))

    for table in PUBLISHED_TABLES:
        target = SUPPLEMENTARY_DIR / table["source_file"]
        candidates = [path for path in source_candidates_for_table(table) if path.exists()]
        if not candidates:
            raise FileNotFoundError(
                f"No source found for official {table['source_file']} "
                f"(legacy={table.get('legacy_source_file', table['source_file'])})"
            )
        source = candidates[0]
        if source.resolve() != target.resolve():
            copy_submission_table(source, target)
        else:
            copy_submission_table(source, target.with_suffix(".tmp.tsv"))
            target.with_suffix(".tmp.tsv").replace(target)


def fill_table_template() -> str:
    text = TABLE_TEMPLATE.read_text(encoding="utf-8")
    for table in PUBLISHED_TABLES:
        placeholder = "{{SUPPLEMENTARY_TABLE_" + table["pub_no"] + "_BODY}}"
        if placeholder not in text:
            raise ValueError(f"Missing placeholder in table template: {placeholder}")
        source_path = SUPPLEMENTARY_DIR / table["source_file"]
        if not source_path.exists():
            raise FileNotFoundError(source_path)
        text = text.replace(placeholder, render_markdown_table(table["source_file"]))

    leftovers = re.findall(r"\{\{SUPPLEMENTARY_TABLE_\d+_BODY\}\}", text)
    if leftovers:
        raise ValueError("Unfilled supplementary table placeholders: " + ", ".join(leftovers))
    return text.rstrip() + "\n"


def figure_markdown_without_old_table_note() -> str:
    text = FIGURE_TEMPLATE.read_text(encoding="utf-8").rstrip()
    marker = "\n## Supplementary Tables and Data Files"
    if marker in text:
        text = text.split(marker, 1)[0].rstrip()
    return text + "\n"


def disposition_for_unpublished(table_number: int, row_count: int) -> str:
    if table_number in FIGURE_SOURCE_OR_DUPLICATIVE:
        return "figure-source, sensitivity, or secondary audit ledger already summarized in Supplementary Figures or main-text analyses; full TSV retained for reproducibility"
    return "machine-readable model, sensitivity, or secondary audit ledger retained for reproducibility; not repeated because it is not a reader-facing summary table"


def write_disposition_manifest() -> None:
    selected = {table["source_file"]: table for table in PUBLISHED_TABLES}
    source_files = sorted(SUPPLEMENTARY_DIR.glob("Supplementary_Table_*.tsv"), key=extract_table_number)

    fieldnames = [
        "source_table_number",
        "source_file",
        "source_data_rows",
        "source_data_columns",
        "publication_disposition",
        "si_pdf_presentation",
        "citation_label",
        "published_si_table_number",
        "published_si_title",
        "rationale",
    ]
    with DISPOSITION_OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for path in source_files:
            table_number = extract_table_number(path)
            row_count = count_data_rows(path)
            col_count = count_columns(path)
            if path.name in selected:
                table = selected[path.name]
                si_pdf_presentation = "optimized_markdown_table" if path.name in OPTIMIZED_COLUMN_SUBSETS else "rendered_markdown_table"
                row = {
                    "source_table_number": table_number,
                    "source_file": path.name,
                    "source_data_rows": row_count,
                    "source_data_columns": col_count,
                    "publication_disposition": "supplementary_table_tsv_with_si_pdf_entry",
                    "si_pdf_presentation": si_pdf_presentation,
                    "citation_label": f"Supplementary Table {table_number}",
                    "published_si_table_number": table_number,
                    "published_si_title": table["title"],
                    "rationale": table["rationale"],
                }
            else:
                row = {
                    "source_table_number": table_number,
                    "source_file": path.name,
                    "source_data_rows": row_count,
                    "source_data_columns": col_count,
                    "publication_disposition": "supplementary_table_tsv_only",
                    "si_pdf_presentation": "not_repeated_in_si_pdf",
                    "citation_label": f"Supplementary Table {table_number}",
                    "published_si_table_number": table_number,
                    "published_si_title": title_from_table_file(path),
                    "rationale": disposition_for_unpublished(table_number, row_count),
                }
            writer.writerow(row)


def assert_unique_sources(tables: Iterable[Dict[str, str]]) -> None:
    seen = set()
    for table in tables:
        source_file = table["source_file"]
        if source_file in seen:
            raise ValueError(f"Duplicate published source file: {source_file}")
        seen.add(source_file)


def main() -> None:
    assert_unique_sources(PUBLISHED_TABLES)
    prepare_curated_table_namespace()
    tables_md = fill_table_template()

    final_md = (
        figure_markdown_without_old_table_note().rstrip()
        + '\n\n<div style="page-break-after: always;"></div>\n\n'
        + tables_md.rstrip()
        + "\n"
    )
    FINAL_OUT.write_text(final_md, encoding="utf-8")
    write_disposition_manifest()

    print(f"Wrote {FINAL_OUT.relative_to(ROOT)}")
    print(f"Wrote {DISPOSITION_OUT.relative_to(ROOT)}")
    print(f"Published SI markdown tables: {len(PUBLISHED_TABLES)}")


if __name__ == "__main__":
    main()
