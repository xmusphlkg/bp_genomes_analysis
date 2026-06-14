#!/usr/bin/env python3
"""Build a submission-facing source-data manifest for final figure panels."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "manuscript" / "submission_data" / "source_data" / "final_source_data_manifest.tsv"


def row(
    figure_label: str,
    panel_label: str,
    panel_title: str,
    final_placement: str,
    frozen_input_files: list[str],
) -> dict[str, str]:
    return {
        "figure_label": figure_label,
        "panel_label": panel_label,
        "panel_title": panel_title,
        "final_placement": final_placement,
        "frozen_input_files": ";".join(frozen_input_files),
    }


ROWS = [
    row(
        "Figure 1",
        "A",
        "Public-genome sampling atlas",
        "main",
        ["manuscript/figure_data/fig01_prn_country_year_summary.tsv"],
    ),
    row(
        "Figure 1",
        "B",
        "Country-year genome density",
        "main",
        ["manuscript/figure_data/fig01_prn_country_year_summary.tsv"],
    ),
    row(
        "Figure 1",
        "C",
        "Recoverable-locus trajectory",
        "main",
        ["manuscript/figure_data/fig01_prn_country_year_summary.tsv"],
    ),
    row(
        "Figure 1",
        "D",
        "Country-readiness scatter",
        "main",
        ["manuscript/figure_data/selected_country/country_selection_scorecard.tsv"],
    ),
    row(
        "Supplementary Figure 13",
        "A-C",
        "Figure 1 context panels",
        "supplementary",
        [
            "manuscript/figure_data/fig01_prn_country_year_summary.tsv",
            "manuscript/supplementary/Supplementary_Table_4_Cohort_Flow_and_Tree_Selection.tsv",
        ],
    ),
    row(
        "Figure 2",
        "A",
        "Structural design-space schematic",
        "main",
        [
            "manuscript/figure_data/prn_event_evidence_manifest.tsv",
            "manuscript/figure_data/is481_target_site_accessibility.tsv",
        ],
    ),
    row(
        "Figure 2",
        "B",
        "Recurrent event catalogue",
        "main",
        [
            "manuscript/figure_data/selected_country/selected_country_structure_reuse.tsv",
            "manuscript/figure_data/prn_junction_confidence_matrix.tsv",
        ],
    ),
    row(
        "Figure 2",
        "C",
        "Country architecture reuse matrix",
        "main",
        ["manuscript/figure_data/selected_country/selected_country_structure_reuse.tsv"],
    ),
    row(
        "Figure 2",
        "D",
        "Structural-null concentration",
        "main",
        ["manuscript/figure_data/structural_event_concentration.tsv"],
    ),
    row(
        "Figure 2",
        "E",
        "Study-block attenuation",
        "main",
        [
            "manuscript/figure_data/study_dependence/structure_reuse_block_reweighted.tsv",
            "manuscript/figure_data/structural_grammar_evidence.tsv",
        ],
    ),
    row(
        "Figure 3",
        "A",
        "Composition-pruned tree with state tracks",
        "main",
        [
            "manuscript/figure_data/figure3_workflow_tree_nodes.tsv",
            "manuscript/figure_data/figure3_workflow_tree_segments.tsv",
            "manuscript/figure_data/fig02_prn_mechanism_calls.tsv",
        ],
    ),
    row(
        "Figure 3",
        "B",
        "Registered ASR frames",
        "main",
        ["manuscript/figure_data/asr_scenario_registry.tsv"],
    ),
    row(
        "Figure 3",
        "C",
        "Stochastic mapping",
        "main",
        ["manuscript/figure_data/asr_stochastic_mapping_summary.tsv"],
    ),
    row(
        "Figure 3",
        "D",
        "Scenario-level stress map",
        "main",
        ["manuscript/figure_data/asr_scenario_registry.tsv"],
    ),
    row(
        "Figure 4",
        "A",
        "Aligned country archive-context timelines",
        "main",
        [
            "manuscript/figure_data/selected_country/country_program_history_manifest.tsv",
            "manuscript/figure_data/fig01_prn_country_year_summary.tsv",
            "manuscript/figure_data/selected_country/selected_country_relative_year_plot_data.tsv",
        ],
    ),
    row(
        "Figure 4",
        "B",
        "Focal epoch IPW fractions and missingness bounds",
        "main",
        [
            "manuscript/figure_data/selected_country/country_epoch_prn_prevalence.tsv",
            "manuscript/figure_data/selected_country/country_epoch_bounds.tsv",
        ],
    ),
    row(
        "Figure 4",
        "C",
        "Relative-year disrupted-fraction trajectories",
        "main",
        [
            "manuscript/figure_data/selected_country/selected_country_relative_year_plot_data.tsv",
            "manuscript/figure_data/figure4_event_centered_pooled.tsv",
        ],
    ),
    row(
        "Figure 4",
        "D",
        "Architecture turnover between archive epochs",
        "main",
        ["manuscript/figure_data/selected_country/country_epoch_architecture_turnover_summary.tsv"],
    ),
    row(
        "Figure 5",
        "A",
        "Nested evidence chain and denominator flow",
        "main",
        ["manuscript/figure_data/fig05_evidence_chain_summary.tsv"],
    ),
    row(
        "Figure 5",
        "B",
        "prn- and comparator-locus specificity",
        "main",
        ["manuscript/figure_data/selected_country/prn_specificity_negative_control.tsv"],
    ),
    row(
        "Figure 5",
        "C",
        "Published PRN phenotype bridge and route-family tiers",
        "main",
        [
            "manuscript/figure_data/biology_bridge_external_context.tsv",
            "manuscript/supplementary/Supplementary_Table_10_Event_Class_Phenotype_Evidence_Tiers.tsv",
        ],
    ),
    row(
        "Figure 5",
        "D",
        "Validation evidence stack",
        "main",
        ["manuscript/figure_data/caller_validation_sensitivity_summary.tsv"],
    ),
    row(
        "Supplementary Figure 1",
        "A-C",
        "Sample-frame reconciliation and country-readiness audit",
        "supplementary",
        [
            "manuscript/submission_data/audit_ledgers/supplementary_table_sources/Supplementary_Table_33_Sample_Frame_Reconciliation.tsv",
            "manuscript/figure_data/selected_country/country_selection_scorecard.tsv",
            "manuscript/figure_data/selected_country/country_epoch_eligibility.tsv",
        ],
    ),
    row(
        "Supplementary Figure 2",
        "A-H",
        "Selected-country sensitivity, block-dependence, missingness diagnostics and read-linked transportability",
        "supplementary",
        [
            "manuscript/figure_data/selected_country/selected_country_year_sensitivity_summary.tsv",
            "manuscript/figure_data/selected_country/prn_interpretability_model.tsv",
            "manuscript/figure_data/selected_country/selected_country_dr_missingness_summary.tsv",
            "manuscript/figure_data/selected_country/selected_country_missingness_tipping_summary.tsv",
            "manuscript/figure_data/study_dependence/selected_country_block_bootstrap.tsv",
            "manuscript/figure_data/selected_country/selected_country_read_linked_transportability_ledger.tsv",
        ],
    ),
    row(
        "Supplementary Figure 3",
        "A-C",
        "Architecture-to-origin anchors and package-level rerun audit",
        "supplementary",
        [
            "manuscript/figure_data/selected_country/selected_country_validation_matrix.tsv",
            "manuscript/figure_data/selected_country/selected_country_origin_package_summary.tsv",
            "manuscript/figure_data/local_rooted_package_tree_summary.tsv",
        ],
    ),
    row(
        "Supplementary Figure 4",
        "A-C",
        "Country-by-epoch architecture turnover",
        "supplementary",
        [
            "manuscript/figure_data/selected_country/country_epoch_architecture_turnover.tsv",
            "manuscript/figure_data/selected_country/country_epoch_architecture_turnover_summary.tsv",
        ],
    ),
    row(
        "Supplementary Figure 5",
        "A-C",
        "Representativeness of tree subsets used for ancestral-state inference",
        "supplementary",
        ["manuscript/figure_data/asr_representativeness_adjustment_summary.tsv"],
    ),
    row(
        "Supplementary Figure 6",
        "A-D",
        "ASR robustness across reconstruction rules, stochastic mapping and balanced resampling",
        "supplementary",
        [
            "manuscript/figure_data/asr_scenario_registry.tsv",
            "manuscript/figure_data/asr_one_global_clone_summary.tsv",
            "manuscript/figure_data/asr_stochastic_mapping_summary.tsv",
            "manuscript/figure_data/study_dependence/asr_study_block_resampling.tsv",
        ],
    ),
    row(
        "Supplementary Figure 7",
        "A-C",
        "Heterogeneity of amplification after local origin or detection",
        "supplementary",
        [
            "manuscript/figure_data/selected_country/origin_burden_prevalence_shift.tsv",
            "manuscript/figure_data/selected_country/selected_country_relative_year_plot_data.tsv",
        ],
    ),
    row(
        "Supplementary Figure 8",
        "A-E",
        "Product-aware descriptive programme-context sensitivity analyses",
        "supplementary",
        [
            "manuscript/figure_data/figure5_association_model_panels.tsv",
            "manuscript/figure_data/figure5_leave_one_country_out_summary.tsv",
            "manuscript/figure_data/figure5_formulation_coverage.tsv",
            "manuscript/figure_data/supplementary_programme_class_summary.tsv",
        ],
    ),
    row(
        "Supplementary Figure 9",
        "A-E",
        "Focal-country readiness audit and descriptive controls",
        "supplementary",
        [
            "manuscript/figure_data/dynamic_identifiability_report.tsv",
            "manuscript/figure_data/dynamic_transmission_advantage_summary.tsv",
            "manuscript/figure_data/dynamic_transmission_advantage_predictions.tsv",
            "manuscript/figure_data/dynamic_counterfactual_summary.tsv",
        ],
    ),
    row(
        "Supplementary Figure 10",
        "A-C",
        "prn-locus structural signal specificity audit",
        "supplementary",
        ["manuscript/figure_data/selected_country/prn_specificity_negative_control.tsv"],
    ),
    row(
        "Supplementary Figure 14",
        "A-D",
        "Structural recurrence evidence compendium",
        "supplementary",
        [
            "manuscript/figure_data/event_definition_hierarchy_sensitivity.tsv",
            "manuscript/figure_data/structural_grammar_evidence.tsv",
            "manuscript/figure_data/event_specific_acquisition_summary.tsv",
            "manuscript/figure_data/prn_junction_confidence_matrix.tsv",
        ],
    ),
    row(
        "Supplementary Figure 15",
        "A-D",
        "Validation and caller-sensitivity compendium",
        "supplementary",
        [
            "manuscript/figure_data/published_overlap_concordance.tsv",
            "manuscript/figure_data/caller_validation_sensitivity_summary.tsv",
            "manuscript/figure_data/prn_threshold_grid_full.tsv",
            "manuscript/figure_data/biology_bridge_external_context.tsv",
        ],
    ),
]


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["figure_label", "panel_label", "panel_title", "final_placement", "frozen_input_files"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(ROWS)


if __name__ == "__main__":
    main()
