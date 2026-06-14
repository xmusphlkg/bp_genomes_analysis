#!/usr/bin/env python3
"""Build manuscript-facing cohort decision logs for the submission scaffold."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


DATA_HOME = Path(
    os.environ.get(
        "PERTUSSIS_PROJECT_DATA_ROOT",
        str(repo_root() / "pertussis_data" / "pertussis_gene"),
    )
)
STEP1_OUTPUTS = DATA_HOME / "step1_ingest" / "outputs"
STEP4_OUTPUTS = DATA_HOME / "step4_prn_validation" / "outputs"
WORKFLOW_MANIFEST = repo_root() / "outputs" / "workflow" / "manifest" / "manifest.tsv"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def count_values(df: pd.DataFrame, column: str, value: str) -> int:
    if column not in df.columns:
        return 0
    return int(df[column].astype(str).str.lower().eq(value.lower()).sum())


def count_in(df: pd.DataFrame, column: str, values: set[str]) -> int:
    if column not in df.columns:
        return 0
    normalized = {value.lower() for value in values}
    return int(df[column].astype(str).str.lower().isin(normalized).sum())


def fasta_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def rel_source(path: Path) -> str:
    for base in (repo_root(), DATA_HOME):
        try:
            return str(path.relative_to(base))
        except ValueError:
            continue
    return str(path)


def prn_interpretability_group(mechanism: str) -> str:
    if mechanism in {
        "intact",
        "coding_disrupted_is481",
        "coding_disrupted_inversion_or_rearrangement",
        "coding_disrupted_other",
    }:
        return "interpretable"
    if mechanism == "uncertain_fragmented_assembly":
        return "uncertain_fragmented"
    return "insufficient"


def final_contract_status(row: pd.Series) -> str:
    if row["in_primary_asr_tree"] == "yes":
        return "primary_asr_tree"
    if row["excluded_pre_gubbins_missingness"] == "yes":
        return "excluded_pre_gubbins_missingness"
    if row["in_core_full_alignment"] == "yes":
        return "core_alignment_only"
    if row["include_in_snippy_ctg"] == "yes":
        return "snippy_eligible_but_not_in_current_alignment"
    return "excluded_before_snippy_alignment"


def final_contract_reason(row: pd.Series) -> str:
    if row["in_primary_asr_tree"] == "yes":
        return "retained_to_rooted_ml_tree_and_asr"
    if row["excluded_pre_gubbins_missingness"] == "yes":
        return f"missingness_fraction={row['missing_fraction_pre_gubbins']}"
    if row["in_core_full_alignment"] == "yes":
        return "present_in_core_alignment_but_absent_from_tip_state_export"
    if row["include_in_snippy_ctg"] == "yes":
        return (
            "operational_gap_between_snippy_eligible_set_and_current_core_alignment;"
            "selection_logic_not_yet_consolidated_as_manuscript_contract"
        )
    return row["snippy_exclusion_reason"] or "not_selected_before_snippy"


def build_decision_log(
    *,
    public_manifest: Path,
    qc_manifest: Path,
    mechanism_calls: Path,
    validation_subset: Path,
    validation_results: Path,
    snippy_plan: Path,
    core_full_alignment: Path,
    pre_gubbins_missingness: Path,
    unpruned_tip_states: Path,
    tip_states: Path,
    workflow_manifest: Path | None,
    out_log: Path,
    out_summary: Path,
    out_supp_flow: Path | None,
) -> None:
    public_df = read_tsv(public_manifest)
    qc_df = read_tsv(qc_manifest)
    mech_df = read_tsv(mechanism_calls)
    validation_df = read_tsv(validation_subset)
    validation_results_df = read_tsv(validation_results)
    snippy_df = read_tsv(snippy_plan)
    missingness_df = read_tsv(pre_gubbins_missingness)
    unpruned_tip_df = read_tsv(unpruned_tip_states)
    tip_df = read_tsv(tip_states)
    workflow_df = read_tsv(workflow_manifest) if workflow_manifest and workflow_manifest.exists() else pd.DataFrame()
    core_ids = set(fasta_ids(core_full_alignment))

    base_cols = [
        "sample_id_canonical",
        "assembly_accession",
        "country",
        "year",
        "source_database",
        "assembly_level",
        "n_contigs",
        "contig_n50",
        "total_sequence_length",
        "raw_reads_available",
        "raw_read_run_count",
        "raw_read_link_status",
        "record_decision",
    ]
    base = qc_df[base_cols].copy()
    base["retained_in_manuscript_qc_cohort"] = "yes"

    mech_cols = [
        "sample_id_canonical",
        "prn_mechanism_call",
        "prn_call_confidence",
        "prn_event_id",
    ]
    base = base.merge(mech_df[mech_cols], on="sample_id_canonical", how="left")
    base["prn_interpretability_group"] = base["prn_mechanism_call"].map(prn_interpretability_group)
    base["prn_interpretable"] = base["prn_interpretability_group"].eq("interpretable").map(yes_no)
    base["prn_disrupted"] = base["prn_mechanism_call"].isin(
        [
            "coding_disrupted_is481",
            "coding_disrupted_inversion_or_rearrangement",
            "coding_disrupted_other",
        ]
    ).map(yes_no)

    validation_cols = [
        "sample_id_canonical",
        "selection_stratum",
        "selection_reason",
    ]
    validation_merge = validation_df[validation_cols].drop_duplicates("sample_id_canonical")
    base = base.merge(validation_merge, on="sample_id_canonical", how="left")
    base["in_validation_subset"] = base["selection_stratum"].notna().map(yes_no)

    snippy_cols = [
        "sample_id_canonical",
        "qc_status",
        "qc_reasons",
        "include_in_snippy_ctg",
        "exclusion_reason",
    ]
    snippy_merge = snippy_df[snippy_cols].drop_duplicates("sample_id_canonical")
    base = base.merge(snippy_merge, on="sample_id_canonical", how="left")
    base["snippy_qc_status"] = base["qc_status"].fillna("")
    base["snippy_qc_reasons"] = base["qc_reasons"].fillna("")
    base["include_in_snippy_ctg"] = (
        base["include_in_snippy_ctg"].fillna("False").str.lower().eq("true").map(yes_no)
    )
    base["snippy_exclusion_reason"] = base["exclusion_reason"].fillna("")
    base = base.drop(columns=["qc_status", "qc_reasons", "exclusion_reason"])

    base["in_core_full_alignment"] = base["assembly_accession"].isin(core_ids).map(yes_no)

    missingness_keep = {
        row["sequence_id"]: row["keep"].strip().lower() == "true"
        for _, row in missingness_df.iterrows()
    }
    missingness_frac = {
        row["sequence_id"]: row["missing_fraction"]
        for _, row in missingness_df.iterrows()
    }
    base["present_in_pre_gubbins_missingness_table"] = base["assembly_accession"].isin(missingness_keep).map(yes_no)
    base["excluded_pre_gubbins_missingness"] = base["assembly_accession"].map(
        lambda acc: yes_no(acc in missingness_keep and not missingness_keep[acc])
    )
    base["missing_fraction_pre_gubbins"] = base["assembly_accession"].map(missingness_frac).fillna("")

    unpruned_nonref_tip_df = unpruned_tip_df[
        unpruned_tip_df["is_reference"].fillna("False").str.lower().ne("true")
    ].copy()
    nonref_tip_df = tip_df[tip_df["is_reference"].fillna("False").str.lower().ne("true")].copy()
    tip_cols = [
        "assembly_accession",
        "prn_state",
        "country_iso3",
        "year",
    ]
    base = base.merge(
        nonref_tip_df[tip_cols].rename(
            columns={
                "prn_state": "asr_tip_prn_state",
                "country_iso3": "asr_country_iso3",
                "year": "asr_year",
            }
        ),
        on="assembly_accession",
        how="left",
    )
    base["in_primary_asr_tree"] = base["asr_tip_prn_state"].notna().map(yes_no)

    base["final_contract_status"] = base.apply(final_contract_status, axis=1)
    base["final_contract_reason"] = base.apply(final_contract_reason, axis=1)

    ordered_cols = [
        "sample_id_canonical",
        "assembly_accession",
        "country",
        "year",
        "source_database",
        "assembly_level",
        "n_contigs",
        "contig_n50",
        "total_sequence_length",
        "record_decision",
        "retained_in_manuscript_qc_cohort",
        "prn_mechanism_call",
        "prn_call_confidence",
        "prn_event_id",
        "prn_interpretability_group",
        "prn_interpretable",
        "prn_disrupted",
        "raw_reads_available",
        "raw_read_run_count",
        "raw_read_link_status",
        "in_validation_subset",
        "selection_stratum",
        "selection_reason",
        "snippy_qc_status",
        "snippy_qc_reasons",
        "include_in_snippy_ctg",
        "snippy_exclusion_reason",
        "in_core_full_alignment",
        "present_in_pre_gubbins_missingness_table",
        "excluded_pre_gubbins_missingness",
        "missing_fraction_pre_gubbins",
        "in_primary_asr_tree",
        "asr_tip_prn_state",
        "asr_country_iso3",
        "asr_year",
        "final_contract_status",
        "final_contract_reason",
    ]
    decision_log = base[ordered_cols].sort_values(
        by=["country", "year", "sample_id_canonical", "assembly_accession"],
        na_position="last",
    )

    public_count = len(public_df)
    qc_anchor_count = len(qc_df)
    combined_manifest_count = len(workflow_df) if not workflow_df.empty else len(mech_df)
    public_assembly_count = count_values(workflow_df, "data_origin", "public_genome_assembly")
    public_read_rescue_count = count_values(workflow_df, "data_origin", "public_read_rescue")
    raw_read_assembly_count = count_values(workflow_df, "data_origin", "public_raw_read_assembly")
    frozen_evidence_chain = repo_root() / "manuscript" / "figure_data" / "fig05_evidence_chain_summary.tsv"
    interpretable_count = count_in(
        mech_df,
        "prn_mechanism_call",
        {
            "intact",
            "coding_disrupted_is481",
            "coding_disrupted_inversion_or_rearrangement",
            "coding_disrupted_other",
        },
    )
    insufficient_count = count_values(mech_df, "prn_mechanism_call", "insufficient_data")
    uncertain_count = count_values(mech_df, "prn_mechanism_call", "uncertain_fragmented_assembly")
    insufficient_or_uncertain_count = insufficient_count + uncertain_count
    mechanism_summary_source = rel_source(mechanism_calls)
    interpretable_definition = "Mechanism call in {intact, coding_disrupted_is481, coding_disrupted_inversion_or_rearrangement, coding_disrupted_other}."
    interpretable_notes = "Current mechanism-call table interpretable subset used for broad mechanism summaries."
    noninterpretable_definition = "Mechanism call in {insufficient_data, uncertain_fragmented_assembly}."
    noninterpretable_notes = f"{insufficient_count} insufficient plus {uncertain_count} uncertain fragmented assemblies."
    if frozen_evidence_chain.exists():
        chain = read_tsv(frozen_evidence_chain)
        chain["n_numeric"] = pd.to_numeric(chain["n"], errors="coerce").fillna(0).astype(int)

        def chain_count(stage_id: str) -> int:
            return int(chain.loc[chain["stage_id"].eq(stage_id), "n_numeric"].sum())

        intact_count = chain_count("intact_boundary")
        disrupted_count = chain_count("structurally_resolved_disrupted")
        interpretable_count = chain_count("interpretable_event_phenotype")
        insufficient_or_uncertain_count = chain_count("noninterpretable_uncertain")
        mechanism_summary_source = rel_source(frozen_evidence_chain)
        interpretable_definition = (
            "Frozen manuscript-facing event/phenotype frame: intact-locus boundary calls plus structurally "
            "resolved disrupted event-class calls."
        )
        interpretable_notes = (
            f"Frozen frame comprises {intact_count} intact-locus boundary calls and "
            f"{disrupted_count} structurally resolved disrupted event-class calls."
        )
        noninterpretable_definition = "Retained records outside the frozen manuscript-facing event/phenotype frame."
        noninterpretable_notes = (
            "Frozen non-interpretable or uncertain records enter missingness bounds rather than event-concentration denominators."
        )
    raw_read_linked_count = int(qc_df["raw_reads_available"].fillna("False").str.lower().eq("true").sum())
    raw_read_linked_count += raw_read_assembly_count
    read_validation_table_count = len(validation_results_df)
    snippy_eligible_count = int(decision_log["include_in_snippy_ctg"].eq("yes").sum())
    core_alignment_total_count = len(core_ids)
    core_alignment_nonref_count = len([seq_id for seq_id in core_ids if seq_id != "Reference"])
    missingness_excluded_count = int((missingness_df["keep"].str.lower() == "false").sum())
    unpruned_asr_tip_total_count = len(unpruned_tip_df)
    unpruned_asr_tip_nonref_count = int(unpruned_nonref_tip_df.shape[0])
    unpruned_asr_intact_nonref_count = int(unpruned_nonref_tip_df["prn_state"].eq("intact").sum())
    unpruned_asr_disrupted_nonref_count = int(unpruned_nonref_tip_df["prn_state"].eq("disrupted").sum())
    unpruned_asr_insufficient_nonref_count = int(unpruned_nonref_tip_df["prn_state"].eq("insufficient_data").sum())
    primary_asr_tip_total_count = len(tip_df)
    primary_asr_tip_nonref_count = int(nonref_tip_df.shape[0])
    primary_asr_intact_nonref_count = int(nonref_tip_df["prn_state"].eq("intact").sum())
    primary_asr_disrupted_nonref_count = int(nonref_tip_df["prn_state"].eq("disrupted").sum())
    primary_asr_insufficient_nonref_count = int(nonref_tip_df["prn_state"].eq("insufficient_data").sum())

    summary_rows = [
        {
            "stage_id": "S00",
            "stage_name": "public_assembly_universe",
            "n_rows": public_count,
            "delta_from_previous": "",
            "source_file": rel_source(public_manifest),
            "rule_or_definition": "All downloaded public assemblies before manuscript-stage retention.",
            "notes": "Universe count used for the top of the funnel.",
        },
        {
            "stage_id": "S01",
            "stage_name": "combined_manuscript_manifest",
            "n_rows": combined_manifest_count,
            "delta_from_previous": combined_manifest_count - public_count,
            "source_file": rel_source(workflow_manifest) if workflow_manifest and workflow_manifest.exists() else rel_source(mechanism_calls),
            "rule_or_definition": "Current combined manifest after public-assembly QC, read-rescue provenance assignment, and de novo public-read augmentation.",
            "notes": (
                f"Comprises {public_assembly_count} public genome assemblies, "
                f"{public_read_rescue_count} public-read rescue records, and "
                f"{raw_read_assembly_count} de novo public_raw_read_assembly records; "
                f"the retained public-assembly QC anchor remains {qc_anchor_count} records."
            ),
        },
        {
            "stage_id": "S02",
            "stage_name": "prn_interpretable_subset",
            "n_rows": interpretable_count,
            "delta_from_previous": interpretable_count - combined_manifest_count,
            "source_file": mechanism_summary_source,
            "rule_or_definition": interpretable_definition,
            "notes": interpretable_notes,
        },
        {
            "stage_id": "S03",
            "stage_name": "prn_insufficient_or_uncertain_subset",
            "n_rows": insufficient_or_uncertain_count,
            "delta_from_previous": "",
            "source_file": mechanism_summary_source,
            "rule_or_definition": noninterpretable_definition,
            "notes": noninterpretable_notes,
        },
        {
            "stage_id": "S04",
            "stage_name": "raw_read_linked_subset",
            "n_rows": raw_read_linked_count,
            "delta_from_previous": "",
            "source_file": rel_source(qc_manifest),
            "rule_or_definition": "QC-retained public assemblies with raw_reads_available == True plus de novo public_raw_read_assembly records.",
            "notes": f"Includes {raw_read_assembly_count} genomes assembled directly from public reads.",
        },
        {
            "stage_id": "S05",
            "stage_name": "read_validation_table",
            "n_rows": read_validation_table_count,
            "delta_from_previous": "",
            "source_file": rel_source(validation_results),
            "rule_or_definition": "Targeted 60-row read-validation table after the 2026-04-12 recovery merge, spanning disrupted, intact, and insufficient states.",
            "notes": "Validation is deliberately partial and functions as both corroboration and rescue/negative-control evidence.",
        },
        {
            "stage_id": "S06",
            "stage_name": "snippy_ctg_qc_pass",
            "n_rows": snippy_eligible_count,
            "delta_from_previous": "",
            "source_file": rel_source(snippy_plan),
            "rule_or_definition": "include_in_snippy_ctg == True in the current Snippy contig-mode plan.",
            "notes": "This is an operational ML-tree eligibility layer, not yet a manuscript-facing representativeness contract.",
        },
        {
            "stage_id": "S07",
            "stage_name": "core_full_alignment_total",
            "n_rows": core_alignment_total_count,
            "delta_from_previous": "",
            "source_file": rel_source(core_full_alignment),
            "rule_or_definition": "All FASTA entries in core.full.aln, including the Reference sequence.",
            "notes": f"Contains {core_alignment_nonref_count} study genomes plus Reference.",
        },
        {
            "stage_id": "S08",
            "stage_name": "excluded_pre_gubbins_missingness",
            "n_rows": missingness_excluded_count,
            "delta_from_previous": "",
            "source_file": rel_source(pre_gubbins_missingness),
            "rule_or_definition": "keep == False in the pre-Gubbins missingness filter.",
            "notes": "Current rerun excludes only GCA_000212975.1 at missing_fraction=0.889548.",
        },
        {
            "stage_id": "S09",
            "stage_name": "unpruned_asr_tree_total_tips",
            "n_rows": unpruned_asr_tip_total_count,
            "delta_from_previous": "",
            "source_file": rel_source(unpruned_tip_states),
            "rule_or_definition": "All tip rows exported into the rooted ML-tree ASR package before composition pruning, including Reference.",
            "notes": f"Contains {unpruned_asr_tip_nonref_count} study genomes plus Reference; retained as the unpruned comparability frame.",
        },
        {
            "stage_id": "S10",
            "stage_name": "unpruned_asr_tree_nonreference_state_breakdown",
            "n_rows": unpruned_asr_tip_nonref_count,
            "delta_from_previous": "",
            "source_file": rel_source(unpruned_tip_states),
            "rule_or_definition": "Non-reference study genomes in the unpruned rooted ML-tree ASR package.",
            "notes": (
                f"Non-reference state counts: intact={unpruned_asr_intact_nonref_count}, "
                f"disrupted={unpruned_asr_disrupted_nonref_count}, "
                f"insufficient_data={unpruned_asr_insufficient_nonref_count}."
            ),
        },
        {
            "stage_id": "S11",
            "stage_name": "composition_pruned_primary_asr_tree_total_tips",
            "n_rows": primary_asr_tip_total_count,
            "delta_from_previous": "",
            "source_file": rel_source(tip_states),
            "rule_or_definition": "Primary ASR quality frame after pruning the 33 nonreference IQ-TREE composition-failed tips, retaining Reference.",
            "notes": f"Contains {primary_asr_tip_nonref_count} study genomes plus Reference.",
        },
        {
            "stage_id": "S12",
            "stage_name": "composition_pruned_primary_asr_tree_nonreference_state_breakdown",
            "n_rows": primary_asr_tip_nonref_count,
            "delta_from_previous": "",
            "source_file": rel_source(tip_states),
            "rule_or_definition": "Non-reference study genomes in the composition-pruned primary ASR quality frame.",
            "notes": (
                f"Non-reference state counts: intact={primary_asr_intact_nonref_count}, "
                f"disrupted={primary_asr_disrupted_nonref_count}, "
                f"insufficient_data={primary_asr_insufficient_nonref_count}."
            ),
        },
    ]
    summary_df = pd.DataFrame(summary_rows)

    out_log.parent.mkdir(parents=True, exist_ok=True)
    decision_log.to_csv(out_log, sep="\t", index=False)
    summary_df.to_csv(out_summary, sep="\t", index=False)
    if out_supp_flow is not None:
        out_supp_flow.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(out_supp_flow, sep="\t", index=False)


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--public-manifest",
        type=Path,
        default=STEP1_OUTPUTS / "bp_public_genome_manifest.tsv",
    )
    parser.add_argument(
        "--qc-manifest",
        type=Path,
        default=STEP1_OUTPUTS / "bp_public_genome_qc_manifest.tsv",
    )
    parser.add_argument(
        "--mechanism-calls",
        type=Path,
        default=STEP4_OUTPUTS / "bp_prn_mechanism_calls.tsv",
    )
    parser.add_argument(
        "--validation-subset",
        type=Path,
        default=STEP4_OUTPUTS / "bp_prn_validation_subset.tsv",
    )
    parser.add_argument(
        "--validation-results",
        type=Path,
        default=STEP4_OUTPUTS / "bp_prn_read_validation.tsv",
    )
    parser.add_argument(
        "--snippy-plan",
        type=Path,
        default=root / "outputs" / "workflow" / "snippy_ctg" / "snippy_ctg_plan.tsv",
    )
    parser.add_argument(
        "--core-full-alignment",
        type=Path,
        default=root / "outputs" / "workflow" / "phylo" / "core.full.aln",
    )
    parser.add_argument(
        "--pre-gubbins-missingness",
        type=Path,
        default=root / "outputs" / "workflow" / "phylo" / "pre_gubbins_missingness.tsv",
    )
    parser.add_argument(
        "--unpruned-tip-states",
        type=Path,
        default=root / "outputs" / "workflow" / "asr" / "tip_states.tsv",
    )
    parser.add_argument(
        "--tip-states",
        type=Path,
        default=root / "outputs" / "workflow" / "asr_sensitivity" / "composition_filtered" / "tip_states.tsv",
    )
    parser.add_argument(
        "--workflow-manifest",
        type=Path,
        default=WORKFLOW_MANIFEST,
    )
    parser.add_argument(
        "--out-log",
        type=Path,
        default=root / "manuscript" / "submission_data" / "cohort" / "master_cohort_decision_log.tsv",
    )
    parser.add_argument(
        "--out-summary",
        type=Path,
        default=root / "manuscript" / "submission_data" / "cohort" / "master_cohort_flow_summary.tsv",
    )
    parser.add_argument(
        "--out-supp-flow",
        type=Path,
        default=root / "manuscript" / "supplementary" / "Supplementary_Table_4_Cohort_Flow_and_Tree_Selection.tsv",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    build_decision_log(
        public_manifest=args.public_manifest,
        qc_manifest=args.qc_manifest,
        mechanism_calls=args.mechanism_calls,
        validation_subset=args.validation_subset,
        validation_results=args.validation_results,
        snippy_plan=args.snippy_plan,
        core_full_alignment=args.core_full_alignment,
        pre_gubbins_missingness=args.pre_gubbins_missingness,
        unpruned_tip_states=args.unpruned_tip_states,
        tip_states=args.tip_states,
        workflow_manifest=args.workflow_manifest,
        out_log=args.out_log,
        out_summary=args.out_summary,
        out_supp_flow=args.out_supp_flow,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
