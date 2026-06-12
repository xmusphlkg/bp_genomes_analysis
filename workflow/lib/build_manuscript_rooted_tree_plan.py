#!/usr/bin/env python3
"""Build the manuscript-contract Snippy plan for the rooted 191-study-genome tree."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DECISION_LOG = ROOT / "manuscript" / "submission_data" / "cohort" / "master_cohort_decision_log.tsv"
DEFAULT_BASE_PLAN = ROOT / "outputs" / "workflow" / "snippy_ctg" / "snippy_ctg_plan.tsv"
DEFAULT_COMPLETED_ROOT = ROOT / "outputs" / "workflow" / "snippy_ctg"
DEFAULT_WORKDIR = ROOT / "outputs" / "workflow" / "manuscript_rooted_tree"

KEEP_STATUSES = {
    "primary_asr_tree",
    "core_alignment_only",
    "excluded_pre_gubbins_missingness",
}

OUTPUT_COLUMNS = [
    "sample_id_canonical",
    "assembly_accession",
    "current_accession",
    "selection_present",
    "phylogeny_manifest_type",
    "phylogeny_tree_role",
    "phylogeny_selection_rule_id",
    "phylogeny_selection_reason",
    "final_contract_status",
    "final_contract_reason",
    "fasta_path",
    "assembly_exists",
    "qc_status",
    "qc_reasons",
    "has_reads",
    "prn_interpretable",
    "prn_call_confidence",
    "evidence_tier",
    "preferred_snippy_mode",
    "planned_snippy_mode",
    "include_in_snippy_ctg",
    "exclusion_reason",
    "snippy_ctg_completed",
]


def load_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def truthy(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "t"}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--decision-log",
        type=Path,
        default=DEFAULT_DECISION_LOG,
        help="Submission decision log defining the manuscript rooted-tree contract.",
    )
    parser.add_argument(
        "--base-plan",
        type=Path,
        default=DEFAULT_BASE_PLAN,
        help="Canonical Snippy contig-mode plan built by M3.",
    )
    parser.add_argument(
        "--completed-root",
        type=Path,
        default=DEFAULT_COMPLETED_ROOT,
        help="Root containing per-accession Snippy contig-mode output directories.",
    )
    parser.add_argument(
        "--out-plan",
        type=Path,
        default=DEFAULT_WORKDIR / "snippy_ctg_plan.tsv",
        help="Output manuscript-contract Snippy plan TSV.",
    )
    parser.add_argument(
        "--out-summary",
        type=Path,
        default=DEFAULT_WORKDIR / "snippy_ctg_plan_summary.json",
        help="Output summary JSON.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    decision_rows = load_tsv_rows(args.decision_log)
    base_rows = load_tsv_rows(args.base_plan)

    base_by_accession = {row.get("assembly_accession", "").strip(): row for row in base_rows}
    output_rows: list[dict[str, str]] = []
    status_counts = Counter()
    exclusion_counts = Counter()
    qc_counts = Counter()
    summary = Counter()

    selected_rows = [
        row for row in decision_rows if row.get("final_contract_status", "").strip() in KEEP_STATUSES
    ]

    seen_accessions: set[str] = set()
    for selection_row in selected_rows:
        accession = selection_row.get("assembly_accession", "").strip()
        if not accession:
            raise ValueError(
                f"Encountered empty assembly_accession in decision log for sample "
                f"{selection_row.get('sample_id_canonical', '')}"
            )
        if accession in seen_accessions:
            raise ValueError(f"Duplicate contract accession in decision log: {accession}")
        seen_accessions.add(accession)

        base_row = base_by_accession.get(accession)
        if base_row is None:
            raise ValueError(f"Contract accession missing from base Snippy plan: {accession}")

        include = base_row.get("include_in_snippy_ctg", "False")
        if include != "True":
            raise ValueError(
                f"Contract accession is not eligible in the base Snippy plan: {accession} "
                f"(include_in_snippy_ctg={include!r})"
            )

        completed = args.completed_root / accession / "snps.aligned.fa"
        final_status = selection_row.get("final_contract_status", "").strip()
        final_reason = selection_row.get("final_contract_reason", "").strip()
        selection_reason = final_reason or final_status or "manuscript_rooted_tree_contract"

        output_rows.append(
            {
                "sample_id_canonical": selection_row.get("sample_id_canonical", ""),
                "assembly_accession": accession,
                "current_accession": accession,
                "selection_present": "True",
                "phylogeny_manifest_type": "submission_contract",
                "phylogeny_tree_role": "manuscript_rooted_ml_tree",
                "phylogeny_selection_rule_id": "submission_contract_rooted_tree_v1",
                "phylogeny_selection_reason": selection_reason,
                "final_contract_status": final_status,
                "final_contract_reason": final_reason,
                "fasta_path": base_row.get("fasta_path", ""),
                "assembly_exists": base_row.get("assembly_exists", "False"),
                "qc_status": base_row.get("qc_status", ""),
                "qc_reasons": base_row.get("qc_reasons", ""),
                "has_reads": base_row.get("has_reads", "False"),
                "prn_interpretable": base_row.get("prn_interpretable", "False"),
                "prn_call_confidence": base_row.get("prn_call_confidence", ""),
                "evidence_tier": base_row.get("evidence_tier", ""),
                "preferred_snippy_mode": base_row.get("preferred_snippy_mode", ""),
                "planned_snippy_mode": base_row.get("planned_snippy_mode", "contigs"),
                "include_in_snippy_ctg": include,
                "exclusion_reason": "",
                "snippy_ctg_completed": "True" if completed.exists() else "False",
            }
        )

        summary["selected_rows"] += 1
        status_counts[final_status] += 1
        qc_counts[base_row.get("qc_status", "") or "missing"] += 1
        if completed.exists():
            summary["completed_rows"] += 1
        else:
            summary["pending_rows"] += 1
            exclusion_counts["missing_completed_snippy_output"] += 1
        if truthy(base_row.get("has_reads", "False")):
            summary["selected_rows_with_reads"] += 1
        else:
            summary["selected_rows_contigs_only"] += 1

    args.out_plan.parent.mkdir(parents=True, exist_ok=True)
    with args.out_plan.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    payload = {
        "decision_log": str(args.decision_log),
        "base_plan": str(args.base_plan),
        "completed_root": str(args.completed_root),
        "selection_rule_id": "submission_contract_rooted_tree_v1",
        "selected_rows": summary["selected_rows"],
        "completed_rows": summary["completed_rows"],
        "pending_rows": summary["pending_rows"],
        "selected_rows_with_reads": summary["selected_rows_with_reads"],
        "selected_rows_contigs_only": summary["selected_rows_contigs_only"],
        "final_contract_status_counts": dict(sorted(status_counts.items())),
        "qc_status_counts": dict(sorted(qc_counts.items())),
        "exclusion_reason_counts": dict(sorted(exclusion_counts.items())),
    }
    args.out_summary.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
