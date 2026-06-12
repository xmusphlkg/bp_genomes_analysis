#!/usr/bin/env python3
"""Build a manifest-scoped Snippy contig-mode plan for an expanded ML phylogeny."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.lib.project_paths import project_module_data_root, project_workflow_root


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SELECTION_MANIFEST = project_module_data_root("step5_phylogeny_asr") / "outputs" / "bp_phylogeny_manifest_balanced.tsv"
DEFAULT_BASE_PLAN = project_workflow_root() / "snippy_ctg" / "snippy_ctg_plan.tsv"
DEFAULT_COMPLETED_ROOT = project_workflow_root() / "snippy_ctg"
DEFAULT_WORKDIR = project_module_data_root("step5_phylogeny_asr") / "work" / "balanced_ml_phylogeny"

OUTPUT_COLUMNS = [
    "sample_id_canonical",
    "assembly_accession",
    "current_accession",
    "selection_present",
    "phylogeny_manifest_type",
    "phylogeny_tree_role",
    "phylogeny_selection_rule_id",
    "phylogeny_selection_reason",
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


def select_accession(row: dict[str, str]) -> str:
    return (row.get("current_accession") or row.get("assembly_accession") or "").strip()


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a custom Snippy contig-mode plan scoped to a Step5 phylogeny manifest."
    )
    parser.add_argument(
        "--selection-manifest",
        type=Path,
        default=DEFAULT_SELECTION_MANIFEST,
        help="Step5 phylogeny manifest defining the cohort to include.",
    )
    parser.add_argument(
        "--base-plan",
        type=Path,
        default=DEFAULT_BASE_PLAN,
        help="Canonical Snippy contig-mode plan built by the M3 bootstrap layer.",
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
        help="Output plan TSV.",
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
    selection_rows = load_tsv_rows(args.selection_manifest)
    base_rows = load_tsv_rows(args.base_plan)

    base_by_accession = {row.get("assembly_accession", "").strip(): row for row in base_rows}
    output_rows: list[dict[str, str]] = []
    exclusion_counts = Counter()
    qc_counts = Counter()
    summary = Counter()

    for selection_row in selection_rows:
        accession = select_accession(selection_row)
        base_row = base_by_accession.get(accession, {})
        include = base_row.get("include_in_snippy_ctg", "False") if base_row else "False"
        exclusion_reason = ""
        if base_row:
            exclusion_reason = base_row.get("exclusion_reason", "")
        else:
            exclusion_reason = "missing_from_base_plan"
        completed = args.completed_root / accession / "snps.aligned.fa"
        output_rows.append(
            {
                "sample_id_canonical": selection_row.get("sample_id_canonical", ""),
                "assembly_accession": selection_row.get("assembly_accession", ""),
                "current_accession": selection_row.get("current_accession", ""),
                "selection_present": "True",
                "phylogeny_manifest_type": selection_row.get("phylogeny_manifest_type", ""),
                "phylogeny_tree_role": selection_row.get("phylogeny_tree_role", ""),
                "phylogeny_selection_rule_id": selection_row.get("phylogeny_selection_rule_id", ""),
                "phylogeny_selection_reason": selection_row.get("phylogeny_selection_reason", ""),
                "fasta_path": base_row.get("fasta_path", ""),
                "assembly_exists": base_row.get("assembly_exists", "False") if base_row else "False",
                "qc_status": base_row.get("qc_status", "") if base_row else "",
                "qc_reasons": base_row.get("qc_reasons", "") if base_row else "",
                "has_reads": base_row.get("has_reads", "False") if base_row else "False",
                "prn_interpretable": base_row.get("prn_interpretable", "False") if base_row else "False",
                "prn_call_confidence": base_row.get("prn_call_confidence", "") if base_row else "",
                "evidence_tier": base_row.get("evidence_tier", "") if base_row else "",
                "preferred_snippy_mode": base_row.get("preferred_snippy_mode", "") if base_row else "",
                "planned_snippy_mode": base_row.get("planned_snippy_mode", "contigs") if base_row else "contigs",
                "include_in_snippy_ctg": include,
                "exclusion_reason": exclusion_reason,
                "snippy_ctg_completed": "True" if completed.exists() else "False",
            }
        )

        summary["selected_rows"] += 1
        if include == "True":
            summary["eligible_rows"] += 1
            if completed.exists():
                summary["completed_rows"] += 1
            else:
                summary["pending_rows"] += 1
            if truthy(base_row.get("has_reads", "False")):
                summary["eligible_rows_with_reads"] += 1
            else:
                summary["eligible_rows_contigs_only"] += 1
        else:
            summary["excluded_rows"] += 1
            exclusion_counts[exclusion_reason or "unknown"] += 1
        qc_counts[(base_row.get("qc_status", "") if base_row else "") or "missing"] += 1

    args.out_plan.parent.mkdir(parents=True, exist_ok=True)
    with args.out_plan.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    payload = {
        "selection_manifest": str(args.selection_manifest),
        "base_plan": str(args.base_plan),
        "completed_root": str(args.completed_root),
        "selected_rows": summary["selected_rows"],
        "eligible_rows": summary["eligible_rows"],
        "completed_rows": summary["completed_rows"],
        "pending_rows": summary["pending_rows"],
        "excluded_rows": summary["excluded_rows"],
        "eligible_rows_with_reads": summary["eligible_rows_with_reads"],
        "eligible_rows_contigs_only": summary["eligible_rows_contigs_only"],
        "qc_status_counts": dict(sorted(qc_counts.items())),
        "exclusion_reason_counts": dict(sorted(exclusion_counts.items())),
    }
    args.out_summary.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
