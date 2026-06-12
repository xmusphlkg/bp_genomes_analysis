#!/usr/bin/env python3
"""Build the canonical Snippy contig-mode execution plan for M3.

This plan is derived from the repository single source of truth
(`state/manifest/manifest.tsv`) plus assembly QC results
(`workflow/assembly_qc/assembly_qc_stats.tsv`).

The current M3 bootstrap strategy is conservative:
  - require a present local FASTA
  - require assembly QC status PASS
  - record whether reads exist, but still use contig mode for the complete-cohort seed tree

Outputs:
  - workflow/snippy_ctg/snippy_ctg_plan.tsv
  - workflow/snippy_ctg/snippy_ctg_plan_report.json
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from workflow.lib.project_paths import project_workflow_root

DEFAULT_MANIFEST = ROOT / "state" / "manifest" / "manifest.tsv"
DEFAULT_QC = project_workflow_root() / "assembly_qc" / "assembly_qc_stats.tsv"
DEFAULT_ASSEMBLY_DIR = ROOT / "pertussis_data" / "bp_genomes_qc" / "assemblies"
DEFAULT_OUTDIR = project_workflow_root() / "snippy_ctg"


TRUTHY = {"1", "true", "yes", "y", "t"}


def as_bool(value: str) -> bool:
    return str(value).strip().lower() in TRUTHY


def load_qc_map(qc_path: Path) -> dict[str, dict[str, str]]:
    qc_map: dict[str, dict[str, str]] = {}
    if not qc_path.exists():
        return qc_map

    with qc_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            qc_map[row["sample_id_canonical"]] = row
    return qc_map


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Snippy contig-mode sample plan")
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--qc", default=str(DEFAULT_QC))
    parser.add_argument("--assembly-dir", default=str(DEFAULT_ASSEMBLY_DIR))
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    qc_path = Path(args.qc)
    assembly_dir = Path(args.assembly_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not manifest_path.exists():
        raise SystemExit(f"Manifest not found: {manifest_path}")

    qc_map = load_qc_map(qc_path)

    plan_rows: list[dict[str, str]] = []
    counts = Counter()
    qc_reason_counts = Counter()
    evidence_tier_counts = Counter()

    with manifest_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample_id = row["sample_id_canonical"]
            assembly_accession = row["assembly_accession"]
            fasta_path = assembly_dir / f"{assembly_accession}.fasta"
            qc_row = qc_map.get(sample_id, {})
            qc_status = qc_row.get("qc_status", "NO_QC")
            qc_reasons = qc_row.get("qc_reasons", "")
            assembly_exists = fasta_path.is_file() and fasta_path.stat().st_size > 0
            has_reads = as_bool(row.get("has_reads", ""))
            prn_interpretable = as_bool(row.get("prn_interpretable", ""))
            evidence_tier = row.get("evidence_tier", "") or "unknown"

            include_in_snippy_ctg = assembly_exists and qc_status == "PASS"
            if include_in_snippy_ctg:
                exclusion_reason = ""
            elif not assembly_exists:
                exclusion_reason = "missing_fasta"
            else:
                exclusion_reason = f"qc_{qc_status.lower()}"
                if qc_reasons:
                    exclusion_reason = f"{exclusion_reason}:{qc_reasons}"

            preferred_snippy_mode = "reads" if has_reads else "contigs"

            plan_rows.append(
                {
                    "sample_id_canonical": sample_id,
                    "assembly_accession": assembly_accession,
                    "fasta_path": str(fasta_path),
                    "assembly_exists": "True" if assembly_exists else "False",
                    "qc_status": qc_status,
                    "qc_reasons": qc_reasons,
                    "has_reads": "True" if has_reads else "False",
                    "prn_interpretable": "True" if prn_interpretable else "False",
                    "prn_call_confidence": row.get("prn_call_confidence", ""),
                    "evidence_tier": evidence_tier,
                    "preferred_snippy_mode": preferred_snippy_mode,
                    "planned_snippy_mode": "contigs",
                    "include_in_snippy_ctg": "True" if include_in_snippy_ctg else "False",
                    "exclusion_reason": exclusion_reason,
                }
            )

            counts["manifest_total"] += 1
            counts["assembly_present" if assembly_exists else "assembly_missing"] += 1
            counts[f"qc_{qc_status.lower()}"] += 1
            evidence_tier_counts[evidence_tier] += 1

            if include_in_snippy_ctg:
                counts["include_in_snippy_ctg"] += 1
                if has_reads:
                    counts["include_with_reads"] += 1
                else:
                    counts["include_contigs_only"] += 1
            else:
                counts[f"exclude_{exclusion_reason.split(':', 1)[0]}"] += 1

            if qc_reasons:
                for reason in qc_reasons.split(";"):
                    if reason:
                        qc_reason_counts[reason] += 1

    fieldnames = [
        "sample_id_canonical",
        "assembly_accession",
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
    ]

    plan_path = outdir / "snippy_ctg_plan.tsv"
    with plan_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(plan_rows)

    report = {
        "build_date": __import__("datetime").datetime.now().isoformat(),
        "manifest_total": counts["manifest_total"],
        "assembly_present": counts["assembly_present"],
        "assembly_missing": counts["assembly_missing"],
        "include_in_snippy_ctg": counts["include_in_snippy_ctg"],
        "include_with_reads": counts["include_with_reads"],
        "include_contigs_only": counts["include_contigs_only"],
        "qc_status_counts": {k[3:]: v for k, v in counts.items() if k.startswith("qc_")},
        "evidence_tier_counts": dict(evidence_tier_counts),
        "qc_reason_counts": dict(qc_reason_counts),
        "selection_rule": "include only rows with local FASTA present and assembly QC status PASS",
        "notes": [
            "Current M3 bootstrap uses contig-mode Snippy to seed the full-cohort phylogeny inputs.",
            "Rows with reads are flagged for later read-backed refinement but remain in contig mode here.",
        ],
    }

    report_path = outdir / "snippy_ctg_plan_report.json"
    with report_path.open("w") as handle:
        json.dump(report, handle, indent=2, default=str)

    print("=== Snippy-ctg Plan Summary ===")
    print(f"Manifest rows:           {counts['manifest_total']}")
    print(f"Assemblies present:      {counts['assembly_present']}")
    print(f"Assemblies missing:      {counts['assembly_missing']}")
    print(f"Included for Snippy-ctg: {counts['include_in_snippy_ctg']}")
    print(f"Included with reads:     {counts['include_with_reads']}")
    print(f"Included contigs only:   {counts['include_contigs_only']}")
    print(f"Plan:                    {plan_path}")
    print(f"Report:                  {report_path}")


if __name__ == "__main__":
    main()
