#!/usr/bin/env python3
"""Build batch and blocker manifests for Step4 raw-read IS validation."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from step4_02_scan_prn_mechanisms import load_tsv_rows, normalize_text, repo_root, write_tsv


PROJECT_ROOT = Path(
    os.environ.get(
        "PERTUSSIS_PROJECT_DATA_ROOT",
        str(repo_root() / "pertussis_data" / "pertussis_gene"),
    )
)


BATCH_COLUMNS = [
    "sample_id_canonical",
    "sra_run_accession",
    "ena_run_accession",
    "read_accession_primary",
    "prn_event_id",
    "prn_mechanism_call",
    "prn_call_confidence",
    "raw_read_link_status",
    "resolved_identifier",
    "reads_1_path",
    "reads_2_path",
    "snippy_bam_path",
    "missing_inputs",
    "batch_status",
    "batch_note",
]

MISSING_COLUMNS = [
    "sample_id_canonical",
    "sra_run_accession",
    "ena_run_accession",
    "read_accession_primary",
    "prn_event_id",
    "prn_mechanism_call",
    "prn_call_confidence",
    "raw_read_link_status",
    "resolved_identifier",
    "reads_1_path",
    "reads_2_path",
    "snippy_bam_path",
    "missing_inputs",
    "batch_status",
    "batch_note",
]


def split_accessions(value: str) -> list[str]:
    tokens: list[str] = []
    for token in normalize_text(value).split(";"):
        cleaned = normalize_text(token)
        if cleaned:
            tokens.append(cleaned)
    return tokens


def unique_preserve_order(values: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def candidate_identifiers(row: dict[str, str], sample_id: str) -> list[str]:
    candidates: list[str] = [sample_id]
    for field in (
        "read_accession_primary",
        "sra_sample_accession",
        "ena_sample_accession",
        "sra_run_accession",
        "ena_run_accession",
    ):
        candidates.extend(split_accessions(row.get(field, "")))
    return unique_preserve_order([item for item in candidates if item])


def resolve_paths(row: dict[str, str], sample_id: str, reads_root: Path, snippy_root: Path) -> dict[str, object]:
    best: dict[str, object] = {
        "resolved_identifier": sample_id,
        "reads_1_path": reads_root / f"{sample_id}_1.fastq.gz",
        "reads_2_path": reads_root / f"{sample_id}_2.fastq.gz",
        "snippy_bam_path": snippy_root / sample_id / "snps.bam",
        "has_reads": False,
        "has_bam": False,
        "score": -1,
    }

    for rank, identifier in enumerate(candidate_identifiers(row, sample_id)):
        reads_1 = reads_root / f"{identifier}_1.fastq.gz"
        reads_2 = reads_root / f"{identifier}_2.fastq.gz"
        snippy_bam = snippy_root / identifier / "snps.bam"
        has_reads = reads_1.exists() and reads_2.exists()
        has_bam = snippy_bam.exists()
        score = (2 if has_reads else 0) + (1 if has_bam else 0)
        if score > int(best["score"]):
            best = {
                "resolved_identifier": identifier,
                "reads_1_path": reads_1,
                "reads_2_path": reads_2,
                "snippy_bam_path": snippy_bam,
                "has_reads": has_reads,
                "has_bam": has_bam,
                "score": score,
                "rank": rank,
            }

    return best


def missing_inputs(has_reads: bool, has_bam: bool) -> str:
    missing: list[str] = []
    if not has_reads:
        missing.extend(["reads_1_fastq", "reads_2_fastq"])
    if not has_bam:
        missing.append("snippy_bam")
    return ",".join(missing)


def to_manifest_row(row: dict[str, str], resolved: dict[str, object], *, status: str, note: str) -> dict[str, str]:
    return {
        "sample_id_canonical": normalize_text(row.get("sample_id_canonical", "")),
        "sra_run_accession": normalize_text(row.get("sra_run_accession", "")),
        "ena_run_accession": normalize_text(row.get("ena_run_accession", "")),
        "read_accession_primary": normalize_text(row.get("read_accession_primary", "")),
        "prn_event_id": normalize_text(row.get("prn_event_id", "")),
        "prn_mechanism_call": normalize_text(row.get("prn_mechanism_call", "")),
        "prn_call_confidence": normalize_text(row.get("prn_call_confidence", "")),
        "raw_read_link_status": normalize_text(row.get("raw_read_link_status", "")),
        "resolved_identifier": normalize_text(str(resolved.get("resolved_identifier", ""))),
        "reads_1_path": str(resolved.get("reads_1_path", "")),
        "reads_2_path": str(resolved.get("reads_2_path", "")),
        "snippy_bam_path": str(resolved.get("snippy_bam_path", "")),
        "missing_inputs": missing_inputs(bool(resolved.get("has_reads", False)), bool(resolved.get("has_bam", False))),
        "batch_status": status,
        "batch_note": note,
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build selected-sample and blocker manifests for Step4 read validation."
    )
    parser.add_argument(
        "--subset",
        type=Path,
        default=PROJECT_ROOT / "step4_prn_validation" / "outputs" / "bp_prn_validation_subset.tsv",
        help="Validation subset manifest.",
    )
    parser.add_argument(
        "--reads-root",
        type=Path,
        default=PROJECT_ROOT / "workflow" / "reads_clean",
        help="Directory containing cleaned paired FASTQ files.",
    )
    parser.add_argument(
        "--snippy-root",
        type=Path,
        default=PROJECT_ROOT / "workflow" / "snippy",
        help="Directory containing read-mode Snippy outputs with snps.bam per sample.",
    )
    parser.add_argument(
        "--batch-label",
        default="current",
        help="Batch label used to determine default output locations.",
    )
    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Skip the first N eligible samples before selecting the batch.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Maximum number of eligible samples to include. Zero means no limit.",
    )
    parser.add_argument(
        "--out-batch",
        type=Path,
        default=None,
        help="Output batch manifest TSV with selected samples.",
    )
    parser.add_argument(
        "--out-missing",
        type=Path,
        default=None,
        help="Output blocker manifest TSV for rows that were not selected.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()

    if args.out_batch is None:
        args.out_batch = (
            PROJECT_ROOT
            / "step4_prn_validation"
            / "work/read_validation"
            / args.batch_label
            / "bp_prn_read_validation_batch.tsv"
        )
    if args.out_missing is None:
        args.out_missing = (
            PROJECT_ROOT
            / "step4_prn_validation"
            / "work/read_validation"
            / args.batch_label
            / "bp_prn_read_validation_missing_inputs.tsv"
        )

    subset_rows = load_tsv_rows(args.subset)

    selected_candidates: list[dict[str, str]] = []
    missing_rows: list[dict[str, str]] = []

    for row in subset_rows:
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        if not sample_id:
            continue

        link_status = normalize_text(row.get("raw_read_link_status", ""))
        resolved = resolve_paths(row, sample_id, args.reads_root, args.snippy_root)

        if link_status != "linked":
            missing_rows.append(
                to_manifest_row(
                    row,
                    resolved,
                    status="excluded_unlinked",
                    note="raw_read_link_status_not_linked",
                )
            )
            continue

        has_reads = bool(resolved.get("has_reads", False))
        has_bam = bool(resolved.get("has_bam", False))

        if has_reads and has_bam:
            selected_candidates.append(
                to_manifest_row(
                    row,
                    resolved,
                    status="selected",
                    note="linked_reads_qc_and_snippy_available",
                )
            )
            continue

        missing_rows.append(
            to_manifest_row(
                row,
                resolved,
                status="blocked_missing_local_inputs",
                note="linked_but_missing_reads_or_snippy_artifacts",
            )
        )

    selected_rows = selected_candidates
    if args.offset:
        for deferred in selected_rows[: args.offset]:
            deferred_row = dict(deferred)
            deferred_row["batch_status"] = "deferred_by_offset_limit"
            deferred_row["batch_note"] = "deferred_by_offset"
            missing_rows.append(deferred_row)
        selected_rows = selected_rows[args.offset :]

    if args.limit > 0 and len(selected_rows) > args.limit:
        for deferred in selected_rows[args.limit :]:
            deferred_row = dict(deferred)
            deferred_row["batch_status"] = "deferred_by_offset_limit"
            deferred_row["batch_note"] = "deferred_by_limit"
            missing_rows.append(deferred_row)
        selected_rows = selected_rows[: args.limit]

    write_tsv(args.out_batch, BATCH_COLUMNS, selected_rows)
    write_tsv(args.out_missing, MISSING_COLUMNS, missing_rows)

    blocked_count = sum(1 for row in missing_rows if row["batch_status"] == "blocked_missing_local_inputs")
    print(f"Wrote batch manifest: {args.out_batch}")
    print(f"Wrote blocker manifest: {args.out_missing}")
    print(f"Selected samples: {len(selected_rows)}")
    print(f"Blocked samples (missing local inputs): {blocked_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
