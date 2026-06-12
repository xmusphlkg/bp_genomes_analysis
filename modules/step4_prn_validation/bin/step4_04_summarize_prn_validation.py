#!/usr/bin/env python3
"""Summarize PRN read-validation outcomes by mechanism class and confidence tier."""

from __future__ import annotations

import argparse
import os
from collections import Counter, defaultdict
from pathlib import Path

from step4_02_scan_prn_mechanisms import (
    load_tsv_rows,
    normalize_text,
    project_module_data_root,
    write_tsv,
)


STEP4_DATA_ROOT = project_module_data_root("step4_prn_validation")


OUTPUT_COLUMNS = [
    "summary_level",
    "prn_mechanism_call",
    "prn_call_confidence",
    "read_validation_status",
    "n_samples",
    "fraction_within_group",
    "major_read_support_class",
    "n_major_read_support_class",
    "major_targeted_locus_assembly_status",
    "n_major_targeted_locus_assembly_status",
    "notes",
]


def summarize_group(
    rows: list[dict[str, str]],
    *,
    summary_level: str,
    mechanism: str,
    confidence: str,
) -> list[dict[str, str]]:
    status_counts = Counter(normalize_text(row.get("read_validation_status", "")) for row in rows)
    support_counts = Counter(normalize_text(row.get("read_support_class", "")) for row in rows)
    locus_counts = Counter(normalize_text(row.get("targeted_locus_assembly_status", "")) for row in rows)
    total = len(rows)
    major_support, major_support_n = support_counts.most_common(1)[0]
    major_locus, major_locus_n = locus_counts.most_common(1)[0]

    output_rows: list[dict[str, str]] = []
    for status, n_samples in sorted(status_counts.items()):
        output_rows.append(
            {
                "summary_level": summary_level,
                "prn_mechanism_call": mechanism,
                "prn_call_confidence": confidence,
                "read_validation_status": status,
                "n_samples": str(n_samples),
                "fraction_within_group": f"{n_samples / total:.6f}",
                "major_read_support_class": major_support,
                "n_major_read_support_class": str(major_support_n),
                "major_targeted_locus_assembly_status": major_locus,
                "n_major_targeted_locus_assembly_status": str(major_locus_n),
                "notes": f"group_total={total}",
            }
        )
    return output_rows


def build_summary_rows(
    validation_rows: list[dict[str, str]],
    mechanism_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    mechanism_by_sample = {
        normalize_text(row.get("sample_id_canonical", "")): row for row in mechanism_rows
    }

    joined_rows: list[dict[str, str]] = []
    for validation_row in validation_rows:
        sample_id = normalize_text(validation_row.get("sample_id_canonical", ""))
        mechanism_row = mechanism_by_sample.get(sample_id, {})
        joined_rows.append(
            {
                **validation_row,
                "prn_mechanism_call": normalize_text(
                    validation_row.get("prn_mechanism_call", "") or mechanism_row.get("prn_mechanism_call", "")
                ),
                "prn_call_confidence": normalize_text(mechanism_row.get("prn_call_confidence", "")),
            }
        )

    by_mechanism: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_mechanism_confidence: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in joined_rows:
        mechanism = normalize_text(row.get("prn_mechanism_call", "")) or "unknown"
        confidence = normalize_text(row.get("prn_call_confidence", "")) or "unknown"
        by_mechanism[mechanism].append(row)
        by_mechanism_confidence[(mechanism, confidence)].append(row)

    output_rows: list[dict[str, str]] = []
    for mechanism in sorted(by_mechanism):
        output_rows.extend(
            summarize_group(
                by_mechanism[mechanism],
                summary_level="mechanism",
                mechanism=mechanism,
                confidence="all",
            )
        )
    for mechanism, confidence in sorted(by_mechanism_confidence):
        output_rows.extend(
            summarize_group(
                by_mechanism_confidence[(mechanism, confidence)],
                summary_level="mechanism_confidence",
                mechanism=mechanism,
                confidence=confidence,
            )
        )
    return output_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Summarize PRN read-validation outcomes by mechanism class and confidence tier."
    )
    parser.add_argument(
        "--validation",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_read_validation.tsv",
        help="Raw validation table from VAL-02.",
    )
    parser.add_argument(
        "--mechanism-calls",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_mechanism_calls.tsv",
        help="Assembly-side PRN mechanism call table.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_validation_summary.tsv",
        help="Validation summary output TSV.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    validation_rows = load_tsv_rows(args.validation)
    mechanism_rows = load_tsv_rows(args.mechanism_calls)
    summary_rows = build_summary_rows(validation_rows, mechanism_rows)
    write_tsv(args.out, OUTPUT_COLUMNS, summary_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
