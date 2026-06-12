#!/usr/bin/env python3
"""Filter an alignment by per-sequence missingness before downstream phylogeny."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REPORT_COLUMNS = [
    "sequence_id",
    "alignment_length",
    "missing_count",
    "missing_fraction",
    "keep",
    "reason",
]

MISSING_CHARS = set("-Nn?")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence_parts: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"): 
                if header is not None:
                    records.append((header, "".join(sequence_parts)))
                header = line[1:].strip()
                sequence_parts = []
            else:
                sequence_parts.append(line)
    if header is not None:
        records.append((header, "".join(sequence_parts)))
    if not records:
        raise ValueError(f"no FASTA records found in {path}")
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index : index + 80] + "\n")


def write_report(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Filter alignment records by gap/N missingness.")
    parser.add_argument("--alignment", type=Path, required=True, help="Input FASTA alignment.")
    parser.add_argument("--out-alignment", type=Path, required=True, help="Filtered FASTA alignment.")
    parser.add_argument("--out-report", type=Path, required=True, help="Per-sequence missingness report TSV.")
    parser.add_argument(
        "--max-missing-fraction",
        type=float,
        default=0.25,
        help="Maximum allowed missing fraction before excluding a sequence.",
    )
    parser.add_argument(
        "--always-keep",
        action="append",
        default=[],
        help="Sequence IDs to retain even if they exceed the missingness threshold.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    keep_ids = set(args.always_keep)
    input_records = read_fasta(args.alignment)

    kept_records: list[tuple[str, str]] = []
    report_rows: list[dict[str, str]] = []

    for sequence_id, sequence in input_records:
        alignment_length = len(sequence)
        if alignment_length == 0:
            raise ValueError(f"empty sequence found for {sequence_id}")
        missing_count = sum(1 for char in sequence if char in MISSING_CHARS)
        missing_fraction = missing_count / alignment_length
        keep = sequence_id in keep_ids or missing_fraction <= args.max_missing_fraction
        if keep:
            kept_records.append((sequence_id, sequence))

        reason = "always_keep" if sequence_id in keep_ids else ("within_threshold" if keep else "excluded_missingness")
        report_rows.append(
            {
                "sequence_id": sequence_id,
                "alignment_length": str(alignment_length),
                "missing_count": str(missing_count),
                "missing_fraction": f"{missing_fraction:.6f}",
                "keep": "True" if keep else "False",
                "reason": reason,
            }
        )

    if not kept_records:
        raise ValueError("missingness filter removed all sequences")

    write_fasta(args.out_alignment, kept_records)
    write_report(args.out_report, report_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())