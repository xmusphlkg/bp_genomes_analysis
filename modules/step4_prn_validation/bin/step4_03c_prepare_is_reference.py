#!/usr/bin/env python3
"""Build a shell-safe IS reference FASTA for Step4 read validation tools."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from step4_02_scan_prn_mechanisms import STEP4_DATA_ROOT, load_reference_fasta, write_tsv


MAP_COLUMNS = [
    "reference_id",
    "is_element_name",
    "source_accession",
    "original_header",
    "sanitized_header",
    "sequence_length_bp",
]


def sanitize_header(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    value = re.sub(r"_+", "_", value).strip("_")
    if not value:
        raise ValueError("reference header sanitized to an empty string")
    return value


def wrap_fasta(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            handle.write(f"{wrap_fasta(sequence)}\n")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Materialize a shell-safe IS reference FASTA for ISMapper-compatible Step4 read validation."
    )
    parser.add_argument(
        "--in-fasta",
        type=Path,
        default=STEP4_DATA_ROOT / "references" / "is_elements" / "bp_is_reference.fasta",
        help="Curated IS reference FASTA from step4_01_build_is_reference.py.",
    )
    parser.add_argument(
        "--out-fasta",
        type=Path,
        default=STEP4_DATA_ROOT / "references" / "is_elements" / "bp_is_reference.shell_safe.fasta",
        help="Shell-safe FASTA with sanitized headers for ISMapper.",
    )
    parser.add_argument(
        "--out-map",
        type=Path,
        default=STEP4_DATA_ROOT / "references" / "is_elements" / "bp_is_reference_shell_safe_map.tsv",
        help="Header translation table between canonical and sanitized FASTA identifiers.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    reference_records = load_reference_fasta(args.in_fasta)

    fasta_records: list[tuple[str, str]] = []
    map_rows: list[dict[str, str]] = []
    seen_headers: set[str] = set()

    for record in reference_records:
        original_header = "|".join(
            [
                record["reference_id"],
                record["is_element_name"],
                record["source_accession"],
            ]
        )
        sanitized_header = sanitize_header(original_header)
        if sanitized_header in seen_headers:
            raise ValueError(f"duplicate sanitized FASTA header generated: {sanitized_header}")
        seen_headers.add(sanitized_header)
        fasta_records.append((sanitized_header, record["sequence"]))
        map_rows.append(
            {
                "reference_id": record["reference_id"],
                "is_element_name": record["is_element_name"],
                "source_accession": record["source_accession"],
                "original_header": original_header,
                "sanitized_header": sanitized_header,
                "sequence_length_bp": str(len(record["sequence"])),
            }
        )

    write_fasta(args.out_fasta, fasta_records)
    write_tsv(args.out_map, MAP_COLUMNS, map_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
