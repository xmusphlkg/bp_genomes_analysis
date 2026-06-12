#!/usr/bin/env python3
"""Extract the per-sequence composition chi-square table from an IQ-TREE log."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


OUTPUT_COLUMNS = [
    "row_index",
    "sequence_id",
    "gap_fraction_pct",
    "composition_test_status",
    "p_value_pct",
    "is_failed",
]

ROW_PATTERN = re.compile(r"^\s*(\d+)\s+(\S+)\s+([0-9.]+%)\s+(passed|failed)\s+([0-9.]+%)\s*$")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Extract the IQ-TREE composition chi-square result table.")
    parser.add_argument("--log", type=Path, required=True, help="IQ-TREE .log file path.")
    parser.add_argument("--out", type=Path, required=True, help="Output TSV path.")
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    rows: list[dict[str, str]] = []
    with args.log.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            match = ROW_PATTERN.match(raw_line.rstrip("\r\n"))
            if match is None:
                continue
            row_index, sequence_id, gap_fraction_pct, status, p_value_pct = match.groups()
            rows.append(
                {
                    "row_index": row_index,
                    "sequence_id": sequence_id,
                    "gap_fraction_pct": gap_fraction_pct,
                    "composition_test_status": status,
                    "p_value_pct": p_value_pct,
                    "is_failed": "True" if status == "failed" else "False",
                }
            )

    if not rows:
        raise ValueError(f"no composition chi-square rows found in {args.log}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())