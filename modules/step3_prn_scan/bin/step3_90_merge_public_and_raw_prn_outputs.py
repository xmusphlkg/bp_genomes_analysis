#!/usr/bin/env python3
"""Merge public-assembly and raw-read-assembly Step3 outputs for downstream Step4 calling."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from tempfile import NamedTemporaryFile


def normalize_text(value: str | None) -> str:
    return (value or "").strip()


def read_tsv(path: Path, *, allow_empty: bool = False) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        if allow_empty:
            return [], []
        raise FileNotFoundError(path)
    if path.stat().st_size == 0:
        if allow_empty:
            return [], []
        raise ValueError(f"empty file: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    if not fieldnames:
        if allow_empty:
            return [], []
        raise ValueError(f"no header found in {path}")
    return fieldnames, rows


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", newline="", encoding="utf-8", delete=False, dir=path.parent) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
        temp_path = Path(handle.name)
    temp_path.replace(path)


def merge_rows(
    public_rows: list[dict[str, str]],
    raw_rows: list[dict[str, str]],
    *,
    key_column: str,
) -> list[dict[str, str]]:
    merged: dict[str, dict[str, str]] = {}
    order: list[str] = []

    for source_rows in (public_rows, raw_rows):
        for row in source_rows:
            key = normalize_text(row.get(key_column, ""))
            if not key:
                continue
            if key not in merged:
                merged[key] = dict(row)
                order.append(key)
            else:
                merged[key].update({k: v for k, v in row.items() if normalize_text(v)})
    return [merged[key] for key in order]


def merge_tsv(
    public_path: Path,
    raw_path: Path,
    output_path: Path,
    *,
    key_column: str,
) -> int:
    public_fields, public_rows = read_tsv(public_path)
    raw_fields, raw_rows = read_tsv(raw_path, allow_empty=True)
    fieldnames = list(public_fields)
    for field in raw_fields:
        if field not in fieldnames:
            fieldnames.append(field)
    merged_rows = merge_rows(public_rows, raw_rows, key_column=key_column)
    normalized_rows = [{field: normalize_text(row.get(field, "")) for field in fieldnames} for row in merged_rows]
    write_tsv(output_path, fieldnames, normalized_rows)
    return len(normalized_rows)


def fasta_iter(path: Path):
    header = None
    seq_chunks: list[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", encoding="utf-8", delete=False, dir=path.parent) as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for idx in range(0, len(sequence), 80):
                handle.write(sequence[idx : idx + 80] + "\n")
        temp_path = Path(handle.name)
    temp_path.replace(path)


def merge_fasta(public_path: Path, raw_path: Path, output_path: Path) -> int:
    merged: dict[str, str] = {}
    order: list[str] = []
    for path in (public_path, raw_path):
        if not path.exists():
            continue
        for header, sequence in fasta_iter(path):
            key = normalize_text(header.split("|", 1)[0])
            if not key:
                continue
            if key not in merged:
                order.append(key)
            merged[key] = (header, sequence)
    records = [merged[key] for key in order]
    write_fasta(output_path, records)
    return len(records)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--public-calls", type=Path, required=True)
    parser.add_argument("--raw-calls", type=Path, required=True)
    parser.add_argument("--public-breakpoints", type=Path, required=True)
    parser.add_argument("--raw-breakpoints", type=Path, required=True)
    parser.add_argument("--public-gap-tsv", type=Path, required=True)
    parser.add_argument("--raw-gap-tsv", type=Path, required=True)
    parser.add_argument("--public-gap-fasta", type=Path, required=True)
    parser.add_argument("--raw-gap-fasta", type=Path, required=True)
    parser.add_argument("--out-calls", type=Path, required=True)
    parser.add_argument("--out-breakpoints", type=Path, required=True)
    parser.add_argument("--out-gap-tsv", type=Path, required=True)
    parser.add_argument("--out-gap-fasta", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    counts = {
        "calls": merge_tsv(
            args.public_calls,
            args.raw_calls,
            args.out_calls,
            key_column="genome_resolved_accession",
        ),
        "breakpoints": merge_tsv(
            args.public_breakpoints,
            args.raw_breakpoints,
            args.out_breakpoints,
            key_column="genome_resolved_accession",
        ),
        "gap_tsv": merge_tsv(
            args.public_gap_tsv,
            args.raw_gap_tsv,
            args.out_gap_tsv,
            key_column="genome_resolved_accession",
        ),
        "gap_fasta": merge_fasta(
            args.public_gap_fasta,
            args.raw_gap_fasta,
            args.out_gap_fasta,
        ),
    }
    for label, count in counts.items():
        print(f"{label}: {count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
