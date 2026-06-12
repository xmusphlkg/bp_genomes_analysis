#!/usr/bin/env python3
"""Remove recombinant columns from a full alignment using Gubbins GFF output.

This script converts a `core.full.aln` alignment plus a Gubbins recombinant-region
GFF into a downstream tree-ready `recomb_filtered.aln` that preserves invariant sites
outside recombinant intervals.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    with path.open("w") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for idx in range(0, len(sequence), 80):
                handle.write(sequence[idx : idx + 80] + "\n")


def load_intervals(gff_path: Path) -> list[tuple[int, int]]:
    intervals: list[tuple[int, int]] = []
    with gff_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue
            if end < start:
                start, end = end, start
            intervals.append((start, end))
    return intervals


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end + 1:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def mask_alignment_columns(alignment_path: Path, intervals: list[tuple[int, int]], output_path: Path) -> dict:
    alignment = read_fasta(alignment_path)
    if not alignment:
        raise ValueError(f"No sequences found in alignment: {alignment_path}")

    aln_len = len(alignment[0][1])
    for _, sequence in alignment:
        if len(sequence) != aln_len:
            raise ValueError("Alignment sequences are not all the same length")

    keep = [True] * aln_len

    masked_columns = 0
    for start, end in intervals:
        # GFF is 1-based inclusive.
        start_idx = max(start - 1, 0)
        end_idx = min(end, aln_len)
        for idx in range(start_idx, end_idx):
            if keep[idx]:
                keep[idx] = False
                masked_columns += 1

    kept_indices = [idx for idx, keep_col in enumerate(keep) if keep_col]
    filtered_records = []
    for header, sequence in alignment:
        filtered_sequence = "".join(sequence[idx] for idx in kept_indices)
        filtered_records.append((header, filtered_sequence))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(filtered_records, output_path)

    return {
        "alignment_length_before": aln_len,
        "alignment_length_after": len(kept_indices),
        "masked_columns": masked_columns,
        "n_sequences": len(filtered_records),
    }


def run(alignment_path: Path, gff_path: Path, output_path: Path, summary_path: Path | None = None) -> dict:
    merged_intervals = merge_intervals(load_intervals(gff_path))
    summary = mask_alignment_columns(alignment_path, merged_intervals, output_path)
    summary["gff_path"] = str(gff_path)
    summary["alignment_path"] = str(alignment_path)
    summary["output_path"] = str(output_path)
    summary["n_intervals_raw"] = len(load_intervals(gff_path))
    summary["n_intervals_merged"] = len(merged_intervals)
    summary["masked_fraction"] = round(
        summary["masked_columns"] / max(summary["alignment_length_before"], 1), 6
    )
    if summary_path is not None:
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with summary_path.open("w") as handle:
            json.dump(summary, handle, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Mask recombinant columns from a full alignment")
    parser.add_argument("--alignment", required=True, help="Input core.full.aln FASTA")
    parser.add_argument("--gff", required=True, help="Gubbins recombinant-region GFF")
    parser.add_argument("--output", required=True, help="Filtered alignment path")
    parser.add_argument("--summary", default="", help="Optional summary JSON path")
    args = parser.parse_args()

    summary = run(
        alignment_path=Path(args.alignment),
        gff_path=Path(args.gff),
        output_path=Path(args.output),
        summary_path=Path(args.summary) if args.summary else None,
    )
    print(json.dumps(summary, indent=2))


if "snakemake" in dir():
    summary_path = Path(str(snakemake.output.filtered_aln)).with_suffix(".mask_summary.json")
    run(
        alignment_path=Path(str(snakemake.input.full_aln)),
        gff_path=Path(str(snakemake.input.gubbins_gff)),
        output_path=Path(str(snakemake.output.filtered_aln)),
        summary_path=summary_path,
    )


if __name__ == "__main__":
    main()