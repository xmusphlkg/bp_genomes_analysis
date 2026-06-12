#!/usr/bin/env python3
"""Aggregate per-sample QUAST reports into a single workflow QC table."""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd


OUTPUT_COLUMNS = [
    "sample_id_canonical",
    "assembly_name",
    "total_length",
    "n_contigs",
    "contig_n50",
    "gc_percent",
    "report_path",
]


def normalize_metric_name(value: str) -> str:
    return (
        (value or "")
        .strip()
        .lower()
        .replace("# ", "")
        .replace("(", "")
        .replace(")", "")
        .replace("%", "pct")
        .replace(">= 0 bp", "")
        .replace("  ", " ")
    )


def parse_quast_report(report_path: Path) -> dict[str, str]:
    with report_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty QUAST report: {report_path}")

    sample_id = report_path.parent.name
    metrics: dict[str, str] = {
        "sample_id_canonical": sample_id,
        "assembly_name": "",
        "total_length": "",
        "n_contigs": "",
        "contig_n50": "",
        "gc_percent": "",
        "report_path": str(report_path),
    }

    if rows and len(rows[0]) >= 2 and rows[0][0] == "Assembly":
        metrics["assembly_name"] = rows[0][1]

    for row in rows:
        if len(row) < 2:
            continue
        metric_name = normalize_metric_name(row[0])
        value = row[1].strip()
        if metric_name.startswith("assembly") and not metrics["assembly_name"]:
            metrics["assembly_name"] = value
        elif metric_name.startswith("total length"):
            metrics["total_length"] = value
        elif metric_name.startswith("contigs"):
            metrics["n_contigs"] = value
        elif metric_name.startswith("n50"):
            metrics["contig_n50"] = value
        elif metric_name.startswith("gc"):
            metrics["gc_percent"] = value

    return metrics


def aggregate_reports(report_paths: list[str]) -> pd.DataFrame:
    records = [parse_quast_report(Path(path)) for path in report_paths]
    frame = pd.DataFrame.from_records(records, columns=OUTPUT_COLUMNS)
    return frame.sort_values(["sample_id_canonical"]).reset_index(drop=True)


if "snakemake" in globals():
    output_table = Path(snakemake.output.table)
    aggregated = aggregate_reports(list(snakemake.input.quast))
    output_table.parent.mkdir(parents=True, exist_ok=True)
    aggregated.to_csv(output_table, sep="\t", index=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Aggregate QUAST reports")
    parser.add_argument("reports", nargs="+", help="Per-sample QUAST report.tsv files")
    parser.add_argument("--out-table", required=True, help="Output TSV path")
    arguments = parser.parse_args()

    aggregated = aggregate_reports(arguments.reports)
    output_table = Path(arguments.out_table)
    output_table.parent.mkdir(parents=True, exist_ok=True)
    aggregated.to_csv(output_table, sep="\t", index=False)
