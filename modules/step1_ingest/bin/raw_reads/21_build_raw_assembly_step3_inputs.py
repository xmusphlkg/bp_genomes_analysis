#!/usr/bin/env python3
"""Materialize Step3 input tables for QC-passed de novo raw-read assemblies.

Runtime environment:
  PROJECT_ENV_KEY: bio_tools
  PROJECT_ENV_NAME: pertussis-bio-tools
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from tempfile import NamedTemporaryFile


ROOT = Path(__file__).resolve().parents[4]
DEFAULT_RAW_QC_PASS = ROOT / "pipelines" / "bp_step1" / "outputs" / "bp_raw_read_assembly_qc_pass.tsv"
DEFAULT_TABLE_OUT = ROOT / "pipelines" / "bp_step1" / "outputs" / "bp_raw_read_step3_table.tsv"
DEFAULT_GENOME_PATHS_OUT = ROOT / "pipelines" / "bp_step1" / "outputs" / "bp_raw_read_step3_genome_paths.tsv"


TABLE_COLUMNS = [
    "genome_resolved_accession",
    "genome_status",
    "year",
    "country",
    "mlst_st",
]

GENOME_PATH_COLUMNS = [
    "resolved_accession",
    "fasta_path",
    "status",
]


def normalize_text(value: str | None) -> str:
    return (value or "").strip()


def pseudo_accession(run_accession: str) -> str:
    return f"RRASM_{normalize_text(run_accession)}"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", newline="", encoding="utf-8", delete=False, dir=path.parent) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
        temp_path = Path(handle.name)
    temp_path.replace(path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build Step3 inputs for QC-passed raw-read assemblies.")
    parser.add_argument("--raw-qc-pass", type=Path, default=DEFAULT_RAW_QC_PASS, help="QC-passed raw assembly TSV.")
    parser.add_argument("--table-out", type=Path, default=DEFAULT_TABLE_OUT, help="Step3 table TSV.")
    parser.add_argument(
        "--genome-paths-out",
        type=Path,
        default=DEFAULT_GENOME_PATHS_OUT,
        help="Step3 genome paths TSV.",
    )
    args = parser.parse_args(argv)

    rows = read_tsv(args.raw_qc_pass)
    table_rows: list[dict[str, str]] = []
    genome_path_rows: list[dict[str, str]] = []

    for row in rows:
        if normalize_text(row.get("qc_decision", "")) != "pass":
            continue
        run_accession = normalize_text(row.get("run_accession", ""))
        contigs_fasta = normalize_text(row.get("contigs_fasta", ""))
        if not run_accession or not contigs_fasta:
            continue
        accession = pseudo_accession(run_accession)
        table_rows.append(
            {
                "genome_resolved_accession": accession,
                "genome_status": "ok",
                "year": normalize_text(row.get("year", "")),
                "country": normalize_text(row.get("country", "")),
                "mlst_st": "",
            }
        )
        genome_path_rows.append(
            {
                "resolved_accession": accession,
                "fasta_path": contigs_fasta,
                "status": "ok",
            }
        )

    write_tsv(args.table_out, TABLE_COLUMNS, table_rows)
    write_tsv(args.genome_paths_out, GENOME_PATH_COLUMNS, genome_path_rows)
    print(f"Wrote Step3 table: {args.table_out}")
    print(f"Wrote genome paths: {args.genome_paths_out}")
    print(f"QC-passed raw assemblies exported: {len(table_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
