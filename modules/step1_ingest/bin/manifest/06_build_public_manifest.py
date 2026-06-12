#!/usr/bin/env python3
"""Build the pre-dedup public genome manifest for step1.

This script reads the cleaned step1 metadata table and writes a stable
candidate manifest with one row per assembly record before duplicate collapse.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Iterable

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from workflow.lib.project_paths import project_module_data_root


MISSING_TOKENS = {"", "missing", "unknown", "not applicable", "n/a"}

OUTPUT_COLUMNS = [
    "sample_id_canonical",
    "assembly_accession",
    "assembly_accession_root",
    "current_accession",
    "biosample_accession",
    "bioproject_accession",
    "sra_run_accession",
    "country",
    "year",
    "month",
    "week_key",
    "date_resolution",
    "source_database",
    "assembly_name",
    "assembly_level",
    "assembly_status",
    "assembly_refseq_category",
    "assembly_release_date",
    "collection_date_raw",
    "geo_raw",
    "host",
    "host_disease",
    "isolation_source",
    "strain",
    "isolate",
    "sequencing_tech",
    "total_sequence_length",
    "gc_percent",
    "n_contigs",
    "contig_n50",
]


def normalize_missing(value: str) -> str:
    return value.strip()


def is_missing(value: str) -> bool:
    return normalize_missing(value).lower() in MISSING_TOKENS


def accession_root(accession: str) -> str:
    value = normalize_missing(accession)
    if not value:
        return ""
    root = re.sub(r"\.\d+$", "", value)
    root = re.sub(r"^(GCA|GCF)_", "", root)
    return root


def canonical_sample_id(biosample_accession: str, assembly_accession: str) -> str:
    biosample = normalize_missing(biosample_accession)
    if biosample and not is_missing(biosample):
        return biosample
    return normalize_missing(assembly_accession)


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        raise ValueError(f"no rows found in {path}")
    return rows


def build_manifest_rows(rows: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    manifest_rows: list[dict[str, str]] = []
    for row in rows:
        assembly_accession = normalize_missing(row["Assembly Accession"])
        biosample_accession = normalize_missing(row["Assembly BioSample Accession"])
        manifest_rows.append(
            {
                "sample_id_canonical": canonical_sample_id(biosample_accession, assembly_accession),
                "assembly_accession": assembly_accession,
                "assembly_accession_root": accession_root(assembly_accession),
                "current_accession": normalize_missing(row["Current Accession"]),
                "biosample_accession": biosample_accession,
                "bioproject_accession": normalize_missing(row["Assembly BioProject Accession"]),
                "sra_run_accession": "",
                "country": normalize_missing(row["country"]),
                "year": normalize_missing(row["year"]),
                "month": normalize_missing(row["month"]),
                "week_key": normalize_missing(row["week_key"]),
                "date_resolution": normalize_missing(row["date_resolution"]),
                "source_database": normalize_missing(row["Source Database"]),
                "assembly_name": normalize_missing(row["Assembly Name"]),
                "assembly_level": normalize_missing(row["Assembly Level"]),
                "assembly_status": normalize_missing(row["Assembly Status"]),
                "assembly_refseq_category": normalize_missing(row["Assembly Refseq Category"]),
                "assembly_release_date": normalize_missing(row["Assembly Release Date"]),
                "collection_date_raw": normalize_missing(row["Assembly BioSample Collection date"]),
                "geo_raw": normalize_missing(row["Assembly BioSample Geographic location"]),
                "host": normalize_missing(row["Assembly BioSample Host"]),
                "host_disease": normalize_missing(row["Assembly BioSample Host disease"]),
                "isolation_source": normalize_missing(row["Assembly BioSample Isolation source"]),
                "strain": normalize_missing(row["Assembly BioSample Strain"]),
                "isolate": normalize_missing(row["Assembly BioSample Isolate"]),
                "sequencing_tech": normalize_missing(row["Assembly Sequencing Tech"]),
                "total_sequence_length": normalize_missing(row["Assembly Stats Total Sequence Length"]),
                "gc_percent": normalize_missing(row["Assembly Stats GC Percent"]),
                "n_contigs": normalize_missing(row["Assembly Stats Number of Contigs"]),
                "contig_n50": normalize_missing(row["Assembly Stats Contig N50"]),
            }
        )

    manifest_rows.sort(key=lambda item: (item["assembly_accession"], item["current_accession"]))
    return manifest_rows


def write_tsv(path: Path, rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_arg_parser() -> argparse.ArgumentParser:
    STEP1_DATA_ROOT = project_module_data_root("step1_ingest")
    parser = argparse.ArgumentParser(
        description="Build the pre-dedup public genome manifest from step1 cleaned metadata."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=STEP1_DATA_ROOT / "bp_metadata_clean.csv",
        help="Input cleaned metadata CSV (default: step1_ingest/bp_metadata_clean.csv).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_public_genome_manifest.tsv",
        help="Output manifest TSV (default: step1_ingest/outputs/bp_public_genome_manifest.tsv).",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    rows = load_rows(args.input)
    manifest_rows = build_manifest_rows(rows)
    write_tsv(args.output, manifest_rows)

    print(f"Wrote {args.output}")
    print(f"Rows: {len(manifest_rows)}")
    print("Header columns:", len(OUTPUT_COLUMNS))
    print("Note: duplicate records were intentionally preserved for later GC-03 resolution.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
