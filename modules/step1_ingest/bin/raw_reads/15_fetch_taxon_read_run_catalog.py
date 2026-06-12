#!/usr/bin/env python3
"""Fetch and refresh a taxon-level ENA read_run catalog for Bordetella pertussis.

Runtime environment:
  PROJECT_ENV_KEY: bio_tools
  PROJECT_ENV_NAME: pertussis-bio-tools

The refreshed catalog is now treated as a versioned snapshot:
- the full current live snapshot is written to ``--output``
- new/changed/removed rows relative to the previous snapshot are written to
  ``--delta-output``
- the previous snapshot is archived before replacement

This keeps the discovery layer reproducible while avoiding silent metadata loss.
"""

from __future__ import annotations

import argparse
import shutil
from datetime import datetime, timezone
from pathlib import Path

from raw_read_utils import (
    default_config_path,
    fetch_ena_rows,
    load_external_reads_config,
    normalize_text,
    read_tsv_if_exists as read_tsv,
    project_module_data_root,
    repo_root,
    write_json,
    write_tsv,
)


STEP1_DATA_ROOT = project_module_data_root("step1_ingest")


OUTPUT_COLUMNS = [
    "run_accession",
    "sample_accession",
    "secondary_sample_accession",
    "study_accession",
    "scientific_name",
    "collection_date",
    "country",
    "location",
    "library_source",
    "library_strategy",
    "library_layout",
    "instrument_platform",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "submitted_ftp",
]

DELTA_COLUMNS = ["change_status"] + OUTPUT_COLUMNS
def utc_now_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def archive_existing_output(output_path: Path, archive_dir: Path) -> Path | None:
    if not output_path.exists():
        return None
    archive_dir.mkdir(parents=True, exist_ok=True)
    archive_path = archive_dir / f"{output_path.stem}.{utc_now_stamp()}{output_path.suffix}"
    shutil.copy2(output_path, archive_path)
    return archive_path


def row_changed(previous: dict[str, str], current: dict[str, str]) -> bool:
    for column in OUTPUT_COLUMNS:
        if normalize_text(previous.get(column, "")) != normalize_text(current.get(column, "")):
            return True
    return False


def build_delta_rows(
    previous_rows: list[dict[str, str]],
    current_rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], dict[str, int]]:
    previous_by_run = {
        normalize_text(row.get("run_accession", "")): {column: normalize_text(row.get(column, "")) for column in OUTPUT_COLUMNS}
        for row in previous_rows
        if normalize_text(row.get("run_accession", ""))
    }
    current_by_run = {
        normalize_text(row.get("run_accession", "")): {column: normalize_text(row.get(column, "")) for column in OUTPUT_COLUMNS}
        for row in current_rows
        if normalize_text(row.get("run_accession", ""))
    }

    delta_rows: list[dict[str, str]] = []
    stats = {
        "new_runs": 0,
        "changed_runs": 0,
        "removed_runs": 0,
        "unchanged_runs": 0,
    }

    for run_accession in sorted(current_by_run):
        current = current_by_run[run_accession]
        previous = previous_by_run.get(run_accession)
        if previous is None:
            stats["new_runs"] += 1
            delta_rows.append({"change_status": "new", **current})
            continue
        if row_changed(previous, current):
            stats["changed_runs"] += 1
            delta_rows.append({"change_status": "changed", **current})
            continue
        stats["unchanged_runs"] += 1

    for run_accession in sorted(set(previous_by_run) - set(current_by_run)):
        stats["removed_runs"] += 1
        delta_rows.append({"change_status": "removed_from_live_query", **previous_by_run[run_accession]})

    return delta_rows, stats
def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fetch the ENA taxon-level read_run catalog for Bordetella pertussis with snapshot/delta tracking."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=default_config_path(),
        help="TOML config file controlling ENA query and refresh behavior.",
    )
    parser.add_argument(
        "--query",
        default="",
        help="Optional ENA query override. Falls back to the config value when omitted.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_ena_taxon_read_run_catalog.tsv",
        help="Full current snapshot TSV.",
    )
    parser.add_argument(
        "--delta-output",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_ena_taxon_read_run_catalog_delta.tsv",
        help="New/changed/removed rows relative to the previous snapshot.",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_ena_taxon_read_run_catalog_refresh_summary.json",
        help="JSON summary of the refresh/delta state.",
    )
    parser.add_argument(
        "--archive-dir",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "catalog_snapshots",
        help="Where to archive the previous full snapshot before replacement.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=None,
        help="Optional network timeout override for the ENA response stream.",
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=None,
        help="Optional retry-count override.",
    )
    parser.add_argument(
        "--sleep-seconds",
        type=float,
        default=None,
        help="Optional retry backoff override.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    config = load_external_reads_config(args.config)
    ena = config.get("ena", {})
    query = args.query or str(ena.get("taxon_query", "tax_tree(520)"))
    timeout_seconds = args.timeout_seconds or int(ena.get("timeout_seconds", 180))
    max_retries = args.max_retries or int(ena.get("max_retries", 3))
    sleep_seconds = args.sleep_seconds if args.sleep_seconds is not None else float(ena.get("sleep_seconds", 1.0))
    previous_rows = read_tsv(args.output)

    current_rows = fetch_ena_rows(
        result="read_run",
        fields=OUTPUT_COLUMNS,
        query=query,
        timeout_seconds=timeout_seconds,
        max_retries=max_retries,
        sleep_seconds=sleep_seconds,
    )
    current_rows = [row for row in current_rows if normalize_text(row.get("run_accession", ""))]
    current_rows.sort(key=lambda row: (row["study_accession"], row["sample_accession"], row["run_accession"]))
    if not current_rows:
        raise ValueError(f"ENA read_run query returned no rows for query: {args.query}")

    delta_rows, delta_stats = build_delta_rows(previous_rows, current_rows)
    archive_path = None
    if previous_rows and (
        delta_stats["new_runs"] or delta_stats["changed_runs"] or delta_stats["removed_runs"]
    ):
        archive_path = archive_existing_output(args.output, args.archive_dir)

    write_tsv(args.output, OUTPUT_COLUMNS, current_rows)
    write_tsv(args.delta_output, DELTA_COLUMNS, delta_rows)

    summary = {
        "config_path": str(args.config),
        "query": query,
        "refreshed_at_utc": datetime.now(timezone.utc).isoformat(),
        "output_path": str(args.output),
        "delta_output_path": str(args.delta_output),
        "archive_path": str(archive_path) if archive_path is not None else "",
        "row_count_current": len(current_rows),
        "row_count_previous": len(previous_rows),
        "unique_samples_current": len({row["sample_accession"] for row in current_rows if row["sample_accession"]}),
        "unique_studies_current": len({row["study_accession"] for row in current_rows if row["study_accession"]}),
        **delta_stats,
    }
    write_json(args.summary_output, summary)

    print(f"Wrote catalog: {args.output}")
    print(f"Wrote delta: {args.delta_output}")
    print(f"Wrote summary: {args.summary_output}")
    if archive_path is not None:
        print(f"Archived previous snapshot: {archive_path}")
    print(f"Query: {query}")
    print(f"Rows: {len(current_rows)}")
    print(f"Unique samples: {summary['unique_samples_current']}")
    print(
        "Delta:"
        f" +{delta_stats['new_runs']} new,"
        f" {delta_stats['changed_runs']} changed,"
        f" {delta_stats['removed_runs']} removed,"
        f" {delta_stats['unchanged_runs']} unchanged"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
