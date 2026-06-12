#!/usr/bin/env python3
"""Build a run-level raw-read download plan prioritized for fast cohort expansion."""

from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path

from raw_read_utils import (
    chunked,
    default_config_path,
    fetch_ena_rows,
    is_supported_short_read_platform,
    load_external_reads_config,
    normalize_text,
    parse_semicolon_ints,
    project_module_data_root,
    repo_root,
)


STEP1_DATA_ROOT = project_module_data_root("step1_ingest")
STEP4_DATA_ROOT = project_module_data_root("step4_prn_validation")

OUTPUT_COLUMNS = [
    "plan_row_id",
    "priority_tier",
    "priority_reason",
    "source_manifest",
    "analysis_cohort_id",
    "sample_id_canonical",
    "biosample_accession",
    "assembly_accession",
    "country",
    "year",
    "run_accession",
    "run_source",
    "raw_read_link_status",
    "raw_read_run_count",
    "ena_library_layout",
    "ena_instrument_platform",
    "ena_fastq_ftp",
    "ena_fastq_md5",
    "ena_fastq_bytes",
    "ena_submitted_ftp",
    "estimated_total_bytes",
    "run_metadata_status",
    "run_compatibility",
    "download_strategy",
]

RUN_METADATA_FIELDS = [
    "run_accession",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "submitted_ftp",
    "library_layout",
    "instrument_platform",
]

def parse_bool(value: str | None) -> bool:
    return normalize_text(value).lower() in {"1", "true", "yes", "y"}


def parse_runs(value: str | None) -> list[str]:
    raw = normalize_text(value)
    if not raw:
        return []
    return [part.strip() for part in raw.split(";") if part.strip()]


def parse_year(value: str | None) -> int:
    text = normalize_text(value)
    if not text:
        return 999999
    try:
        return int(text)
    except ValueError:
        return 999999


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fetch_run_metadata_rows(
    run_accessions: list[str],
    *,
    timeout_seconds: int,
    max_retries: int,
    sleep_seconds: float,
) -> list[dict[str, str]]:
    query = " OR ".join(f'run_accession="{run_accession}"' for run_accession in run_accessions)
    return fetch_ena_rows(
        result="read_run",
        fields=RUN_METADATA_FIELDS,
        query=query,
        timeout_seconds=timeout_seconds,
        max_retries=max_retries,
        sleep_seconds=sleep_seconds,
    )


def build_run_metadata_lookup(
    run_accessions: list[str],
    *,
    batch_size: int,
    timeout_seconds: int,
    max_retries: int,
    sleep_seconds: float,
) -> tuple[dict[str, dict[str, str]], set[str]]:
    lookup: dict[str, dict[str, str]] = {}
    failed_runs: set[str] = set()

    for batch in chunked(run_accessions, batch_size):
        try:
            rows = fetch_run_metadata_rows(
                batch,
                timeout_seconds=timeout_seconds,
                max_retries=max_retries,
                sleep_seconds=sleep_seconds,
            )
        except RuntimeError:
            failed_runs.update(batch)
            continue

        for row in rows:
            run_accession = normalize_text(row.get("run_accession", ""))
            if not run_accession:
                continue
            lookup[run_accession] = {
                "run_accession": run_accession,
                "fastq_ftp": normalize_text(row.get("fastq_ftp", "")),
                "fastq_md5": normalize_text(row.get("fastq_md5", "")),
                "fastq_bytes": normalize_text(row.get("fastq_bytes", "")),
                "submitted_ftp": normalize_text(row.get("submitted_ftp", "")),
                "library_layout": normalize_text(row.get("library_layout", "")),
                "instrument_platform": normalize_text(row.get("instrument_platform", "")),
            }

        time.sleep(sleep_seconds)

    return lookup, failed_runs
def classify_run_download(metadata: dict[str, str] | None, run_source: str, config: dict) -> tuple[str, str, str]:
    if metadata is None:
        if run_source == "SRA":
            return "batch_lookup_failed", "metadata_unresolved", "sra_toolkit_fallback"
        return "batch_lookup_failed", "metadata_unresolved", "skip_incompatible"

    fastq_urls = [item for item in metadata.get("fastq_ftp", "").split(";") if item.strip()]
    layout = metadata.get("library_layout", "").upper()
    platform = metadata.get("instrument_platform", "").upper()

    if not fastq_urls:
        if run_source == "SRA":
            return "resolved", "no_fastq_ftp", "sra_toolkit_fallback"
        return "resolved", "no_fastq_ftp", "skip_incompatible"
    if layout != "PAIRED":
        return "resolved", "not_paired", "skip_incompatible"
    if not is_supported_short_read_platform(platform, config):
        return "resolved", "not_supported_short_read_platform", "skip_incompatible"
    if len(fastq_urls) != 2:
        return "resolved", "unexpected_fastq_file_count", "skip_incompatible"
    return "resolved", "paired_short_read_fastq", "ena_fastq"


def priority_from_row(row: dict[str, str]) -> tuple[int, str]:
    cohort = normalize_text(row.get("analysis_cohort_id", "")).upper()
    if cohort == "D":
        return 1, "cohort_D_validation_pool"
    if cohort == "C":
        return 2, "cohort_C_country_year_gap_fill"
    if cohort == "B":
        return 3, "cohort_B_trend_support"
    if cohort == "A":
        return 4, "cohort_A_general_expansion"
    return 9, "unclassified_cohort"


def choose_better(existing: dict[str, str], candidate: dict[str, str]) -> dict[str, str]:
    existing_tier = int(existing["priority_tier"])
    candidate_tier = int(candidate["priority_tier"])
    if candidate_tier < existing_tier:
        return candidate
    if candidate_tier > existing_tier:
        return existing

    existing_year = parse_year(existing.get("year", ""))
    candidate_year = parse_year(candidate.get("year", ""))
    if candidate_year < existing_year:
        return candidate
    if candidate_year > existing_year:
        return existing

    existing_sample = existing.get("sample_id_canonical", "")
    candidate_sample = candidate.get("sample_id_canonical", "")
    if candidate_sample < existing_sample:
        return candidate
    return existing


def build_candidates(
    manifest_paths: list[Path],
    *,
    resolve_run_metadata: bool,
    max_runs: int | None,
    ena_batch_size: int,
    ena_timeout_seconds: int,
    ena_max_retries: int,
    ena_sleep_seconds: float,
    config: dict,
) -> tuple[list[dict[str, str]], int]:
    skipped_no_reads = 0
    run_to_best_row: dict[str, dict[str, str]] = {}

    for manifest_path in manifest_paths:
        rows = load_rows(manifest_path)
        for row in rows:
            if not parse_bool(row.get("raw_reads_available", "")):
                skipped_no_reads += 1
                continue

            priority_tier, priority_reason = priority_from_row(row)
            sra_runs = parse_runs(row.get("sra_run_accession", ""))
            ena_runs = parse_runs(row.get("ena_run_accession", ""))
            candidates: list[tuple[str, str]] = [(run, "SRA") for run in sra_runs]
            candidates.extend((run, "ENA") for run in ena_runs)

            for run_accession, run_source in candidates:
                candidate = {
                    "plan_row_id": "",
                    "priority_tier": str(priority_tier),
                    "priority_reason": priority_reason,
                    "source_manifest": manifest_path.name,
                    "analysis_cohort_id": normalize_text(row.get("analysis_cohort_id", "")),
                    "sample_id_canonical": normalize_text(row.get("sample_id_canonical", "")),
                    "biosample_accession": normalize_text(row.get("biosample_accession", "")),
                    "assembly_accession": normalize_text(row.get("assembly_accession", "")),
                    "country": normalize_text(row.get("country", "")),
                    "year": normalize_text(row.get("year", "")),
                    "run_accession": run_accession,
                    "run_source": run_source,
                    "raw_read_link_status": normalize_text(row.get("raw_read_link_status", "")),
                    "raw_read_run_count": normalize_text(row.get("raw_read_run_count", "")),
                }
                existing = run_to_best_row.get(run_accession)
                if existing is None:
                    run_to_best_row[run_accession] = candidate
                else:
                    run_to_best_row[run_accession] = choose_better(existing, candidate)

    candidates = list(run_to_best_row.values())
    candidates.sort(
        key=lambda item: (
            int(item["priority_tier"]),
            parse_year(item.get("year", "")),
            normalize_text(item.get("country", "")),
            normalize_text(item.get("sample_id_canonical", "")),
            normalize_text(item.get("run_accession", "")),
        )
    )

    for index, row in enumerate(candidates, start=1):
        row["plan_row_id"] = f"plan_{index:07d}"

    if max_runs is not None:
        candidates = candidates[:max_runs]

    if not resolve_run_metadata:
        for row in candidates:
            row["ena_library_layout"] = ""
            row["ena_instrument_platform"] = ""
            row["ena_fastq_ftp"] = ""
            row["ena_fastq_md5"] = ""
            row["ena_fastq_bytes"] = ""
            row["ena_submitted_ftp"] = ""
            row["estimated_total_bytes"] = ""
            row["run_metadata_status"] = "not_requested"
            row["run_compatibility"] = "metadata_not_requested"
            row["download_strategy"] = "sra_toolkit_fallback" if row["run_source"] == "SRA" else "skip_incompatible"
        return candidates, skipped_no_reads

    run_accessions = [row["run_accession"] for row in candidates if row.get("run_accession")]
    metadata_lookup, failed_runs = build_run_metadata_lookup(
        run_accessions,
        batch_size=ena_batch_size,
        timeout_seconds=ena_timeout_seconds,
        max_retries=ena_max_retries,
        sleep_seconds=ena_sleep_seconds,
    )

    for row in candidates:
        run_accession = row["run_accession"]
        metadata = metadata_lookup.get(run_accession)
        metadata_status, compatibility, strategy = classify_run_download(metadata, row["run_source"], config)
        if metadata is None and run_accession not in failed_runs:
            metadata_status = "not_found"

        fastq_bytes = metadata.get("fastq_bytes", "") if metadata is not None else ""
        estimated_total_bytes = str(sum(parse_semicolon_ints(fastq_bytes))) if fastq_bytes else ""

        row["ena_library_layout"] = metadata.get("library_layout", "") if metadata is not None else ""
        row["ena_instrument_platform"] = metadata.get("instrument_platform", "") if metadata is not None else ""
        row["ena_fastq_ftp"] = metadata.get("fastq_ftp", "") if metadata is not None else ""
        row["ena_fastq_md5"] = metadata.get("fastq_md5", "") if metadata is not None else ""
        row["ena_fastq_bytes"] = fastq_bytes
        row["ena_submitted_ftp"] = metadata.get("submitted_ftp", "") if metadata is not None else ""
        row["estimated_total_bytes"] = estimated_total_bytes
        row["run_metadata_status"] = metadata_status
        row["run_compatibility"] = compatibility
        row["download_strategy"] = strategy

    return candidates, skipped_no_reads


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def existing_manifest_paths(paths: list[Path]) -> list[Path]:
    existing = [path for path in paths if path.exists()]
    if existing:
        return existing
    missing = "\n".join(str(path) for path in paths)
    raise FileNotFoundError(f"No input manifests found. Checked:\n{missing}")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a run-level raw-read download plan for distributed execution."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=default_config_path(),
        help="TOML config file shared by the external raw-read planning chain.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        action="append",
        default=[],
        help="Input cohort manifest TSV. Repeat for multiple files. If omitted, defaults to cohort D then C.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=STEP4_DATA_ROOT / "inputs" / "bp_raw_reads_download_plan.tsv",
        help="Output run-level download plan TSV.",
    )
    parser.add_argument(
        "--max-runs",
        type=int,
        default=None,
        help="Optional cap on emitted run accessions (for dry runs).",
    )
    parser.add_argument(
        "--skip-run-metadata",
        action="store_true",
        help="Do not enrich runs with ENA fastq/layout metadata.",
    )
    parser.add_argument(
        "--ena-batch-size",
        type=int,
        default=None,
        help="Optional ENA metadata batch-size override.",
    )
    parser.add_argument(
        "--ena-timeout-seconds",
        type=int,
        default=None,
        help="Optional ENA metadata timeout override.",
    )
    parser.add_argument(
        "--ena-max-retries",
        type=int,
        default=None,
        help="Optional ENA metadata retry override.",
    )
    parser.add_argument(
        "--ena-sleep-seconds",
        type=float,
        default=None,
        help="Optional ENA metadata pause override.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    config = load_external_reads_config(args.config)
    ena = config.get("ena", {})

    default_manifests = [
        STEP1_DATA_ROOT / "outputs" / "bp_cohort_D_validation.tsv",
        STEP1_DATA_ROOT / "outputs" / "bp_cohort_C_country_year.tsv",
    ]
    manifest_paths = existing_manifest_paths(args.manifest or default_manifests)

    rows, skipped_no_reads = build_candidates(
        manifest_paths,
        resolve_run_metadata=not args.skip_run_metadata,
        max_runs=args.max_runs,
        ena_batch_size=args.ena_batch_size or int(ena.get("run_metadata_batch_size", 100)),
        ena_timeout_seconds=args.ena_timeout_seconds or int(ena.get("timeout_seconds", 30)),
        ena_max_retries=args.ena_max_retries or int(ena.get("max_retries", 3)),
        ena_sleep_seconds=args.ena_sleep_seconds if args.ena_sleep_seconds is not None else float(ena.get("sleep_seconds", 0.1)),
        config=config,
    )

    write_tsv(args.out, rows)

    unique_samples = len({row["sample_id_canonical"] for row in rows if row["sample_id_canonical"]})
    print(f"Wrote plan: {args.out}")
    print(f"Input manifests: {', '.join(str(path) for path in manifest_paths)}")
    print(f"Planned runs: {len(rows)}")
    print(f"Unique samples represented: {unique_samples}")
    print(f"Rows skipped because raw reads unavailable: {skipped_no_reads}")
    if rows and not args.skip_run_metadata:
        strategy_counts: dict[str, int] = {}
        for row in rows:
            strategy = row["download_strategy"]
            strategy_counts[strategy] = strategy_counts.get(strategy, 0) + 1
        strategy_summary = ", ".join(
            f"{strategy}={count}" for strategy, count in sorted(strategy_counts.items())
        )
        print(f"Download strategies: {strategy_summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
