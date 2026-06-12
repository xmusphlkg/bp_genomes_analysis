#!/usr/bin/env python3
"""Build a targeted external raw-read intake plan for high-value country windows.

Runtime environment:
  PROJECT_ENV_KEY: bio_tools
  PROJECT_ENV_NAME: pertussis-bio-tools

The targeted subset is now derived only from the cleaned external gap-fill plan.
That removes a second ENA query pass and guarantees that the targeted selection
inherits the exact same compatibility, exclusion, and dedup policy as the
upstream planning step.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

from raw_read_utils import (
    default_config_path,
    load_external_reads_config,
    load_targeted_targets,
    normalize_text,
    parse_fastq_bytes,
    parse_year_from_date,
    read_tsv_rows,
    read_tsv_with_header,
    project_module_data_root,
    repo_root,
    write_tsv,
)


ROOT = repo_root()
STEP1_DATA_ROOT = project_module_data_root("step1_ingest")
DEFAULT_PLAN = STEP1_DATA_ROOT / "outputs" / "bp_external_raw_reads_only_plan.tsv"
DEFAULT_SAMPLES = STEP1_DATA_ROOT / "outputs" / "bp_external_raw_reads_only_samples.tsv"
DEFAULT_OUT_PLAN = STEP1_DATA_ROOT / "outputs" / "bp_targeted_external_raw_reads_plan.tsv"
DEFAULT_OUT_RUNS = STEP1_DATA_ROOT / "outputs" / "bp_targeted_external_raw_reads_runs.txt"
DEFAULT_OUT_INVENTORY = STEP1_DATA_ROOT / "outputs" / "bp_targeted_external_raw_reads_inventory.tsv"


def candidate_rows(plan_rows: list[dict[str, str]], target: dict[str, str]) -> list[dict[str, str]]:
    study = normalize_text(target.get("study_accession", ""))
    country_filter = normalize_text(target.get("country_filter", "")).casefold()
    year_min = normalize_text(target.get("year_min", ""))
    year_max = normalize_text(target.get("year_max", ""))
    out: list[dict[str, str]] = []

    for row in plan_rows:
        if normalize_text(row.get("study_accession", "")) != study:
            continue

        country_raw = normalize_text(row.get("ena_country_raw", "")) or normalize_text(row.get("country", ""))
        if country_filter and country_raw.casefold() != country_filter:
            continue

        year = normalize_text(row.get("year", "")) or parse_year_from_date(row.get("collection_date", ""))
        if year_min and year.isdigit() and int(year) < int(year_min):
            continue
        if year_max and year.isdigit() and int(year) > int(year_max):
            continue

        merged = dict(row)
        merged["country"] = normalize_text(row.get("country", "")) or normalize_text(target.get("country_name", ""))
        merged["country_iso3"] = normalize_text(target.get("country_iso3", ""))
        merged["year"] = year
        merged["study_accession"] = study
        merged["target_label"] = normalize_text(target.get("target_label", ""))
        merged["estimated_total_bytes_numeric"] = str(
            int(normalize_text(row.get("estimated_total_bytes", "")) or "0")
            or parse_fastq_bytes(row.get("ena_fastq_bytes", ""))
        )
        out.append(merged)
    return out


def sort_candidates(rows: list[dict[str, str]], sort_mode: str) -> list[dict[str, str]]:
    if sort_mode == "smallest_bytes":
        return sorted(
            rows,
            key=lambda row: (
                int(normalize_text(row.get("estimated_total_bytes_numeric", "")) or "0"),
                normalize_text(row.get("collection_date", "")),
                normalize_text(row.get("run_accession", "")),
            ),
        )
    if sort_mode == "largest_bytes":
        return sorted(
            rows,
            key=lambda row: (
                -int(normalize_text(row.get("estimated_total_bytes_numeric", "")) or "0"),
                normalize_text(row.get("collection_date", "")),
                normalize_text(row.get("run_accession", "")),
            ),
        )
    return sorted(rows, key=lambda row: normalize_text(row.get("run_accession", "")))


def select_target_rows(rows: list[dict[str, str]], target: dict[str, str]) -> list[dict[str, str]]:
    rows = sort_candidates(rows, normalize_text(target.get("sort_mode", "")))
    max_samples = int(str(target.get("max_samples", "") or "0"))
    max_per_year = int(str(target.get("max_per_year", "") or "0"))
    if max_per_year > 0:
        year_counts: dict[str, int] = defaultdict(int)
        selected: list[dict[str, str]] = []
        for row in rows:
            year = normalize_text(row.get("year", "")) or "unknown"
            if year_counts[year] >= max_per_year:
                continue
            selected.append(row)
            year_counts[year] += 1
        rows = selected
    if max_samples > 0:
        rows = rows[:max_samples]
    for idx, row in enumerate(rows, start=1):
        row["priority_tier"] = "1"
        row["priority_reason"] = f"targeted_country_gapfill::{target['target_label']}"
        row["analysis_cohort_id"] = "E_TARGETED"
        row["target_rank_within_label"] = str(idx)
    return rows


def build_inventory_rows(selected_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    inventory: list[dict[str, str]] = []
    for row in selected_rows:
        inventory.append(
            {
                "target_label": normalize_text(row.get("target_label", "")),
                "country_iso3": normalize_text(row.get("country_iso3", "")),
                "study_accession": normalize_text(row.get("study_accession", "")),
                "sample_id_canonical": normalize_text(row.get("sample_id_canonical", "")),
                "biosample_accession": normalize_text(row.get("biosample_accession", "")),
                "run_accession": normalize_text(row.get("run_accession", "")),
                "year": normalize_text(row.get("year", "")),
                "collection_date": normalize_text(row.get("collection_date", "")),
                "estimated_total_bytes": normalize_text(row.get("estimated_total_bytes_numeric", "")),
                "ena_country_raw": normalize_text(row.get("ena_country_raw", "")),
                "ena_location_raw": normalize_text(row.get("ena_location_raw", "")),
            }
        )
    return inventory


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build targeted external gap-fill intake plan.")
    parser.add_argument(
        "--config",
        type=Path,
        default=default_config_path(),
        help="TOML config file containing targeted country-window definitions.",
    )
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN, help="External-only run-level plan TSV.")
    parser.add_argument(
        "--samples",
        type=Path,
        default=DEFAULT_SAMPLES,
        help="External-only sample-level TSV. Retained for interface stability.",
    )
    parser.add_argument("--out-plan", type=Path, default=DEFAULT_OUT_PLAN, help="Targeted run-level plan TSV.")
    parser.add_argument("--out-runs", type=Path, default=DEFAULT_OUT_RUNS, help="One-run-per-line list.")
    parser.add_argument(
        "--out-inventory",
        type=Path,
        default=DEFAULT_OUT_INVENTORY,
        help="Inventory TSV with country/year metadata for selected runs.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    config = load_external_reads_config(args.config)
    targets = load_targeted_targets(config)
    if not targets:
        raise ValueError(f"No targeted targets defined in config: {args.config}")

    plan_fields, plan_rows = read_tsv_with_header(args.plan)
    _ = read_tsv_rows(args.samples)

    selected_rows: list[dict[str, str]] = []
    for target in targets:
        candidates = candidate_rows(plan_rows, target)
        chosen = select_target_rows(candidates, target)
        selected_rows.extend(chosen)
        print(
            f"[target] {target['target_label']}: selected {len(chosen)} of {len(candidates)} candidates",
            file=sys.stderr,
        )

    seen_runs: set[str] = set()
    deduped_rows: list[dict[str, str]] = []
    for row in selected_rows:
        run_accession = normalize_text(row.get("run_accession", ""))
        if not run_accession or run_accession in seen_runs:
            continue
        seen_runs.add(run_accession)
        deduped_rows.append(row)

    output_fields = list(plan_fields)
    for extra in [
        "country_iso3",
        "collection_date",
        "ena_country_raw",
        "ena_location_raw",
        "study_accession",
        "target_label",
        "target_rank_within_label",
        "estimated_total_bytes_numeric",
    ]:
        if extra not in output_fields:
            output_fields.append(extra)

    inventory_rows = build_inventory_rows(deduped_rows)
    write_tsv(args.out_plan, output_fields, deduped_rows)
    write_tsv(
        args.out_inventory,
        [
            "target_label",
            "country_iso3",
            "study_accession",
            "sample_id_canonical",
            "biosample_accession",
            "run_accession",
            "year",
            "collection_date",
            "estimated_total_bytes",
            "ena_country_raw",
            "ena_location_raw",
        ],
        inventory_rows,
    )

    args.out_runs.parent.mkdir(parents=True, exist_ok=True)
    with args.out_runs.open("w", encoding="utf-8") as handle:
        for row in deduped_rows:
            handle.write(f"{normalize_text(row.get('run_accession', ''))}\n")

    print(f"Wrote targeted plan: {args.out_plan}")
    print(f"Wrote targeted inventory: {args.out_inventory}")
    print(f"Wrote run list: {args.out_runs}")
    print(f"Selected runs: {len(deduped_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
