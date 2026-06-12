#!/usr/bin/env python3
"""Build an external raw-read-only gap-fill manifest and incremental run plan.

Runtime environment:
  PROJECT_ENV_KEY: bio_tools
  PROJECT_ENV_NAME: pertussis-bio-tools

This planner standardizes the discovery layer around one shared policy source:
1. compatibility rules come from ``config/external_reads.toml``
2. exclusion of already planned/processed rows comes from the same config
3. incremental deltas are computed against the previous emitted outputs
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path

from raw_read_utils import (
    collect_accessions,
    classify_ena_run,
    coverage_sets_from_rows,
    default_config_path,
    detect_run_source,
    discover_exclusion_paths,
    load_external_reads_config,
    normalize_country_for_output,
    normalize_text,
    parse_fastq_bytes,
    parse_year_from_date,
    read_tsv_if_exists,
    read_tsv_rows,
    project_module_data_root,
    project_workflow_root,
    repo_root,
    write_json,
    write_tsv,
)


STEP1_DATA_ROOT = project_module_data_root("step1_ingest")
WORKFLOW_DATA_ROOT = project_workflow_root()


RUN_PLAN_COLUMNS = [
    "plan_row_id",
    "priority_tier",
    "priority_reason",
    "source_manifest",
    "analysis_cohort_id",
    "sample_id_canonical",
    "biosample_accession",
    "secondary_sample_accession",
    "study_accession",
    "scientific_name",
    "assembly_accession",
    "country",
    "year",
    "collection_date",
    "ena_country_raw",
    "ena_location_raw",
    "run_accession",
    "run_source",
    "raw_read_link_status",
    "raw_read_run_count",
    "ena_library_layout",
    "ena_instrument_platform",
    "ena_library_source",
    "ena_library_strategy",
    "ena_fastq_ftp",
    "ena_fastq_md5",
    "ena_fastq_bytes",
    "ena_submitted_ftp",
    "estimated_total_bytes",
    "run_metadata_status",
    "run_compatibility",
    "download_strategy",
]

SAMPLE_COLUMNS = [
    "sample_id_canonical",
    "biosample_accession",
    "secondary_sample_accession",
    "study_accession",
    "scientific_name",
    "country_values",
    "collection_year_min",
    "collection_year_max",
    "run_count_total",
    "run_count_paired_short_read_fastq",
    "run_count_paired_illumina_fastq",
    "total_estimated_bytes",
    "run_accessions_all",
    "run_accessions_paired_short_read_fastq",
    "run_accessions_paired_illumina_fastq",
    "library_layouts",
    "instrument_platforms",
    "eligibility_status",
    "eligibility_reason",
]

SUMMARY_DEFAULT = STEP1_DATA_ROOT / "outputs" / "bp_external_raw_reads_only_refresh_summary.json"
RUN_COLUMNS = ("run_accession", "sra_run_accession", "ena_run_accession")
BIOSAMPLE_COLUMNS = ("biosample_accession", "sample_accession")
SAMPLE_COLUMNS_GENERIC = ("secondary_sample_accession", "sra_sample_accession", "ena_sample_accession")
SAMPLE_ID_COLUMNS = ("sample_id_canonical",)


def build_gapfill_rows(
    *,
    catalog_rows: list[dict[str, str]],
    manifest_rows: list[dict[str, str]],
    exclusion_paths: list[Path],
    include_incompatible_runs: bool,
    config: dict,
) -> tuple[list[dict[str, str]], list[dict[str, str]], Counter[str], dict[str, dict[str, int]]]:
    manifest_biosamples, manifest_samples, manifest_runs, manifest_sample_ids = coverage_sets_from_rows(
        manifest_rows,
        biosample_columns=BIOSAMPLE_COLUMNS,
        sample_columns=SAMPLE_COLUMNS_GENERIC,
        run_columns=RUN_COLUMNS,
        sample_id_columns=SAMPLE_ID_COLUMNS,
    )

    excluded_biosamples: set[str] = set()
    excluded_samples: set[str] = set()
    excluded_runs: set[str] = set()
    excluded_sample_ids: set[str] = set()
    exclusion_file_stats: dict[str, dict[str, int]] = {}

    for path in exclusion_paths:
        rows = read_tsv_if_exists(path)
        biosamples, samples, runs, sample_ids = coverage_sets_from_rows(
            rows,
            biosample_columns=BIOSAMPLE_COLUMNS,
            sample_columns=SAMPLE_COLUMNS_GENERIC,
            run_columns=RUN_COLUMNS,
            sample_id_columns=SAMPLE_ID_COLUMNS,
        )
        excluded_biosamples.update(biosamples)
        excluded_samples.update(samples)
        excluded_runs.update(runs)
        excluded_sample_ids.update(sample_ids)
        exclusion_file_stats[str(path)] = {
            "rows": len(rows),
            "biosamples": len(biosamples),
            "samples": len(samples),
            "runs": len(runs),
            "sample_ids": len(sample_ids),
        }

    exclusion_stats: Counter[str] = Counter()
    sample_to_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    sample_meta: dict[str, dict[str, object]] = {}
    kept_run_rows: list[dict[str, str]] = []

    for row in catalog_rows:
        run_accession = normalize_text(row.get("run_accession", ""))
        biosample_accession = normalize_text(row.get("sample_accession", ""))
        secondary_sample_accession = normalize_text(row.get("secondary_sample_accession", ""))
        sample_id_canonical = secondary_sample_accession or biosample_accession or run_accession

        if run_accession in manifest_runs:
            exclusion_stats["drop_manifest_run_accession"] += 1
            continue
        if biosample_accession and biosample_accession in manifest_biosamples:
            exclusion_stats["drop_manifest_biosample_accession"] += 1
            continue
        if biosample_accession and biosample_accession in manifest_samples:
            exclusion_stats["drop_manifest_primary_sample_accession"] += 1
            continue
        if secondary_sample_accession and secondary_sample_accession in manifest_samples:
            exclusion_stats["drop_manifest_secondary_sample_accession"] += 1
            continue
        if sample_id_canonical and sample_id_canonical in manifest_sample_ids:
            exclusion_stats["drop_manifest_sample_id"] += 1
            continue

        if run_accession in excluded_runs:
            exclusion_stats["drop_existing_planned_or_processed_run_accession"] += 1
            continue
        if biosample_accession and biosample_accession in excluded_biosamples:
            exclusion_stats["drop_existing_planned_or_processed_biosample"] += 1
            continue
        if biosample_accession and biosample_accession in excluded_samples:
            exclusion_stats["drop_existing_planned_or_processed_primary_sample"] += 1
            continue
        if secondary_sample_accession and secondary_sample_accession in excluded_samples:
            exclusion_stats["drop_existing_planned_or_processed_secondary_sample"] += 1
            continue
        if sample_id_canonical and sample_id_canonical in excluded_sample_ids:
            exclusion_stats["drop_existing_planned_or_processed_sample_id"] += 1
            continue

        run_metadata_status, run_compatibility, download_strategy = classify_ena_run(row, config)
        estimated_total_bytes = parse_fastq_bytes(row.get("fastq_bytes", ""))
        run_source = detect_run_source(run_accession)
        collection_date = normalize_text(row.get("collection_date", ""))
        country_raw = normalize_text(row.get("country", ""))
        location_raw = normalize_text(row.get("location", ""))
        country = normalize_country_for_output(country_raw)
        year = parse_year_from_date(collection_date)

        run_row = {
            "plan_row_id": "",
            "priority_tier": "5" if download_strategy == "ena_fastq" else "8",
            "priority_reason": "external_raw_reads_only_gapfill",
            "source_manifest": "bp_ena_taxon_read_run_catalog.tsv",
            "analysis_cohort_id": "E",
            "sample_id_canonical": sample_id_canonical,
            "biosample_accession": biosample_accession,
            "secondary_sample_accession": secondary_sample_accession,
            "study_accession": normalize_text(row.get("study_accession", "")),
            "scientific_name": normalize_text(row.get("scientific_name", "")) or "Bordetella pertussis",
            "assembly_accession": "",
            "country": country,
            "year": year,
            "collection_date": collection_date,
            "ena_country_raw": country_raw,
            "ena_location_raw": location_raw,
            "run_accession": run_accession,
            "run_source": run_source,
            "raw_read_link_status": "catalog_direct_taxon_read_run",
            "raw_read_run_count": "",
            "ena_library_layout": normalize_text(row.get("library_layout", "")),
            "ena_instrument_platform": normalize_text(row.get("instrument_platform", "")),
            "ena_library_source": normalize_text(row.get("library_source", "")),
            "ena_library_strategy": normalize_text(row.get("library_strategy", "")),
            "ena_fastq_ftp": normalize_text(row.get("fastq_ftp", "")),
            "ena_fastq_md5": normalize_text(row.get("fastq_md5", "")),
            "ena_fastq_bytes": normalize_text(row.get("fastq_bytes", "")),
            "ena_submitted_ftp": normalize_text(row.get("submitted_ftp", "")),
            "estimated_total_bytes": str(estimated_total_bytes) if estimated_total_bytes else "",
            "run_metadata_status": run_metadata_status,
            "run_compatibility": run_compatibility,
            "download_strategy": download_strategy,
        }

        sample_to_rows[sample_id_canonical].append(run_row)
        if download_strategy == "ena_fastq" or include_incompatible_runs:
            kept_run_rows.append(run_row)
            exclusion_stats["keep_run_rows"] += 1
        else:
            exclusion_stats["keep_sample_incompatible_only"] += 1

        meta = sample_meta.setdefault(
            sample_id_canonical,
            {
                "study_accession": "",
                "scientific_name": run_row["scientific_name"],
                "secondary_sample_accession": secondary_sample_accession,
                "countries": set(),
                "years": [],
            },
        )
        if not meta["study_accession"]:
            meta["study_accession"] = run_row["study_accession"]
        if not meta["secondary_sample_accession"]:
            meta["secondary_sample_accession"] = secondary_sample_accession
        if country:
            meta["countries"].add(country)
        if year.isdigit():
            meta["years"].append(int(year))

    sample_rows: list[dict[str, str]] = []
    for sample_id, rows in sorted(sample_to_rows.items(), key=lambda item: item[0]):
        all_runs = sorted(row["run_accession"] for row in rows)
        eligible_runs = sorted(row["run_accession"] for row in rows if row["download_strategy"] == "ena_fastq")
        total_bytes = sum(int(row["estimated_total_bytes"]) for row in rows if row["estimated_total_bytes"])
        layouts = sorted({row["ena_library_layout"] for row in rows if row["ena_library_layout"]})
        platforms = sorted({row["ena_instrument_platform"] for row in rows if row["ena_instrument_platform"]})
        first = rows[0]
        meta = sample_meta.get(sample_id, {})
        years = sorted(meta.get("years", []))
        sample_rows.append(
            {
                "sample_id_canonical": sample_id,
                "biosample_accession": first["biosample_accession"],
                "secondary_sample_accession": str(meta.get("secondary_sample_accession", "")),
                "study_accession": str(meta.get("study_accession", "")),
                "scientific_name": str(meta.get("scientific_name", "Bordetella pertussis")),
                "country_values": ";".join(sorted(meta.get("countries", set()))),
                "collection_year_min": str(years[0]) if years else "",
                "collection_year_max": str(years[-1]) if years else "",
                "run_count_total": str(len(all_runs)),
                "run_count_paired_short_read_fastq": str(len(eligible_runs)),
                "run_count_paired_illumina_fastq": str(len(eligible_runs)),
                "total_estimated_bytes": str(total_bytes) if total_bytes else "",
                "run_accessions_all": ";".join(all_runs),
                "run_accessions_paired_short_read_fastq": ";".join(eligible_runs),
                "run_accessions_paired_illumina_fastq": ";".join(eligible_runs),
                "library_layouts": ";".join(layouts),
                "instrument_platforms": ";".join(platforms),
                "eligibility_status": "eligible" if eligible_runs else "incompatible_only",
                "eligibility_reason": (
                    "paired_short_read_fastq_present" if eligible_runs else "no_supported_paired_short_read_fastq"
                ),
            }
        )

    eligible_run_counts_by_sample = Counter(
        row["sample_id_canonical"] for row in kept_run_rows if row["download_strategy"] == "ena_fastq"
    )
    for index, row in enumerate(
        sorted(
            kept_run_rows,
            key=lambda current: (
                int(current["priority_tier"]),
                -int(current["estimated_total_bytes"] or "0"),
                current["study_accession"],
                current["sample_id_canonical"],
                current["run_accession"],
            ),
        ),
        start=1,
    ):
        row["plan_row_id"] = f"gapfill_{index:07d}"
        row["raw_read_run_count"] = str(eligible_run_counts_by_sample[row["sample_id_canonical"]])

    return sample_rows, kept_run_rows, exclusion_stats, exclusion_file_stats


def new_rows_by_key(
    rows: list[dict[str, str]],
    previous_rows: list[dict[str, str]],
    key_columns: tuple[str, ...],
) -> list[dict[str, str]]:
    previous_keys = {
        tuple(normalize_text(row.get(column, "")) for column in key_columns)
        for row in previous_rows
    }
    out: list[dict[str, str]] = []
    for row in rows:
        key = tuple(normalize_text(row.get(column, "")) for column in key_columns)
        if key not in previous_keys:
            out.append(row)
    return out


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build an external raw-read-only sample manifest and compatible run-level gap-fill plan "
            "by subtracting the canonical retained manifest plus existing plan/processing states "
            "from the ENA taxon read_run catalog."
        )
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=default_config_path(),
        help="TOML config file controlling compatibility and exclusion behavior.",
    )
    parser.add_argument(
        "--catalog",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_ena_taxon_read_run_catalog.tsv",
        help="Input ENA taxon read_run catalog TSV.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=WORKFLOW_DATA_ROOT / "manifest" / "manifest.tsv",
        help="Canonical retained public genome manifest TSV.",
    )
    parser.add_argument(
        "--sample-out",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_external_raw_reads_only_samples.tsv",
        help="Output sample-level external raw-read-only manifest TSV.",
    )
    parser.add_argument(
        "--plan-out",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_external_raw_reads_only_plan.tsv",
        help="Output run-level external raw-read-only plan TSV.",
    )
    parser.add_argument(
        "--sample-delta-out",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_external_raw_reads_only_samples_incremental.tsv",
        help="Sample-level incremental TSV containing only newly surfaced candidates.",
    )
    parser.add_argument(
        "--plan-delta-out",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_external_raw_reads_only_plan_incremental.tsv",
        help="Run-level incremental TSV containing only newly surfaced candidates.",
    )
    parser.add_argument(
        "--summary-out",
        type=Path,
        default=repo_root() / SUMMARY_DEFAULT,
        help="JSON summary of retained, excluded, and incremental counts.",
    )
    parser.add_argument(
        "--exclude-tsv",
        type=Path,
        action="append",
        default=[],
        help="Additional TSVs whose run/sample accessions should be treated as already planned or processed.",
    )
    parser.add_argument(
        "--include-incompatible-runs",
        action="store_true",
        help="Emit incompatible runs in the run-level plan instead of sample-only summaries.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    config = load_external_reads_config(args.config)

    previous_sample_rows = read_tsv_if_exists(args.sample_out)
    previous_plan_rows = read_tsv_if_exists(args.plan_out)

    catalog_rows = read_tsv_rows(args.catalog)
    manifest_rows = read_tsv_rows(args.manifest)
    exclusion_paths = discover_exclusion_paths(repo_root(), args.exclude_tsv, config)

    sample_rows, plan_rows, exclusion_stats, exclusion_file_stats = build_gapfill_rows(
        catalog_rows=catalog_rows,
        manifest_rows=manifest_rows,
        exclusion_paths=exclusion_paths,
        include_incompatible_runs=args.include_incompatible_runs,
        config=config,
    )

    sample_delta_rows = new_rows_by_key(sample_rows, previous_sample_rows, ("sample_id_canonical",))
    plan_delta_rows = new_rows_by_key(plan_rows, previous_plan_rows, ("run_accession",))

    write_tsv(args.sample_out, SAMPLE_COLUMNS, sample_rows)
    write_tsv(args.plan_out, RUN_PLAN_COLUMNS, plan_rows)
    write_tsv(args.sample_delta_out, SAMPLE_COLUMNS, sample_delta_rows)
    write_tsv(args.plan_delta_out, RUN_PLAN_COLUMNS, plan_delta_rows)

    eligible_sample_count = sum(1 for row in sample_rows if row["eligibility_status"] == "eligible")
    summary = {
        "config_path": str(args.config),
        "catalog_path": str(args.catalog),
        "manifest_path": str(args.manifest),
        "sample_output_path": str(args.sample_out),
        "plan_output_path": str(args.plan_out),
        "sample_delta_output_path": str(args.sample_delta_out),
        "plan_delta_output_path": str(args.plan_delta_out),
        "compatibility_policy": config.get("compatibility", {}),
        "exclusion_paths": [str(path) for path in exclusion_paths],
        "exclusion_file_stats": exclusion_file_stats,
        "external_only_samples": len(sample_rows),
        "eligible_paired_short_read_samples": eligible_sample_count,
        "run_plan_rows_emitted": len(plan_rows),
        "incremental_sample_rows": len(sample_delta_rows),
        "incremental_run_rows": len(plan_delta_rows),
        "exclusion_stats": dict(sorted(exclusion_stats.items())),
        "top_studies_by_sample_count": sorted(
            Counter(row["study_accession"] for row in sample_rows if row["study_accession"]).items(),
            key=lambda item: (-item[1], item[0]),
        )[:25],
        "unique_country_annotations": len(collect_accessions(sample_rows, ("country_values",))),
    }
    write_json(args.summary_out, summary)

    print(f"Wrote sample manifest: {args.sample_out}")
    print(f"Wrote run plan: {args.plan_out}")
    print(f"Wrote incremental sample manifest: {args.sample_delta_out}")
    print(f"Wrote incremental run plan: {args.plan_delta_out}")
    print(f"Wrote summary: {args.summary_out}")
    print(f"External-only samples: {len(sample_rows)}")
    print(f"Eligible paired short-read samples: {eligible_sample_count}")
    print(f"Run-plan rows emitted: {len(plan_rows)}")
    print(f"Incremental new sample rows: {len(sample_delta_rows)}")
    print(f"Incremental new run rows: {len(plan_delta_rows)}")
    print("Exclusion/retention summary:")
    for key, value in sorted(exclusion_stats.items()):
        print(f"  {key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
