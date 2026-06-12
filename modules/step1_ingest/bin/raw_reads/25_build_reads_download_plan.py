#!/usr/bin/env python3
"""Build a normalized reads download plan from the retained manifest.

Runtime environment:
  PROJECT_ENV_KEY: bio_tools
  PROJECT_ENV_NAME: pertussis-bio-tools

This finalizer reuses the same compatibility outputs produced upstream and makes
the per-sample run selection policy explicit in ``config/external_reads.toml``.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from raw_read_utils import (
    default_config_path,
    download_selection_metric,
    load_external_reads_config,
    metric_value,
    normalize_text,
    read_tsv_rows,
    project_module_data_root,
    project_workflow_root,
    repo_root,
    write_tsv,
)


ROOT = repo_root()
STEP4_DATA_ROOT = project_module_data_root("step4_prn_validation")
MANIFEST = project_workflow_root() / "manifest" / "manifest.tsv"
EXISTING_PLAN = STEP4_DATA_ROOT / "inputs" / "bp_raw_reads_download_plan.tsv"
OUTPUT_DIR = project_workflow_root() / "reads_plan"
COMPATIBLE_LABELS = {"paired_short_read_fastq", "paired_illumina_fastq"}


def choose_better(existing: dict[str, str], candidate: dict[str, str], metric: str) -> dict[str, str]:
    existing_metric = metric_value(existing, metric)
    candidate_metric = metric_value(candidate, metric)
    if metric == "largest_estimated_total_bytes":
        if candidate_metric > existing_metric:
            return candidate
        if candidate_metric < existing_metric:
            return existing
    elif metric == "smallest_estimated_total_bytes":
        if candidate_metric < existing_metric:
            return candidate
        if candidate_metric > existing_metric:
            return existing

    existing_tier = int(normalize_text(existing.get("priority_tier", "")) or "999")
    candidate_tier = int(normalize_text(candidate.get("priority_tier", "")) or "999")
    if candidate_tier < existing_tier:
        return candidate
    if candidate_tier > existing_tier:
        return existing

    existing_run = normalize_text(existing.get("run_accession", ""))
    candidate_run = normalize_text(candidate.get("run_accession", ""))
    if candidate_run < existing_run:
        return candidate
    return existing


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a normalized per-sample reads download plan from the current workflow manifest."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=default_config_path(),
        help="TOML config file containing download-plan selection policy.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=MANIFEST,
        help="Workflow manifest whose retained samples define the download universe.",
    )
    parser.add_argument(
        "--plan",
        type=Path,
        default=EXISTING_PLAN,
        help="Run-level raw-read plan TSV to normalize into one chosen run per sample.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Directory where the normalized plan and summary will be written.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    config = load_external_reads_config(args.config)
    selection_metric = download_selection_metric(config)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = read_tsv_rows(args.manifest)
    manifest_samples = {
        normalize_text(row.get("sample_id_canonical", "")): row
        for row in manifest_rows
        if normalize_text(row.get("sample_id_canonical", ""))
    }
    plan_rows = read_tsv_rows(args.plan)

    skipped_reasons = Counter()
    compatible: list[dict[str, str]] = []
    seen_samples: set[str] = set()

    for row in plan_rows:
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        compatibility = normalize_text(row.get("run_compatibility", ""))
        strategy = normalize_text(row.get("download_strategy", ""))

        if sample_id not in manifest_samples:
            skipped_reasons["not_in_manifest"] += 1
            continue
        if compatibility not in COMPATIBLE_LABELS:
            skipped_reasons[f"incompatible:{compatibility or 'missing'}"] += 1
            continue
        if strategy != "ena_fastq":
            skipped_reasons[f"non_downloadable_strategy:{strategy or 'missing'}"] += 1
            continue

        compatible.append(row)
        seen_samples.add(sample_id)

    best_per_sample: dict[str, dict[str, str]] = {}
    for row in compatible:
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        if sample_id not in best_per_sample:
            best_per_sample[sample_id] = row
            continue
        best_per_sample[sample_id] = choose_better(best_per_sample[sample_id], row, selection_metric)

    out_cols = [
        "sample_id_canonical",
        "assembly_accession",
        "run_accession",
        "ena_library_layout",
        "ena_instrument_platform",
        "ena_fastq_ftp",
        "ena_fastq_md5",
        "estimated_total_bytes",
        "download_strategy",
        "priority_tier",
    ]

    plan_path = args.output_dir / "reads_download_plan.tsv"
    output_rows = [
        {column: row.get(column, "") for column in out_cols}
        for sample_id, row in sorted(best_per_sample.items())
    ]
    write_tsv(
        plan_path,
        out_cols,
        output_rows,
    )

    total_bytes = sum(int(normalize_text(row.get("estimated_total_bytes", "")) or "0") for row in best_per_sample.values())
    samples_no_reads = set(manifest_samples) - seen_samples

    summary_lines = [
        "=== Reads Download Plan Summary ===",
        f"Config: {args.config}",
        f"Selection metric: {selection_metric}",
        f"Manifest samples: {len(manifest_samples)}",
        f"Samples with compatible paired short-read reads: {len(best_per_sample)}",
        f"Samples without compatible reads: {len(samples_no_reads)}",
        f"Total estimated download size: {total_bytes / (1024**3):.1f} GB",
        f"Average per-sample: {total_bytes / max(len(best_per_sample), 1) / (1024**2):.1f} MB",
        "",
        "=== Skip Reasons ===",
    ]
    for reason, count in sorted(skipped_reasons.items(), key=lambda item: (-item[1], item[0])):
        summary_lines.append(f"  {reason}: {count}")

    summary_path = args.output_dir / "reads_plan_summary.txt"
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    print("\n".join(summary_lines))
    print(f"\nOutput: {plan_path}")
    print(f"Summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
