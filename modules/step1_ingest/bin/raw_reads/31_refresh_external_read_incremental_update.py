#!/usr/bin/env python3
"""Refresh the external public-read discovery layer as an incremental update.

Runtime environment:
  PROJECT_ENV_KEY: bio_tools
  PROJECT_ENV_NAME: pertussis-bio-tools

This wrapper standardizes the refresh sequence:
1. pull a fresh ENA taxon catalog snapshot
2. rebuild the external-only candidate sample/run plans against the canonical
   manifest and existing planned/processed states
3. optionally rebuild the targeted country gap-fill subset

The wrapper does not download FASTQs. It only refreshes the discovery and
planning layer so repeated runs stay incremental and non-destructive.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from raw_read_utils import default_config_path, load_external_reads_config
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from workflow.lib.project_paths import project_module_data_root, project_workflow_root


ROOT = Path(__file__).resolve().parents[4]
SCRIPT_DIR = ROOT / "modules" / "step1_ingest" / "bin" / "raw_reads"
STEP1_DATA_ROOT = project_module_data_root("step1_ingest")


def run_command(command: list[str]) -> None:
    print("[run]", " ".join(str(part) for part in command))
    subprocess.run(command, check=True)


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Refresh the external public-read acquisition/cleaning plan as an incremental update."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=default_config_path(),
        help="TOML config file shared by the external read refresh chain.",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python interpreter to use for child refresh scripts.",
    )
    parser.add_argument(
        "--query",
        default="",
        help="Optional ENA taxon query override. Falls back to the config value when omitted.",
    )
    parser.add_argument(
        "--skip-targeted",
        action="store_true",
        help="Refresh only the catalog and external-only plan; skip targeted country subset regeneration.",
    )
    parser.add_argument(
        "--exclude-tsv",
        type=Path,
        action="append",
        default=[],
        help="Additional TSVs to treat as already planned or processed when rebuilding the external-only plan.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    config = load_external_reads_config(args.config)
    query = args.query or str(config.get("ena", {}).get("taxon_query", "tax_tree(520)"))

    catalog_output = STEP1_DATA_ROOT / "outputs" / "bp_ena_taxon_read_run_catalog.tsv"
    catalog_delta = STEP1_DATA_ROOT / "outputs" / "bp_ena_taxon_read_run_catalog_delta.tsv"
    catalog_summary = STEP1_DATA_ROOT / "outputs" / "bp_ena_taxon_read_run_catalog_refresh_summary.json"
    sample_output = STEP1_DATA_ROOT / "outputs" / "bp_external_raw_reads_only_samples.tsv"
    plan_output = STEP1_DATA_ROOT / "outputs" / "bp_external_raw_reads_only_plan.tsv"
    sample_delta = STEP1_DATA_ROOT / "outputs" / "bp_external_raw_reads_only_samples_incremental.tsv"
    plan_delta = STEP1_DATA_ROOT / "outputs" / "bp_external_raw_reads_only_plan_incremental.tsv"
    gapfill_summary = STEP1_DATA_ROOT / "outputs" / "bp_external_raw_reads_only_refresh_summary.json"

    fetch_cmd = [
        args.python,
        str(SCRIPT_DIR / "15_fetch_taxon_read_run_catalog.py"),
        "--config",
        str(args.config),
        "--query",
        query,
        "--output",
        str(catalog_output),
        "--delta-output",
        str(catalog_delta),
        "--summary-output",
        str(catalog_summary),
    ]
    run_command(fetch_cmd)

    gapfill_cmd = [
        args.python,
        str(SCRIPT_DIR / "16_build_external_gapfill.py"),
        "--config",
        str(args.config),
        "--catalog",
        str(catalog_output),
        "--manifest",
        str(project_workflow_root() / "manifest" / "manifest.tsv"),
        "--sample-out",
        str(sample_output),
        "--plan-out",
        str(plan_output),
        "--sample-delta-out",
        str(sample_delta),
        "--plan-delta-out",
        str(plan_delta),
        "--summary-out",
        str(gapfill_summary),
    ]
    for path in args.exclude_tsv:
        gapfill_cmd.extend(["--exclude-tsv", str(path)])
    run_command(gapfill_cmd)

    if not args.skip_targeted:
        targeted_cmd = [
            args.python,
            str(SCRIPT_DIR / "20_build_targeted_country_gapfill_subset.py"),
            "--config",
            str(args.config),
            "--plan",
            str(plan_output),
            "--samples",
            str(sample_output),
            "--out-plan",
            str(STEP1_DATA_ROOT / "outputs" / "bp_targeted_external_raw_reads_plan.tsv"),
            "--out-runs",
            str(STEP1_DATA_ROOT / "outputs" / "bp_targeted_external_raw_reads_runs.txt"),
            "--out-inventory",
            str(STEP1_DATA_ROOT / "outputs" / "bp_targeted_external_raw_reads_inventory.tsv"),
        ]
        run_command(targeted_cmd)

    catalog_payload = read_json(catalog_summary)
    gapfill_payload = read_json(gapfill_summary)

    print("\n=== Incremental Refresh Summary ===")
    if catalog_payload:
        print(
            "Catalog:"
            f" {catalog_payload.get('row_count_current', 0)} rows,"
            f" +{catalog_payload.get('new_runs', 0)} new,"
            f" {catalog_payload.get('changed_runs', 0)} changed,"
            f" {catalog_payload.get('removed_runs', 0)} removed"
        )
    if gapfill_payload:
        print(
            "Gapfill:"
            f" {gapfill_payload.get('external_only_samples', 0)} samples,"
            f" {gapfill_payload.get('run_plan_rows_emitted', 0)} run rows,"
            f" {gapfill_payload.get('incremental_sample_rows', 0)} new samples,"
            f" {gapfill_payload.get('incremental_run_rows', 0)} new runs"
        )
    print(f"Config: {args.config}")
    print(f"Catalog delta: {catalog_delta}")
    print(f"Plan delta: {plan_delta}")
    print(f"Summary JSON: {gapfill_summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
