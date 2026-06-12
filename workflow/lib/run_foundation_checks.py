#!/usr/bin/env python3
"""Run the foundation checks and produce a consolidated readiness report."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from workflow.lib.project_paths import project_module_data_root, project_workflow_root


def run_json_command(command: list[str], output_path: Path, cwd: Path) -> dict:
    subprocess.run(command, cwd=str(cwd), check=True)
    with output_path.open(encoding="utf-8") as handle:
        return json.load(handle)


def main() -> None:
    base = Path(__file__).resolve().parents[1]

    parser = argparse.ArgumentParser(
        description="Run reads availability, vaccine-variable coverage, and validation-feasibility checks"
    )
    parser.add_argument(
        "--manifest",
        default=str(project_workflow_root() / "manifest" / "manifest.tsv"),
        help="Unified manifest TSV for the reads-availability check.",
    )
    parser.add_argument(
        "--outdir",
        default=str(project_workflow_root() / "checkpoints"),
        help="Directory for readiness JSON reports.",
    )
    parser.add_argument(
        "--out-runs",
        default=str(project_workflow_root() / "manifest" / "runs.tsv"),
        help="Output run table for the reads-availability check.",
    )
    parser.add_argument(
        "--prn-curation",
        default=str(project_module_data_root("public_health") / "inputs" / "curation" / "vaccine_formulation_curation.tsv"),
        help="Optional TSV with curated prn-in-vaccine rows for the vaccine-variable coverage check; accepts the richer formulation curation file.",
    )
    parser.add_argument("--min-reads-pct", type=float, default=30.0)
    parser.add_argument(
        "--skip-reads-availability",
        dest="skip_reads_availability",
        action="store_true",
        help="Skip the reads-availability check when it has already been run upstream.",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict] = {}

    if not args.skip_reads_availability:
        reads_report_path = outdir / "reads_availability_report.json"
        print("Running reads availability check...")
        try:
            results["reads_availability"] = run_json_command(
                [
                    sys.executable,
                    str(base.parent / "modules/step1_ingest/bin/raw_reads/24_trace_reads_availability.py"),
                    "--manifest",
                    args.manifest,
                    "--min-pct",
                    str(args.min_reads_pct),
                    "--out-runs",
                    args.out_runs,
                    "--out-report",
                    str(reads_report_path),
                ],
                reads_report_path,
                base,
            )
        except Exception as error:
            results["reads_availability"] = {"error": str(error)}
            print(f"  Reads availability error: {error}")

    vaccine_report_path = outdir / "vaccine_variable_coverage_report.json"
    print("\nRunning vaccine-variable coverage check...")
    try:
        results["vaccine_variable_coverage"] = run_json_command(
            [
                sys.executable,
                str(base.parent / "modules/public_health/bin/ph_09_assess_vaccine_variable_coverage.py"),
                "--prn-curation",
                args.prn_curation,
                "--out",
                str(vaccine_report_path),
            ],
            vaccine_report_path,
            base,
        )
    except Exception as error:
        results["vaccine_variable_coverage"] = {"error": str(error)}
        print(f"  Vaccine-variable coverage error: {error}")

    validation_report_path = outdir / "validation_feasibility_report.json"
    print("\nRunning validation-feasibility check...")
    try:
        results["validation_feasibility"] = run_json_command(
            [
                sys.executable,
                str(base.parent / "modules/step4_prn_validation/bin/step4_03b_assess_validation_feasibility.py"),
                "--out",
                str(validation_report_path),
            ],
            validation_report_path,
            base,
        )
    except Exception as error:
        results["validation_feasibility"] = {"error": str(error)}
        print(f"  Validation-feasibility error: {error}")

    out_path = outdir / "foundation_checks_report.json"
    with out_path.open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2, default=str)

    print(f"\n{'=' * 60}")
    print("FOUNDATION CHECK REPORT")
    print(f"{'=' * 60}")
    for check_name, check_data in results.items():
        if "error" in check_data:
            print(f"  {check_name}: ERROR — {check_data['error']}")
        else:
            print(f"  {check_name}: {check_data.get('decision', 'UNKNOWN')}")
    print(f"\nFull report: {out_path}")


if __name__ == "__main__":
    main()
