#!/usr/bin/env python3
"""Audit reporting-era country coverage for direct, blocked, and global-backlog status."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from ph_utils import normalize_text, repo_root, write_tsv
from ph_utils import project_workflow_root


OUTPUT_COLUMNS = [
    "country_iso3",
    "region_who",
    "current_match_type",
    "target_status",
    "coverage_status",
    "blocker_reason",
    "priority_group",
]


WORKFLOW_DATA_ROOT = project_workflow_root()


def load_master_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_country: dict[str, dict[str, str]] = {}
    for row in rows:
        iso3 = normalize_text(row.get("country_iso3", "")).upper()
        if iso3 and iso3 not in by_country:
            by_country[iso3] = row
    return [by_country[iso3] for iso3 in sorted(by_country)]


def audit_row(row: dict[str, str]) -> dict[str, str]:
    iso3 = normalize_text(row.get("country_iso3", "")).upper()
    region = normalize_text(row.get("region_who", "")).upper()
    match_type = normalize_text(row.get("reporting_era_match_type", ""))

    if match_type == "country_direct":
        coverage_status = "complete"
        target_status = "country_direct"
        blocker_reason = ""
    elif match_type == "country_alias_proxy":
        coverage_status = "blocked"
        target_status = "country_direct_required"
        blocker_reason = "country alias fallback still in use; country-level row not yet curated"
    elif match_type == "regional_proxy":
        coverage_status = "blocked"
        target_status = "country_direct_required"
        blocker_reason = "regional proxy retained because no country-level national source row was curated in this pass"
    elif match_type == "global_proxy":
        coverage_status = "global_backlog"
        target_status = "global_backlog"
        blocker_reason = "outside regional-first scope; queued for later country-level curation"
    else:
        coverage_status = "blocked"
        target_status = "country_direct_required"
        blocker_reason = "no reporting-era match resolved from current curation"

    if coverage_status == "global_backlog":
        priority_group = "global_backlog"
    elif region in {"AMRO", "EURO"}:
        priority_group = "regional_first"
    elif coverage_status == "complete":
        priority_group = "country_direct_existing"
    else:
        priority_group = "regional_or_manual_blocker"

    return {
        "country_iso3": iso3,
        "region_who": region,
        "current_match_type": match_type,
        "target_status": target_status,
        "coverage_status": coverage_status,
        "blocker_reason": blocker_reason,
        "priority_group": priority_group,
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Audit reporting-era coverage status from master country-year output.")
    parser.add_argument(
        "--master",
        type=Path,
        default=WORKFLOW_DATA_ROOT / "epi" / "ap_exposure_index.tsv",
        help="Master country-year public-health table",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=WORKFLOW_DATA_ROOT / "epi" / "ph_reporting_era_coverage_audit.tsv",
        help="Coverage audit output TSV",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    rows = [audit_row(row) for row in load_master_rows(args.master)]
    write_tsv(args.out, OUTPUT_COLUMNS, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
