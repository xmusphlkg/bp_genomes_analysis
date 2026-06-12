#!/usr/bin/env python3
"""Build a reporting-era resolution worklist for blocked and backlog countries."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from ph_12_audit_reporting_era_coverage import audit_row
from ph_utils import (
    index_reporting_era_rows,
    match_reporting_era_row,
    normalize_text,
    read_delimited_rows,
    repo_root,
    project_workflow_root,
    write_tsv,
)


OUTPUT_COLUMNS = [
    "country_iso3",
    "country_name",
    "region_who",
    "current_match_type",
    "coverage_status",
    "target_status",
    "priority_group",
    "country_row_present",
    "country_row_confidence",
    "country_row_has_milestones",
    "country_row_has_interim_proxy_dates",
    "country_reporting_registry_sources",
    "country_all_registry_sources",
    "proxy_anchor_scope_type",
    "proxy_anchor_iso3",
    "proxy_primary_source_url",
    "proxy_secondary_source_url",
    "suggested_action",
    "suggested_search_hint",
]

STATUS_ORDER = {"blocked": 0, "global_backlog": 1, "complete": 2}


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


def build_country_row_lookup(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        if normalize_text(row.get("scope_type", "")).lower() != "country":
            continue
        iso3 = normalize_text(row.get("iso3", "")).upper()
        if iso3:
            lookup[iso3] = row
    return lookup


def build_registry_counters(rows: list[dict[str, str]]) -> tuple[Counter[str], Counter[str]]:
    reporting_counts: Counter[str] = Counter()
    all_counts: Counter[str] = Counter()
    for row in rows:
        iso3 = normalize_text(row.get("country_iso3", "")).upper()
        if not iso3:
            continue
        all_counts[iso3] += 1
        if normalize_text(row.get("source_domain", "")) == "reporting_era":
            reporting_counts[iso3] += 1
    return reporting_counts, all_counts


def row_has_milestones(row: dict[str, str]) -> bool:
    milestone_columns = (
        "pcr_lab_guideline_year",
        "reporting_case_definition_change_year",
        "surveillance_platform_change_year",
    )
    return any(normalize_text(row.get(column, "")) for column in milestone_columns)


def row_has_interim_proxy_dates(row: dict[str, str]) -> bool:
    return "interim" in normalize_text(row.get("coverage_note", "")).lower()


def suggested_action(
    *,
    coverage_status: str,
    country_row_present: bool,
    country_row_has_milestones: bool,
    country_row_has_interim_proxy_dates: bool,
) -> str:
    if coverage_status == "blocked":
        return "curate_national_reporting_row"
    if coverage_status == "global_backlog":
        return "queue_country_row_after_regional_pass"
    if not country_row_present or not country_row_has_milestones or country_row_has_interim_proxy_dates:
        return "extract_country_specific_milestones"
    return "maintain_country_row"


def suggested_search_hint(country_name: str, region_who: str, coverage_status: str) -> str:
    if coverage_status == "complete":
        return ""
    if region_who == "AMRO":
        return (
            f"{country_name} pertussis OR tos ferina OR coqueluche "
            "official ministry of health surveillance case definition notification"
        )
    if region_who == "EURO":
        return (
            f"{country_name} pertussis OR whooping cough "
            "official public health institute case definition notification"
        )
    return f"{country_name} pertussis official surveillance case definition notification"


def build_worklist_rows(
    master_rows: list[dict[str, str]],
    reporting_rows: list[dict[str, str]],
    registry_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    indexed_rows = index_reporting_era_rows(reporting_rows)
    country_rows = build_country_row_lookup(reporting_rows)
    reporting_counts, all_counts = build_registry_counters(registry_rows)

    output: list[dict[str, str]] = []
    for master_row in master_rows:
        iso3 = normalize_text(master_row.get("country_iso3", "")).upper()
        country_name = normalize_text(master_row.get("country_name", ""))
        region_who = normalize_text(master_row.get("region_who", "")).upper()
        matched_row, match_type = match_reporting_era_row(iso3, region_who, indexed_rows)
        coverage = audit_row(
            {
                "country_iso3": iso3,
                "region_who": region_who,
                "reporting_era_match_type": match_type,
            }
        )
        country_row = country_rows.get(iso3, {})
        country_row_present = bool(country_row)
        has_milestones = row_has_milestones(country_row)
        has_interim_proxy_dates = row_has_interim_proxy_dates(country_row)

        output.append(
            {
                "country_iso3": iso3,
                "country_name": country_name,
                "region_who": region_who,
                "current_match_type": match_type,
                "coverage_status": coverage["coverage_status"],
                "target_status": coverage["target_status"],
                "priority_group": coverage["priority_group"],
                "country_row_present": "1" if country_row_present else "0",
                "country_row_confidence": normalize_text(country_row.get("confidence", "")),
                "country_row_has_milestones": "1" if has_milestones else "0",
                "country_row_has_interim_proxy_dates": "1" if has_interim_proxy_dates else "0",
                "country_reporting_registry_sources": str(reporting_counts.get(iso3, 0)),
                "country_all_registry_sources": str(all_counts.get(iso3, 0)),
                "proxy_anchor_scope_type": normalize_text(matched_row.get("scope_type", "")),
                "proxy_anchor_iso3": normalize_text(matched_row.get("iso3", "")).upper(),
                "proxy_primary_source_url": normalize_text(matched_row.get("primary_source_url", "")),
                "proxy_secondary_source_url": normalize_text(matched_row.get("secondary_source_url", "")),
                "suggested_action": suggested_action(
                    coverage_status=coverage["coverage_status"],
                    country_row_present=country_row_present,
                    country_row_has_milestones=has_milestones,
                    country_row_has_interim_proxy_dates=has_interim_proxy_dates,
                ),
                "suggested_search_hint": suggested_search_hint(country_name, region_who, coverage["coverage_status"]),
            }
        )

    return sorted(
        output,
        key=lambda row: (
            STATUS_ORDER.get(row["coverage_status"], 99),
            normalize_text(row["region_who"]),
            normalize_text(row["country_iso3"]),
        ),
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a reporting-era resolution worklist from master outputs, raw curation, and source registry."
    )
    parser.add_argument(
        "--master",
        type=Path,
        default=WORKFLOW_DATA_ROOT / "epi" / "ap_exposure_index.tsv",
        help="Master country-year public-health table",
    )
    parser.add_argument(
        "--reporting-era",
        type=Path,
        default=repo_root() / "modules/public_health/inputs/raw/report_cases/pertussis_diagnosis_reporting_era_indicators.csv",
        help="Raw reporting-era curation CSV",
    )
    parser.add_argument(
        "--source-registry",
        type=Path,
        default=repo_root() / "modules/public_health/inputs/curation/public_health_source_registry.tsv",
        help="Canonical source registry TSV",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=WORKFLOW_DATA_ROOT / "epi" / "ph_reporting_era_resolution_worklist.tsv",
        help="Worklist TSV output path",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    worklist_rows = build_worklist_rows(
        load_master_rows(args.master),
        read_delimited_rows(args.reporting_era, delimiter=","),
        read_delimited_rows(args.source_registry, delimiter="\t"),
    )
    write_tsv(args.out, OUTPUT_COLUMNS, worklist_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
