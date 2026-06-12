#!/usr/bin/env python3
"""Standardize pertussis diagnosis/reporting era indicator curation."""

from __future__ import annotations

import argparse
from pathlib import Path

from ph_utils import (
    current_freeze_date,
    extract_min_year,
    normalize_text,
    read_delimited_rows,
    repo_root,
    project_module_data_root,
    split_source_tokens,
    write_tsv,
)


OUTPUT_COLUMNS = [
    "scope_type",
    "country_or_region",
    "iso3",
    "pcr_lab_guideline_year",
    "pcr_lab_guideline_year_min",
    "pcr_lab_guideline_exact_date",
    "reporting_case_definition_change_year",
    "reporting_case_definition_change_year_min",
    "reporting_case_definition_change_exact_date",
    "surveillance_platform_change_year",
    "surveillance_platform_change_year_min",
    "surveillance_platform_change_exact_date",
    "era_indicator_summary",
    "primary_source_title",
    "primary_source_url",
    "secondary_source_title",
    "secondary_source_url",
    "evidence_note",
    "confidence",
    "coverage_note",
    "source_file",
    "data_freeze_date",
]


PUBLIC_HEALTH_DATA_ROOT = project_module_data_root("public_health")


def load_rows(path: Path) -> list[dict[str, str]]:
    rows = read_delimited_rows(path, encoding="utf-8-sig", delimiter=",")
    output: list[dict[str, str]] = []
    for row in rows:
        pcr_year_min = extract_min_year(row.get("pcr_lab_guideline_year", ""))
        reporting_year_min = extract_min_year(row.get("reporting_case_definition_change_year", ""))
        platform_year_min = extract_min_year(row.get("surveillance_platform_change_year", ""))
        output.append(
            {
                "scope_type": normalize_text(row.get("scope_type", "")).lower(),
                "country_or_region": normalize_text(row.get("country_or_region", "")),
                "iso3": normalize_text(row.get("iso3", "")).upper(),
                "pcr_lab_guideline_year": normalize_text(row.get("pcr_lab_guideline_year", "")),
                "pcr_lab_guideline_year_min": "" if pcr_year_min is None else str(pcr_year_min),
                "pcr_lab_guideline_exact_date": normalize_text(row.get("pcr_lab_guideline_exact_date", "")),
                "reporting_case_definition_change_year": normalize_text(
                    row.get("reporting_case_definition_change_year", "")
                ),
                "reporting_case_definition_change_year_min": (
                    "" if reporting_year_min is None else str(reporting_year_min)
                ),
                "reporting_case_definition_change_exact_date": normalize_text(
                    row.get("reporting_case_definition_change_exact_date", "")
                ),
                "surveillance_platform_change_year": normalize_text(row.get("surveillance_platform_change_year", "")),
                "surveillance_platform_change_year_min": "" if platform_year_min is None else str(platform_year_min),
                "surveillance_platform_change_exact_date": normalize_text(
                    row.get("surveillance_platform_change_exact_date", "")
                ),
                "era_indicator_summary": normalize_text(row.get("era_indicator_summary", "")),
                "primary_source_title": normalize_text(row.get("primary_source_title", "")),
                "primary_source_url": normalize_text(row.get("primary_source_url", "")),
                "secondary_source_title": normalize_text(row.get("secondary_source_title", "")),
                "secondary_source_url": normalize_text(row.get("secondary_source_url", "")),
                "evidence_note": normalize_text(row.get("evidence_note", "")),
                "confidence": normalize_text(row.get("confidence", "")).lower(),
                "coverage_note": normalize_text(row.get("coverage_note", "")),
                "source_file": path.name,
                "data_freeze_date": current_freeze_date(),
            }
        )
    output.sort(key=lambda row: (row["scope_type"], row["iso3"], row["country_or_region"]))
    return output


def validate_registry_urls(rows: list[dict[str, str]], registry_path: Path) -> None:
    if not registry_path.exists():
        raise SystemExit(f"ERROR: canonical source registry not found: {registry_path}")
    registry_rows = read_delimited_rows(registry_path, delimiter="\t")
    registry_urls = {
        normalize_text(row.get("source_url", ""))
        for row in registry_rows
        if normalize_text(row.get("source_url", ""))
    }
    missing: list[str] = []
    for row in rows:
        for field in ("primary_source_url", "secondary_source_url"):
            for token in split_source_tokens(row.get(field, "")):
                if token not in registry_urls:
                    missing.append(token)
    if missing:
        missing_text = "\n".join(sorted(dict.fromkeys(missing)))
        raise SystemExit(f"ERROR: reporting-era source URLs missing from canonical registry:\n{missing_text}")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Standardize pertussis diagnosis/reporting era indicators.")
    parser.add_argument(
        "--input",
        type=Path,
        default=repo_root()
        / "modules/public_health/inputs/raw/report_cases/pertussis_diagnosis_reporting_era_indicators.csv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_reporting_era_indicators.tsv",
    )
    parser.add_argument(
        "--source-registry",
        type=Path,
        default=repo_root() / "modules/public_health/inputs/curation/public_health_source_registry.tsv",
        help="Canonical registry used to validate cited URLs",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    rows = load_rows(args.input)
    validate_registry_urls(rows, args.source_registry)
    write_tsv(args.out, OUTPUT_COLUMNS, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
