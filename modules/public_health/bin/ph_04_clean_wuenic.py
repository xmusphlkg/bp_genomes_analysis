#!/usr/bin/env python3
"""Clean WUENIC DTP coverage export into a country-year table."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from ph_utils import (
    current_freeze_date,
    find_first_input,
    load_country_name_map,
    normalize_country,
    normalize_text,
    parse_export_date_from_name,
    parse_float,
    read_xlsx_sheet_rows,
    repo_root,
    project_module_data_root,
    rows_to_dicts,
    write_tsv,
)


OUTPUT_COLUMNS = [
    "country_iso3",
    "country_name",
    "year",
    "dtp1_coverage_admin",
    "dtp1_coverage_official",
    "dtp3_coverage_admin",
    "dtp3_coverage_official",
    "dtp3_coverage",
    "booster_antigen",
    "booster_coverage",
    "source_name",
    "source_url",
    "source_release_date",
    "data_freeze_date",
    "notes",
]


PUBLIC_HEALTH_DATA_ROOT = project_module_data_root("public_health")


def default_input_path() -> Path:
    return find_first_input(
        repo_root() / "modules/public_health" / "inputs" / "raw" / "wuenic",
        (".xlsx",),
    )


def build_output_rows(input_path: Path, country_map_path: Path) -> list[dict[str, str]]:
    country_map = load_country_name_map(country_map_path)
    raw_rows = rows_to_dicts(read_xlsx_sheet_rows(input_path, sheet_name="Sheet1"))
    grouped: dict[tuple[str, str], dict[str, str]] = defaultdict(dict)

    for raw_row in raw_rows:
        if normalize_text(raw_row.get("GROUP", "")) != "COUNTRIES":
            continue
        raw_country = normalize_text(raw_row.get("NAME", ""))
        country_row = normalize_country(raw_country, country_map)
        if country_row.get("match_status") != "normalized":
            continue

        year = normalize_text(raw_row.get("YEAR", ""))
        antigen = normalize_text(raw_row.get("ANTIGEN", ""))
        category = normalize_text(raw_row.get("COVERAGE_CATEGORY", "")).casefold()
        coverage = parse_float(raw_row.get("COVERAGE", ""))
        if not year or coverage is None:
            continue

        key = (country_row["country_iso3"], year)
        grouped[key].update(
            {
                "country_iso3": country_row["country_iso3"],
                "country_name": country_row["normalized_country_name"],
                "year": year,
                "source_name": "WHO/UNICEF WUENIC",
                "source_url": "https://immunizationdata.who.int/dashboard",
                "source_release_date": parse_export_date_from_name(input_path),
                "data_freeze_date": current_freeze_date(),
            }
        )

        if antigen == "DTPCV1":
            grouped[key][f"dtp1_coverage_{category}"] = f"{coverage:.2f}"
        elif antigen == "DTPCV3":
            grouped[key][f"dtp3_coverage_{category}"] = f"{coverage:.2f}"
        else:
            grouped[key]["booster_antigen"] = antigen
            grouped[key]["booster_coverage"] = f"{coverage:.2f}"

    output_rows: list[dict[str, str]] = []
    for key in sorted(grouped):
        row = grouped[key]
        dtp3 = (
            normalize_text(row.get("dtp3_coverage_official", ""))
            or normalize_text(row.get("dtp3_coverage_wuenic", ""))
            or normalize_text(row.get("dtp3_coverage_admin", ""))
        )
        output_rows.append(
            {
                "country_iso3": row["country_iso3"],
                "country_name": row["country_name"],
                "year": row["year"],
                "dtp1_coverage_admin": normalize_text(row.get("dtp1_coverage_admin", "")),
                "dtp1_coverage_official": normalize_text(row.get("dtp1_coverage_official", "")),
                "dtp3_coverage_admin": normalize_text(row.get("dtp3_coverage_admin", "")),
                "dtp3_coverage_official": normalize_text(row.get("dtp3_coverage_official", "")),
                "dtp3_coverage": dtp3,
                "booster_antigen": normalize_text(row.get("booster_antigen", "")),
                "booster_coverage": normalize_text(row.get("booster_coverage", "")),
                "source_name": row["source_name"],
                "source_url": row["source_url"],
                "source_release_date": row["source_release_date"],
                "data_freeze_date": row["data_freeze_date"],
                "notes": "dtp3_preference=official_then_wuenic_then_admin;booster_fields_blank_if_not_present",
            }
        )
    return output_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Clean WUENIC DTP coverage export.")
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--country-map",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_country_name_map.tsv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_wuenic_clean.tsv",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    if args.input is None:
        args.input = default_input_path()
    output_rows = build_output_rows(args.input, args.country_map)
    write_tsv(args.out, OUTPUT_COLUMNS, output_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
