#!/usr/bin/env python3
"""Clean WHO pertussis reported-case exports into a country-year table."""

from __future__ import annotations

import argparse
from pathlib import Path

from ph_utils import (
    current_freeze_date,
    find_first_input,
    load_country_name_map,
    normalize_country,
    normalize_text,
    parse_export_date_from_name,
    parse_int,
    read_xlsx_sheet_rows,
    repo_root,
    project_module_data_root,
    write_tsv,
)


OUTPUT_COLUMNS = [
    "country_iso3",
    "country_name",
    "year",
    "reported_cases",
    "incidence_per_100k",
    "source_name",
    "source_url",
    "source_release_date",
    "data_freeze_date",
    "raw_country_string",
    "notes",
]


PUBLIC_HEALTH_DATA_ROOT = project_module_data_root("public_health")


def default_input_path() -> Path:
    input_dir = repo_root() / "modules/public_health" / "inputs" / "raw" / "who_cases"
    preferred = sorted(
        path
        for path in input_dir.iterdir()
        if path.is_file()
        and path.suffix.lower() == ".xlsx"
        and "reported cases" in path.name.lower()
    )
    if preferred:
        return preferred[0]
    return find_first_input(input_dir, (".xlsx",))


def build_output_rows(input_path: Path, country_map_path: Path) -> list[dict[str, str]]:
    country_map = load_country_name_map(country_map_path)
    rows = read_xlsx_sheet_rows(input_path, sheet_name="Sheet1")
    header = rows[0]
    year_columns = [column for column in header[2:] if normalize_text(column).isdigit()]
    output_rows: list[dict[str, str]] = []

    for row in rows[1:]:
        padded = row + [""] * (len(header) - len(row))
        raw_country = normalize_text(padded[0])
        disease = normalize_text(padded[1])
        if not raw_country or disease != "Pertussis":
            continue

        country_row = normalize_country(raw_country, country_map)
        if country_row.get("match_status") != "normalized":
            continue

        for year in year_columns:
            value = normalize_text(padded[header.index(year)])
            reported_cases = parse_int(value)
            if reported_cases is None:
                continue
            output_rows.append(
                {
                    "country_iso3": country_row["country_iso3"],
                    "country_name": country_row["normalized_country_name"],
                    "year": year,
                    "reported_cases": str(reported_cases),
                    "incidence_per_100k": "",
                    "source_name": "WHO Immunization Data portal",
                    "source_url": "https://immunizationdata.who.int/dashboard",
                    "source_release_date": parse_export_date_from_name(input_path),
                    "data_freeze_date": current_freeze_date(),
                    "raw_country_string": raw_country,
                    "notes": "country_year_reported_cases_only;incidence_not_present_in_country_rows",
                }
            )

    output_rows.sort(key=lambda row: (row["country_iso3"], int(row["year"])))
    return output_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Clean WHO pertussis reported-case export.")
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
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_who_cases_clean.tsv",
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
