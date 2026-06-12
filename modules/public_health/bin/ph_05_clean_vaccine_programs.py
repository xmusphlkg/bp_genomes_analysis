#!/usr/bin/env python3
"""Build interval-based vaccine-program metadata from WHO and curated inputs."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from ph_utils import (
    current_freeze_date,
    find_first_input,
    normalize_text,
    parse_export_date_from_name,
    parse_float,
    parse_int,
    read_xlsx_sheet_rows,
    read_vaccine_program_rows,
    repo_root,
    project_module_data_root,
    rows_to_dicts,
    write_tsv,
)


OUTPUT_COLUMNS = [
    "country_iso3",
    "country_name",
    "year_start",
    "year_end",
    "vaccine_program_type",
    "acellular_vs_whole_cell",
    "prn_in_vaccine",
    "primary_series_schedule",
    "booster_schedule",
    "program_change_year",
    "source_name",
    "source_url",
    "source_release_date",
    "curation_confidence",
    "notes",
]


PUBLIC_HEALTH_DATA_ROOT = project_module_data_root("public_health")


def introduction_status_category(value: str) -> str:
    value = normalize_text(value)
    if value == "Yes":
        return "routine"
    if value in {"Yes (P)", "Yes (R)", "Yes (A)", "Yes (O)", "High risk"}:
        return "targeted"
    return "none"


def load_ap_introduction_rows(path: Path) -> dict[str, dict[str, str]]:
    raw_rows = rows_to_dicts(read_xlsx_sheet_rows(path, sheet_name="Sheet1"))
    by_country: dict[str, dict[str, str]] = {}
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in raw_rows:
        iso3 = normalize_text(row.get("ISO_3_CODE", ""))
        year = normalize_text(row.get("YEAR", ""))
        if iso3 and year.isdigit():
            grouped[iso3].append(row)

    for iso3, rows in grouped.items():
        rows.sort(key=lambda row: int(normalize_text(row.get("YEAR", "0")) or "0"))
        any_ap_rows = [row for row in rows if introduction_status_category(row.get("INTRO", "")) in {"routine", "targeted"}]
        routine_rows = [row for row in rows if introduction_status_category(row.get("INTRO", "")) == "routine"]
        by_country[iso3] = {
            "country_name": normalize_text(rows[0].get("COUNTRYNAME", "")),
            "min_year": normalize_text(rows[0].get("YEAR", "")),
            "max_year": normalize_text(rows[-1].get("YEAR", "")),
            "first_any_ap_year": "" if not any_ap_rows else normalize_text(any_ap_rows[0].get("YEAR", "")),
            "first_any_ap_status": "" if not any_ap_rows else normalize_text(any_ap_rows[0].get("INTRO", "")),
            "first_routine_ap_year": "" if not routine_rows else normalize_text(routine_rows[0].get("YEAR", "")),
            "first_routine_ap_status": "" if not routine_rows else normalize_text(routine_rows[0].get("INTRO", "")),
            "latest_ap_status": next(
                (
                    normalize_text(row.get("INTRO", ""))
                    for row in reversed(rows)
                    if normalize_text(row.get("INTRO", ""))
                ),
                "",
            ),
        }
    return by_country


def load_program_csv(path: Path) -> dict[str, dict[str, str]]:
    rows = read_vaccine_program_rows(path)
    return {normalize_text(row.get("CODE", "")): row for row in rows if normalize_text(row.get("CODE", ""))}


def schedule_summary(program_row: dict[str, str]) -> tuple[str, str]:
    first_month = parse_float(program_row.get("TimeFirstShot", ""))
    last_month = parse_float(program_row.get("TimeLastShot", ""))
    dose_count = parse_int(program_row.get("VaccineDose", ""))
    primary = ""
    if first_month is not None or dose_count is not None:
        primary = f"first_dose_month={'' if first_month is None else f'{first_month:.2f}'};routine_dose_count={'' if dose_count is None else dose_count}"

    booster_parts = []
    if last_month is not None:
        booster_parts.append(f"last_routine_dose_month={last_month:.2f}")
    booster_parts.append(f"pregnant={'yes' if normalize_text(program_row.get('VaccinePregnant', '')) == '1' else 'no'}")
    booster_parts.append(f"adult={'yes' if normalize_text(program_row.get('VaccineAdult', '')) == '1' else 'no'}")
    booster_parts.append(f"risk={'yes' if normalize_text(program_row.get('VaccineRisk', '')) == '1' else 'no'}")
    if normalize_text(program_row.get("VaccinePregnantTime", "")):
        booster_parts.append(f"pregnant_timing={normalize_text(program_row['VaccinePregnantTime'])}")
    if normalize_text(program_row.get("VaccinePregnantIntroYear", "")):
        booster_parts.append(f"pregnant_intro_year={normalize_text(program_row['VaccinePregnantIntroYear'])}")
    if normalize_text(program_row.get("VaccinePregnantIntroDate", "")):
        booster_parts.append(f"pregnant_intro_date={normalize_text(program_row['VaccinePregnantIntroDate'])}")
    return primary, ";".join(booster_parts)


def append_interval(
    output_rows: list[dict[str, str]],
    *,
    iso3: str,
    country_name: str,
    year_start: int,
    year_end: int,
    vaccine_program_type: str,
    acellular_vs_whole_cell: str,
    program_change_year: str,
    primary_series_schedule: str,
    booster_schedule: str,
    source_release_date: str,
    notes: str,
) -> None:
    output_rows.append(
        {
            "country_iso3": iso3,
            "country_name": country_name,
            "year_start": str(year_start),
            "year_end": str(year_end),
            "vaccine_program_type": vaccine_program_type,
            "acellular_vs_whole_cell": acellular_vs_whole_cell,
            "prn_in_vaccine": "unknown",
            "primary_series_schedule": primary_series_schedule,
            "booster_schedule": booster_schedule,
            "program_change_year": program_change_year,
            "source_name": "WHO aP introduction export; supplemental vaccine_program.csv",
            "source_url": "https://immunizationdata.who.int/dashboard",
            "source_release_date": source_release_date,
            "curation_confidence": "medium",
            "notes": notes,
        }
    )


def build_output_rows(ap_path: Path, csv_path: Path) -> list[dict[str, str]]:
    ap_by_country = load_ap_introduction_rows(ap_path)
    program_by_country = load_program_csv(csv_path)
    source_release_date = parse_export_date_from_name(ap_path) or current_freeze_date()

    output_rows: list[dict[str, str]] = []
    for iso3 in sorted(ap_by_country):
        ap_row = ap_by_country[iso3]
        country_name = ap_row["country_name"]
        min_year = int(ap_row["min_year"])
        max_year = int(ap_row["max_year"])
        first_any_ap_year = parse_int(ap_row.get("first_any_ap_year", ""))
        first_routine_ap_year = parse_int(ap_row.get("first_routine_ap_year", ""))
        program_row = program_by_country.get(iso3, {})
        primary_series_schedule, booster_schedule = schedule_summary(program_row)
        supplemental_note = ""
        if program_row:
            supplemental_note = "supplemental_program_csv_present"
        else:
            supplemental_note = "supplemental_program_csv_missing"

        common_status_note = (
            f"first_any_ap_year={ap_row.get('first_any_ap_year', '') or 'na'};"
            f"first_any_ap_status={ap_row.get('first_any_ap_status', '') or 'na'};"
            f"first_routine_ap_year={ap_row.get('first_routine_ap_year', '') or 'na'};"
            f"first_routine_ap_status={ap_row.get('first_routine_ap_status', '') or 'na'};"
            f"latest_ap_status={ap_row.get('latest_ap_status', '') or 'na'}"
        )

        if first_any_ap_year is None:
            append_interval(
                output_rows,
                iso3=iso3,
                country_name=country_name,
                year_start=min_year,
                year_end=max_year,
                vaccine_program_type="whole_cell_or_unknown",
                acellular_vs_whole_cell="whole_cell_or_unknown",
                program_change_year="",
                primary_series_schedule=primary_series_schedule,
                booster_schedule=booster_schedule,
                source_release_date=source_release_date,
                notes=f"no_ap_intro_in_who_export;{common_status_note};{supplemental_note}",
            )
            continue

        if first_any_ap_year > min_year:
            append_interval(
                output_rows,
                iso3=iso3,
                country_name=country_name,
                year_start=min_year,
                year_end=first_any_ap_year - 1,
                vaccine_program_type="pre_ap_introduction_or_whole_cell",
                acellular_vs_whole_cell="whole_cell_or_unknown",
                program_change_year=str(first_any_ap_year),
                primary_series_schedule=primary_series_schedule,
                booster_schedule=booster_schedule,
                source_release_date=source_release_date,
                notes=f"before_first_ap_intro;{common_status_note};{supplemental_note}",
            )

        if first_routine_ap_year is None:
            append_interval(
                output_rows,
                iso3=iso3,
                country_name=country_name,
                year_start=first_any_ap_year,
                year_end=max_year,
                vaccine_program_type="ap_targeted_or_partial_nonroutine",
                acellular_vs_whole_cell="targeted_or_partial_ap",
                program_change_year=str(first_any_ap_year),
                primary_series_schedule=primary_series_schedule,
                booster_schedule=booster_schedule,
                source_release_date=source_release_date,
                notes=f"nonroutine_or_partial_ap_only;{common_status_note};{supplemental_note}",
            )
            continue

        if first_any_ap_year < first_routine_ap_year:
            append_interval(
                output_rows,
                iso3=iso3,
                country_name=country_name,
                year_start=first_any_ap_year,
                year_end=first_routine_ap_year - 1,
                vaccine_program_type="ap_targeted_or_partial_nonroutine",
                acellular_vs_whole_cell="targeted_or_partial_ap",
                program_change_year=str(first_any_ap_year),
                primary_series_schedule=primary_series_schedule,
                booster_schedule=booster_schedule,
                source_release_date=source_release_date,
                notes=f"before_first_routine_ap_intro;{common_status_note};{supplemental_note}",
            )

        append_interval(
            output_rows,
            iso3=iso3,
            country_name=country_name,
            year_start=first_routine_ap_year,
            year_end=max_year,
            vaccine_program_type="ap_introduced_routine_or_mixed",
            acellular_vs_whole_cell="mixed_or_acellular",
            program_change_year=str(first_routine_ap_year),
            primary_series_schedule=primary_series_schedule,
            booster_schedule=booster_schedule,
            source_release_date=source_release_date,
            notes=f"routine_ap_intro_recorded;{common_status_note};{supplemental_note}",
        )

    return output_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build vaccine-program metadata table.")
    parser.add_argument(
        "--ap-input",
        type=Path,
        default=repo_root() / "modules/public_health" / "inputs" / "raw" / "vaccine_program_docs" / "Introduction of aP (acellular pertussis) vaccine 2026-21-03 23-11 UTC.xlsx",
    )
    parser.add_argument(
        "--supplemental-csv",
        type=Path,
        default=repo_root() / "modules/public_health" / "inputs" / "raw" / "vaccine_program_docs" / "vaccine_program.csv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_country_program_metadata.tsv",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    output_rows = build_output_rows(args.ap_input, args.supplemental_csv)
    write_tsv(args.out, OUTPUT_COLUMNS, output_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
