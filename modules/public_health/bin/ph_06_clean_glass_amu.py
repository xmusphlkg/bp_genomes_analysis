#!/usr/bin/env python3
"""Clean WHO GLASS antimicrobial-use exports into a country-year AMU table."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

from ph_utils import (
    current_freeze_date,
    find_first_input,
    load_source_meta,
    normalize_date_string,
    normalize_text,
    parse_float,
    parse_export_date_from_name,
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
    "macrolide_use_ddd_per_1000_per_day",
    "total_antibiotic_use_ddd_per_1000_per_day",
    "source_name",
    "source_url",
    "source_release_date",
    "data_freeze_date",
    "notes",
]

SOURCE_URL = "https://www.who.int/news/item/25-09-2025-updated-who-dashboard-offers-new-insights-on-antimicrobial-resistance-and-use"
SOURCE_NAME = "WHO GLASS AMU"


PUBLIC_HEALTH_DATA_ROOT = project_module_data_root("public_health")


def find_glass_workbook(input_dir: Path) -> Path:
    candidates = sorted(
        path
        for path in input_dir.rglob("*")
        if path.is_file()
        and path.suffix.lower() == ".xlsx"
        and not path.name.startswith("~$")
        and "glass" in path.name.casefold()
    )
    if candidates:
        return candidates[0]
    return find_first_input(input_dir, (".xlsx",))


def parse_release_date_from_intro(workbook_path: Path) -> str:
    intro_rows = read_xlsx_sheet_rows(workbook_path, sheet_name="Introduction")
    for row in intro_rows:
        text = normalize_text(" ".join(row))
        match = re.search(r"Data released on:\s*([0-9]{2}-[0-9]{2}-[0-9]{4})", text)
        if match:
            return normalize_date_string(match.group(1))
    return ""


def load_note_map(workbook_path: Path) -> dict[str, str]:
    rows = rows_to_dicts(read_xlsx_sheet_rows(workbook_path, sheet_name="Notes"))
    note_map: dict[str, str] = {}
    for row in rows:
        code = normalize_text(row.get("Notes", "")).casefold()
        explanation = normalize_text(row.get("Explenation", "")) or normalize_text(row.get("Explanation", ""))
        if code and explanation:
            note_map[code] = explanation
    return note_map


def parse_note_codes(raw_value: str, note_map: dict[str, str]) -> list[str]:
    text = normalize_text(raw_value)
    if not text:
        return []
    codes: list[str] = []
    for token in re.findall(r"[A-Za-z]+", text):
        token = token.casefold()
        if token in note_map:
            codes.append(token)
            continue
        if len(token) > 1:
            for char in token:
                if char in note_map:
                    codes.append(char)
    return sorted(dict.fromkeys(codes))


def format_metric(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6f}"


def build_output_rows(
    input_dir: Path,
    source_meta_path: Path,
) -> list[dict[str, str]]:
    workbook_path = find_glass_workbook(input_dir)
    source_meta = load_source_meta(source_meta_path).get(workbook_path.name, {})
    note_map = load_note_map(workbook_path)
    release_date = (
        parse_release_date_from_intro(workbook_path)
        or normalize_date_string(source_meta.get("export_date", ""))
        or parse_export_date_from_name(workbook_path)
    )
    source_name = normalize_text(source_meta.get("source_name", "")) or SOURCE_NAME
    source_url = normalize_text(source_meta.get("source_url", "")) or SOURCE_URL

    context_rows = rows_to_dicts(read_xlsx_sheet_rows(workbook_path, sheet_name="Data contextual info"))
    aware_rows = rows_to_dicts(read_xlsx_sheet_rows(workbook_path, sheet_name="Antibiotic_Use_AWaRe"))
    atc_rows = rows_to_dicts(read_xlsx_sheet_rows(workbook_path, sheet_name="Antimicrobial_Use_ATC4"))

    context_by_key: dict[tuple[str, str], dict[str, str]] = {}
    for row in context_rows:
        key = (normalize_text(row.get("CountryIso3", "")), normalize_text(row.get("Year", "")))
        if all(key) and key not in context_by_key:
            context_by_key[key] = row

    totals_by_key: dict[tuple[str, str], float] = defaultdict(float)
    total_notes_by_key: dict[tuple[str, str], list[str]] = defaultdict(list)
    total_country_names: dict[tuple[str, str], str] = {}
    for row in aware_rows:
        key = (normalize_text(row.get("COUNTRY", "")), normalize_text(row.get("YEAR", "")))
        did = parse_float(row.get("DID", ""))
        if not all(key) or did is None:
            continue
        totals_by_key[key] += did
        total_country_names.setdefault(key, normalize_text(row.get("CountryTerritoryArea", "")))
        total_notes_by_key[key].extend(parse_note_codes(row.get("Notes", ""), note_map))

    macrolide_by_key: dict[tuple[str, str], float] = defaultdict(float)
    macrolide_notes_by_key: dict[tuple[str, str], list[str]] = defaultdict(list)
    macrolide_country_names: dict[tuple[str, str], str] = {}
    for row in atc_rows:
        if normalize_text(row.get("ATC4", "")) != "J01FA":
            continue
        key = (normalize_text(row.get("CountryIso3", "")), normalize_text(row.get("Year", "")))
        did = parse_float(row.get("DID", ""))
        if not all(key) or did is None:
            continue
        macrolide_by_key[key] += did
        macrolide_country_names.setdefault(key, normalize_text(row.get("CountryTerritoryArea", "")))
        macrolide_notes_by_key[key].extend(parse_note_codes(row.get("Notes", ""), note_map))

    output_rows: list[dict[str, str]] = []
    for key in sorted(set(context_by_key) | set(totals_by_key) | set(macrolide_by_key)):
        country_iso3, year = key
        context = context_by_key.get(key, {})
        note_codes = sorted(
            dict.fromkeys(total_notes_by_key.get(key, []) + macrolide_notes_by_key.get(key, []))
        )
        note_summaries = [note_map[code] for code in note_codes if code in note_map]
        notes = [
            "unit=DDD_per_1000_inhabitants_per_day",
            "total_antibiotic_metric=aware_sum_did",
            "macrolide_metric=atc4_j01fa_sum_did",
        ]
        population_coverage = normalize_text(context.get("Population coverage (%)", ""))
        health_level = normalize_text(context.get("HealthLevel", ""))
        health_sector = normalize_text(context.get("HealthSector", ""))
        type_of_data = normalize_text(context.get("TypeOfData", ""))
        comparability = normalize_text(context.get("ComparabilityWithPreviousYearData", ""))
        disclaimer = normalize_text(context.get("Disclaimer", ""))
        if population_coverage:
            notes.append(f"population_coverage_pct={population_coverage}")
        if health_level:
            notes.append(f"glass_health_level={health_level}")
        if health_sector:
            notes.append(f"glass_health_sector={health_sector}")
        if type_of_data:
            notes.append(f"glass_type_of_data={type_of_data}")
        if comparability:
            notes.append(f"comparability_with_previous_year={comparability}")
        if note_codes:
            notes.append(f"glass_note_codes={','.join(note_codes)}")
        if note_summaries:
            notes.append("glass_note_summary=" + " | ".join(note_summaries))
        if disclaimer:
            notes.append(f"glass_disclaimer={disclaimer}")

        country_name = (
            normalize_text(context.get("CountryTerritoryArea", ""))
            or total_country_names.get(key, "")
            or macrolide_country_names.get(key, "")
        )
        output_rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": country_name,
                "year": year,
                "macrolide_use_ddd_per_1000_per_day": format_metric(macrolide_by_key.get(key)),
                "total_antibiotic_use_ddd_per_1000_per_day": format_metric(totals_by_key.get(key)),
                "source_name": source_name,
                "source_url": source_url,
                "source_release_date": release_date,
                "data_freeze_date": current_freeze_date(),
                "notes": ";".join(note for note in notes if note),
            }
        )
    return output_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Clean WHO GLASS AMU export.")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=repo_root() / "modules/public_health" / "inputs" / "raw" / "glass_amu",
    )
    parser.add_argument(
        "--source-meta",
        type=Path,
        default=repo_root() / "modules/public_health" / "inputs" / "raw" / "glass_amu" / "source_meta.tsv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_glass_amu_clean.tsv",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    output_rows = build_output_rows(args.input_dir, args.source_meta)
    write_tsv(args.out, OUTPUT_COLUMNS, output_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
