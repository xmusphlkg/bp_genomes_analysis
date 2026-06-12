#!/usr/bin/env python3
"""Build the master public-health country-year table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from ph_utils import (
    current_freeze_date,
    index_reporting_era_rows,
    match_reporting_era_row,
    normalize_text,
    parse_int,
    reporting_era_post_flag,
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
    "region_who",
    "reported_cases",
    "incidence_per_100k",
    "dtp3_coverage",
    "booster_coverage",
    "vaccine_program_type",
    "prn_in_vaccine",
    "acellular_vs_whole_cell",
    "macrolide_use_ddd_per_1000_per_day",
    "total_antibiotic_use_ddd_per_1000_per_day",
    "post_covid_period",
    "reporting_era_record_iso3",
    "reporting_era_scope_type",
    "reporting_era_match_type",
    "reporting_era_confidence",
    "pcr_lab_guideline_year",
    "reporting_case_definition_change_year",
    "surveillance_platform_change_year",
    "post_pcr_lab_guideline_era",
    "post_reporting_case_definition_change_era",
    "post_surveillance_platform_change_era",
    "genomes_count",
    "n_genomes_prn_interpretable",
    "n_prn_disrupted",
    "frac_prn_disrupted",
    "n_read_supported_prn_disrupted",
    "n_mr_marked",
    "frac_23s_A2047G",
    "genomes_per_case",
    "surveillance_source",
    "vaccine_source",
    "amu_source",
    "source_release_date",
    "cases_last_updated",
    "vaccine_last_updated",
    "amu_last_updated",
    "data_freeze_date",
    "notes",
]


PUBLIC_HEALTH_DATA_ROOT = project_module_data_root("public_health")


def load_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_who_region_map(path: Path) -> dict[str, str]:
    rows = rows_to_dicts(read_xlsx_sheet_rows(path, sheet_name="Sheet1"))
    region_map: dict[str, str] = {}
    for row in rows:
        iso3 = normalize_text(row.get("ISO_3_CODE", ""))
        region = normalize_text(row.get("WHO_REGION", ""))
        if iso3 and region and iso3 not in region_map:
            region_map[iso3] = region
    return region_map


def latest_date(values: list[str]) -> str:
    dates = sorted(value for value in values if normalize_text(value))
    return dates[-1] if dates else ""


def program_row_for_year(program_rows: list[dict[str, str]], year: int) -> dict[str, str] | None:
    for row in program_rows:
        start = parse_int(row.get("year_start", ""))
        end = parse_int(row.get("year_end", ""))
        if start is None or end is None:
            continue
        if start <= year <= end:
            return row
    return None


def build_output_rows(
    cases_rows: list[dict[str, str]],
    vaccine_rows: list[dict[str, str]],
    program_rows: list[dict[str, str]],
    glass_rows: list[dict[str, str]],
    esac_rows: list[dict[str, str]],
    reporting_era_rows: list[dict[str, str]],
    region_map: dict[str, str],
) -> list[dict[str, str]]:
    cases_by_key = {
        (normalize_text(row["country_iso3"]), normalize_text(row["year"])): row
        for row in cases_rows
    }
    vaccine_by_key = {
        (normalize_text(row["country_iso3"]), normalize_text(row["year"])): row
        for row in vaccine_rows
    }
    glass_by_key = {
        (normalize_text(row["country_iso3"]), normalize_text(row["year"])): row
        for row in glass_rows
    }
    esac_by_key = {
        (normalize_text(row["country_iso3"]), normalize_text(row["year"])): row
        for row in esac_rows
    }
    program_by_country: dict[str, list[dict[str, str]]] = {}
    for row in program_rows:
        program_by_country.setdefault(normalize_text(row["country_iso3"]), []).append(row)
    reporting_era_index = index_reporting_era_rows(reporting_era_rows)

    all_keys = sorted(set(cases_by_key) | set(vaccine_by_key) | set(glass_by_key) | set(esac_by_key))
    output_rows: list[dict[str, str]] = []
    for country_iso3, year_text in all_keys:
        year = int(year_text)
        case_row = cases_by_key.get((country_iso3, year_text), {})
        vaccine_row = vaccine_by_key.get((country_iso3, year_text), {})
        glass_row = glass_by_key.get((country_iso3, year_text), {})
        esac_row = esac_by_key.get((country_iso3, year_text), {})
        program_row = program_row_for_year(program_by_country.get(country_iso3, []), year) or {}
        region_who = region_map.get(country_iso3, "")
        reporting_era_row, reporting_era_match_type = match_reporting_era_row(
            country_iso3,
            region_who,
            reporting_era_index,
        )

        country_name = (
            normalize_text(case_row.get("country_name", ""))
            or normalize_text(vaccine_row.get("country_name", ""))
            or normalize_text(program_row.get("country_name", ""))
            or normalize_text(esac_row.get("country_name", ""))
            or normalize_text(glass_row.get("country_name", ""))
        )
        macrolide_value = normalize_text(esac_row.get("macrolide_use_ddd_per_1000_per_day", "")) or normalize_text(
            glass_row.get("macrolide_use_ddd_per_1000_per_day", "")
        )
        total_antibiotic_value = normalize_text(
            glass_row.get("total_antibiotic_use_ddd_per_1000_per_day", "")
        ) or normalize_text(esac_row.get("total_antibiotic_use_ddd_per_1000_per_day", ""))
        amu_source_parts = [
            normalize_text(esac_row.get("source_name", "")),
            normalize_text(glass_row.get("source_name", "")),
        ]
        amu_source = "; ".join(part for part in dict.fromkeys(amu_source_parts) if part)

        notes = [
            "post_covid_period_defined_as_year_ge_2024_for_country_year_resolution",
        ]
        if not macrolide_value and not total_antibiotic_value:
            notes.append("amu_missing_country_level_data")
        elif not macrolide_value:
            notes.append("macrolide_use_missing_country_year")
        elif not total_antibiotic_value:
            notes.append("total_antibiotic_use_missing_country_year")
        if normalize_text(esac_row.get("macrolide_use_ddd_per_1000_per_day", "")):
            notes.append("macrolide_metric_preferred_from_esacnet")
        elif normalize_text(glass_row.get("macrolide_use_ddd_per_1000_per_day", "")):
            notes.append("macrolide_metric_fallback_from_glass")
        if normalize_text(glass_row.get("total_antibiotic_use_ddd_per_1000_per_day", "")):
            notes.append("total_antibiotic_metric_preferred_from_glass")
        elif normalize_text(esac_row.get("total_antibiotic_use_ddd_per_1000_per_day", "")):
            notes.append("total_antibiotic_metric_fallback_from_esacnet")
        if not normalize_text(vaccine_row.get("booster_coverage", "")):
            notes.append("booster_coverage_not_present_in_current_wuenic_export")

        source_release_dates = [
            normalize_text(case_row.get("source_release_date", "")),
            normalize_text(vaccine_row.get("source_release_date", "")),
            normalize_text(program_row.get("source_release_date", "")),
            normalize_text(glass_row.get("source_release_date", "")),
            normalize_text(esac_row.get("source_release_date", "")),
        ]
        vaccine_source_parts = []
        if vaccine_row:
            vaccine_source_parts.append(normalize_text(vaccine_row.get("source_name", "")))
        if program_row:
            vaccine_source_parts.append(normalize_text(program_row.get("source_name", "")))

        output_rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": country_name,
                "year": year_text,
                "region_who": region_who,
                "reported_cases": normalize_text(case_row.get("reported_cases", "")),
                "incidence_per_100k": normalize_text(case_row.get("incidence_per_100k", "")),
                "dtp3_coverage": normalize_text(vaccine_row.get("dtp3_coverage", "")),
                "booster_coverage": normalize_text(vaccine_row.get("booster_coverage", "")),
                "vaccine_program_type": normalize_text(program_row.get("vaccine_program_type", "")),
                "prn_in_vaccine": normalize_text(program_row.get("prn_in_vaccine", "")),
                "acellular_vs_whole_cell": normalize_text(program_row.get("acellular_vs_whole_cell", "")),
                "macrolide_use_ddd_per_1000_per_day": macrolide_value,
                "total_antibiotic_use_ddd_per_1000_per_day": total_antibiotic_value,
                "post_covid_period": "1" if year >= 2024 else "0",
                "reporting_era_record_iso3": normalize_text(reporting_era_row.get("iso3", "")),
                "reporting_era_scope_type": normalize_text(reporting_era_row.get("scope_type", "")),
                "reporting_era_match_type": reporting_era_match_type,
                "reporting_era_confidence": normalize_text(reporting_era_row.get("confidence", "")),
                "pcr_lab_guideline_year": normalize_text(reporting_era_row.get("pcr_lab_guideline_year_min", "")),
                "reporting_case_definition_change_year": normalize_text(
                    reporting_era_row.get("reporting_case_definition_change_year_min", "")
                ),
                "surveillance_platform_change_year": normalize_text(
                    reporting_era_row.get("surveillance_platform_change_year_min", "")
                ),
                "post_pcr_lab_guideline_era": reporting_era_post_flag(
                    year,
                    normalize_text(reporting_era_row.get("pcr_lab_guideline_year_min", "")),
                ),
                "post_reporting_case_definition_change_era": reporting_era_post_flag(
                    year,
                    normalize_text(reporting_era_row.get("reporting_case_definition_change_year_min", "")),
                ),
                "post_surveillance_platform_change_era": reporting_era_post_flag(
                    year,
                    normalize_text(reporting_era_row.get("surveillance_platform_change_year_min", "")),
                ),
                "genomes_count": "",
                "n_genomes_prn_interpretable": "",
                "n_prn_disrupted": "",
                "frac_prn_disrupted": "",
                "n_read_supported_prn_disrupted": "",
                "n_mr_marked": "",
                "frac_23s_A2047G": "",
                "genomes_per_case": "",
                "surveillance_source": normalize_text(case_row.get("source_name", "")),
                "vaccine_source": "; ".join(part for part in vaccine_source_parts if part),
                "amu_source": amu_source,
                "source_release_date": latest_date(source_release_dates),
                "cases_last_updated": normalize_text(case_row.get("source_release_date", "")),
                "vaccine_last_updated": latest_date(
                    [
                        normalize_text(vaccine_row.get("source_release_date", "")),
                        normalize_text(program_row.get("source_release_date", "")),
                    ]
                ),
                "amu_last_updated": latest_date(
                    [
                        normalize_text(glass_row.get("source_release_date", "")),
                        normalize_text(esac_row.get("source_release_date", "")),
                    ]
                ),
                "data_freeze_date": current_freeze_date(),
                "notes": ";".join(notes),
            }
        )

    return output_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build the master public-health country-year table.")
    parser.add_argument(
        "--cases",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_who_cases_clean.tsv",
    )
    parser.add_argument(
        "--vaccine",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_wuenic_clean.tsv",
    )
    parser.add_argument(
        "--program",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_country_program_metadata.tsv",
    )
    parser.add_argument(
        "--glass",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_glass_amu_clean.tsv",
    )
    parser.add_argument(
        "--esacnet",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_esacnet_amu_clean.tsv",
    )
    parser.add_argument(
        "--who-region-source",
        type=Path,
        default=repo_root() / "modules" / "public_health" / "inputs" / "raw" / "vaccine_program_docs" / "Introduction of aP (acellular pertussis) vaccine 2026-21-03 23-11 UTC.xlsx",
    )
    parser.add_argument(
        "--reporting-era",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_reporting_era_indicators.tsv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_country_year_master.tsv",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    output_rows = build_output_rows(
        cases_rows=load_tsv_rows(args.cases),
        vaccine_rows=load_tsv_rows(args.vaccine),
        program_rows=load_tsv_rows(args.program),
        glass_rows=load_tsv_rows(args.glass),
        esac_rows=load_tsv_rows(args.esacnet),
        reporting_era_rows=load_tsv_rows(args.reporting_era) if args.reporting_era.exists() else [],
        region_map=load_who_region_map(args.who_region_source),
    )
    write_tsv(args.out, OUTPUT_COLUMNS, output_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
