#!/usr/bin/env python3
"""Clean ECDC ESAC-Net J01FA and J01 exports into a harmonized Europe-focused table."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from ph_utils import (
    current_freeze_date,
    load_country_name_map,
    load_source_meta,
    normalize_country,
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

SOURCE_NAME = "ECDC ESAC-Net"
SOURCE_URL = "https://qap.ecdc.europa.eu/public/extensions/AMC2_Dashboard/AMC2_Dashboard.html"


PUBLIC_HEALTH_DATA_ROOT = project_module_data_root("public_health")


def discover_workbooks(input_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in input_dir.iterdir()
        if path.is_file() and path.suffix.lower() == ".xlsx" and not path.name.startswith("~$")
    )


def infer_care_setting(path: Path, source_meta_row: dict[str, str]) -> str:
    name = " ".join(
        [
            path.name.casefold(),
            normalize_text(source_meta_row.get("release_note", "")).casefold(),
            normalize_text(source_meta_row.get("notes", "")).casefold(),
        ]
    )
    if "total care" in name:
        return "total_care"
    if "community" in name:
        return "community"
    return "unknown"


def infer_metric_group(path: Path, source_meta_row: dict[str, str]) -> str:
    text = " ".join(
        [
            path.name.casefold(),
            normalize_text(source_meta_row.get("release_note", "")).casefold(),
            normalize_text(source_meta_row.get("notes", "")).casefold(),
        ]
    )
    if "j01fa" in text or "macrolide" in text:
        return "macrolide"
    if re.search(r"\bj01\b", text) or "antibacterials for systemic use" in text:
        return "total_antibiotic"
    return "unknown"


def format_metric(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6f}"


def build_output_rows(
    input_dir: Path,
    source_meta_path: Path,
    country_map_path: Path,
) -> list[dict[str, str]]:
    country_map = load_country_name_map(country_map_path)
    source_meta = load_source_meta(source_meta_path)

    metrics_by_key: dict[tuple[str, str], dict[str, dict[str, float]]] = {}
    country_name_by_key: dict[tuple[str, str], str] = {}
    source_dates_by_key: dict[tuple[str, str], list[str]] = {}
    source_names_by_key: dict[tuple[str, str], list[str]] = {}
    source_urls_by_key: dict[tuple[str, str], list[str]] = {}
    unresolved_countries: set[str] = set()

    for workbook_path in discover_workbooks(input_dir):
        meta_row = source_meta.get(workbook_path.name, {})
        care_setting = infer_care_setting(workbook_path, meta_row)
        metric_group = infer_metric_group(workbook_path, meta_row)
        release_date = (
            normalize_date_string(meta_row.get("export_date", ""))
            or parse_export_date_from_name(workbook_path)
        )
        source_name = normalize_text(meta_row.get("source_name", "")) or SOURCE_NAME
        source_url = normalize_text(meta_row.get("source_url", "")) or SOURCE_URL

        rows = rows_to_dicts(read_xlsx_sheet_rows(workbook_path, sheet_name="Sheet1"))
        for row in rows:
            raw_country = normalize_text(row.get("Country", ""))
            if not raw_country or "EU/EEA crude population-weighted mean" in raw_country:
                continue
            mapped = normalize_country(raw_country, country_map)
            if not normalize_text(mapped.get("country_iso3", "")):
                unresolved_countries.add(raw_country)
                continue
            year = normalize_text(row.get("Year", ""))
            value = parse_float(row.get("DDD per 1000 inhabitants per day", ""))
            if not year or value is None:
                continue
            key = (normalize_text(mapped["country_iso3"]), year)
            metrics_by_key.setdefault(key, {})
            metrics_by_key[key].setdefault(metric_group, {})
            metrics_by_key[key][metric_group][care_setting] = value
            country_name_by_key[key] = normalize_text(mapped.get("normalized_country_name", ""))
            source_dates_by_key.setdefault(key, [])
            source_names_by_key.setdefault(key, [])
            source_urls_by_key.setdefault(key, [])
            if release_date:
                source_dates_by_key[key].append(release_date)
            if source_name:
                source_names_by_key[key].append(source_name)
            if source_url:
                source_urls_by_key[key].append(source_url)

    if unresolved_countries:
        unresolved = ", ".join(sorted(unresolved_countries))
        raise ValueError(f"unresolved ESAC-Net countries: {unresolved}")

    output_rows: list[dict[str, str]] = []
    for key in sorted(metrics_by_key):
        country_iso3, year = key
        metrics = metrics_by_key[key]
        macrolide_metrics = metrics.get("macrolide", {})
        total_antibiotic_metrics = metrics.get("total_antibiotic", {})
        macrolide_total_care = macrolide_metrics.get("total_care")
        macrolide_community = macrolide_metrics.get("community")
        total_antibiotic_total_care = total_antibiotic_metrics.get("total_care")
        total_antibiotic_community = total_antibiotic_metrics.get("community")
        chosen_macrolide = (
            macrolide_total_care if macrolide_total_care is not None else macrolide_community
        )
        chosen_total_antibiotic = (
            total_antibiotic_total_care
            if total_antibiotic_total_care is not None
            else total_antibiotic_community
        )
        source_date = sorted(source_dates_by_key.get(key, []))[-1] if source_dates_by_key.get(key) else ""
        source_name = "; ".join(sorted(dict.fromkeys(source_names_by_key.get(key, []))))
        source_url = "; ".join(sorted(dict.fromkeys(source_urls_by_key.get(key, []))))

        notes = [
            "unit=DDD_per_1000_inhabitants_per_day",
        ]
        if macrolide_metrics:
            notes.append("ecdc_macrolide_atc_group=J01FA_macrolides")
            notes.append(
                "macrolide_preferred_sector="
                + ("total_care" if macrolide_total_care is not None else "community")
            )
            if macrolide_community is not None:
                notes.append(f"macrolide_community_sector_value={macrolide_community:.6f}")
            if macrolide_total_care is not None:
                notes.append(f"macrolide_total_care_sector_value={macrolide_total_care:.6f}")
        if total_antibiotic_metrics:
            notes.append("ecdc_total_antibiotic_atc_group=J01_systemic_antibacterials")
            notes.append(
                "total_antibiotic_preferred_sector="
                + ("total_care" if total_antibiotic_total_care is not None else "community")
            )
            if total_antibiotic_community is not None:
                notes.append(
                    f"total_antibiotic_community_sector_value={total_antibiotic_community:.6f}"
                )
            if total_antibiotic_total_care is not None:
                notes.append(
                    f"total_antibiotic_total_care_sector_value={total_antibiotic_total_care:.6f}"
                )
        # Unknown files are preserved in notes rather than silently discarded so
        # newer frozen exports remain auditable until the parser is expanded.
        if "unknown" in metrics:
            notes.append("ecdc_metric_group_unknown_for_one_or_more_files")

        output_rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": country_name_by_key.get(key, ""),
                "year": year,
                "macrolide_use_ddd_per_1000_per_day": format_metric(chosen_macrolide),
                "total_antibiotic_use_ddd_per_1000_per_day": format_metric(chosen_total_antibiotic),
                "source_name": source_name,
                "source_url": source_url,
                "source_release_date": source_date,
                "data_freeze_date": current_freeze_date(),
                "notes": ";".join(note for note in notes if note),
            }
        )
    return output_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Clean ECDC ESAC-Net AMU export.")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=repo_root() / "modules/public_health" / "inputs" / "raw" / "esacnet_amu",
    )
    parser.add_argument(
        "--source-meta",
        type=Path,
        default=repo_root() / "modules/public_health" / "inputs" / "raw" / "esacnet_amu" / "source_meta.tsv",
    )
    parser.add_argument(
        "--country-map",
        type=Path,
        default=repo_root() / "outputs" / "workflow" / "epi" / "ph_country_name_map.tsv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_esacnet_amu_clean.tsv",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    output_rows = build_output_rows(args.input_dir, args.source_meta, args.country_map)
    write_tsv(args.out, OUTPUT_COLUMNS, output_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
