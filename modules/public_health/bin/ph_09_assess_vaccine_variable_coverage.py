#!/usr/bin/env python3
"""Assess vaccine-variable coverage for downstream ecology modeling."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from ph_utils import project_workflow_root


REPO_ROOT = Path(__file__).resolve().parents[3]
WORKFLOW_DATA_ROOT = project_workflow_root()
UNKNOWN_TOKENS = {"", "unknown", "na", "n/a", "none", "missing", "nan"}
LEGACY_CURATION_REQUIRED_COLUMNS = {"country_iso3", "prn_in_vaccine"}
RICH_CURATION_REQUIRED_COLUMNS = {"country_iso3", "prn_in_vaccine_curated", "year_start", "year_end"}


def normalize_text(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def is_known_series(series: pd.Series) -> pd.Series:
    normalized = series.fillna("").astype(str).str.strip().str.casefold()
    return ~normalized.isin(UNKNOWN_TOKENS)


def parse_year(value: object) -> int | None:
    text = normalize_text(value)
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def load_prn_curation(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["country_iso3", "prn_in_vaccine", "year_start", "year_end", "notes"])

    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    column_names = set(table.columns)

    if LEGACY_CURATION_REQUIRED_COLUMNS.issubset(column_names):
        normalized = table.copy()
    elif RICH_CURATION_REQUIRED_COLUMNS.issubset(column_names):
        normalized = pd.DataFrame(
            {
                "country_iso3": table["country_iso3"],
                "prn_in_vaccine": table["prn_in_vaccine_curated"],
                "year_start": table["year_start"],
                "year_end": table["year_end"],
                "notes": table.get("notes", ""),
            }
        )
    else:
        expected = ", ".join(sorted(LEGACY_CURATION_REQUIRED_COLUMNS | RICH_CURATION_REQUIRED_COLUMNS))
        raise ValueError(
            "Prn curation file is missing required columns for either the legacy overlay "
            f"or the richer formulation curation format. Expected a subset of: {expected}"
        )

    for optional_col in ("year_start", "year_end", "notes"):
        if optional_col not in normalized.columns:
            normalized[optional_col] = ""
    return normalized


def build_curation_mask(ph: pd.DataFrame, curation: pd.DataFrame) -> tuple[pd.Series, int]:
    if curation.empty or "country_iso3" not in ph.columns:
        return pd.Series(False, index=ph.index), 0

    iso_series = ph["country_iso3"].fillna("").astype(str).str.strip().str.upper()
    has_year = "year" in ph.columns
    year_series = pd.to_numeric(ph["year"], errors="coerce") if has_year else None

    curation_mask = pd.Series(False, index=ph.index)
    applied_rows = 0

    for row in curation.to_dict(orient="records"):
        iso3 = normalize_text(row.get("country_iso3", "")).upper()
        prn_value = normalize_text(row.get("prn_in_vaccine", "")).casefold()
        if not iso3 or prn_value in UNKNOWN_TOKENS:
            continue

        row_mask = iso_series.eq(iso3)
        if has_year:
            start_year = parse_year(row.get("year_start", ""))
            end_year = parse_year(row.get("year_end", ""))
            if start_year is not None:
                row_mask &= year_series.ge(start_year).fillna(False)
            if end_year is not None:
                row_mask &= year_series.le(end_year).fillna(False)

        if row_mask.any():
            curation_mask |= row_mask
            applied_rows += 1

    return curation_mask, applied_rows


def assess_vaccine_variable_coverage(
    ph_master_path: str,
    min_countries_ap: int = 15,
    min_countries_prn_form: int = 5,
    prn_curation_path: str = "",
) -> dict:
    """Assess whether the public panel is ready for the ecology exposure index."""

    ph = pd.read_csv(ph_master_path, sep="\t", dtype=str)

    has_program_type = (
        is_known_series(ph["vaccine_program_type"])
        if "vaccine_program_type" in ph.columns
        else pd.Series(False, index=ph.index)
    )
    has_acellular = (
        is_known_series(ph["acellular_vs_whole_cell"])
        if "acellular_vs_whole_cell" in ph.columns
        else pd.Series(False, index=ph.index)
    )
    has_ap_info = has_program_type | has_acellular
    countries_with_ap = ph.loc[has_ap_info, "country_iso3"].nunique() if "country_iso3" in ph.columns else 0

    has_prn_form_master = (
        is_known_series(ph["prn_in_vaccine"])
        if "prn_in_vaccine" in ph.columns
        else pd.Series(False, index=ph.index)
    )
    countries_with_prn_form_master = (
        ph.loc[has_prn_form_master, "country_iso3"].nunique() if "country_iso3" in ph.columns else 0
    )

    curation_file = Path(prn_curation_path) if prn_curation_path else Path()
    curation_rows = load_prn_curation(curation_file) if prn_curation_path else pd.DataFrame()
    has_prn_form_curation, curation_rows_applied = build_curation_mask(ph, curation_rows)

    has_prn_form = has_prn_form_master | has_prn_form_curation
    countries_with_prn_form = ph.loc[has_prn_form, "country_iso3"].nunique() if "country_iso3" in ph.columns else 0
    countries_added_by_curation = max(0, countries_with_prn_form - countries_with_prn_form_master)

    has_dtp3 = is_known_series(ph["dtp3_coverage"]) if "dtp3_coverage" in ph.columns else pd.Series(False, index=ph.index)
    countries_with_dtp3 = ph.loc[has_dtp3, "country_iso3"].nunique() if "country_iso3" in ph.columns else 0

    ap_info_ready = countries_with_ap >= min_countries_ap
    prn_formulation_ready = countries_with_prn_form >= min_countries_prn_form

    if ap_info_ready and prn_formulation_ready:
        decision = "PROCEED"
        recommendation = (
            f"aP intro data available for {countries_with_ap} countries (>={min_countries_ap}). "
            f"Prn formulation data available for {countries_with_prn_form} countries (>={min_countries_prn_form}). "
            "Can build both V1 and V2 exposure indices."
        )
        if countries_added_by_curation:
            recommendation += f" Curation contributed {countries_added_by_curation} country entries."
    elif ap_info_ready:
        decision = "PROCEED_V1_ONLY"
        recommendation = (
            f"aP intro data sufficient ({countries_with_ap} countries). "
            f"Prn formulation data insufficient ({countries_with_prn_form}/{min_countries_prn_form}). "
            "Proceed with V1 exposure index only. Prioritize Prn formulation curation for V2."
        )
    else:
        decision = "NEEDS_DOWNGRADE"
        recommendation = (
            f"aP intro data insufficient ({countries_with_ap}/{min_countries_ap} countries). "
            "Downgrade epidemiological claims or invest in manual WHO data curation."
        )

    country_summary: dict[str, dict[str, object]] = {}
    if "country_iso3" in ph.columns:
        for iso3 in ph["country_iso3"].dropna().astype(str).str.strip().unique():
            rows = ph[ph["country_iso3"].astype(str).str.strip() == iso3]
            indexer = rows.index
            has_prn_master_country = bool(has_prn_form_master[indexer].any())
            has_prn_curation_country = bool(has_prn_form_curation[indexer].any())
            country_summary[iso3] = {
                "has_ap_info": bool(has_ap_info[indexer].any()),
                "has_prn_form": has_prn_master_country or has_prn_curation_country,
                "has_prn_form_master": has_prn_master_country,
                "has_prn_form_curation": has_prn_curation_country,
                "has_dtp3": bool(has_dtp3[indexer].any()),
                "n_country_years": len(rows),
            }

    return {
        "assessment": "vaccine_variable_coverage",
        "date": pd.Timestamp.now().isoformat(),
        "countries_with_ap_info": countries_with_ap,
        "countries_with_prn_formulation": countries_with_prn_form,
        "countries_with_prn_formulation_from_master": countries_with_prn_form_master,
        "countries_with_prn_formulation_added_by_curation": countries_added_by_curation,
        "countries_with_dtp3": countries_with_dtp3,
        "threshold_ap": min_countries_ap,
        "threshold_prn_form": min_countries_prn_form,
        "ap_info_ready": ap_info_ready,
        "prn_formulation_ready": prn_formulation_ready,
        "decision": decision,
        "recommendation": recommendation,
        "prn_curation_file": str(curation_file) if prn_curation_path else "",
        "prn_curation_file_exists": bool(prn_curation_path) and curation_file.exists(),
        "prn_curation_rows_loaded": int(len(curation_rows)),
        "prn_curation_rows_applied": int(curation_rows_applied),
        "country_detail": country_summary,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Assess vaccine-variable coverage")
    parser.add_argument("--ph-master", default=str(WORKFLOW_DATA_ROOT / "epi" / "ap_exposure_index.tsv"))
    parser.add_argument("--min-ap", type=int, default=15)
    parser.add_argument("--min-prn-form", type=int, default=5)
    parser.add_argument(
        "--prn-curation",
        default=str(REPO_ROOT / "modules" / "public_health" / "inputs" / "curation" / "vaccine_formulation_curation.tsv"),
        help=(
            "Optional TSV with curated prn-in-vaccine rows. Accepts either the legacy "
            "overlay format (country_iso3 + prn_in_vaccine) or the richer "
            "vaccine_formulation_curation.tsv format."
        ),
    )
    parser.add_argument(
        "--out",
        default=str(WORKFLOW_DATA_ROOT / "checkpoints" / "vaccine_variable_coverage_report.json"),
    )
    args = parser.parse_args()

    report = assess_vaccine_variable_coverage(
        args.ph_master,
        args.min_ap,
        args.min_prn_form,
        prn_curation_path=args.prn_curation,
    )

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, default=str)

    print(f"\nVaccine-variable coverage — {report['decision']}")
    print(f"  Countries with aP info: {report['countries_with_ap_info']}/{report['threshold_ap']}")
    print(
        f"  Countries with Prn formulation: "
        f"{report['countries_with_prn_formulation']}/{report['threshold_prn_form']} "
        f"(master={report['countries_with_prn_formulation_from_master']}, "
        f"+curation={report['countries_with_prn_formulation_added_by_curation']})"
    )
    print(f"  Countries with DTP3: {report['countries_with_dtp3']}")
    print(f"  Recommendation: {report['recommendation']}")
