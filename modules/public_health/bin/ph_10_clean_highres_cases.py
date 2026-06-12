#!/usr/bin/env python3
"""Standardize monthly/weekly pertussis surveillance extracts for targeted validation analyses."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from ph_utils import (
    current_freeze_date,
    index_reporting_era_rows,
    match_reporting_era_row,
    normalize_text,
    read_xlsx_sheet_rows,
    reporting_era_post_flag,
    repo_root,
    project_module_data_root,
    project_workflow_root,
    rows_to_dicts,
)


SHEET_CONFIG = {
    "AU": {"country_iso3": "AUS", "country_name": "Australia"},
    "CN": {"country_iso3": "CHN", "country_name": "China"},
    "JP": {"country_iso3": "JPN", "country_name": "Japan"},
    "NZ": {"country_iso3": "NZL", "country_name": "New Zealand"},
    "SE": {"country_iso3": "SWE", "country_name": "Sweden"},
    "US": {"country_iso3": "USA", "country_name": "United States"},
}


OUTPUT_COLUMNS = [
    "country_iso3",
    "country_name",
    "source_sheet",
    "time_resolution",
    "date",
    "year",
    "month",
    "week",
    "interval_index",
    "fractional_year",
    "cases",
    "annual_cases",
    "share_of_annual_cases",
    "source_url",
    "source_file",
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
    "data_freeze_date",
    "notes",
]


PUBLIC_HEALTH_DATA_ROOT = project_module_data_root("public_health")
WORKFLOW_DATA_ROOT = project_workflow_root()


OVERLAP_COLUMNS = [
    "country_iso3",
    "country_name",
    "time_resolution",
    "first_year",
    "last_year",
    "n_rows",
    "n_country_years",
    "genomic_overlap_year_count",
    "genomic_overlap_years",
    "first_prn_detection_year",
    "first_local_origin_year",
    "reporting_era_record_iso3",
    "reporting_era_scope_type",
    "reporting_era_match_type",
    "reporting_era_confidence",
    "pcr_lab_guideline_year",
    "reporting_case_definition_change_year",
    "surveillance_platform_change_year",
    "notes",
]


def default_input_path() -> Path:
    return repo_root() / "modules/public_health/inputs/raw/who_cases/Pertussis incidence.xlsx"


def load_who_region_map(path: Path) -> dict[str, str]:
    rows = rows_to_dicts(read_xlsx_sheet_rows(path, sheet_name="Sheet1"))
    region_map: dict[str, str] = {}
    for row in rows:
        iso3 = normalize_text(row.get("ISO_3_CODE", "")).upper()
        region = normalize_text(row.get("WHO_REGION", "")).upper()
        if iso3 and region and iso3 not in region_map:
            region_map[iso3] = region
    return region_map


def annotate_reporting_era(
    combined: pd.DataFrame,
    reporting_era_rows: list[dict[str, str]],
    region_map: dict[str, str],
) -> pd.DataFrame:
    output = combined.copy()
    reporting_era_index = index_reporting_era_rows(reporting_era_rows)
    output["region_who"] = output["country_iso3"].map(lambda code: region_map.get(normalize_text(code).upper(), ""))

    matched_rows = output.apply(
        lambda row: match_reporting_era_row(
            normalize_text(row.get("country_iso3", "")),
            normalize_text(row.get("region_who", "")),
            reporting_era_index,
        ),
        axis=1,
    )
    output["reporting_era_record"] = matched_rows.map(lambda item: item[0])
    output["reporting_era_match_type"] = matched_rows.map(lambda item: item[1])
    output["reporting_era_record_iso3"] = output["reporting_era_record"].map(
        lambda row: normalize_text(row.get("iso3", ""))
    )
    output["reporting_era_scope_type"] = output["reporting_era_record"].map(
        lambda row: normalize_text(row.get("scope_type", ""))
    )
    output["reporting_era_confidence"] = output["reporting_era_record"].map(
        lambda row: normalize_text(row.get("confidence", ""))
    )
    output["pcr_lab_guideline_year"] = output["reporting_era_record"].map(
        lambda row: normalize_text(row.get("pcr_lab_guideline_year_min", ""))
    )
    output["reporting_case_definition_change_year"] = output["reporting_era_record"].map(
        lambda row: normalize_text(row.get("reporting_case_definition_change_year_min", ""))
    )
    output["surveillance_platform_change_year"] = output["reporting_era_record"].map(
        lambda row: normalize_text(row.get("surveillance_platform_change_year_min", ""))
    )
    output["post_pcr_lab_guideline_era"] = output.apply(
        lambda row: reporting_era_post_flag(
            int(row["year"]),
            row["pcr_lab_guideline_year"],
        ),
        axis=1,
    )
    output["post_reporting_case_definition_change_era"] = output.apply(
        lambda row: reporting_era_post_flag(
            int(row["year"]),
            row["reporting_case_definition_change_year"],
        ),
        axis=1,
    )
    output["post_surveillance_platform_change_era"] = output.apply(
        lambda row: reporting_era_post_flag(
            int(row["year"]),
            row["surveillance_platform_change_year"],
        ),
        axis=1,
    )
    return output.drop(columns=["region_who", "reporting_era_record"])


def load_highres_cases(
    input_path: Path,
    reporting_era_rows: list[dict[str, str]],
    region_map: dict[str, str],
) -> pd.DataFrame:
    workbook = pd.ExcelFile(input_path)
    rows = []
    for sheet_name, meta in SHEET_CONFIG.items():
        df = pd.read_excel(input_path, sheet_name=sheet_name)
        df = df.rename(columns=str.lower).copy()
        df["date"] = pd.to_datetime(df["date"], errors="coerce")
        df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
        df["month"] = pd.to_numeric(df.get("month"), errors="coerce").astype("Int64")
        df["week"] = pd.to_numeric(df.get("week"), errors="coerce").astype("Int64")
        df["cases"] = pd.to_numeric(df["cases"], errors="coerce")
        df["url"] = df["url"].fillna("")
        df = df[df["cases"].notna() & df["year"].notna()].copy()

        resolution = df["week"].notna().map({True: "weekly", False: "monthly"})
        df["time_resolution"] = resolution
        df["interval_index"] = df["week"].where(df["week"].notna(), df["month"]).astype("Int64")
        df["fractional_year"] = df["year"].astype(float) + (
            (df["week"].astype(float) - 0.5) / 52.0
        ).where(df["week"].notna(), (df["month"].astype(float) - 0.5) / 12.0)
        df["country_iso3"] = meta["country_iso3"]
        df["country_name"] = meta["country_name"]
        df["source_sheet"] = sheet_name
        df["source_file"] = input_path.name
        df["data_freeze_date"] = current_freeze_date()
        df["notes"] = df["time_resolution"].map(
            {
                "weekly": "high_resolution_national_surveillance_weekly",
                "monthly": "high_resolution_national_surveillance_monthly",
            }
        )
        rows.append(df)

    combined = pd.concat(rows, ignore_index=True)
    annual = (
        combined.groupby(["country_iso3", "year"], dropna=False)["cases"]
        .sum()
        .rename("annual_cases")
        .reset_index()
    )
    combined = combined.merge(annual, on=["country_iso3", "year"], how="left")
    combined["share_of_annual_cases"] = combined["cases"] / combined["annual_cases"]
    combined["date"] = combined["date"].dt.date.astype(str)
    combined = annotate_reporting_era(combined, reporting_era_rows, region_map)
    combined = combined.rename(columns={"url": "source_url"})
    return combined[OUTPUT_COLUMNS].sort_values(["country_iso3", "date", "interval_index"])


def build_overlap_summary(highres_df: pd.DataFrame) -> pd.DataFrame:
    ipw = pd.read_csv(WORKFLOW_DATA_ROOT / "epi" / "ipw_prevalence.tsv", sep="\t")
    origin_events = pd.read_csv(WORKFLOW_DATA_ROOT / "asr" / "origin_events.tsv", sep="\t")

    first_detection = (
        ipw.loc[ipw["n_prn_disrupted"].fillna(0) > 0, ["country_iso3", "year"]]
        .groupby("country_iso3", dropna=False)["year"]
        .min()
        .rename("first_prn_detection_year")
        .reset_index()
    )

    origin_country_rows = []
    subtree_dir = WORKFLOW_DATA_ROOT / "asr" / "event_subtrees"
    for path in sorted(subtree_dir.glob("origin_*.descendant_tips.tsv")):
        descendants = pd.read_csv(path, sep="\t")
        disrupted = descendants[descendants["observed_prn_state"] == "disrupted"].copy()
        disrupted = disrupted[disrupted["country_iso3"].notna() & disrupted["year"].notna()].copy()
        if disrupted.empty:
            continue
        origin_country_rows.append(
            disrupted.groupby("country_iso3", dropna=False)["year"].min().reset_index().assign(origin_id=path.stem.split(".")[0])
        )

    if origin_country_rows:
        origin_country = pd.concat(origin_country_rows, ignore_index=True)
        origin_country = origin_country.rename(columns={"year": "origin_first_year"})
        first_local_origin = (
            origin_country.groupby("country_iso3", dropna=False)["origin_first_year"]
            .min()
            .rename("first_local_origin_year")
            .reset_index()
        )
    else:
        first_local_origin = pd.DataFrame(columns=["country_iso3", "first_local_origin_year"])

    genomic_overlap = (
        ipw.loc[ipw["n_genomes_prn_interpretable"].fillna(0) > 0, ["country_iso3", "year"]]
        .drop_duplicates()
        .groupby("country_iso3", dropna=False)["year"]
        .agg(
            genomic_overlap_year_count="count",
            genomic_overlap_years=lambda values: ",".join(str(int(value)) for value in sorted(values)),
        )
        .reset_index()
    )

    summary = (
        highres_df.groupby(["country_iso3", "country_name", "time_resolution"], dropna=False)
        .agg(
            first_year=("year", "min"),
            last_year=("year", "max"),
            n_rows=("date", "count"),
            n_country_years=("year", "nunique"),
            reporting_era_record_iso3=("reporting_era_record_iso3", "first"),
            reporting_era_scope_type=("reporting_era_scope_type", "first"),
            reporting_era_match_type=("reporting_era_match_type", "first"),
            reporting_era_confidence=("reporting_era_confidence", "first"),
            pcr_lab_guideline_year=("pcr_lab_guideline_year", "first"),
            reporting_case_definition_change_year=("reporting_case_definition_change_year", "first"),
            surveillance_platform_change_year=("surveillance_platform_change_year", "first"),
        )
        .reset_index()
        .merge(genomic_overlap, on="country_iso3", how="left")
        .merge(first_detection, on="country_iso3", how="left")
        .merge(first_local_origin, on="country_iso3", how="left")
    )

    summary["notes"] = summary.apply(
        lambda row: (
            "usable_for_event_centered_validation"
            if (
                pd.notna(row["first_prn_detection_year"])
                and row["first_year"] <= row["first_prn_detection_year"] <= row["last_year"]
                and pd.notna(row["genomic_overlap_year_count"])
                and row["genomic_overlap_year_count"] >= 2
            )
            else (
                "usable_for_overlap_validation_only"
                if pd.notna(row["genomic_overlap_year_count"]) and row["genomic_overlap_year_count"] >= 2
                else "limited_genomic_overlap"
            )
        ),
        axis=1,
    )
    return summary[OVERLAP_COLUMNS].sort_values("country_iso3")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Standardize high-resolution pertussis surveillance workbook.")
    parser.add_argument("--input", type=Path, default=default_input_path())
    parser.add_argument(
        "--era-indicators",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_reporting_era_indicators.tsv",
    )
    parser.add_argument(
        "--who-region-source",
        type=Path,
        default=repo_root()
        / "modules/public_health/inputs/raw/vaccine_program_docs/Introduction of aP (acellular pertussis) vaccine 2026-21-03 23-11 UTC.xlsx",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_highres_cases.tsv",
    )
    parser.add_argument(
        "--overlap-out",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_highres_overlap_summary.tsv",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    reporting_era_rows = []
    if args.era_indicators.exists():
        reporting_era_rows = pd.read_csv(args.era_indicators, sep="\t", dtype=str).fillna("").to_dict("records")
    region_map = load_who_region_map(args.who_region_source)
    highres_df = load_highres_cases(args.input, reporting_era_rows, region_map)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    highres_df.to_csv(args.out, sep="\t", index=False)

    overlap_df = build_overlap_summary(highres_df)
    overlap_df.to_csv(args.overlap_out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
