#!/usr/bin/env python3
"""Construct workflow-native exposure-index outputs for country-year ecology models."""

from __future__ import annotations

import json
import re
import sys
from itertools import product
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from workflow.lib.project_paths import project_module_data_root

NOTE_YEAR_PATTERN = re.compile(r"(?P<field>first_any_ap_year|first_routine_ap_year)=(?P<value>na|\d{4})")
DEFAULT_FORMULATION_CURATION_PATH = project_module_data_root("public_health") / "inputs" / "curation" / "vaccine_formulation_curation.tsv"
DEFAULT_PRODUCT_METADATA_PATH = project_module_data_root("public_health") / "inputs" / "curation" / "vaccine_product_metadata.tsv"
FORMULATION_REQUIRED_COLUMNS = {
    "country_iso3",
    "year_start",
    "year_end",
    "ap_timing_anchor_year",
    "primary_series_formulation",
    "booster_formulation",
    "prn_in_vaccine_curated",
    "prn_in_vaccine_source_class",
    "formulation_confidence",
    "source_name",
    "source_url",
    "source_release_date",
    "notes",
}
PRODUCT_METADATA_REQUIRED_COLUMNS = {
    "country_iso3",
    "country_name",
    "year_start",
    "year_end",
    "exposure_role",
    "region_scope",
    "product_name",
    "manufacturer",
    "product_platform",
    "ap_prn_positive_fraction",
    "population_share",
    "share_basis",
    "evidence_confidence",
    "source_name",
    "source_url",
    "source_release_date",
    "notes",
}
PRN_VALUE_MAP = {"yes": 1.0, "mixed": 0.5, "no": 0.0}
AP_POSITIVE_PRIMARY = {
    "ap_prn_positive",
    "ap_prn_negative",
    "ap_or_mixed_unknown_prn",
    "mixed_brand_heterogeneous",
    "mixed_or_partial_ap",
}
PRODUCT_ROLE_PREFIX = {
    "routine_primary": "routine_primary",
    "routine_booster": "routine_booster",
    "maternal": "maternal",
}
FORMULATION_AP_SHARE_MAP = {
    "wp_only": 0.0,
    "wp_or_unknown": 0.0,
    "none_recorded": 0.0,
    "ap_prn_positive": 1.0,
    "ap_prn_negative": 1.0,
    "dtap_or_tdap_prn_positive": 1.0,
    "mixed_brand_heterogeneous": 1.0,
    "ap_or_mixed_unknown_prn": 1.0,
    "mixed_or_partial_ap": 0.5,
}
FORMULATION_PRN_SHARE_MAP = {
    "wp_only": 0.0,
    "wp_or_unknown": 0.0,
    "none_recorded": 0.0,
    "ap_prn_positive": 1.0,
    "ap_prn_negative": 0.0,
    "dtap_or_tdap_prn_positive": 1.0,
    "mixed_brand_heterogeneous": 0.5,
    "ap_or_mixed_unknown_prn": 0.5,
    "mixed_or_partial_ap": 0.5,
}
PRODUCT_PLATFORM_AP_SHARE_MAP = {
    "wp": 0.0,
    "ap_prn_positive": 1.0,
    "ap_prn_negative": 1.0,
    "ap_mixed": 1.0,
}
PRODUCT_CONFIDENCE_WEIGHT_MAP = {
    "high": 1.0,
    "medium": 0.75,
    "low": 0.5,
    "unknown": 0.25,
    "": 0.25,
}


def normalize_text(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def coerce_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def resolve_optional_input_path(path: str | None, default_path: Path) -> str:
    text = normalize_text(path)
    if text:
        return text
    return str(default_path)


def canonical_country_name(values: pd.Series) -> str:
    clean = values.dropna().astype(str).str.strip()
    clean = clean[clean.str.len().gt(0)]
    if clean.empty:
        return ""
    counts = clean.value_counts()
    top_count = int(counts.max())
    candidates = counts[counts.eq(top_count)].index.tolist()
    return min(candidates, key=lambda value: (len(value), value))


def parse_year_from_notes(notes: str, field_name: str) -> float:
    if not isinstance(notes, str):
        return np.nan
    for match in NOTE_YEAR_PATTERN.finditer(notes):
        if match.group("field") == field_name:
            value = match.group("value")
            return np.nan if value == "na" else float(value)
    return np.nan


def expand_program_metadata(program_metadata_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = pd.read_csv(program_metadata_path, sep="\t", dtype=str)
    metadata["year_start"] = coerce_numeric(metadata.get("year_start", pd.Series(dtype=str)))
    metadata["year_end"] = coerce_numeric(metadata.get("year_end", pd.Series(dtype=str)))
    metadata["first_any_ap_year"] = metadata.get("notes", pd.Series(dtype=str)).map(
        lambda value: parse_year_from_notes(value, "first_any_ap_year")
    )
    metadata["first_routine_ap_year"] = metadata.get("notes", pd.Series(dtype=str)).map(
        lambda value: parse_year_from_notes(value, "first_routine_ap_year")
    )
    metadata["booster_flag"] = metadata.get("booster_schedule", pd.Series(dtype=str)).fillna("").str.strip().ne("").astype(int)

    country_summary = (
        metadata.groupby("country_iso3", dropna=False)
        .agg(
            country_name_program_canonical=("country_name", canonical_country_name),
            first_any_ap_year=("first_any_ap_year", "min"),
            first_routine_ap_year=("first_routine_ap_year", "min"),
        )
        .reset_index()
    )

    annual_rows: list[dict[str, object]] = []
    for row in metadata.itertuples(index=False):
        if np.isnan(row.year_start) or np.isnan(row.year_end):
            continue
        for year_value in range(int(row.year_start), int(row.year_end) + 1):
            annual_rows.append(
                {
                    "country_iso3": row.country_iso3,
                    "country_name_program": row.country_name,
                    "year": year_value,
                    "program_metadata_vaccine_program_type": row.vaccine_program_type,
                    "program_metadata_acellular_vs_whole_cell": row.acellular_vs_whole_cell,
                    "program_metadata_prn_in_vaccine": row.prn_in_vaccine,
                    "program_metadata_booster_flag": row.booster_flag,
                    "first_any_ap_year": row.first_any_ap_year,
                    "first_routine_ap_year": row.first_routine_ap_year,
                }
            )
    annual_frame = pd.DataFrame.from_records(annual_rows)
    return annual_frame, country_summary


def load_formulation_curation(path: str | None) -> pd.DataFrame:
    columns = sorted(FORMULATION_REQUIRED_COLUMNS)
    if not path:
        return pd.DataFrame(columns=columns)
    curation_path = Path(path)
    if not curation_path.exists():
        return pd.DataFrame(columns=columns)

    frame = pd.read_csv(curation_path, sep="\t", dtype=str).fillna("")
    missing = FORMULATION_REQUIRED_COLUMNS - set(frame.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"Formulation curation file is missing required columns: {missing_text}")

    frame["country_iso3"] = frame["country_iso3"].map(lambda value: normalize_text(value).upper())
    frame["year_start"] = coerce_numeric(frame["year_start"])
    frame["year_end"] = coerce_numeric(frame["year_end"])
    frame["ap_timing_anchor_year"] = coerce_numeric(frame["ap_timing_anchor_year"])
    return frame


def load_product_metadata(path: str | None) -> pd.DataFrame:
    columns = sorted(PRODUCT_METADATA_REQUIRED_COLUMNS)
    if not path:
        return pd.DataFrame(columns=columns)
    metadata_path = Path(path)
    if not metadata_path.exists():
        return pd.DataFrame(columns=columns)

    frame = pd.read_csv(metadata_path, sep="\t", dtype=str).fillna("")
    missing = PRODUCT_METADATA_REQUIRED_COLUMNS - set(frame.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"Product metadata file is missing required columns: {missing_text}")

    frame["country_iso3"] = frame["country_iso3"].map(lambda value: normalize_text(value).upper())
    frame["year_start"] = coerce_numeric(frame["year_start"])
    frame["year_end"] = coerce_numeric(frame["year_end"])
    frame["ap_prn_positive_fraction"] = coerce_numeric(frame["ap_prn_positive_fraction"])
    frame["population_share"] = coerce_numeric(frame["population_share"])
    frame["exposure_role"] = frame["exposure_role"].map(normalize_text)
    frame["product_platform"] = frame["product_platform"].map(normalize_text)
    frame["evidence_confidence"] = frame["evidence_confidence"].map(lambda value: normalize_text(value).lower())
    return frame


def infer_primary_series_formulation(vaccine_program_type: str, acellular_vs_whole_cell: str) -> str:
    vaccine_program_type = normalize_text(vaccine_program_type)
    acellular_vs_whole_cell = normalize_text(acellular_vs_whole_cell)
    if vaccine_program_type == "pre_ap_introduction_or_whole_cell" or acellular_vs_whole_cell == "whole_cell_or_unknown":
        return "wp_only"
    if vaccine_program_type == "ap_targeted_or_partial_nonroutine" or acellular_vs_whole_cell == "targeted_or_partial_ap":
        return "mixed_or_partial_ap"
    if vaccine_program_type == "ap_introduced_routine_or_mixed" or acellular_vs_whole_cell == "mixed_or_acellular":
        return "ap_or_mixed_unknown_prn"
    return "unknown"


def infer_booster_formulation(primary_formulation: str, booster_flag: float) -> str:
    if int(booster_flag or 0) == 0:
        return "none_recorded"
    if primary_formulation == "wp_only":
        return "wp_or_unknown"
    if primary_formulation == "mixed_or_partial_ap":
        return "mixed_or_partial_ap"
    if primary_formulation == "ap_or_mixed_unknown_prn":
        return "ap_or_mixed_unknown_prn"
    return "unknown"


def infer_prn_value(primary_formulation: str) -> str:
    if primary_formulation == "wp_only":
        return "no"
    if primary_formulation == "mixed_or_partial_ap":
        return "mixed"
    return "unknown"


def program_phase_supports_ap(vaccine_program_type: object) -> bool:
    if pd.isna(vaccine_program_type):
        return np.nan
    text = normalize_text(vaccine_program_type)
    if not text:
        return np.nan
    return text in {
        "ap_introduced_routine_or_mixed",
        "ap_targeted_or_partial_nonroutine",
    }


def formulation_class_supports_ap(program_formulation_class: object) -> bool:
    return normalize_text(program_formulation_class) in {
        "routine_ap_prn_positive",
        "routine_ap_prn_negative",
        "routine_ap_mixed",
        "routine_ap_unknown",
    }


def classify_program_formulation(
    primary_series_formulation: object,
    prn_in_vaccine_curated: object,
    vaccine_program_type_effective: object,
) -> str:
    primary = normalize_text(primary_series_formulation)
    prn_value = normalize_text(prn_in_vaccine_curated).lower()
    vaccine_program_type = normalize_text(vaccine_program_type_effective)

    if primary == "wp_only":
        if prn_value == "mixed":
            return "routine_ap_mixed"
        return "wp_only_or_pre_ap"
    if prn_value == "yes":
        return "routine_ap_prn_positive"
    if prn_value == "no":
        return "routine_ap_prn_negative"
    if prn_value == "mixed":
        return "routine_ap_mixed"
    if primary in AP_POSITIVE_PRIMARY:
        return "routine_ap_unknown"
    if vaccine_program_type == "pre_ap_introduction_or_whole_cell":
        return "wp_only_or_pre_ap"
    if vaccine_program_type in {
        "ap_introduced_routine_or_mixed",
        "ap_targeted_or_partial_nonroutine",
    }:
        return "routine_ap_unknown"
    return "unknown"


def determine_precedence_rule(
    prn_in_vaccine_source_class: object,
    program_formulation_conflict: object,
    ap_exposure_v2_available: object,
) -> str:
    source_class = normalize_text(prn_in_vaccine_source_class)
    has_v2 = bool(ap_exposure_v2_available)
    if bool(program_formulation_conflict):
        return "curated_formulation_preferred_over_program_phase"
    if source_class == "program_phase_inferred":
        return "program_phase_inferred_fallback"
    if has_v2:
        return "curated_formulation_applied"
    return "program_metadata_only_fallback"


def formulation_to_ap_share(formulation: object) -> float:
    text = normalize_text(formulation)
    value = FORMULATION_AP_SHARE_MAP.get(text)
    return np.nan if value is None else float(value)


def formulation_to_prn_share(formulation: object) -> float:
    text = normalize_text(formulation)
    value = FORMULATION_PRN_SHARE_MAP.get(text)
    return np.nan if value is None else float(value)


def prn_value_to_fraction(value: object) -> float:
    text = normalize_text(value).lower()
    mapped = PRN_VALUE_MAP.get(text)
    return np.nan if mapped is None else float(mapped)


def product_platform_to_ap_share(platform: object) -> float:
    text = normalize_text(platform)
    value = PRODUCT_PLATFORM_AP_SHARE_MAP.get(text)
    return np.nan if value is None else float(value)


def product_platform_to_prn_positive_share(platform: object, ap_prn_positive_fraction: object) -> float:
    text = normalize_text(platform)
    if text == "ap_mixed":
        return float(ap_prn_positive_fraction)
    if text == "ap_prn_positive":
        return 1.0
    if text == "ap_prn_negative":
        return 0.0
    if text == "wp":
        return 0.0
    return np.nan


def confidence_to_weight(value: object) -> float:
    text = normalize_text(value).lower()
    return float(PRODUCT_CONFIDENCE_WEIGHT_MAP.get(text, PRODUCT_CONFIDENCE_WEIGHT_MAP["unknown"]))


def join_unique_text(values: pd.Series) -> str:
    seen: list[str] = []
    for value in values.fillna("").astype(str):
        text = value.strip()
        if not text or text in seen:
            continue
        seen.append(text)
    return "; ".join(seen)


def blend_raw_with_fallback(raw_value: object, fallback_value: object, blend_weight: object) -> float:
    raw = np.nan if pd.isna(raw_value) else float(raw_value)
    fallback = np.nan if pd.isna(fallback_value) else float(fallback_value)
    weight = 0.0 if pd.isna(blend_weight) else float(blend_weight)
    weight = min(1.0, max(0.0, weight))
    if np.isnan(raw) and np.isnan(fallback):
        return np.nan
    if np.isnan(raw):
        return fallback
    if np.isnan(fallback):
        return raw
    return weight * raw + (1.0 - weight) * fallback


def apply_formulation_curation(merged: pd.DataFrame, formulation_curation: pd.DataFrame) -> pd.DataFrame:
    output = merged.copy()
    output["country_iso3"] = output.get("country_iso3", pd.Series(dtype=str)).fillna("").astype(str).str.upper()
    output["year"] = coerce_numeric(output.get("year", pd.Series(dtype=str)))

    curated_columns = [
        "ap_timing_anchor_year",
        "primary_series_formulation",
        "booster_formulation",
        "prn_in_vaccine_curated",
        "prn_in_vaccine_source_class",
        "formulation_confidence",
        "formulation_source_name",
        "formulation_source_url",
        "formulation_source_release_date",
        "formulation_notes",
    ]
    for column in curated_columns:
        output[column] = ""

    if not formulation_curation.empty:
        for row in formulation_curation.itertuples(index=False):
            row_mask = output["country_iso3"].eq(row.country_iso3)
            if not np.isnan(row.year_start):
                row_mask &= output["year"].ge(float(row.year_start))
            if not np.isnan(row.year_end):
                row_mask &= output["year"].le(float(row.year_end))
            if not bool(row_mask.any()):
                continue

            output.loc[row_mask, "ap_timing_anchor_year"] = "" if np.isnan(row.ap_timing_anchor_year) else f"{int(row.ap_timing_anchor_year)}"
            output.loc[row_mask, "primary_series_formulation"] = normalize_text(row.primary_series_formulation)
            output.loc[row_mask, "booster_formulation"] = normalize_text(row.booster_formulation)
            output.loc[row_mask, "prn_in_vaccine_curated"] = normalize_text(row.prn_in_vaccine_curated).lower()
            output.loc[row_mask, "prn_in_vaccine_source_class"] = normalize_text(row.prn_in_vaccine_source_class)
            output.loc[row_mask, "formulation_confidence"] = normalize_text(row.formulation_confidence)
            output.loc[row_mask, "formulation_source_name"] = normalize_text(row.source_name)
            output.loc[row_mask, "formulation_source_url"] = normalize_text(row.source_url)
            output.loc[row_mask, "formulation_source_release_date"] = normalize_text(row.source_release_date)
            output.loc[row_mask, "formulation_notes"] = normalize_text(row.notes)

    output["acellular_vs_whole_cell_effective"] = output.get(
        "program_metadata_acellular_vs_whole_cell",
        pd.Series(index=output.index, dtype=str),
    ).fillna(output.get("acellular_vs_whole_cell", pd.Series(index=output.index, dtype=str)))

    fallback_primary = output.apply(
        lambda row: infer_primary_series_formulation(
            row.get("vaccine_program_type_effective", ""),
            row.get("acellular_vs_whole_cell_effective", ""),
        ),
        axis=1,
    )
    output["primary_series_formulation"] = output["primary_series_formulation"].where(
        output["primary_series_formulation"].astype(str).str.len().gt(0),
        fallback_primary,
    )

    fallback_booster = output.apply(
        lambda row: infer_booster_formulation(
            normalize_text(row.get("primary_series_formulation", "")),
            float(row.get("booster_flag", 0) or 0),
        ),
        axis=1,
    )
    output["booster_formulation"] = output["booster_formulation"].where(
        output["booster_formulation"].astype(str).str.len().gt(0),
        fallback_booster,
    )

    fallback_prn = output.apply(
        lambda row: infer_prn_value(normalize_text(row.get("primary_series_formulation", ""))),
        axis=1,
    )
    output["prn_in_vaccine_curated"] = output["prn_in_vaccine_curated"].where(
        output["prn_in_vaccine_curated"].astype(str).str.len().gt(0),
        fallback_prn,
    )

    output["prn_in_vaccine_source_class"] = output["prn_in_vaccine_source_class"].where(
        output["prn_in_vaccine_source_class"].astype(str).str.len().gt(0),
        np.where(output["prn_in_vaccine_curated"].isin(["yes", "mixed", "no"]), "program_phase_inferred", "unknown"),
    )
    output["formulation_confidence"] = output["formulation_confidence"].where(
        output["formulation_confidence"].astype(str).str.len().gt(0),
        np.where(output["prn_in_vaccine_source_class"].eq("program_phase_inferred"), "medium", "unknown"),
    )

    output["ap_timing_anchor_year"] = coerce_numeric(output["ap_timing_anchor_year"])
    output["ap_timing_anchor_year_effective"] = output["ap_timing_anchor_year"].fillna(
        output["first_routine_ap_year"].fillna(output["first_any_ap_year"])
    )
    output["prn_formulation_component"] = output["prn_in_vaccine_curated"].map(PRN_VALUE_MAP)
    output["ap_exposure_v2_available"] = output["prn_in_vaccine_curated"].isin(["yes", "mixed", "no"])
    output["program_supports_ap"] = output["vaccine_program_type_effective"].map(program_phase_supports_ap)
    output["program_formulation_class"] = output.apply(
        lambda row: classify_program_formulation(
            row.get("primary_series_formulation", ""),
            row.get("prn_in_vaccine_curated", ""),
            row.get("vaccine_program_type_effective", ""),
        ),
        axis=1,
    )
    output["resolved_supports_ap"] = output["program_formulation_class"].map(formulation_class_supports_ap)
    output["program_formulation_conflict"] = (
        output["program_supports_ap"].notna()
        & (
            output["program_supports_ap"].astype(bool)
            != output["resolved_supports_ap"].astype(bool)
        )
    ) & output["program_formulation_class"].isin(
        [
            "wp_only_or_pre_ap",
            "routine_ap_prn_positive",
            "routine_ap_prn_negative",
            "routine_ap_mixed",
            "routine_ap_unknown",
        ]
    )
    output["exposure_precedence_rule"] = output.apply(
        lambda row: determine_precedence_rule(
            row.get("prn_in_vaccine_source_class", ""),
            row.get("program_formulation_conflict", False),
            row.get("ap_exposure_v2_available", False),
        ),
        axis=1,
    )
    return output


def fallback_role_ap_share(row: pd.Series, role: str) -> float:
    if role == "routine_primary":
        return formulation_to_ap_share(row.get("primary_series_formulation", ""))
    if role == "routine_booster":
        return formulation_to_ap_share(row.get("booster_formulation", ""))
    return np.nan


def fallback_role_prn_positive_share(row: pd.Series, role: str) -> float:
    if role == "routine_primary":
        ap_share = fallback_role_ap_share(row, role)
        prn_share = prn_value_to_fraction(row.get("prn_in_vaccine_curated", ""))
    elif role == "routine_booster":
        ap_share = fallback_role_ap_share(row, role)
        prn_share = formulation_to_prn_share(row.get("booster_formulation", ""))
    else:
        return np.nan
    if np.isnan(ap_share) or np.isnan(prn_share):
        return np.nan
    return float(ap_share) * float(prn_share)


def apply_product_metadata(merged: pd.DataFrame, product_metadata: pd.DataFrame) -> pd.DataFrame:
    output = merged.copy()

    for prefix in PRODUCT_ROLE_PREFIX.values():
        numeric_defaults = {
            f"product_{prefix}_observed_share": np.nan,
            f"product_{prefix}_weighted_confidence": np.nan,
            f"product_{prefix}_bridge_weight": np.nan,
            f"product_{prefix}_wp_share_raw": np.nan,
            f"product_{prefix}_ap_share_raw": np.nan,
            f"product_{prefix}_ap_prn_positive_share_raw": np.nan,
            f"product_{prefix}_ap_prn_negative_share_raw": np.nan,
            f"product_{prefix}_ap_share_conf_weighted": np.nan,
            f"product_{prefix}_ap_prn_positive_share_conf_weighted": np.nan,
            f"product_{prefix}_ap_prn_negative_share_conf_weighted": np.nan,
        }
        for column, default_value in numeric_defaults.items():
            output[column] = default_value
        output[f"product_{prefix}_product_list"] = ""
        output[f"product_{prefix}_share_basis"] = ""
        output[f"product_{prefix}_source_name"] = ""

    if not product_metadata.empty:
        annual_rows: list[dict[str, object]] = []
        for row in product_metadata.itertuples(index=False):
            if np.isnan(row.year_start) or np.isnan(row.year_end):
                continue
            ap_share = product_platform_to_ap_share(row.product_platform)
            prn_positive_share = product_platform_to_prn_positive_share(
                row.product_platform,
                row.ap_prn_positive_fraction,
            )
            if np.isnan(ap_share) or np.isnan(prn_positive_share):
                continue
            prn_negative_share = max(0.0, float(ap_share) - float(prn_positive_share))
            confidence_weight = confidence_to_weight(row.evidence_confidence)
            for year_value in range(int(row.year_start), int(row.year_end) + 1):
                annual_rows.append(
                    {
                        "country_iso3": row.country_iso3,
                        "year": year_value,
                        "exposure_role": normalize_text(row.exposure_role),
                        "observed_share": float(row.population_share),
                        "confidence_weight": confidence_weight,
                        "wp_share_raw": float(row.population_share) * (1.0 - float(ap_share)),
                        "ap_share_raw": float(row.population_share) * float(ap_share),
                        "ap_prn_positive_share_raw": float(row.population_share) * float(prn_positive_share),
                        "ap_prn_negative_share_raw": float(row.population_share) * float(prn_negative_share),
                        "product_name": normalize_text(row.product_name),
                        "share_basis": normalize_text(row.share_basis),
                        "source_name": normalize_text(row.source_name),
                    }
                )

        annual = pd.DataFrame.from_records(annual_rows)
        if not annual.empty:
            grouped = (
                annual.groupby(["country_iso3", "year", "exposure_role"], dropna=False)
                .agg(
                    observed_share=("observed_share", "sum"),
                    weighted_confidence_numerator=("confidence_weight", lambda values: np.nan),
                    wp_share_raw=("wp_share_raw", "sum"),
                    ap_share_raw=("ap_share_raw", "sum"),
                    ap_prn_positive_share_raw=("ap_prn_positive_share_raw", "sum"),
                    ap_prn_negative_share_raw=("ap_prn_negative_share_raw", "sum"),
                    product_list=("product_name", join_unique_text),
                    share_basis=("share_basis", join_unique_text),
                    source_name=("source_name", join_unique_text),
                )
                .reset_index()
            )
            weighted_confidence = (
                annual.assign(weighted_confidence_component=annual["observed_share"] * annual["confidence_weight"])
                .groupby(["country_iso3", "year", "exposure_role"], dropna=False)[
                    ["observed_share", "weighted_confidence_component"]
                ]
                .sum()
                .reset_index()
            )
            weighted_confidence["weighted_confidence"] = np.where(
                weighted_confidence["observed_share"].gt(0),
                weighted_confidence["weighted_confidence_component"] / weighted_confidence["observed_share"],
                np.nan,
            )
            grouped = grouped.drop(columns=["weighted_confidence_numerator"]).merge(
                weighted_confidence[["country_iso3", "year", "exposure_role", "weighted_confidence"]],
                on=["country_iso3", "year", "exposure_role"],
                how="left",
            )

            for role, prefix in PRODUCT_ROLE_PREFIX.items():
                role_frame = grouped.loc[grouped["exposure_role"].eq(role)].copy()
                if role_frame.empty:
                    continue
                role_frame = role_frame.rename(
                    columns={
                        "observed_share": f"product_{prefix}_observed_share",
                        "weighted_confidence": f"product_{prefix}_weighted_confidence",
                        "wp_share_raw": f"product_{prefix}_wp_share_raw",
                        "ap_share_raw": f"product_{prefix}_ap_share_raw",
                        "ap_prn_positive_share_raw": f"product_{prefix}_ap_prn_positive_share_raw",
                        "ap_prn_negative_share_raw": f"product_{prefix}_ap_prn_negative_share_raw",
                        "product_list": f"product_{prefix}_product_list",
                        "share_basis": f"product_{prefix}_share_basis",
                        "source_name": f"product_{prefix}_source_name",
                    }
                )
                role_frame = role_frame.drop(columns=["exposure_role"]).set_index(["country_iso3", "year"])
                output = output.set_index(["country_iso3", "year"])
                shared_index = output.index.intersection(role_frame.index)
                if shared_index.empty:
                    output = output.reset_index()
                    continue
                for column in role_frame.columns:
                    output.loc[shared_index, column] = role_frame.loc[shared_index, column]
                output = output.reset_index()

    for role, prefix in PRODUCT_ROLE_PREFIX.items():
        observed_share_column = f"product_{prefix}_observed_share"
        weighted_confidence_column = f"product_{prefix}_weighted_confidence"
        bridge_weight_column = f"product_{prefix}_bridge_weight"
        fallback_ap_share = output.apply(lambda row: fallback_role_ap_share(row, role), axis=1)
        fallback_prn_positive_share = output.apply(
            lambda row: fallback_role_prn_positive_share(row, role),
            axis=1,
        )
        output[bridge_weight_column] = (
            output[observed_share_column].fillna(0).astype(float)
            * output[weighted_confidence_column].fillna(0).astype(float)
        ).clip(lower=0.0, upper=1.0)
        output[f"product_{prefix}_ap_share_conf_weighted"] = output.apply(
            lambda row: blend_raw_with_fallback(
                row.get(f"product_{prefix}_ap_share_raw"),
                fallback_ap_share.loc[row.name],
                row.get(bridge_weight_column),
            ),
            axis=1,
        )
        output[f"product_{prefix}_ap_prn_positive_share_conf_weighted"] = output.apply(
            lambda row: blend_raw_with_fallback(
                row.get(f"product_{prefix}_ap_prn_positive_share_raw"),
                fallback_prn_positive_share.loc[row.name],
                row.get(bridge_weight_column),
            ),
            axis=1,
        )
        output[f"product_{prefix}_ap_prn_negative_share_conf_weighted"] = np.where(
            output[f"product_{prefix}_ap_share_conf_weighted"].notna()
            & output[f"product_{prefix}_ap_prn_positive_share_conf_weighted"].notna(),
            np.maximum(
                0.0,
                output[f"product_{prefix}_ap_share_conf_weighted"]
                - output[f"product_{prefix}_ap_prn_positive_share_conf_weighted"],
            ),
            np.nan,
        )

    output["routine_primary_ap_prn_positive_coverage_proxy"] = (
        output["dtp3_coverage"] / 100.0
    ) * output["product_routine_primary_ap_prn_positive_share_conf_weighted"]
    output["routine_primary_ap_prn_negative_coverage_proxy"] = (
        output["dtp3_coverage"] / 100.0
    ) * output["product_routine_primary_ap_prn_negative_share_conf_weighted"]
    output["routine_primary_wp_coverage_proxy"] = (
        output["dtp3_coverage"] / 100.0
    ) * np.maximum(0.0, 1.0 - output["product_routine_primary_ap_share_conf_weighted"].fillna(0.0))
    return output


def select_primary_value(values: list[float]) -> float:
    ordered = sorted(float(value) for value in values)
    return ordered[len(ordered) // 2]


def standardize_series(series: pd.Series) -> pd.Series:
    valid = series.dropna()
    if valid.empty:
        return pd.Series(np.nan, index=series.index)
    standard_deviation = valid.std(ddof=0)
    if standard_deviation == 0 or np.isnan(standard_deviation):
        return pd.Series(0.0, index=series.index)
    return (series - valid.mean()) / standard_deviation


def load_vaccine_variable_coverage_decision(ph_master_path: str) -> tuple[str, str]:
    default = ("UNKNOWN", "Vaccine-variable coverage report not found")
    coverage_path = (
        Path(ph_master_path).resolve().parents[3]
        / "outputs"
        / "workflow"
        / "checkpoints"
        / "vaccine_variable_coverage_report.json"
    )
    if not coverage_path.exists():
        return default
    with coverage_path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    recommendation = payload.get("recommendation", payload.get("summary", ""))
    return str(payload.get("decision", "UNKNOWN")), str(recommendation)


def write_formulation_outputs(base_frame: pd.DataFrame, output_index_path: str) -> tuple[Path, Path]:
    output_path = Path(output_index_path)
    country_year_path = output_path.with_name("formulation_curation_country_year.tsv")
    summary_path = output_path.with_name("formulation_curation_summary.tsv")

    country_year_columns = [
        "country_iso3",
        "country_name",
        "year",
        "vaccine_program_type_effective",
        "primary_series_formulation",
        "booster_formulation",
        "prn_in_vaccine_curated",
        "prn_in_vaccine_source_class",
        "formulation_confidence",
        "ap_timing_anchor_year_effective",
        "program_formulation_class",
        "program_formulation_conflict",
        "exposure_precedence_rule",
        "ap_exposure_v2_available",
        "formulation_source_name",
        "formulation_source_url",
        "formulation_source_release_date",
        "formulation_notes",
    ]
    base_frame[country_year_columns].to_csv(country_year_path, sep="\t", index=False)

    summary = (
        base_frame.groupby("country_iso3", dropna=False)
        .agg(
            country_name=("country_name", canonical_country_name),
            n_country_years=("year", "count"),
            n_with_known_prn=("ap_exposure_v2_available", "sum"),
            first_year=("year", "min"),
            last_year=("year", "max"),
            dominant_primary_series_formulation=("primary_series_formulation", lambda values: values.mode().iloc[0] if not values.mode().empty else ""),
            dominant_booster_formulation=("booster_formulation", lambda values: values.mode().iloc[0] if not values.mode().empty else ""),
            dominant_prn_in_vaccine_curated=("prn_in_vaccine_curated", lambda values: values.mode().iloc[0] if not values.mode().empty else ""),
            dominant_formulation_confidence=("formulation_confidence", lambda values: values.mode().iloc[0] if not values.mode().empty else ""),
        )
        .reset_index()
    )
    summary["known_prn_fraction"] = np.where(
        summary["n_country_years"].gt(0),
        summary["n_with_known_prn"] / summary["n_country_years"],
        np.nan,
    )
    summary.to_csv(summary_path, sep="\t", index=False)
    return country_year_path, summary_path


def write_product_metadata_outputs(base_frame: pd.DataFrame, output_index_path: str) -> tuple[Path, Path]:
    output_path = Path(output_index_path)
    country_year_path = output_path.with_name("product_metadata_country_year.tsv")
    summary_path = output_path.with_name("product_metadata_summary.tsv")

    country_year_columns = [
        "country_iso3",
        "country_name",
        "year",
        "primary_series_formulation",
        "booster_formulation",
        "prn_in_vaccine_curated",
        "product_routine_primary_observed_share",
        "product_routine_primary_weighted_confidence",
        "product_routine_primary_bridge_weight",
        "product_routine_primary_wp_share_raw",
        "product_routine_primary_ap_share_raw",
        "product_routine_primary_ap_prn_positive_share_raw",
        "product_routine_primary_ap_prn_negative_share_raw",
        "product_routine_primary_ap_share_conf_weighted",
        "product_routine_primary_ap_prn_positive_share_conf_weighted",
        "product_routine_primary_ap_prn_negative_share_conf_weighted",
        "product_routine_primary_product_list",
        "product_routine_primary_share_basis",
        "product_routine_booster_observed_share",
        "product_routine_booster_weighted_confidence",
        "product_routine_booster_bridge_weight",
        "product_routine_booster_ap_share_conf_weighted",
        "product_routine_booster_ap_prn_positive_share_conf_weighted",
        "product_routine_booster_ap_prn_negative_share_conf_weighted",
        "product_routine_booster_product_list",
        "product_routine_booster_share_basis",
        "product_maternal_observed_share",
        "product_maternal_weighted_confidence",
        "product_maternal_bridge_weight",
        "product_maternal_ap_share_conf_weighted",
        "product_maternal_ap_prn_positive_share_conf_weighted",
        "product_maternal_ap_prn_negative_share_conf_weighted",
        "product_maternal_product_list",
        "product_maternal_share_basis",
        "routine_primary_ap_prn_positive_coverage_proxy",
        "routine_primary_ap_prn_negative_coverage_proxy",
        "routine_primary_wp_coverage_proxy",
    ]
    available_country_year_columns = [column for column in country_year_columns if column in base_frame.columns]
    base_frame[available_country_year_columns].to_csv(country_year_path, sep="\t", index=False)

    summary = (
        base_frame.groupby("country_iso3", dropna=False)
        .agg(
            country_name=("country_name", canonical_country_name),
            n_country_years=("year", "count"),
            n_years_with_primary_product_metadata=("product_routine_primary_observed_share", lambda values: int(values.fillna(0).gt(0).sum())),
            n_years_with_booster_product_metadata=("product_routine_booster_observed_share", lambda values: int(values.fillna(0).gt(0).sum())),
            n_years_with_maternal_product_metadata=("product_maternal_observed_share", lambda values: int(values.fillna(0).gt(0).sum())),
            mean_primary_prn_positive_share=("product_routine_primary_ap_prn_positive_share_conf_weighted", "mean"),
            mean_primary_prn_negative_share=("product_routine_primary_ap_prn_negative_share_conf_weighted", "mean"),
            mean_primary_ap_share=("product_routine_primary_ap_share_conf_weighted", "mean"),
            mean_booster_prn_positive_share=("product_routine_booster_ap_prn_positive_share_conf_weighted", "mean"),
            mean_maternal_prn_positive_share=("product_maternal_ap_prn_positive_share_conf_weighted", "mean"),
            dominant_primary_products=("product_routine_primary_product_list", lambda values: values.mode().iloc[0] if not values.mode().empty else ""),
            dominant_primary_share_basis=("product_routine_primary_share_basis", lambda values: values.mode().iloc[0] if not values.mode().empty else ""),
        )
        .reset_index()
    )
    summary.to_csv(summary_path, sep="\t", index=False)
    return country_year_path, summary_path


def write_precedence_audit(base_frame: pd.DataFrame, output_index_path: str) -> Path:
    output_path = Path(output_index_path)
    audit_path = output_path.with_name("formulation_precedence_audit.tsv")
    audit_columns = [
        "country_iso3",
        "country_name",
        "year",
        "vaccine_program_type_effective",
        "primary_series_formulation",
        "booster_formulation",
        "prn_in_vaccine_curated",
        "formulation_confidence",
        "program_supports_ap",
        "program_formulation_class",
        "resolved_supports_ap",
        "program_formulation_conflict",
        "exposure_precedence_rule",
        "prn_in_vaccine_source_class",
        "formulation_source_name",
        "formulation_source_url",
        "formulation_notes",
    ]
    base_frame[audit_columns].sort_values(["program_formulation_conflict", "country_iso3", "year"], ascending=[False, True, True]).to_csv(
        audit_path,
        sep="\t",
        index=False,
    )
    return audit_path


def build_index(
    ph_master_path: str,
    program_metadata_path: str,
    formulation_curation_path: str | None,
    product_metadata_path: str | None,
    output_index_path: str,
    output_figure_path: str,
    exposure_version: str,
    lambda_range: list[float],
    gamma_range: list[float],
) -> pd.DataFrame:
    ph_master = pd.read_csv(ph_master_path, sep="\t", dtype=str)
    ph_master["year"] = coerce_numeric(ph_master.get("year", pd.Series(dtype=str)))
    ph_master["dtp3_coverage"] = coerce_numeric(ph_master.get("dtp3_coverage", pd.Series(dtype=str)))
    ph_master["reported_cases"] = coerce_numeric(ph_master.get("reported_cases", pd.Series(dtype=str))).fillna(0)
    ph_master["genomes_per_case"] = coerce_numeric(ph_master.get("genomes_per_case", pd.Series(dtype=str)))
    ph_master["post_covid_period"] = coerce_numeric(ph_master.get("post_covid_period", pd.Series(dtype=str))).fillna(0)

    annual_program, country_summary = expand_program_metadata(program_metadata_path)
    merged = ph_master.merge(
        country_summary.drop(columns=["country_name_program_canonical"], errors="ignore"),
        on="country_iso3",
        how="left",
    )
    if not annual_program.empty:
        merged = merged.merge(annual_program, on=["country_iso3", "year"], how="left", suffixes=("", "_program"))

    merged["first_any_ap_year"] = merged["first_any_ap_year"].fillna(merged.get("first_any_ap_year_program"))
    merged["first_routine_ap_year"] = merged["first_routine_ap_year"].fillna(merged.get("first_routine_ap_year_program"))
    merged["booster_flag"] = coerce_numeric(
        merged.get("program_metadata_booster_flag", pd.Series(index=merged.index, dtype=float))
    ).fillna(0)
    merged["vaccine_program_type_effective"] = merged.get(
        "program_metadata_vaccine_program_type",
        pd.Series(index=merged.index, dtype=str),
    ).fillna(merged.get("vaccine_program_type", pd.Series(index=merged.index, dtype=str)))

    formulation_curation = load_formulation_curation(
        resolve_optional_input_path(formulation_curation_path, DEFAULT_FORMULATION_CURATION_PATH)
    )
    merged = apply_formulation_curation(merged, formulation_curation)
    product_metadata = load_product_metadata(
        resolve_optional_input_path(product_metadata_path, DEFAULT_PRODUCT_METADATA_PATH)
    )
    merged = apply_product_metadata(merged, product_metadata)

    merged["years_since_any_ap_intro"] = np.where(
        merged["first_any_ap_year"].notna(),
        np.maximum(0, merged["year"] - merged["first_any_ap_year"]),
        0,
    )
    merged["years_since_routine_ap_intro"] = np.where(
        merged["first_routine_ap_year"].notna(),
        np.maximum(0, merged["year"] - merged["first_routine_ap_year"]),
        0,
    )
    merged["years_since_ap_timing_anchor"] = np.where(
        merged["ap_timing_anchor_year_effective"].notna(),
        np.maximum(0, merged["year"] - merged["ap_timing_anchor_year_effective"]),
        0,
    )
    merged["ap_policy_phase_score"] = np.select(
        [
            merged["vaccine_program_type_effective"].eq("ap_introduced_routine_or_mixed"),
            merged["vaccine_program_type_effective"].eq("ap_targeted_or_partial_nonroutine"),
        ],
        [1.0, 0.5],
        default=0.0,
    )

    merged["dtp3_component_z"] = standardize_series(merged["dtp3_coverage"])
    merged["years_since_routine_ap_intro_z"] = standardize_series(merged["years_since_routine_ap_intro"])
    merged["years_since_ap_timing_anchor_z"] = standardize_series(merged["years_since_ap_timing_anchor"])
    merged["prn_formulation_component_z"] = standardize_series(merged["prn_formulation_component"])
    merged["booster_component"] = merged["booster_flag"].astype(float)
    merged["product_routine_primary_ap_prn_positive_share_z"] = standardize_series(
        merged["product_routine_primary_ap_prn_positive_share_conf_weighted"]
    )
    merged["product_routine_booster_ap_prn_positive_share_z"] = standardize_series(
        merged["product_routine_booster_ap_prn_positive_share_conf_weighted"]
    )
    merged["ap_exposure_v3_primary_component_available"] = (
        merged["product_routine_primary_observed_share"].notna()
    )
    merged["ap_exposure_v3_booster_component_required"] = merged["booster_flag"].fillna(0).astype(float).gt(0)
    merged["ap_exposure_v3_booster_component_available"] = (
        merged["product_routine_booster_observed_share"].notna()
    )
    merged["ap_exposure_v3_available"] = (
        merged["ap_exposure_v3_primary_component_available"]
        & (
            ~merged["ap_exposure_v3_booster_component_required"]
            | merged["ap_exposure_v3_booster_component_available"]
        )
    )
    merged["ap_exposure_v3_component_status"] = np.select(
        [
            ~merged["ap_exposure_v3_primary_component_available"]
            & merged["ap_exposure_v3_booster_component_required"]
            & ~merged["ap_exposure_v3_booster_component_available"],
            ~merged["ap_exposure_v3_primary_component_available"],
            merged["ap_exposure_v3_booster_component_required"]
            & ~merged["ap_exposure_v3_booster_component_available"],
        ],
        [
            "missing_primary_and_booster_role_specific_product_components",
            "missing_primary_role_specific_product_component",
            "missing_booster_role_specific_product_component",
        ],
        default="complete_role_specific_product_components",
    )

    vaccine_variable_decision, vaccine_variable_recommendation = load_vaccine_variable_coverage_decision(ph_master_path)
    primary_lambda = select_primary_value(lambda_range)
    primary_gamma = select_primary_value(gamma_range)
    delta_prn = 1.0

    records: list[pd.DataFrame] = []
    for lambda_value, gamma_value in product(lambda_range, gamma_range):
        parameter_frame = merged.copy()
        parameter_frame["exposure_version"] = exposure_version
        parameter_frame["exposure_lambda_years"] = float(lambda_value)
        parameter_frame["exposure_gamma_booster"] = float(gamma_value)
        parameter_frame["exposure_delta_prn"] = float(delta_prn)
        parameter_frame["exposure_formula_id"] = f"{exposure_version}_lambda_{lambda_value:g}_gamma_{gamma_value:g}"
        parameter_frame["is_primary_parameterization"] = (
            (float(lambda_value) == primary_lambda) and (float(gamma_value) == primary_gamma)
        )
        parameter_frame["ap_exposure_v1_score"] = (
            parameter_frame["dtp3_component_z"]
            + float(lambda_value) * parameter_frame["years_since_routine_ap_intro_z"]
            + float(gamma_value) * parameter_frame["booster_component"]
        )
        parameter_frame["ap_exposure_v2_score"] = (
            parameter_frame["dtp3_component_z"]
            + float(lambda_value) * parameter_frame["years_since_ap_timing_anchor_z"]
            + float(gamma_value) * parameter_frame["booster_component"]
            + float(delta_prn) * parameter_frame["prn_formulation_component_z"]
        )
        parameter_frame["ap_exposure_v3_score"] = (
            parameter_frame["dtp3_component_z"]
            + float(lambda_value) * parameter_frame["years_since_ap_timing_anchor_z"]
            + float(delta_prn) * parameter_frame["product_routine_primary_ap_prn_positive_share_z"]
            + float(gamma_value) * parameter_frame["product_routine_booster_ap_prn_positive_share_z"]
        )
        parameter_frame["dtp3_only_score"] = parameter_frame["dtp3_coverage"]
        parameter_frame["exposure_version_effective"] = np.where(
            parameter_frame["ap_exposure_v2_available"],
            "v2_curated",
            "v1_fallback",
        )
        parameter_frame["vaccine_variable_decision"] = vaccine_variable_decision
        parameter_frame["vaccine_variable_recommendation"] = vaccine_variable_recommendation
        parameter_frame["exposure_score_interpretation"] = (
            "heuristic_global_z_score_composite_not_absolute_biologic_scale"
        )
        parameter_frame["exposure_notes"] = (
            "V1 uses DTP3 plus years since routine/mixed aP intro plus booster_flag; "
            "V2 replaces WHO-only timing with curated formulation-aware timing and adds PRN-in-formulation status; "
            "V3 swaps the coarse PRN term for role-specific product metadata on routine primary and booster exposure; "
            "all score versions are heuristic global z-score composites intended for relative ranking, not absolute biologic effect scales."
        )
        records.append(parameter_frame)

    output_frame = pd.concat(records, ignore_index=True)
    output_path = Path(output_index_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_frame.to_csv(output_path, sep="\t", index=False)

    primary_rows = output_frame.loc[output_frame["is_primary_parameterization"]].copy()
    formulation_country_year_path, formulation_summary_path = write_formulation_outputs(primary_rows, output_index_path)
    product_country_year_path, product_summary_path = write_product_metadata_outputs(primary_rows, output_index_path)
    precedence_audit_path = write_precedence_audit(primary_rows, output_index_path)

    figure_path = Path(output_figure_path)
    figure_path.parent.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({"figure.figsize": (12.0, 4.8)})
    figure, axes = plt.subplots(1, 3)

    program_summary = (
        primary_rows.groupby(["year", "vaccine_program_type_effective"], dropna=False)[["ap_exposure_v1_score", "ap_exposure_v2_score"]]
        .median()
        .reset_index()
    )
    for program_type, group in program_summary.groupby("vaccine_program_type_effective", dropna=False):
        label = program_type if isinstance(program_type, str) and program_type else "missing_program_type"
        axes[0].plot(group["year"], group["ap_exposure_v1_score"], label=label)
    axes[0].set_title("Primary V1 score by program phase")
    axes[0].set_xlabel("Year")
    axes[0].set_ylabel("Median exposure score")
    axes[0].legend(fontsize=6, loc="best")

    v2_summary = (
        primary_rows.groupby(["year", "prn_in_vaccine_curated"], dropna=False)["ap_exposure_v2_score"]
        .median()
        .reset_index()
    )
    for prn_value, group in v2_summary.groupby("prn_in_vaccine_curated", dropna=False):
        label = prn_value if isinstance(prn_value, str) and prn_value else "unknown"
        axes[1].plot(group["year"], group["ap_exposure_v2_score"], label=label)
    axes[1].set_title("Primary V2 score by PRN status")
    axes[1].set_xlabel("Year")
    axes[1].set_ylabel("Median exposure score")
    axes[1].legend(fontsize=6, loc="best")

    heatmap_rows = []
    for lambda_value, gamma_value in product(sorted(lambda_range), sorted(gamma_range)):
        alternative = output_frame.loc[
            output_frame["exposure_formula_id"].eq(f"{exposure_version}_lambda_{lambda_value:g}_gamma_{gamma_value:g}")
        ]
        joined = primary_rows[["country_iso3", "year", "ap_exposure_v2_score"]].merge(
            alternative[["country_iso3", "year", "ap_exposure_v2_score"]],
            on=["country_iso3", "year"],
            suffixes=("_primary", "_alt"),
        )
        correlation = joined[["ap_exposure_v2_score_primary", "ap_exposure_v2_score_alt"]].corr().iloc[0, 1]
        heatmap_rows.append({"lambda": lambda_value, "gamma": gamma_value, "correlation": correlation})
    heatmap = pd.DataFrame(heatmap_rows).pivot(index="gamma", columns="lambda", values="correlation")
    image = axes[2].imshow(heatmap.values, aspect="auto", origin="lower", vmin=0, vmax=1, cmap="viridis")
    axes[2].set_xticks(range(len(heatmap.columns)), [f"{value:g}" for value in heatmap.columns])
    axes[2].set_yticks(range(len(heatmap.index)), [f"{value:g}" for value in heatmap.index])
    axes[2].set_xlabel("lambda years")
    axes[2].set_ylabel("gamma booster")
    axes[2].set_title("Correlation with primary V2 score")
    figure.colorbar(image, ax=axes[2], fraction=0.046, pad=0.04)
    figure.tight_layout()
    figure.savefig(figure_path)
    plt.close(figure)

    print(f"Wrote exposure index: {output_path}")
    print(f"Wrote formulation country-year table: {formulation_country_year_path}")
    print(f"Wrote formulation summary table: {formulation_summary_path}")
    print(f"Wrote product metadata country-year table: {product_country_year_path}")
    print(f"Wrote product metadata summary table: {product_summary_path}")
    print(f"Wrote formulation precedence audit: {precedence_audit_path}")
    return output_frame


if "snakemake" in globals():
    build_index(
        ph_master_path=snakemake.input.ph_master,
        program_metadata_path=snakemake.input.program_metadata,
        formulation_curation_path=snakemake.input.get("formulation_curation", ""),
        product_metadata_path=snakemake.input.get("product_metadata", ""),
        output_index_path=snakemake.output.index,
        output_figure_path=snakemake.output.sensitivity,
        exposure_version=snakemake.params.version,
        lambda_range=[float(value) for value in snakemake.params.lambda_range],
        gamma_range=[float(value) for value in snakemake.params.gamma_range],
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Build workflow exposure indices")
    parser.add_argument("--ph-master", required=True, help="Public health country-year master TSV")
    parser.add_argument("--program-metadata", required=True, help="Country program metadata TSV")
    parser.add_argument(
        "--formulation-curation",
        default="modules/public_health/inputs/curation/vaccine_formulation_curation.tsv",
        help="Detailed vaccine formulation curation TSV",
    )
    parser.add_argument(
        "--product-metadata",
        default="modules/public_health/inputs/curation/vaccine_product_metadata.tsv",
        help="Role-specific vaccine product metadata TSV",
    )
    parser.add_argument("--index-out", required=True, help="Output exposure index TSV")
    parser.add_argument("--sensitivity-out", required=True, help="Output sensitivity PDF")
    parser.add_argument("--version", default="v1", help="Exposure index version label")
    parser.add_argument("--lambda-range", nargs="+", type=float, default=[0.5, 1.0, 2.0])
    parser.add_argument("--gamma-range", nargs="+", type=float, default=[0.0, 0.5, 1.0])
    arguments = parser.parse_args()

    build_index(
        ph_master_path=arguments.ph_master,
        program_metadata_path=arguments.program_metadata,
        formulation_curation_path=arguments.formulation_curation,
        product_metadata_path=arguments.product_metadata,
        output_index_path=arguments.index_out,
        output_figure_path=arguments.sensitivity_out,
        exposure_version=arguments.version,
        lambda_range=arguments.lambda_range,
        gamma_range=arguments.gamma_range,
    )
