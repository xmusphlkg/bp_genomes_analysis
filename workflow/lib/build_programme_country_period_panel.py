#!/usr/bin/env python3
"""Build programme-surveillance country-period panels from the annual ecology layer."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from workflow.lib.project_paths import project_workflow_root


ELIGIBLE_PRIMARY_CLASSES = {
    "wp_only_or_pre_ap",
    "routine_ap_prn_positive",
    "routine_ap_prn_negative",
    "routine_ap_mixed",
}
AP_SUPPORTING_CLASSES = {
    "routine_ap_prn_positive",
    "routine_ap_prn_negative",
    "routine_ap_mixed",
    "routine_ap_unknown",
}
CONFIDENCE_RANK = {"high": 3, "medium": 2, "low": 1, "unknown": 0, "": 0}
RANK_TO_CONFIDENCE = {value: key for key, value in CONFIDENCE_RANK.items()}
ANALYSIS_PANEL_PREFIX = "programme_country_period"
DEFAULT_ORIGIN_DESCENDANTS_DIR = project_workflow_root() / "asr" / "event_subtrees"


def normalize_text(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def coerce_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def parse_bool(value: object) -> bool:
    text = normalize_text(value).lower()
    return text in {"true", "1", "yes", "y", "t"}


def mode_or_empty(series: pd.Series) -> str:
    clean = series.dropna().astype(str)
    clean = clean[clean.str.len() > 0]
    if clean.empty:
        return ""
    return clean.mode().iloc[0]


def dominant_with_latest_tiebreak(series: pd.Series, years: pd.Series | None = None) -> str:
    clean = series.fillna("").astype(str).str.strip()
    clean = clean[clean.str.len() > 0]
    if clean.empty:
        return ""

    counts = clean.value_counts()
    top_count = int(counts.max())
    candidates = set(counts[counts.eq(top_count)].index.tolist())
    if len(candidates) == 1:
        return next(iter(candidates))

    if years is None:
        return sorted(candidates)[0]

    year_frame = pd.DataFrame({"value": series, "year": years})
    year_frame["value"] = year_frame["value"].fillna("").astype(str).str.strip()
    year_frame = year_frame.loc[year_frame["value"].isin(candidates)].copy()
    year_frame["year"] = coerce_numeric(year_frame["year"])
    year_frame = year_frame.sort_values(["year", "value"], ascending=[True, True])
    if year_frame.empty:
        return sorted(candidates)[0]
    return normalize_text(year_frame.iloc[-1]["value"])


def conservative_confidence(series: pd.Series) -> str:
    values = [CONFIDENCE_RANK.get(normalize_text(value).lower(), 0) for value in series]
    if not values:
        return ""
    return RANK_TO_CONFIDENCE[min(values)]


def period_class_for_group(series: pd.Series, years: pd.Series | None = None) -> str:
    clean = series.fillna("").astype(str).str.strip()
    clean = clean[(clean.str.len() > 0) & clean.ne("unknown")]
    if clean.empty:
        return "unknown"
    unique = sorted(set(clean.tolist()))
    if len(unique) == 1:
        return unique[0]
    if set(unique).issubset(AP_SUPPORTING_CLASSES):
        return "routine_ap_mixed"
    if years is None:
        return dominant_with_latest_tiebreak(clean)
    return dominant_with_latest_tiebreak(clean, years.loc[clean.index])


def transition_details(frame: pd.DataFrame) -> tuple[bool, str, str, str, int]:
    clean = frame["program_formulation_class"].fillna("").astype(str).str.strip()
    clean = clean[(clean.str.len() > 0) & clean.ne("unknown")]
    unique_classes = sorted(set(clean.tolist()))
    if not unique_classes:
        return False, "", "unknown", "", 0

    dominant_class = dominant_with_latest_tiebreak(clean, frame.loc[clean.index, "year"])
    transition_flag = len(unique_classes) > 1
    if not transition_flag:
        return False, "", dominant_class, unique_classes[0], 1

    includes_wp = "wp_only_or_pre_ap" in unique_classes
    includes_ap = any(value in AP_SUPPORTING_CLASSES for value in unique_classes)
    if includes_wp and includes_ap:
        transition_type = "ap_vs_wp_transition"
    else:
        transition_type = "within_ap_transition"
    return True, transition_type, dominant_class, ",".join(unique_classes), len(unique_classes)


def numeric_min_or_nan(series: pd.Series) -> float:
    valid = coerce_numeric(series).dropna()
    if valid.empty:
        return np.nan
    return float(valid.min())


def numeric_max_or_nan(series: pd.Series) -> float:
    valid = coerce_numeric(series).dropna()
    if valid.empty:
        return np.nan
    return float(valid.max())


def build_origin_bridge(prevalence: pd.DataFrame, origin_descendants_dir: str | None) -> pd.DataFrame:
    annual = prevalence[["country_iso3", "year"]].copy()
    annual["country_iso3"] = annual["country_iso3"].fillna("").astype(str).str.upper()
    annual["year"] = coerce_numeric(annual["year"])
    annual = annual.loc[annual["country_iso3"].ne("") & annual["year"].notna()].copy()
    annual["year"] = annual["year"].astype(int)
    annual = annual.drop_duplicates(subset=["country_iso3", "year"]).copy()

    annual["n_origin_clades_active"] = 0
    annual["n_new_origins_detected"] = 0
    annual["first_local_origin_year"] = np.nan

    if not origin_descendants_dir:
        descendants = pd.DataFrame()
    else:
        rows: list[pd.DataFrame] = []
        subtree_dir = Path(origin_descendants_dir)
        for path in sorted(subtree_dir.glob("origin_*.descendant_tips.tsv")):
            frame = pd.read_csv(path, sep="\t", dtype=str)
            frame["origin_id"] = path.stem.split(".")[0]
            rows.append(frame)
        descendants = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()

    if not descendants.empty:
        descendants["country_iso3"] = descendants.get("country_iso3", pd.Series(dtype=str)).fillna("").astype(str).str.upper()
        descendants["year"] = coerce_numeric(descendants.get("year", pd.Series(dtype=str)))
        descendants = descendants.loc[
            descendants.get("observed_prn_state", pd.Series(dtype=str)).fillna("").eq("disrupted")
            & descendants["country_iso3"].ne("")
            & descendants["year"].notna()
        ].copy()
        descendants["year"] = descendants["year"].astype(int)

        active = (
            descendants.groupby(["country_iso3", "year"], dropna=False)["origin_id"]
            .nunique()
            .rename("n_origin_clades_active")
            .reset_index()
        )
        annual = annual.merge(active, on=["country_iso3", "year"], how="left", suffixes=("", "_bridge"))
        annual["n_origin_clades_active"] = (
            coerce_numeric(annual.get("n_origin_clades_active_bridge", pd.Series(index=annual.index, dtype=float)))
            .fillna(annual["n_origin_clades_active"])
            .fillna(0)
            .astype(int)
        )
        annual = annual.drop(columns=[column for column in ["n_origin_clades_active_bridge"] if column in annual.columns])

        local_first = (
            descendants.groupby(["country_iso3", "origin_id"], dropna=False)["year"]
            .min()
            .rename("first_local_origin_year")
            .reset_index()
        )
        new_origins = (
            local_first.groupby(["country_iso3", "first_local_origin_year"], dropna=False)
            .size()
            .rename("n_new_origins_detected")
            .reset_index()
            .rename(columns={"first_local_origin_year": "year"})
        )
        annual = annual.merge(new_origins, on=["country_iso3", "year"], how="left", suffixes=("", "_bridge"))
        annual["n_new_origins_detected"] = (
            coerce_numeric(annual.get("n_new_origins_detected_bridge", pd.Series(index=annual.index, dtype=float)))
            .fillna(annual["n_new_origins_detected"])
            .fillna(0)
            .astype(int)
        )
        annual = annual.drop(columns=[column for column in ["n_new_origins_detected_bridge"] if column in annual.columns])

        country_first_local = (
            local_first.groupby("country_iso3", dropna=False)["first_local_origin_year"]
            .min()
            .rename("first_local_origin_year")
            .reset_index()
        )
        annual = annual.drop(columns=["first_local_origin_year"], errors="ignore").merge(
            country_first_local,
            on="country_iso3",
            how="left",
        )

    detection_rows = prevalence.copy()
    detection_rows["country_iso3"] = detection_rows.get("country_iso3", pd.Series(dtype=str)).fillna("").astype(str).str.upper()
    detection_rows["year"] = coerce_numeric(detection_rows.get("year", pd.Series(dtype=str)))
    detection_rows["n_prn_disrupted"] = coerce_numeric(detection_rows.get("n_prn_disrupted", pd.Series(dtype=str))).fillna(0)
    detection_rows = detection_rows.loc[
        detection_rows["country_iso3"].ne("") & detection_rows["year"].notna() & detection_rows["n_prn_disrupted"].gt(0)
    ].copy()
    detection_rows["year"] = detection_rows["year"].astype(int)
    first_detection = (
        detection_rows.groupby("country_iso3", dropna=False)["year"]
        .min()
        .rename("first_prn_detection_year")
        .reset_index()
    )

    annual = annual.merge(first_detection, on="country_iso3", how="left")
    return annual[
        [
            "country_iso3",
            "year",
            "n_origin_clades_active",
            "n_new_origins_detected",
            "first_local_origin_year",
            "first_prn_detection_year",
        ]
    ]


def exclusion_reason(
    row: pd.Series,
    *,
    min_interpretable: int,
) -> str:
    if pd.isna(row["response_ipw_prevalence"]) or float(row["response_ipw_weight_total"] or 0) <= 0:
        return "missing_ipw_response"
    if pd.isna(row["reported_cases_period"]) or float(row["reported_cases_period"] or 0) <= 0:
        return "missing_reported_cases"
    if pd.isna(row["genomes_per_case_effective"]) or float(row["genomes_per_case_effective"] or 0) <= 0:
        return "missing_genomes_per_case"
    if int(row["response_n_genomes_prn_interpretable"] or 0) < min_interpretable:
        return f"below_interpretable_threshold_{min_interpretable}"
    if normalize_text(row["program_formulation_class_concurrent"]) not in ELIGIBLE_PRIMARY_CLASSES:
        return f"class_excluded_{normalize_text(row['program_formulation_class_concurrent']) or 'unknown'}"
    return ""


def prepare_annual_dataset_from_frames(
    exposure: pd.DataFrame,
    prevalence: pd.DataFrame,
    origin_descendants_dir: str | None = None,
) -> pd.DataFrame:
    exposure = exposure.copy()
    prevalence = prevalence.copy()

    for frame in [exposure, prevalence]:
        frame["country_iso3"] = frame.get("country_iso3", pd.Series(dtype=str)).fillna("").astype(str).str.upper()
        frame["year"] = coerce_numeric(frame.get("year", pd.Series(dtype=str)))

    exposure["reported_cases"] = coerce_numeric(exposure.get("reported_cases", pd.Series(dtype=str)))
    exposure["genomes_per_case"] = coerce_numeric(exposure.get("genomes_per_case", pd.Series(dtype=str)))
    exposure["post_covid_period"] = coerce_numeric(exposure.get("post_covid_period", pd.Series(dtype=str))).fillna(0)
    exposure["program_formulation_conflict"] = exposure.get(
        "program_formulation_conflict",
        pd.Series(index=exposure.index, dtype=str),
    ).map(parse_bool)
    exposure["years_since_any_ap_use"] = coerce_numeric(
        exposure.get("years_since_any_ap_use", exposure.get("years_since_any_ap_intro", pd.Series(dtype=str)))
    )

    prevalence["n_genomes_total"] = coerce_numeric(prevalence.get("n_genomes_total", pd.Series(dtype=str))).fillna(0)
    prevalence["n_genomes_prn_interpretable"] = coerce_numeric(
        prevalence.get("n_genomes_prn_interpretable", pd.Series(dtype=str))
    ).fillna(0)
    prevalence["n_prn_disrupted"] = coerce_numeric(prevalence.get("n_prn_disrupted", pd.Series(dtype=str))).fillna(0)
    prevalence["ipw_weight_total"] = coerce_numeric(prevalence.get("ipw_weight_total", pd.Series(dtype=str)))
    prevalence["ipw_prevalence"] = coerce_numeric(prevalence.get("ipw_prevalence", pd.Series(dtype=str)))
    prevalence["naive_prevalence"] = coerce_numeric(prevalence.get("naive_prevalence", pd.Series(dtype=str)))
    prevalence["boundary_lower_prevalence"] = coerce_numeric(
        prevalence.get("boundary_lower_prevalence", pd.Series(dtype=str))
    )
    prevalence["boundary_upper_prevalence"] = coerce_numeric(
        prevalence.get("boundary_upper_prevalence", pd.Series(dtype=str))
    )
    prevalence["n_missing_outcomes"] = coerce_numeric(prevalence.get("n_missing_outcomes", pd.Series(dtype=str)))

    exposure = exposure.loc[
        exposure["country_iso3"].ne("")
        & exposure["year"].notna()
        & exposure.get("is_primary_parameterization", pd.Series(index=exposure.index, dtype=str))
        .fillna("")
        .astype(str)
        .str.lower()
        .eq("true")
    ].copy()
    exposure["year"] = exposure["year"].astype(int)

    prevalence_subset = prevalence[
        [
            "country_iso3",
            "year",
            "n_genomes_total",
            "n_genomes_prn_interpretable",
            "n_prn_disrupted",
            "ipw_weight_total",
            "ipw_prevalence",
            "naive_prevalence",
            "boundary_lower_prevalence",
            "boundary_upper_prevalence",
            "n_missing_outcomes",
        ]
    ].rename(
        columns={
            "n_genomes_total": "response_n_genomes_total",
            "n_genomes_prn_interpretable": "response_n_genomes_prn_interpretable",
            "n_prn_disrupted": "response_n_prn_disrupted",
            "ipw_weight_total": "response_ipw_weight_total",
            "ipw_prevalence": "response_ipw_prevalence",
            "naive_prevalence": "response_naive_prevalence",
            "boundary_lower_prevalence": "response_boundary_lower_prevalence",
            "boundary_upper_prevalence": "response_boundary_upper_prevalence",
            "n_missing_outcomes": "response_n_missing_outcomes",
        }
    )

    dataset = exposure.merge(prevalence_subset, on=["country_iso3", "year"], how="left")
    for column in [
        "response_n_genomes_total",
        "response_n_genomes_prn_interpretable",
        "response_n_prn_disrupted",
    ]:
        dataset[column] = coerce_numeric(dataset.get(column, pd.Series(index=dataset.index, dtype=str))).fillna(0)
    for column in [
        "response_ipw_weight_total",
        "response_ipw_prevalence",
        "response_naive_prevalence",
        "response_boundary_lower_prevalence",
        "response_boundary_upper_prevalence",
        "response_n_missing_outcomes",
    ]:
        dataset[column] = coerce_numeric(dataset.get(column, pd.Series(index=dataset.index, dtype=str)))
    dataset["response_n_missing_outcomes"] = dataset["response_n_missing_outcomes"].fillna(
        dataset["response_n_genomes_total"] - dataset["response_n_genomes_prn_interpretable"]
    )
    boundary_lower = pd.Series(
        np.where(
            dataset["response_n_genomes_total"].gt(0),
            dataset["response_n_prn_disrupted"] / dataset["response_n_genomes_total"],
            np.nan,
        ),
        index=dataset.index,
    )
    boundary_upper = pd.Series(
        np.where(
            dataset["response_n_genomes_total"].gt(0),
            (
                dataset["response_n_prn_disrupted"]
                + dataset["response_n_missing_outcomes"]
            ) / dataset["response_n_genomes_total"],
            np.nan,
        ),
        index=dataset.index,
    )
    dataset["response_boundary_lower_prevalence"] = dataset["response_boundary_lower_prevalence"].fillna(boundary_lower)
    dataset["response_boundary_upper_prevalence"] = dataset["response_boundary_upper_prevalence"].fillna(boundary_upper)
    dataset["workflow_genomes_per_case"] = np.where(
        dataset["reported_cases"].fillna(0).gt(0),
        dataset["response_n_genomes_total"] / dataset["reported_cases"],
        np.nan,
    )
    dataset["genomes_per_case_effective"] = dataset["genomes_per_case"].fillna(dataset["workflow_genomes_per_case"])

    origin_bridge = build_origin_bridge(prevalence, origin_descendants_dir)
    dataset = dataset.merge(origin_bridge, on=["country_iso3", "year"], how="left")
    dataset["n_origin_clades_active"] = coerce_numeric(
        dataset.get("n_origin_clades_active", pd.Series(index=dataset.index, dtype=str))
    ).fillna(0)
    dataset["n_new_origins_detected"] = coerce_numeric(
        dataset.get("n_new_origins_detected", pd.Series(index=dataset.index, dtype=str))
    ).fillna(0)
    dataset["first_local_origin_year"] = coerce_numeric(
        dataset.get("first_local_origin_year", pd.Series(index=dataset.index, dtype=str))
    )
    dataset["first_prn_detection_year"] = coerce_numeric(
        dataset.get("first_prn_detection_year", pd.Series(index=dataset.index, dtype=str))
    )

    first_prn_positive_year = (
        dataset.loc[dataset["program_formulation_class"].eq("routine_ap_prn_positive"), ["country_iso3", "year"]]
        .groupby("country_iso3", dropna=False)["year"]
        .min()
        .rename("first_prn_positive_routine_use_year")
        .reset_index()
    )
    dataset = dataset.merge(first_prn_positive_year, on="country_iso3", how="left")
    dataset["first_prn_positive_routine_use_year"] = coerce_numeric(
        dataset.get("first_prn_positive_routine_use_year", pd.Series(index=dataset.index, dtype=str))
    )
    dataset["years_since_prn_positive_routine_use"] = np.where(
        dataset["first_prn_positive_routine_use_year"].notna()
        & dataset["year"].ge(dataset["first_prn_positive_routine_use_year"]),
        dataset["year"] - dataset["first_prn_positive_routine_use_year"],
        np.nan,
    )
    dataset["years_since_first_local_origin"] = np.where(
        dataset["first_local_origin_year"].notna() & dataset["year"].ge(dataset["first_local_origin_year"]),
        dataset["year"] - dataset["first_local_origin_year"],
        np.nan,
    )
    dataset["years_since_first_prn_detection"] = np.where(
        dataset["first_prn_detection_year"].notna() & dataset["year"].ge(dataset["first_prn_detection_year"]),
        dataset["year"] - dataset["first_prn_detection_year"],
        np.nan,
    )
    return dataset


def prepare_annual_dataset(
    exposure_path: str,
    prevalence_path: str,
    origin_descendants_dir: str | None = None,
) -> pd.DataFrame:
    exposure = pd.read_csv(exposure_path, sep="\t", dtype=str)
    prevalence = pd.read_csv(prevalence_path, sep="\t", dtype=str)
    return prepare_annual_dataset_from_frames(
        exposure,
        prevalence,
        origin_descendants_dir=origin_descendants_dir,
    )


def build_period_panel(
    annual: pd.DataFrame,
    *,
    bin_size: int,
    min_interpretable: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    annual = annual.copy()
    if "response_n_missing_outcomes" not in annual.columns:
        annual["response_n_missing_outcomes"] = (
            annual.get("response_n_genomes_total", pd.Series(index=annual.index, dtype=float)).fillna(0)
            - annual.get("response_n_genomes_prn_interpretable", pd.Series(index=annual.index, dtype=float)).fillna(0)
        )
    for column in [
        "n_new_origins_detected",
        "n_origin_clades_active",
        "first_local_origin_year",
        "first_prn_detection_year",
        "years_since_prn_positive_routine_use",
        "years_since_any_ap_use",
    ]:
        if column not in annual.columns:
            annual[column] = np.nan if "year" in column else 0
    annual["period_start"] = (annual["year"] // bin_size) * bin_size
    annual["period_end"] = annual["period_start"] + (bin_size - 1)
    annual["weighted_ipw_successes"] = annual["response_ipw_prevalence"] * annual["response_ipw_weight_total"]

    for (country_iso3, period_start), frame in annual.groupby(["country_iso3", "period_start"], dropna=False):
        period_end = int(frame["period_end"].iloc[0])
        total_weight = float(frame["response_ipw_weight_total"].fillna(0).sum())
        weighted_ipw_successes = float(frame["weighted_ipw_successes"].fillna(0).sum())
        reported_cases_period = float(frame["reported_cases"].fillna(0).sum())
        total_genomes = float(frame["response_n_genomes_total"].fillna(0).sum())
        total_interpretable = float(frame["response_n_genomes_prn_interpretable"].fillna(0).sum())
        total_disrupted = float(frame["response_n_prn_disrupted"].fillna(0).sum())
        total_missing = float(frame["response_n_missing_outcomes"].fillna(0).sum())

        transition_flag, transition_type, concurrent_class, unique_classes, n_unique_classes = transition_details(frame)
        notes: list[str] = []
        if bool(frame["program_formulation_conflict"].fillna(False).any()):
            notes.append("contains_curated_vs_program_phase_conflict_resolved_by_curated_precedence")
        if transition_flag:
            notes.append(f"transition_period_flag={transition_type}")
        if frame["response_n_genomes_total"].fillna(0).eq(0).any() and reported_cases_period > 0:
            notes.append("reported_cases_aggregated_across_full_period_including_zero_genome_years")

        first_local_origin_year = numeric_min_or_nan(frame["first_local_origin_year"])
        first_prn_detection_year = numeric_min_or_nan(frame["first_prn_detection_year"])
        years_since_first_local_origin = (
            float(period_end - first_local_origin_year)
            if not np.isnan(first_local_origin_year) and period_end >= first_local_origin_year
            else np.nan
        )
        years_since_first_prn_detection = (
            float(period_end - first_prn_detection_year)
            if not np.isnan(first_prn_detection_year) and period_end >= first_prn_detection_year
            else np.nan
        )

        row = {
            "analysis_panel": f"{ANALYSIS_PANEL_PREFIX}_{bin_size}y",
            "country_iso3": country_iso3,
            "country_name": dominant_with_latest_tiebreak(frame["country_name"], frame["year"]),
            "period_start": int(period_start),
            "period_end": int(period_end),
            "period_label": f"{int(period_start)}-{int(period_end)}",
            "years_observed": ",".join(str(value) for value in sorted(frame["year"].astype(int).tolist())),
            "n_years_observed": int(frame["year"].nunique()),
            "n_annual_rows": int(len(frame)),
            "reported_cases_period": reported_cases_period,
            "response_n_genomes_total": int(total_genomes),
            "response_n_genomes_prn_interpretable": int(total_interpretable),
            "response_n_prn_disrupted": int(total_disrupted),
            "response_n_missing_outcomes": int(total_missing),
            "response_ipw_weight_total": total_weight,
            "response_ipw_successes_est": weighted_ipw_successes,
            "response_ipw_prevalence": np.nan if total_weight <= 0 else weighted_ipw_successes / total_weight,
            "response_naive_prevalence": np.nan if total_interpretable <= 0 else total_disrupted / total_interpretable,
            "response_boundary_lower_prevalence": (
                np.nan if total_genomes <= 0 else total_disrupted / total_genomes
            ),
            "response_boundary_upper_prevalence": (
                np.nan if total_genomes <= 0 else (total_disrupted + total_missing) / total_genomes
            ),
            "genomes_per_case_effective": np.nan if reported_cases_period <= 0 else total_genomes / reported_cases_period,
            "post_covid_period": int(frame["post_covid_period"].fillna(0).max()),
            "program_formulation_class_concurrent": concurrent_class,
            "program_formulation_class": concurrent_class,
            "program_formulation_class_lagged": "",
            "lagged_class_available": False,
            "lagged_from_period_label": "",
            "transition_period_flag": bool(transition_flag),
            "transition_period_type": transition_type,
            "unique_program_classes_within_period": unique_classes,
            "n_unique_program_classes_within_period": int(n_unique_classes),
            "share_years_prn_positive_within_period": float(
                frame["program_formulation_class"].fillna("").astype(str).eq("routine_ap_prn_positive").mean()
            ),
            "years_since_prn_positive_routine_use": numeric_max_or_nan(frame["years_since_prn_positive_routine_use"]),
            "years_since_any_ap_use": numeric_max_or_nan(frame["years_since_any_ap_use"]),
            "n_new_origins_detected_period": int(frame["n_new_origins_detected"].fillna(0).sum()),
            "n_active_origin_clades_period": int(frame["n_origin_clades_active"].fillna(0).max()),
            "first_local_origin_year": first_local_origin_year,
            "first_prn_detection_year": first_prn_detection_year,
            "years_since_first_local_origin": years_since_first_local_origin,
            "years_since_first_prn_detection": years_since_first_prn_detection,
            "has_local_origin_by_period_end": bool(not np.isnan(first_local_origin_year)),
            "has_prn_detection_by_period_end": bool(not np.isnan(first_prn_detection_year)),
            "dominant_prn_in_vaccine_curated": dominant_with_latest_tiebreak(frame["prn_in_vaccine_curated"], frame["year"]),
            "dominant_primary_series_formulation": dominant_with_latest_tiebreak(frame["primary_series_formulation"], frame["year"]),
            "dominant_booster_formulation": dominant_with_latest_tiebreak(frame["booster_formulation"], frame["year"]),
            "formulation_confidence_period": conservative_confidence(frame["formulation_confidence"]),
            "period_contains_conflict": bool(frame["program_formulation_conflict"].fillna(False).any()),
            "dominant_exposure_precedence_rule": dominant_with_latest_tiebreak(frame["exposure_precedence_rule"], frame["year"]),
            "country_row_share": np.nan,
            "primary_panel_eligible": False,
            "exclusion_reason": "",
            "notes": ";".join(notes),
        }
        row["exclusion_reason"] = exclusion_reason(pd.Series(row), min_interpretable=min_interpretable)
        row["primary_panel_eligible"] = row["exclusion_reason"] == ""
        rows.append(row)

    panel = pd.DataFrame.from_records(rows)
    if panel.empty:
        return panel

    panel = panel.sort_values(["country_iso3", "period_start"]).reset_index(drop=True)
    panel["program_formulation_class_lagged"] = (
        panel.groupby("country_iso3", dropna=False)["program_formulation_class_concurrent"].shift(1).fillna("")
    )
    panel["lagged_class_available"] = panel["program_formulation_class_lagged"].astype(str).str.len().gt(0)
    panel["lagged_from_period_label"] = (
        panel.groupby("country_iso3", dropna=False)["period_label"].shift(1).fillna("")
    )

    eligible = panel.loc[panel["primary_panel_eligible"]].copy()
    n_eligible = len(eligible)
    share_lookup = {}
    if n_eligible:
        share_lookup = eligible.groupby("country_iso3", dropna=False).size().div(n_eligible).to_dict()
    panel["country_row_share"] = panel["country_iso3"].map(share_lookup).fillna(0.0)
    return panel.sort_values(["primary_panel_eligible", "country_iso3", "period_start"], ascending=[False, True, True])


def write_coverage_report(panel: pd.DataFrame, output_path: Path, *, min_interpretable: int, bin_size: int) -> Path:
    coverage_path = output_path.with_name(f"{output_path.stem}_coverage.tsv")
    eligible = panel.loc[panel["primary_panel_eligible"]].copy()
    coverage = pd.DataFrame(
        [
            {
                "metric": "analysis_panel",
                "value": f"{ANALYSIS_PANEL_PREFIX}_{bin_size}y",
                "notes": "primary programme-surveillance period panel label",
            },
            {
                "metric": "period_bin_size_years",
                "value": bin_size,
                "notes": "years per country-period bin",
            },
            {
                "metric": "min_interpretable_genomes_per_period",
                "value": min_interpretable,
                "notes": "eligibility threshold",
            },
            {
                "metric": "eligible_rows",
                "value": int(len(eligible)),
                "notes": "period rows eligible for the primary programme-surveillance model",
            },
            {
                "metric": "eligible_countries",
                "value": int(eligible["country_iso3"].nunique()) if not eligible.empty else 0,
                "notes": "unique countries contributing eligible period rows",
            },
            {
                "metric": "max_country_row_share",
                "value": "" if eligible.empty else f"{eligible['country_row_share'].max():.6f}",
                "notes": "largest eligible-row share contributed by a single country",
            },
            {
                "metric": "eligible_rows_high_confidence",
                "value": int(eligible["formulation_confidence_period"].astype(str).str.lower().eq("high").sum()),
                "notes": "eligible rows with high-confidence formulation curation",
            },
            {
                "metric": "eligible_rows_transition_flagged",
                "value": int(eligible["transition_period_flag"].fillna(False).sum()),
                "notes": "eligible rows spanning within-period class transitions",
            },
            {
                "metric": "eligible_rows_with_lagged_class",
                "value": int(eligible["lagged_class_available"].fillna(False).sum()),
                "notes": "eligible rows with a prior-period class available for lagged exposure sensitivity",
            },
        ]
    )
    coverage.to_csv(coverage_path, sep="\t", index=False)
    return coverage_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Build programme-surveillance country-period panel.")
    parser.add_argument("--exposure", required=True, help="Workflow annual exposure index TSV")
    parser.add_argument("--prevalence", required=True, help="Workflow IPW prevalence TSV")
    parser.add_argument("--out", required=True, help="Output period panel TSV")
    parser.add_argument(
        "--origin-descendants-dir",
        default=str(DEFAULT_ORIGIN_DESCENDANTS_DIR),
        help="Directory containing origin_*.descendant_tips.tsv bridge files",
    )
    parser.add_argument("--bin-size", type=int, default=5, help="Country-period bin size in years")
    parser.add_argument(
        "--min-interpretable",
        type=int,
        default=5,
        help="Minimum interpretable genomes required for a primary-panel period row",
    )
    args = parser.parse_args()

    annual = prepare_annual_dataset(
        args.exposure,
        args.prevalence,
        origin_descendants_dir=args.origin_descendants_dir,
    )
    panel = build_period_panel(
        annual,
        bin_size=args.bin_size,
        min_interpretable=args.min_interpretable,
    )

    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    panel.to_csv(output_path, sep="\t", index=False)
    coverage_path = write_coverage_report(
        panel,
        output_path,
        min_interpretable=args.min_interpretable,
        bin_size=args.bin_size,
    )
    print(f"Wrote programme-surveillance country-period panel: {output_path}")
    print(f"Wrote programme-surveillance coverage report: {coverage_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
