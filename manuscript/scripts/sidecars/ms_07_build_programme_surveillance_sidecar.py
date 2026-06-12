#!/usr/bin/env python3
"""Build manuscript-facing programme-surveillance sidecar extracts."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


PRIMARY_MODEL_ID = "programme_country_period_5y_primary_programme_only_ipw"
PRIMARY_BRIDGE_MODEL_ID = "programme_country_period_5y_primary_programme_plus_bridge_ipw"
EXCLUDE_USA_MODEL_ID = "programme_country_period_5y_exclude_USA_programme_only_ipw"
HIGH_CONFIDENCE_MODEL_ID = "programme_country_period_5y_high_confidence_programme_only_ipw"
TRANSITION_EXCLUDED_MODEL_ID = "programme_country_period_5y_transition_excluded_programme_only_ipw"
LAGGED_MODEL_ID = "programme_country_period_5y_lagged_programme_only_ipw"
SUPPORTING_MODEL_ID = "programme_country_period_3y_supporting_programme_only_ipw"
PROGRAM_CLASS_ORDER = [
    "wp_only_or_pre_ap",
    "routine_ap_prn_negative",
    "routine_ap_mixed",
    "routine_ap_prn_positive",
    "transition_mixed_within_period",
    "routine_ap_unknown",
    "unknown",
]


def normalize_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and np.isnan(value):
        return ""
    return str(value).strip()


def parse_bool(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1", "yes", "y", "t"})
    )


def coerce_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def program_class_label(value: object) -> str:
    return normalize_text(value) or "unknown"


def class_sort_key(value: object) -> int:
    text = normalize_text(value)
    try:
        return PROGRAM_CLASS_ORDER.index(text)
    except ValueError:
        return len(PROGRAM_CLASS_ORDER)


def parse_program_class_term(term: object) -> str:
    text = normalize_text(term)
    if "[T." not in text or not text.endswith("]"):
        return ""
    return text.split("[T.", 1)[1][:-1]


def prepare_panel(path: str) -> pd.DataFrame:
    panel = pd.read_csv(path, sep="\t", dtype=str)
    numeric_columns = [
        "period_start",
        "period_end",
        "reported_cases_period",
        "response_n_genomes_total",
        "response_n_genomes_prn_interpretable",
        "response_n_prn_disrupted",
        "response_n_missing_outcomes",
        "response_ipw_weight_total",
        "response_ipw_successes_est",
        "response_ipw_prevalence",
        "response_naive_prevalence",
        "response_boundary_lower_prevalence",
        "response_boundary_upper_prevalence",
        "genomes_per_case_effective",
        "post_covid_period",
        "country_row_share",
        "share_years_prn_positive_within_period",
        "years_since_prn_positive_routine_use",
        "years_since_any_ap_use",
        "n_new_origins_detected_period",
        "n_active_origin_clades_period",
        "first_local_origin_year",
        "first_prn_detection_year",
        "years_since_first_local_origin",
        "years_since_first_prn_detection",
    ]
    for column in numeric_columns:
        panel[column] = coerce_numeric(panel.get(column, pd.Series(dtype=str)))
    bool_columns = [
        "primary_panel_eligible",
        "period_contains_conflict",
        "transition_period_flag",
        "lagged_class_available",
        "has_local_origin_by_period_end",
        "has_prn_detection_by_period_end",
    ]
    for column in bool_columns:
        panel[column] = parse_bool(panel.get(column, pd.Series(dtype=str)))
    return panel


def prepare_model_results(path: str) -> pd.DataFrame:
    results = pd.read_csv(path, sep="\t", dtype=str)
    numeric_columns = [
        "effect_estimate",
        "ci_lower",
        "ci_upper",
        "p_value",
        "q_value",
        "n_rows",
        "n_countries",
    ]
    for column in numeric_columns:
        results[column] = coerce_numeric(results.get(column, pd.Series(dtype=str)))
    return results


def prepare_model_diagnostics(path: str) -> pd.DataFrame:
    diagnostics = pd.read_csv(path, sep="\t", dtype=str)
    diagnostics["converged"] = parse_bool(diagnostics.get("converged", pd.Series(dtype=str)))
    diagnostics["n_rows"] = coerce_numeric(diagnostics.get("n_rows", pd.Series(dtype=str)))
    diagnostics["n_countries"] = coerce_numeric(diagnostics.get("n_countries", pd.Series(dtype=str)))
    return diagnostics


def prepare_audit(path: str) -> pd.DataFrame:
    audit = pd.read_csv(path, sep="\t", dtype=str)
    audit["year"] = coerce_numeric(audit.get("year", pd.Series(dtype=str)))
    audit["program_formulation_conflict"] = parse_bool(
        audit.get("program_formulation_conflict", pd.Series(dtype=str))
    )
    return audit


def load_coverage(path: str) -> dict[str, str]:
    coverage = pd.read_csv(path, sep="\t", dtype=str)
    return {
        normalize_text(row.metric): normalize_text(row.value)
        for row in coverage.itertuples(index=False)
    }


def eligible_panel(panel: pd.DataFrame) -> pd.DataFrame:
    return panel.loc[panel["primary_panel_eligible"]].copy()


def programme_class_column(frame: pd.DataFrame) -> str:
    if "program_formulation_class_concurrent" in frame.columns:
        return "program_formulation_class_concurrent"
    return "program_formulation_class"


def build_panel_eligibility(panel: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "analysis_panel",
        "country_iso3",
        "country_name",
        "period_label",
        "period_start",
        "period_end",
        "n_years_observed",
        "response_n_genomes_prn_interpretable",
        "response_n_prn_disrupted",
        "response_ipw_prevalence",
        "response_naive_prevalence",
        "response_boundary_lower_prevalence",
        "response_boundary_upper_prevalence",
        "reported_cases_period",
        "genomes_per_case_effective",
        "program_formulation_class_concurrent",
        "program_formulation_class_lagged",
        "lagged_class_available",
        "transition_period_flag",
        "transition_period_type",
        "unique_program_classes_within_period",
        "share_years_prn_positive_within_period",
        "n_new_origins_detected_period",
        "n_active_origin_clades_period",
        "years_since_first_local_origin",
        "years_since_first_prn_detection",
        "formulation_confidence_period",
        "period_contains_conflict",
        "dominant_exposure_precedence_rule",
        "country_row_share",
        "primary_panel_eligible",
        "exclusion_reason",
        "notes",
    ]
    available = [column for column in columns if column in panel.columns]
    return panel[available].sort_values(
        ["primary_panel_eligible", "country_iso3", "period_start"],
        ascending=[False, True, True],
    )


def build_programme_class_summary(panel: pd.DataFrame) -> pd.DataFrame:
    eligible = eligible_panel(panel)
    if eligible.empty:
        return pd.DataFrame()
    class_column = programme_class_column(eligible)

    rows: list[dict[str, object]] = []
    for program_class, frame in eligible.groupby(class_column, dropna=False):
        total_weight = float(frame["response_ipw_weight_total"].fillna(0).sum())
        total_successes = float(frame["response_ipw_successes_est"].fillna(0).sum())
        total_interpretable = float(frame["response_n_genomes_prn_interpretable"].fillna(0).sum())
        total_disrupted = float(frame["response_n_prn_disrupted"].fillna(0).sum())
        rows.append(
            {
                "program_formulation_class": normalize_text(program_class),
                "n_period_rows": int(len(frame)),
                "n_countries": int(frame["country_iso3"].nunique()),
                "country_iso3_list": ",".join(sorted(frame["country_iso3"].dropna().astype(str).unique())),
                "total_interpretable_genomes": int(total_interpretable),
                "total_prn_disrupted": int(total_disrupted),
                "pooled_naive_prevalence": (
                    np.nan if total_interpretable <= 0 else total_disrupted / total_interpretable
                ),
                "pooled_ipw_prevalence": np.nan if total_weight <= 0 else total_successes / total_weight,
                "mean_period_ipw_prevalence": frame["response_ipw_prevalence"].mean(),
                "median_period_ipw_prevalence": frame["response_ipw_prevalence"].median(),
                "high_confidence_rows": int(
                    frame["formulation_confidence_period"].fillna("").astype(str).str.lower().eq("high").sum()
                ),
                "rows_with_precedence_conflict": int(frame["period_contains_conflict"].fillna(False).sum()),
                "transition_flagged_rows": int(frame.get("transition_period_flag", pd.Series(False, index=frame.index)).fillna(False).sum()),
                "lagged_class_available_rows": int(frame.get("lagged_class_available", pd.Series(False, index=frame.index)).fillna(False).sum()),
                "year_window": (
                    f"{int(frame['period_start'].min())}-{int(frame['period_end'].max())}"
                    if not frame.empty
                    else ""
                ),
            }
        )
    summary = pd.DataFrame.from_records(rows)
    return summary.sort_values(
        by=["program_formulation_class"],
        key=lambda values: values.map(class_sort_key),
    ).reset_index(drop=True)


def build_country_influence(panel: pd.DataFrame) -> pd.DataFrame:
    eligible = eligible_panel(panel)
    if eligible.empty:
        return pd.DataFrame()
    class_column = programme_class_column(eligible)

    rows: list[dict[str, object]] = []
    max_share = float(eligible["country_row_share"].max())
    for country_iso3, frame in eligible.groupby("country_iso3", dropna=False):
        total_weight = float(frame["response_ipw_weight_total"].fillna(0).sum())
        total_successes = float(frame["response_ipw_successes_est"].fillna(0).sum())
        total_interpretable = float(frame["response_n_genomes_prn_interpretable"].fillna(0).sum())
        total_disrupted = float(frame["response_n_prn_disrupted"].fillna(0).sum())
        rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": normalize_text(frame["country_name"].mode().iloc[0]) if not frame["country_name"].dropna().empty else "",
                "eligible_rows": int(len(frame)),
                "country_row_share": float(frame["country_row_share"].iloc[0]),
                "program_formulation_classes": ",".join(
                    sorted(
                        {
                            normalize_text(value)
                            for value in frame[class_column].tolist()
                            if normalize_text(value)
                        }
                    )
                ),
                "high_confidence_rows": int(
                    frame["formulation_confidence_period"].fillna("").astype(str).str.lower().eq("high").sum()
                ),
                "rows_with_precedence_conflict": int(frame["period_contains_conflict"].fillna(False).sum()),
                "transition_flagged_rows": int(frame.get("transition_period_flag", pd.Series(False, index=frame.index)).fillna(False).sum()),
                "lagged_class_available_rows": int(frame.get("lagged_class_available", pd.Series(False, index=frame.index)).fillna(False).sum()),
                "total_interpretable_genomes": int(total_interpretable),
                "total_prn_disrupted": int(total_disrupted),
                "pooled_naive_prevalence": (
                    np.nan if total_interpretable <= 0 else total_disrupted / total_interpretable
                ),
                "pooled_ipw_prevalence": np.nan if total_weight <= 0 else total_successes / total_weight,
                "mean_period_ipw_prevalence": frame["response_ipw_prevalence"].mean(),
                "period_labels": ",".join(frame["period_label"].astype(str).tolist()),
                "dominance_flag": (
                    "largest_country_contributor"
                    if float(frame["country_row_share"].iloc[0]) == max_share
                    else ""
                ),
            }
        )
    return pd.DataFrame.from_records(rows).sort_values(
        ["country_row_share", "country_iso3"], ascending=[False, True]
    )


def classify_effect_stability(row: pd.Series) -> str:
    estimate = row["effect_estimate"]
    lower = row["ci_lower"]
    upper = row["ci_upper"]
    if pd.isna(estimate) or pd.isna(lower) or pd.isna(upper):
        return "missing"
    if max(abs(estimate), abs(lower), abs(upper)) >= 10:
        return "extreme_or_quasi_separated"
    if lower <= 0 <= upper:
        return "imprecise_crosses_null"
    return "stable_directional_signal"


def prevalence_stability(row: pd.Series) -> str:
    estimate = row["effect_estimate"]
    lower = row["ci_lower"]
    upper = row["ci_upper"]
    if pd.isna(estimate) or pd.isna(lower) or pd.isna(upper):
        return "missing"
    if lower <= 0 <= upper:
        return "imprecise_crosses_null"
    return "stable_directional_signal"


def safe_odds_ratio(value: float) -> float:
    if pd.isna(value) or abs(float(value)) >= 10:
        return np.nan
    return float(np.exp(value))


def build_model_terms(results: pd.DataFrame, diagnostics: pd.DataFrame) -> pd.DataFrame:
    coeffs = results.loc[
        results["result_type"].eq("coefficient")
        & results["estimate_term"].astype(str).str.contains("program_formulation_class")
    ].copy()
    if coeffs.empty:
        return pd.DataFrame()

    coeffs["program_class"] = coeffs["estimate_term"].map(parse_program_class_term)
    coeffs["term_label"] = coeffs["program_class"].map(lambda value: f"{value} vs wp_only_or_pre_ap")
    coeffs["estimate_stability"] = coeffs.apply(classify_effect_stability, axis=1)
    coeffs["odds_ratio"] = coeffs["effect_estimate"].map(safe_odds_ratio)
    coeffs["odds_ratio_ci_lower"] = coeffs["ci_lower"].map(safe_odds_ratio)
    coeffs["odds_ratio_ci_upper"] = coeffs["ci_upper"].map(safe_odds_ratio)

    diagnostic_subset = diagnostics[
        ["model_id", "converged", "covariance_type", "notes"]
    ].rename(columns={"notes": "diagnostic_notes"})
    coeffs = coeffs.merge(diagnostic_subset, on="model_id", how="left")

    model_order = {
        PRIMARY_MODEL_ID: 0,
        PRIMARY_BRIDGE_MODEL_ID: 1,
        EXCLUDE_USA_MODEL_ID: 2,
        HIGH_CONFIDENCE_MODEL_ID: 3,
        TRANSITION_EXCLUDED_MODEL_ID: 4,
        LAGGED_MODEL_ID: 5,
        SUPPORTING_MODEL_ID: 6,
    }
    return coeffs.sort_values(
        by=["model_id", "term_label"],
        key=lambda values: values.map(model_order).fillna(99) if values.name == "model_id" else values,
    )


def build_adjusted_prevalence_summary(results: pd.DataFrame) -> pd.DataFrame:
    adjusted = results.loc[
        results["result_type"].isin(["adjusted_prevalence", "adjusted_prevalence_difference"])
        & results["response_track"].eq("ipw")
    ].copy()
    if adjusted.empty:
        return adjusted
    adjusted["estimate_stability"] = adjusted.apply(prevalence_stability, axis=1)
    return adjusted.sort_values(
        ["model_spec", "sensitivity_label", "result_type", "program_class", "excluded_country_iso3"],
        ascending=[True, True, True, True, True],
    )


def effect_direction(value: object) -> str:
    numeric = float(value)
    if numeric > 0:
        return "positive"
    if numeric < 0:
        return "negative"
    return "null"


def build_leave_one_country_out_summary(results: pd.DataFrame) -> pd.DataFrame:
    required = {"sensitivity_label", "model_spec", "response_track", "result_type"}
    if not required.issubset(results.columns):
        return pd.DataFrame()
    loo = results.loc[
        results["sensitivity_label"].eq("leave_one_country_out")
        & results["model_spec"].eq("programme_only")
        & results["response_track"].eq("ipw")
        & results["result_type"].eq("adjusted_prevalence_difference")
    ].copy()
    if loo.empty:
        return loo

    primary = results.loc[
        results["model_id"].eq(PRIMARY_MODEL_ID)
        & results["result_type"].eq("adjusted_prevalence_difference")
    ][["estimate_term", "effect_estimate"]].rename(columns={"effect_estimate": "primary_effect_estimate"})
    loo = loo.merge(primary, on="estimate_term", how="left")
    loo["same_direction_as_primary"] = np.where(
        loo["primary_effect_estimate"].notna(),
        loo["effect_estimate"].map(effect_direction) == loo["primary_effect_estimate"].map(effect_direction),
        np.nan,
    )
    return loo.sort_values(["program_class", "excluded_country_iso3"]).reset_index(drop=True)


def build_precedence_conflicts(audit: pd.DataFrame) -> pd.DataFrame:
    conflicts = audit.loc[audit["program_formulation_conflict"]].copy()
    if conflicts.empty:
        return conflicts
    columns = [
        "country_iso3",
        "country_name",
        "year",
        "vaccine_program_type_effective",
        "primary_series_formulation",
        "booster_formulation",
        "prn_in_vaccine_curated",
        "formulation_confidence",
        "program_formulation_class",
        "exposure_precedence_rule",
        "prn_in_vaccine_source_class",
        "formulation_source_name",
        "formulation_source_url",
        "formulation_notes",
    ]
    available = [column for column in columns if column in conflicts.columns]
    return conflicts[available].sort_values(["country_iso3", "year"]).reset_index(drop=True)


def build_curation_priorities(panel: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for country_iso3, frame in panel.groupby("country_iso3", dropna=False):
        unknown_class_rows = frame["exclusion_reason"].fillna("").astype(str).str.startswith("class_excluded_")
        threshold_ready_unknown_rows = unknown_class_rows & frame["response_n_genomes_prn_interpretable"].fillna(0).ge(5)
        rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": normalize_text(frame["country_name"].mode().iloc[0]) if not frame["country_name"].dropna().empty else "",
                "unknown_class_rows": int(unknown_class_rows.sum()),
                "threshold_ready_unknown_rows": int(threshold_ready_unknown_rows.sum()),
                "unknown_class_interpretable_genomes": int(
                    frame.loc[unknown_class_rows, "response_n_genomes_prn_interpretable"].fillna(0).sum()
                ),
                "periods_unknown_class": ",".join(frame.loc[unknown_class_rows, "period_label"].astype(str).tolist()),
                "missing_cases_rows": int(frame["exclusion_reason"].fillna("").eq("missing_reported_cases").sum()),
                "below_threshold_rows": int(
                    frame["exclusion_reason"].fillna("").astype(str).str.startswith("below_interpretable_threshold_").sum()
                ),
                "priority_status": (
                    "would_expand_primary_panel_if_curated"
                    if int(threshold_ready_unknown_rows.sum()) > 0
                    else "count_limited_or_nonpanel_blocker"
                ),
            }
        )
    summary = pd.DataFrame.from_records(rows)
    summary = summary.loc[summary["unknown_class_rows"].gt(0)].copy()
    if summary.empty:
        return summary
    return summary.sort_values(
        ["threshold_ready_unknown_rows", "unknown_class_interpretable_genomes", "country_iso3"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


def high_confidence_strategy(diagnostics: pd.DataFrame) -> tuple[bool, str, str]:
    high_confidence_row = diagnostics.loc[diagnostics["model_id"].eq(HIGH_CONFIDENCE_MODEL_ID)].head(1)
    estimable = bool(not high_confidence_row.empty and bool(high_confidence_row["converged"].iloc[0]))
    if estimable:
        return True, "converged_high_confidence_subset_model", "adjusted_model_available"

    note = normalize_text(high_confidence_row["notes"].iloc[0]) if not high_confidence_row.empty else "not_run"
    lowered = note.lower()
    if "reference_class_all_zero_prevalence" in lowered or "singular matrix" in lowered:
        return False, note, "supplementary_descriptive_fallback_only"
    return False, note, "descriptive_only_no_adjusted_fallback"


def build_high_confidence_fallback(panel: pd.DataFrame, diagnostics: pd.DataFrame) -> pd.DataFrame:
    high_confidence = eligible_panel(panel)
    high_confidence = high_confidence.loc[
        high_confidence["formulation_confidence_period"].fillna("").astype(str).str.lower().eq("high")
    ].copy()
    if high_confidence.empty:
        return pd.DataFrame()

    _, note, strategy = high_confidence_strategy(diagnostics)
    fallback_panel = high_confidence.copy()
    fallback_panel["primary_panel_eligible"] = True
    summary = build_programme_class_summary(fallback_panel)
    if summary.empty:
        return summary
    summary.insert(0, "subset_label", "high_confidence_only")
    summary.insert(1, "high_confidence_strategy", strategy)
    summary.insert(2, "high_confidence_note", note)
    return summary


def build_two_stage_uncertainty_summary(path: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    numeric_columns = [
        "n_replicates",
        "bootstrap_mean",
        "bootstrap_median",
        "bootstrap_ci_lower",
        "bootstrap_ci_upper",
        "bootstrap_interval_width",
        "main_model_point_estimate",
        "main_model_ci_lower",
        "main_model_ci_upper",
        "main_model_ci_width",
        "propagated_ci_lower",
        "propagated_ci_upper",
        "propagated_interval_width",
    ]
    for column in numeric_columns:
        frame[column] = coerce_numeric(frame.get(column, pd.Series(dtype=str)))
    if "interval_narrower_than_single_model" in frame.columns:
        frame["interval_narrower_than_single_model"] = parse_bool(
            frame["interval_narrower_than_single_model"]
        )
    return frame.sort_values(["model_spec", "result_type", "estimate_term"]).reset_index(drop=True)


def mechanism_group(value: object) -> str:
    text = normalize_text(value).lower()
    if not text:
        return "Unknown"
    if text == "intact":
        return "Intact"
    if "is481" in text:
        return "IS481 insertion"
    if "inversion" in text or "rearrangement" in text:
        return "Inversion / rearrangement"
    if "insufficient" in text:
        return "Insufficient data"
    return "Other disruptions"


def build_representativeness_audit(manifest_path: str, tip_states_path: str) -> pd.DataFrame:
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    tip_states = pd.read_csv(tip_states_path, sep="\t", dtype=str)

    manifest["prn_interpretable"] = parse_bool(manifest.get("prn_interpretable", pd.Series(dtype=str)))
    manifest["year"] = coerce_numeric(manifest.get("year", pd.Series(dtype=str)))
    full = manifest.loc[manifest["prn_interpretable"]].copy()
    full["mechanism_group"] = full["prn_mechanism_call"].map(mechanism_group)
    full["year_band"] = np.where(
        full["year"].notna(),
        (full["year"] // 5 * 5).astype("Int64").astype(str) + "s",
        "unknown",
    )
    full["lineage_proxy"] = full.get("phylo_lineage", pd.Series(dtype=str)).fillna("")
    full["lineage_proxy"] = full["lineage_proxy"].where(
        full["lineage_proxy"].astype(str).str.len().gt(0),
        full.get("mlst_st", pd.Series(dtype=str)).fillna("").map(lambda value: f"MLST_{value}" if normalize_text(value) else "unknown"),
    )

    tree = tip_states.loc[
        tip_states.get("is_reference", pd.Series(dtype=str)).fillna("").astype(str).str.lower().ne("true")
        & tip_states["sample_id_canonical"].fillna("").astype(str).str.len().gt(0)
        & tip_states["prn_state"].fillna("").isin(["intact", "disrupted"])
    ].copy()
    tree["year"] = coerce_numeric(tree.get("year", pd.Series(dtype=str)))
    tree["mechanism_group"] = tree["observed_prn_mechanism_call"].map(mechanism_group)
    tree["year_band"] = np.where(
        tree["year"].notna(),
        (tree["year"] // 5 * 5).astype("Int64").astype(str) + "s",
        "unknown",
    )
    tree["lineage_proxy"] = tree.get("phylo_lineage", pd.Series(dtype=str)).fillna("")
    tree["lineage_proxy"] = tree["lineage_proxy"].where(
        tree["lineage_proxy"].astype(str).str.len().gt(0),
        tree.get("mlst_st", pd.Series(dtype=str)).fillna("").map(lambda value: f"MLST_{value}" if normalize_text(value) else "unknown"),
    )

    rows: list[dict[str, object]] = []
    dimensions = {
        "country_iso3": "country",
        "year_band": "year_band",
        "lineage_proxy": "lineage_proxy",
        "mechanism_group": "mechanism_group",
    }
    for column, label in dimensions.items():
        full_counts = full[column].fillna("unknown").astype(str).value_counts(dropna=False)
        tree_counts = tree[column].fillna("unknown").astype(str).value_counts(dropna=False)
        total_full = int(full_counts.sum())
        total_tree = int(tree_counts.sum())
        categories = sorted(set(full_counts.index).union(set(tree_counts.index)))
        for category in categories:
            full_count = int(full_counts.get(category, 0))
            tree_count = int(tree_counts.get(category, 0))
            full_fraction = float(full_count / total_full) if total_full else np.nan
            tree_fraction = float(tree_count / total_tree) if total_tree else np.nan
            rows.append(
                {
                    "comparison_dimension": label,
                    "category": category,
                    "full_interpretable_count": full_count,
                    "full_interpretable_fraction": full_fraction,
                    "tree_subset_count": tree_count,
                    "tree_subset_fraction": tree_fraction,
                    "absolute_fraction_gap": np.nan if np.isnan(full_fraction) or np.isnan(tree_fraction) else abs(tree_fraction - full_fraction),
                    "fraction_ratio_tree_to_full": (
                        np.nan if np.isnan(full_fraction) or full_fraction == 0 else tree_fraction / full_fraction
                    ),
                }
            )
    return pd.DataFrame.from_records(rows).sort_values(
        ["comparison_dimension", "absolute_fraction_gap", "category"],
        ascending=[True, False, True],
    )


def build_validation_summary(
    mechanism_calls_path: str,
    read_validation_path: str,
    validation_evidence_path: str,
) -> pd.DataFrame:
    mechanism_calls = pd.read_csv(mechanism_calls_path, sep="\t", dtype=str)
    read_validation = pd.read_csv(read_validation_path, sep="\t", dtype=str)
    validation_evidence = pd.read_csv(validation_evidence_path, sep="\t", dtype=str)

    disrupted = mechanism_calls.loc[
        ~mechanism_calls["prn_mechanism_call"].fillna("").isin(["intact", "insufficient_data"])
    ].copy()
    disrupted["mechanism_group"] = disrupted["prn_mechanism_call"].map(mechanism_group)
    validation_lookup = read_validation[["sample_id_canonical", "read_validation_status"]].drop_duplicates()
    disrupted = disrupted.merge(
        validation_lookup,
        on="sample_id_canonical",
        how="left",
        suffixes=("", "_from_reads"),
    )
    if "read_validation_status_from_reads" in disrupted.columns:
        existing_status = disrupted.get("read_validation_status", pd.Series(dtype=str)).fillna("").astype(str)
        replace_mask = existing_status.str.strip().str.lower().isin({"", "not_run", "tool_output_missing"})
        disrupted["read_validation_status"] = existing_status.where(
            ~replace_mask,
            disrupted["read_validation_status_from_reads"].fillna(existing_status),
        )

    supported_statuses = {"supported", "supported_concordant", "supported_candidate"}
    no_signal_statuses = {"no_prn_is_signal_detected"}
    unresolved_statuses = {"unresolved", "tool_output_missing", "not_run", ""}

    evidence_lookup = (
        validation_evidence.sort_values("validation_level")
        .drop_duplicates(subset=["mechanism_group"], keep="first")
        .set_index("mechanism_group")
    )

    rows: list[dict[str, object]] = []
    for mechanism_group_name, frame in disrupted.groupby("mechanism_group", dropna=False):
        statuses = frame["read_validation_status"].fillna("").astype(str).str.strip().str.lower()
        evidence_row = evidence_lookup.loc[mechanism_group_name] if mechanism_group_name in evidence_lookup.index else None
        rows.append(
            {
                "mechanism_group": mechanism_group_name,
                "n_disrupted_samples": int(len(frame)),
                "n_with_validation_rows": int(statuses.ne("").sum()),
                "n_supported_like": int(statuses.isin(supported_statuses).sum()),
                "n_no_signal": int(statuses.isin(no_signal_statuses).sum()),
                "n_unresolved_or_not_run": int(statuses.isin(unresolved_statuses).sum()),
                "representative_validation_level": "" if evidence_row is None else normalize_text(evidence_row["validation_level"]),
                "representative_sample_id_canonical": "" if evidence_row is None else normalize_text(evidence_row["sample_id_canonical"]),
                "representative_supporting_evidence": "" if evidence_row is None else normalize_text(evidence_row["supporting_read_or_public_longread"]),
            }
        )
    return pd.DataFrame.from_records(rows).sort_values(
        by=["mechanism_group"],
        key=lambda values: values.map(class_sort_key) if values.name == "mechanism_group" else values,
    )


def difference_sign_lookup(frame: pd.DataFrame) -> dict[tuple[str, str], str]:
    lookup = {}
    for row in frame.itertuples(index=False):
        if pd.isna(row.effect_estimate):
            continue
        lookup[(normalize_text(row.program_class), normalize_text(row.reference_class))] = (
            "positive" if row.effect_estimate > 0 else "negative" if row.effect_estimate < 0 else "null"
        )
    return lookup


def bridge_attenuation_detected(adjusted_prevalence: pd.DataFrame) -> bool:
    primary = adjusted_prevalence.loc[
        adjusted_prevalence["model_id"].eq(PRIMARY_MODEL_ID)
        & adjusted_prevalence["result_type"].eq("adjusted_prevalence_difference")
    ][["program_class", "reference_class", "effect_estimate"]].rename(columns={"effect_estimate": "primary_effect"})
    bridge = adjusted_prevalence.loc[
        adjusted_prevalence["model_id"].eq(PRIMARY_BRIDGE_MODEL_ID)
        & adjusted_prevalence["result_type"].eq("adjusted_prevalence_difference")
    ][["program_class", "reference_class", "effect_estimate"]].rename(columns={"effect_estimate": "bridge_effect"})
    merged = primary.merge(bridge, on=["program_class", "reference_class"], how="inner")
    if merged.empty:
        return False
    return bool(
        (
            merged["primary_effect"].fillna(0).map(np.sign) != merged["bridge_effect"].fillna(0).map(np.sign)
        ).any()
        or (
            merged["bridge_effect"].abs() < (0.75 * merged["primary_effect"].abs())
        ).any()
    )


def build_readiness_summary(
    panel_5y: pd.DataFrame,
    coverage_5y: dict[str, str],
    coverage_3y: dict[str, str],
    model_terms: pd.DataFrame,
    adjusted_prevalence: pd.DataFrame,
    diagnostics: pd.DataFrame,
    uncertainty_summary: pd.DataFrame | None = None,
    curation_priorities: pd.DataFrame | None = None,
) -> pd.DataFrame:
    if uncertainty_summary is None:
        uncertainty_summary = pd.DataFrame()
    if curation_priorities is None:
        curation_priorities = pd.DataFrame()
    eligible_countries_5y = int(float(coverage_5y.get("eligible_countries", "0") or 0))
    eligible_rows_5y = int(float(coverage_5y.get("eligible_rows", "0") or 0))
    max_country_share_5y = float(coverage_5y.get("max_country_row_share", "0") or 0)
    eligible_rows_high_confidence_5y = int(float(coverage_5y.get("eligible_rows_high_confidence", "0") or 0))
    eligible_countries_3y = int(float(coverage_3y.get("eligible_countries", "0") or 0))
    coverage_target_met = eligible_countries_5y >= 10 and max_country_share_5y <= 0.35

    primary_terms = model_terms.loc[model_terms["model_id"].eq(PRIMARY_MODEL_ID)].copy()
    primary_extreme = bool(primary_terms["estimate_stability"].eq("extreme_or_quasi_separated").any())

    difference_rows = adjusted_prevalence.loc[
        adjusted_prevalence["result_type"].eq("adjusted_prevalence_difference")
    ].copy()
    primary_signs = difference_sign_lookup(difference_rows.loc[difference_rows["model_id"].eq(PRIMARY_MODEL_ID)])
    exclude_usa_signs = difference_sign_lookup(difference_rows.loc[difference_rows["model_id"].eq(EXCLUDE_USA_MODEL_ID)])
    transition_signs = difference_sign_lookup(difference_rows.loc[difference_rows["model_id"].eq(TRANSITION_EXCLUDED_MODEL_ID)])
    lagged_signs = difference_sign_lookup(difference_rows.loc[difference_rows["model_id"].eq(LAGGED_MODEL_ID)])

    comparable_keys = set(primary_signs).intersection(exclude_usa_signs)
    exclude_usa_sign_flip = any(primary_signs[key] != exclude_usa_signs[key] for key in comparable_keys)
    comparable_transition = set(primary_signs).intersection(transition_signs)
    transition_sign_flip = any(primary_signs[key] != transition_signs[key] for key in comparable_transition)
    comparable_lagged = set(primary_signs).intersection(lagged_signs)
    lagged_sign_flip = any(primary_signs[key] != lagged_signs[key] for key in comparable_lagged)

    loo = build_leave_one_country_out_summary(adjusted_prevalence)
    leave_one_country_out_flip_count = int((loo["same_direction_as_primary"] == False).sum()) if not loo.empty else 0

    bootstrap_narrower = bool(
        "interval_narrower_than_single_model" in uncertainty_summary.columns
        and uncertainty_summary["interval_narrower_than_single_model"].fillna(False).astype(bool).any()
    )
    bridge_attenuation = bridge_attenuation_detected(adjusted_prevalence)

    high_confidence_estimable, high_confidence_note, high_confidence_strategy_label = high_confidence_strategy(
        diagnostics
    )

    eligible_5y = eligible_panel(panel_5y)
    conflicts = int(eligible_5y["period_contains_conflict"].fillna(False).sum()) if not eligible_5y.empty else 0
    transition_rows = int(
        eligible_5y.get("transition_period_flag", pd.Series(False, index=eligible_5y.index)).fillna(False).sum()
    ) if not eligible_5y.empty else 0
    lagged_available_rows = int(
        eligible_5y.get("lagged_class_available", pd.Series(False, index=eligible_5y.index)).fillna(False).sum()
    ) if not eligible_5y.empty else 0

    top_curation_target = ""
    top_curation_priority = ""
    if not curation_priorities.empty:
        top_curation_target = normalize_text(curation_priorities.iloc[0]["country_iso3"])
        top_curation_priority = normalize_text(curation_priorities.iloc[0]["priority_status"])

    if not coverage_target_met:
        headline_recommendation = "expand_panel_before_programme_ecology_headline"
        recommendation_note = "5-year panel still misses the minimum country-diversity target."
    elif (
        primary_extreme
        or exclude_usa_sign_flip
        or transition_sign_flip
        or lagged_sign_flip
        or leave_one_country_out_flip_count > 0
    ):
        headline_recommendation = "country_dependent_or_design_sensitive_ecological_signal_only"
        recommendation_note = (
            "Primary contrasts remain sensitive to country composition or exposure design and should not anchor a policy-facing headline."
        )
    elif not high_confidence_estimable:
        headline_recommendation = "country_dependent_or_design_sensitive_ecological_signal_only"
        recommendation_note = "High-confidence subset remains support-only."
    else:
        headline_recommendation = "cautious_programme_heterogeneity_summary"
        recommendation_note = "Coverage target is met and upgraded sensitivity checks remain directionally aligned."

    summary = pd.DataFrame(
        [
            {
                "headline_recommendation": headline_recommendation,
                "recommendation_note": recommendation_note,
                "eligible_rows_5y": eligible_rows_5y,
                "eligible_countries_5y": eligible_countries_5y,
                "eligible_rows_high_confidence_5y": eligible_rows_high_confidence_5y,
                "eligible_countries_3y": eligible_countries_3y,
                "max_country_row_share_5y": max_country_share_5y,
                "coverage_target_met": coverage_target_met,
                "primary_extreme_coefficients": primary_extreme,
                "exclude_usa_sign_flip": exclude_usa_sign_flip,
                "transition_excluded_sign_flip": transition_sign_flip,
                "lagged_exposure_sign_flip": lagged_sign_flip,
                "leave_one_country_out_flip_count": leave_one_country_out_flip_count,
                "bridge_attenuation_detected": bridge_attenuation,
                "bootstrap_narrower_than_single_model": bootstrap_narrower,
                "high_confidence_estimable": high_confidence_estimable,
                "high_confidence_note": high_confidence_note,
                "high_confidence_strategy": high_confidence_strategy_label,
                "eligible_rows_with_precedence_conflict_5y": conflicts,
                "eligible_rows_transition_flagged_5y": transition_rows,
                "eligible_rows_with_lagged_class_5y": lagged_available_rows,
                "top_curation_target": top_curation_target,
                "top_curation_priority": top_curation_priority,
            }
        ]
    )
    return summary


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build manuscript-facing programme-surveillance sidecar extracts.")
    parser.add_argument("--panel-5y", default="outputs/workflow/epi/programme_country_period_panel_5y.tsv")
    parser.add_argument("--panel-3y", default="outputs/workflow/epi/programme_country_period_panel_3y.tsv")
    parser.add_argument("--coverage-5y", default="outputs/workflow/epi/programme_country_period_panel_5y_coverage.tsv")
    parser.add_argument("--coverage-3y", default="outputs/workflow/epi/programme_country_period_panel_3y_coverage.tsv")
    parser.add_argument("--model-results", default="outputs/workflow/epi/programme_program_model_results.tsv")
    parser.add_argument("--model-diagnostics", default="outputs/workflow/epi/programme_program_model_diagnostics.tsv")
    parser.add_argument("--precedence-audit", default="outputs/workflow/epi/formulation_precedence_audit.tsv")
    parser.add_argument("--two-stage-uncertainty", default="outputs/workflow/epi/programme_two_stage_uncertainty_summary.tsv")
    parser.add_argument("--manifest", default="outputs/workflow/manifest/manifest.tsv")
    parser.add_argument("--tip-states", default="outputs/workflow/asr/tip_states.tsv")
    parser.add_argument("--mechanism-calls", default="manuscript/figure_data/fig02_prn_mechanism_calls.tsv")
    parser.add_argument("--read-validation", default="manuscript/figure_data/figure6_read_validation.tsv")
    parser.add_argument("--validation-evidence", default="manuscript/figure_data/validation_evidence.tsv")
    parser.add_argument("--outdir", default="manuscript/figure_data")
    args = parser.parse_args()

    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)

    panel_5y = prepare_panel(args.panel_5y)
    panel_3y = prepare_panel(args.panel_3y)
    model_results = prepare_model_results(args.model_results)
    model_diagnostics = prepare_model_diagnostics(args.model_diagnostics)
    audit = prepare_audit(args.precedence_audit)
    coverage_5y = load_coverage(args.coverage_5y)
    coverage_3y = load_coverage(args.coverage_3y)
    uncertainty_summary = build_two_stage_uncertainty_summary(args.two_stage_uncertainty)

    panel_eligibility = build_panel_eligibility(panel_5y)
    programme_class_summary = build_programme_class_summary(panel_5y)
    country_influence = build_country_influence(panel_5y)
    model_terms = build_model_terms(model_results, model_diagnostics)
    adjusted_prevalence = build_adjusted_prevalence_summary(model_results)
    leave_one_country_out = build_leave_one_country_out_summary(model_results)
    precedence_conflicts = build_precedence_conflicts(audit)
    curation_priorities = build_curation_priorities(panel_5y)
    high_confidence_fallback = build_high_confidence_fallback(panel_5y, model_diagnostics)
    representativeness_audit = build_representativeness_audit(args.manifest, args.tip_states)
    validation_summary = build_validation_summary(
        args.mechanism_calls,
        args.read_validation,
        args.validation_evidence,
    )
    readiness_summary = build_readiness_summary(
        panel_5y,
        coverage_5y,
        coverage_3y,
        model_terms,
        adjusted_prevalence,
        model_diagnostics,
        uncertainty_summary,
        curation_priorities,
    )

    write_tsv(panel_eligibility, output_dir / "supplementary_programme_panel_eligibility.tsv")
    write_tsv(programme_class_summary, output_dir / "supplementary_programme_class_summary.tsv")
    write_tsv(country_influence, output_dir / "supplementary_programme_country_influence.tsv")
    write_tsv(model_terms, output_dir / "supplementary_programme_model_terms.tsv")
    write_tsv(adjusted_prevalence, output_dir / "supplementary_programme_adjusted_prevalence.tsv")
    write_tsv(leave_one_country_out, output_dir / "supplementary_programme_leave_one_country_out.tsv")
    write_tsv(uncertainty_summary, output_dir / "supplementary_programme_two_stage_uncertainty.tsv")
    write_tsv(precedence_conflicts, output_dir / "supplementary_programme_precedence_conflicts.tsv")
    write_tsv(curation_priorities, output_dir / "supplementary_programme_curation_priorities.tsv")
    write_tsv(high_confidence_fallback, output_dir / "supplementary_programme_high_confidence_fallback.tsv")
    write_tsv(representativeness_audit, output_dir / "supplementary_programme_representativeness_audit.tsv")
    write_tsv(validation_summary, output_dir / "supplementary_programme_validation_summary.tsv")
    write_tsv(readiness_summary, output_dir / "supplementary_programme_readiness.tsv")

    print(f"Wrote programme-surveillance sidecar extracts to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
