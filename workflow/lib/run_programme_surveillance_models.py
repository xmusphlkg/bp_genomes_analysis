#!/usr/bin/env python3
"""Fit programme-surveillance ecological models on country-period panels."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import os
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from patsy import build_design_matrices


CLASS_ORDER = [
    "wp_only_or_pre_ap",
    "routine_ap_prn_negative",
    "routine_ap_mixed",
    "routine_ap_prn_positive",
]
PROGRAMME_MODELS_MAX_WORKERS_ENV = "PROGRAMME_MODELS_MAX_WORKERS"
PROGRAMME_MODELS_DEFAULT_MAX_WORKERS = 32


def resolve_parallel_workers(
    task_count: int | None = None,
    *,
    requested_max_workers: int | None = None,
    cpu_count: int | None = None,
) -> int:
    available_cpu = max(1, cpu_count or os.cpu_count() or 1)
    workers = (
        requested_max_workers
        if requested_max_workers is not None
        else min(available_cpu, PROGRAMME_MODELS_DEFAULT_MAX_WORKERS)
    )
    workers = max(1, min(int(workers), available_cpu))
    if task_count is not None:
        workers = min(workers, max(1, int(task_count)))
    return workers


RESULT_COLUMNS = [
    "model_id",
    "analysis_panel",
    "sensitivity_label",
    "model_spec",
    "response_track",
    "exposure_column",
    "excluded_country_iso3",
    "country_filter",
    "n_rows",
    "n_countries",
    "year_window",
    "model_family",
    "weighting_scheme",
    "result_type",
    "estimate_term",
    "program_class",
    "reference_class",
    "effect_scale",
    "effect_estimate",
    "ci_lower",
    "ci_upper",
    "p_value",
    "q_value",
    "notes",
]

DIAGNOSTIC_COLUMNS = [
    "model_id",
    "analysis_panel",
    "sensitivity_label",
    "model_spec",
    "response_track",
    "exposure_column",
    "excluded_country_iso3",
    "n_rows",
    "n_countries",
    "year_window",
    "converged",
    "covariance_type",
    "log_likelihood",
    "notes",
]


def normalize_text(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def coerce_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    if not p_values:
        return []
    ordered = sorted(enumerate(p_values), key=lambda item: item[1])
    adjusted = [0.0] * len(p_values)
    running_min = 1.0
    total = len(p_values)
    for reverse_rank, (original_index, p_value) in enumerate(reversed(ordered), start=1):
        rank = total - reverse_rank + 1
        scaled = min(1.0, p_value * total / rank)
        running_min = min(running_min, scaled)
        adjusted[original_index] = running_min
    return adjusted


def standardize(series: pd.Series, *, fill_missing_with_zero: bool = False) -> pd.Series:
    valid = series.dropna()
    if valid.empty:
        output = pd.Series(np.nan, index=series.index)
    else:
        deviation = valid.std(ddof=0)
        if deviation == 0 or np.isnan(deviation):
            output = pd.Series(0.0, index=series.index)
        else:
            output = (series - valid.mean()) / deviation
    if fill_missing_with_zero:
        return output.fillna(0.0)
    return output


def safe_token(value: object) -> str:
    return normalize_text(value).lower().replace(" ", "_")


def resolve_env_max_workers() -> int | None:
    env_text = str(os.environ.get(PROGRAMME_MODELS_MAX_WORKERS_ENV, "")).strip()
    if not env_text:
        return None
    try:
        return int(env_text)
    except ValueError as exc:
        raise SystemExit(
            f"ERROR: {PROGRAMME_MODELS_MAX_WORKERS_ENV} must be an integer, got {env_text!r}"
        ) from exc


def response_track_label(response_track: str) -> str:
    return {
        "ipw": "ipw_fractional_response_no_pseudotrials",
        "naive": "naive_grouped_binomial_interpretable_trials",
        "boundary_lower": "boundary_lower_total_genome_trials",
        "boundary_upper": "boundary_upper_total_genome_trials",
    }.get(response_track, response_track)


def build_formula(exposure_column: str, *, include_bridge: bool, response_column: str = "analysis_response") -> str:
    base = (
        f"{response_column} ~ "
        f"C({exposure_column}, Treatment(reference='wp_only_or_pre_ap')) + "
        "log1p_reported_cases_z + genomes_per_case_z + post_covid_period"
    )
    if not include_bridge:
        return base
    return (
        base
        + " + n_new_origins_detected_period_z"
        + " + n_active_origin_clades_period_z"
        + " + years_since_first_local_origin_z"
        + " + years_since_first_prn_detection_z"
        + " + has_local_origin_by_period_end"
        + " + has_prn_detection_by_period_end"
    )


def prevalence_boundary_state(frame: pd.DataFrame) -> str:
    total_weight = float(frame["model_weight"].fillna(0).sum())
    total_successes = float(frame["model_successes"].fillna(0).sum())
    if total_weight <= 0:
        return "no_weight"
    if np.isclose(total_successes, 0.0):
        return "all_zero_prevalence"
    if np.isclose(total_successes, total_weight):
        return "all_one_prevalence"
    return "nonboundary"


def structural_diagnostic_tokens(frame: pd.DataFrame, exposure_column: str) -> list[str]:
    tokens: list[str] = []
    class_values = sorted(
        {
            normalize_text(value)
            for value in frame[exposure_column].tolist()
            if normalize_text(value)
        }
    )
    if len(class_values) <= 1:
        tokens.append("single_program_class_subset")
    reference_rows = frame.loc[frame[exposure_column].eq("wp_only_or_pre_ap")].copy()
    if not reference_rows.empty:
        reference_state = prevalence_boundary_state(reference_rows)
        if reference_state != "nonboundary":
            tokens.append(f"reference_class_{reference_state}")
    for program_class, class_frame in frame.groupby(exposure_column, dropna=False):
        class_label = safe_token(program_class or "unknown")
        if len(class_frame) == 1:
            tokens.append(f"class_{class_label}_single_row")
        class_state = prevalence_boundary_state(class_frame)
        if class_state != "nonboundary":
            tokens.append(f"class_{class_label}_{class_state}")
    return list(dict.fromkeys(tokens))


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
        panel[column] = (
            panel.get(column, pd.Series(dtype=str))
            .fillna("")
            .astype(str)
            .str.lower()
            .isin(["true", "1", "yes"])
        )
    panel["period_start"] = panel["period_start"].astype("Int64")
    panel["period_end"] = panel["period_end"].astype("Int64")
    return panel


def attach_response_track(frame: pd.DataFrame, response_track: str) -> pd.DataFrame:
    output = frame.copy()
    if response_track == "ipw":
        output["model_weight"] = output["response_ipw_weight_total"]
        output["model_successes"] = output["response_ipw_successes_est"]
        output["model_prevalence"] = output["response_ipw_prevalence"]
        output["analysis_response"] = output["model_prevalence"]
    elif response_track == "naive":
        output["model_weight"] = output["response_n_genomes_prn_interpretable"]
        output["model_successes"] = output["response_n_prn_disrupted"]
        output["model_prevalence"] = output["response_naive_prevalence"]
    elif response_track == "boundary_lower":
        output["model_weight"] = output["response_n_genomes_total"]
        output["model_successes"] = output["response_n_prn_disrupted"]
        output["model_prevalence"] = output["response_boundary_lower_prevalence"]
    elif response_track == "boundary_upper":
        output["model_weight"] = output["response_n_genomes_total"]
        output["model_successes"] = (
            output["response_n_prn_disrupted"].fillna(0) + output["response_n_missing_outcomes"].fillna(0)
        )
        output["model_prevalence"] = output["response_boundary_upper_prevalence"]
    else:
        raise ValueError(f"Unsupported response track: {response_track}")

    if response_track != "ipw":
        output["analysis_response"] = np.where(
            output["model_weight"].fillna(0).gt(0),
            output["model_successes"] / output["model_weight"],
            np.nan,
        )
    output["response_prevalence_cc"] = (
        output["model_successes"].fillna(0.0) + 0.5
    ) / (output["model_weight"] + 1.0)
    return output


def model_frame(panel: pd.DataFrame, *, exposure_column: str, response_track: str) -> pd.DataFrame:
    frame = panel.loc[panel["primary_panel_eligible"]].copy()
    frame = frame.loc[frame[exposure_column].fillna("").astype(str).str.len().gt(0)].copy()
    frame = attach_response_track(frame, response_track)
    frame = frame.loc[
        frame["model_prevalence"].notna()
        & frame["model_weight"].notna()
        & frame["model_weight"].gt(0)
        & frame["reported_cases_period"].notna()
        & frame["reported_cases_period"].gt(0)
        & frame["genomes_per_case_effective"].notna()
        & frame["genomes_per_case_effective"].gt(0)
    ].copy()
    frame["log1p_reported_cases_z"] = standardize(np.log1p(frame["reported_cases_period"]))
    frame["genomes_per_case_z"] = standardize(frame["genomes_per_case_effective"])
    frame["n_new_origins_detected_period_z"] = standardize(
        frame["n_new_origins_detected_period"].fillna(0),
        fill_missing_with_zero=True,
    )
    frame["n_active_origin_clades_period_z"] = standardize(
        frame["n_active_origin_clades_period"].fillna(0),
        fill_missing_with_zero=True,
    )
    frame["years_since_first_local_origin_z"] = standardize(
        frame["years_since_first_local_origin"],
        fill_missing_with_zero=True,
    )
    frame["years_since_first_prn_detection_z"] = standardize(
        frame["years_since_first_prn_detection"],
        fill_missing_with_zero=True,
    )
    frame["has_local_origin_by_period_end"] = frame["has_local_origin_by_period_end"].astype(int)
    frame["has_prn_detection_by_period_end"] = frame["has_prn_detection_by_period_end"].astype(int)
    return frame


def diagnostic_row(
    *,
    model_id: str,
    analysis_panel: str,
    sensitivity_label: str,
    model_spec: str,
    response_track: str,
    exposure_column: str,
    excluded_country_iso3: str,
    frame: pd.DataFrame,
    converged: bool,
    covariance_type: str,
    log_likelihood: float | None,
    notes: str,
) -> dict[str, object]:
    if frame.empty:
        year_window = "NA"
    else:
        year_window = f"{int(frame['period_start'].min())}-{int(frame['period_end'].max())}"
    return {
        "model_id": model_id,
        "analysis_panel": analysis_panel,
        "sensitivity_label": sensitivity_label,
        "model_spec": model_spec,
        "response_track": response_track,
        "exposure_column": exposure_column,
        "excluded_country_iso3": excluded_country_iso3,
        "n_rows": int(len(frame)),
        "n_countries": int(frame["country_iso3"].nunique()) if not frame.empty else 0,
        "year_window": year_window,
        "converged": converged,
        "covariance_type": covariance_type,
        "log_likelihood": log_likelihood,
        "notes": notes,
    }


def extract_coefficient_results(
    result,
    frame: pd.DataFrame,
    *,
    model_id: str,
    analysis_panel: str,
    sensitivity_label: str,
    model_spec: str,
    response_track: str,
    exposure_column: str,
    excluded_country_iso3: str,
    country_filter: str,
    model_family: str,
    notes: str,
    regularized: bool = False,
) -> list[dict[str, object]]:
    year_window = f"{int(frame['period_start'].min())}-{int(frame['period_end'].max())}"
    terms = list(result.params.index)
    p_values = []
    for term in terms:
        try:
            p_values.append(float(result.pvalues[term]))
        except Exception:
            p_values.append(np.nan)
    valid_p = [value for value in p_values if not np.isnan(value)]
    q_lookup = {}
    if valid_p:
        adjusted = benjamini_hochberg(valid_p)
        iterator = iter(adjusted)
        for term, p_value in zip(terms, p_values):
            q_lookup[term] = np.nan if np.isnan(p_value) else next(iterator)
    rows: list[dict[str, object]] = []
    for term, p_value in zip(terms, p_values):
        estimate = float(result.params[term])
        if regularized:
            lower = np.nan
            upper = np.nan
        else:
            standard_error = float(result.bse[term])
            lower = estimate - 1.96 * standard_error
            upper = estimate + 1.96 * standard_error
        note_text = (
            f"{notes};q_value_scope=within_model_term_family_bh_not_manuscript_wide_fdr"
            if notes
            else "q_value_scope=within_model_term_family_bh_not_manuscript_wide_fdr"
        )
        rows.append(
            {
                "model_id": model_id,
                "analysis_panel": analysis_panel,
                "sensitivity_label": sensitivity_label,
                "model_spec": model_spec,
                "response_track": response_track,
                "exposure_column": exposure_column,
                "excluded_country_iso3": excluded_country_iso3,
                "country_filter": country_filter,
                "n_rows": int(len(frame)),
                "n_countries": int(frame["country_iso3"].nunique()),
                "year_window": year_window,
                "model_family": model_family,
                "weighting_scheme": response_track_label(response_track),
                "result_type": "coefficient",
                "estimate_term": term,
                "program_class": "",
                "reference_class": "wp_only_or_pre_ap",
                "effect_scale": "log_odds",
                "effect_estimate": estimate,
                "ci_lower": lower,
                "ci_upper": upper,
                "p_value": p_value,
                "q_value": q_lookup.get(term, np.nan),
                "notes": note_text,
            }
        )
    return rows


def parameter_draws(result, n_draws: int = 400) -> np.ndarray | None:
    try:
        covariance = np.asarray(result.cov_params(), dtype=float)
        params = np.asarray(result.params, dtype=float)
    except Exception:
        return None


def fit_response_model(
    frame: pd.DataFrame,
    *,
    formula: str,
    response_track: str,
    regularized: bool = False,
):
    n_countries = int(frame["country_iso3"].nunique())
    if response_track == "ipw":
        model = smf.glm(
            formula=formula,
            data=frame,
            family=sm.families.Binomial(),
        )
        if regularized:
            fitted = model.fit_regularized(alpha=0.1, L1_wt=0.0, maxiter=200)
            covariance_type = "regularized_ridge"
        else:
            if n_countries >= 3:
                try:
                    fitted = model.fit(cov_type="cluster", cov_kwds={"groups": frame["country_iso3"]}, maxiter=200)
                    covariance_type = "cluster"
                except Exception:
                    fitted = model.fit(cov_type="HC1", maxiter=200)
                    covariance_type = "hc1_fallback_from_cluster"
            else:
                fitted = model.fit(cov_type="HC1", maxiter=200)
                covariance_type = "hc1"
        model_family = f"statsmodels_fractional_logit_{covariance_type}"
        return fitted, covariance_type, model_family

    model = smf.glm(
        formula=formula,
        data=frame,
        family=sm.families.Binomial(),
        freq_weights=frame["model_weight"],
    )
    if regularized:
        fitted = model.fit_regularized(alpha=0.1, L1_wt=0.0, maxiter=200)
        covariance_type = "regularized_ridge"
    else:
        if n_countries >= 3:
            try:
                fitted = model.fit(cov_type="cluster", cov_kwds={"groups": frame["country_iso3"]}, maxiter=200)
                covariance_type = "cluster"
            except Exception:
                fitted = model.fit(cov_type="HC1", maxiter=200)
                covariance_type = "hc1_fallback_from_cluster"
        else:
            fitted = model.fit(cov_type="HC1", maxiter=200)
            covariance_type = "hc1"
    model_family = f"statsmodels_glm_binomial_grouped_{covariance_type}"
    return fitted, covariance_type, model_family
    if covariance.shape[0] != covariance.shape[1] or covariance.shape[0] != params.shape[0]:
        return None
    covariance = covariance + np.eye(covariance.shape[0]) * 1e-10
    try:
        return np.random.default_rng(0).multivariate_normal(params, covariance, size=n_draws)
    except Exception:
        return None


def marginal_prediction_rows(
    result,
    frame: pd.DataFrame,
    *,
    model_id: str,
    analysis_panel: str,
    sensitivity_label: str,
    model_spec: str,
    response_track: str,
    exposure_column: str,
    excluded_country_iso3: str,
    country_filter: str,
    model_family: str,
    notes: str,
    regularized: bool = False,
) -> list[dict[str, object]]:
    present_classes = [value for value in CLASS_ORDER if value in set(frame[exposure_column])]
    if not present_classes:
        return []

    draws = None if regularized else parameter_draws(result)
    design_info = None if regularized else result.model.data.design_info
    year_window = f"{int(frame['period_start'].min())}-{int(frame['period_end'].max())}"
    rows: list[dict[str, object]] = []
    reference_class = "wp_only_or_pre_ap" if "wp_only_or_pre_ap" in present_classes else present_classes[0]
    reference_draws = None
    reference_point = None

    for program_class in present_classes:
        counterfactual = frame.copy()
        counterfactual[exposure_column] = program_class
        point_value = float(np.mean(result.predict(counterfactual)))

        lower = np.nan
        upper = np.nan
        class_draws = None
        if draws is not None and design_info is not None:
            try:
                design_matrix = build_design_matrices([design_info], counterfactual, return_type="dataframe")[0]
                linear_predictor = np.asarray(design_matrix, dtype=float) @ draws.T
                probabilities = 1.0 / (1.0 + np.exp(-linear_predictor))
                class_draws = probabilities.mean(axis=0)
                lower = float(np.quantile(class_draws, 0.025))
                upper = float(np.quantile(class_draws, 0.975))
            except Exception:
                class_draws = None

        if program_class == reference_class:
            reference_draws = class_draws
            reference_point = point_value

        rows.append(
            {
                "model_id": model_id,
                "analysis_panel": analysis_panel,
                "sensitivity_label": sensitivity_label,
                "model_spec": model_spec,
                "response_track": response_track,
                "exposure_column": exposure_column,
                "excluded_country_iso3": excluded_country_iso3,
                "country_filter": country_filter,
                "n_rows": int(len(frame)),
                "n_countries": int(frame["country_iso3"].nunique()),
                "year_window": year_window,
                "model_family": model_family,
                "weighting_scheme": response_track_label(response_track),
                "result_type": "adjusted_prevalence",
                "estimate_term": f"adjusted_prevalence::{program_class}",
                "program_class": program_class,
                "reference_class": reference_class,
                "effect_scale": "probability",
                "effect_estimate": point_value,
                "ci_lower": lower,
                "ci_upper": upper,
                "p_value": np.nan,
                "q_value": np.nan,
                "notes": notes,
            }
        )

        if program_class != reference_class:
            difference = np.nan if reference_point is None else point_value - reference_point
            diff_lower = np.nan
            diff_upper = np.nan
            if class_draws is not None and reference_draws is not None:
                diff_draws = class_draws - reference_draws
                diff_lower = float(np.quantile(diff_draws, 0.025))
                diff_upper = float(np.quantile(diff_draws, 0.975))
            rows.append(
                {
                    "model_id": model_id,
                    "analysis_panel": analysis_panel,
                    "sensitivity_label": sensitivity_label,
                    "model_spec": model_spec,
                    "response_track": response_track,
                    "exposure_column": exposure_column,
                    "excluded_country_iso3": excluded_country_iso3,
                    "country_filter": country_filter,
                    "n_rows": int(len(frame)),
                    "n_countries": int(frame["country_iso3"].nunique()),
                    "year_window": year_window,
                    "model_family": model_family,
                    "weighting_scheme": response_track_label(response_track),
                    "result_type": "adjusted_prevalence_difference",
                    "estimate_term": f"adjusted_prevalence_difference::{program_class}::{reference_class}",
                    "program_class": program_class,
                    "reference_class": reference_class,
                    "effect_scale": "probability_difference",
                    "effect_estimate": difference,
                    "ci_lower": diff_lower,
                    "ci_upper": diff_upper,
                    "p_value": np.nan,
                    "q_value": np.nan,
                    "notes": notes,
                }
            )
    return rows


def fit_single_model(
    results: list[dict[str, object]],
    diagnostics: list[dict[str, object]],
    frame: pd.DataFrame,
    *,
    model_id: str,
    analysis_panel: str,
    sensitivity_label: str,
    model_spec: str,
    response_track: str,
    exposure_column: str,
    country_filter: str,
    excluded_country_iso3: str = "",
    notes: str,
    regularized: bool = False,
) -> None:
    if frame.empty:
        diagnostics.append(
            diagnostic_row(
                model_id=model_id,
                analysis_panel=analysis_panel,
                sensitivity_label=sensitivity_label,
                model_spec=model_spec,
                response_track=response_track,
                exposure_column=exposure_column,
                excluded_country_iso3=excluded_country_iso3,
                frame=frame,
                converged=False,
                covariance_type="NA",
                log_likelihood=None,
                notes="empty_subset",
            )
        )
        return

    if "wp_only_or_pre_ap" not in set(frame[exposure_column]):
        diagnostics.append(
            diagnostic_row(
                model_id=model_id,
                analysis_panel=analysis_panel,
                sensitivity_label=sensitivity_label,
                model_spec=model_spec,
                response_track=response_track,
                exposure_column=exposure_column,
                excluded_country_iso3=excluded_country_iso3,
                frame=frame,
                converged=False,
                covariance_type="NA",
                log_likelihood=None,
                notes="missing_reference_class_wp_only_or_pre_ap",
            )
        )
        return

    include_bridge = model_spec == "programme_plus_bridge"
    formula = build_formula(exposure_column, include_bridge=include_bridge)
    covariance_type = "NA"
    try:
        fitted, covariance_type, model_family = fit_response_model(
            frame,
            formula=formula,
            response_track=response_track,
            regularized=regularized,
        )
        results.extend(
            extract_coefficient_results(
                fitted,
                frame,
                model_id=model_id,
                analysis_panel=analysis_panel,
                sensitivity_label=sensitivity_label,
                model_spec=model_spec,
                response_track=response_track,
                exposure_column=exposure_column,
                excluded_country_iso3=excluded_country_iso3,
                country_filter=country_filter,
                model_family=model_family,
                notes=notes,
                regularized=regularized,
            )
        )
        results.extend(
            marginal_prediction_rows(
                fitted,
                frame,
                model_id=model_id,
                analysis_panel=analysis_panel,
                sensitivity_label=sensitivity_label,
                model_spec=model_spec,
                response_track=response_track,
                exposure_column=exposure_column,
                excluded_country_iso3=excluded_country_iso3,
                country_filter=country_filter,
                model_family=model_family,
                notes=notes,
                regularized=regularized,
            )
        )
        log_likelihood = None if regularized else float(fitted.llf)
        diagnostics.append(
            diagnostic_row(
                model_id=model_id,
                analysis_panel=analysis_panel,
                sensitivity_label=sensitivity_label,
                model_spec=model_spec,
                response_track=response_track,
                exposure_column=exposure_column,
                excluded_country_iso3=excluded_country_iso3,
                frame=frame,
                converged=True,
                covariance_type=covariance_type,
                log_likelihood=log_likelihood,
                notes=notes,
            )
        )
    except Exception as error:  # pragma: no cover - diagnostic path
        tokens = structural_diagnostic_tokens(frame, exposure_column)
        failure_note = f"fit_failed:{error}"
        if tokens:
            failure_note = ";".join([failure_note, *tokens])
        diagnostics.append(
            diagnostic_row(
                model_id=model_id,
                analysis_panel=analysis_panel,
                sensitivity_label=sensitivity_label,
                model_spec=model_spec,
                response_track=response_track,
                exposure_column=exposure_column,
                excluded_country_iso3=excluded_country_iso3,
                frame=frame,
                converged=False,
                covariance_type=covariance_type,
                log_likelihood=None,
                notes=failure_note,
            )
        )


def run_single_model_task(task: dict[str, object]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    results: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []
    fit_single_model(
        results,
        diagnostics,
        task["frame"],
        model_id=str(task["model_id"]),
        analysis_panel=str(task["analysis_panel"]),
        sensitivity_label=str(task["sensitivity_label"]),
        model_spec=str(task["model_spec"]),
        response_track=str(task["response_track"]),
        exposure_column=str(task["exposure_column"]),
        country_filter=str(task["country_filter"]),
        excluded_country_iso3=str(task.get("excluded_country_iso3", "")),
        notes=str(task["notes"]),
        regularized=bool(task.get("regularized", False)),
    )
    return results, diagnostics


def run_country_leave_one_out(
    results: list[dict[str, object]],
    diagnostics: list[dict[str, object]],
    frame: pd.DataFrame,
    *,
    analysis_panel: str,
    notes: str,
    max_workers: int | None = None,
) -> None:
    tasks: list[dict[str, object]] = []
    for excluded_country in sorted(frame["country_iso3"].dropna().astype(str).unique()):
        subset = frame.loc[frame["country_iso3"].ne(excluded_country)].copy()
        tasks.append(
            {
                "frame": subset,
                "model_id": f"{analysis_panel}_leave_one_out_{excluded_country}",
                "analysis_panel": analysis_panel,
                "sensitivity_label": "leave_one_country_out",
                "model_spec": "programme_only",
                "response_track": "ipw",
                "exposure_column": "program_formulation_class_concurrent",
                "country_filter": f"exclude_country={excluded_country}",
                "excluded_country_iso3": excluded_country,
                "notes": notes + f";leave_one_country_out={excluded_country}",
            }
        )

    worker_count = resolve_parallel_workers(len(tasks), requested_max_workers=max_workers)
    if worker_count <= 1:
        task_outputs = [run_single_model_task(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=worker_count) as pool:
            task_outputs = list(pool.map(run_single_model_task, tasks))

    for task_results, task_diagnostics in task_outputs:
        results.extend(task_results)
        diagnostics.extend(task_diagnostics)


def run_all_models(
    primary_panel: pd.DataFrame,
    secondary_panel: pd.DataFrame | None,
    *,
    max_workers: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    results: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []

    primary_label = (
        normalize_text(primary_panel["analysis_panel"].iloc[0])
        if not primary_panel.empty
        else "programme_country_period_primary"
    )
    shared_notes = (
        "programme_surveillance_design_upgrade;"
        "headline_focus=surveillance_heterogeneity_not_causal_policy_effect;"
        "sampling_covariates_are_diagnostics_not_policy_mechanisms"
    )

    primary_ipw_concurrent = model_frame(
        primary_panel,
        exposure_column="program_formulation_class_concurrent",
        response_track="ipw",
    )
    fit_single_model(
        results,
        diagnostics,
        primary_ipw_concurrent,
        model_id=f"{primary_label}_primary_programme_only_ipw",
        analysis_panel=primary_label,
        sensitivity_label="primary_period_panel",
        model_spec="programme_only",
        response_track="ipw",
        exposure_column="program_formulation_class_concurrent",
        country_filter="all_primary_eligible_period_rows",
        notes=shared_notes + ";response=ipw;exposure=concurrent_programme_class",
    )
    fit_single_model(
        results,
        diagnostics,
        primary_ipw_concurrent,
        model_id=f"{primary_label}_primary_programme_plus_bridge_ipw",
        analysis_panel=primary_label,
        sensitivity_label="primary_period_panel",
        model_spec="programme_plus_bridge",
        response_track="ipw",
        exposure_column="program_formulation_class_concurrent",
        country_filter="all_primary_eligible_period_rows",
        notes=shared_notes + ";response=ipw;exposure=concurrent_programme_class;bridge=origin_and_detection_context",
    )

    for response_track in ["naive", "boundary_lower", "boundary_upper"]:
        fit_single_model(
            results,
            diagnostics,
            model_frame(primary_panel, exposure_column="program_formulation_class_concurrent", response_track=response_track),
            model_id=f"{primary_label}_primary_programme_only_{response_track}",
            analysis_panel=primary_label,
            sensitivity_label="response_track_sensitivity",
            model_spec="programme_only",
            response_track=response_track,
            exposure_column="program_formulation_class_concurrent",
            country_filter="all_primary_eligible_period_rows",
            notes=shared_notes + f";response={response_track};exposure=concurrent_programme_class",
        )

    if "USA" in set(primary_ipw_concurrent["country_iso3"]):
        fit_single_model(
            results,
            diagnostics,
            primary_ipw_concurrent.loc[primary_ipw_concurrent["country_iso3"].ne("USA")].copy(),
            model_id=f"{primary_label}_exclude_USA_programme_only_ipw",
            analysis_panel=primary_label,
            sensitivity_label="exclude_usa",
            model_spec="programme_only",
            response_track="ipw",
            exposure_column="program_formulation_class_concurrent",
            country_filter="exclude_country=USA",
            excluded_country_iso3="USA",
            notes=shared_notes + ";response=ipw;exposure=concurrent_programme_class;country_dependence_check=USA_excluded",
        )

    high_confidence_panel = primary_panel.loc[
        primary_panel["formulation_confidence_period"].astype(str).str.lower().eq("high")
    ].copy()
    fit_single_model(
        results,
        diagnostics,
        model_frame(high_confidence_panel, exposure_column="program_formulation_class_concurrent", response_track="ipw"),
        model_id=f"{primary_label}_high_confidence_programme_only_ipw",
        analysis_panel=primary_label,
        sensitivity_label="high_confidence_only",
        model_spec="programme_only",
        response_track="ipw",
        exposure_column="program_formulation_class_concurrent",
        country_filter="formulation_confidence_period=high",
        notes=shared_notes + ";response=ipw;row_subset=high_confidence_only",
    )

    transition_excluded = primary_panel.loc[~primary_panel["transition_period_flag"].fillna(False)].copy()
    fit_single_model(
        results,
        diagnostics,
        model_frame(transition_excluded, exposure_column="program_formulation_class_concurrent", response_track="ipw"),
        model_id=f"{primary_label}_transition_excluded_programme_only_ipw",
        analysis_panel=primary_label,
        sensitivity_label="transition_excluded",
        model_spec="programme_only",
        response_track="ipw",
        exposure_column="program_formulation_class_concurrent",
        country_filter="transition_period_flag=false",
        notes=shared_notes + ";response=ipw;transition_rows_excluded=true",
    )

    lagged_subset = primary_panel.loc[primary_panel["lagged_class_available"].fillna(False)].copy()
    fit_single_model(
        results,
        diagnostics,
        model_frame(lagged_subset, exposure_column="program_formulation_class_lagged", response_track="ipw"),
        model_id=f"{primary_label}_lagged_programme_only_ipw",
        analysis_panel=primary_label,
        sensitivity_label="lagged_exposure",
        model_spec="programme_only",
        response_track="ipw",
        exposure_column="program_formulation_class_lagged",
        country_filter="lagged_class_available=true",
        notes=shared_notes + ";response=ipw;exposure=lagged_programme_class",
    )

    fit_single_model(
        results,
        diagnostics,
        primary_ipw_concurrent,
        model_id=f"{primary_label}_ridge_penalized_programme_only_ipw",
        analysis_panel=primary_label,
        sensitivity_label="ridge_penalized_support",
        model_spec="programme_only",
        response_track="ipw",
        exposure_column="program_formulation_class_concurrent",
        country_filter="all_primary_eligible_period_rows",
        notes=shared_notes + ";response=ipw;penalized_support=true",
        regularized=True,
    )

    run_country_leave_one_out(
        results,
        diagnostics,
        primary_ipw_concurrent,
        analysis_panel=primary_label,
        notes=shared_notes + ";response=ipw;exposure=concurrent_programme_class",
        max_workers=max_workers,
    )

    if secondary_panel is not None and not secondary_panel.empty:
        secondary_label = normalize_text(secondary_panel["analysis_panel"].iloc[0]) or "programme_country_period_secondary"
        fit_single_model(
            results,
            diagnostics,
            model_frame(secondary_panel, exposure_column="program_formulation_class_concurrent", response_track="ipw"),
            model_id=f"{secondary_label}_supporting_programme_only_ipw",
            analysis_panel=secondary_label,
            sensitivity_label="supporting_period_panel",
            model_spec="programme_only",
            response_track="ipw",
            exposure_column="program_formulation_class_concurrent",
            country_filter="all_primary_eligible_period_rows",
            notes=shared_notes + ";response=ipw;supporting_panel=yes",
        )

    return (
        pd.DataFrame(results, columns=RESULT_COLUMNS),
        pd.DataFrame(diagnostics, columns=DIAGNOSTIC_COLUMNS),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Fit programme-surveillance models.")
    parser.add_argument("--primary-panel", required=True, help="Primary country-period panel TSV")
    parser.add_argument("--secondary-panel", default="", help="Optional supporting panel TSV")
    parser.add_argument("--results-out", required=True, help="Model results TSV")
    parser.add_argument("--diagnostics-out", required=True, help="Model diagnostics TSV")
    parser.add_argument(
        "--max-workers",
        type=int,
        default=None,
        help=(
            "Maximum process workers for leave-one-country-out fits. "
            f"Defaults to the {PROGRAMME_MODELS_MAX_WORKERS_ENV} env var or an auto cap of "
            f"{PROGRAMME_MODELS_DEFAULT_MAX_WORKERS}."
        ),
    )
    args = parser.parse_args()

    primary_panel = prepare_panel(args.primary_panel)
    secondary_panel = prepare_panel(args.secondary_panel) if normalize_text(args.secondary_panel) else None
    requested_max_workers = args.max_workers
    if requested_max_workers is None:
        requested_max_workers = resolve_env_max_workers()

    results_frame, diagnostics_frame = run_all_models(
        primary_panel,
        secondary_panel,
        max_workers=requested_max_workers,
    )

    results_path = Path(args.results_out)
    diagnostics_path = Path(args.diagnostics_out)
    results_path.parent.mkdir(parents=True, exist_ok=True)
    diagnostics_path.parent.mkdir(parents=True, exist_ok=True)
    results_frame.to_csv(results_path, sep="\t", index=False)
    diagnostics_frame.to_csv(diagnostics_path, sep="\t", index=False)
    print(f"Wrote programme-surveillance model results: {results_path}")
    print(f"Wrote programme-surveillance model diagnostics: {diagnostics_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
