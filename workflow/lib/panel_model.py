#!/usr/bin/env python3
"""Fit workflow-native country-level association models for manuscript ecology analyses."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm


PANEL_MAX_WORKERS_ENV = "PANEL_MAX_WORKERS"
PANEL_DEFAULT_MAX_WORKERS = 32


RESULT_COLUMNS = [
    "model_id",
    "analysis_cohort",
    "response_variable",
    "model_family",
    "country_filter",
    "year_window",
    "n_country_year_cells",
    "n_countries",
    "covariates",
    "random_effects",
    "weighting_scheme",
    "estimate_term",
    "effect_scale",
    "effect_estimate",
    "ci_lower",
    "ci_upper",
    "p_value",
    "q_value",
    "sensitivity_label",
    "focal_exposure_family",
    "excluded_country_iso3",
    "notes",
    "exposure_formula_id",
    "exposure_lambda",
    "exposure_gamma",
    "exposure_delta_prn",
    "uses_ipw_response",
]

DIAGNOSTIC_COLUMNS = [
    "model_id",
    "sensitivity_label",
    "focal_exposure_family",
    "excluded_country_iso3",
    "n_country_year_cells",
    "n_countries",
    "year_window",
    "converged",
    "log_likelihood",
    "notes",
]


def resolve_parallel_workers(
    task_count: int | None = None,
    *,
    requested_max_workers: int | None = None,
    cpu_count: int | None = None,
) -> int:
    available_cpu = max(1, cpu_count or os.cpu_count() or 1)
    workers = requested_max_workers if requested_max_workers is not None else min(available_cpu, PANEL_DEFAULT_MAX_WORKERS)
    workers = max(1, min(int(workers), available_cpu))
    if task_count is not None:
        workers = min(workers, max(1, int(task_count)))
    return workers


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


def standardize(series: pd.Series) -> pd.Series:
    valid = series.dropna()
    if valid.empty:
        return pd.Series(np.nan, index=series.index)
    deviation = valid.std(ddof=0)
    if deviation == 0 or np.isnan(deviation):
        return pd.Series(0.0, index=series.index)
    return (series - valid.mean()) / deviation


def prepare_panel_dataset(exposure_path: str, prevalence_path: str) -> pd.DataFrame:
    exposure = pd.read_csv(exposure_path, sep="\t", dtype=str)
    prevalence = pd.read_csv(prevalence_path, sep="\t", dtype=str)

    for frame in [exposure, prevalence]:
        frame["year"] = coerce_numeric(frame.get("year", pd.Series(dtype=str)))
        frame["country_iso3"] = frame.get("country_iso3", pd.Series(dtype=str)).fillna("")

    exposure["ap_exposure_v1_score"] = coerce_numeric(exposure.get("ap_exposure_v1_score", pd.Series(dtype=str)))
    exposure["ap_exposure_v2_score"] = coerce_numeric(exposure.get("ap_exposure_v2_score", pd.Series(dtype=str)))
    exposure["ap_exposure_v3_score"] = coerce_numeric(exposure.get("ap_exposure_v3_score", pd.Series(dtype=str)))
    exposure["dtp3_only_score"] = coerce_numeric(exposure.get("dtp3_only_score", pd.Series(dtype=str)))
    exposure["dtp3_coverage"] = coerce_numeric(exposure.get("dtp3_coverage", pd.Series(dtype=str)))
    exposure["reported_cases"] = coerce_numeric(exposure.get("reported_cases", pd.Series(dtype=str))).fillna(0)
    exposure["genomes_per_case"] = coerce_numeric(exposure.get("genomes_per_case", pd.Series(dtype=str)))
    exposure["post_covid_period"] = coerce_numeric(exposure.get("post_covid_period", pd.Series(dtype=str))).fillna(0)
    exposure["exposure_lambda_years"] = coerce_numeric(exposure.get("exposure_lambda_years", pd.Series(dtype=str)))
    exposure["exposure_gamma_booster"] = coerce_numeric(exposure.get("exposure_gamma_booster", pd.Series(dtype=str)))
    exposure["exposure_delta_prn"] = coerce_numeric(exposure.get("exposure_delta_prn", pd.Series(dtype=str)))
    exposure["ap_exposure_v2_available"] = (
        exposure.get("ap_exposure_v2_available", pd.Series(dtype=str))
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .isin(["true", "1", "yes"])
    )
    exposure["ap_exposure_v3_available"] = (
        exposure.get("ap_exposure_v3_available", pd.Series(dtype=str))
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        .isin(["true", "1", "yes"])
    )

    prevalence["n_genomes_prn_interpretable"] = coerce_numeric(prevalence.get("n_genomes_prn_interpretable", pd.Series(dtype=str))).fillna(0)
    prevalence["n_genomes_total"] = coerce_numeric(prevalence.get("n_genomes_total", pd.Series(dtype=str))).fillna(0)
    prevalence["ipw_weight_total"] = coerce_numeric(prevalence.get("ipw_weight_total", pd.Series(dtype=str)))
    prevalence["ipw_prevalence"] = coerce_numeric(prevalence.get("ipw_prevalence", pd.Series(dtype=str)))
    prevalence["naive_prevalence"] = coerce_numeric(prevalence.get("naive_prevalence", pd.Series(dtype=str)))

    prevalence_subset = prevalence[
        [
            "country_iso3",
            "year",
            "n_genomes_total",
            "n_genomes_prn_interpretable",
            "ipw_weight_total",
            "ipw_prevalence",
            "naive_prevalence",
        ]
    ].rename(
        columns={
            "n_genomes_total": "response_n_genomes_total",
            "n_genomes_prn_interpretable": "response_n_genomes_prn_interpretable",
            "ipw_weight_total": "response_ipw_weight_total",
            "ipw_prevalence": "response_ipw_prevalence",
            "naive_prevalence": "response_naive_prevalence",
        }
    )

    dataset = exposure.merge(prevalence_subset, on=["country_iso3", "year"], how="inner")
    dataset = dataset.loc[dataset["country_iso3"].ne("") & dataset["year"].notna()].copy()
    dataset["year"] = dataset["year"].astype(int)
    dataset["workflow_genomes_per_case"] = np.where(
        dataset["reported_cases"].fillna(0).gt(0),
        dataset["response_n_genomes_total"] / dataset["reported_cases"],
        np.nan,
    )
    dataset["genomes_per_case_effective"] = dataset["genomes_per_case"].fillna(dataset["workflow_genomes_per_case"])
    return dataset


def build_design_matrix(frame: pd.DataFrame, exposure_column: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    design = pd.DataFrame(index=frame.index)
    design["Intercept"] = 1.0
    design[exposure_column] = standardize(frame[exposure_column])
    design["log1p_reported_cases_z"] = standardize(np.log1p(frame["reported_cases"]))
    design["post_covid_period"] = frame["post_covid_period"].astype(float)
    design["genomes_per_case_z"] = standardize(frame["genomes_per_case_effective"])

    response = pd.DataFrame(index=frame.index)
    response["prevalence"] = frame["response_ipw_prevalence"]
    response["weights"] = frame["response_ipw_weight_total"]
    return design, response


def fit_glm(frame: pd.DataFrame, exposure_column: str, cluster_robust: bool = False):
    design, response = build_design_matrix(frame, exposure_column)
    model = sm.GLM(
        response["prevalence"],
        design,
        family=sm.families.Binomial(),
    )
    if cluster_robust or frame["country_iso3"].nunique() >= 3:
        return model.fit(cov_type="cluster", cov_kwds={"groups": frame["country_iso3"]}, maxiter=200)
    return model.fit(cov_type="HC1", maxiter=200)


def extract_result_rows(
    result,
    frame: pd.DataFrame,
    *,
    model_id: str,
    sensitivity_label: str,
    focal_exposure_family: str,
    excluded_country_iso3: str,
    covariates: str,
    country_filter: str,
    notes: str,
    exposure_formula_id: str,
    exposure_lambda: float | None,
    exposure_gamma: float | None,
    exposure_delta_prn: float | None,
    exposure_term: str,
    model_family: str,
) -> list[dict[str, object]]:
    terms_to_report = [exposure_term, "log1p_reported_cases_z", "post_covid_period", "genomes_per_case_z"]
    p_values = [float(result.pvalues[term]) for term in terms_to_report]
    q_values = benjamini_hochberg(p_values)
    year_window = f"{int(frame['year'].min())}-{int(frame['year'].max())}"

    rows: list[dict[str, object]] = []
    for term, q_value in zip(terms_to_report, q_values):
        coefficient = float(result.params[term])
        standard_error = float(result.bse[term])
        note_text = (
            f"{notes};q_value_scope=within_model_reported_terms_bh_not_manuscript_wide_fdr"
            if notes
            else "q_value_scope=within_model_reported_terms_bh_not_manuscript_wide_fdr"
        )
        rows.append(
            {
                "model_id": model_id,
                "analysis_cohort": "C_IPW",
                "response_variable": "ipw_prn_disrupted_prevalence",
                "model_family": model_family,
                "country_filter": country_filter,
                "year_window": year_window,
                "n_country_year_cells": int(len(frame)),
                "n_countries": int(frame["country_iso3"].nunique()),
                "covariates": covariates,
                "random_effects": "not_fit",
                "weighting_scheme": "ipw_fractional_response_no_pseudotrials",
                "estimate_term": term,
                "effect_scale": "log_odds",
                "effect_estimate": coefficient,
                "ci_lower": coefficient - 1.96 * standard_error,
                "ci_upper": coefficient + 1.96 * standard_error,
                "p_value": float(result.pvalues[term]),
                "q_value": q_value,
                "sensitivity_label": sensitivity_label,
                "focal_exposure_family": focal_exposure_family,
                "excluded_country_iso3": excluded_country_iso3,
                "notes": note_text,
                "exposure_formula_id": exposure_formula_id,
                "exposure_lambda": exposure_lambda,
                "exposure_gamma": exposure_gamma,
                "exposure_delta_prn": exposure_delta_prn,
                "uses_ipw_response": True,
            }
        )
    return rows


def diagnostic_row(
    model_id: str,
    sensitivity_label: str,
    focal_exposure_family: str,
    excluded_country_iso3: str,
    frame: pd.DataFrame,
    converged: bool,
    log_likelihood: float | None,
    notes: str,
) -> dict[str, object]:
    if frame.empty:
        year_window = "NA"
    else:
        year_window = f"{int(frame['year'].min())}-{int(frame['year'].max())}"
    return {
        "model_id": model_id,
        "sensitivity_label": sensitivity_label,
        "focal_exposure_family": focal_exposure_family,
        "excluded_country_iso3": excluded_country_iso3,
        "n_country_year_cells": int(len(frame)),
        "n_countries": int(frame["country_iso3"].nunique()) if not frame.empty else 0,
        "year_window": year_window,
        "converged": converged,
        "log_likelihood": log_likelihood,
        "notes": notes,
    }


def run_single_model(
    results: list[dict[str, object]],
    diagnostics: list[dict[str, object]],
    frame: pd.DataFrame,
    *,
    exposure_column: str,
    model_id: str,
    sensitivity_label: str,
    focal_exposure_family: str,
    excluded_country_iso3: str = "",
    cluster_robust: bool = False,
    country_filter: str,
    notes: str,
) -> None:
    if frame.empty:
        diagnostics.append(
            diagnostic_row(
                model_id,
                sensitivity_label,
                focal_exposure_family,
                excluded_country_iso3,
                frame,
                False,
                None,
                "empty_subset",
            )
        )
        return

    exposure_formula_id = str(frame.get("exposure_formula_id", pd.Series([""])).iloc[0]) if "exposure_formula_id" in frame else ""
    exposure_lambda = float(frame["exposure_lambda_years"].iloc[0]) if "exposure_lambda_years" in frame and pd.notna(frame["exposure_lambda_years"].iloc[0]) else None
    exposure_gamma = float(frame["exposure_gamma_booster"].iloc[0]) if "exposure_gamma_booster" in frame and pd.notna(frame["exposure_gamma_booster"].iloc[0]) else None
    exposure_delta_prn = float(frame["exposure_delta_prn"].iloc[0]) if "exposure_delta_prn" in frame and pd.notna(frame["exposure_delta_prn"].iloc[0]) else None

    covariates = f"{exposure_column}_z,log1p_reported_cases_z,post_covid_period,genomes_per_case_z"
    effective_cluster = cluster_robust or frame["country_iso3"].nunique() >= 3
    model_family = "statsmodels_fractional_logit_cluster" if effective_cluster else "statsmodels_fractional_logit_hc1"

    try:
        result = fit_glm(frame, exposure_column, cluster_robust=cluster_robust)
        results.extend(
            extract_result_rows(
                result,
                frame,
                model_id=model_id,
                sensitivity_label=sensitivity_label,
                focal_exposure_family=focal_exposure_family,
                excluded_country_iso3=excluded_country_iso3,
                covariates=covariates,
                country_filter=country_filter,
                notes=notes,
                exposure_formula_id=exposure_formula_id,
                exposure_lambda=exposure_lambda,
                exposure_gamma=exposure_gamma,
                exposure_delta_prn=exposure_delta_prn,
                exposure_term=exposure_column,
                model_family=model_family,
            )
        )
        diagnostics.append(
            diagnostic_row(
                model_id,
                sensitivity_label,
                focal_exposure_family,
                excluded_country_iso3,
                frame,
                True,
                float(result.llf),
                notes,
            )
        )
    except Exception as error:  # pragma: no cover - diagnostic path
        diagnostics.append(
            diagnostic_row(
                model_id,
                sensitivity_label,
                focal_exposure_family,
                excluded_country_iso3,
                frame,
                False,
                None,
                f"fit_failed:{error}",
            )
        )


def run_single_model_task(task: dict[str, object]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    results: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []
    run_single_model(
        results,
        diagnostics,
        task["frame"],
        exposure_column=str(task["exposure_column"]),
        model_id=str(task["model_id"]),
        sensitivity_label=str(task["sensitivity_label"]),
        focal_exposure_family=str(task["focal_exposure_family"]),
        excluded_country_iso3=str(task.get("excluded_country_iso3", "")),
        cluster_robust=bool(task.get("cluster_robust", False)),
        country_filter=str(task["country_filter"]),
        notes=str(task["notes"]),
    )
    return results, diagnostics


def write_leave_one_out_summary(results_frame: pd.DataFrame, output_dir: Path) -> Path:
    focal_rows = results_frame.loc[
        results_frame["sensitivity_label"].isin(
            [
                "leave_one_country_out_ap_exposure_v2",
                "leave_one_country_out_ap_exposure_v3",
                "leave_one_country_out_ap_exposure_v1",
                "leave_one_country_out_legacy_dtp3",
            ]
        )
        & (
            ((results_frame["focal_exposure_family"] == "v3") & (results_frame["estimate_term"] == "ap_exposure_v3_score"))
            | ((results_frame["focal_exposure_family"] == "v2") & (results_frame["estimate_term"] == "ap_exposure_v2_score"))
            | ((results_frame["focal_exposure_family"] == "v1") & (results_frame["estimate_term"] == "ap_exposure_v1_score"))
            | ((results_frame["focal_exposure_family"] == "dtp3") & (results_frame["estimate_term"] == "dtp3_coverage"))
        )
    ].copy()

    primary_lookup = (
        results_frame.loc[
            (
                ((results_frame["focal_exposure_family"] == "v3") & (results_frame["sensitivity_label"] == "primary_ap_exposure_v3"))
                | ((results_frame["focal_exposure_family"] == "v2") & (results_frame["sensitivity_label"] == "primary_ap_exposure_v2"))
                | ((results_frame["focal_exposure_family"] == "v1") & (results_frame["sensitivity_label"] == "primary_ap_exposure_v1"))
                | ((results_frame["focal_exposure_family"] == "dtp3") & (results_frame["sensitivity_label"] == "legacy_dtp3_proxy"))
            )
            & (
                ((results_frame["focal_exposure_family"] == "v3") & (results_frame["estimate_term"] == "ap_exposure_v3_score"))
                | ((results_frame["focal_exposure_family"] == "v2") & (results_frame["estimate_term"] == "ap_exposure_v2_score"))
                | ((results_frame["focal_exposure_family"] == "v1") & (results_frame["estimate_term"] == "ap_exposure_v1_score"))
                | ((results_frame["focal_exposure_family"] == "dtp3") & (results_frame["estimate_term"] == "dtp3_coverage"))
            ),
            ["focal_exposure_family", "effect_estimate"],
        ]
        .drop_duplicates(subset=["focal_exposure_family"])
        .set_index("focal_exposure_family")["effect_estimate"]
        .to_dict()
    )
    focal_rows["primary_effect_estimate"] = focal_rows["focal_exposure_family"].map(primary_lookup)
    focal_rows["same_direction_as_primary"] = np.where(
        focal_rows["primary_effect_estimate"].notna(),
        np.sign(focal_rows["effect_estimate"]) == np.sign(focal_rows["primary_effect_estimate"]),
        np.nan,
    )
    leave_one_out_path = output_dir / "panel_model_leave_one_country_out.tsv"
    focal_rows.to_csv(leave_one_out_path, sep="\t", index=False)
    return leave_one_out_path


def write_coverage_report(base_common: pd.DataFrame, output_dir: Path) -> Path:
    n_total = len(base_common)
    report = pd.DataFrame(
        [
            {
                "metric": "panel_rows_primary_parameterization",
                "value": n_total,
                "notes": "rows with IPW outcome, >=5 interpretable genomes, reported cases, and effective genomes_per_case",
            },
            {
                "metric": "panel_countries_primary_parameterization",
                "value": int(base_common["country_iso3"].nunique()) if n_total else 0,
                "notes": "unique countries in the primary manuscript ecology panel",
            },
            {
                "metric": "v2_rows_with_known_prn",
                "value": int(base_common["ap_exposure_v2_available"].sum()) if n_total else 0,
                "notes": "primary-panel rows with curated yes/mixed/no Prn formulation status",
            },
            {
                "metric": "v3_rows_with_product_metadata",
                "value": int(base_common["ap_exposure_v3_available"].sum()) if n_total else 0,
                "notes": "primary-panel rows with non-missing product-aware V3 exposure inputs",
            },
            {
                "metric": "v2_known_fraction",
                "value": "" if n_total == 0 else f"{base_common['ap_exposure_v2_available'].mean():.6f}",
                "notes": "fraction of primary-panel rows eligible for V2",
            },
            {
                "metric": "v3_known_fraction",
                "value": "" if n_total == 0 else f"{base_common['ap_exposure_v3_available'].mean():.6f}",
                "notes": "fraction of primary-panel rows eligible for V3",
            },
        ]
    )
    coverage_path = output_dir / "panel_model_coverage_report.tsv"
    report.to_csv(coverage_path, sep="\t", index=False)
    return coverage_path


def run_models(
    dataset: pd.DataFrame,
    results_out: str,
    diagnostics_out: str,
    *,
    max_workers: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    results: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []

    primary_dataset = dataset.loc[
        dataset["is_primary_parameterization"].astype(str).str.lower().eq("true")
    ].copy()
    base_common = primary_dataset.loc[
        primary_dataset["reported_cases"].notna()
        & primary_dataset["genomes_per_case_effective"].notna()
        & primary_dataset["response_ipw_prevalence"].notna()
        & primary_dataset["response_ipw_weight_total"].gt(0)
        & primary_dataset["response_n_genomes_prn_interpretable"].ge(5)
    ].copy()

    base_v3 = base_common.loc[base_common["ap_exposure_v3_score"].notna()].copy()
    base_v2 = base_common.loc[base_common["ap_exposure_v2_score"].notna()].copy()
    base_v1 = base_common.loc[base_common["ap_exposure_v1_score"].notna()].copy()
    base_dtp3 = base_common.loc[base_common["dtp3_coverage"].notna()].copy()
    v3_known_fraction = float(base_common["ap_exposure_v3_available"].mean()) if not base_common.empty else float("nan")
    v2_known_fraction = float(base_common["ap_exposure_v2_available"].mean()) if not base_common.empty else float("nan")

    v3_notes = (
        "workflow_native_m7_primary_v3;response_uses_ipw_prevalence;"
        "aPExposure_V3_uses_DTP3_plus_curated_ap_timing_plus_role_specific_product_metadata_on_routine_primary_and_booster_exposure;"
        f"v3_known_fraction={v3_known_fraction:.3f}"
    )
    v2_notes = (
        "workflow_native_m7_primary_v2;response_uses_ipw_prevalence;"
        "aPExposure_V2_uses_DTP3_plus_curated_ap_timing_plus_booster_flag_plus_PRN_formulation_status;"
        f"v2_known_fraction={v2_known_fraction:.3f}"
    )
    v1_notes = (
        "workflow_native_m7_primary_v1;response_uses_ipw_prevalence;"
        "aPExposure_V1_uses_DTP3_plus_years_since_first_routine_ap_intro_plus_booster_flag"
    )
    dtp3_notes = "legacy_proxy_model_for_comparison_only;DTP3_not_interpreted_as_mechanistic_selection_pressure"

    run_single_model(
        results,
        diagnostics,
        base_v3,
        exposure_column="ap_exposure_v3_score",
        model_id="m7_primary_ipw_ap_exposure_v3_glm",
        sensitivity_label="primary_ap_exposure_v3",
        focal_exposure_family="v3",
        country_filter="all_available_country_year_cells_with_complete_v3_covariates",
        notes=v3_notes,
    )
    run_single_model(
        results,
        diagnostics,
        base_v2,
        exposure_column="ap_exposure_v2_score",
        model_id="m7_primary_ipw_ap_exposure_v2_glm",
        sensitivity_label="primary_ap_exposure_v2",
        focal_exposure_family="v2",
        country_filter="all_available_country_year_cells_with_complete_v2_covariates",
        notes=v2_notes,
    )
    run_single_model(
        results,
        diagnostics,
        base_v1,
        exposure_column="ap_exposure_v1_score",
        model_id="m7_primary_ipw_ap_exposure_v1_glm",
        sensitivity_label="primary_ap_exposure_v1",
        focal_exposure_family="v1",
        country_filter="all_available_country_year_cells_with_complete_v1_covariates",
        notes=v1_notes,
    )
    run_single_model(
        results,
        diagnostics,
        base_dtp3,
        exposure_column="dtp3_coverage",
        model_id="m7_legacy_dtp3_proxy_glm",
        sensitivity_label="legacy_dtp3_proxy",
        focal_exposure_family="dtp3",
        country_filter="same_country_year_cells_as_primary_when_dtp3_present",
        notes=dtp3_notes,
    )
    if base_v3["country_iso3"].nunique() >= 3:
        run_single_model(
            results,
            diagnostics,
            base_v3,
            exposure_column="ap_exposure_v3_score",
            model_id="m7_primary_ipw_ap_exposure_v3_cluster_robust",
            sensitivity_label="cluster_robust_ap_exposure_v3",
            focal_exposure_family="v3",
            cluster_robust=True,
            country_filter="all_available_country_year_cells_with_complete_v3_covariates",
            notes=v3_notes + ";cluster_robust_se_by_country",
        )
    if base_v2["country_iso3"].nunique() >= 3:
        run_single_model(
            results,
            diagnostics,
            base_v2,
            exposure_column="ap_exposure_v2_score",
            model_id="m7_primary_ipw_ap_exposure_v2_cluster_robust",
            sensitivity_label="cluster_robust_ap_exposure_v2",
            focal_exposure_family="v2",
            cluster_robust=True,
            country_filter="all_available_country_year_cells_with_complete_v2_covariates",
            notes=v2_notes + ";cluster_robust_se_by_country",
        )

    special_sensitivity_specs = [
        (
            "v3",
            "ap_exposure_v3_score",
            base_v3,
            "exclude_usa_ap_exposure_v3",
            "m7_exclude_usa_ap_exposure_v3_glm",
            base_v3.loc[base_v3["country_iso3"].ne("USA")].copy(),
            "exclude_country=USA",
            v3_notes + ";explicit_exclusion=USA",
        ),
        (
            "v2",
            "ap_exposure_v2_score",
            base_v2,
            "exclude_usa_ap_exposure_v2",
            "m7_exclude_usa_ap_exposure_v2_glm",
            base_v2.loc[base_v2["country_iso3"].ne("USA")].copy(),
            "exclude_country=USA",
            v2_notes + ";explicit_exclusion=USA",
        ),
        (
            "v3",
            "ap_exposure_v3_score",
            base_v3,
            "exclude_post2020_china_ap_exposure_v3",
            "m7_exclude_post2020_china_ap_exposure_v3_glm",
            base_v3.loc[~(base_v3["country_iso3"].eq("CHN") & base_v3["year"].ge(2021))].copy(),
            "exclude_country_year_subset=CHN_year_ge_2021",
            v3_notes + ";explicit_exclusion=proxy_for_mr_dominant_post2020_china_rows",
        ),
        (
            "v1",
            "ap_exposure_v1_score",
            base_v1,
            "exclude_usa_ap_exposure_v1",
            "m7_exclude_usa_ap_exposure_v1_glm",
            base_v1.loc[base_v1["country_iso3"].ne("USA")].copy(),
            "exclude_country=USA",
            v1_notes + ";explicit_exclusion=USA",
        ),
        (
            "dtp3",
            "dtp3_coverage",
            base_dtp3,
            "exclude_usa_legacy_dtp3",
            "m7_exclude_usa_legacy_dtp3_glm",
            base_dtp3.loc[base_dtp3["country_iso3"].ne("USA")].copy(),
            "exclude_country=USA",
            dtp3_notes + ";explicit_exclusion=USA",
        ),
        (
            "v2",
            "ap_exposure_v2_score",
            base_v2,
            "exclude_post2020_china_ap_exposure_v2",
            "m7_exclude_post2020_china_ap_exposure_v2_glm",
            base_v2.loc[~(base_v2["country_iso3"].eq("CHN") & base_v2["year"].ge(2021))].copy(),
            "exclude_country_year_subset=CHN_year_ge_2021",
            v2_notes + ";explicit_exclusion=proxy_for_mr_dominant_post2020_china_rows",
        ),
        (
            "v1",
            "ap_exposure_v1_score",
            base_v1,
            "exclude_post2020_china_ap_exposure_v1",
            "m7_exclude_post2020_china_ap_exposure_v1_glm",
            base_v1.loc[~(base_v1["country_iso3"].eq("CHN") & base_v1["year"].ge(2021))].copy(),
            "exclude_country_year_subset=CHN_year_ge_2021",
            v1_notes + ";explicit_exclusion=proxy_for_mr_dominant_post2020_china_rows",
        ),
        (
            "dtp3",
            "dtp3_coverage",
            base_dtp3,
            "exclude_post2020_china_legacy_dtp3",
            "m7_exclude_post2020_china_legacy_dtp3_glm",
            base_dtp3.loc[~(base_dtp3["country_iso3"].eq("CHN") & base_dtp3["year"].ge(2021))].copy(),
            "exclude_country_year_subset=CHN_year_ge_2021",
            dtp3_notes + ";explicit_exclusion=proxy_for_mr_dominant_post2020_china_rows",
        ),
    ]
    for exposure_family, exposure_column, _, sensitivity_label, model_id, subset, country_filter, notes in special_sensitivity_specs:
        run_single_model(
            results,
            diagnostics,
            subset,
            exposure_column=exposure_column,
            model_id=model_id,
            sensitivity_label=sensitivity_label,
            focal_exposure_family=exposure_family,
            country_filter=country_filter,
            notes=notes,
        )

    leave_one_out_tasks: list[dict[str, object]] = []
    for exposure_family, exposure_column, base_frame, label, base_notes in [
        ("v3", "ap_exposure_v3_score", base_v3, "leave_one_country_out_ap_exposure_v3", v3_notes),
        ("v2", "ap_exposure_v2_score", base_v2, "leave_one_country_out_ap_exposure_v2", v2_notes),
        ("v1", "ap_exposure_v1_score", base_v1, "leave_one_country_out_ap_exposure_v1", v1_notes),
        ("dtp3", "dtp3_coverage", base_dtp3, "leave_one_country_out_legacy_dtp3", dtp3_notes),
    ]:
        countries = sorted(country for country in base_frame["country_iso3"].dropna().unique() if country)
        for excluded_country in countries:
            subset = base_frame.loc[base_frame["country_iso3"].ne(excluded_country)].copy()
            leave_one_out_tasks.append(
                {
                    "frame": subset,
                    "exposure_column": exposure_column,
                    "model_id": f"m7_{label}_exclude_{excluded_country}_glm",
                    "sensitivity_label": label,
                    "focal_exposure_family": exposure_family,
                    "excluded_country_iso3": excluded_country,
                    "country_filter": f"exclude_country={excluded_country}",
                    "notes": base_notes + f";leave_one_country_out={excluded_country}",
                }
            )

    loo_workers = resolve_parallel_workers(len(leave_one_out_tasks), requested_max_workers=max_workers)
    if loo_workers <= 1:
        task_outputs = [run_single_model_task(task) for task in leave_one_out_tasks]
    else:
        with ProcessPoolExecutor(max_workers=loo_workers) as pool:
            task_outputs = list(pool.map(run_single_model_task, leave_one_out_tasks))
    for task_results, task_diagnostics in task_outputs:
        results.extend(task_results)
        diagnostics.extend(task_diagnostics)

    results_frame = pd.DataFrame(results, columns=RESULT_COLUMNS)
    diagnostics_frame = pd.DataFrame(diagnostics, columns=DIAGNOSTIC_COLUMNS)

    results_path = Path(results_out)
    diagnostics_path = Path(diagnostics_out)
    results_path.parent.mkdir(parents=True, exist_ok=True)
    diagnostics_path.parent.mkdir(parents=True, exist_ok=True)
    results_frame.to_csv(results_path, sep="\t", index=False)

    diagnostics_table_path = diagnostics_path.with_suffix(".tsv")
    diagnostics_frame.to_csv(diagnostics_table_path, sep="\t", index=False)

    dataset_path = results_path.with_name("panel_model_country_year_dataset.tsv")
    base_common.to_csv(dataset_path, sep="\t", index=False)
    leave_one_out_path = write_leave_one_out_summary(results_frame, results_path.parent)
    coverage_path = write_coverage_report(base_common, results_path.parent)

    plt.rcParams.update({"figure.figsize": (11, 4.5)})
    figure, axes = plt.subplots(1, 2)

    primary_focal = results_frame.loc[
        (
            ((results_frame["focal_exposure_family"] == "v3") & (results_frame["sensitivity_label"] == "primary_ap_exposure_v3") & (results_frame["estimate_term"] == "ap_exposure_v3_score"))
            | ((results_frame["focal_exposure_family"] == "v2") & (results_frame["sensitivity_label"] == "primary_ap_exposure_v2") & (results_frame["estimate_term"] == "ap_exposure_v2_score"))
            | ((results_frame["focal_exposure_family"] == "v1") & (results_frame["sensitivity_label"] == "primary_ap_exposure_v1") & (results_frame["estimate_term"] == "ap_exposure_v1_score"))
            | ((results_frame["focal_exposure_family"] == "dtp3") & (results_frame["sensitivity_label"] == "legacy_dtp3_proxy") & (results_frame["estimate_term"] == "dtp3_coverage"))
        )
    ].copy()
    if not primary_focal.empty:
        label_map = {"v3": "aPExposure V3", "v2": "aPExposure V2", "v1": "aPExposure V1", "dtp3": "DTP3 only"}
        primary_focal["label"] = primary_focal["focal_exposure_family"].map(label_map)
        primary_focal = primary_focal.sort_values("effect_estimate")
        positions = np.arange(len(primary_focal))
        axes[0].errorbar(
            primary_focal["effect_estimate"],
            positions,
            xerr=[
                primary_focal["effect_estimate"] - primary_focal["ci_lower"],
                primary_focal["ci_upper"] - primary_focal["effect_estimate"],
            ],
            fmt="o",
            color="#0a5c91",
            ecolor="#8aa1b4",
            capsize=3,
        )
        axes[0].axvline(0, color="0.7", linewidth=1)
        axes[0].set_yticks(positions, primary_focal["label"])
        axes[0].set_xlabel("Log-odds estimate")
        axes[0].set_title("Primary exposure comparison")
    else:
        axes[0].text(0.5, 0.5, "No fitted primary coefficients", ha="center", va="center")
        axes[0].set_axis_off()

    loo_v2 = pd.read_csv(leave_one_out_path, sep="\t")
    loo_v2 = loo_v2.loc[loo_v2["focal_exposure_family"].eq("v2")].copy()
    if not loo_v2.empty:
        loo_v2 = loo_v2.sort_values("effect_estimate")
        positions = np.arange(len(loo_v2))
        axes[1].errorbar(
            loo_v2["effect_estimate"],
            positions,
            xerr=[
                loo_v2["effect_estimate"] - loo_v2["ci_lower"],
                loo_v2["ci_upper"] - loo_v2["effect_estimate"],
            ],
            fmt="o",
            color="#b24745",
            ecolor="#d9a6a4",
            capsize=3,
        )
        axes[1].axvline(0, color="0.7", linewidth=1)
        axes[1].set_yticks(positions, loo_v2["excluded_country_iso3"])
        axes[1].set_xlabel("Log-odds estimate")
        axes[1].set_title("Leave-one-country-out V2")
    else:
        axes[1].text(0.5, 0.5, "No leave-one-country-out V2 fits", ha="center", va="center")
        axes[1].set_axis_off()

    figure.tight_layout()
    figure.savefig(diagnostics_path)
    plt.close(figure)

    print(f"Wrote model results: {results_path}")
    print(f"Wrote primary dataset: {dataset_path}")
    print(f"Wrote leave-one-country-out summary: {leave_one_out_path}")
    print(f"Wrote coverage report: {coverage_path}")
    return results_frame, diagnostics_frame


if "snakemake" in globals():
    dataset = prepare_panel_dataset(snakemake.input.exposure, snakemake.input.prevalence)
    run_models(dataset, snakemake.output.results, snakemake.output.diagnostics)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Fit workflow panel association models")
    parser.add_argument("--exposure", required=True, help="Exposure index TSV")
    parser.add_argument("--prevalence", required=True, help="IPW prevalence TSV")
    parser.add_argument("--results-out", required=True, help="Panel model results TSV")
    parser.add_argument("--diagnostics-out", required=True, help="Diagnostics PDF path")
    parser.add_argument(
        "--max-workers",
        type=int,
        default=None,
        help=(
            "Maximum process workers for leave-one-country-out fits. "
            f"Defaults to the {PANEL_MAX_WORKERS_ENV} env var or an auto cap of {PANEL_DEFAULT_MAX_WORKERS}."
        ),
    )
    arguments = parser.parse_args()

    requested_max_workers = arguments.max_workers
    if requested_max_workers is None:
        env_text = str(os.environ.get(PANEL_MAX_WORKERS_ENV, "")).strip()
        if env_text:
            try:
                requested_max_workers = int(env_text)
            except ValueError as exc:
                raise SystemExit(f"ERROR: {PANEL_MAX_WORKERS_ENV} must be an integer, got {env_text!r}") from exc

    dataset = prepare_panel_dataset(arguments.exposure, arguments.prevalence)
    run_models(dataset, arguments.results_out, arguments.diagnostics_out, max_workers=requested_max_workers)
