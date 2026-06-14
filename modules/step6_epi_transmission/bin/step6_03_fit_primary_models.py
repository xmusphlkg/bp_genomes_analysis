#!/usr/bin/env python3
"""Fit primary ecological models from the country-year analysis input table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import numpy as np
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.lib.project_paths import project_module_data_root


MODEL_OUTPUT_COLUMNS = [
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
    "q_value_scope",
    "sensitivity_label",
    "notes",
]

DIAGNOSTIC_COLUMNS = [
    "model_id",
    "n_obs",
    "converged",
    "n_iter",
    "design_rank",
    "log_likelihood",
    "notes",
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def normalize_text(value: object) -> str:
    return str(value or "").strip()


def load_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_float(value: str) -> float | None:
    value = normalize_text(value).replace(",", "")
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def bh_adjust(p_values: list[float]) -> list[float]:
    n = len(p_values)
    order = np.argsort(p_values)
    adjusted = [0.0] * n
    running = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        true_rank = n - rank + 1
        value = p_values[idx] * n / true_rank
        running = min(running, value)
        adjusted[idx] = min(max(running, 0.0), 1.0)
    return adjusted


def z_scores(values: list[float]) -> np.ndarray:
    arr = np.array(values, dtype=float)
    std = arr.std(ddof=0)
    if std == 0:
        return arr * 0.0
    return (arr - arr.mean()) / std


def summarize_input_coverage(rows: list[dict[str, str]]) -> dict[str, int]:
    counts = {
        "core_complete": 0,
        "macrolide_complete": 0,
        "total_antibiotic_complete": 0,
        "non_ap_program_cells": 0,
    }
    for row in rows:
        trials = parse_float(row.get("n_genomes_prn_interpretable", ""))
        successes = parse_float(row.get("n_prn_disrupted", ""))
        reported_cases = parse_float(row.get("reported_cases", ""))
        dtp3 = parse_float(row.get("dtp3_coverage", ""))
        genomes_per_case = parse_float(row.get("genomes_per_case", ""))
        vaccine_program_type = normalize_text(row.get("vaccine_program_type", ""))
        core_complete = None not in {trials, successes, reported_cases, dtp3, genomes_per_case}
        if core_complete and trials and trials > 0 and reported_cases and reported_cases > 0 and successes <= trials:
            counts["core_complete"] += 1
            if vaccine_program_type != "ap_introduced_routine_or_mixed":
                counts["non_ap_program_cells"] += 1
            if parse_float(row.get("macrolide_use_ddd_per_1000_per_day", "")) is not None:
                counts["macrolide_complete"] += 1
            if parse_float(row.get("total_antibiotic_use_ddd_per_1000_per_day", "")) is not None:
                counts["total_antibiotic_complete"] += 1
    return counts


def build_model_dataset(rows: list[dict[str, str]]) -> tuple[np.ndarray, np.ndarray, list[str], dict[str, str]]:
    prepared: list[dict[str, float | str]] = []
    for row in rows:
        trials = parse_float(row.get("n_genomes_prn_interpretable", ""))
        successes = parse_float(row.get("n_prn_disrupted", ""))
        reported_cases = parse_float(row.get("reported_cases", ""))
        dtp3 = parse_float(row.get("dtp3_coverage", ""))
        genomes_per_case = parse_float(row.get("genomes_per_case", ""))
        post_covid_period = parse_float(row.get("post_covid_period", "")) or 0.0
        if None in {trials, successes, reported_cases, dtp3, genomes_per_case}:
            continue
        if trials <= 0 or reported_cases <= 0 or successes > trials:
            continue
        prepared.append(
            {
                "country_iso3": normalize_text(row.get("country_iso3", "")),
                "year": normalize_text(row.get("year", "")),
                "successes": successes,
                "trials": trials,
                "dtp3": dtp3,
                "log_cases": float(np.log1p(reported_cases)),
                "post_covid_period": post_covid_period,
                "genomes_per_case": genomes_per_case,
            }
        )

    if len(prepared) < 5:
        raise ValueError("not enough country-year cells with complete core covariates for primary model fitting")

    dtp3_z = z_scores([float(row["dtp3"]) for row in prepared])
    log_cases_z = z_scores([float(row["log_cases"]) for row in prepared])
    genomes_per_case_z = z_scores([float(row["genomes_per_case"]) for row in prepared])
    X = np.column_stack(
        [
            np.ones(len(prepared)),
            dtp3_z,
            log_cases_z,
            np.array([float(row["post_covid_period"]) for row in prepared], dtype=float),
            genomes_per_case_z,
        ]
    )
    y = np.column_stack(
        [
            np.array([float(row["successes"]) for row in prepared], dtype=float),
            np.array([float(row["trials"]) - float(row["successes"]) for row in prepared], dtype=float),
        ]
    )
    term_names = [
        "Intercept",
        "dtp3_coverage_z",
        "log1p_reported_cases_z",
        "post_covid_period",
        "genomes_per_case_z",
    ]
    metadata = {
        "n_country_year_cells": str(len(prepared)),
        "n_countries": str(len({normalize_text(row["country_iso3"]) for row in prepared})),
        "year_window": f"{min(int(normalize_text(row['year'])) for row in prepared)}-{max(int(normalize_text(row['year'])) for row in prepared)}",
        "country_groups": [normalize_text(row["country_iso3"]) for row in prepared],
    }
    return X, y, term_names, metadata


def fit_primary_glm(
    X: np.ndarray,
    y: np.ndarray,
    country_groups: list[str],
) -> tuple[sm.GLM, str, list[str]]:
    model = sm.GLM(y, X, family=sm.families.Binomial())
    warning_names: list[str] = []
    country_array = np.asarray(country_groups, dtype=object)
    n_clusters = len({group for group in country_groups if group})
    if n_clusters >= 3:
        try:
            result = model.fit(
                maxiter=200,
                disp=0,
                cov_type="cluster",
                cov_kwds={"groups": country_array},
            )
            return result, "country_cluster", warning_names
        except Exception as exc:
            warning_names.append(f"cluster_covariance_fallback={type(exc).__name__}")
    result = model.fit(maxiter=200, disp=0, cov_type="HC1")
    return result, "hc1", warning_names


def build_model_outputs(rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    coverage = summarize_input_coverage(rows)
    X, y, term_names, metadata = build_model_dataset(rows)
    country_groups = [str(value) for value in metadata.pop("country_groups")]
    result, covariance_type, fit_notes = fit_primary_glm(X, y, country_groups)
    conf_int = result.conf_int()
    p_values = [float(value) for value in result.pvalues]
    q_values = bh_adjust(p_values)
    n_clusters = len({group for group in country_groups if group})

    model_id = "int03_primary_binomial_glm_statsmodels_v2"
    notes = (
        f"statsmodels_glm_binomial_{covariance_type}_covariance;"
        "country_random_intercepts_not_fit_due_sparse_country_year_panel;"
        f"country_cluster_robust_covariance={covariance_type};"
        f"country_cluster_count={n_clusters};"
        "q_value_scope=within_primary_model_reported_terms_bh_not_analysis_wide_fdr;"
        "incidence_excluded_because_current_cases_table_lacks_country_level_incidence;"
        f"vaccine_program_type_excluded_due_sparse_nonroutine_or_unknown_contrast_n={coverage['non_ap_program_cells']};"
        f"macrolide_complete_rows_available_but_excluded_from_primary_due_sparse_overlap_n={coverage['macrolide_complete']};"
        f"total_antibiotic_complete_rows_available_but_excluded_from_primary_due_sparse_overlap_n={coverage['total_antibiotic_complete']}"
    )
    if fit_notes:
        notes = f"{notes};{';'.join(fit_notes)}"
    covariates = ",".join(term_names[1:])

    model_rows: list[dict[str, str]] = []
    for index, term_name in enumerate(term_names):
        model_rows.append(
            {
                "model_id": model_id,
                "analysis_cohort": "C",
                "response_variable": "n_prn_disrupted / n_genomes_prn_interpretable",
                "model_family": f"statsmodels_glm_binomial_{covariance_type}_covariance",
                "country_filter": "all_available_country_year_cells_with_complete_core_covariates",
                "year_window": metadata["year_window"],
                "n_country_year_cells": metadata["n_country_year_cells"],
                "n_countries": metadata["n_countries"],
                "covariates": covariates,
                "random_effects": "not_fit",
                "weighting_scheme": "grouped_binomial_trials",
                "estimate_term": term_name,
                "effect_scale": "log_odds",
                "effect_estimate": f"{float(result.params[index]):.6f}",
                "ci_lower": f"{float(conf_int[index, 0]):.6f}",
                "ci_upper": f"{float(conf_int[index, 1]):.6f}",
                "p_value": f"{p_values[index]:.6g}",
                "q_value": f"{q_values[index]:.6g}",
                "q_value_scope": "within_primary_model_reported_terms_bh_not_analysis_wide_fdr",
                "sensitivity_label": "primary_core",
                "notes": notes,
            }
        )

    iteration = result.fit_history.get("iteration", "")
    diagnostic_rows = [
        {
            "model_id": model_id,
            "n_obs": metadata["n_country_year_cells"],
            "converged": "true" if bool(result.converged) else "false",
            "n_iter": str(iteration),
            "design_rank": str(int(np.linalg.matrix_rank(X))),
            "log_likelihood": f"{float(result.llf):.6f}",
            "notes": notes,
        }
    ]
    return model_rows, diagnostic_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Fit primary ecological models.")
    parser.add_argument(
        "--input",
        type=Path,
        default=project_module_data_root("step6_epi_transmission") / "outputs" / "bp_country_year_analysis_input.tsv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=project_module_data_root("step6_epi_transmission") / "outputs" / "bp_country_year_association_models.tsv",
    )
    parser.add_argument(
        "--diagnostics-out",
        type=Path,
        default=project_module_data_root("step6_epi_transmission") / "outputs" / "bp_country_year_model_diagnostics.tsv",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    rows = load_tsv_rows(args.input)
    model_rows, diagnostic_rows = build_model_outputs(rows)
    write_tsv(args.out, MODEL_OUTPUT_COLUMNS, model_rows)
    write_tsv(args.diagnostics_out, DIAGNOSTIC_COLUMNS, diagnostic_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
