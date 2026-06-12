#!/usr/bin/env python3
"""Run AMU-specific exploratory sensitivity analyses on exact-overlap country-year subsets."""

from __future__ import annotations

import argparse
import csv
import math
import warnings
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
    "sensitivity_label",
    "notes",
]

DIAGNOSTIC_COLUMNS = [
    "model_group_id",
    "source_table",
    "analysis_cohort",
    "amu_metric",
    "country_filter",
    "year_window",
    "n_obs",
    "n_countries",
    "standard_glm_converged",
    "standard_glm_warning_types",
    "standard_glm_log_likelihood",
    "dropped_covariates",
    "ridge_alphas_fit",
    "run_status",
    "notes",
]

OVERLAP_MANIFEST_COLUMNS = [
    "analysis_cohort",
    "source_table",
    "amu_metric",
    "country_iso3",
    "year",
    "n_genomes_prn_interpretable",
    "n_prn_disrupted",
    "reported_cases",
    "dtp3_coverage",
    "amu_value",
    "country_filter",
    "year_window",
    "run_status",
    "notes",
]

RIDGE_ALPHAS = [0.01, 0.1, 1.0]
MIN_EXPLORATORY_OBS = 8


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def normalize_text(value: object) -> str:
    return str(value or "").strip()


def parse_float(value: str) -> float | None:
    text = normalize_text(value).replace(",", "")
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_int(value: str) -> int | None:
    text = normalize_text(value)
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def load_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def z_scores(values: list[float]) -> np.ndarray:
    arr = np.array(values, dtype=float)
    std = arr.std(ddof=0)
    if std == 0:
        return arr * 0.0
    return (arr - arr.mean()) / std


def format_alpha(alpha: float) -> str:
    return str(alpha).replace(".", "p")


def format_year_window(records: list[dict[str, object]]) -> str:
    years = [int(record["year"]) for record in records]
    return f"{min(years)}-{max(years)}"


def build_exact_overlap_records(rows: list[dict[str, str]], amu_key: str) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for row in rows:
        trials = parse_float(row.get("n_genomes_prn_interpretable", ""))
        successes = parse_float(row.get("n_prn_disrupted", ""))
        reported_cases = parse_float(row.get("reported_cases", ""))
        dtp3 = parse_float(row.get("dtp3_coverage", ""))
        amu_value = parse_float(row.get(amu_key, ""))
        year = parse_int(row.get("year", ""))
        if None in {trials, successes, reported_cases, dtp3, amu_value, year}:
            continue
        if trials <= 0 or reported_cases <= 0 or successes > trials:
            continue
        records.append(
            {
                "country_iso3": normalize_text(row.get("country_iso3", "")),
                "year": year,
                "successes": successes,
                "trials": trials,
                "reported_cases": reported_cases,
                "dtp3": dtp3,
                "log_cases": math.log1p(reported_cases),
                "amu": amu_value,
            }
        )
    return records


def design_matrix(records: list[dict[str, object]]) -> tuple[np.ndarray, list[str], list[str]]:
    X_columns: list[np.ndarray] = [np.ones(len(records), dtype=float)]
    term_names = ["Intercept"]
    dropped_covariates: list[str] = []
    covariate_specs = [
        ("dtp3_coverage_z", "dtp3"),
        ("log1p_reported_cases_z", "log_cases"),
        ("amu_z", "amu"),
    ]
    for term_name, key in covariate_specs:
        values = np.array([float(record[key]) for record in records], dtype=float)
        if values.std(ddof=0) == 0:
            dropped_covariates.append(term_name)
            continue
        X_columns.append(z_scores(values.tolist()))
        term_names.append(term_name)
    return np.column_stack(X_columns), term_names, dropped_covariates


def grouped_binomial_response(records: list[dict[str, object]]) -> np.ndarray:
    return np.column_stack(
        [
            np.array([float(record["successes"]) for record in records], dtype=float),
            np.array([float(record["trials"]) - float(record["successes"]) for record in records], dtype=float),
        ]
    )


def fit_standard_glm(X: np.ndarray, y: np.ndarray) -> tuple[sm.GLM, list[str]]:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = sm.GLM(y, X, family=sm.families.Binomial()).fit(maxiter=200, disp=0)
    warning_types = sorted({type(item.message).__name__ for item in caught})
    return result, warning_types


def fit_ridge_glm(X: np.ndarray, y: np.ndarray, alpha: float):
    return sm.GLM(y, X, family=sm.families.Binomial()).fit_regularized(
        alpha=alpha,
        L1_wt=0.0,
        maxiter=1000,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run exploratory AMU-only sensitivity analyses.")
    step6_root = project_module_data_root("step6_epi_transmission")
    parser.add_argument(
        "--core-input",
        type=Path,
        default=step6_root / "outputs" / "bp_country_year_analysis_input.tsv",
    )
    parser.add_argument(
        "--balanced-input",
        type=Path,
        default=step6_root / "outputs" / "bp_country_year_analysis_input_phylo_balanced.tsv",
    )
    parser.add_argument(
        "--full-input",
        type=Path,
        default=step6_root / "outputs" / "bp_country_year_analysis_input_phylo_full.tsv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=step6_root / "outputs" / "bp_country_year_amu_exploratory_models.tsv",
    )
    parser.add_argument(
        "--diagnostics-out",
        type=Path,
        default=step6_root / "outputs" / "bp_country_year_amu_exploratory_diagnostics.tsv",
    )
    parser.add_argument(
        "--overlap-manifest-out",
        type=Path,
        default=step6_root / "outputs" / "bp_country_year_amu_exploratory_overlap_manifest.tsv",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()

    datasets = [
        {
            "source_table": args.core_input.name,
            "analysis_cohort": "C",
            "rows": load_tsv_rows(args.core_input),
        },
        {
            "source_table": args.balanced_input.name,
            "analysis_cohort": "A",
            "rows": load_tsv_rows(args.balanced_input),
        },
        {
            "source_table": args.full_input.name,
            "analysis_cohort": "A",
            "rows": load_tsv_rows(args.full_input),
        },
    ]
    amu_metrics = [
        {
            "key": "macrolide_use_ddd_per_1000_per_day",
            "label": "macrolide_use",
        },
        {
            "key": "total_antibiotic_use_ddd_per_1000_per_day",
            "label": "total_antibiotic_use",
        },
    ]

    model_rows: list[dict[str, str]] = []
    diagnostic_rows: list[dict[str, str]] = []
    overlap_manifest_rows: list[dict[str, str]] = []

    for dataset in datasets:
        for metric in amu_metrics:
            records = build_exact_overlap_records(dataset["rows"], metric["key"])
            n_obs = len(records)
            n_countries = len({normalize_text(record["country_iso3"]) for record in records})
            year_window = format_year_window(records) if records else ""
            country_filter = f"exact_overlap_rows_with_{metric['key']}_and_complete_core_covariates"
            model_group_id = (
                f"int04_amu_exploratory_{Path(dataset['source_table']).stem}_{metric['label']}_ridge_path_v1"
            )
            notes = [
                "exploratory_only_not_for_primary_inference",
                "standard_glm_attempted_before_ridge_path",
                f"min_obs_required_for_ridge={MIN_EXPLORATORY_OBS}",
            ]

            standard_glm_converged = ""
            standard_glm_warning_types = ""
            standard_glm_log_likelihood = ""
            dropped_covariates = ""
            ridge_alphas_fit = ""
            run_status = "not_fit_due_no_overlap_rows"

            if records:
                X, term_names, dropped = design_matrix(records)
                y = grouped_binomial_response(records)
                dropped_covariates = ",".join(dropped)
                if n_obs < MIN_EXPLORATORY_OBS:
                    run_status = f"not_fit_due_n_obs_lt_{MIN_EXPLORATORY_OBS}"
                    notes.append("standard_glm_skipped_below_min_exploratory_obs_threshold")
                else:
                    standard_result, warning_types = fit_standard_glm(X, y)
                    standard_glm_converged = "true" if bool(standard_result.converged) else "false"
                    standard_glm_warning_types = ",".join(warning_types)
                    standard_glm_log_likelihood = f"{float(standard_result.llf):.6f}"
                    run_status = "ridge_path_fit"
                    ridge_alphas_fit = ",".join(str(alpha) for alpha in RIDGE_ALPHAS)
                    amu_index = term_names.index("amu_z") if "amu_z" in term_names else None
                    if amu_index is not None:
                        for alpha in RIDGE_ALPHAS:
                            ridge_result = fit_ridge_glm(X, y, alpha)
                            alpha_text = format_alpha(alpha)
                            model_id = f"{model_group_id}_alpha_{alpha_text}"
                            sensitivity_label = f"exploratory_amu_{metric['label']}_ridge_alpha_{alpha_text}"
                            row_notes = notes + [
                                f"source_table={dataset['source_table']}",
                                f"penalty_alpha={alpha}",
                                f"standard_glm_warning_types={standard_glm_warning_types or 'none'}",
                            ]
                            if dropped_covariates:
                                row_notes.append(f"dropped_covariates={dropped_covariates}")
                            model_rows.append(
                                {
                                    "model_id": model_id,
                                    "analysis_cohort": dataset["analysis_cohort"],
                                    "response_variable": "n_prn_disrupted / n_genomes_prn_interpretable",
                                    "model_family": "statsmodels_glm_binomial_ridge_no_random_effects",
                                    "country_filter": country_filter,
                                    "year_window": year_window,
                                    "n_country_year_cells": str(n_obs),
                                    "n_countries": str(n_countries),
                                    "covariates": ",".join(term_names[1:]),
                                    "random_effects": "not_fit",
                                    "weighting_scheme": "grouped_binomial_trials",
                                    "estimate_term": "amu_z",
                                    "effect_scale": "ridge_penalized_log_odds",
                                    "effect_estimate": f"{float(ridge_result.params[amu_index]):.6f}",
                                    "ci_lower": "",
                                    "ci_upper": "",
                                    "p_value": "",
                                    "q_value": "",
                                    "sensitivity_label": sensitivity_label,
                                    "notes": ";".join(row_notes),
                                }
                            )
                    else:
                        run_status = "not_fit_due_amu_zero_variance"
            diagnostic_rows.append(
                {
                    "model_group_id": model_group_id,
                    "source_table": dataset["source_table"],
                    "analysis_cohort": dataset["analysis_cohort"],
                    "amu_metric": metric["label"],
                    "country_filter": country_filter,
                    "year_window": year_window,
                    "n_obs": str(n_obs),
                    "n_countries": str(n_countries),
                    "standard_glm_converged": standard_glm_converged,
                    "standard_glm_warning_types": standard_glm_warning_types,
                    "standard_glm_log_likelihood": standard_glm_log_likelihood,
                    "dropped_covariates": dropped_covariates,
                    "ridge_alphas_fit": ridge_alphas_fit,
                    "run_status": run_status,
                    "notes": ";".join(notes),
                }
            )

            for record in sorted(records, key=lambda item: (normalize_text(item["country_iso3"]), int(item["year"]))):
                overlap_manifest_rows.append(
                    {
                        "analysis_cohort": dataset["analysis_cohort"],
                        "source_table": dataset["source_table"],
                        "amu_metric": metric["label"],
                        "country_iso3": normalize_text(record["country_iso3"]),
                        "year": str(int(record["year"])),
                        "n_genomes_prn_interpretable": f"{float(record['trials']):.0f}",
                        "n_prn_disrupted": f"{float(record['successes']):.0f}",
                        "reported_cases": f"{float(record['reported_cases']):.6f}",
                        "dtp3_coverage": f"{float(record['dtp3']):.6f}",
                        "amu_value": f"{float(record['amu']):.6f}",
                        "country_filter": country_filter,
                        "year_window": year_window,
                        "run_status": run_status,
                        "notes": "exact_overlap_row_with_complete_core_covariates_and_non_null_amu",
                    }
                )

    write_tsv(args.out, MODEL_OUTPUT_COLUMNS, model_rows)
    write_tsv(args.diagnostics_out, DIAGNOSTIC_COLUMNS, diagnostic_rows)
    write_tsv(args.overlap_manifest_out, OVERLAP_MANIFEST_COLUMNS, overlap_manifest_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
