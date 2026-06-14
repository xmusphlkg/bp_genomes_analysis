#!/usr/bin/env python3
"""Run labeled sensitivity ecological models without overwriting primary outputs."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import warnings
from collections import defaultdict
from pathlib import Path
from types import ModuleType
from typing import Callable
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

SENSITIVITY_Q_VALUE_SCOPE = "within_sensitivity_model_reported_terms_bh_not_analysis_wide_fdr"

DIAGNOSTIC_COLUMNS = [
    "model_id",
    "n_obs",
    "converged",
    "n_iter",
    "design_rank",
    "log_likelihood",
    "notes",
]

BASE_COVARIATES = [
    {"source_key": "dtp3", "term": "dtp3_coverage_z", "kind": "z"},
    {"source_key": "log_cases", "term": "log1p_reported_cases_z", "kind": "z"},
    {"source_key": "post_covid_period", "term": "post_covid_period", "kind": "raw"},
    {"source_key": "genomes_per_case", "term": "genomes_per_case_z", "kind": "z"},
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


def z_scores(values: list[float]) -> np.ndarray:
    arr = np.array(values, dtype=float)
    std = arr.std(ddof=0)
    if std == 0:
        return arr * 0.0
    return (arr - arr.mean()) / std


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


def load_module(module_name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"unable to load module at {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_phylogeny_analysis_input(
    manifest_path: Path,
    manifest_label: str,
    schema_path: Path,
    mechanism_calls_path: Path,
    validation_path: Path,
    marker_table_path: Path,
    country_map_path: Path,
    public_health_path: Path,
    out_path: Path,
) -> list[dict[str, str]]:
    module_bin = Path(__file__).resolve().parent
    step6_01 = load_module(
        f"step6_01_{manifest_label}",
        module_bin / "step6_01_build_country_year_genomic_summaries.py",
    )
    step6_02 = load_module(
        f"step6_02_{manifest_label}",
        module_bin / "step6_02_join_public_health.py",
    )

    manifest_rows = step6_01.load_tsv_rows(manifest_path)
    selected_rows = [
        row
        for row in manifest_rows
        if normalize_text(row.get("analysis_cohort_id", "")) == "A"
        and normalize_text(row.get("phylogeny_selected_for_tree", "")).casefold() == "true"
    ]
    genomic_rows = step6_01.build_output_rows(
        cohort_rows=selected_rows,
        mechanism_rows=step6_01.load_tsv_rows(mechanism_calls_path),
        validation_rows=step6_01.load_tsv_rows(validation_path),
        marker_rows=step6_01.load_tsv_rows(marker_table_path),
        country_map=step6_01.load_country_map(country_map_path),
    )
    for row in genomic_rows:
        notes = normalize_text(row.get("notes", ""))
        row["notes"] = notes.replace(
            "country_year_genomic_summary_from_cohort_C",
            "country_year_genomic_summary_from_cohort_A",
        )

    fieldnames = step6_02.load_schema(schema_path)
    analysis_rows = step6_02.build_output_rows(
        genomic_rows=genomic_rows,
        public_health_rows=step6_02.load_tsv_rows(public_health_path),
        fieldnames=fieldnames,
    )
    for row in analysis_rows:
        notes = [normalize_text(row.get("notes", "")), f"analysis_input_rebuilt_from={manifest_path.name}"]
        row["notes"] = ";".join(note for note in notes if note)

    write_tsv(out_path, fieldnames, analysis_rows)
    return analysis_rows


def is_core_complete(row: dict[str, str]) -> bool:
    trials = parse_float(row.get("n_genomes_prn_interpretable", ""))
    successes = parse_float(row.get("n_prn_disrupted", ""))
    reported_cases = parse_float(row.get("reported_cases", ""))
    dtp3 = parse_float(row.get("dtp3_coverage", ""))
    genomes_per_case = parse_float(row.get("genomes_per_case", ""))
    if None in {trials, successes, reported_cases, dtp3, genomes_per_case}:
        return False
    return bool(trials and trials > 0 and reported_cases and reported_cases > 0 and successes <= trials)


def load_first_routine_ap_years(path: Path) -> dict[str, int]:
    first_years: dict[str, int] = {}
    for row in load_tsv_rows(path):
        if normalize_text(row.get("vaccine_program_type", "")) != "ap_introduced_routine_or_mixed":
            continue
        year = parse_int(row.get("year_start", "")) or parse_int(row.get("program_change_year", ""))
        country_iso3 = normalize_text(row.get("country_iso3", ""))
        if not country_iso3 or year is None:
            continue
        first_years[country_iso3] = min(year, first_years.get(country_iso3, year))
    return first_years


def country_totals(rows: list[dict[str, str]]) -> dict[str, float]:
    totals: dict[str, float] = defaultdict(float)
    for row in rows:
        country_iso3 = normalize_text(row.get("country_iso3", ""))
        genomes_total = parse_float(row.get("n_genomes_total", "")) or 0.0
        if country_iso3:
            totals[country_iso3] += genomes_total
    return totals


def build_prepared_rows(
    rows: list[dict[str, str]],
    row_filter: Callable[[dict[str, str]], bool],
    extra_value_builders: dict[str, Callable[[dict[str, str]], float | None]] | None = None,
) -> list[dict[str, object]]:
    prepared: list[dict[str, object]] = []
    for row in rows:
        if not row_filter(row):
            continue
        trials = parse_float(row.get("n_genomes_prn_interpretable", ""))
        successes = parse_float(row.get("n_prn_disrupted", ""))
        reported_cases = parse_float(row.get("reported_cases", ""))
        dtp3 = parse_float(row.get("dtp3_coverage", ""))
        genomes_per_case = parse_float(row.get("genomes_per_case", ""))
        year = parse_int(row.get("year", ""))
        if None in {trials, successes, reported_cases, dtp3, genomes_per_case, year}:
            continue
        if trials <= 0 or reported_cases <= 0 or successes > trials:
            continue
        record: dict[str, object] = {
            "country_iso3": normalize_text(row.get("country_iso3", "")),
            "year": year,
            "successes": successes,
            "trials": trials,
            "dtp3": dtp3,
            "log_cases": float(np.log1p(reported_cases)),
            "post_covid_period": parse_float(row.get("post_covid_period", "")) or 0.0,
            "genomes_per_case": genomes_per_case,
        }
        for key, builder in (extra_value_builders or {}).items():
            record[key] = builder(row)
        prepared.append(record)
    return prepared


def format_year_window(prepared_rows: list[dict[str, object]]) -> str:
    years = [int(row["year"]) for row in prepared_rows]
    return f"{min(years)}-{max(years)}"


def fit_grouped_binomial_with_country_robust_covariance(
    X: np.ndarray,
    y: np.ndarray,
    country_groups: list[str],
) -> tuple[object, str, list[str]]:
    model = sm.GLM(y, X, family=sm.families.Binomial())
    fit_notes: list[str] = []
    n_clusters = len({group for group in country_groups if group})
    if n_clusters >= 3:
        try:
            result = model.fit(
                maxiter=200,
                disp=0,
                cov_type="cluster",
                cov_kwds={"groups": np.asarray(country_groups, dtype=object)},
            )
            return result, "country_cluster", fit_notes
        except Exception as exc:
            fit_notes.append(f"cluster_covariance_fallback={type(exc).__name__}")
    result = model.fit(maxiter=200, disp=0, cov_type="HC1")
    return result, "hc1", fit_notes


def fit_labeled_glm(
    rows: list[dict[str, str]],
    model_id: str,
    analysis_cohort: str,
    sensitivity_label: str,
    country_filter: str,
    row_filter: Callable[[dict[str, str]], bool],
    base_notes: list[str],
    extra_covariates: list[dict[str, str]] | None = None,
    extra_value_builders: dict[str, Callable[[dict[str, str]], float | None]] | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    covariate_specs = BASE_COVARIATES + (extra_covariates or [])
    prepared_rows = build_prepared_rows(rows, row_filter=row_filter, extra_value_builders=extra_value_builders)
    prepared_rows = [
        row
        for row in prepared_rows
        if all(row.get(spec["source_key"]) is not None for spec in covariate_specs)
    ]
    if len(prepared_rows) < 5:
        raise ValueError(f"{model_id} has fewer than 5 complete country-year cells")

    X_columns: list[np.ndarray] = [np.ones(len(prepared_rows), dtype=float)]
    term_names = ["Intercept"]
    dropped_covariates: list[str] = []
    for spec in covariate_specs:
        values = np.array([float(row[spec["source_key"]]) for row in prepared_rows], dtype=float)
        if values.std(ddof=0) == 0:
            dropped_covariates.append(spec["term"])
            continue
        if spec["kind"] == "z":
            values = z_scores(values.tolist())
        X_columns.append(values)
        term_names.append(spec["term"])

    if len(term_names) < 2:
        raise ValueError(f"{model_id} has no variable covariates after filtering")

    X = np.column_stack(X_columns)
    y = np.column_stack(
        [
            np.array([float(row["successes"]) for row in prepared_rows], dtype=float),
            np.array([float(row["trials"]) - float(row["successes"]) for row in prepared_rows], dtype=float),
        ]
    )
    country_groups = [normalize_text(row["country_iso3"]) for row in prepared_rows]
    n_clusters = len({group for group in country_groups if group})

    warning_names: list[str] = []
    covariance_notes: list[str] = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result, covariance_type, covariance_notes = fit_grouped_binomial_with_country_robust_covariance(
            X,
            y,
            country_groups,
        )
        warning_names = sorted({type(item.message).__name__ for item in caught})

    conf_int = result.conf_int()
    p_values = [float(value) for value in result.pvalues]
    q_values = bh_adjust(p_values)
    covariates = ",".join(term_names[1:])

    notes = list(base_notes)
    if dropped_covariates:
        notes.append(f"dropped_zero_variance_covariates={','.join(dropped_covariates)}")
    if warning_names:
        notes.append(f"fit_warning_types={','.join(warning_names)}")
    if covariance_notes:
        notes.extend(covariance_notes)
    notes.append(f"country_cluster_robust_covariance={covariance_type}")
    notes.append(f"country_cluster_count={n_clusters}")
    notes.append(f"sensitivity_label={sensitivity_label}")
    notes.append(f"q_value_scope={SENSITIVITY_Q_VALUE_SCOPE}")
    note_text = ";".join(note for note in notes if note)

    model_rows: list[dict[str, str]] = []
    for index, term_name in enumerate(term_names):
        model_rows.append(
            {
                "model_id": model_id,
                "analysis_cohort": analysis_cohort,
                "response_variable": "n_prn_disrupted / n_genomes_prn_interpretable",
                "model_family": f"statsmodels_glm_binomial_{covariance_type}_covariance",
                "country_filter": country_filter,
                "year_window": format_year_window(prepared_rows),
                "n_country_year_cells": str(len(prepared_rows)),
                "n_countries": str(len({normalize_text(row["country_iso3"]) for row in prepared_rows})),
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
                "q_value_scope": SENSITIVITY_Q_VALUE_SCOPE,
                "sensitivity_label": sensitivity_label,
                "notes": note_text,
            }
        )

    diagnostic_rows = [
        {
            "model_id": model_id,
            "n_obs": str(len(prepared_rows)),
            "converged": "true" if bool(result.converged) else "false",
            "n_iter": str(result.fit_history.get("iteration", "")),
            "design_rank": str(int(np.linalg.matrix_rank(X))),
            "log_likelihood": f"{float(result.llf):.6f}",
            "notes": note_text,
        }
    ]
    return model_rows, diagnostic_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run labeled sensitivity ecological models.")
    step6_root = project_module_data_root("step6_epi_transmission")
    step5_root = project_module_data_root("step5_phylogeny_asr")
    step4_root = project_module_data_root("step4_prn_validation")
    step2_root = project_module_data_root("step2_typing")
    public_health_root = project_module_data_root("public_health")
    parser.add_argument(
        "--analysis-input",
        type=Path,
        default=step6_root / "outputs" / "bp_country_year_analysis_input.tsv",
    )
    parser.add_argument(
        "--primary-models",
        type=Path,
        default=step6_root / "outputs" / "bp_country_year_association_models.tsv",
    )
    parser.add_argument(
        "--primary-diagnostics",
        type=Path,
        default=step6_root / "outputs" / "bp_country_year_model_diagnostics.tsv",
    )
    parser.add_argument(
        "--schema",
        type=Path,
        default=step6_root / "outputs" / "bp_country_year_analysis_input.schema.tsv",
    )
    parser.add_argument(
        "--balanced-manifest",
        type=Path,
        default=step5_root / "outputs" / "bp_phylogeny_manifest_balanced.tsv",
    )
    parser.add_argument(
        "--full-manifest",
        type=Path,
        default=step5_root / "outputs" / "bp_phylogeny_manifest_full.tsv",
    )
    parser.add_argument(
        "--mechanism-calls",
        type=Path,
        default=step4_root / "outputs" / "bp_prn_mechanism_calls.tsv",
    )
    parser.add_argument(
        "--validation",
        type=Path,
        default=step4_root / "outputs" / "bp_prn_read_validation.tsv",
    )
    parser.add_argument(
        "--marker-table",
        type=Path,
        default=step2_root / "outputs" / "bp_qc_merged_mlst_markers.tsv",
    )
    parser.add_argument(
        "--country-map",
        type=Path,
        default=public_health_root / "outputs" / "ph_country_name_map.tsv",
    )
    parser.add_argument(
        "--public-health",
        type=Path,
        default=public_health_root / "outputs" / "ph_country_year_master.tsv",
    )
    parser.add_argument(
        "--program-metadata",
        type=Path,
        default=public_health_root / "outputs" / "ph_country_program_metadata.tsv",
    )
    parser.add_argument(
        "--balanced-analysis-out",
        type=Path,
        default=step6_root / "outputs" / "bp_country_year_analysis_input_phylo_balanced.tsv",
    )
    parser.add_argument(
        "--full-analysis-out",
        type=Path,
        default=step6_root / "outputs" / "bp_country_year_analysis_input_phylo_full.tsv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=step6_root / "outputs" / "bp_country_year_association_models_with_sensitivity.tsv",
    )
    parser.add_argument(
        "--diagnostics-out",
        type=Path,
        default=step6_root / "outputs" / "bp_country_year_model_diagnostics_with_sensitivity.tsv",
    )
    parser.add_argument(
        "--rebuild-phylogeny-analysis-inputs",
        action="store_true",
        help=(
            "Explicitly rebuild the balanced/full phylogeny analysis-input tables before fitting sensitivities. "
            "Without this flag, existing immutable inputs are loaded fail-closed."
        ),
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()

    primary_model_rows = load_tsv_rows(args.primary_models)
    primary_diagnostic_rows = load_tsv_rows(args.primary_diagnostics)
    core_rows = load_tsv_rows(args.analysis_input)
    if args.rebuild_phylogeny_analysis_inputs:
        balanced_rows = build_phylogeny_analysis_input(
            manifest_path=args.balanced_manifest,
            manifest_label="balanced",
            schema_path=args.schema,
            mechanism_calls_path=args.mechanism_calls,
            validation_path=args.validation,
            marker_table_path=args.marker_table,
            country_map_path=args.country_map,
            public_health_path=args.public_health,
            out_path=args.balanced_analysis_out,
        )
        full_rows = build_phylogeny_analysis_input(
            manifest_path=args.full_manifest,
            manifest_label="full",
            schema_path=args.schema,
            mechanism_calls_path=args.mechanism_calls,
            validation_path=args.validation,
            marker_table_path=args.marker_table,
            country_map_path=args.country_map,
            public_health_path=args.public_health,
            out_path=args.full_analysis_out,
        )
        phylo_input_note = "phylogeny_analysis_input_rebuilt_explicitly"
    else:
        if not args.balanced_analysis_out.exists() or not args.full_analysis_out.exists():
            raise FileNotFoundError(
                "Balanced/full phylogeny analysis inputs are missing. "
                "Re-run with --rebuild-phylogeny-analysis-inputs to rebuild them explicitly."
            )
        balanced_rows = load_tsv_rows(args.balanced_analysis_out)
        full_rows = load_tsv_rows(args.full_analysis_out)
        phylo_input_note = "phylogeny_analysis_input_loaded_from_existing_outputs"

    first_routine_ap_years = load_first_routine_ap_years(args.program_metadata)
    all_country_totals = country_totals(core_rows)
    countries_with_total_genomes_ge10 = {
        country_iso3 for country_iso3, total in all_country_totals.items() if total >= 10
    }

    def with_known_routine_ap(row: dict[str, str]) -> float | None:
        country_iso3 = normalize_text(row.get("country_iso3", ""))
        year = parse_int(row.get("year", ""))
        first_year = first_routine_ap_years.get(country_iso3)
        if year is None or first_year is None:
            return None
        return float(year - first_year)

    sensitivity_configs = [
        {
            "rows": balanced_rows,
            "model_id": "int04_sensitivity_cohort_a_balanced_binomial_glm_v1",
            "analysis_cohort": "A",
            "sensitivity_label": "cohort_A_phylogeny_balanced",
            "country_filter": "phylogeny_selected_country_year_cells_from_balanced_manifest",
            "row_filter": lambda row: is_core_complete(row),
            "base_notes": [
                "statsmodels_glm_binomial_country_cluster_or_hc1_covariance",
                "int04_sensitivity_analysis",
                phylo_input_note,
                f"analysis_input_source={args.balanced_analysis_out.name}",
            ],
        },
        {
            "rows": full_rows,
            "model_id": "int04_sensitivity_cohort_a_full_binomial_glm_v1",
            "analysis_cohort": "A",
            "sensitivity_label": "cohort_A_phylogeny_full",
            "country_filter": "phylogeny_selected_country_year_cells_from_full_manifest",
            "row_filter": lambda row: is_core_complete(row),
            "base_notes": [
                "statsmodels_glm_binomial_country_cluster_or_hc1_covariance",
                "int04_sensitivity_analysis",
                phylo_input_note,
                f"analysis_input_source={args.full_analysis_out.name}",
            ],
        },
        {
            "rows": core_rows,
            "model_id": "int04_sensitivity_pre2020_binomial_glm_v1",
            "analysis_cohort": "C",
            "sensitivity_label": "year_window_pre2020",
            "country_filter": "all_available_country_year_cells_with_complete_core_covariates",
            "row_filter": lambda row: is_core_complete(row) and (parse_int(row.get("year", "")) or 0) <= 2019,
            "base_notes": [
                "statsmodels_glm_binomial_country_cluster_or_hc1_covariance",
                "int04_sensitivity_analysis",
                "year_window_restricted_to_pre2020",
            ],
        },
        {
            "rows": core_rows,
            "model_id": "int04_sensitivity_country_total_genomes_ge10_binomial_glm_v1",
            "analysis_cohort": "C",
            "sensitivity_label": "country_filter_total_genomes_ge10",
            "country_filter": "countries_with_total_n_genomes_ge10_in_primary_analysis_input",
            "row_filter": lambda row: is_core_complete(row)
            and normalize_text(row.get("country_iso3", "")) in countries_with_total_genomes_ge10,
            "base_notes": [
                "statsmodels_glm_binomial_country_cluster_or_hc1_covariance",
                "int04_sensitivity_analysis",
                "country_filter_total_n_genomes_ge10",
            ],
        },
        {
            "rows": core_rows,
            "model_id": "int04_sensitivity_exclude_unclear_vaccine_binomial_glm_v1",
            "analysis_cohort": "C",
            "sensitivity_label": "country_filter_exclude_unclear_vaccine",
            "country_filter": "exclude_rows_with_whole_cell_or_unknown_vaccine_program_type",
            "row_filter": lambda row: is_core_complete(row)
            and normalize_text(row.get("vaccine_program_type", "")) != "whole_cell_or_unknown",
            "base_notes": [
                "statsmodels_glm_binomial_country_cluster_or_hc1_covariance",
                "int04_sensitivity_analysis",
                "exclude_unclear_vaccine_program_rows",
            ],
        },
        {
            "rows": core_rows,
            "model_id": "int04_sensitivity_years_since_routine_ap_binomial_glm_v1",
            "analysis_cohort": "C",
            "sensitivity_label": "vaccine_timing_years_since_routine_ap",
            "country_filter": "all_available_country_year_cells_with_complete_core_covariates_and_known_routine_ap_intro_year",
            "row_filter": lambda row: is_core_complete(row),
            "base_notes": [
                "statsmodels_glm_binomial_country_cluster_or_hc1_covariance",
                "int04_sensitivity_analysis",
                "negative_years_since_routine_ap_intro_indicate_pre_intro_country_years",
            ],
            "extra_covariates": [
                {
                    "source_key": "years_since_routine_ap_intro",
                    "term": "years_since_routine_ap_intro_z",
                    "kind": "z",
                }
            ],
            "extra_value_builders": {
                "years_since_routine_ap_intro": with_known_routine_ap,
            },
        },
    ]

    sensitivity_model_rows: list[dict[str, str]] = []
    sensitivity_diagnostic_rows: list[dict[str, str]] = []
    for config in sensitivity_configs:
        model_rows, diagnostic_rows = fit_labeled_glm(
            rows=config["rows"],
            model_id=config["model_id"],
            analysis_cohort=config["analysis_cohort"],
            sensitivity_label=config["sensitivity_label"],
            country_filter=config["country_filter"],
            row_filter=config["row_filter"],
            base_notes=config["base_notes"],
            extra_covariates=config.get("extra_covariates"),
            extra_value_builders=config.get("extra_value_builders"),
        )
        sensitivity_model_rows.extend(model_rows)
        sensitivity_diagnostic_rows.extend(diagnostic_rows)

    write_tsv(args.out, MODEL_OUTPUT_COLUMNS, primary_model_rows + sensitivity_model_rows)
    write_tsv(args.diagnostics_out, DIAGNOSTIC_COLUMNS, primary_diagnostic_rows + sensitivity_diagnostic_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
