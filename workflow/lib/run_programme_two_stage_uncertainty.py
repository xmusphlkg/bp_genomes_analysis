#!/usr/bin/env python3
"""Propagate missingness-model uncertainty into programme-surveillance summaries."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from patsy import build_design_matrices

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from build_programme_country_period_panel import build_period_panel, prepare_annual_dataset_from_frames
from missingness_model import (
    FULL_FEATURES,
    REDUCED_FEATURES,
    fit_logistic_model,
    load_inputs,
    score_rows,
)
from run_programme_surveillance_models import (
    CLASS_ORDER,
    build_formula,
    fit_response_model,
    model_frame,
    parameter_draws,
)
from workflow.lib.project_paths import project_workflow_root


PROGRAMME_UNCERTAINTY_MAX_WORKERS_ENV = "PROGRAMME_UNCERTAINTY_MAX_WORKERS"
PROGRAMME_UNCERTAINTY_DEFAULT_MAX_WORKERS = 32
_BOOTSTRAP_CONTEXT: dict[str, object] = {}
DEFAULT_WORKFLOW_ROOT = project_workflow_root()


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
        else min(available_cpu, PROGRAMME_UNCERTAINTY_DEFAULT_MAX_WORKERS)
    )
    workers = max(1, min(int(workers), available_cpu))
    if task_count is not None:
        workers = min(workers, max(1, int(task_count)))
    return workers


def resolve_env_max_workers() -> int | None:
    env_text = str(os.environ.get(PROGRAMME_UNCERTAINTY_MAX_WORKERS_ENV, "")).strip()
    if not env_text:
        return None
    try:
        return int(env_text)
    except ValueError as exc:
        raise SystemExit(
            f"ERROR: {PROGRAMME_UNCERTAINTY_MAX_WORKERS_ENV} must be an integer, got {env_text!r}"
        ) from exc


def build_replicate_seeds(n_replicates: int, base_seed: int) -> list[int]:
    if n_replicates <= 0:
        return []
    seed_sequence = np.random.SeedSequence(int(base_seed))
    child_sequences = seed_sequence.spawn(int(n_replicates))
    return [int(sequence.generate_state(1, dtype=np.uint64)[0]) for sequence in child_sequences]


def set_bootstrap_context(
    manifest_frame: pd.DataFrame,
    exposure_frame: pd.DataFrame,
    *,
    origin_descendants_dir: str,
    weight_truncation: float,
    bin_size: int,
    min_interpretable: int,
) -> None:
    global _BOOTSTRAP_CONTEXT
    _BOOTSTRAP_CONTEXT = {
        "manifest_frame": manifest_frame,
        "exposure_frame": exposure_frame,
        "origin_descendants_dir": origin_descendants_dir,
        "weight_truncation": weight_truncation,
        "bin_size": bin_size,
        "min_interpretable": min_interpretable,
    }


def _init_bootstrap_worker(
    manifest_frame: pd.DataFrame,
    exposure_frame: pd.DataFrame,
    origin_descendants_dir: str,
    weight_truncation: float,
    bin_size: int,
    min_interpretable: int,
) -> None:
    set_bootstrap_context(
        manifest_frame,
        exposure_frame,
        origin_descendants_dir=origin_descendants_dir,
        weight_truncation=weight_truncation,
        bin_size=bin_size,
        min_interpretable=min_interpretable,
    )


def parse_bool(value: object) -> bool:
    text = str(value).strip().lower()
    return text in {"true", "1", "yes"}


def coerce_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def bootstrap_prediction_frame(manifest_frame: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    labeled = manifest_frame.dropna(subset=["prn_interpretable_numeric"]).copy()
    sampled = labeled.iloc[rng.integers(0, len(labeled), size=len(labeled))].copy()

    full_rows = sampled.dropna(subset=FULL_FEATURES).copy()
    reduced_rows = sampled.dropna(subset=REDUCED_FEATURES).copy()
    full_model, full_means, full_scales, _ = fit_logistic_model(full_rows, FULL_FEATURES)
    reduced_model, reduced_means, reduced_scales, _ = fit_logistic_model(reduced_rows, REDUCED_FEATURES)

    predictions = manifest_frame[["sample_id_canonical"]].copy()
    predictions["prob_interpretable"] = np.nan

    full_mask = manifest_frame[FULL_FEATURES].notna().all(axis=1)
    if bool(full_mask.any()):
        predictions.loc[full_mask, "prob_interpretable"] = score_rows(
            manifest_frame.loc[full_mask, FULL_FEATURES],
            full_model,
            full_means,
            full_scales,
            FULL_FEATURES,
        )

    reduced_mask = predictions["prob_interpretable"].isna() & manifest_frame[REDUCED_FEATURES].notna().all(axis=1)
    if bool(reduced_mask.any()):
        predictions.loc[reduced_mask, "prob_interpretable"] = score_rows(
            manifest_frame.loc[reduced_mask, REDUCED_FEATURES],
            reduced_model,
            reduced_means,
            reduced_scales,
            REDUCED_FEATURES,
        )

    fallback_probability = float(sampled["prn_interpretable_numeric"].astype(float).mean())
    predictions["prob_interpretable"] = predictions["prob_interpretable"].fillna(fallback_probability)
    return predictions


def build_prevalence_from_predictions(
    manifest_frame: pd.DataFrame,
    prediction_frame: pd.DataFrame,
    *,
    weight_truncation: float,
) -> pd.DataFrame:
    frame = manifest_frame.merge(prediction_frame, on="sample_id_canonical", how="left")
    frame["country_iso3"] = frame["country_iso3"].fillna("").astype(str).str.upper()
    frame["year"] = coerce_numeric(frame["year"])
    frame["prn_disrupted_flag"] = frame.get("prn_disrupted", pd.Series(index=frame.index, dtype=str)).map(parse_bool)
    frame = frame.loc[frame["country_iso3"].ne("") & frame["year"].notna()].copy()
    frame["year"] = frame["year"].astype(int)
    frame["prob_interpretable"] = coerce_numeric(frame["prob_interpretable"]).clip(lower=1e-6, upper=1.0)
    frame["ipw_weight"] = (1.0 / frame["prob_interpretable"]).clip(upper=float(weight_truncation))
    frame["weighted_interpretable"] = np.where(frame["prn_interpretable_numeric"], frame["ipw_weight"], 0.0)
    frame["weighted_disrupted"] = np.where(
        frame["prn_interpretable_numeric"] & frame["prn_disrupted_flag"],
        frame["ipw_weight"],
        0.0,
    )

    aggregated = (
        frame.groupby(["country_iso3", "year"], dropna=False)
        .agg(
            n_genomes_total=("sample_id_canonical", "count"),
            n_genomes_prn_interpretable=("prn_interpretable_numeric", "sum"),
            n_prn_disrupted=("prn_disrupted_flag", "sum"),
            ipw_weight_total=("weighted_interpretable", "sum"),
            ipw_weighted_disrupted=("weighted_disrupted", "sum"),
        )
        .reset_index()
    )
    aggregated["naive_prevalence"] = np.where(
        aggregated["n_genomes_prn_interpretable"] > 0,
        aggregated["n_prn_disrupted"] / aggregated["n_genomes_prn_interpretable"],
        np.nan,
    )
    aggregated["ipw_prevalence"] = np.where(
        aggregated["ipw_weight_total"] > 0,
        aggregated["ipw_weighted_disrupted"] / aggregated["ipw_weight_total"],
        np.nan,
    )
    aggregated["boundary_lower_prevalence"] = np.where(
        aggregated["n_genomes_total"] > 0,
        aggregated["n_prn_disrupted"] / aggregated["n_genomes_total"],
        np.nan,
    )
    aggregated["n_missing_outcomes"] = (
        aggregated["n_genomes_total"] - aggregated["n_genomes_prn_interpretable"]
    )
    aggregated["boundary_upper_prevalence"] = np.where(
        aggregated["n_genomes_total"] > 0,
        (
            aggregated["n_prn_disrupted"] + aggregated["n_missing_outcomes"]
        ) / aggregated["n_genomes_total"],
        np.nan,
    )
    return aggregated


def logistic_probability(linear_predictor: np.ndarray) -> np.ndarray:
    clipped = np.clip(linear_predictor, -30, 30)
    return 1.0 / (1.0 + np.exp(-clipped))


def extract_programme_estimates(
    result,
    frame: pd.DataFrame,
    model_spec: str,
    *,
    n_parameter_draws: int = 200,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    draw_matrix = parameter_draws(result, n_draws=n_parameter_draws)
    term_positions = {term: index for index, term in enumerate(result.params.index)}
    for term in result.params.index:
        if "C(program_formulation_class_concurrent" not in term:
            continue
        if draw_matrix is not None and term in term_positions:
            term_values = draw_matrix[:, term_positions[term]]
        else:
            term_values = np.asarray([float(result.params[term])], dtype=float)
        for value in term_values:
            rows.append(
                {
                    "model_spec": model_spec,
                    "result_type": "coefficient",
                    "estimate_term": term,
                    "program_class": "",
                    "reference_class": "wp_only_or_pre_ap",
                    "effect_estimate": float(value),
                }
            )

    present_classes = [value for value in CLASS_ORDER if value in set(frame["program_formulation_class_concurrent"])]
    reference_class = "wp_only_or_pre_ap" if "wp_only_or_pre_ap" in present_classes else present_classes[0]
    design_info = getattr(result.model.data, "design_info", None)
    reference_draws = None
    reference_point = None
    for program_class in present_classes:
        counterfactual = frame.copy()
        counterfactual["program_formulation_class_concurrent"] = program_class
        value = float(np.mean(result.predict(counterfactual)))
        if program_class == reference_class:
            reference_point = value
            if draw_matrix is not None and design_info is not None:
                try:
                    design_matrix = build_design_matrices([design_info], counterfactual, return_type="dataframe")[0]
                    reference_linear = np.asarray(design_matrix, dtype=float) @ draw_matrix.T
                    reference_draws = logistic_probability(reference_linear).mean(axis=0)
                except Exception:
                    reference_draws = None
        if program_class == reference_class or reference_point is None:
            continue

        class_draws = None
        if draw_matrix is not None and design_info is not None:
            try:
                design_matrix = build_design_matrices([design_info], counterfactual, return_type="dataframe")[0]
                class_linear = np.asarray(design_matrix, dtype=float) @ draw_matrix.T
                class_draws = logistic_probability(class_linear).mean(axis=0)
            except Exception:
                class_draws = None

        if class_draws is not None and reference_draws is not None:
            effect_values = class_draws - reference_draws
        else:
            effect_values = np.asarray([value - reference_point], dtype=float)
        for effect_value in np.asarray(effect_values, dtype=float):
            rows.append(
                {
                    "model_spec": model_spec,
                    "result_type": "adjusted_prevalence_difference",
                    "estimate_term": f"adjusted_prevalence_difference::{program_class}::{reference_class}",
                    "program_class": program_class,
                    "reference_class": reference_class,
                    "effect_estimate": float(effect_value),
                }
            )
    return rows


def fit_bootstrap_period_models(panel: pd.DataFrame) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    frame = model_frame(
        panel,
        exposure_column="program_formulation_class_concurrent",
        response_track="ipw",
    )
    if frame.empty or "wp_only_or_pre_ap" not in set(frame["program_formulation_class_concurrent"]):
        return output

    for model_spec in ["programme_only", "programme_plus_bridge"]:
        formula = build_formula(
            "program_formulation_class_concurrent",
            include_bridge=model_spec == "programme_plus_bridge",
        )
        try:
            result, _covariance_type, _model_family = fit_response_model(
                frame,
                formula=formula,
                response_track="ipw",
                regularized=False,
            )
        except Exception:
            result, _covariance_type, _model_family = fit_response_model(
                frame,
                formula=formula,
                response_track="ipw",
                regularized=True,
            )
        output.extend(extract_programme_estimates(result, frame, model_spec))
    return output


def run_bootstrap_replicate_task(task: tuple[int, int]) -> list[dict[str, object]]:
    replicate_id, replicate_seed = task
    context = _BOOTSTRAP_CONTEXT
    if not context:
        raise RuntimeError("bootstrap context is not initialized")

    manifest_frame = context["manifest_frame"]
    exposure_frame = context["exposure_frame"]
    rng = np.random.default_rng(int(replicate_seed))
    prediction_frame = bootstrap_prediction_frame(manifest_frame, rng)
    prevalence_frame = build_prevalence_from_predictions(
        manifest_frame,
        prediction_frame,
        weight_truncation=float(context["weight_truncation"]),
    )
    annual_frame = prepare_annual_dataset_from_frames(
        exposure_frame,
        prevalence_frame,
        origin_descendants_dir=str(context["origin_descendants_dir"]),
    )
    panel = build_period_panel(
        annual_frame,
        bin_size=int(context["bin_size"]),
        min_interpretable=int(context["min_interpretable"]),
    )
    return [{"replicate_id": int(replicate_id), **row} for row in fit_bootstrap_period_models(panel)]


def summarize_bootstrap_results(replicates: pd.DataFrame, main_results: pd.DataFrame) -> pd.DataFrame:
    main_subset = main_results.loc[
        main_results["sensitivity_label"].eq("primary_period_panel")
        & main_results["response_track"].eq("ipw")
        & main_results["exposure_column"].eq("program_formulation_class_concurrent")
        & main_results["result_type"].isin(["coefficient", "adjusted_prevalence_difference"])
        & main_results["model_spec"].isin(["programme_only", "programme_plus_bridge"])
    ].copy()
    main_subset["main_model_ci_width"] = main_subset["ci_upper"] - main_subset["ci_lower"]
    join_columns = ["model_spec", "result_type", "estimate_term", "program_class", "reference_class"]
    if replicates.empty:
        summary = main_subset[join_columns + ["effect_estimate", "ci_lower", "ci_upper", "main_model_ci_width", "notes"]].rename(
            columns={
                "effect_estimate": "main_model_point_estimate",
                "ci_lower": "main_model_ci_lower",
                "ci_upper": "main_model_ci_upper",
                "notes": "main_model_notes",
            }
        ).copy()
        summary["n_replicates"] = 0
        summary["bootstrap_mean"] = np.nan
        summary["bootstrap_median"] = np.nan
        summary["bootstrap_ci_lower"] = np.nan
        summary["bootstrap_ci_upper"] = np.nan
        summary["bootstrap_interval_width"] = np.nan
        summary["propagated_ci_lower"] = np.nan
        summary["propagated_ci_upper"] = np.nan
        summary["propagated_interval_width"] = np.nan
        summary["interval_narrower_than_single_model"] = np.nan
        summary["uncertainty_summary_method"] = "bootstrap_unavailable_no_successful_replicates"
        return summary

    grouped = (
        replicates.groupby(["model_spec", "result_type", "estimate_term", "program_class", "reference_class"], dropna=False)
        .agg(
            n_replicates=("replicate_id", "nunique"),
            bootstrap_mean=("effect_estimate", "mean"),
            bootstrap_median=("effect_estimate", "median"),
            bootstrap_ci_lower=("effect_estimate", lambda values: float(np.quantile(values, 0.025))),
            bootstrap_ci_upper=("effect_estimate", lambda values: float(np.quantile(values, 0.975))),
        )
        .reset_index()
    )
    grouped["bootstrap_interval_width"] = grouped["bootstrap_ci_upper"] - grouped["bootstrap_ci_lower"]

    summary = grouped.merge(
        main_subset[join_columns + ["effect_estimate", "ci_lower", "ci_upper", "main_model_ci_width", "notes"]].rename(
            columns={
                "effect_estimate": "main_model_point_estimate",
                "ci_lower": "main_model_ci_lower",
                "ci_upper": "main_model_ci_upper",
                "notes": "main_model_notes",
            }
        ),
        on=join_columns,
        how="left",
    )
    summary["propagated_ci_lower"] = summary["bootstrap_ci_lower"]
    summary["propagated_ci_upper"] = summary["bootstrap_ci_upper"]
    summary["propagated_interval_width"] = summary["bootstrap_interval_width"]
    summary["interval_narrower_than_single_model"] = np.where(
        summary["main_model_ci_width"].notna(),
        summary["bootstrap_interval_width"] < summary["main_model_ci_width"],
        np.nan,
    )
    summary["uncertainty_summary_method"] = "bootstrap_percentile_interval"
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Run programme-surveillance two-stage uncertainty propagation.")
    parser.add_argument("--manifest", default=str(DEFAULT_WORKFLOW_ROOT / "manifest" / "manifest.tsv"))
    parser.add_argument("--assembly-qc", default=str(DEFAULT_WORKFLOW_ROOT / "assembly_qc" / "assembly_qc_stats.tsv"))
    parser.add_argument("--exposure", default=str(DEFAULT_WORKFLOW_ROOT / "epi" / "ap_exposure_index.tsv"))
    parser.add_argument(
        "--main-results",
        default=str(DEFAULT_WORKFLOW_ROOT / "epi" / "programme_program_model_results.tsv"),
    )
    parser.add_argument(
        "--origin-descendants-dir",
        default=str(DEFAULT_WORKFLOW_ROOT / "asr" / "event_subtrees"),
    )
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--replicates-out", required=True)
    parser.add_argument("--n-replicates", type=int, default=80)
    parser.add_argument("--weight-truncation", type=float, default=20.0)
    parser.add_argument("--bin-size", type=int, default=5)
    parser.add_argument("--min-interpretable", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260409)
    parser.add_argument(
        "--max-workers",
        type=int,
        default=None,
        help=(
            "Maximum process workers for bootstrap replicates. "
            f"Defaults to the {PROGRAMME_UNCERTAINTY_MAX_WORKERS_ENV} env var or an auto cap of "
            f"{PROGRAMME_UNCERTAINTY_DEFAULT_MAX_WORKERS}."
        ),
    )
    args = parser.parse_args()

    manifest_frame = load_inputs(args.manifest, args.assembly_qc)
    exposure_frame = pd.read_csv(args.exposure, sep="\t", dtype=str)
    main_results = pd.read_csv(args.main_results, sep="\t", dtype=str)
    for column in ["ci_lower", "ci_upper", "effect_estimate"]:
        main_results[column] = coerce_numeric(main_results.get(column, pd.Series(dtype=str)))

    requested_max_workers = args.max_workers
    if requested_max_workers is None:
        requested_max_workers = resolve_env_max_workers()

    set_bootstrap_context(
        manifest_frame,
        exposure_frame,
        origin_descendants_dir=args.origin_descendants_dir,
        weight_truncation=args.weight_truncation,
        bin_size=args.bin_size,
        min_interpretable=args.min_interpretable,
    )
    replicate_tasks = list(enumerate(build_replicate_seeds(args.n_replicates, args.seed), start=1))
    worker_count = resolve_parallel_workers(len(replicate_tasks), requested_max_workers=requested_max_workers)

    replicate_rows: list[dict[str, object]] = []
    if worker_count <= 1:
        for task in replicate_tasks:
            replicate_rows.extend(run_bootstrap_replicate_task(task))
    elif replicate_tasks:
        with ProcessPoolExecutor(
            max_workers=worker_count,
            initializer=_init_bootstrap_worker,
            initargs=(
                manifest_frame,
                exposure_frame,
                args.origin_descendants_dir,
                args.weight_truncation,
                args.bin_size,
                args.min_interpretable,
            ),
        ) as pool:
            for task_rows in pool.map(run_bootstrap_replicate_task, replicate_tasks):
                replicate_rows.extend(task_rows)
    else:
        replicate_rows = []

    replicates = pd.DataFrame.from_records(replicate_rows)
    summary = summarize_bootstrap_results(replicates, main_results)

    summary_path = Path(args.summary_out)
    replicates_path = Path(args.replicates_out)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    replicates_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, sep="\t", index=False)
    replicates.to_csv(replicates_path, sep="\t", index=False)
    print(f"Wrote programme-surveillance two-stage uncertainty summary: {summary_path}")
    print(f"Wrote programme-surveillance two-stage uncertainty replicates: {replicates_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
