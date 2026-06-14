#!/usr/bin/env python3
"""Fit a workflow-native missingness model for prn interpretability."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm


FULL_FEATURES = ["year_numeric", "has_reads_numeric", "log_total_length", "log_n_contigs"]
REDUCED_FEATURES = ["year_numeric", "has_reads_numeric"]


def parse_bool(value: object) -> bool:
    text = str(value).strip().lower()
    return text in {"true", "1", "yes"}


def coerce_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def normalize_qc_columns(qc_frame: pd.DataFrame) -> pd.DataFrame:
    rename_map = {
        "sample": "sample_id_canonical",
        "sample_id": "sample_id_canonical",
        "total length": "total_length",
        "total_length": "total_length",
        "n_contigs": "n_contigs",
        "contigs": "n_contigs",
        "contig_n50": "contig_n50",
        "n50": "contig_n50",
    }
    available = {column: rename_map[column] for column in qc_frame.columns if column in rename_map}
    return qc_frame.rename(columns=available)


def roc_auc_score_manual(target: pd.Series, probabilities: np.ndarray) -> float:
    target_array = target.astype(int).to_numpy()
    positive = int(target_array.sum())
    negative = int(len(target_array) - positive)
    if positive == 0 or negative == 0:
        return float("nan")
    ranks = pd.Series(probabilities).rank(method="average").to_numpy()
    auc = (ranks[target_array == 1].sum() - positive * (positive + 1) / 2) / (positive * negative)
    return float(auc)


def brier_score_manual(target: pd.Series, probabilities: np.ndarray) -> float:
    target_array = target.astype(int).to_numpy(dtype=float)
    probability_array = np.asarray(probabilities, dtype=float)
    if len(target_array) == 0 or len(target_array) != len(probability_array):
        return float("nan")
    return float(np.mean(np.square(target_array - probability_array)))


def make_stratified_folds(target: pd.Series, n_splits: int, random_seed: int = 20260506) -> list[tuple[np.ndarray, np.ndarray]]:
    if n_splits < 2:
        return []
    target_array = target.astype(int).to_numpy()
    rng = np.random.default_rng(random_seed)
    fold_ids = np.full(len(target_array), -1, dtype=int)
    for class_value in np.unique(target_array):
        class_idx = np.flatnonzero(target_array == class_value)
        rng.shuffle(class_idx)
        for offset, row_idx in enumerate(class_idx):
            fold_ids[row_idx] = offset % n_splits
    all_idx = np.arange(len(target_array))
    folds: list[tuple[np.ndarray, np.ndarray]] = []
    for fold_id in range(n_splits):
        test_idx = all_idx[fold_ids == fold_id]
        train_idx = all_idx[fold_ids != fold_id]
        if len(test_idx) == 0 or len(train_idx) == 0:
            continue
        folds.append((train_idx, test_idx))
    return folds


def load_inputs(manifest_path: str, assembly_qc_path: str | None = None) -> pd.DataFrame:
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    manifest["sample_id_canonical"] = manifest["sample_id_canonical"].astype(str)
    manifest["year_numeric"] = coerce_numeric(manifest.get("year", pd.Series(index=manifest.index, dtype=str)))
    manifest["has_reads_numeric"] = manifest.get("has_reads", pd.Series(index=manifest.index, dtype=str)).map(parse_bool).astype(int)
    manifest["prn_interpretable_numeric"] = (
        manifest.get("prn_interpretable", pd.Series(index=manifest.index, dtype=str)).map(parse_bool)
    )

    manifest["total_length"] = coerce_numeric(
        manifest.get("total_sequence_length", pd.Series(index=manifest.index, dtype=str))
    )
    manifest["n_contigs_numeric"] = coerce_numeric(
        manifest.get("n_contigs", pd.Series(index=manifest.index, dtype=str))
    )

    if assembly_qc_path and Path(assembly_qc_path).exists():
        qc_frame = pd.read_csv(assembly_qc_path, sep="\t", dtype=str)
        qc_frame = normalize_qc_columns(qc_frame)
        if "sample_id_canonical" in qc_frame.columns:
            keep_columns = [
                column
                for column in ["sample_id_canonical", "total_length", "n_contigs", "contig_n50"]
                if column in qc_frame.columns
            ]
            qc_frame = qc_frame[keep_columns].drop_duplicates(subset=["sample_id_canonical"]).copy()
            rename_qc = {}
            if "total_length" in qc_frame.columns:
                rename_qc["total_length"] = "total_length_qc"
            if "n_contigs" in qc_frame.columns:
                rename_qc["n_contigs"] = "n_contigs_qc"
            qc_frame = qc_frame.rename(columns=rename_qc)
            qc_frame["total_length_qc"] = coerce_numeric(qc_frame.get("total_length_qc", pd.Series(dtype=str)))
            qc_frame["n_contigs_qc"] = coerce_numeric(qc_frame.get("n_contigs_qc", pd.Series(dtype=str)))
            manifest = manifest.merge(qc_frame, on="sample_id_canonical", how="left")
            manifest["total_length"] = manifest["total_length"].fillna(manifest.get("total_length_qc"))
            manifest["n_contigs_numeric"] = manifest["n_contigs_numeric"].fillna(manifest.get("n_contigs_qc"))

    manifest["log_total_length"] = np.log1p(manifest["total_length"])
    manifest["log_n_contigs"] = np.log1p(manifest["n_contigs_numeric"])
    return manifest


def fit_logistic_model(frame: pd.DataFrame, feature_columns: list[str]) -> tuple[object, pd.Series, pd.Series, dict[str, float]]:
    feature_frame = frame[feature_columns].astype(float)
    feature_means = feature_frame.mean()
    feature_scales = feature_frame.std(ddof=0).replace(0, 1.0).fillna(1.0)
    standardized = (feature_frame - feature_means) / feature_scales
    design = sm.add_constant(standardized, has_constant="add")
    target = frame["prn_interpretable_numeric"].astype(int)
    model = sm.GLM(target, design, family=sm.families.Binomial())
    fit_method = "glm"
    try:
        fitted = model.fit(maxiter=200)
    except Exception:
        fitted = model.fit_regularized(alpha=1e-6, L1_wt=0.0, maxiter=200)
        fit_method = "glm_regularized"
    probabilities = np.asarray(fitted.predict(design), dtype=float)
    predictions = (probabilities >= 0.5).astype(int)
    metrics = {
        "n_rows": int(len(frame)),
        "training_accuracy": float((predictions == target.to_numpy()).mean()),
        "training_auc": roc_auc_score_manual(target, probabilities),
        "training_brier": brier_score_manual(target, probabilities),
        "mean_probability": float(np.mean(probabilities)),
        "fit_method": fit_method,
    }
    return fitted, feature_means, feature_scales, metrics


def cross_validated_probabilities(frame: pd.DataFrame, feature_columns: list[str]) -> np.ndarray:
    target = frame["prn_interpretable_numeric"].astype(int)
    if len(frame) < 4 or target.nunique() < 2:
        return np.repeat(np.nan, len(frame))
    min_class = int(target.value_counts().min())
    n_splits = min(5, min_class)
    if n_splits < 2:
        return np.repeat(np.nan, len(frame))
    output = np.repeat(np.nan, len(frame))
    for train_idx, test_idx in make_stratified_folds(target, n_splits=n_splits):
        train_frame = frame.iloc[train_idx].copy()
        test_frame = frame.iloc[test_idx].copy()
        try:
            fitted, feature_means, feature_scales, _metrics = fit_logistic_model(train_frame, feature_columns)
            output[test_idx] = score_rows(test_frame, fitted, feature_means, feature_scales, feature_columns)
        except Exception:
            return np.repeat(np.nan, len(frame))
    return output


def summarise_probability_metrics(target: pd.Series, probabilities: np.ndarray) -> dict[str, float]:
    target_series = pd.Series(target, dtype=float)
    probability_array = np.asarray(probabilities, dtype=float)
    mask = target_series.notna() & np.isfinite(probability_array)
    if not mask.any():
        return {
            "n_rows": 0,
            "accuracy": float("nan"),
            "auc": float("nan"),
            "brier": float("nan"),
        }
    observed = target_series.loc[mask].astype(int)
    predicted = probability_array[mask.to_numpy()]
    labels = (predicted >= 0.5).astype(int)
    return {
        "n_rows": int(mask.sum()),
        "accuracy": float((labels == observed.to_numpy()).mean()),
        "auc": roc_auc_score_manual(observed, predicted),
        "brier": brier_score_manual(observed, predicted),
    }


def serialize_model(
    model,
    feature_means: pd.Series,
    feature_scales: pd.Series,
    feature_columns: list[str],
    metrics: dict[str, float],
    out_of_fold_metrics: dict[str, float],
) -> dict[str, object]:
    coefficients = []
    for feature_name in feature_columns:
        coefficients.append(
            {
                "feature": feature_name,
                "coefficient": float(model.params[feature_name]),
                "feature_mean": float(feature_means[feature_name]),
                "feature_scale": float(feature_scales[feature_name]),
            }
        )
    return {
        "features": feature_columns,
        "intercept": float(model.params["const"]),
        "metrics": metrics,
        "out_of_fold_metrics": out_of_fold_metrics,
        "coefficients": coefficients,
    }


def score_rows(frame: pd.DataFrame, model, feature_means: pd.Series, feature_scales: pd.Series, feature_columns: list[str]) -> np.ndarray:
    standardized = (frame[feature_columns].astype(float) - feature_means) / feature_scales
    design = sm.add_constant(standardized, has_constant="add")
    return np.asarray(model.predict(design), dtype=float)


def build_predictions(frame: pd.DataFrame) -> tuple[dict[str, object], pd.DataFrame]:
    labeled = frame.dropna(subset=["prn_interpretable_numeric"]).copy()
    full_rows = labeled.dropna(subset=FULL_FEATURES).copy()
    reduced_rows = labeled.dropna(subset=REDUCED_FEATURES).copy()

    full_model, full_means, full_scales, full_metrics = fit_logistic_model(full_rows, FULL_FEATURES)
    reduced_model, reduced_means, reduced_scales, reduced_metrics = fit_logistic_model(reduced_rows, REDUCED_FEATURES)
    full_oof = cross_validated_probabilities(full_rows, FULL_FEATURES)
    reduced_oof = cross_validated_probabilities(reduced_rows, REDUCED_FEATURES)
    full_oof_metrics = summarise_probability_metrics(full_rows["prn_interpretable_numeric"], full_oof)
    reduced_oof_metrics = summarise_probability_metrics(reduced_rows["prn_interpretable_numeric"], reduced_oof)

    predictions = labeled[["sample_id_canonical", "prn_interpretable_numeric"]].copy()
    predictions["prob_interpretable"] = np.nan
    predictions["prob_interpretable_oof"] = np.nan
    predictions["prediction_source"] = "fallback"
    predictions["oof_prediction_source"] = "not_available"

    full_index = full_rows.index
    predictions.loc[full_index, "prob_interpretable"] = score_rows(full_rows, full_model, full_means, full_scales, FULL_FEATURES)
    predictions.loc[full_index, "prob_interpretable_oof"] = full_oof
    predictions.loc[full_index, "prediction_source"] = "full_model"
    predictions.loc[full_index, "oof_prediction_source"] = "full_model_oof"

    reduced_only = predictions["prob_interpretable"].isna() & labeled[REDUCED_FEATURES].notna().all(axis=1)
    reduced_index = labeled.index[reduced_only]
    if len(reduced_index) > 0:
        reduced_subset = labeled.loc[reduced_index, REDUCED_FEATURES]
        predictions.loc[reduced_index, "prob_interpretable"] = score_rows(
            reduced_subset,
            reduced_model,
            reduced_means,
            reduced_scales,
            REDUCED_FEATURES,
        )
        predictions.loc[reduced_index, "prediction_source"] = "reduced_model"
        reduced_oof_series = pd.Series(reduced_oof, index=reduced_rows.index, dtype=float)
        predictions.loc[reduced_index, "prob_interpretable_oof"] = reduced_oof_series.reindex(reduced_index).to_numpy()
        predictions.loc[reduced_index, "oof_prediction_source"] = "reduced_model_oof"

    fallback_probability = float(labeled["prn_interpretable_numeric"].astype(float).mean())
    predictions["prob_interpretable"] = predictions["prob_interpretable"].fillna(fallback_probability)
    predictions["prediction_source"] = predictions["prediction_source"].fillna("fallback")
    predictions["oof_prediction_source"] = predictions["oof_prediction_source"].fillna("not_available")
    predictions["pred"] = (predictions["prob_interpretable"] >= 0.5).astype(int)
    predictions["pred_oof"] = np.where(
        predictions["prob_interpretable_oof"].notna(),
        (predictions["prob_interpretable_oof"] >= 0.5).astype(int),
        np.nan,
    )
    predictions["prediction_probability_scope"] = np.where(
        predictions["prediction_source"].eq("fallback"),
        "fallback_marginal_probability",
        "training_scope_fitted_probability",
    )
    predictions["oof_probability_scope"] = np.where(
        predictions["prob_interpretable_oof"].notna(),
        "out_of_fold_probability",
        "not_available",
    )
    predictions["in_model"] = predictions["prediction_source"].isin(["full_model", "reduced_model"])
    predictions["y_actual"] = predictions["prn_interpretable_numeric"].astype(int)

    metadata = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "model_version": "workflow_missingness_v2",
        "n_total_manifest_rows": int(len(frame)),
        "n_labeled_rows": int(len(labeled)),
        "fallback_probability": fallback_probability,
        "full_model": serialize_model(
            full_model,
            full_means,
            full_scales,
            FULL_FEATURES,
            full_metrics,
            full_oof_metrics,
        ),
        "reduced_model": serialize_model(
            reduced_model,
            reduced_means,
            reduced_scales,
            REDUCED_FEATURES,
            reduced_metrics,
            reduced_oof_metrics,
        ),
        "prediction_summary": {
            "n_full_model": int((predictions["prediction_source"] == "full_model").sum()),
            "n_reduced_model": int((predictions["prediction_source"] == "reduced_model").sum()),
            "n_fallback": int((predictions["prediction_source"] == "fallback").sum()),
            "n_out_of_fold_predictions": int(predictions["prob_interpretable_oof"].notna().sum()),
            "mean_probability": float(predictions["prob_interpretable"].mean()),
            "min_probability": float(predictions["prob_interpretable"].min()),
            "max_probability": float(predictions["prob_interpretable"].max()),
            "training_probability_scope": "training_scope_fitted_probability",
            "out_of_fold_probability_scope": "out_of_fold_probability_when_available",
        },
    }
    return metadata, predictions


def build_report_html(metadata: dict[str, object], predictions: pd.DataFrame) -> str:
    summary_rows = pd.DataFrame(
        [
            {
                "section": "full_model",
                "n_rows": metadata["full_model"]["metrics"]["n_rows"],
                "training_accuracy": metadata["full_model"]["metrics"]["training_accuracy"],
                "training_auc": metadata["full_model"]["metrics"]["training_auc"],
                "oof_accuracy": metadata["full_model"]["out_of_fold_metrics"]["accuracy"],
                "oof_auc": metadata["full_model"]["out_of_fold_metrics"]["auc"],
            },
            {
                "section": "reduced_model",
                "n_rows": metadata["reduced_model"]["metrics"]["n_rows"],
                "training_accuracy": metadata["reduced_model"]["metrics"]["training_accuracy"],
                "training_auc": metadata["reduced_model"]["metrics"]["training_auc"],
                "oof_accuracy": metadata["reduced_model"]["out_of_fold_metrics"]["accuracy"],
                "oof_auc": metadata["reduced_model"]["out_of_fold_metrics"]["auc"],
            },
        ]
    )
    coefficient_rows = []
    for model_name in ["full_model", "reduced_model"]:
        for row in metadata[model_name]["coefficients"]:
            coefficient_rows.append({"model": model_name, **row})
    coefficient_frame = pd.DataFrame(coefficient_rows)
    top_predictions = predictions.sort_values("prob_interpretable", ascending=False).head(20)

    return "".join(
        [
            "<html><head><meta charset='utf-8'><title>Missingness diagnostics</title></head><body>",
            "<h1>Workflow Missingness Diagnostics</h1>",
            f"<p>Generated at {metadata['generated_at']}</p>",
            pd.DataFrame([metadata["prediction_summary"]]).to_html(index=False),
            "<h2>Model metrics</h2>",
            summary_rows.to_html(index=False),
            "<h2>Coefficients</h2>",
            coefficient_frame.to_html(index=False),
            "<h2>Highest predicted interpretability</h2>",
            top_predictions.to_html(index=False),
            "</body></html>",
        ]
    )


def write_sidecar_outputs(model_path: Path, metadata: dict[str, object], predictions: pd.DataFrame) -> None:
    output_dir = model_path.parent
    predictions_out = output_dir / "missingness_model_predictions.tsv"
    coefficients_out = output_dir / "missingness_model_coefficients.tsv"
    summary_out = output_dir / "missingness_model_summary.txt"

    prediction_frame = predictions[
        [
            "sample_id_canonical",
            "y_actual",
            "prob_interpretable",
            "prob_interpretable_oof",
            "pred",
            "pred_oof",
            "in_model",
            "prediction_source",
            "oof_prediction_source",
            "prediction_probability_scope",
            "oof_probability_scope",
        ]
    ].copy()
    prediction_frame.to_csv(predictions_out, sep="\t", index=False)

    coefficient_rows = []
    for model_name in ["full_model", "reduced_model"]:
        model_payload = metadata[model_name]
        for row in model_payload["coefficients"]:
            coefficient_rows.append(
                {
                    "model": model_name,
                    "feature": row["feature"],
                    "coefficient": row["coefficient"],
                    "feature_mean": row["feature_mean"],
                    "feature_scale": row["feature_scale"],
                    "fit_method": model_payload["metrics"]["fit_method"],
                    "training_metric_provenance": "training_scope_fitted_probabilities",
                    "training_accuracy": model_payload["metrics"]["training_accuracy"],
                    "training_auc": model_payload["metrics"]["training_auc"],
                    "training_brier": model_payload["metrics"]["training_brier"],
                    "oof_metric_provenance": "stratified_out_of_fold_probabilities",
                    "oof_accuracy": model_payload["out_of_fold_metrics"]["accuracy"],
                    "oof_auc": model_payload["out_of_fold_metrics"]["auc"],
                    "oof_brier": model_payload["out_of_fold_metrics"]["brier"],
                }
            )
    pd.DataFrame(coefficient_rows).to_csv(coefficients_out, sep="\t", index=False)

    summary_lines = [
        "=== Workflow Missingness Model Summary ===",
        f"Total samples with labels: {metadata['n_labeled_rows']}",
        f"Model samples (complete features): {metadata['full_model']['metrics']['n_rows']}",
        f"Accuracy: {metadata['full_model']['metrics']['training_accuracy']:.4f}",
        f"AUC: {metadata['full_model']['metrics']['training_auc']:.4f}",
        f"Out-of-fold Accuracy: {metadata['full_model']['out_of_fold_metrics']['accuracy']:.4f}"
        if np.isfinite(metadata["full_model"]["out_of_fold_metrics"]["accuracy"])
        else "Out-of-fold Accuracy: nan",
        f"Out-of-fold AUC: {metadata['full_model']['out_of_fold_metrics']['auc']:.4f}"
        if np.isfinite(metadata["full_model"]["out_of_fold_metrics"]["auc"])
        else "Out-of-fold AUC: nan",
    ]
    summary_out.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")


def run_model(manifest_path: str, model_out: str, report_out: str, assembly_qc_path: str | None = None) -> None:
    frame = load_inputs(manifest_path, assembly_qc_path)
    metadata, predictions = build_predictions(frame)
    metadata["predictions"] = predictions.to_dict(orient="records")

    model_path = Path(model_out)
    report_path = Path(report_out)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    with model_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
    report_path.write_text(build_report_html(metadata, predictions), encoding="utf-8")
    write_sidecar_outputs(model_path, metadata, predictions)


if "snakemake" in globals():
    run_model(
        manifest_path=snakemake.input.manifest,
        assembly_qc_path=snakemake.input.get("assembly_qc"),
        model_out=snakemake.output.model,
        report_out=snakemake.output.report,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Fit workflow missingness model")
    parser.add_argument("--manifest", required=True, help="Manifest TSV")
    parser.add_argument("--assembly-qc", default="", help="Optional aggregated assembly QC TSV")
    parser.add_argument("--model-out", required=True, help="Output model JSON")
    parser.add_argument("--report-out", required=True, help="Output HTML report")
    arguments = parser.parse_args()

    run_model(
        manifest_path=arguments.manifest,
        assembly_qc_path=arguments.assembly_qc or None,
        model_out=arguments.model_out,
        report_out=arguments.report_out,
    )
