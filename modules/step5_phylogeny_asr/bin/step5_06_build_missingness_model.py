#!/usr/bin/env python3
"""T10 — PRN Missingness Model: logistic model predicting prn interpretability.

Builds a logistic regression predicting prn_interpretable ~ f(assembly quality, metadata).
Used to:
  1. Justify non-random missingness assumption in ASR
  2. Generate imputation probabilities for "model_impute" ASR scenario
  3. Identify borderline samples for re-examination

Inputs:
    workflow/manifest/manifest.tsv
    workflow/assembly_qc/assembly_qc_stats.tsv (optional, enriches model)

Outputs:
    workflow/missingness_model/missingness_model_coefficients.tsv
    workflow/missingness_model/missingness_model_predictions.tsv
    workflow/missingness_model/missingness_model_summary.txt
    workflow/missingness_model/missingness_model_diagnostics.png

Usage:
    python modules/step5_phylogeny_asr/bin/step5_06_build_missingness_model.py
"""

import csv
import os
import sys
import warnings
from pathlib import Path
from collections import Counter

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from workflow.lib.project_paths import project_workflow_root


WORKFLOW_ROOT = project_workflow_root()
MANIFEST = Path(__file__).resolve().parents[3] / "state" / "manifest" / "manifest.tsv"
QC_STATS = WORKFLOW_ROOT / "assembly_qc" / "assembly_qc_stats.tsv"
OUTPUT_DIR = WORKFLOW_ROOT / "missingness_model"


def safe_float(val, default=np.nan):
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def safe_int(val, default=0):
    try:
        return int(val)
    except (ValueError, TypeError):
        return default


def load_data():
    """Load manifest + optional QC data, return structured arrays.
    
    Uses only universally-available features (year, has_reads, assembly size,
    n_contigs) since prn-specific features are by definition absent for
    non-interpretable samples.
    """
    # Load manifest
    with open(MANIFEST, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        manifest = list(reader)

    # Load QC if available — enriches total_length / n_contigs
    qc_map = {}
    if QC_STATS.exists():
        with open(QC_STATS, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row.get("qc_status") not in ("NO_FASTA", "FAIL"):
                    qc_map[row["sample_id_canonical"]] = row

    records = []
    for row in manifest:
        sid = row["sample_id_canonical"]
        qc = qc_map.get(sid, {})

        # Response variable
        interp_raw = row.get("prn_interpretable", "")
        if interp_raw not in ("True", "False", "1", "0", "true", "false"):
            continue
        y = 1 if interp_raw in ("True", "1", "true") else 0

        # Universal predictors (available for both groups)
        year = safe_float(row.get("year", ""))
        has_reads = 1 if row.get("has_reads", "") in ("True", "1", "true") else 0
        
        # Assembly stats: prefer QC-computed, fallback to manifest metadata
        total_length = safe_float(qc.get("total_length", row.get("total_sequence_length", "")))
        n_contigs = safe_float(qc.get("n_contigs", row.get("n_contigs", "")))

        records.append({
            "sample_id_canonical": sid,
            "y": y,
            "year": year,
            "has_reads": has_reads,
            "log_total_length": np.log10(total_length) if not np.isnan(total_length) and total_length > 0 else np.nan,
            "log_n_contigs": np.log10(n_contigs) if not np.isnan(n_contigs) and n_contigs > 0 else np.nan,
        })

    return records


def fit_logistic(X, y, max_iter=100, lr=0.01):
    """Simple logistic regression via gradient descent (no sklearn dependency).
    
    X: (n, p) array with intercept column
    y: (n,) binary array
    Returns: coefficients (p,), convergence info
    """
    n, p = X.shape
    beta = np.zeros(p)

    for iteration in range(max_iter):
        z = X @ beta
        # Clip for numerical stability
        z = np.clip(z, -20, 20)
        prob = 1.0 / (1.0 + np.exp(-z))
        grad = X.T @ (prob - y) / n
        beta -= lr * grad

        if np.max(np.abs(grad)) < 1e-6:
            break

    return beta, iteration + 1


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    records = load_data()
    print(f"Loaded {len(records)} samples with interpretability labels")

    # Select features — tiered: full model if assembly stats available, else reduced
    full_features = ["year", "has_reads", "log_total_length", "log_n_contigs"]
    reduced_features = ["year", "has_reads"]

    # Try full model first
    feature_names = full_features
    valid = []
    for r in records:
        vals = [r[f] for f in feature_names]
        if not any(np.isnan(v) for v in vals):
            valid.append(r)

    # Check both classes present
    y_check = set(r["y"] for r in valid)
    if len(valid) < 50 or len(y_check) < 2:
        print(f"Full model: {len(valid)} samples, {len(y_check)} classes — falling back to reduced.")
        feature_names = reduced_features
        valid = []
        for r in records:
            vals = [r.get(f, np.nan) for f in feature_names]
            if not any(np.isnan(v) for v in vals):
                valid.append(r)
        y_check = set(r["y"] for r in valid)
        print(f"Reduced model: {len(valid)} samples, {len(y_check)} classes")

    if len(valid) < 20 or len(y_check) < 2:
        print("ERROR: Insufficient data for model fitting (need >=20 samples with both classes)", file=sys.stderr)
        sys.exit(1)

    y = np.array([r["y"] for r in valid])
    X_raw = np.array([[r[f] for f in feature_names] for r in valid])

    # Standardize features
    means = np.nanmean(X_raw, axis=0)
    stds = np.nanstd(X_raw, axis=0)
    stds[stds == 0] = 1
    X_std = (X_raw - means) / stds

    # Add intercept
    X = np.column_stack([np.ones(len(valid)), X_std])
    names_with_intercept = ["intercept"] + feature_names

    # Fit
    beta, n_iter = fit_logistic(X, y, max_iter=500, lr=0.5)

    # Predictions
    z = X @ beta
    z = np.clip(z, -20, 20)
    probs = 1.0 / (1.0 + np.exp(-z))

    # Accuracy
    preds = (probs >= 0.5).astype(int)
    accuracy = np.mean(preds == y)
    baseline = max(np.mean(y), 1 - np.mean(y))

    # AUC (simple trapezoidal)
    sorted_idx = np.argsort(-probs)
    y_sorted = y[sorted_idx]
    n_pos = np.sum(y)
    n_neg = np.sum(1 - y)
    if n_pos > 0 and n_neg > 0:
        tpr_list = np.cumsum(y_sorted) / n_pos
        fpr_list = np.cumsum(1 - y_sorted) / n_neg
        auc = float(np.trapezoid(tpr_list, fpr_list)) if hasattr(np, 'trapezoid') else float(np.trapz(tpr_list, fpr_list))
    else:
        auc = float("nan")

    # Output coefficients
    coef_path = OUTPUT_DIR / "missingness_model_coefficients.tsv"
    with open(coef_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["feature", "coefficient_std", "mean", "std"])
        for i, name in enumerate(names_with_intercept):
            m = means[i - 1] if i > 0 else 0
            s = stds[i - 1] if i > 0 else 1
            writer.writerow([name, f"{beta[i]:.4f}", f"{m:.4f}", f"{s:.4f}"])

    # Output predictions for all records (including invalid ones with NaN prob)
    pred_path = OUTPUT_DIR / "missingness_model_predictions.tsv"
    valid_set = {r["sample_id_canonical"] for r in valid}
    prob_map = {valid[i]["sample_id_canonical"]: probs[i] for i in range(len(valid))}

    with open(pred_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id_canonical", "y_actual", "prob_interpretable", "pred", "in_model"])
        for r in records:
            sid = r["sample_id_canonical"]
            p = prob_map.get(sid, "")
            pred = ""
            if p != "":
                pred = 1 if p >= 0.5 else 0
                p = f"{p:.4f}"
            writer.writerow([sid, r["y"], p, pred, sid in valid_set])

    # Summary
    summary_lines = [
        "=== PRN Missingness Model Summary ===",
        f"Total samples with labels: {len(records)}",
        f"Interpretable (y=1): {sum(1 for r in records if r['y']==1)}",
        f"Not interpretable (y=0): {sum(1 for r in records if r['y']==0)}",
        f"",
        f"Model samples (complete features): {len(valid)}",
        f"Features: {', '.join(feature_names)}",
        f"Iterations to converge: {n_iter}",
        f"",
        f"=== Performance ===",
        f"Accuracy: {accuracy:.3f}  (baseline: {baseline:.3f})",
        f"AUC: {auc:.3f}",
        f"",
        f"=== Coefficients (standardized) ===",
    ]
    for i, name in enumerate(names_with_intercept):
        summary_lines.append(f"  {name:25s}  {beta[i]:+.4f}")

    summary_lines.extend([
        "",
        "=== Interpretation ===",
        "Positive coefficients → higher prob of interpretable call",
        "Model used for ASR 'model_impute' scenario: P(interpretable) per sample",
    ])

    summary_path = OUTPUT_DIR / "missingness_model_summary.txt"
    with open(summary_path, "w") as f:
        f.write("\n".join(summary_lines) + "\n")

    print("\n".join(summary_lines))
    print(f"\nOutputs:")
    print(f"  Coefficients: {coef_path}")
    print(f"  Predictions: {pred_path}")
    print(f"  Summary: {summary_path}")


if __name__ == "__main__":
    main()
