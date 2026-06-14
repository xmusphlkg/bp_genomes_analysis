#!/usr/bin/env python3
"""Build manuscript-facing sidecars for the April 13/14 analysis upgrade pass.

This script adds reviewer-facing analysis layers that can be cited directly in
the manuscript package without re-running the full workflow:
1. Selected-country year-composition and leave-one-year-out sensitivity.
2. Country-by-epoch *prn*-disruption architecture turnover.
3. Local origin burden versus prevalence-shift bridge.
4. PRN-locus structural signal specificity audit across same-frame antigens and pseudo-control loci.
5. Genome-level PRN interpretability model summary table.
6. Quality-restricted selected-country sensitivity repeats.
7. Missingness tipping-point stress-test summaries.
8. Threshold-robust selected-country readiness screening summaries.

The goal is not to replace the main workflow, but to freeze manuscript-facing TSVs
that can be rendered into Extended Data figures and cited directly in the draft.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import statsmodels.api as sm


REPO_ROOT = Path(__file__).resolve().parents[3]
FIGURE_DATA_DIR = REPO_ROOT / "manuscript" / "figure_data"
SELECTED_DIR = FIGURE_DATA_DIR / "selected_country"
SUPP_DIR = REPO_ROOT / "manuscript" / "supplementary"
SELECTED_COUNTRY_CURATION_DIR = REPO_ROOT / "manuscript" / "curation" / "selected_country"
AUDIT_SUPP_TABLE_DIR = REPO_ROOT / "manuscript" / "submission_data" / "audit_ledgers" / "supplementary_table_sources"

HISTORY_PATH = SELECTED_DIR / "country_program_history_manifest.tsv"
EPOCH_PREVALENCE_PATH = SELECTED_DIR / "country_epoch_prn_prevalence.tsv"
EPOCH_ELIGIBILITY_PATH = SELECTED_DIR / "country_epoch_eligibility.tsv"
SELECTION_SCORECARD_PATH = SELECTED_DIR / "country_selection_scorecard.tsv"
YEAR_COMPOSITION_PATH = SUPP_DIR / "Supplementary_Table_34_Selected_Country_Epoch_Year_Composition.tsv"
YEAR_COMPOSITION_FALLBACK_PATH = AUDIT_SUPP_TABLE_DIR / "Supplementary_Table_34_Selected_Country_Epoch_Year_Composition.tsv"
EPOCH_CONTRAST_PATH = SELECTED_DIR / "country_epoch_contrast_summary.tsv"
ORIGIN_PACKAGE_PATH = SELECTED_DIR / "selected_country_origin_package_summary.tsv"
ORIGIN_SHIFT_PATH = SELECTED_DIR / "selected_country_origin_amplification.tsv"
DETECTION_SHIFT_PATH = SELECTED_DIR / "selected_country_detection_amplification.tsv"
STRUCTURE_REUSE_PATH = SELECTED_DIR / "selected_country_structure_reuse.tsv"
EVIDENCE_GRID_PATH = SELECTED_DIR / "cross_country_evidence_grid.tsv"

IPW_PATH = REPO_ROOT / "outputs" / "workflow" / "epi" / "ipw_prevalence.tsv"
MANIFEST_PATH = REPO_ROOT / "outputs" / "workflow" / "manifest" / "manifest.tsv"
ASSEMBLY_QC_PATH = REPO_ROOT / "outputs" / "workflow" / "assembly_qc" / "assembly_qc_stats.tsv"
STEP2_MARKERS_PATH = FIGURE_DATA_DIR / "pseudo_control_loci" / "pseudo_control_marker_status.tsv"
PSEUDO_CONTROL_STATUS_PATH = FIGURE_DATA_DIR / "pseudo_control_loci" / "pseudo_control_marker_status.tsv"
RESCUED_MANIFEST_OVERRIDES_PATH = SELECTED_COUNTRY_CURATION_DIR / "rescued_prn_overrides.tsv"
MISSINGNESS_MODEL_SUMMARY_PATH = REPO_ROOT / "outputs" / "workflow" / "missingness_model" / "missingness_model_summary.txt"
MISSINGNESS_MODEL_COEFFICIENTS_PATH = REPO_ROOT / "outputs" / "workflow" / "missingness_model" / "missingness_model_coefficients.tsv"
ASR_ROOTING_SENSITIVITY_PATH = REPO_ROOT / "manuscript" / "figure_data" / "asr_rooting_sensitivity.tsv"
ASR_SENSITIVITY_PATH = REPO_ROOT / "manuscript" / "figure_data" / "figure3_workflow_asr_sensitivity.tsv"
ASR_MK_PATH = REPO_ROOT / "manuscript" / "figure_data" / "asr_mk_origin_uncertainty.tsv"
ASR_RESAMPLING_DIR = REPO_ROOT / "outputs" / "workflow" / "asr_resampling"
PRIMARY_ORIGIN_EVENTS_PATH = FIGURE_DATA_DIR / "figure3_workflow_origin_events.tsv"

YEAR_SENS_SUMMARY_PATH = SELECTED_DIR / "selected_country_year_sensitivity_summary.tsv"
YEAR_SENS_LOYO_PATH = SELECTED_DIR / "selected_country_leave_one_year_out.tsv"
TURNOVER_LONG_PATH = SELECTED_DIR / "country_epoch_architecture_turnover.tsv"
TURNOVER_SUMMARY_PATH = SELECTED_DIR / "country_epoch_architecture_turnover_summary.tsv"
ORIGIN_BRIDGE_PATH = SELECTED_DIR / "origin_burden_prevalence_shift.tsv"
ANTIGEN_COMBINED_PATH = SELECTED_DIR / "prn_specificity_negative_control.tsv"
MISSINGNESS_INTERPRETABILITY_PATH = SELECTED_DIR / "prn_interpretability_model.tsv"
QUALITY_SENSITIVITY_PATH = SELECTED_DIR / "selected_country_quality_restricted_sensitivity.tsv"
MISSINGNESS_TIPPING_GRID_PATH = SELECTED_DIR / "selected_country_missingness_tipping_grid.tsv"
READINESS_THRESHOLD_GRID_PATH = SELECTED_DIR / "selected_country_threshold_robustness_grid.tsv"
MISSINGNESS_DR_SUMMARY_PATH = SELECTED_DIR / "selected_country_dr_missingness_summary.tsv"
ASR_SCENARIO_REGISTRY_PATH = REPO_ROOT / "manuscript" / "figure_data" / "asr_scenario_registry.tsv"
ASR_ONE_GLOBAL_CLONE_SUMMARY_PATH = REPO_ROOT / "manuscript" / "figure_data" / "asr_one_global_clone_summary.tsv"
FIG03_ORIGINS_PATH = FIGURE_DATA_DIR / "fig03_independent_origins.tsv"
ARCHITECTURE_ORIGIN_VALIDATION_PATH = FIGURE_DATA_DIR / "architecture_origin_validation_matrix.tsv"
ORIGIN_CONFIDENCE_TIER_PATH = FIGURE_DATA_DIR / "origin_confidence_tier_table.tsv"


YEAR_SENS_SUPP_PATH = SUPP_DIR / "Supplementary_Table_35_Selected_Country_Leave_One_Year_Out.tsv"
TURNOVER_SUPP_PATH = SUPP_DIR / "Supplementary_Table_36_Country_Epoch_Architecture_Turnover.tsv"
ORIGIN_BRIDGE_SUPP_PATH = SUPP_DIR / "Supplementary_Table_37_Origin_Burden_Prevalence_Shift.tsv"
ANTIGEN_SUPP_PATH = SUPP_DIR / "Supplementary_Table_38_PRN_Specificity_Negative_Control.tsv"
MISSINGNESS_INTERPRETABILITY_SUPP_PATH = SUPP_DIR / "Supplementary_Table_39_prn_Interpretability_Model.tsv"
QUALITY_SENSITIVITY_SUPP_PATH = SUPP_DIR / "Supplementary_Table_40_Quality_Restricted_Selected_Country_Sensitivity.tsv"
MISSINGNESS_TIPPING_SUPP_PATH = SUPP_DIR / "Supplementary_Table_41_Missingness_Tipping_Point_Summary.tsv"
READINESS_THRESHOLD_SUPP_PATH = SUPP_DIR / "Supplementary_Table_42_Selected_Country_Threshold_Robustness.tsv"
MISSINGNESS_DR_SUPP_PATH = SUPP_DIR / "Supplementary_Table_43_Selected_Country_DR_Missingness_Summary.tsv"
ASR_SCENARIO_REGISTRY_SUPP_PATH = SUPP_DIR / "Supplementary_Table_44_ASR_Scenario_Registry.tsv"
ASR_SCENARIO_REGISTRY_READER_SUPP_PATH = SUPP_DIR / "Supplementary_Table_6_ASR_Scenario_Registry.tsv"
ORIGIN_CONFIDENCE_TIER_SUPP_PATH = SUPP_DIR / "Supplementary_Table_45_Origin_Confidence_Tiers.tsv"


PRIMARY_COUNTRY_PAIRS = {
    "USA": ("usa_wp_only", "usa_ap_prn_background"),
    "NZL": ("nzl_wp_only", "nzl_ap_with_prn"),
    "AUS": ("aus_wp_only", "aus_ap_with_prn"),
    "JPN": ("jpn_pre2012_mixed_ap", "jpn_ap_without_prn"),
}

QUALITY_SUBSET_RULES = [
    {
        "subset_id": "all_interpretable",
        "subset_label": "All interpretable genomes",
        "requires_reads": False,
        "max_n_contigs": np.nan,
        "min_n50": np.nan,
    },
    {
        "subset_id": "has_reads",
        "subset_label": "Read-linked interpretable genomes",
        "requires_reads": True,
        "max_n_contigs": np.nan,
        "min_n50": np.nan,
    },
    {
        "subset_id": "reads_low_fragmentation",
        "subset_label": "Read-linked interpretable genomes with <=100 contigs",
        "requires_reads": True,
        "max_n_contigs": 100,
        "min_n50": np.nan,
    },
    {
        "subset_id": "reads_high_contiguity",
        "subset_label": "Read-linked interpretable genomes with N50 >= 50 kb",
        "requires_reads": True,
        "max_n_contigs": np.nan,
        "min_n50": 50_000,
    },
    {
        "subset_id": "reads_qc_compact",
        "subset_label": "Read-linked interpretable genomes with <=100 contigs and N50 >= 50 kb",
        "requires_reads": True,
        "max_n_contigs": 100,
        "min_n50": 50_000,
    },
]

TIPPING_DELTAS = np.round(np.arange(-1.0, 1.001, 0.05), 2)
ODDS_MULTIPLIERS = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
PROBABILITY_EPSILON = 1e-3

CONTRAST_CORE_COLUMNS = [
    "country_iso3",
    "country_name",
    "previous_epoch_id",
    "next_epoch_id",
    "previous_epoch_label",
    "next_epoch_label",
    "contrast_type",
    "previous_prn_in_formulation",
    "next_prn_in_formulation",
    "previous_n_prn_interpretable",
    "next_n_prn_interpretable",
    "previous_naive_prevalence",
    "next_naive_prevalence",
    "previous_ipw_prevalence",
    "next_ipw_prevalence",
    "delta_naive_prevalence",
    "delta_ipw_prevalence",
    "previous_bound_lower",
    "previous_bound_upper",
    "next_bound_lower",
    "next_bound_upper",
    "bounds_direction",
    "contrast_eligible",
    "notes",
]

BOUND_WIDTH_POLICIES = {
    "wide": 1.0,
    "medium": 0.67,
    "strict": 0.5,
}

LOCUS_CONFIG = [
    {
        "locus": "prn",
        "label": "PRN",
        "category": "acellular_antigen",
        "signal_definition": "step4_coding_disrupted",
        "locus_length_bp": 2733,
    },
    {
        "locus": "fim2",
        "label": "Fim2",
        "category": "acellular_antigen",
        "signal_definition": "step2_marker_below_threshold",
        "locus_length_bp": 1266,
    },
    {
        "locus": "fim3",
        "label": "Fim3",
        "category": "acellular_antigen",
        "signal_definition": "step2_marker_below_threshold",
        "locus_length_bp": 1266,
    },
    {
        "locus": "brkA",
        "label": "BrkA",
        "category": "structure_matched_autotransporter",
        "signal_definition": "pseudo_control_marker_below_threshold",
        "locus_length_bp": 3033,
    },
    {
        "locus": "tcfA",
        "label": "TcfA",
        "category": "structure_matched_autotransporter",
        "signal_definition": "pseudo_control_marker_below_threshold",
        "locus_length_bp": 1944,
    },
    {
        "locus": "vag8",
        "label": "Vag8",
        "category": "structure_matched_autotransporter",
        "signal_definition": "pseudo_control_marker_below_threshold",
        "locus_length_bp": 2748,
    },
    {
        "locus": "sphB1",
        "label": "SphB1",
        "category": "structure_matched_autotransporter",
        "signal_definition": "pseudo_control_marker_below_threshold",
        "locus_length_bp": 2796,
    },
    {
        "locus": "phg",
        "label": "Phg",
        "category": "pertactin_homologous_autotransporter",
        "signal_definition": "pseudo_control_marker_below_threshold",
        "locus_length_bp": 1257,
    },
    {
        "locus": "bapC",
        "label": "BapC",
        "category": "secondary_autotransporter_reference_pseudogene_caveat",
        "signal_definition": "pseudo_control_marker_below_threshold",
        "locus_length_bp": 2987,
    },
    {
        "locus": "fhaB",
        "label": "FHA/FhaB",
        "category": "secondary_large_vaccine_adhesin",
        "signal_definition": "pseudo_control_marker_below_threshold",
        "locus_length_bp": 10773,
    },
]


def parse_bool(series: pd.Series) -> pd.Series:
    if series is None:
        return pd.Series(dtype=bool)
    normalized = (
        series.astype(str)
        .str.strip()
        .str.lower()
        .replace({"": np.nan, "nan": np.nan, "none": np.nan})
    )
    return normalized.isin({"true", "1", "yes", "y", "t"})


def normalize_text(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .where(~series.isna(), "")
        .str.strip()
        .replace({"nan": "", "None": "", "NA": ""})
    )


def event_label(event_id: str) -> str:
    event = str(event_id or "")
    if "gap1043" in event:
        return "IS481 gap1043"
    if "gap1042" in event and "is481" in event:
        return "IS481 gap1042"
    if "gap1041" in event and "is481" in event:
        return "IS481 gap1041"
    if "cov58" in event:
        return "Rearrangement cov58"
    if "cov91" in event:
        return "Rearrangement cov91"
    if "cov94" in event:
        return "Rearrangement cov94"
    if "other" in event and "gap1042" in event:
        return "Other gap1042"
    if "other" in event and "gap1041" in event:
        return "Other gap1041"
    if "other" in event and "gap1040" in event:
        return "Other gap1040"
    if "insufficient" in event:
        return "Insufficient"
    if not event:
        return "Unknown"
    return event


def safe_divide(num: float, den: float) -> float:
    if den is None or pd.isna(den) or den == 0:
        return np.nan
    return float(num) / float(den)


def continuity_corrected_prevalence(n_disrupted: float, n_interpretable: float) -> float:
    """Reference prevalence for odds-multiplier sensitivity when observed cells are 0 or all."""
    if pd.isna(n_disrupted) or pd.isna(n_interpretable) or n_interpretable < 0:
        return np.nan
    return float(n_disrupted + 0.5) / float(n_interpretable + 1.0)


def probability_from_odds_multiplier(reference_prevalence: float, odds_multiplier: float) -> float:
    if pd.isna(reference_prevalence) or pd.isna(odds_multiplier):
        return np.nan
    reference_prevalence = float(np.clip(reference_prevalence, PROBABILITY_EPSILON, 1.0 - PROBABILITY_EPSILON))
    odds = reference_prevalence / (1.0 - reference_prevalence)
    adjusted_odds = odds * float(odds_multiplier)
    return float(adjusted_odds / (1.0 + adjusted_odds))


def choose_sign(delta: float) -> str:
    if pd.isna(delta):
        return "not_estimable"
    if delta > 0:
        return "increase"
    if delta < 0:
        return "decrease"
    return "no_change"


def total_variation_distance(a: Iterable[float], b: Iterable[float]) -> float:
    a_arr = np.asarray(list(a), dtype=float)
    b_arr = np.asarray(list(b), dtype=float)
    if a_arr.size == 0 and b_arr.size == 0:
        return np.nan
    a_arr = np.nan_to_num(a_arr, nan=0.0)
    b_arr = np.nan_to_num(b_arr, nan=0.0)
    return 0.5 * np.abs(a_arr - b_arr).sum()


def first_nonmissing(*series: pd.Series) -> pd.Series:
    output = pd.Series(np.nan, index=series[0].index if series else None, dtype=float)
    for item in series:
        output = output.where(output.notna(), item)
    return output


def load_rescued_manifest_overrides() -> pd.DataFrame:
    if not RESCUED_MANIFEST_OVERRIDES_PATH.exists():
        return pd.DataFrame(columns=["assembly_accession", "prn_interpretable", "prn_disrupted"])
    frame = pd.read_csv(RESCUED_MANIFEST_OVERRIDES_PATH, sep="\t", dtype=str)
    frame["assembly_accession"] = normalize_text(frame["assembly_accession"])
    frame["prn_interpretable"] = parse_bool(frame["prn_interpretable"])
    frame["prn_disrupted"] = parse_bool(frame["prn_disrupted"])
    return frame.loc[frame["assembly_accession"] != "", ["assembly_accession", "prn_interpretable", "prn_disrupted"]]


def parse_missingness_model_summary(path: Path) -> dict[str, float]:
    metrics: dict[str, float] = {}
    if not path.exists():
        return metrics
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or ":" not in line:
            continue
        key, value = [part.strip() for part in line.split(":", 1)]
        if key == "Total samples with labels":
            metrics["labeled_samples_total"] = float(value)
        elif key == "Model samples (complete features)":
            metrics["model_samples_complete"] = float(value)
        elif key == "Out-of-fold Accuracy":
            metrics["model_accuracy_oof"] = float(value.split()[0])
        elif key == "Out-of-fold AUC":
            metrics["model_auc_oof"] = float(value.split()[0])
        elif key == "Accuracy":
            metrics["model_accuracy"] = float(value.split()[0])
        elif key == "AUC":
            metrics["model_auc"] = float(value.split()[0])
    if "model_accuracy_oof" in metrics:
        metrics["model_accuracy"] = metrics["model_accuracy_oof"]
    if "model_auc_oof" in metrics:
        metrics["model_auc"] = metrics["model_auc_oof"]
    return metrics


def clip_probability(values: pd.Series | np.ndarray, lower: float = PROBABILITY_EPSILON) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    upper = 1.0 - lower
    return np.clip(array, lower, upper)


def safe_quantile(values: Iterable[float], q: float) -> float:
    series = pd.Series(list(values), dtype=float).dropna()
    if series.empty:
        return np.nan
    return float(series.quantile(q))


def build_missingness_feature_frame(
    frame: pd.DataFrame,
    include_epoch_indicator: bool = False,
    epoch_reference: str | None = None,
) -> pd.DataFrame:
    output = pd.DataFrame(index=frame.index)
    output["year"] = pd.to_numeric(frame.get("year", pd.Series(index=frame.index, dtype=float)), errors="coerce")
    output["has_reads"] = parse_bool(frame.get("has_reads", pd.Series(index=frame.index, dtype=str))).astype(float)

    total_length = pd.to_numeric(frame.get("qc_total_length", pd.Series(index=frame.index, dtype=float)), errors="coerce")
    n_contigs = pd.to_numeric(frame.get("qc_n_contigs", pd.Series(index=frame.index, dtype=float)), errors="coerce")
    output["log_total_length"] = np.log1p(total_length.clip(lower=0))
    output["log_n_contigs"] = np.log1p(n_contigs.clip(lower=0))
    output["missing_total_length"] = total_length.isna().astype(float)
    output["missing_n_contigs"] = n_contigs.isna().astype(float)

    if include_epoch_indicator:
        epoch_id = normalize_text(frame.get("epoch_id", pd.Series(index=frame.index, dtype=str)))
        reference = epoch_reference or (epoch_id.iloc[0] if not epoch_id.empty else "")
        output["epoch_is_nonreference"] = epoch_id.ne(reference).astype(float)

    medians = output.median(axis=0, numeric_only=True)
    imputed = output.fillna(medians).fillna(0.0)
    return imputed.astype(float)


def fit_probability_model(
    features: pd.DataFrame,
    outcome: pd.Series,
    sample_weight: pd.Series | np.ndarray | None = None,
) -> dict[str, object]:
    y = pd.to_numeric(outcome, errors="coerce")
    mask = y.notna()
    if not mask.any():
        return {"model_type": "constant", "constant_probability": 0.5}

    x_use = features.loc[mask]
    y_use = y.loc[mask].astype(int)
    if y_use.nunique() < 2:
        return {"model_type": "constant", "constant_probability": float(y_use.mean())}

    weight_use = None
    if sample_weight is not None:
        weight_series = pd.Series(sample_weight, index=features.index, dtype=float)
        weight_use = weight_series.loc[mask].to_numpy(dtype=float)

    design = sm.add_constant(x_use.to_numpy(dtype=float), has_constant="add")
    try:
        model = sm.GLM(
            y_use.to_numpy(dtype=float),
            design,
            family=sm.families.Binomial(),
            freq_weights=weight_use,
        ).fit()
    except Exception:
        model = sm.GLM(
            y_use.to_numpy(dtype=float),
            design,
            family=sm.families.Binomial(),
            freq_weights=weight_use,
        ).fit_regularized(alpha=1e-6, L1_wt=0.0)
    return {"model_type": "logistic", "model": model}


def predict_probability(model_bundle: dict[str, object], features: pd.DataFrame) -> np.ndarray:
    if model_bundle["model_type"] == "constant":
        return np.repeat(float(model_bundle["constant_probability"]), len(features))
    model = model_bundle["model"]
    design = sm.add_constant(features.to_numpy(dtype=float), has_constant="add")
    return np.asarray(model.predict(design), dtype=float)


def make_stratified_folds(y: pd.Series, n_splits: int, random_state: int = 20260415) -> list[tuple[np.ndarray, np.ndarray]]:
    if n_splits < 2:
        return []
    rng = np.random.default_rng(random_state)
    y_array = pd.Series(y, dtype=int).to_numpy()
    fold_ids = np.full(len(y_array), -1, dtype=int)
    for class_value in np.unique(y_array):
        class_idx = np.flatnonzero(y_array == class_value)
        rng.shuffle(class_idx)
        for offset, idx in enumerate(class_idx):
            fold_ids[idx] = offset % n_splits
    folds: list[tuple[np.ndarray, np.ndarray]] = []
    all_idx = np.arange(len(y_array))
    for fold_id in range(n_splits):
        test_idx = all_idx[fold_ids == fold_id]
        train_idx = all_idx[fold_ids != fold_id]
        if len(test_idx) == 0 or len(train_idx) == 0:
            continue
        folds.append((train_idx, test_idx))
    return folds


def compute_accuracy(observed: pd.Series, predicted: pd.Series) -> float:
    obs = pd.Series(observed, dtype=float)
    pred = pd.Series(predicted, index=obs.index, dtype=float)
    mask = obs.notna() & pred.notna()
    if not mask.any():
        return np.nan
    return float((obs.loc[mask].astype(int) == (pred.loc[mask] >= 0.5).astype(int)).mean())


def compute_brier_score(observed: pd.Series, predicted: pd.Series) -> float:
    obs = pd.Series(observed, dtype=float)
    pred = pd.Series(predicted, index=obs.index, dtype=float)
    mask = obs.notna() & pred.notna()
    if not mask.any():
        return np.nan
    residual = obs.loc[mask].astype(float) - pred.loc[mask].astype(float)
    return float(np.mean(np.square(residual)))


def compute_roc_auc(observed: pd.Series, predicted: pd.Series) -> float:
    obs = pd.Series(observed, dtype=float)
    pred = pd.Series(predicted, index=obs.index, dtype=float)
    mask = obs.notna() & pred.notna()
    obs = obs.loc[mask].astype(int)
    pred = pred.loc[mask].astype(float)
    if obs.empty or obs.nunique() < 2:
        return np.nan
    n_pos = int(obs.sum())
    n_neg = int((1 - obs).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = pred.rank(method="average")
    rank_sum_pos = float(ranks.loc[obs == 1].sum())
    auc = (rank_sum_pos - (n_pos * (n_pos + 1) / 2.0)) / (n_pos * n_neg)
    return float(auc)


def build_brier_decomposition(observed: pd.Series, predicted: np.ndarray, n_bins: int = 10) -> dict[str, float]:
    obs = pd.Series(observed, dtype=float)
    pred = pd.Series(predicted, index=obs.index, dtype=float)
    mask = obs.notna() & pred.notna()
    obs = obs.loc[mask]
    pred = pred.loc[mask]
    if obs.empty:
        return {
            "brier_reliability": np.nan,
            "brier_resolution": np.nan,
            "brier_uncertainty": np.nan,
        }

    binned = pd.DataFrame({"observed": obs, "predicted": pred})
    try:
        binned["calibration_bin"] = pd.qcut(
            binned["predicted"],
            q=min(n_bins, binned["predicted"].nunique()),
            labels=False,
            duplicates="drop",
        )
    except ValueError:
        binned["calibration_bin"] = 0
    grouped = (
        binned.groupby("calibration_bin", dropna=False)
        .agg(bin_n=("observed", "size"), observed_mean=("observed", "mean"), predicted_mean=("predicted", "mean"))
        .reset_index(drop=True)
    )
    overall = float(obs.mean())
    grouped["bin_fraction"] = grouped["bin_n"] / grouped["bin_n"].sum()
    reliability = float((grouped["bin_fraction"] * (grouped["predicted_mean"] - grouped["observed_mean"]) ** 2).sum())
    resolution = float((grouped["bin_fraction"] * (grouped["observed_mean"] - overall) ** 2).sum())
    uncertainty = float(overall * (1.0 - overall))
    return {
        "brier_reliability": reliability,
        "brier_resolution": resolution,
        "brier_uncertainty": uncertainty,
    }


def build_calibration_curve_rows(
    observed: pd.Series,
    predicted: np.ndarray,
    model_variant: str,
    country_iso3: str = "",
    country_name: str = "",
) -> list[dict[str, object]]:
    obs = pd.Series(observed, dtype=float)
    pred = pd.Series(predicted, index=obs.index, dtype=float)
    mask = obs.notna() & pred.notna()
    obs = obs.loc[mask]
    pred = pred.loc[mask]
    if obs.empty:
        return []

    calibration = pd.DataFrame({"observed": obs, "predicted": pred})
    try:
        calibration["calibration_bin"] = pd.qcut(
            calibration["predicted"],
            q=min(10, calibration["predicted"].nunique()),
            labels=False,
            duplicates="drop",
        )
    except ValueError:
        calibration["calibration_bin"] = 0

    rows: list[dict[str, object]] = []
    for item in (
        calibration.groupby("calibration_bin", dropna=False)
        .agg(bin_n=("observed", "size"), observed_mean=("observed", "mean"), predicted_mean=("predicted", "mean"))
        .reset_index()
        .itertuples(index=False)
    ):
        rows.append(
            {
                "row_type": "calibration_bin",
                "model_variant": model_variant,
                "country_iso3": country_iso3,
                "country_name": country_name,
                "feature": "",
                "feature_definition": "",
                "calibration_bin": int(item.calibration_bin) + 1,
                "bin_n": int(item.bin_n),
                "bin_predicted_mean": float(item.predicted_mean),
                "bin_observed_mean": float(item.observed_mean),
            }
        )
    return rows


def estimate_calibration_slope_intercept(observed: pd.Series, predicted: np.ndarray) -> tuple[float, float]:
    obs = pd.Series(observed, dtype=float)
    pred = pd.Series(predicted, index=obs.index, dtype=float)
    mask = obs.notna() & pred.notna()
    obs = obs.loc[mask]
    pred = clip_probability(pred.loc[mask].to_numpy(dtype=float))
    if obs.empty or pd.Series(obs).nunique() < 2:
        return (np.nan, np.nan)
    try:
        design = sm.add_constant(np.log(pred / (1.0 - pred)))
        model = sm.GLM(obs.to_numpy(dtype=float), design, family=sm.families.Binomial()).fit()
    except Exception:
        return (np.nan, np.nan)
    params = model.params
    return (float(params[0]), float(params[1]) if len(params) > 1 else np.nan)


def cross_validated_probabilities(features: pd.DataFrame, outcome: pd.Series) -> np.ndarray:
    y = pd.to_numeric(outcome, errors="coerce")
    mask = y.notna()
    x_use = features.loc[mask]
    y_use = y.loc[mask].astype(int)
    output = np.repeat(np.nan, len(features))
    if len(y_use) < 4 or y_use.nunique() < 2:
        return output
    min_class = int(y_use.value_counts().min())
    n_splits = min(5, min_class)
    if n_splits < 2:
        return output

    fold_pred = np.repeat(np.nan, len(y_use))
    try:
        for train_idx, test_idx in make_stratified_folds(y_use, n_splits=n_splits):
            bundle = fit_probability_model(x_use.iloc[train_idx], y_use.iloc[train_idx])
            fold_pred[test_idx] = predict_probability(bundle, x_use.iloc[test_idx])
    except Exception:
        return output
    output[np.flatnonzero(mask.to_numpy())] = fold_pred
    return output


def summarise_probability_predictions(observed: pd.Series, predicted: np.ndarray) -> dict[str, float]:
    obs = pd.Series(observed, dtype=float)
    pred = pd.Series(predicted, index=obs.index, dtype=float)
    mask = obs.notna() & pred.notna()
    obs = obs.loc[mask]
    pred = pred.loc[mask]
    if obs.empty:
        return {"accuracy": np.nan, "auc": np.nan, "brier": np.nan}

    accuracy = compute_accuracy(obs, pred)
    auc = compute_roc_auc(obs, pred)
    brier = compute_brier_score(obs, pred)
    return {"accuracy": accuracy, "auc": auc, "brier": brier}


def collect_epoch_profiles(epoch_eligibility: pd.DataFrame, history: pd.DataFrame) -> pd.DataFrame:
    merged = epoch_eligibility.copy()
    if "confidence_level" not in merged.columns:
        confidence = (
            history[["country_iso3", "epoch_id", "confidence_level"]]
            .drop_duplicates(subset=["country_iso3", "epoch_id"])
            .copy()
        )
        confidence["confidence_level"] = normalize_text(confidence["confidence_level"]).replace("", "unknown")
        merged = merged.merge(confidence, on=["country_iso3", "epoch_id"], how="left")
    merged["confidence_level"] = normalize_text(merged["confidence_level"]).replace("", "unknown")
    rows: list[dict[str, object]] = []
    for country_iso3, frame in merged.groupby("country_iso3", dropna=False):
        ordered = frame.sort_values(["start_year", "epoch_id"]).copy()
        rows.append(
            {
                "country_iso3": country_iso3,
                "epoch_interpretable_profile": "; ".join(
                    f"{row.epoch_id}={int(pd.to_numeric(pd.Series([row.n_prn_interpretable]), errors='coerce').fillna(0).iloc[0])}"
                    for row in ordered.itertuples(index=False)
                ),
                "epoch_bound_width_profile": "; ".join(
                    f"{row.epoch_id}={pd.to_numeric(pd.Series([row.bound_width]), errors='coerce').iloc[0]:.3f}"
                    if pd.notna(pd.to_numeric(pd.Series([row.bound_width]), errors="coerce").iloc[0])
                    else f"{row.epoch_id}=NA"
                    for row in ordered.itertuples(index=False)
                ),
                "epoch_confidence_profile": "; ".join(
                    f"{row.epoch_id}={row.confidence_level or 'unknown'}" for row in ordered.itertuples(index=False)
                ),
            }
        )
    return pd.DataFrame(rows)


def load_inputs() -> dict[str, pd.DataFrame]:
    history = pd.read_csv(HISTORY_PATH, sep="\t")
    epochs = pd.read_csv(EPOCH_PREVALENCE_PATH, sep="\t")
    epoch_eligibility = pd.read_csv(EPOCH_ELIGIBILITY_PATH, sep="\t")
    scorecard = pd.read_csv(SELECTION_SCORECARD_PATH, sep="\t")
    year_comp = pd.read_csv(
        YEAR_COMPOSITION_PATH if YEAR_COMPOSITION_PATH.exists() else YEAR_COMPOSITION_FALLBACK_PATH,
        sep="\t",
    )
    contrasts = pd.read_csv(EPOCH_CONTRAST_PATH, sep="\t")
    origin_packages = pd.read_csv(ORIGIN_PACKAGE_PATH, sep="\t")
    origin_shift = pd.read_csv(ORIGIN_SHIFT_PATH, sep="\t")
    detection_shift = pd.read_csv(DETECTION_SHIFT_PATH, sep="\t")
    structure_reuse = pd.read_csv(STRUCTURE_REUSE_PATH, sep="\t")
    evidence_grid = pd.read_csv(EVIDENCE_GRID_PATH, sep="\t")
    ipw = pd.read_csv(IPW_PATH, sep="\t")
    manifest = pd.read_csv(MANIFEST_PATH, sep="\t", low_memory=False)
    assembly_qc = pd.read_csv(ASSEMBLY_QC_PATH, sep="\t", low_memory=False)
    step2 = pd.read_csv(STEP2_MARKERS_PATH, sep="\t", low_memory=False)
    step2 = step2.rename(
        columns={
            "biosample_accession": "Assembly BioSample Accession",
            "current_accession": "Current Accession",
        }
    )
    rescued_overrides = load_rescued_manifest_overrides()

    for df, year_col in [
        (history, "start_year"),
        (history, "end_year"),
        (epochs, "start_year"),
        (epochs, "end_year"),
        (epoch_eligibility, "start_year"),
        (epoch_eligibility, "end_year"),
        (year_comp, "year"),
        (year_comp, "epoch_start_year"),
        (year_comp, "epoch_end_year"),
        (ipw, "year"),
        (manifest, "year"),
    ]:
        if year_col in df.columns:
            df[year_col] = pd.to_numeric(df[year_col], errors="coerce")

    manifest["assembly_accession"] = normalize_text(manifest.get("assembly_accession", pd.Series(dtype=str)))
    manifest["prn_interpretable"] = parse_bool(manifest.get("prn_interpretable", pd.Series(dtype=str)))
    manifest["prn_disrupted"] = parse_bool(manifest.get("prn_disrupted", pd.Series(dtype=str)))
    manifest["country_iso3"] = normalize_text(manifest.get("country_iso3", pd.Series(dtype=str)))
    manifest["has_reads"] = parse_bool(manifest.get("has_reads", pd.Series(dtype=str)))
    manifest["n_contigs"] = pd.to_numeric(manifest.get("n_contigs", pd.Series(dtype=float)), errors="coerce")
    manifest["contig_n50"] = pd.to_numeric(manifest.get("contig_n50", pd.Series(dtype=float)), errors="coerce")
    manifest["total_sequence_length"] = pd.to_numeric(
        manifest.get("total_sequence_length", pd.Series(dtype=float)), errors="coerce"
    )

    if not rescued_overrides.empty:
        manifest = manifest.merge(rescued_overrides, on="assembly_accession", how="left", suffixes=("", "_override"))
        for column in ["prn_interpretable", "prn_disrupted"]:
            override_column = f"{column}_override"
            manifest[column] = manifest[override_column].where(manifest[override_column].notna(), manifest[column])
            manifest = manifest.drop(columns=[override_column])

    assembly_qc["assembly_accession"] = normalize_text(assembly_qc.get("assembly_accession", pd.Series(dtype=str)))
    assembly_qc["n_contigs"] = pd.to_numeric(assembly_qc.get("n_contigs", pd.Series(dtype=float)), errors="coerce")
    assembly_qc["n50"] = pd.to_numeric(assembly_qc.get("n50", pd.Series(dtype=float)), errors="coerce")
    assembly_qc["total_length"] = pd.to_numeric(assembly_qc.get("total_length", pd.Series(dtype=float)), errors="coerce")
    manifest = manifest.merge(
        assembly_qc[["assembly_accession", "n_contigs", "n50", "total_length"]].rename(
            columns={"n_contigs": "qc_n_contigs", "n50": "qc_n50", "total_length": "qc_total_length"}
        ),
        on="assembly_accession",
        how="left",
    )
    manifest["qc_n_contigs"] = first_nonmissing(manifest["n_contigs"], manifest["qc_n_contigs"])
    manifest["qc_n50"] = first_nonmissing(manifest["contig_n50"], manifest["qc_n50"])
    manifest["qc_total_length"] = first_nonmissing(manifest["total_sequence_length"], manifest["qc_total_length"])
    manifest = manifest.loc[manifest["country_iso3"].ne("") & manifest["year"].notna()].copy()

    epochs["country_iso3"] = normalize_text(epochs["country_iso3"])
    history["country_iso3"] = normalize_text(history["country_iso3"])
    epoch_eligibility["country_iso3"] = normalize_text(epoch_eligibility["country_iso3"])
    epoch_eligibility["country_name"] = normalize_text(epoch_eligibility["country_name"])
    for column in [
        "n_retained_genomes",
        "n_prn_interpretable",
        "n_prn_disrupted",
        "n_prn_uninterpretable_or_uncertain",
        "bound_width",
        "n_package_level_hard_anchors",
    ]:
        if column in epoch_eligibility.columns:
            epoch_eligibility[column] = pd.to_numeric(epoch_eligibility[column], errors="coerce")
    scorecard["country_iso3"] = normalize_text(scorecard["country_iso3"])
    scorecard["country_name"] = normalize_text(scorecard["country_name"])
    for column in [
        "n_package_level_hard_anchors",
        "n_comparable_epochs",
        "n_informative_epochs",
        "max_epoch_bound_width",
    ]:
        if column in scorecard.columns:
            scorecard[column] = pd.to_numeric(scorecard[column], errors="coerce")
    year_comp["country_iso3"] = normalize_text(year_comp["country_iso3"])
    ipw["country_iso3"] = normalize_text(ipw["country_iso3"])
    origin_packages["origin_country_iso3"] = normalize_text(origin_packages.get("origin_country_iso3", pd.Series(dtype=str)))
    evidence_grid["country_iso3"] = normalize_text(evidence_grid["country_iso3"])
    structure_reuse["country_iso3"] = normalize_text(structure_reuse["country_iso3"])
    structure_reuse["prn_event_id"] = normalize_text(structure_reuse["prn_event_id"])

    return {
        "history": history,
        "epochs": epochs,
        "epoch_eligibility": epoch_eligibility,
        "scorecard": scorecard,
        "year_comp": year_comp,
        "contrasts": contrasts,
        "origin_packages": origin_packages,
        "origin_shift": origin_shift,
        "detection_shift": detection_shift,
        "structure_reuse": structure_reuse,
        "evidence_grid": evidence_grid,
        "ipw": ipw,
        "manifest": manifest,
        "step2": step2,
    }


def build_interpretability_model_table(manifest: pd.DataFrame) -> pd.DataFrame:
    metrics = parse_missingness_model_summary(MISSINGNESS_MODEL_SUMMARY_PATH)
    label = parse_bool(manifest.get("prn_interpretable", pd.Series(index=manifest.index, dtype=str))).astype(int)
    feature_frame = build_missingness_feature_frame(manifest)
    feature_definitions = {
        "year": "Collection year of the genome record",
        "year_numeric": "Collection year of the genome record",
        "has_reads": "Raw-read linkage available for the assembly",
        "has_reads_numeric": "Raw-read linkage available for the assembly",
        "log_total_length": "Log-transformed assembly total length",
        "log_n_contigs": "Log-transformed contig count as a fragmentation proxy",
        "missing_total_length": "Indicator that assembly total length was unavailable",
        "missing_n_contigs": "Indicator that contig count was unavailable",
    }
    variant_features = {
        "year_only": ["year"],
        "reads_only": ["has_reads"],
        "assembly_only": ["log_total_length", "log_n_contigs", "missing_total_length", "missing_n_contigs"],
        "full_minus_reads": ["year", "log_total_length", "log_n_contigs", "missing_total_length", "missing_n_contigs"],
        "full_model": [
            "year",
            "has_reads",
            "log_total_length",
            "log_n_contigs",
            "missing_total_length",
            "missing_n_contigs",
        ],
    }
    variant_rows: list[dict[str, object]] = []
    calibration_rows: list[dict[str, object]] = []
    full_model_probabilities: np.ndarray | None = None

    for variant_name, columns in variant_features.items():
        features = feature_frame[columns].copy()
        bundle = fit_probability_model(features, label)
        probabilities = clip_probability(predict_probability(bundle, features))
        summary = summarise_probability_predictions(label, probabilities)
        cv_probabilities = cross_validated_probabilities(features, label)
        cv_summary = summarise_probability_predictions(label, cv_probabilities)
        calibration_intercept, calibration_slope = estimate_calibration_slope_intercept(label, probabilities)
        brier_parts = build_brier_decomposition(label, probabilities)

        variant_rows.append(
            {
                "row_type": "model_diagnostic",
                "model_variant": variant_name,
                "country_iso3": "",
                "country_name": "",
                "feature": "",
                "feature_definition": "",
                "calibration_bin": np.nan,
                "bin_n": np.nan,
                "bin_predicted_mean": np.nan,
                "bin_observed_mean": np.nan,
                "labeled_samples_total": float(len(label)),
                "model_samples_complete": float(len(label)),
                "model_accuracy": summary["accuracy"],
                "model_auc": summary["auc"],
                "model_brier": summary["brier"],
                "cv_model_accuracy": cv_summary["accuracy"],
                "cv_model_auc": cv_summary["auc"],
                "cv_model_brier": cv_summary["brier"],
                "calibration_intercept": calibration_intercept,
                "calibration_slope": calibration_slope,
                "brier_reliability": brier_parts["brier_reliability"],
                "brier_resolution": brier_parts["brier_resolution"],
                "brier_uncertainty": brier_parts["brier_uncertainty"],
                "standardized_coefficient": np.nan,
                "direction_of_association_with_interpretable_call": "",
                "mean": np.nan,
                "std": np.nan,
                "notes": (
                    "In-sample metrics on the labeled manifest after median imputation and missingness indicators."
                ),
            }
        )
        calibration_rows.extend(build_calibration_curve_rows(label, probabilities, variant_name))
        if variant_name == "full_model":
            full_model_probabilities = probabilities

    coeff_rows: list[dict[str, object]] = []
    if MISSINGNESS_MODEL_COEFFICIENTS_PATH.exists():
        coeff = pd.read_csv(MISSINGNESS_MODEL_COEFFICIENTS_PATH, sep="\t")
        coeff["feature"] = normalize_text(coeff["feature"])
        coeff["feature"] = coeff["feature"].replace({"year_numeric": "year", "has_reads_numeric": "has_reads"})
        coeff = coeff.loc[coeff["feature"].ne("intercept")].copy()
        coefficient_column = "coefficient_std" if "coefficient_std" in coeff.columns else "coefficient"
        mean_column = "mean" if "mean" in coeff.columns else "feature_mean"
        std_column = "std" if "std" in coeff.columns else "feature_scale"
        if "model" not in coeff.columns:
            coeff["model"] = "full_model"
        coeff["standardized_coefficient"] = pd.to_numeric(coeff[coefficient_column], errors="coerce")
        coeff["direction_of_association_with_interpretable_call"] = coeff["standardized_coefficient"].map(
            lambda value: "Positive" if pd.notna(value) and value >= 0 else "Negative"
        )
        coeff["feature_definition"] = coeff["feature"].map(feature_definitions).fillna("")
        coeff["mean_numeric"] = pd.to_numeric(coeff.get(mean_column, pd.Series(index=coeff.index)), errors="coerce")
        coeff["std_numeric"] = pd.to_numeric(coeff.get(std_column, pd.Series(index=coeff.index)), errors="coerce")
        coeff["coefficient_notes"] = np.where(
            coeff.get("oof_metric_provenance", pd.Series("", index=coeff.index)).fillna("").astype(str).str.len().gt(0),
            "Workflow missingness-model coefficient; training metrics are in-sample and OOF metrics are reported separately.",
            "Workflow-frozen standardized coefficient from the legacy missingness model export.",
        )
        for row in coeff.to_dict("records"):
            coeff_rows.append(
                {
                    "row_type": "feature_coefficient",
                    "model_variant": row.get("model", "full_model"),
                    "country_iso3": "",
                    "country_name": "",
                    "feature": row.get("feature", ""),
                    "feature_definition": row.get("feature_definition", ""),
                    "calibration_bin": np.nan,
                    "bin_n": np.nan,
                    "bin_predicted_mean": np.nan,
                    "bin_observed_mean": np.nan,
                    "labeled_samples_total": metrics.get("labeled_samples_total", float(len(label))),
                    "model_samples_complete": metrics.get("model_samples_complete", float(len(label))),
                    "model_accuracy": metrics.get("model_accuracy", np.nan),
                    "model_auc": metrics.get("model_auc", np.nan),
                    "model_brier": np.nan,
                    "cv_model_accuracy": np.nan,
                    "cv_model_auc": np.nan,
                    "cv_model_brier": np.nan,
                    "calibration_intercept": np.nan,
                    "calibration_slope": np.nan,
                    "brier_reliability": np.nan,
                    "brier_resolution": np.nan,
                    "brier_uncertainty": np.nan,
                    "standardized_coefficient": row.get("standardized_coefficient", np.nan),
                    "direction_of_association_with_interpretable_call": row.get(
                        "direction_of_association_with_interpretable_call",
                        "",
                    ),
                    "mean": row.get("mean_numeric", np.nan),
                    "std": row.get("std_numeric", np.nan),
                    "notes": row.get("coefficient_notes", ""),
                }
            )

    prevalence_rows: list[dict[str, object]] = []
    if full_model_probabilities is not None:
        prevalence_rows.append(
            {
                "row_type": "overall_summary",
                "model_variant": "full_model",
                "country_iso3": "",
                "country_name": "",
                "feature": "observed_interpretable_fraction",
                "feature_definition": "Observed fraction of genomes with interpretable PRN calls in the labeled manifest.",
                "calibration_bin": np.nan,
                "bin_n": np.nan,
                "bin_predicted_mean": np.nan,
                "bin_observed_mean": np.nan,
                "labeled_samples_total": float(len(label)),
                "model_samples_complete": float(len(label)),
                "model_accuracy": np.nan,
                "model_auc": np.nan,
                "model_brier": np.nan,
                "cv_model_accuracy": np.nan,
                "cv_model_auc": np.nan,
                "cv_model_brier": np.nan,
                "calibration_intercept": np.nan,
                "calibration_slope": np.nan,
                "brier_reliability": np.nan,
                "brier_resolution": np.nan,
                "brier_uncertainty": np.nan,
                "standardized_coefficient": float(label.mean()),
                "direction_of_association_with_interpretable_call": "",
                "mean": float(label.mean()),
                "std": float(label.std(ddof=0)),
                "notes": "Included so the supplementary table carries the manifest-wide interpretability base rate.",
            }
        )

    output = pd.DataFrame(coeff_rows + variant_rows + calibration_rows + prevalence_rows)
    feature_order = {
        "year": 0,
        "has_reads": 1,
        "log_total_length": 2,
        "log_n_contigs": 3,
        "missing_total_length": 4,
        "missing_n_contigs": 5,
        "observed_interpretable_fraction": 6,
        "": 99,
    }
    row_order = {
        "feature_coefficient": 0,
        "model_diagnostic": 1,
        "calibration_bin": 2,
        "overall_summary": 3,
    }
    output["row_order"] = output["row_type"].map(row_order).fillna(99)
    output["feature_order"] = output["feature"].map(feature_order).fillna(99)
    return output.sort_values(["row_order", "model_variant", "feature_order", "calibration_bin"]).drop(
        columns=["row_order", "feature_order"]
    ).reset_index(drop=True)


def build_quality_restricted_sensitivity(
    manifest: pd.DataFrame,
    history: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    def filter_epoch(country_iso3: str, epoch_id: str, subset_rule: dict[str, object]) -> pd.DataFrame:
        epoch_rows = history.loc[
            (history["country_iso3"] == country_iso3) & (history["epoch_id"] == epoch_id),
            ["start_year", "end_year"],
        ]
        if epoch_rows.empty:
            return pd.DataFrame(columns=manifest.columns)
        bounds = epoch_rows.iloc[0]
        frame = manifest.loc[
            (manifest["country_iso3"] == country_iso3)
            & manifest["year"].between(bounds["start_year"], bounds["end_year"], inclusive="both")
            & manifest["prn_interpretable"]
        ].copy()
        if subset_rule["requires_reads"]:
            frame = frame.loc[frame["has_reads"]]
        max_n_contigs = subset_rule["max_n_contigs"]
        if pd.notna(max_n_contigs):
            frame = frame.loc[frame["qc_n_contigs"].notna() & (frame["qc_n_contigs"] <= max_n_contigs)]
        min_n50 = subset_rule["min_n50"]
        if pd.notna(min_n50):
            frame = frame.loc[frame["qc_n50"].notna() & (frame["qc_n50"] >= min_n50)]
        return frame

    for country_iso3, (pre_epoch, post_epoch) in PRIMARY_COUNTRY_PAIRS.items():
        country_name_rows = history.loc[history["country_iso3"] == country_iso3, "country_name"].dropna()
        country_name = country_name_rows.iloc[0] if not country_name_rows.empty else country_iso3
        baseline_pre = filter_epoch(country_iso3, pre_epoch, QUALITY_SUBSET_RULES[0])
        baseline_post = filter_epoch(country_iso3, post_epoch, QUALITY_SUBSET_RULES[0])
        baseline_delta = safe_divide(baseline_post["prn_disrupted"].sum(), len(baseline_post)) - safe_divide(
            baseline_pre["prn_disrupted"].sum(), len(baseline_pre)
        )
        baseline_sign = choose_sign(baseline_delta)

        for subset_rule in QUALITY_SUBSET_RULES:
            pre = filter_epoch(country_iso3, pre_epoch, subset_rule)
            post = filter_epoch(country_iso3, post_epoch, subset_rule)
            pre_n = int(len(pre))
            post_n = int(len(post))
            pre_disrupted = int(pre["prn_disrupted"].sum()) if pre_n else 0
            post_disrupted = int(post["prn_disrupted"].sum()) if post_n else 0
            pre_prev = safe_divide(pre_disrupted, pre_n)
            post_prev = safe_divide(post_disrupted, post_n)
            delta = np.nan if pd.isna(pre_prev) or pd.isna(post_prev) else post_prev - pre_prev
            sign = choose_sign(delta)
            subset_estimable = pre_n > 0 and post_n > 0
            rows.append(
                {
                    "country_iso3": country_iso3,
                    "country_name": country_name,
                    "comparison_pre_epoch_id": pre_epoch,
                    "comparison_post_epoch_id": post_epoch,
                    "subset_id": subset_rule["subset_id"],
                    "subset_label": subset_rule["subset_label"],
                    "requires_reads": subset_rule["requires_reads"],
                    "max_n_contigs": subset_rule["max_n_contigs"],
                    "min_n50": subset_rule["min_n50"],
                    "pre_n_interpretable_subset": pre_n,
                    "pre_n_disrupted_subset": pre_disrupted,
                    "post_n_interpretable_subset": post_n,
                    "post_n_disrupted_subset": post_disrupted,
                    "pre_naive_prevalence_subset": pre_prev,
                    "post_naive_prevalence_subset": post_prev,
                    "subset_naive_delta": delta,
                    "subset_naive_sign": sign,
                    "subset_estimable": subset_estimable,
                    "passes_main_epoch_floor_8": pre_n >= 8 and post_n >= 8,
                    "baseline_naive_delta_all_interpretable": baseline_delta,
                    "baseline_naive_sign_all_interpretable": baseline_sign,
                    "direction_matches_all_interpretable": subset_estimable and sign == baseline_sign,
                    "subset_note": (
                        "one_or_both_epochs_empty_after_quality_filter"
                        if not subset_estimable
                        else ""
                    ),
                }
            )

    return pd.DataFrame(rows).sort_values(
        ["country_iso3", "subset_id"]
    ).reset_index(drop=True)


def build_missingness_tipping_summary(
    epoch_eligibility: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    grid_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for country_iso3, (pre_epoch, post_epoch) in PRIMARY_COUNTRY_PAIRS.items():
        frame = epoch_eligibility.loc[
            (epoch_eligibility["country_iso3"] == country_iso3)
            & (epoch_eligibility["epoch_id"].isin([pre_epoch, post_epoch]))
        ].copy()
        if len(frame) != 2:
            continue
        pre_row = frame.loc[frame["epoch_id"] == pre_epoch].iloc[0]
        post_row = frame.loc[frame["epoch_id"] == post_epoch].iloc[0]
        country_name = pre_row.get("country_name") or post_row.get("country_name") or country_iso3

        pre_i = float(pre_row["n_prn_interpretable"])
        pre_d = float(pre_row["n_prn_disrupted"])
        pre_m = float(pre_row["n_prn_uninterpretable_or_uncertain"])
        post_i = float(post_row["n_prn_interpretable"])
        post_d = float(post_row["n_prn_disrupted"])
        post_m = float(post_row["n_prn_uninterpretable_or_uncertain"])

        pre_obs = safe_divide(pre_d, pre_i)
        post_obs = safe_divide(post_d, post_i)
        pre_reference = continuity_corrected_prevalence(pre_d, pre_i)
        post_reference = continuity_corrected_prevalence(post_d, post_i)
        baseline_delta = np.nan if pd.isna(pre_obs) or pd.isna(post_obs) else post_obs - pre_obs
        baseline_sign = choose_sign(baseline_delta)

        for delta_shift in TIPPING_DELTAS:
            pre_missing_prob = np.clip(pre_obs + delta_shift, 0.0, 1.0) if pd.notna(pre_obs) else np.nan
            post_missing_prob = np.clip(post_obs + delta_shift, 0.0, 1.0) if pd.notna(post_obs) else np.nan
            pre_aug = safe_divide(pre_d + pre_m * pre_missing_prob, pre_i + pre_m)
            post_aug = safe_divide(post_d + post_m * post_missing_prob, post_i + post_m)
            adjusted_delta = np.nan if pd.isna(pre_aug) or pd.isna(post_aug) else post_aug - pre_aug
            adjusted_sign = choose_sign(adjusted_delta)
            sign_changed = adjusted_sign != baseline_sign
            full_reversal = (
                baseline_sign in {"increase", "decrease"}
                and adjusted_sign in {"increase", "decrease"}
                and adjusted_sign != baseline_sign
            )
            grid_rows.append(
                {
                    "country_iso3": country_iso3,
                    "country_name": country_name,
                    "comparison_pre_epoch_id": pre_epoch,
                    "comparison_post_epoch_id": post_epoch,
                    "sensitivity_model": "additive_probability_offset",
                    "delta_shift_applied_to_missing_probability": delta_shift,
                    "pre_missing_odds_multiplier": np.nan,
                    "post_missing_odds_multiplier": np.nan,
                    "odds_multiplier_log2_distance": np.nan,
                    "pre_observed_interpretable_prevalence": pre_obs,
                    "post_observed_interpretable_prevalence": post_obs,
                    "pre_continuity_corrected_reference_prevalence": pre_reference,
                    "post_continuity_corrected_reference_prevalence": post_reference,
                    "pre_missing_probability_assumed": pre_missing_prob,
                    "post_missing_probability_assumed": post_missing_prob,
                    "pre_adjusted_total_prevalence": pre_aug,
                    "post_adjusted_total_prevalence": post_aug,
                    "baseline_delta_interpretable_only": baseline_delta,
                    "adjusted_delta_total": adjusted_delta,
                    "baseline_sign": baseline_sign,
                    "adjusted_sign": adjusted_sign,
                    "sign_changed": sign_changed,
                    "full_reversal": full_reversal,
                }
            )

        for pre_multiplier in ODDS_MULTIPLIERS:
            for post_multiplier in ODDS_MULTIPLIERS:
                pre_missing_prob = probability_from_odds_multiplier(pre_reference, pre_multiplier)
                post_missing_prob = probability_from_odds_multiplier(post_reference, post_multiplier)
                pre_aug = safe_divide(pre_d + pre_m * pre_missing_prob, pre_i + pre_m)
                post_aug = safe_divide(post_d + post_m * post_missing_prob, post_i + post_m)
                adjusted_delta = np.nan if pd.isna(pre_aug) or pd.isna(post_aug) else post_aug - pre_aug
                adjusted_sign = choose_sign(adjusted_delta)
                sign_changed = adjusted_sign != baseline_sign
                full_reversal = (
                    baseline_sign in {"increase", "decrease"}
                    and adjusted_sign in {"increase", "decrease"}
                    and adjusted_sign != baseline_sign
                )
                log2_distance = max(abs(np.log2(pre_multiplier)), abs(np.log2(post_multiplier)))
                grid_rows.append(
                    {
                        "country_iso3": country_iso3,
                        "country_name": country_name,
                        "comparison_pre_epoch_id": pre_epoch,
                        "comparison_post_epoch_id": post_epoch,
                        "sensitivity_model": "odds_multiplier_pattern_mixture",
                        "delta_shift_applied_to_missing_probability": np.nan,
                        "pre_missing_odds_multiplier": pre_multiplier,
                        "post_missing_odds_multiplier": post_multiplier,
                        "odds_multiplier_log2_distance": log2_distance,
                        "pre_observed_interpretable_prevalence": pre_obs,
                        "post_observed_interpretable_prevalence": post_obs,
                        "pre_continuity_corrected_reference_prevalence": pre_reference,
                        "post_continuity_corrected_reference_prevalence": post_reference,
                        "pre_missing_probability_assumed": pre_missing_prob,
                        "post_missing_probability_assumed": post_missing_prob,
                        "pre_adjusted_total_prevalence": pre_aug,
                        "post_adjusted_total_prevalence": post_aug,
                        "baseline_delta_interpretable_only": baseline_delta,
                        "adjusted_delta_total": adjusted_delta,
                        "baseline_sign": baseline_sign,
                        "adjusted_sign": adjusted_sign,
                        "sign_changed": sign_changed,
                        "full_reversal": full_reversal,
                    }
                )

        country_grid = pd.DataFrame([row for row in grid_rows if row["country_iso3"] == country_iso3])
        additive_grid = country_grid.loc[country_grid["sensitivity_model"].eq("additive_probability_offset")]
        odds_grid = country_grid.loc[country_grid["sensitivity_model"].eq("odds_multiplier_pattern_mixture")]
        sign_change = additive_grid.loc[additive_grid["sign_changed"]]
        full_reversal = additive_grid.loc[additive_grid["full_reversal"]]
        odds_sign_change = odds_grid.loc[odds_grid["sign_changed"]].sort_values(
            ["odds_multiplier_log2_distance", "pre_missing_odds_multiplier", "post_missing_odds_multiplier"]
        )
        odds_full_reversal = odds_grid.loc[odds_grid["full_reversal"]].sort_values(
            ["odds_multiplier_log2_distance", "pre_missing_odds_multiplier", "post_missing_odds_multiplier"]
        )
        odds_sign_row = odds_sign_change.iloc[0] if not odds_sign_change.empty else None
        odds_reversal_row = odds_full_reversal.iloc[0] if not odds_full_reversal.empty else None
        summary_rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": country_name,
                "comparison_pre_epoch_id": pre_epoch,
                "comparison_post_epoch_id": post_epoch,
                "pre_n_interpretable": pre_i,
                "pre_n_missing": pre_m,
                "post_n_interpretable": post_i,
                "post_n_missing": post_m,
                "pre_missing_fraction": safe_divide(pre_m, pre_i + pre_m),
                "post_missing_fraction": safe_divide(post_m, post_i + post_m),
                "baseline_delta_interpretable_only": baseline_delta,
                "baseline_sign": baseline_sign,
                "sign_change_observed_within_grid": not sign_change.empty,
                "full_reversal_observed_within_grid": not full_reversal.empty,
                "min_delta_for_sign_change": (
                    np.nan if sign_change.empty else sign_change.iloc[0]["delta_shift_applied_to_missing_probability"]
                ),
                "min_abs_delta_for_sign_change": (
                    np.nan
                    if sign_change.empty
                    else sign_change["delta_shift_applied_to_missing_probability"].abs().min()
                ),
                "min_delta_for_full_reversal": (
                    np.nan if full_reversal.empty else full_reversal.iloc[0]["delta_shift_applied_to_missing_probability"]
                ),
                "min_abs_delta_for_full_reversal": (
                    np.nan
                    if full_reversal.empty
                    else full_reversal["delta_shift_applied_to_missing_probability"].abs().min()
                ),
                "odds_multiplier_grid_values": ",".join(str(value).rstrip("0").rstrip(".") for value in ODDS_MULTIPLIERS),
                "odds_multiplier_sign_change_observed_within_grid": not odds_sign_change.empty,
                "odds_multiplier_full_reversal_observed_within_grid": not odds_full_reversal.empty,
                "min_log2_odds_multiplier_distance_for_sign_change": (
                    np.nan if odds_sign_row is None else odds_sign_row["odds_multiplier_log2_distance"]
                ),
                "pre_odds_multiplier_at_min_sign_change": (
                    np.nan if odds_sign_row is None else odds_sign_row["pre_missing_odds_multiplier"]
                ),
                "post_odds_multiplier_at_min_sign_change": (
                    np.nan if odds_sign_row is None else odds_sign_row["post_missing_odds_multiplier"]
                ),
                "min_log2_odds_multiplier_distance_for_full_reversal": (
                    np.nan if odds_reversal_row is None else odds_reversal_row["odds_multiplier_log2_distance"]
                ),
                "pre_odds_multiplier_at_min_full_reversal": (
                    np.nan if odds_reversal_row is None else odds_reversal_row["pre_missing_odds_multiplier"]
                ),
                "post_odds_multiplier_at_min_full_reversal": (
                    np.nan if odds_reversal_row is None else odds_reversal_row["post_missing_odds_multiplier"]
                ),
            }
        )

    grid = pd.DataFrame(grid_rows).sort_values(
        [
            "country_iso3",
            "sensitivity_model",
            "delta_shift_applied_to_missing_probability",
            "pre_missing_odds_multiplier",
            "post_missing_odds_multiplier",
        ]
    ).reset_index(drop=True)
    summary = pd.DataFrame(summary_rows).sort_values("country_iso3").reset_index(drop=True)
    return summary, grid


def build_threshold_robustness(
    epoch_eligibility: pd.DataFrame,
    scorecard: pd.DataFrame,
    evidence_grid: pd.DataFrame,
    history: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    scorecard_use = scorecard[
        [
            "country_iso3",
            "country_name",
            "country_mechanism_anchor",
            "n_package_level_hard_anchors",
        ]
    ].drop_duplicates(subset=["country_iso3"]).copy()
    evidence_direction = (
        evidence_grid[["country_iso3", "prevalence_direction", "final_interpretation_tier"]]
        .drop_duplicates(subset=["country_iso3"])
        .set_index("country_iso3")
        .to_dict(orient="index")
    )
    epoch_use = epoch_eligibility.copy()
    if "confidence_level" not in epoch_use.columns:
        confidence = (
            history[["country_iso3", "epoch_id", "confidence_level"]]
            .drop_duplicates(subset=["country_iso3", "epoch_id"])
            .copy()
        )
        confidence["confidence_level"] = normalize_text(confidence["confidence_level"]).replace("", "unknown")
        epoch_use = epoch_use.merge(confidence, on=["country_iso3", "epoch_id"], how="left")
    epoch_use["confidence_level"] = normalize_text(epoch_use["confidence_level"]).replace("", "unknown")

    country_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    combo_id = 0
    for min_interpretable in [6, 8, 10, 12]:
        for min_epochs in [2, 3]:
            for bound_policy, bound_width_max in BOUND_WIDTH_POLICIES.items():
                for confidence_policy, allowed_confidence in [
                    ("high_or_medium", {"high", "medium"}),
                    ("high_only", {"high"}),
                ]:
                    for mechanism_anchor_required in [False, True]:
                        for package_level_anchor_required in [False, True]:
                            if package_level_anchor_required and not mechanism_anchor_required:
                                continue
                            combo_id += 1
                            stage1_countries: list[str] = []
                            triangulated_countries: list[str] = []
                            stage1_only_countries: list[str] = []
                            fails_stage1_countries: list[str] = []
                            focal_status: dict[str, str] = {}
                            for item in scorecard_use.itertuples(index=False):
                                country_epochs = epoch_use.loc[epoch_use["country_iso3"] == item.country_iso3].copy()
                                eligible_epochs = country_epochs.loc[
                                    (country_epochs["n_prn_interpretable"].fillna(0) >= min_interpretable)
                                    & (country_epochs["bound_width"].fillna(np.inf) <= bound_width_max)
                                    & country_epochs["confidence_level"].isin(allowed_confidence)
                                ].copy()
                                n_eligible_epochs = int(len(eligible_epochs))
                                stage1_eligible = n_eligible_epochs >= min_epochs
                                any_mechanism_anchor = normalize_text(pd.Series([item.country_mechanism_anchor])).iloc[0] != "none"
                                local_package_anchor = (
                                    pd.to_numeric(pd.Series([item.n_package_level_hard_anchors]), errors="coerce")
                                    .fillna(0)
                                    .iloc[0]
                                    > 0
                                )
                                stage1_countries.extend([item.country_iso3] if stage1_eligible else [])
                                if package_level_anchor_required:
                                    stage2_eligible = bool(stage1_eligible and local_package_anchor)
                                else:
                                    stage2_eligible = bool(stage1_eligible and any_mechanism_anchor)

                                if not stage1_eligible:
                                    selection_state = "fails_stage1"
                                    fails_stage1_countries.append(item.country_iso3)
                                elif stage2_eligible:
                                    selection_state = "triangulated"
                                    triangulated_countries.append(item.country_iso3)
                                else:
                                    selection_state = "stage1_only"
                                    stage1_only_countries.append(item.country_iso3)

                                direction = evidence_direction.get(item.country_iso3, {}).get("prevalence_direction", "")
                                focal_status[item.country_iso3] = (
                                    f"{selection_state}_{direction}"
                                    if selection_state != "fails_stage1" and direction
                                    else selection_state
                                )
                                country_rows.append(
                                    {
                                        "combo_id": combo_id,
                                        "min_interpretable_per_epoch": min_interpretable,
                                        "min_eligible_epochs": min_epochs,
                                        "bound_width_policy": bound_policy,
                                        "bound_width_max": bound_width_max,
                                        "confidence_policy": confidence_policy,
                                        "mechanism_anchor_required": mechanism_anchor_required,
                                        "package_level_anchor_required": package_level_anchor_required,
                                        "country_iso3": item.country_iso3,
                                        "country_name": item.country_name,
                                        "country_mechanism_anchor": item.country_mechanism_anchor,
                                        "n_package_level_hard_anchors": item.n_package_level_hard_anchors,
                                        "n_eligible_epochs": n_eligible_epochs,
                                        "stage1_epidemiologic_eligible": stage1_eligible,
                                        "stage2_mechanistic_eligible": stage2_eligible,
                                        "selection_state": selection_state,
                                        "included_in_stage1_under_combo": stage1_eligible,
                                        "included_in_stage2_under_combo": stage2_eligible,
                                        "prevalence_direction": direction,
                                        "final_interpretation_tier": evidence_direction.get(item.country_iso3, {}).get(
                                            "final_interpretation_tier", ""
                                        ),
                                    }
                                )
                            summary_rows.append(
                                {
                                    "combo_id": combo_id,
                                    "min_interpretable_per_epoch": min_interpretable,
                                    "min_eligible_epochs": min_epochs,
                                    "bound_width_policy": bound_policy,
                                    "bound_width_max": bound_width_max,
                                    "confidence_policy": confidence_policy,
                                    "mechanism_anchor_required": mechanism_anchor_required,
                                    "package_level_anchor_required": package_level_anchor_required,
                                    "n_stage1_countries": len(stage1_countries),
                                    "stage1_countries": ";".join(sorted(stage1_countries)),
                                    "n_triangulated_countries": len(triangulated_countries),
                                    "triangulated_countries": ";".join(sorted(triangulated_countries)),
                                    "n_stage1_only_countries": len(stage1_only_countries),
                                    "stage1_only_countries": ";".join(sorted(stage1_only_countries)),
                                    "n_fails_stage1_countries": len(fails_stage1_countries),
                                    "fails_stage1_countries": ";".join(sorted(fails_stage1_countries)),
                                    "usa_status": focal_status.get("USA", "fails_stage1"),
                                    "nzl_status": focal_status.get("NZL", "fails_stage1"),
                                    "jpn_status": focal_status.get("JPN", "fails_stage1"),
                                    "aus_status": focal_status.get("AUS", "fails_stage1"),
                                }
                            )

    country_grid = pd.DataFrame(country_rows).sort_values(
        ["combo_id", "country_iso3"]
    ).reset_index(drop=True)
    summary = pd.DataFrame(summary_rows).sort_values("combo_id").reset_index(drop=True)
    return summary, country_grid


def lookup_threshold_state(
    country_grid: pd.DataFrame,
    country_iso3: str,
    *,
    min_interpretable_per_epoch: int,
    min_eligible_epochs: int,
    bound_width_policy: str,
    confidence_policy: str,
    mechanism_anchor_required: bool,
    package_level_anchor_required: bool,
) -> str:
    row = country_grid.loc[
        (country_grid["country_iso3"] == country_iso3)
        & (country_grid["min_interpretable_per_epoch"] == min_interpretable_per_epoch)
        & (country_grid["min_eligible_epochs"] == min_eligible_epochs)
        & (country_grid["bound_width_policy"] == bound_width_policy)
        & (country_grid["confidence_policy"] == confidence_policy)
        & (country_grid["mechanism_anchor_required"] == mechanism_anchor_required)
        & (country_grid["package_level_anchor_required"] == package_level_anchor_required)
    ]
    if row.empty:
        return "missing"
    return str(row.iloc[0]["selection_state"])


def determine_first_failure_rule(country_iso3: str, country_grid: pd.DataFrame) -> str:
    baseline_state = lookup_threshold_state(
        country_grid,
        country_iso3,
        min_interpretable_per_epoch=6,
        min_eligible_epochs=2,
        bound_width_policy="wide",
        confidence_policy="high_or_medium",
        mechanism_anchor_required=False,
        package_level_anchor_required=False,
    )
    if baseline_state == "fails_stage1":
        return "fewer_than_two_liberal_stage1_epochs"
    if baseline_state == "stage1_only":
        return "no_country_level_mechanism_anchor"

    rule_checks = [
        (
            "high_confidence_only",
            {
                "min_interpretable_per_epoch": 6,
                "min_eligible_epochs": 2,
                "bound_width_policy": "wide",
                "confidence_policy": "high_only",
                "mechanism_anchor_required": False,
                "package_level_anchor_required": False,
            },
        ),
        (
            "min_interpretable_ge_8",
            {
                "min_interpretable_per_epoch": 8,
                "min_eligible_epochs": 2,
                "bound_width_policy": "wide",
                "confidence_policy": "high_or_medium",
                "mechanism_anchor_required": False,
                "package_level_anchor_required": False,
            },
        ),
        (
            "three_epoch_requirement",
            {
                "min_interpretable_per_epoch": 6,
                "min_eligible_epochs": 3,
                "bound_width_policy": "wide",
                "confidence_policy": "high_or_medium",
                "mechanism_anchor_required": False,
                "package_level_anchor_required": False,
            },
        ),
        (
            "medium_bound_width_requirement",
            {
                "min_interpretable_per_epoch": 6,
                "min_eligible_epochs": 2,
                "bound_width_policy": "medium",
                "confidence_policy": "high_or_medium",
                "mechanism_anchor_required": False,
                "package_level_anchor_required": False,
            },
        ),
        (
            "strict_bound_width_requirement",
            {
                "min_interpretable_per_epoch": 6,
                "min_eligible_epochs": 2,
                "bound_width_policy": "strict",
                "confidence_policy": "high_or_medium",
                "mechanism_anchor_required": False,
                "package_level_anchor_required": False,
            },
        ),
        (
            "min_interpretable_ge_10",
            {
                "min_interpretable_per_epoch": 10,
                "min_eligible_epochs": 2,
                "bound_width_policy": "wide",
                "confidence_policy": "high_or_medium",
                "mechanism_anchor_required": False,
                "package_level_anchor_required": False,
            },
        ),
        (
            "local_origin_package_anchor_requirement",
            {
                "min_interpretable_per_epoch": 6,
                "min_eligible_epochs": 2,
                "bound_width_policy": "wide",
                "confidence_policy": "high_or_medium",
                "mechanism_anchor_required": True,
                "package_level_anchor_required": True,
            },
        ),
        (
            "min_interpretable_ge_12",
            {
                "min_interpretable_per_epoch": 12,
                "min_eligible_epochs": 2,
                "bound_width_policy": "wide",
                "confidence_policy": "high_or_medium",
                "mechanism_anchor_required": False,
                "package_level_anchor_required": False,
            },
        ),
    ]
    for label, params in rule_checks:
        if lookup_threshold_state(country_grid, country_iso3, **params) != "triangulated":
            return label
    return "stable_across_tested_rules"


def build_selection_scorecard_multiverse_summary(
    scorecard: pd.DataFrame,
    epoch_eligibility: pd.DataFrame,
    history: pd.DataFrame,
    country_grid: pd.DataFrame,
) -> pd.DataFrame:
    scorecard_base = scorecard.drop(
        columns=[
            "selection_frequency_stage1",
            "selection_frequency_stage2",
            "selection_frequency_primary_only",
            "selection_frequency_context",
            "n_multiverse_combos",
            "n_stage1_combos",
            "n_stage2_combos",
            "n_primary_only_combos",
            "n_context_combos",
            "first_failure_rule",
            "epoch_interpretable_profile",
            "epoch_bound_width_profile",
            "epoch_confidence_profile",
            "screening_exclusion_primary_reason",
            "multiverse_screening_note",
        ],
        errors="ignore",
    )
    profile = collect_epoch_profiles(epoch_eligibility, history)
    summary = (
        country_grid.groupby(["country_iso3", "country_name"], dropna=False)
        .agg(
            selection_frequency_stage1=("stage1_epidemiologic_eligible", "mean"),
            selection_frequency_stage2=("stage2_mechanistic_eligible", "mean"),
            selection_frequency_primary_only=("selection_state", lambda s: float((s == "stage1_only").mean())),
            selection_frequency_context=("selection_state", lambda s: float((s == "fails_stage1").mean())),
            n_multiverse_combos=("combo_id", "nunique"),
            n_stage1_combos=("stage1_epidemiologic_eligible", lambda s: int(pd.Series(s).fillna(False).sum())),
            n_stage2_combos=("stage2_mechanistic_eligible", lambda s: int(pd.Series(s).fillna(False).sum())),
            n_primary_only_combos=("selection_state", lambda s: int((s == "stage1_only").sum())),
            n_context_combos=("selection_state", lambda s: int((s == "fails_stage1").sum())),
        )
        .reset_index()
    )
    summary["first_failure_rule"] = summary["country_iso3"].map(
        lambda iso3: determine_first_failure_rule(iso3, country_grid)
    )
    output = scorecard_base.merge(summary, on=["country_iso3", "country_name"], how="left").merge(
        profile,
        on="country_iso3",
        how="left",
    )
    output["screening_exclusion_primary_reason"] = output["first_failure_rule"].fillna("not_classified")
    output["multiverse_screening_note"] = output.apply(
        lambda row: (
            "Stable Stage 2 triangulated country across tested Stage 1 thresholds and anchor policies."
            if pd.to_numeric(pd.Series([row.get("selection_frequency_stage2")]), errors="coerce").fillna(0).iloc[0] >= 0.75
            else "Country is Stage 1 estimable but remains primary-only because triangulation is absent or sensitive to stricter anchor rules."
            if (
                pd.to_numeric(pd.Series([row.get("selection_frequency_stage1")]), errors="coerce").fillna(0).iloc[0] > 0
                and pd.to_numeric(pd.Series([row.get("selection_frequency_stage2")]), errors="coerce").fillna(0).iloc[0] == 0
            )
            else "Country remains Stage 1 estimable in part of the multiverse, but triangulation is not stable."
            if pd.to_numeric(pd.Series([row.get("selection_frequency_stage1")]), errors="coerce").fillna(0).iloc[0] > 0
            else "Country fails the liberal Stage 1 screen and is retained only as context."
        ),
        axis=1,
    )
    return output


def summarise_weight_distribution(weights: np.ndarray) -> dict[str, float]:
    series = pd.Series(weights, dtype=float).replace([np.inf, -np.inf], np.nan).dropna()
    if series.empty:
        return {
            "weight_mean": np.nan,
            "weight_max": np.nan,
            "weight_q95": np.nan,
            "weight_q99": np.nan,
            "weight_fraction_gt_10": np.nan,
            "weight_fraction_gt_20": np.nan,
        }
    return {
        "weight_mean": float(series.mean()),
        "weight_max": float(series.max()),
        "weight_q95": float(series.quantile(0.95)),
        "weight_q99": float(series.quantile(0.99)),
        "weight_fraction_gt_10": float((series > 10).mean()),
        "weight_fraction_gt_20": float((series > 20).mean()),
    }


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    value_arr = np.asarray(values, dtype=float)
    weight_arr = np.asarray(weights, dtype=float)
    mask = np.isfinite(value_arr) & np.isfinite(weight_arr)
    if mask.sum() == 0:
        return np.nan
    value_arr = value_arr[mask]
    weight_arr = weight_arr[mask]
    if weight_arr.sum() <= 0:
        return np.nan
    return float(np.sum(value_arr * weight_arr) / np.sum(weight_arr))


def probability_or_nan(value: float) -> float:
    if not np.isfinite(value):
        return np.nan
    if value < 0.0 or value > 1.0:
        return np.nan
    return float(value)


def estimator_status(values: Iterable[float]) -> str:
    finite_values = [float(value) for value in values if np.isfinite(float(value))]
    if not finite_values:
        return "not_estimable"
    if any(value < 0.0 or value > 1.0 for value in finite_values):
        return "failed_out_of_bounds"
    return "ok"


def classify_identifiability_tier(row: pd.Series) -> str:
    if not bool(row.get("sign_stable_across_estimators", False)):
        return "fragile"
    reversal = pd.to_numeric(pd.Series([row.get("min_abs_delta_for_full_reversal")]), errors="coerce").iloc[0]
    odds_reversal = pd.to_numeric(
        pd.Series([row.get("min_log2_odds_multiplier_distance_for_full_reversal")]), errors="coerce"
    ).iloc[0]
    bounds_direction = str(row.get("bounds_direction", ""))
    if (
        ("nonoverlapping" in bounds_direction)
        and (pd.isna(reversal) or reversal >= 0.5)
        and (pd.isna(odds_reversal) or odds_reversal >= 2.0)
    ):
        return "stable"
    if (pd.isna(reversal) or reversal >= 0.25) and (pd.isna(odds_reversal) or odds_reversal >= 1.0):
        return "bounded"
    return "fragile"


def build_missingness_dr_summary(
    manifest: pd.DataFrame,
    history: pd.DataFrame,
    contrasts: pd.DataFrame,
    tipping_summary: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    epoch_manifest = assign_epochs_to_manifest(manifest, history)
    if epoch_manifest.empty:
        return pd.DataFrame(), contrasts.copy()

    contrast_rows: list[dict[str, object]] = []
    for country_iso3, (pre_epoch, post_epoch) in PRIMARY_COUNTRY_PAIRS.items():
        country_frame = epoch_manifest.loc[
            (epoch_manifest["country_iso3"] == country_iso3)
            & (epoch_manifest["epoch_id"].isin([pre_epoch, post_epoch]))
        ].copy()
        if country_frame.empty:
            continue
        country_name = (
            country_frame["country_name"].dropna().iloc[0]
            if "country_name" in country_frame.columns and country_frame["country_name"].notna().any()
            else history.loc[history["country_iso3"] == country_iso3, "country_name"].dropna().iloc[0]
        )
        country_frame["epoch_id"] = normalize_text(country_frame["epoch_id"])
        country_frame["prn_interpretable"] = parse_bool(country_frame["prn_interpretable"])
        country_frame["prn_disrupted"] = parse_bool(country_frame["prn_disrupted"])

        features = build_missingness_feature_frame(
            country_frame,
            include_epoch_indicator=True,
            epoch_reference=pre_epoch,
        )
        observed = country_frame["prn_interpretable"].astype(int)
        propensity_bundle = fit_probability_model(features, observed)
        propensity_in_sample = clip_probability(predict_probability(propensity_bundle, features))
        propensity_cross_fit = cross_validated_probabilities(features, observed)
        propensity = clip_probability(
            np.where(
                np.isfinite(propensity_cross_fit),
                propensity_cross_fit,
                propensity_in_sample,
            )
        )
        inverse_probability = 1.0 / propensity

        observed_outcome = country_frame["prn_disrupted"].where(country_frame["prn_interpretable"]).astype(float)
        observed_feature_frame = features.loc[country_frame["prn_interpretable"]].copy()
        observed_outcome_frame = observed_outcome.loc[country_frame["prn_interpretable"]].copy()
        outcome_bundle = fit_probability_model(observed_feature_frame, observed_outcome_frame)
        outcome_pred = clip_probability(predict_probability(outcome_bundle, features))
        outcome_cross_fit = cross_validated_probabilities(observed_feature_frame, observed_outcome_frame)
        if np.isfinite(outcome_cross_fit).any():
            observed_idx = np.flatnonzero(country_frame["prn_interpretable"].to_numpy(dtype=bool))
            replace_mask = np.isfinite(outcome_cross_fit)
            outcome_pred[observed_idx[replace_mask]] = clip_probability(outcome_cross_fit[replace_mask])

        observed_mask = country_frame["prn_interpretable"].to_numpy(dtype=bool)
        y_filled = observed_outcome.fillna(0.0).to_numpy(dtype=float)
        weights_cap20 = np.minimum(inverse_probability, 20.0)
        observed_weights = inverse_probability[observed_mask]
        p01 = safe_quantile(observed_weights, 0.01)
        p99 = safe_quantile(observed_weights, 0.99)
        p05 = safe_quantile(observed_weights, 0.05)
        p95 = safe_quantile(observed_weights, 0.95)
        weights_p01_p99 = np.clip(inverse_probability, p01, p99) if pd.notna(p01) and pd.notna(p99) else inverse_probability
        weights_p05_p95 = np.clip(inverse_probability, p05, p95) if pd.notna(p05) and pd.notna(p95) else inverse_probability
        overlap_weights = 1.0 - propensity

        epoch_estimates: dict[str, dict[str, float]] = {}
        aipw_raw_estimates: dict[str, float] = {}
        for epoch_id in [pre_epoch, post_epoch]:
            epoch_mask = country_frame["epoch_id"].eq(epoch_id).to_numpy(dtype=bool)
            observed_epoch = epoch_mask & observed_mask
            raw_aipw = (
                float(
                    np.mean(
                        outcome_pred[epoch_mask]
                        + observed_mask[epoch_mask]
                        * (y_filled[epoch_mask] - outcome_pred[epoch_mask])
                        / propensity[epoch_mask]
                    )
                )
                if epoch_mask.sum() > 0
                else np.nan
            )
            aipw_raw_estimates[epoch_id] = raw_aipw
            epoch_estimates[epoch_id] = {
                "naive": safe_divide(y_filled[observed_epoch].sum(), observed_epoch.sum()),
                "ipw_cap20": weighted_mean(y_filled[observed_epoch], weights_cap20[observed_epoch]),
                "ipw_untruncated": weighted_mean(y_filled[observed_epoch], inverse_probability[observed_epoch]),
                "ipw_p01_p99": weighted_mean(y_filled[observed_epoch], weights_p01_p99[observed_epoch]),
                "ipw_p05_p95": weighted_mean(y_filled[observed_epoch], weights_p05_p95[observed_epoch]),
                "overlap_observed": weighted_mean(y_filled[observed_epoch], overlap_weights[observed_epoch]),
                "aipw": probability_or_nan(raw_aipw),
            }

        delta_map = {
            f"delta_{key}_prevalence": (
                np.nan
                if pd.isna(epoch_estimates[pre_epoch][key]) or pd.isna(epoch_estimates[post_epoch][key])
                else epoch_estimates[post_epoch][key] - epoch_estimates[pre_epoch][key]
            )
            for key in epoch_estimates[pre_epoch]
        }
        sign_set = {
            choose_sign(value)
            for value in delta_map.values()
            if pd.notna(value)
        } - {"not_estimable"}
        sign_set_no_zero = sign_set - {"no_change"}
        sign_stable = len(sign_set_no_zero) <= 1

        diagnostics = summarise_weight_distribution(observed_weights)
        interpretable_propensity = propensity[observed_mask]
        missing_propensity = propensity[~observed_mask]
        aipw_status = estimator_status(aipw_raw_estimates.values())
        contrast_rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": country_name,
                "previous_epoch_id": pre_epoch,
                "next_epoch_id": post_epoch,
                "n_total_previous_epoch": int(country_frame["epoch_id"].eq(pre_epoch).sum()),
                "n_total_next_epoch": int(country_frame["epoch_id"].eq(post_epoch).sum()),
                "n_observed_previous_epoch": int((country_frame["epoch_id"].eq(pre_epoch) & country_frame["prn_interpretable"]).sum()),
                "n_observed_next_epoch": int((country_frame["epoch_id"].eq(post_epoch) & country_frame["prn_interpretable"]).sum()),
                "propensity_median_interpretable": safe_quantile(interpretable_propensity, 0.5),
                "propensity_q05_interpretable": safe_quantile(interpretable_propensity, 0.05),
                "propensity_q95_interpretable": safe_quantile(interpretable_propensity, 0.95),
                "propensity_median_missing": safe_quantile(missing_propensity, 0.5),
                "propensity_q05_missing": safe_quantile(missing_propensity, 0.05),
                "propensity_q95_missing": safe_quantile(missing_propensity, 0.95),
                "propensity_crossfit_fraction": float(np.isfinite(propensity_cross_fit).mean()),
                "outcome_crossfit_fraction_observed": float(np.isfinite(outcome_cross_fit).mean()),
                "p01_weight_truncation": p01,
                "p99_weight_truncation": p99,
                "p05_weight_truncation": p05,
                "p95_weight_truncation": p95,
                "dr_estimator_label": "aipw_crossfit_bounded_fail_closed",
                "aipw_estimator_status": aipw_status,
                "aipw_out_of_bounds": aipw_status == "failed_out_of_bounds",
                "previous_aipw_raw_prevalence": aipw_raw_estimates[pre_epoch],
                "next_aipw_raw_prevalence": aipw_raw_estimates[post_epoch],
                "sign_stable_across_estimators": sign_stable,
                "consistent_nonzero_delta_sign": next(iter(sign_set_no_zero)) if len(sign_set_no_zero) == 1 else "mixed_or_zero",
                **diagnostics,
                **{
                    f"previous_{key}_prevalence": value
                    for key, value in epoch_estimates[pre_epoch].items()
                },
                **{
                    f"next_{key}_prevalence": value
                    for key, value in epoch_estimates[post_epoch].items()
                },
                **delta_map,
            }
        )

    summary = pd.DataFrame(contrast_rows).sort_values("country_iso3").reset_index(drop=True)
    core_columns = [column for column in CONTRAST_CORE_COLUMNS if column in contrasts.columns]
    augmented = contrasts.loc[:, core_columns].copy() if core_columns else contrasts.copy()
    merge_columns = [
        "country_iso3",
        "previous_epoch_id",
        "next_epoch_id",
    ]
    summary_for_merge = summary.drop(columns=["country_name"], errors="ignore")
    collision_renames = {
        column: f"{column}_dr"
        for column in summary_for_merge.columns
        if column not in merge_columns and column in augmented.columns
    }
    summary_for_merge = summary_for_merge.rename(columns=collision_renames)
    augmented = augmented.merge(summary_for_merge, on=merge_columns, how="left")
    if tipping_summary is not None and not tipping_summary.empty:
        tipping_merge_columns = [
            column
            for column in [
                "country_iso3",
                "comparison_pre_epoch_id",
                "comparison_post_epoch_id",
                "min_abs_delta_for_full_reversal",
                "min_abs_delta_for_sign_change",
                "odds_multiplier_grid_values",
                "odds_multiplier_sign_change_observed_within_grid",
                "odds_multiplier_full_reversal_observed_within_grid",
                "min_log2_odds_multiplier_distance_for_sign_change",
                "pre_odds_multiplier_at_min_sign_change",
                "post_odds_multiplier_at_min_sign_change",
                "min_log2_odds_multiplier_distance_for_full_reversal",
                "pre_odds_multiplier_at_min_full_reversal",
                "post_odds_multiplier_at_min_full_reversal",
            ]
            if column in tipping_summary.columns
        ]
        augmented = augmented.merge(
            tipping_summary[tipping_merge_columns].rename(
                columns={
                    "comparison_pre_epoch_id": "previous_epoch_id",
                    "comparison_post_epoch_id": "next_epoch_id",
                }
            ),
            on=["country_iso3", "previous_epoch_id", "next_epoch_id"],
            how="left",
        )
    augmented["identifiability_tier"] = augmented.apply(classify_identifiability_tier, axis=1)
    return summary, augmented


def summarize_asr_output_dir(output_dir: Path | str | float) -> dict[str, float]:
    if output_dir is None or (isinstance(output_dir, float) and pd.isna(output_dir)):
        return {
            "largest_disrupted_clade_share": np.nan,
            "max_disrupted_tips_per_origin": np.nan,
        }
    output_path = REPO_ROOT / str(output_dir)
    origin_path = output_path / "origin_events.tsv"
    tip_path = output_path / "tip_states.tsv"
    if not origin_path.exists() or not tip_path.exists():
        return {
            "largest_disrupted_clade_share": np.nan,
            "max_disrupted_tips_per_origin": np.nan,
        }
    origin = pd.read_csv(origin_path, sep="\t")
    tips = pd.read_csv(tip_path, sep="\t")
    disrupted_tip_count = int(normalize_text(tips.get("prn_state", pd.Series(dtype=str))).eq("disrupted").sum())
    max_disrupted = pd.to_numeric(origin.get("n_tips_disrupted", pd.Series(dtype=float)), errors="coerce").max()
    return {
        "largest_disrupted_clade_share": safe_divide(max_disrupted, disrupted_tip_count),
        "max_disrupted_tips_per_origin": max_disrupted,
    }


def count_true_like(values: Iterable[object]) -> int:
    series = pd.Series(list(values), dtype=object)
    return int(series.map(lambda value: bool(value) if pd.notna(value) else False).sum())


def build_asr_scenario_registry() -> tuple[pd.DataFrame, pd.DataFrame]:
    rooting = pd.read_csv(ASR_ROOTING_SENSITIVITY_PATH, sep="\t")
    rooting["scenario"] = normalize_text(rooting["scenario"])
    support = pd.read_csv(ASR_SENSITIVITY_PATH, sep="\t")
    support["scenario"] = normalize_text(support["scenario"])
    support_lookup = support.drop_duplicates(subset=["scenario"]).set_index("scenario").to_dict(orient="index")
    mk = pd.read_csv(ASR_MK_PATH, sep="\t")
    mk["scenario"] = normalize_text(mk["scenario"])
    mk_lookup = mk.drop_duplicates(subset=["scenario"]).set_index("scenario").to_dict(orient="index")

    registry_rows: list[dict[str, object]] = []
    for row in rooting.itertuples(index=False):
        scenario = row.scenario
        mk_row = mk_lookup.get(scenario, {})
        output_summary = summarize_asr_output_dir(getattr(row, "output_dir", np.nan))
        fitch_origin_events = pd.to_numeric(
            pd.Series([getattr(row, "fitch_origin_events", np.nan)]), errors="coerce"
        ).iloc[0]
        tip_count = pd.to_numeric(pd.Series([getattr(row, "tip_count", np.nan)]), errors="coerce").iloc[0]
        disrupted_tip_count = pd.to_numeric(
            pd.Series([getattr(row, "disrupted_tip_count", np.nan)]), errors="coerce"
        ).iloc[0]
        notes = getattr(row, "notes", "")
        support_scenario = (
            "composition_pruned_primary_quality_frame"
            if scenario == "composition_filtered_reference_rooted_primary"
            else scenario
        )
        support_row = support_lookup.get(support_scenario, {})
        if support_row:
            tip_count = pd.to_numeric(pd.Series([support_row.get("tip_count", np.nan)]), errors="coerce").iloc[0]
            disrupted_tip_count = pd.to_numeric(
                pd.Series([support_row.get("disrupted_tip_count", np.nan)]), errors="coerce"
            ).iloc[0]
            fitch_origin_events = pd.to_numeric(
                pd.Series([support_row.get("fitch_origin_events", np.nan)]), errors="coerce"
            ).iloc[0]
            pastml_origin_events = pd.to_numeric(
                pd.Series([support_row.get("pastml_origin_events", np.nan)]), errors="coerce"
            ).iloc[0]
            notes = support_row.get("notes", notes)
        else:
            pastml_origin_events = pd.to_numeric(
                pd.Series([getattr(row, "pastml_origin_events", np.nan)]), errors="coerce"
            ).iloc[0]
        if scenario == "composition_filtered_reference_rooted_primary" and PRIMARY_ORIGIN_EVENTS_PATH.exists():
            primary_origins = pd.read_csv(PRIMARY_ORIGIN_EVENTS_PATH, sep="\t")
            if not primary_origins.empty:
                primary_origin_count = int(len(primary_origins))
                if pd.isna(fitch_origin_events):
                    fitch_origin_events = primary_origin_count
                elif int(fitch_origin_events) != primary_origin_count:
                    fitch_origin_events = primary_origin_count
                tree_ids = normalize_text(primary_origins.get("phylo_tree_id", pd.Series(dtype=str)))
                if tree_ids.str.contains("chn_rescue", case=False, na=False).any():
                    rescue_note = (
                        "Primary ASR quality frame with targeted CHN legacy-gap PRN rescue; "
                        "33 nonreference IQ-TREE composition-failed tips pruned."
                    )
                    notes = f"{notes} {rescue_note}".strip()
        registry_rows.append(
            {
                "scenario_id": scenario,
                "scenario_source": "rooting_sensitivity",
                "scenario_class": (
                    "primary"
                    if scenario == "composition_filtered_reference_rooted_primary"
                    else "rooting"
                ),
                "analysis_frame": getattr(row, "analysis_frame", ""),
                "rooting_mode": getattr(row, "rooting_mode", ""),
                "tip_count": tip_count,
                "disrupted_tip_count": disrupted_tip_count,
                "fitch_origin_events": fitch_origin_events,
                "pastml_origin_events": pastml_origin_events,
                "mk_origin_count_mean": pd.to_numeric(
                    pd.Series([mk_row.get("mk_origin_count_mean", np.nan)]), errors="coerce"
                ).iloc[0],
                "mk_origin_count_lower_95": pd.to_numeric(
                    pd.Series([mk_row.get("mk_origin_count_lower_95", np.nan)]), errors="coerce"
                ).iloc[0],
                "mk_origin_count_upper_95": pd.to_numeric(
                    pd.Series([mk_row.get("mk_origin_count_upper_95", np.nan)]), errors="coerce"
                ).iloc[0],
                "largest_disrupted_clade_share": output_summary["largest_disrupted_clade_share"],
                "max_disrupted_tips_per_origin": output_summary["max_disrupted_tips_per_origin"],
                "rejects_one_global_clone_fitch": bool(
                    pd.to_numeric(pd.Series([fitch_origin_events]), errors="coerce").fillna(0).iloc[0] > 1
                ),
                "rejects_one_global_clone_mk95": bool(
                    pd.to_numeric(pd.Series([mk_row.get("mk_origin_count_lower_95", np.nan)]), errors="coerce").fillna(0).iloc[0] > 1
                )
                if pd.notna(pd.to_numeric(pd.Series([mk_row.get("mk_origin_count_lower_95", np.nan)]), errors="coerce").iloc[0])
                else np.nan,
                "notes": notes,
            }
        )

    support_output_dir = {
        "unpruned_support_ge_70": "outputs/workflow/asr_sensitivity/support_70",
        "unpruned_support_ge_90": "outputs/workflow/asr_sensitivity/support_90",
    }
    for row in support.loc[support["scenario"].str.contains("support_ge_", na=False)].itertuples(index=False):
        output_summary = summarize_asr_output_dir(support_output_dir.get(row.scenario, np.nan))
        registry_rows.append(
            {
                "scenario_id": row.scenario,
                "scenario_source": "support_threshold_sensitivity",
                "scenario_class": "support_threshold",
                "analysis_frame": "unpruned_comparability_frame",
                "rooting_mode": "reference_support_threshold",
                "tip_count": pd.to_numeric(pd.Series([row.tip_count]), errors="coerce").iloc[0],
                "disrupted_tip_count": pd.to_numeric(pd.Series([row.disrupted_tip_count]), errors="coerce").iloc[0],
                "fitch_origin_events": pd.to_numeric(pd.Series([row.fitch_origin_events]), errors="coerce").iloc[0],
                "pastml_origin_events": pd.to_numeric(pd.Series([row.pastml_origin_events]), errors="coerce").iloc[0],
                "mk_origin_count_mean": np.nan,
                "mk_origin_count_lower_95": np.nan,
                "mk_origin_count_upper_95": np.nan,
                "largest_disrupted_clade_share": output_summary["largest_disrupted_clade_share"],
                "max_disrupted_tips_per_origin": output_summary["max_disrupted_tips_per_origin"],
                "rejects_one_global_clone_fitch": bool(
                    pd.to_numeric(pd.Series([row.fitch_origin_events]), errors="coerce").fillna(0).iloc[0] > 1
                ),
                "rejects_one_global_clone_mk95": np.nan,
                "notes": row.notes,
            }
        )

    for scheme_dir in sorted(ASR_RESAMPLING_DIR.glob("*_balanced")):
        for replicate_dir in sorted(scheme_dir.glob("replicate_*")):
            origin_path = replicate_dir / "origin_events.tsv"
            tip_path = replicate_dir / "tip_states.tsv"
            if not origin_path.exists() or not tip_path.exists():
                continue
            origin = pd.read_csv(origin_path, sep="\t")
            tips = pd.read_csv(tip_path, sep="\t")
            disrupted_tip_count = int(normalize_text(tips.get("prn_state", pd.Series(dtype=str))).eq("disrupted").sum())
            fitch_origin_events = int(len(origin))
            max_disrupted = pd.to_numeric(origin.get("n_tips_disrupted", pd.Series(dtype=float)), errors="coerce").max()
            registry_rows.append(
                {
                    "scenario_id": f"{scheme_dir.name}_{replicate_dir.name}",
                    "scenario_source": "resampling",
                    "scenario_class": f"resampling_{scheme_dir.name}",
                    "analysis_frame": scheme_dir.name,
                    "rooting_mode": "reference",
                    "tip_count": int(len(tips)),
                    "disrupted_tip_count": disrupted_tip_count,
                    "fitch_origin_events": fitch_origin_events,
                    "pastml_origin_events": np.nan,
                    "mk_origin_count_mean": np.nan,
                    "mk_origin_count_lower_95": np.nan,
                    "mk_origin_count_upper_95": np.nan,
                    "largest_disrupted_clade_share": safe_divide(max_disrupted, disrupted_tip_count),
                    "max_disrupted_tips_per_origin": max_disrupted,
                    "rejects_one_global_clone_fitch": bool(fitch_origin_events > 1),
                    "rejects_one_global_clone_mk95": np.nan,
                    "notes": f"Resampled balanced tree replicate from {scheme_dir.name}.",
                }
            )

    registry = pd.DataFrame(registry_rows).sort_values(["scenario_class", "scenario_id"]).reset_index(drop=True)
    summary = (
        registry.groupby("scenario_class", dropna=False)
        .agg(
            n_scenarios=("scenario_id", "nunique"),
            n_reject_one_global_clone_fitch=("rejects_one_global_clone_fitch", count_true_like),
            min_fitch_origin_events=("fitch_origin_events", "min"),
            median_fitch_origin_events=("fitch_origin_events", "median"),
            max_fitch_origin_events=("fitch_origin_events", "max"),
            median_largest_disrupted_clade_share=("largest_disrupted_clade_share", "median"),
            min_mk_lower_95=("mk_origin_count_lower_95", "min"),
            n_reject_one_global_clone_mk95=("rejects_one_global_clone_mk95", count_true_like),
        )
        .reset_index()
    )
    overall = pd.DataFrame(
        [
            {
                "scenario_class": "overall",
                "n_scenarios": int(registry["scenario_id"].nunique()),
                "n_reject_one_global_clone_fitch": count_true_like(registry["rejects_one_global_clone_fitch"]),
                "min_fitch_origin_events": registry["fitch_origin_events"].min(),
                "median_fitch_origin_events": registry["fitch_origin_events"].median(),
                "max_fitch_origin_events": registry["fitch_origin_events"].max(),
                "median_largest_disrupted_clade_share": registry["largest_disrupted_clade_share"].median(),
                "min_mk_lower_95": registry["mk_origin_count_lower_95"].min(),
                "n_reject_one_global_clone_mk95": count_true_like(registry["rejects_one_global_clone_mk95"]),
            }
        ]
    )
    summary = pd.concat([overall, summary], ignore_index=True, sort=False)
    return registry, summary


def clean_scalar(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "na"} else text


def truthy_scalar(value: object) -> bool:
    return clean_scalar(value).lower() in {"true", "1", "yes", "y", "t"}


def classify_origin_confidence_tier(row: pd.Series) -> tuple[int, str, str, str]:
    n_disrupted = pd.to_numeric(pd.Series([row.get("n_tips_disrupted")]), errors="coerce").fillna(0).iloc[0]
    representative_level = clean_scalar(row.get("representative_validation_level"))
    dominant_event_level = clean_scalar(row.get("dominant_event_validation_level"))
    origin_hard_anchor = truthy_scalar(row.get("origin_package_hard_anchor"))
    dominant_event_hard_anchor = truthy_scalar(row.get("dominant_event_hard_anchor"))
    package_anchored = representative_level not in {"", "assembly_only", "read_validation_unresolved", "other"}

    if n_disrupted <= 1:
        return (
            3,
            "Tier 3 singleton exploratory origin",
            "singleton_lower_confidence",
            "exploratory_ledger_not_core_package_support",
        )
    if origin_hard_anchor and package_anchored:
        return (
            1,
            "Tier 1 package-anchored non-singleton origin",
            "read_or_longread_anchored_non_singleton",
            "core_repeated_origin_package_support",
        )
    if origin_hard_anchor or dominant_event_hard_anchor or dominant_event_level not in {"", "assembly_only", "read_validation_unresolved", "other"}:
        return (
            2,
            "Tier 2 bounded non-singleton origin",
            "non_singleton_event_or_partial_anchor",
            "supporting_bounded_context",
        )
    return (
        2,
        "Tier 2 bounded non-singleton origin",
        "non_singleton_tree_only",
        "supporting_bounded_context",
    )


def build_origin_confidence_tier_table() -> pd.DataFrame:
    origins = pd.read_csv(FIG03_ORIGINS_PATH, sep="\t")
    validation = pd.read_csv(ARCHITECTURE_ORIGIN_VALIDATION_PATH, sep="\t")
    origins["origin_id"] = normalize_text(origins["origin_id"])
    validation["origin_id"] = normalize_text(validation["origin_id"])
    merged = origins.merge(
        validation[
            [
                "origin_id",
                "origin_is_major_package",
                "tree_representative_validation_level",
                "representative_validation_level",
                "dominant_event_validation_level",
                "origin_package_hard_anchor",
                "dominant_event_hard_anchor",
                "origin_package_event_anchored_only",
                "exemplar_replacement_applied",
                "validation_priority",
                "evidence_alignment",
                "followup_class",
                "public_data_recovery_status",
            ]
        ],
        on="origin_id",
        how="left",
    )
    merged["tree_representative_validation_level"] = merged["tree_representative_validation_level"].fillna(
        merged["validation_level"]
    )
    merged["representative_validation_level"] = merged["representative_validation_level"].fillna(
        merged["validation_level"]
    )
    merged["dominant_event_validation_level"] = merged["dominant_event_validation_level"].fillna("")
    for column in [
        "n_tips_total",
        "n_tips_disrupted",
        "n_countries",
        "first_year",
        "last_year",
        "branch_support",
        "origin_support_score",
    ]:
        merged[column] = pd.to_numeric(merged[column], errors="coerce")

    tier_rows: list[dict[str, object]] = []
    for row in merged.itertuples(index=False):
        row_series = pd.Series(row._asdict())
        tier_rank, tier_label, evidence_class, primary_use = classify_origin_confidence_tier(row_series)
        representative_level = clean_scalar(row_series.get("representative_validation_level")) or "unknown"
        tree_level = clean_scalar(row_series.get("tree_representative_validation_level")) or "unknown"
        dominant_event_level = clean_scalar(row_series.get("dominant_event_validation_level")) or "unknown"
        tier_rows.append(
            {
                "origin_id": row.origin_id,
                "origin_confidence_tier_rank": tier_rank,
                "origin_confidence_tier": tier_label,
                "origin_evidence_class": evidence_class,
                "primary_use_in_manuscript": primary_use,
                "n_tips_disrupted": row.n_tips_disrupted,
                "n_tips_total": row.n_tips_total,
                "n_countries": row.n_countries,
                "first_year": row.first_year,
                "last_year": row.last_year,
                "major_mlst_st": clean_scalar(row.major_mlst_st),
                "dominant_prn_mechanism": clean_scalar(row.dominant_prn_mechanism),
                "dominant_prn_event_id": clean_scalar(row.dominant_prn_event_id),
                "branch_support": row.branch_support,
                "origin_support_score": row.origin_support_score,
                "tree_representative_validation_level": tree_level,
                "representative_validation_level": representative_level,
                "dominant_event_validation_level": dominant_event_level,
                "origin_package_hard_anchor": truthy_scalar(row_series.get("origin_package_hard_anchor")),
                "dominant_event_hard_anchor": truthy_scalar(row_series.get("dominant_event_hard_anchor")),
                "exemplar_replacement_applied": truthy_scalar(row_series.get("exemplar_replacement_applied")),
                "validation_priority": clean_scalar(row_series.get("validation_priority")),
                "evidence_alignment": clean_scalar(row_series.get("evidence_alignment")),
                "public_data_recovery_status": clean_scalar(row_series.get("public_data_recovery_status")),
                "interpretation_note": (
                    "Non-singleton origin with package-level representative support."
                    if tier_rank == 1
                    else "Singleton origins remain a lower-confidence exploratory ledger."
                    if tier_rank == 3
                    else "Non-singleton origin retained as bounded support but not promoted to strongest package evidence."
                ),
            }
        )
    return pd.DataFrame(tier_rows).sort_values(
        ["origin_confidence_tier_rank", "n_tips_disrupted", "origin_id"],
        ascending=[True, False, True],
    ).reset_index(drop=True)


def build_year_sensitivity(
    epochs: pd.DataFrame,
    year_comp: pd.DataFrame,
    ipw: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    leave_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    ipw_subset = ipw[["country_iso3", "year", "ipw_weight_total", "ipw_weighted_disrupted"]].copy()
    ipw_subset["ipw_weight_total"] = pd.to_numeric(ipw_subset["ipw_weight_total"], errors="coerce").fillna(0.0)
    ipw_subset["ipw_weighted_disrupted"] = pd.to_numeric(
        ipw_subset["ipw_weighted_disrupted"], errors="coerce"
    ).fillna(0.0)

    for country_iso3, (pre_epoch, post_epoch) in PRIMARY_COUNTRY_PAIRS.items():
        country_epochs = epochs.loc[
            (epochs["country_iso3"] == country_iso3)
            & (epochs["epoch_id"].isin([pre_epoch, post_epoch]))
        ].copy()
        if len(country_epochs) != 2:
            continue
        country_name = country_epochs["country_name"].iloc[0]

        yc = year_comp.loc[
            (year_comp["country_iso3"] == country_iso3)
            & (year_comp["epoch_id"].isin([pre_epoch, post_epoch]))
        ].copy()
        yc["n_prn_interpretable"] = pd.to_numeric(yc["n_prn_interpretable"], errors="coerce").fillna(0).astype(int)
        yc["n_prn_disrupted"] = pd.to_numeric(yc["n_prn_disrupted"], errors="coerce").fillna(0).astype(int)
        yc = yc.merge(ipw_subset, on=["country_iso3", "year"], how="left")
        yc["ipw_weight_total"] = yc["ipw_weight_total"].fillna(0.0)
        yc["ipw_weighted_disrupted"] = yc["ipw_weighted_disrupted"].fillna(0.0)

        epoch_stats: dict[str, dict[str, float]] = {}
        top_years: dict[str, int | float] = {}
        for epoch_id in [pre_epoch, post_epoch]:
            epoch_df = yc.loc[yc["epoch_id"] == epoch_id].copy()
            epoch_df = epoch_df.loc[epoch_df["n_prn_interpretable"] > 0].copy()
            total_interpretable = int(epoch_df["n_prn_interpretable"].sum())
            total_disrupted = int(epoch_df["n_prn_disrupted"].sum())
            total_weight = float(epoch_df["ipw_weight_total"].sum())
            total_weighted_disrupted = float(epoch_df["ipw_weighted_disrupted"].sum())
            top_share = (
                float(epoch_df["share_of_epoch_interpretable"].max())
                if "share_of_epoch_interpretable" in epoch_df.columns and not epoch_df.empty
                else np.nan
            )
            top_year = (
                int(epoch_df.sort_values("share_of_epoch_interpretable", ascending=False)["year"].iloc[0])
                if not epoch_df.empty
                else np.nan
            )
            epoch_stats[epoch_id] = {
                "n_years": int(epoch_df["year"].nunique()),
                "n_interpretable": total_interpretable,
                "n_disrupted": total_disrupted,
                "naive": safe_divide(total_disrupted, total_interpretable),
                "ipw": safe_divide(total_weighted_disrupted, total_weight),
                "ipw_total": total_weight,
                "ipw_disrupted": total_weighted_disrupted,
                "top_share": top_share,
                "top_year": top_year,
            }
            top_years[epoch_id] = top_year

        baseline_naive_delta = epoch_stats[post_epoch]["naive"] - epoch_stats[pre_epoch]["naive"]
        baseline_ipw_delta = epoch_stats[post_epoch]["ipw"] - epoch_stats[pre_epoch]["ipw"]
        baseline_naive_sign = choose_sign(baseline_naive_delta)
        baseline_ipw_sign = choose_sign(baseline_ipw_delta)

        for epoch_id in [pre_epoch, post_epoch]:
            epoch_df = yc.loc[(yc["epoch_id"] == epoch_id) & (yc["n_prn_interpretable"] > 0)].copy()
            for row in epoch_df.itertuples(index=False):
                adjusted = {key: dict(value) for key, value in epoch_stats.items()}
                adjusted[epoch_id]["n_interpretable"] -= int(row.n_prn_interpretable)
                adjusted[epoch_id]["n_disrupted"] -= int(row.n_prn_disrupted)
                adjusted[epoch_id]["ipw_total"] -= float(row.ipw_weight_total)
                adjusted[epoch_id]["ipw_disrupted"] -= float(row.ipw_weighted_disrupted)
                adjusted[epoch_id]["naive"] = safe_divide(
                    adjusted[epoch_id]["n_disrupted"], adjusted[epoch_id]["n_interpretable"]
                )
                adjusted[epoch_id]["ipw"] = safe_divide(
                    adjusted[epoch_id]["ipw_disrupted"], adjusted[epoch_id]["ipw_total"]
                )
                leave_naive_delta = adjusted[post_epoch]["naive"] - adjusted[pre_epoch]["naive"]
                leave_ipw_delta = adjusted[post_epoch]["ipw"] - adjusted[pre_epoch]["ipw"]
                leave_naive_sign = choose_sign(leave_naive_delta)
                leave_ipw_sign = choose_sign(leave_ipw_delta)
                leave_rows.append(
                    {
                        "country_iso3": country_iso3,
                        "country_name": country_name,
                        "comparison_pre_epoch_id": pre_epoch,
                        "comparison_post_epoch_id": post_epoch,
                        "dropped_epoch_id": epoch_id,
                        "dropped_year": int(row.year),
                        "dropped_share_of_epoch_interpretable": float(row.share_of_epoch_interpretable),
                        "dropped_share_of_epoch_disrupted": float(row.share_of_epoch_disrupted),
                        "baseline_naive_delta": baseline_naive_delta,
                        "baseline_naive_sign": baseline_naive_sign,
                        "leave_one_year_naive_delta": leave_naive_delta,
                        "leave_one_year_naive_sign": leave_naive_sign,
                        "naive_sign_changed": leave_naive_sign != baseline_naive_sign,
                        "baseline_ipw_delta": baseline_ipw_delta,
                        "baseline_ipw_sign": baseline_ipw_sign,
                        "leave_one_year_ipw_delta": leave_ipw_delta,
                        "leave_one_year_ipw_sign": leave_ipw_sign,
                        "ipw_sign_changed": leave_ipw_sign != baseline_ipw_sign,
                        "absolute_naive_delta_change": (
                            np.nan
                            if pd.isna(leave_naive_delta) or pd.isna(baseline_naive_delta)
                            else abs(leave_naive_delta - baseline_naive_delta)
                        ),
                        "absolute_ipw_delta_change": (
                            np.nan
                            if pd.isna(leave_ipw_delta) or pd.isna(baseline_ipw_delta)
                            else abs(leave_ipw_delta - baseline_ipw_delta)
                        ),
                    }
                )

        country_leave = pd.DataFrame([row for row in leave_rows if row["country_iso3"] == country_iso3])
        top_drop_naive_changed = False
        top_drop_ipw_changed = False
        top_drop_rows: dict[str, pd.Series | None] = {}
        for epoch_id in [pre_epoch, post_epoch]:
            top_year = top_years.get(epoch_id)
            if pd.isna(top_year):
                top_drop_rows[epoch_id] = None
                continue
            top_row = country_leave.loc[
                (country_leave["dropped_epoch_id"] == epoch_id) & (country_leave["dropped_year"] == int(top_year))
            ]
            if not top_row.empty:
                top_drop_rows[epoch_id] = top_row.iloc[0]
                top_drop_naive_changed = top_drop_naive_changed or bool(top_row["naive_sign_changed"].iloc[0])
                top_drop_ipw_changed = top_drop_ipw_changed or bool(top_row["ipw_sign_changed"].iloc[0])
            else:
                top_drop_rows[epoch_id] = None

        pre_top_row = top_drop_rows.get(pre_epoch)
        post_top_row = top_drop_rows.get(post_epoch)

        summary_rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": country_name,
                "comparison_pre_epoch_id": pre_epoch,
                "comparison_post_epoch_id": post_epoch,
                "baseline_naive_delta": baseline_naive_delta,
                "baseline_naive_sign": baseline_naive_sign,
                "baseline_ipw_delta": baseline_ipw_delta,
                "baseline_ipw_sign": baseline_ipw_sign,
                "pre_epoch_years_with_interpretable": epoch_stats[pre_epoch]["n_years"],
                "post_epoch_years_with_interpretable": epoch_stats[post_epoch]["n_years"],
                "pre_epoch_max_single_year_share": epoch_stats[pre_epoch]["top_share"],
                "post_epoch_max_single_year_share": epoch_stats[post_epoch]["top_share"],
                "pre_epoch_top_year": epoch_stats[pre_epoch]["top_year"],
                "post_epoch_top_year": epoch_stats[post_epoch]["top_year"],
                "pre_epoch_top_year_drop_naive_delta": (
                    np.nan if pre_top_row is None else pre_top_row["leave_one_year_naive_delta"]
                ),
                "post_epoch_top_year_drop_naive_delta": (
                    np.nan if post_top_row is None else post_top_row["leave_one_year_naive_delta"]
                ),
                "pre_epoch_top_year_drop_ipw_delta": (
                    np.nan if pre_top_row is None else pre_top_row["leave_one_year_ipw_delta"]
                ),
                "post_epoch_top_year_drop_ipw_delta": (
                    np.nan if post_top_row is None else post_top_row["leave_one_year_ipw_delta"]
                ),
                "pre_epoch_top_year_drop_naive_sign": (
                    "not_estimable" if pre_top_row is None else pre_top_row["leave_one_year_naive_sign"]
                ),
                "post_epoch_top_year_drop_naive_sign": (
                    "not_estimable" if post_top_row is None else post_top_row["leave_one_year_naive_sign"]
                ),
                "pre_epoch_top_year_drop_ipw_sign": (
                    "not_estimable" if pre_top_row is None else pre_top_row["leave_one_year_ipw_sign"]
                ),
                "post_epoch_top_year_drop_ipw_sign": (
                    "not_estimable" if post_top_row is None else post_top_row["leave_one_year_ipw_sign"]
                ),
                "any_leave_one_year_naive_sign_change": bool(country_leave["naive_sign_changed"].any()),
                "any_leave_one_year_ipw_sign_change": bool(country_leave["ipw_sign_changed"].any()),
                "max_leave_one_year_naive_delta_change": country_leave["absolute_naive_delta_change"].max(),
                "max_leave_one_year_ipw_delta_change": country_leave["absolute_ipw_delta_change"].max(),
                "top_year_drop_changes_naive_sign": top_drop_naive_changed,
                "top_year_drop_changes_ipw_sign": top_drop_ipw_changed,
            }
        )

    return (
        pd.DataFrame(summary_rows).sort_values("country_iso3").reset_index(drop=True),
        pd.DataFrame(leave_rows).sort_values(["country_iso3", "dropped_epoch_id", "dropped_year"]).reset_index(drop=True),
    )


def assign_epochs_to_manifest(manifest: pd.DataFrame, history: pd.DataFrame) -> pd.DataFrame:
    selected_history = history.loc[history["country_iso3"].isin(PRIMARY_COUNTRY_PAIRS.keys())].copy()
    frames: list[pd.DataFrame] = []
    for item in selected_history.itertuples(index=False):
        subset = manifest.loc[
            (manifest["country_iso3"] == item.country_iso3)
            & manifest["year"].between(item.start_year, item.end_year, inclusive="both")
        ].copy()
        if subset.empty:
            continue
        subset["epoch_id"] = item.epoch_id
        subset["epoch_label"] = getattr(item, "epoch_label", item.epoch_id)
        subset["epoch_type"] = getattr(item, "epoch_type", "")
        frames.append(subset)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def build_architecture_turnover(
    manifest: pd.DataFrame,
    history: pd.DataFrame,
    structure_reuse: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    epoch_manifest = assign_epochs_to_manifest(manifest, history)
    if epoch_manifest.empty:
        return pd.DataFrame(), pd.DataFrame()
    country_names = history[["country_iso3", "country_name"]].drop_duplicates()
    epoch_totals_all = (
        epoch_manifest.loc[epoch_manifest["prn_disrupted"]]
        .groupby(["country_iso3", "epoch_id"], dropna=False)
        .size()
        .rename("epoch_total_disrupted")
        .reset_index()
    )
    event_meta = (
        structure_reuse[
            [
                "prn_event_id",
                "prn_mechanism_call",
                "major_recurrent_architecture",
                "dominant_after_origin_collapse",
                "validation_level",
                "hard_anchor",
            ]
        ]
        .drop_duplicates(subset=["prn_event_id"])
        .copy()
    )

    disrupted = epoch_manifest.loc[epoch_manifest["prn_disrupted"]].copy()
    disrupted["prn_event_id"] = normalize_text(disrupted["prn_event_id"])
    disrupted["prn_mechanism_call"] = normalize_text(disrupted["prn_mechanism_call"])

    long_table = (
        disrupted.groupby(
            ["country_iso3", "epoch_id", "epoch_label", "epoch_type", "prn_event_id", "prn_mechanism_call"],
            dropna=False,
        )
        .size()
        .rename("epoch_event_count")
        .reset_index()
    )
    long_table = long_table.merge(epoch_totals_all, on=["country_iso3", "epoch_id"], how="left")
    long_table = long_table.merge(country_names, on="country_iso3", how="left")
    long_table["epoch_event_share"] = long_table["epoch_event_count"] / long_table["epoch_total_disrupted"]
    long_table = long_table.merge(event_meta, on=["prn_event_id", "prn_mechanism_call"], how="left")
    long_table["event_label"] = long_table["prn_event_id"].map(event_label)
    long_table["dominant_epoch_architecture"] = (
        long_table.groupby(["country_iso3", "epoch_id"])["epoch_event_count"]
        .transform(lambda s: s == s.max())
        .astype(bool)
    )

    summary_rows: list[dict[str, object]] = []
    for country_iso3, frame in long_table.groupby("country_iso3", dropna=False):
        country_name = frame["country_name"].dropna().iloc[0] if frame["country_name"].notna().any() else ""
        epoch_vectors: dict[str, pd.Series] = {}
        epoch_frames: dict[str, pd.DataFrame] = {}
        for epoch_id, epoch_frame in frame.groupby("epoch_id", dropna=False):
            epoch_frame = epoch_frame.sort_values(["epoch_event_count", "prn_event_id"], ascending=[False, True]).copy()
            epoch_frames[epoch_id] = epoch_frame
            epoch_vectors[epoch_id] = epoch_frame.set_index("prn_event_id")["epoch_event_share"]

        compare_pairs: list[tuple[str, str]] = []
        if country_iso3 in PRIMARY_COUNTRY_PAIRS:
            compare_pairs.append(PRIMARY_COUNTRY_PAIRS[country_iso3])
        epoch_order = (
            history.loc[history["country_iso3"] == country_iso3, ["epoch_id", "start_year"]]
            .drop_duplicates()
            .sort_values(["start_year", "epoch_id"])["epoch_id"]
            .tolist()
        )
        compare_pairs.extend(list(zip(epoch_order[:-1], epoch_order[1:])))
        seen_pairs: set[tuple[str, str]] = set()

        for prev_epoch, next_epoch in compare_pairs:
            if (prev_epoch, next_epoch) in seen_pairs:
                continue
            seen_pairs.add((prev_epoch, next_epoch))
            prev_frame = epoch_frames.get(prev_epoch, pd.DataFrame())
            next_frame = epoch_frames.get(next_epoch, pd.DataFrame())
            prev_vector_raw = epoch_vectors.get(prev_epoch, pd.Series(dtype=float))
            next_vector_raw = epoch_vectors.get(next_epoch, pd.Series(dtype=float))
            all_events = sorted(set(prev_vector_raw.index) | set(next_vector_raw.index))
            if not all_events:
                all_events = [""]
            prev_vector = prev_vector_raw.reindex(all_events, fill_value=0.0)
            next_vector = next_vector_raw.reindex(all_events, fill_value=0.0)
            if prev_frame.empty:
                prev_dom = {
                    "prn_event_id": "none_observed",
                    "event_label": "None observed",
                    "prn_mechanism_call": "none_observed",
                    "epoch_event_share": 0.0,
                }
            else:
                prev_dom = prev_frame.iloc[0]
            if next_frame.empty:
                next_dom = {
                    "prn_event_id": "none_observed",
                    "event_label": "None observed",
                    "prn_mechanism_call": "none_observed",
                    "epoch_event_share": 0.0,
                }
            else:
                next_dom = next_frame.iloc[0]
            prev_total = epoch_totals_all.loc[
                (epoch_totals_all["country_iso3"] == country_iso3) & (epoch_totals_all["epoch_id"] == prev_epoch),
                "epoch_total_disrupted",
            ]
            next_total = epoch_totals_all.loc[
                (epoch_totals_all["country_iso3"] == country_iso3) & (epoch_totals_all["epoch_id"] == next_epoch),
                "epoch_total_disrupted",
            ]
            summary_rows.append(
                {
                    "country_iso3": country_iso3,
                    "country_name": country_name,
                    "previous_epoch_id": prev_epoch,
                    "next_epoch_id": next_epoch,
                    "previous_epoch_label": (
                        prev_frame["epoch_label"].iloc[0] if not prev_frame.empty else ""
                    ),
                    "next_epoch_label": (
                        next_frame["epoch_label"].iloc[0] if not next_frame.empty else ""
                    ),
                    "comparison_type": (
                        "primary_selected_country_pair"
                        if country_iso3 in PRIMARY_COUNTRY_PAIRS and PRIMARY_COUNTRY_PAIRS[country_iso3] == (prev_epoch, next_epoch)
                        else "adjacent_epoch_pair"
                    ),
                    "previous_dominant_event_id": prev_dom["prn_event_id"],
                    "previous_dominant_event_label": prev_dom["event_label"],
                    "previous_dominant_mechanism": prev_dom["prn_mechanism_call"],
                    "previous_dominant_share": prev_dom["epoch_event_share"],
                    "next_dominant_event_id": next_dom["prn_event_id"],
                    "next_dominant_event_label": next_dom["event_label"],
                    "next_dominant_mechanism": next_dom["prn_mechanism_call"],
                    "next_dominant_share": next_dom["epoch_event_share"],
                    "dominant_event_changed": prev_dom["prn_event_id"] != next_dom["prn_event_id"],
                    "dominant_mechanism_changed": prev_dom["prn_mechanism_call"] != next_dom["prn_mechanism_call"],
                    "architecture_total_variation_distance": total_variation_distance(prev_vector, next_vector),
                    "n_events_previous_epoch": int(len(prev_frame)),
                    "n_events_next_epoch": int(len(next_frame)),
                    "epoch_total_disrupted_previous": int(prev_total.iloc[0]) if not prev_total.empty else 0,
                    "epoch_total_disrupted_next": int(next_total.iloc[0]) if not next_total.empty else 0,
                }
            )

    long_out = long_table.sort_values(
        ["country_iso3", "epoch_id", "epoch_event_count"], ascending=[True, True, False]
    ).reset_index(drop=True)
    summary_out = pd.DataFrame(summary_rows)
    if not summary_out.empty:
        summary_out = summary_out.sort_values(["country_iso3", "previous_epoch_id"]).reset_index(drop=True)
    return (long_out, summary_out)


def build_origin_bridge(
    origin_packages: pd.DataFrame,
    origin_shift: pd.DataFrame,
    detection_shift: pd.DataFrame,
    evidence_grid: pd.DataFrame,
) -> pd.DataFrame:
    origin_packages = origin_packages.copy()
    for column in ["n_disrupted_descendants", "n_total_descendants", "follow_up_years", "origin_package_hard_anchor"]:
        if column in origin_packages.columns:
            if column == "origin_package_hard_anchor":
                origin_packages[column] = parse_bool(origin_packages[column])
            else:
                origin_packages[column] = pd.to_numeric(origin_packages[column], errors="coerce")

    origin_burden = (
        origin_packages.groupby("origin_country_iso3", dropna=False)
        .agg(
            n_local_origin_packages=("origin_id", "nunique"),
            n_origin_packages_ge3_descendants=("established_ge3_descendants", lambda s: pd.to_numeric(s, errors="coerce").fillna(0).sum()),
            total_package_descendants=("n_disrupted_descendants", "sum"),
            max_package_descendants=("n_disrupted_descendants", "max"),
            median_package_follow_up_years=("follow_up_years", "median"),
            n_hard_anchor_packages=("origin_package_hard_anchor", "sum"),
        )
        .reset_index()
        .rename(columns={"origin_country_iso3": "country_iso3"})
    )
    scaffold = (
        evidence_grid[
            [
                "country_iso3",
                "country_name",
                "amplification_pattern",
                "final_interpretation_tier",
                "prevalence_direction",
                "bounds_stability",
            ]
        ]
        .drop_duplicates(subset=["country_iso3"])
        .reset_index(drop=True)
        .copy()
    )
    scaffold["country_order"] = np.arange(len(scaffold))

    origin_shift_use = origin_shift[
        [
            "country_iso3",
            "post_minus_pre_ipw_prevalence",
            "amplification_pattern",
            "peak_origin_clade_descendants",
            "peak_origin_clades_active",
        ]
    ].copy()
    origin_shift_use["shift_source"] = "first_local_origin"
    detection_shift_use = detection_shift[
        [
            "country_iso3",
            "post_minus_pre_ipw_prevalence",
            "amplification_pattern",
            "peak_origin_clade_descendants",
            "peak_origin_clades_active",
        ]
    ].copy()
    detection_shift_use["shift_source"] = "first_prn_detection"

    preferred_shift = (
        pd.concat([origin_shift_use, detection_shift_use], ignore_index=True)
        .sort_values(["country_iso3", "shift_source"], key=lambda s: s.map({"first_local_origin": 0, "first_prn_detection": 1}))
        .drop_duplicates(subset=["country_iso3"], keep="first")
    )

    output = scaffold.merge(origin_burden, on="country_iso3", how="left").merge(
        preferred_shift,
        on="country_iso3",
        how="left",
        suffixes=("_evidence_grid", ""),
    )
    for column in [
        "n_local_origin_packages",
        "n_origin_packages_ge3_descendants",
        "total_package_descendants",
        "max_package_descendants",
        "n_hard_anchor_packages",
    ]:
        output[column] = pd.to_numeric(output[column], errors="coerce").fillna(0).astype(int)
    output["has_local_origin_package"] = output["n_local_origin_packages"] > 0
    output["shift_source"] = output["shift_source"].fillna("not_estimable_no_anchor")
    output["amplification_pattern"] = (
        output["amplification_pattern"]
        .fillna(output["amplification_pattern_evidence_grid"])
        .fillna("not_observed")
    )
    output["packages_per_peak_descendant"] = output.apply(
        lambda row: safe_divide(row["n_local_origin_packages"], row["peak_origin_clade_descendants"]), axis=1
    )
    output = output.drop(columns=["amplification_pattern_evidence_grid"])
    return output.sort_values(["country_order", "n_local_origin_packages"], ascending=[True, False]).drop(
        columns=["country_order"]
    ).reset_index(drop=True)


def build_negative_control(
    manifest: pd.DataFrame,
    history: pd.DataFrame,
    step2: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    marker_status_columns = [
        f"marker_status_{item['locus']}"
        for item in LOCUS_CONFIG
        if item["locus"] != "prn"
    ]
    step2_columns = [
        "Assembly BioSample Accession",
        "Current Accession",
        "country",
        "year",
    ] + [column for column in marker_status_columns if column in step2.columns]
    step2_use = step2[step2_columns].copy()
    step2_use = step2_use.rename(
        columns={
            "Assembly BioSample Accession": "biosample_accession",
            "Current Accession": "step2_current_accession",
        }
    )
    step2_use["biosample_accession"] = normalize_text(step2_use["biosample_accession"])
    step2_use["step2_current_accession"] = normalize_text(step2_use["step2_current_accession"])
    step2_use = step2_use.loc[step2_use["biosample_accession"] != ""].drop_duplicates(subset=["biosample_accession"])
    if PSEUDO_CONTROL_STATUS_PATH.exists() and PSEUDO_CONTROL_STATUS_PATH.stat().st_size > 0:
        pseudo_status = pd.read_csv(PSEUDO_CONTROL_STATUS_PATH, sep="\t", dtype=str)
        pseudo_status["biosample_accession"] = normalize_text(pseudo_status.get("biosample_accession", pd.Series(dtype=str)))
        pseudo_status_columns = [column for column in marker_status_columns if column in pseudo_status.columns]
        if pseudo_status_columns:
            pseudo_status = (
                pseudo_status[["biosample_accession"] + pseudo_status_columns]
                .loc[pseudo_status["biosample_accession"] != ""]
                .drop_duplicates(subset=["biosample_accession"])
            )
            step2_use = step2_use.merge(pseudo_status, on="biosample_accession", how="left", suffixes=("", "_pseudo"))
            for column in pseudo_status_columns:
                pseudo_column = f"{column}_pseudo"
                if pseudo_column not in step2_use.columns:
                    continue
                existing = normalize_text(step2_use.get(column, pd.Series(index=step2_use.index, dtype=str)))
                incoming = normalize_text(step2_use[pseudo_column])
                step2_use[column] = existing.where(existing.ne("") & existing.ne("missing"), incoming)
                step2_use = step2_use.drop(columns=[pseudo_column])
    manifest_use = manifest.copy()
    manifest_use["biosample_accession"] = normalize_text(manifest_use["biosample_accession"])
    manifest_use = manifest_use.loc[manifest_use["biosample_accession"] != ""].drop_duplicates(subset=["biosample_accession"])
    for column in marker_status_columns:
        step2_use[column] = normalize_text(step2_use.get(column, pd.Series(index=step2_use.index, dtype=str))).replace("", "missing")

    overlap = manifest_use.merge(step2_use, on="biosample_accession", how="inner")
    if "year_x" in overlap.columns and "year" not in overlap.columns:
        overlap = overlap.rename(columns={"year_x": "year"})
    if "country_iso3_x" in overlap.columns and "country_iso3" not in overlap.columns:
        overlap = overlap.rename(columns={"country_iso3_x": "country_iso3"})
    overlap = overlap.drop_duplicates(subset=["biosample_accession"]).copy()
    country_name_map = (
        history[["country_iso3", "country_name"]].drop_duplicates().set_index("country_iso3")["country_name"].to_dict()
    )

    global_rows: list[dict[str, object]] = []
    epoch_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    overlap_selected = assign_epochs_to_manifest(overlap, history)
    if overlap_selected.empty:
        overlap_selected = pd.DataFrame(columns=overlap.columns.tolist() + ["epoch_id", "epoch_label", "epoch_type"])

    for item in LOCUS_CONFIG:
        locus = item["locus"]
        label = item["label"]
        locus_length_bp = int(item["locus_length_bp"])
        locus_length_kb = locus_length_bp / 1000.0
        if locus == "prn":
            interpretable = parse_bool(overlap["prn_interpretable"])
            signal_positive = parse_bool(overlap["prn_disrupted"])
            noninterpretable = ~interpretable
        else:
            status = normalize_text(overlap.get(f"marker_status_{locus}", pd.Series(dtype=str))).replace("", "missing")
            interpretable = status.eq("ok")
            signal_positive = status.eq("below_threshold")
            noninterpretable = ~interpretable

        signal_frame = overlap.loc[signal_positive].copy()
        global_rows.append(
            {
                "locus": locus,
                "locus_label": label,
                "locus_category": item["category"],
                "signal_definition": item["signal_definition"],
                "locus_length_bp": locus_length_bp,
                "n_overlap_genomes": int(len(overlap)),
                "n_locus_interpretable": int(interpretable.sum()),
                "n_signal_positive": int(signal_positive.sum()),
                "signal_positive_fraction_among_interpretable": safe_divide(signal_positive.sum(), interpretable.sum()),
                "signal_rate_per_kb_interpretable_genome": safe_divide(
                    signal_positive.sum(),
                    interpretable.sum() * locus_length_kb,
                ),
                "signal_rate_per_1000_interpretable_genomes": safe_divide(
                    signal_positive.sum() * 1000,
                    interpretable.sum(),
                ),
                "n_noninterpretable_or_missing": int(noninterpretable.sum()),
                "n_countries_with_signal": int(signal_frame["country_iso3"].nunique()),
                "n_years_with_signal": int(signal_frame["year"].nunique()),
                "recurrent_signal_flag": bool(signal_positive.sum() >= 5 and signal_frame["country_iso3"].nunique() >= 2),
            }
        )

        global_summary = (
            "recurrent_structural_signal_detected"
            if (signal_positive.sum() >= 5 and signal_frame["country_iso3"].nunique() >= 2)
            else "no_recurrent_signal_detected"
        )
        if locus == "prn":
            interpretation_note = "PRN shows the expected recurrent coding-disruption signal on the overlap frame."
        elif item["category"] in {
            "structure_matched_autotransporter",
            "pertactin_homologous_autotransporter",
        }:
            interpretation_note = (
                "A structure-matched pseudo-control shows a recurrent non-full-length marker signal; interpret as a "
                "possible locus-specific artifact/background signal rather than PRN-specific evidence."
                if global_summary == "recurrent_structural_signal_detected"
                else "No PRN-like recurrent non-full-length signal was detected in this structure-matched pseudo-control."
            )
        elif item["category"] == "secondary_autotransporter_reference_pseudogene_caveat":
            interpretation_note = (
                "Secondary autotransporter control; Tohama-I carries a bapC pseudogene caveat, so non-full-length calls "
                "are reported but not treated as primary PRN-locus specificity evidence."
            )
        elif item["category"] == "acellular_antigen":
            interpretation_note = (
                "No PRN-like recurrent non-full-length signal was detected in this same-Step2 antigen control."
            )
        else:
            interpretation_note = (
                "Secondary large vaccine-antigen/adhesin control; reported as a size- and surface-exposure comparator "
                "rather than a direct autotransporter-family match. Any non-full-length calls should be interpreted "
                "against their much lower rate than PRN."
            )
        summary_rows.append(
            {
                "locus": locus,
                "locus_label": label,
                "signal_definition": item["signal_definition"],
                "global_signal_summary": global_summary,
                "interpretation_note": interpretation_note,
            }
        )

        if overlap_selected.empty:
            continue
        if locus == "prn":
            overlap_selected["_interpretable"] = parse_bool(overlap_selected["prn_interpretable"])
            overlap_selected["_signal_positive"] = parse_bool(overlap_selected["prn_disrupted"])
        else:
            status = normalize_text(
                overlap_selected.get(f"marker_status_{locus}", pd.Series(index=overlap_selected.index, dtype=str))
            ).replace("", "missing")
            overlap_selected["_interpretable"] = status.eq("ok")
            overlap_selected["_signal_positive"] = status.eq("below_threshold")
        grouped = (
            overlap_selected.groupby(["country_iso3", "epoch_id", "epoch_label", "epoch_type"], dropna=False)
            .agg(
                n_total_overlap_genomes=("biosample_accession", "nunique"),
                n_locus_interpretable=("_interpretable", "sum"),
                n_signal_positive=("_signal_positive", "sum"),
            )
            .reset_index()
        )
        grouped["n_noninterpretable_or_missing"] = grouped["n_total_overlap_genomes"] - grouped["n_locus_interpretable"]
        grouped["signal_positive_fraction_among_interpretable"] = grouped.apply(
            lambda row: safe_divide(row["n_signal_positive"], row["n_locus_interpretable"]),
            axis=1,
        )
        grouped["country_name"] = grouped["country_iso3"].map(country_name_map)
        grouped["locus"] = locus
        grouped["locus_label"] = label
        grouped["locus_category"] = item["category"]
        grouped["signal_definition"] = item["signal_definition"]
        grouped["locus_length_bp"] = locus_length_bp
        grouped["n_overlap_genomes"] = int(len(overlap))
        grouped["signal_rate_per_kb_interpretable_genome"] = grouped.apply(
            lambda row: safe_divide(
                row["n_signal_positive"],
                row["n_locus_interpretable"] * locus_length_kb,
            ),
            axis=1,
        )
        grouped["signal_rate_per_1000_interpretable_genomes"] = grouped.apply(
            lambda row: safe_divide(row["n_signal_positive"] * 1000, row["n_locus_interpretable"]),
            axis=1,
        )
        grouped["n_countries_with_signal"] = int(signal_frame["country_iso3"].nunique())
        grouped["n_years_with_signal"] = int(signal_frame["year"].nunique())
        grouped["recurrent_signal_flag"] = bool(signal_positive.sum() >= 5 and signal_frame["country_iso3"].nunique() >= 2)
        grouped["global_signal_summary"] = global_summary
        grouped["interpretation_note"] = interpretation_note
        epoch_rows.extend(grouped.to_dict(orient="records"))

    return (
        pd.DataFrame(global_rows).sort_values("locus").reset_index(drop=True),
        pd.DataFrame(epoch_rows).sort_values(["locus", "country_iso3", "epoch_id"]).reset_index(drop=True),
        pd.DataFrame(summary_rows).sort_values("locus").reset_index(drop=True),
    )


def write_dual(df: pd.DataFrame, main_path: Path, supp_path: Path) -> None:
    main_path.parent.mkdir(parents=True, exist_ok=True)
    supp_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(main_path, sep="\t", index=False)
    df.to_csv(supp_path, sep="\t", index=False)


def main() -> None:
    data = load_inputs()

    interpretability_table = build_interpretability_model_table(data["manifest"])
    threshold_summary, threshold_grid = build_threshold_robustness(
        epoch_eligibility=data["epoch_eligibility"],
        scorecard=data["scorecard"],
        evidence_grid=data["evidence_grid"],
        history=data["history"],
    )
    updated_scorecard = build_selection_scorecard_multiverse_summary(
        scorecard=data["scorecard"],
        epoch_eligibility=data["epoch_eligibility"],
        history=data["history"],
        country_grid=threshold_grid,
    )
    tipping_summary, tipping_grid = build_missingness_tipping_summary(
        epoch_eligibility=data["epoch_eligibility"],
    )
    dr_summary, updated_contrasts = build_missingness_dr_summary(
        manifest=data["manifest"],
        history=data["history"],
        contrasts=data["contrasts"],
        tipping_summary=tipping_summary,
    )
    asr_registry, asr_one_global_clone_summary = build_asr_scenario_registry()
    origin_confidence_tiers = build_origin_confidence_tier_table()

    year_summary, leave_one_year = build_year_sensitivity(
        epochs=data["epochs"],
        year_comp=data["year_comp"],
        ipw=data["ipw"],
    )
    turnover_long, turnover_summary = build_architecture_turnover(
        manifest=data["manifest"],
        history=data["history"],
        structure_reuse=data["structure_reuse"],
    )
    origin_bridge = build_origin_bridge(
        origin_packages=data["origin_packages"],
        origin_shift=data["origin_shift"],
        detection_shift=data["detection_shift"],
        evidence_grid=data["evidence_grid"],
    )
    antigen_global, antigen_epoch, antigen_summary = build_negative_control(
        manifest=data["manifest"],
        history=data["history"],
        step2=data["step2"],
    )
    quality_sensitivity = build_quality_restricted_sensitivity(
        manifest=data["manifest"],
        history=data["history"],
    )

    YEAR_SENS_SUMMARY_PATH.parent.mkdir(parents=True, exist_ok=True)
    year_summary.to_csv(YEAR_SENS_SUMMARY_PATH, sep="\t", index=False)
    turnover_summary.to_csv(TURNOVER_SUMMARY_PATH, sep="\t", index=False)
    updated_scorecard.to_csv(SELECTION_SCORECARD_PATH, sep="\t", index=False)
    updated_contrasts.to_csv(EPOCH_CONTRAST_PATH, sep="\t", index=False)
    write_dual(leave_one_year, YEAR_SENS_LOYO_PATH, YEAR_SENS_SUPP_PATH)
    write_dual(turnover_long, TURNOVER_LONG_PATH, TURNOVER_SUPP_PATH)
    write_dual(origin_bridge, ORIGIN_BRIDGE_PATH, ORIGIN_BRIDGE_SUPP_PATH)
    write_dual(
        interpretability_table,
        MISSINGNESS_INTERPRETABILITY_PATH,
        MISSINGNESS_INTERPRETABILITY_SUPP_PATH,
    )
    write_dual(
        quality_sensitivity,
        QUALITY_SENSITIVITY_PATH,
        QUALITY_SENSITIVITY_SUPP_PATH,
    )
    write_dual(
        threshold_summary,
        READINESS_THRESHOLD_GRID_PATH.with_name("selected_country_threshold_robustness_summary.tsv"),
        READINESS_THRESHOLD_SUPP_PATH,
    )
    threshold_grid.to_csv(READINESS_THRESHOLD_GRID_PATH, sep="\t", index=False)
    write_dual(
        dr_summary,
        MISSINGNESS_DR_SUMMARY_PATH,
        MISSINGNESS_DR_SUPP_PATH,
    )
    write_dual(
        tipping_summary,
        MISSINGNESS_TIPPING_GRID_PATH.with_name("selected_country_missingness_tipping_summary.tsv"),
        MISSINGNESS_TIPPING_SUPP_PATH,
    )
    tipping_grid.to_csv(MISSINGNESS_TIPPING_GRID_PATH, sep="\t", index=False)
    write_dual(
        asr_registry,
        ASR_SCENARIO_REGISTRY_PATH,
        ASR_SCENARIO_REGISTRY_SUPP_PATH,
    )
    asr_registry.to_csv(ASR_SCENARIO_REGISTRY_READER_SUPP_PATH, sep="\t", index=False)
    ASR_ONE_GLOBAL_CLONE_SUMMARY_PATH.parent.mkdir(parents=True, exist_ok=True)
    asr_one_global_clone_summary.to_csv(ASR_ONE_GLOBAL_CLONE_SUMMARY_PATH, sep="\t", index=False)
    write_dual(
        origin_confidence_tiers,
        ORIGIN_CONFIDENCE_TIER_PATH,
        ORIGIN_CONFIDENCE_TIER_SUPP_PATH,
    )
    antigen_global_combined = antigen_global.merge(
        antigen_summary.drop(columns=["signal_definition"]),
        on=["locus", "locus_label"],
        how="left",
    )
    antigen_global_combined["table_scope"] = "global_overlap_frame"
    antigen_epoch_supp = antigen_epoch.copy()
    antigen_epoch_supp["table_scope"] = "selected_country_epoch_overlap_frame"
    antigen_supp = pd.concat([antigen_global_combined, antigen_epoch_supp], ignore_index=True, sort=False)
    write_dual(antigen_supp, ANTIGEN_COMBINED_PATH, ANTIGEN_SUPP_PATH)

    print("Built analysis upgrade sidecars:")
    for path in [
        YEAR_SENS_SUMMARY_PATH,
        YEAR_SENS_LOYO_PATH,
        TURNOVER_LONG_PATH,
        TURNOVER_SUMMARY_PATH,
        ORIGIN_BRIDGE_PATH,
        ANTIGEN_COMBINED_PATH,
        MISSINGNESS_INTERPRETABILITY_PATH,
        QUALITY_SENSITIVITY_PATH,
        MISSINGNESS_DR_SUMMARY_PATH,
        MISSINGNESS_TIPPING_GRID_PATH,
        READINESS_THRESHOLD_GRID_PATH,
        ASR_SCENARIO_REGISTRY_PATH,
        ASR_ONE_GLOBAL_CLONE_SUMMARY_PATH,
        ORIGIN_CONFIDENCE_TIER_PATH,
        SELECTION_SCORECARD_PATH,
        EPOCH_CONTRAST_PATH,
    ]:
        print(f" - {path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
