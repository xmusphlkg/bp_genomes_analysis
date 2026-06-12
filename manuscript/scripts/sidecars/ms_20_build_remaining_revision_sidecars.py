#!/usr/bin/env python3
"""Build sidecars that close the remaining revision-plan audit gaps.

The earlier upgrade scripts already create the main country-selection,
missingness, ASR, structural-concentration, and negative-control tables. This
script adds the remaining reviewer-facing ledgers that are mostly
cross-sectional summaries over those frozen outputs:

1. Read-linked subset transportability across country, time, quality and
   mechanism strata.
2. Tree-subset inclusion/IPW diagnostics for representativeness-aware ASR.
3. Local origin-package support profiles, including patristic distance
   summaries where local trees are available.
4. Within-origin structural concentration.
5. Event-definition hierarchy sensitivity for the structural-concentration
   claim.
"""

from __future__ import annotations

import os
import math
import sys
from itertools import combinations
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from Bio import Phylo
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "workflow" / "lib"))

from project_paths import project_module_data_root


DATA_HOME = Path(
    os.environ.get(
        "PERTUSSIS_PROJECT_DATA_ROOT",
        str(REPO_ROOT / "pertussis_data" / "pertussis_gene"),
    )
)
FIGURE_DATA_DIR = REPO_ROOT / "manuscript" / "figure_data"
SELECTED_DIR = FIGURE_DATA_DIR / "selected_country"
SUPP_DIR = REPO_ROOT / "manuscript" / "supplementary"
STEP4_OUTPUTS = project_module_data_root("step4_prn_validation") / "outputs"

MANIFEST_PATH = REPO_ROOT / "outputs" / "workflow" / "manifest" / "manifest.tsv"
ASR_TIP_STATES_PATH = REPO_ROOT / "outputs" / "workflow" / "asr" / "tip_states.tsv"
ASR_SCENARIO_REGISTRY_PATH = FIGURE_DATA_DIR / "asr_scenario_registry.tsv"
ASR_REPRESENTATIVENESS_AUDIT_PATH = FIGURE_DATA_DIR / "supplementary_programme_representativeness_audit.tsv"
ASR_RESAMPLING_SUMMARY_PATH = FIGURE_DATA_DIR / "figure3_workflow_asr_resampling.tsv"
LOCAL_PACKAGE_SUMMARY_PATH = FIGURE_DATA_DIR / "local_rooted_package_tree_summary.tsv"
ORIGIN_PACKAGE_SUMMARY_PATH = SELECTED_DIR / "selected_country_origin_package_summary.tsv"
ORIGIN_PACKAGE_CONTEXT_PATH = FIGURE_DATA_DIR / "origin_package_context.tsv"
EVENT_MANIFEST_PATH = FIGURE_DATA_DIR / "prn_event_evidence_manifest.tsv"
QUALITY_SENSITIVITY_PATH = SELECTED_DIR / "selected_country_quality_restricted_sensitivity.tsv"
DETECTABILITY_SUMMARY_PATH = FIGURE_DATA_DIR / "prn_event_class_detectability.tsv"
DETECTABILITY_DETAIL_PATH = FIGURE_DATA_DIR / "prn_event_class_detectability_detail.tsv"
DETECTABILITY_STRESS_RESULTS_PATH = (
    project_module_data_root("step4_prn_validation")
    / "work"
    / "read_validation"
    / "detectability_stress"
    / "bp_prn_detectability_stress_results.tsv"
)

READ_LINKED_TRANSPORT_PATH = SELECTED_DIR / "selected_country_read_linked_transportability_ledger.tsv"
READ_LINKED_TRANSPORT_SUPP_PATH = SUPP_DIR / "Supplementary_Table_46_Read_Linked_Transportability.tsv"
ASR_REPRESENTATIVENESS_PATH = FIGURE_DATA_DIR / "asr_representativeness_adjustment_summary.tsv"
ASR_REPRESENTATIVENESS_SUPP_PATH = SUPP_DIR / "Supplementary_Table_51_ASR_Representativeness_Adjustment.tsv"
LOCAL_SUPPORT_PATH = FIGURE_DATA_DIR / "local_origin_package_support_profile.tsv"
LOCAL_SUPPORT_SUPP_PATH = SUPP_DIR / "Supplementary_Table_48_Local_Origin_Package_Support_Profile.tsv"
WITHIN_ORIGIN_PATH = FIGURE_DATA_DIR / "within_origin_structural_concentration.tsv"
WITHIN_ORIGIN_SUPP_PATH = SUPP_DIR / "Supplementary_Table_49_Within_Origin_Structural_Concentration.tsv"
EVENT_HIERARCHY_PATH = FIGURE_DATA_DIR / "event_definition_hierarchy_sensitivity.tsv"
EVENT_HIERARCHY_SUPP_PATH = SUPP_DIR / "Supplementary_Table_50_Event_Definition_Hierarchy_Sensitivity.tsv"
MISSINGNESS_VISUAL_PATH = SELECTED_DIR / "selected_country_missingness_visual_summary.tsv"

SELECTED_COUNTRIES = ("USA", "NZL", "AUS", "JPN")
TOP_COUNTRY_MIN_INTERPRETABLE = 10


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def rel_source(path: Path) -> str:
    for base in (REPO_ROOT, DATA_HOME):
        try:
            return str(path.relative_to(base))
        except ValueError:
            continue
    return str(path)


def write_dual(df: pd.DataFrame, main_path: Path, supp_path: Path) -> None:
    write_tsv(df, main_path)
    write_tsv(df, supp_path)


def clean_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = str(value).strip()
    if text.casefold() in {"", "nan", "none", "na", "not available"}:
        return ""
    return text


def to_bool(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.strip()
        .str.casefold()
        .isin({"true", "1", "yes", "y", "t"})
    )


def to_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def first_nonempty(values: Iterable[object]) -> str:
    for value in values:
        text = clean_text(value)
        if text:
            return text
    return ""


def mechanism_group(call: object, disrupted: object | None = None) -> str:
    text = clean_text(call).casefold()
    if disrupted is not None and not bool(disrupted):
        return "intact"
    if not text or text == "intact":
        return "intact"
    if "is481" in text:
        return "IS481"
    if "inversion" in text or "rearrangement" in text:
        return "inversion_or_rearrangement"
    if "insufficient" in text:
        return "insufficient"
    return "other_disruption"


def year_band(value: object) -> str:
    year = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(year):
        return "unknown_year"
    year = int(year)
    if year < 1990:
        return "pre_1990"
    if year < 2000:
        return "1990_1999"
    if year < 2010:
        return "2000_2009"
    if year < 2015:
        return "2010_2014"
    if year < 2020:
        return "2015_2019"
    return "2020_plus"


def assembly_quality(row: pd.Series) -> str:
    has_reads = bool(row.get("has_reads_bool", False))
    n_contigs = row.get("n_contigs_num", np.nan)
    n50 = row.get("contig_n50_num", np.nan)
    if not has_reads:
        return "assembly_only_or_no_raw_reads"
    if pd.notna(n_contigs) and pd.notna(n50) and n_contigs <= 100 and n50 >= 50_000:
        return "read_linked_compact"
    if pd.notna(n_contigs) and n_contigs <= 100:
        return "read_linked_low_fragmentation"
    if pd.notna(n50) and n50 >= 50_000:
        return "read_linked_high_n50"
    return "read_linked_fragmented_or_low_n50"


def prepare_manifest() -> pd.DataFrame:
    manifest = read_tsv(MANIFEST_PATH)
    manifest["prn_interpretable_bool"] = to_bool(manifest["prn_interpretable"])
    manifest["prn_disrupted_bool"] = to_bool(manifest["prn_disrupted"])
    manifest["has_reads_bool"] = to_bool(manifest["has_reads"])
    manifest["year_num"] = to_numeric(manifest["year"])
    manifest["n_contigs_num"] = to_numeric(manifest["n_contigs"])
    manifest["contig_n50_num"] = to_numeric(manifest["contig_n50"])
    manifest["total_sequence_length_num"] = to_numeric(manifest["total_sequence_length"])
    manifest["year_band"] = manifest["year"].map(year_band)
    manifest["mechanism_group"] = manifest.apply(
        lambda row: mechanism_group(row.get("prn_mechanism_call", ""), row.get("prn_disrupted_bool", False)),
        axis=1,
    )
    manifest["assembly_quality_bucket"] = manifest.apply(assembly_quality, axis=1)
    manifest["country_iso3_clean"] = manifest["country_iso3"].replace("", "unknown")
    return manifest


def summarise_group(frame: pd.DataFrame) -> dict[str, object]:
    if frame.empty:
        return {
            "n_interpretable": 0,
            "fraction_of_subset": 0.0,
            "n_disrupted": 0,
            "disrupted_fraction": np.nan,
            "median_n_contigs": np.nan,
            "median_contig_n50": np.nan,
            "median_year": np.nan,
        }
    return {
        "n_interpretable": int(len(frame)),
        "fraction_of_subset": np.nan,
        "n_disrupted": int(frame["prn_disrupted_bool"].sum()),
        "disrupted_fraction": float(frame["prn_disrupted_bool"].mean()),
        "median_n_contigs": float(frame["n_contigs_num"].median(skipna=True)),
        "median_contig_n50": float(frame["contig_n50_num"].median(skipna=True)),
        "median_year": float(frame["year_num"].median(skipna=True)),
    }


def build_transport_dimension(
    manifest: pd.DataFrame,
    dimension: str,
    column: str,
    categories: Iterable[str] | None = None,
) -> list[dict[str, object]]:
    full = manifest.loc[manifest["prn_interpretable_bool"]].copy()
    read_linked = full.loc[full["has_reads_bool"]].copy()
    if categories is None:
        categories = sorted(set(full[column].dropna().astype(str)) | set(read_linked[column].dropna().astype(str)))
    rows: list[dict[str, object]] = []
    total_full = len(full)
    total_read = len(read_linked)
    for category in categories:
        full_group = full.loc[full[column].astype(str) == str(category)]
        read_group = read_linked.loc[read_linked[column].astype(str) == str(category)]
        if dimension == "country" and category not in SELECTED_COUNTRIES and len(full_group) < TOP_COUNTRY_MIN_INTERPRETABLE:
            continue
        full_summary = summarise_group(full_group)
        read_summary = summarise_group(read_group)
        full_fraction = len(full_group) / total_full if total_full else np.nan
        read_fraction = len(read_group) / total_read if total_read else np.nan
        rows.append(
            {
                "comparison_dimension": dimension,
                "category": category,
                "full_interpretable_count": full_summary["n_interpretable"],
                "read_linked_interpretable_count": read_summary["n_interpretable"],
                "full_interpretable_fraction": full_fraction,
                "read_linked_interpretable_fraction": read_fraction,
                "fraction_gap_read_minus_full": read_fraction - full_fraction,
                "read_to_full_fraction_ratio": (
                    read_fraction / full_fraction if full_fraction and not pd.isna(full_fraction) else np.nan
                ),
                "full_disrupted_count": full_summary["n_disrupted"],
                "read_linked_disrupted_count": read_summary["n_disrupted"],
                "full_disrupted_fraction": full_summary["disrupted_fraction"],
                "read_linked_disrupted_fraction": read_summary["disrupted_fraction"],
                "disrupted_fraction_gap_read_minus_full": (
                    read_summary["disrupted_fraction"] - full_summary["disrupted_fraction"]
                    if not pd.isna(read_summary["disrupted_fraction"]) and not pd.isna(full_summary["disrupted_fraction"])
                    else np.nan
                ),
                "full_median_n_contigs": full_summary["median_n_contigs"],
                "read_linked_median_n_contigs": read_summary["median_n_contigs"],
                "full_median_contig_n50": full_summary["median_contig_n50"],
                "read_linked_median_contig_n50": read_summary["median_contig_n50"],
                "full_median_year": full_summary["median_year"],
                "read_linked_median_year": read_summary["median_year"],
                "notes": "Read-linked subset compared with all PRN-interpretable genomes in the retained manifest.",
            }
        )
    return rows


def build_read_linked_transportability(manifest: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    rows.extend(build_transport_dimension(manifest, "country", "country_iso3_clean"))
    rows.extend(build_transport_dimension(manifest, "year_band", "year_band"))
    rows.extend(build_transport_dimension(manifest, "mechanism_group", "mechanism_group"))
    rows.extend(build_transport_dimension(manifest, "assembly_quality", "assembly_quality_bucket"))

    if QUALITY_SENSITIVITY_PATH.exists():
        quality = read_tsv(QUALITY_SENSITIVITY_PATH)
        has_reads_rows = quality.loc[quality["subset_id"] == "has_reads"].copy()
        for row in has_reads_rows.itertuples(index=False):
            rows.append(
                {
                    "comparison_dimension": "selected_country_contrast_direction",
                    "category": f"{row.country_iso3}:{row.comparison_pre_epoch_id}->{row.comparison_post_epoch_id}",
                    "full_interpretable_count": "",
                    "read_linked_interpretable_count": "",
                    "full_interpretable_fraction": "",
                    "read_linked_interpretable_fraction": "",
                    "fraction_gap_read_minus_full": "",
                    "read_to_full_fraction_ratio": "",
                    "full_disrupted_count": "",
                    "read_linked_disrupted_count": "",
                    "full_disrupted_fraction": row.baseline_naive_delta_all_interpretable,
                    "read_linked_disrupted_fraction": row.subset_naive_delta,
                    "disrupted_fraction_gap_read_minus_full": "",
                    "full_median_n_contigs": "",
                    "read_linked_median_n_contigs": "",
                    "full_median_contig_n50": "",
                    "read_linked_median_contig_n50": "",
                    "full_median_year": "",
                    "read_linked_median_year": "",
                    "notes": (
                        "Contrast-level row: full/read-linked fields encode all-interpretable and "
                        "read-linked naive epoch-delta estimates; "
                        f"direction_matches_all_interpretable={row.direction_matches_all_interpretable}; "
                        f"subset_estimable={row.subset_estimable}."
                    ),
                }
            )
    return pd.DataFrame(rows)


def build_missingness_visual_summary() -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    if (SELECTED_DIR / "prn_interpretability_model.tsv").exists():
        interp = read_tsv(SELECTED_DIR / "prn_interpretability_model.tsv")
        calibration = interp.loc[interp["row_type"] == "calibration_bin"].copy()
        if not calibration.empty:
            calibration["display_scope"] = "calibration_curve"
            rows.append(calibration)
    if (SELECTED_DIR / "selected_country_dr_missingness_summary.tsv").exists():
        dr = read_tsv(SELECTED_DIR / "selected_country_dr_missingness_summary.tsv")
        melted = dr.melt(
            id_vars=["country_iso3", "country_name", "previous_epoch_id", "next_epoch_id"],
            value_vars=[
                col
                for col in dr.columns
                if col.startswith("delta_") or col.startswith("weight_") or col.startswith("propensity_")
            ],
            var_name="metric",
            value_name="value",
        )
        melted["row_type"] = "dr_weight_or_delta_metric"
        melted["display_scope"] = "dr_missingness_summary"
        rows.append(melted)
    tipping_grid_path = SELECTED_DIR / "selected_country_missingness_tipping_grid.tsv"
    if tipping_grid_path.exists():
        tipping = read_tsv(tipping_grid_path)
        odds = tipping.loc[tipping["sensitivity_model"] == "odds_multiplier_pattern_mixture"].copy()
        odds["row_type"] = "odds_multiplier_tipping_grid"
        odds["display_scope"] = "tipping_heatmap"
        rows.append(odds)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True, sort=False)


def fit_tree_inclusion_model(manifest: pd.DataFrame) -> pd.DataFrame:
    if not ASR_TIP_STATES_PATH.exists():
        return pd.DataFrame()
    tip_states = read_tsv(ASR_TIP_STATES_PATH)
    tree_assemblies = set(tip_states["assembly_accession"].dropna().astype(str))
    frame = manifest.loc[manifest["prn_interpretable_bool"]].copy()
    frame["in_primary_asr_tree"] = frame["assembly_accession"].isin(tree_assemblies)
    if frame["in_primary_asr_tree"].nunique() < 2:
        return pd.DataFrame()

    top_countries = frame["country_iso3_clean"].value_counts().head(12).index
    frame["country_model"] = np.where(frame["country_iso3_clean"].isin(top_countries), frame["country_iso3_clean"], "other")
    features = [
        "year_num",
        "has_reads_bool",
        "n_contigs_num",
        "contig_n50_num",
        "country_model",
        "year_band",
        "mechanism_group",
    ]
    numeric_features = ["year_num", "n_contigs_num", "contig_n50_num"]
    categorical_features = ["has_reads_bool", "country_model", "year_band", "mechanism_group"]
    preprocessor = ColumnTransformer(
        transformers=[
            ("num", Pipeline([("imputer", SimpleImputer(strategy="median")), ("scale", StandardScaler())]), numeric_features),
            (
                "cat",
                Pipeline(
                    [
                        ("imputer", SimpleImputer(strategy="most_frequent")),
                        ("onehot", OneHotEncoder(handle_unknown="ignore")),
                    ]
                ),
                categorical_features,
            ),
        ]
    )
    model = Pipeline(
        [
            ("preprocessor", preprocessor),
            (
                "classifier",
                LogisticRegression(max_iter=2000, class_weight="balanced", solver="lbfgs", random_state=20260418),
            ),
        ]
    )
    y = frame["in_primary_asr_tree"].astype(int)
    model.fit(frame[features], y)
    predicted = model.predict_proba(frame[features])[:, 1]
    predicted = np.clip(predicted, 0.01, 0.99)
    frame["tree_inclusion_probability"] = predicted
    frame["tree_inclusion_ipw"] = np.where(frame["in_primary_asr_tree"], 1.0 / predicted, np.nan)
    auc = roc_auc_score(y, predicted)

    tree_frame = frame.loc[frame["in_primary_asr_tree"]].copy()
    weight = tree_frame["tree_inclusion_ipw"]
    weight_sum = weight.sum()
    weighted_disrupted = float((tree_frame["prn_disrupted_bool"].astype(float) * weight).sum() / weight_sum)

    rows: list[dict[str, object]] = [
        {
            "row_type": "tree_inclusion_model",
            "comparison_dimension": "overall",
            "category": "primary_asr_tree",
            "full_interpretable_count": int(len(frame)),
            "tree_subset_count": int(tree_frame.shape[0]),
            "tree_subset_fraction": float(tree_frame.shape[0] / len(frame)),
            "model_auc": float(auc),
            "weight_mean": float(weight.mean()),
            "weight_median": float(weight.median()),
            "weight_q95": float(weight.quantile(0.95)),
            "weight_max": float(weight.max()),
            "unweighted_tree_disrupted_fraction": float(tree_frame["prn_disrupted_bool"].mean()),
            "ipw_tree_disrupted_fraction": weighted_disrupted,
            "notes": "Inverse-probability-of-tree-inclusion diagnostic; weights are descriptive and do not re-run ASR transitions.",
        }
    ]

    for dimension, column in [
        ("country", "country_iso3_clean"),
        ("year_band", "year_band"),
        ("mechanism_group", "mechanism_group"),
    ]:
        for category, subset in tree_frame.groupby(column, dropna=False):
            if subset.empty:
                continue
            rows.append(
                {
                    "row_type": "tree_inclusion_weighted_composition",
                    "comparison_dimension": dimension,
                    "category": category,
                    "full_interpretable_count": int((frame[column] == category).sum()),
                    "tree_subset_count": int(len(subset)),
                    "tree_subset_fraction": float(len(subset) / len(tree_frame)),
                    "model_auc": float(auc),
                    "weight_mean": float(subset["tree_inclusion_ipw"].mean()),
                    "weight_median": float(subset["tree_inclusion_ipw"].median()),
                    "weight_q95": float(subset["tree_inclusion_ipw"].quantile(0.95)),
                    "weight_max": float(subset["tree_inclusion_ipw"].max()),
                    "unweighted_tree_disrupted_fraction": float(subset["prn_disrupted_bool"].mean()),
                    "ipw_tree_disrupted_fraction": float(
                        (
                            subset["prn_disrupted_bool"].astype(float) * subset["tree_inclusion_ipw"]
                        ).sum()
                        / subset["tree_inclusion_ipw"].sum()
                    ),
                    "notes": "Composition stratum within the primary ASR tree subset.",
                }
            )
    return pd.DataFrame(rows)


def build_asr_representativeness_summary(manifest: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if ASR_REPRESENTATIVENESS_AUDIT_PATH.exists():
        audit = read_tsv(ASR_REPRESENTATIVENESS_AUDIT_PATH)
        for dimension, group in audit.groupby("comparison_dimension", dropna=False):
            absolute_gap = pd.to_numeric(group["absolute_fraction_gap"], errors="coerce")
            max_idx = absolute_gap.idxmax()
            max_row = group.loc[max_idx]
            rows.append(
                {
                    "row_type": "representativeness_gap_summary",
                    "scenario_or_dimension": dimension,
                    "n_categories": int(group.shape[0]),
                    "max_absolute_fraction_gap": float(absolute_gap.max(skipna=True)),
                    "mean_absolute_fraction_gap": float(absolute_gap.mean(skipna=True)),
                    "largest_gap_category": max_row.get("category", ""),
                    "full_interpretable_fraction_largest_gap": max_row.get("full_interpretable_fraction", ""),
                    "tree_subset_fraction_largest_gap": max_row.get("tree_subset_fraction", ""),
                    "n_scenarios": "",
                    "median_fitch_origin_events": "",
                    "min_fitch_origin_events": "",
                    "max_fitch_origin_events": "",
                    "proportion_rejecting_one_global_clone": "",
                    "notes": "Tree-subset representativeness audit against the full PRN-interpretable cohort.",
                }
            )
    if ASR_RESAMPLING_SUMMARY_PATH.exists():
        resampling = read_tsv(ASR_RESAMPLING_SUMMARY_PATH)
        for row in resampling.itertuples(index=False):
            rows.append(
                {
                    "row_type": "composition_balanced_resampling_summary",
                    "scenario_or_dimension": row.scheme,
                    "n_categories": "",
                    "max_absolute_fraction_gap": "",
                    "mean_absolute_fraction_gap": "",
                    "largest_gap_category": "",
                    "full_interpretable_fraction_largest_gap": "",
                    "tree_subset_fraction_largest_gap": "",
                    "n_scenarios": row.n_replicates,
                    "median_fitch_origin_events": row.fitch_origin_events_median,
                    "min_fitch_origin_events": row.fitch_origin_events_min,
                    "max_fitch_origin_events": row.fitch_origin_events_max,
                    "proportion_rejecting_one_global_clone": "",
                    "notes": "Balanced ASR resampling summary used as the representativeness-aware robustness layer.",
                }
            )
    if ASR_SCENARIO_REGISTRY_PATH.exists():
        registry = read_tsv(ASR_SCENARIO_REGISTRY_PATH)
        for scenario_class, group in registry.groupby("scenario_class", dropna=False):
            fitch = pd.to_numeric(group["fitch_origin_events"], errors="coerce")
            rejects = to_bool(group["rejects_one_global_clone_fitch"])
            rows.append(
                {
                    "row_type": "registered_asr_scenario_distribution",
                    "scenario_or_dimension": scenario_class,
                    "n_categories": "",
                    "max_absolute_fraction_gap": "",
                    "mean_absolute_fraction_gap": "",
                    "largest_gap_category": "",
                    "full_interpretable_fraction_largest_gap": "",
                    "tree_subset_fraction_largest_gap": "",
                    "n_scenarios": int(group["scenario_id"].nunique()),
                    "median_fitch_origin_events": float(fitch.median(skipna=True)) if fitch.notna().any() else "",
                    "min_fitch_origin_events": float(fitch.min(skipna=True)) if fitch.notna().any() else "",
                    "max_fitch_origin_events": float(fitch.max(skipna=True)) if fitch.notna().any() else "",
                    "proportion_rejecting_one_global_clone": float(rejects.mean()) if len(rejects) else "",
                    "notes": "Registered ASR scenario distribution before stochastic-mapping augmentation.",
                }
            )
    ipw = fit_tree_inclusion_model(manifest)
    if not ipw.empty:
        rows.extend(ipw.to_dict(orient="records"))
    return pd.DataFrame(rows)


def pairwise_patristic_summary(tree_path: Path, tip_labels: list[str]) -> dict[str, object]:
    if not tree_path.exists() or len(tip_labels) < 2:
        return {
            "pairwise_distance_pairs": 0,
            "pairwise_distance_min": np.nan,
            "pairwise_distance_median": np.nan,
            "pairwise_distance_max": np.nan,
        }
    tree = Phylo.read(tree_path, "newick")
    available = {terminal.name for terminal in tree.get_terminals()}
    labels = [label for label in tip_labels if label in available]
    distances = []
    for left, right in combinations(labels, 2):
        try:
            distances.append(float(tree.distance(left, right)))
        except Exception:
            continue
    if not distances:
        return {
            "pairwise_distance_pairs": 0,
            "pairwise_distance_min": np.nan,
            "pairwise_distance_median": np.nan,
            "pairwise_distance_max": np.nan,
        }
    arr = np.array(distances, dtype=float)
    return {
        "pairwise_distance_pairs": int(len(arr)),
        "pairwise_distance_min": float(np.min(arr)),
        "pairwise_distance_median": float(np.median(arr)),
        "pairwise_distance_max": float(np.max(arr)),
    }


def build_local_origin_support_profile(manifest: pd.DataFrame) -> pd.DataFrame:
    if not LOCAL_PACKAGE_SUMMARY_PATH.exists():
        return pd.DataFrame()
    local = read_tsv(LOCAL_PACKAGE_SUMMARY_PATH)
    if ORIGIN_PACKAGE_CONTEXT_PATH.exists():
        origin_summary = read_tsv(ORIGIN_PACKAGE_CONTEXT_PATH)
    elif ORIGIN_PACKAGE_SUMMARY_PATH.exists():
        origin_summary = read_tsv(ORIGIN_PACKAGE_SUMMARY_PATH)
    else:
        origin_summary = pd.DataFrame()
    event_lookup = (
        manifest.loc[manifest["prn_disrupted_bool"], ["assembly_accession", "prn_event_id"]]
        .drop_duplicates()
        .set_index("assembly_accession")["prn_event_id"]
        .to_dict()
    )
    rows: list[dict[str, object]] = []
    for row in local.itertuples(index=False):
        origin_id = row.origin_id
        tip_state_path = REPO_ROOT / str(row.local_asr_dir) / "tip_states.tsv"
        tree_path = REPO_ROOT / str(row.local_asr_dir) / "rooted_ml_tree.reference_rooted.nwk"
        if not tip_state_path.exists():
            continue
        tips = read_tsv(tip_state_path)
        tips["is_disrupted"] = tips["prn_state"].astype(str) == "disrupted"
        tips["year_num"] = to_numeric(tips["year"])
        tips["event_id"] = tips["assembly_accession"].map(event_lookup).fillna("")
        disrupted = tips.loc[tips["is_disrupted"]].copy()
        dominant_event = clean_text(getattr(row, "dominant_prn_event_id", ""))
        dominant_event_count = int((disrupted["event_id"] == dominant_event).sum()) if dominant_event else 0
        dist_summary = pairwise_patristic_summary(tree_path, disrupted["tree_tip_label"].dropna().astype(str).tolist())
        origin_row = (
            origin_summary.loc[origin_summary["origin_id"] == origin_id].iloc[0].to_dict()
            if not origin_summary.empty and (origin_summary["origin_id"] == origin_id).any()
            else {}
        )
        rows.append(
            {
                "origin_id": origin_id,
                "dominant_prn_event_id": dominant_event,
                "package_country_iso3": getattr(row, "package_country_iso3", ""),
                "local_tip_count": getattr(row, "local_tip_count", ""),
                "local_disrupted_tip_count": int(disrupted.shape[0]),
                "local_fitch_origin_events": getattr(row, "local_fitch_origin_events", ""),
                "local_pastml_origin_events": getattr(row, "local_pastml_origin_events", ""),
                "local_origin_consistency_status": getattr(row, "local_origin_consistency_status", ""),
                "preserves_single_origin_consistent_package": getattr(row, "preserves_single_origin_consistent_package", ""),
                "representative_covering_branch_support": getattr(row, "representative_covering_branch_support", ""),
                "local_country_count": int(disrupted["country_iso3"].replace("", np.nan).nunique(dropna=True)),
                "local_year_min": int(disrupted["year_num"].min()) if disrupted["year_num"].notna().any() else "",
                "local_year_max": int(disrupted["year_num"].max()) if disrupted["year_num"].notna().any() else "",
                "local_mlst_st_count": int(disrupted["mlst_st"].replace("", np.nan).nunique(dropna=True)),
                "dominant_event_disrupted_tip_count": dominant_event_count,
                "dominant_event_disrupted_tip_share": (
                    dominant_event_count / disrupted.shape[0] if disrupted.shape[0] else np.nan
                ),
                "origin_table_n_disrupted_descendants": origin_row.get("n_disrupted_descendants", ""),
                "origin_table_branch_support": origin_row.get("branch_support", ""),
                "origin_table_confidence": origin_row.get("origin_confidence", ""),
                **dist_summary,
                "notes": "Pairwise distances are patristic distances from the local rooted SNP tree among disrupted local tips.",
            }
        )
    return pd.DataFrame(rows)


def build_within_origin_structural_concentration(manifest: pd.DataFrame) -> pd.DataFrame:
    subtree_dir = REPO_ROOT / "outputs" / "workflow" / "asr" / "event_subtrees"
    if not subtree_dir.exists():
        return pd.DataFrame()
    event_lookup = (
        manifest.loc[:, ["assembly_accession", "prn_event_id", "prn_disrupted_bool", "country_iso3_clean", "year", "mlst_st"]]
        .drop_duplicates(subset=["assembly_accession"])
        .set_index("assembly_accession")
        .to_dict(orient="index")
    )
    rows: list[dict[str, object]] = []
    for path in sorted(subtree_dir.glob("origin_*.descendant_tips.tsv")):
        origin_id = path.name.replace(".descendant_tips.tsv", "")
        tips = read_tsv(path)
        records = []
        for tip in tips.itertuples(index=False):
            accession = clean_text(getattr(tip, "assembly_accession", ""))
            metadata = event_lookup.get(accession, {})
            is_disrupted = bool(metadata.get("prn_disrupted_bool", False)) or clean_text(getattr(tip, "observed_prn_state", "")) == "disrupted"
            if is_disrupted:
                records.append(
                    {
                        "assembly_accession": accession,
                        "prn_event_id": clean_text(metadata.get("prn_event_id", "")),
                        "country_iso3": clean_text(metadata.get("country_iso3_clean", "")),
                        "year": pd.to_numeric(pd.Series([metadata.get("year", "")]), errors="coerce").iloc[0],
                        "mlst_st": clean_text(metadata.get("mlst_st", "")),
                    }
                )
        if not records:
            continue
        df = pd.DataFrame(records)
        counts = df["prn_event_id"].replace("", "unknown_event").value_counts()
        top3 = counts.head(3).sum()
        dominant_event = counts.index[0]
        rows.append(
            {
                "origin_id": origin_id,
                "n_disrupted_descendant_tips": int(df.shape[0]),
                "n_unique_events": int(counts.shape[0]),
                "dominant_prn_event_id": dominant_event,
                "dominant_event_count": int(counts.iloc[0]),
                "dominant_event_share": float(counts.iloc[0] / df.shape[0]),
                "top3_event_share": float(top3 / df.shape[0]),
                "country_count": int(df["country_iso3"].replace("", np.nan).nunique(dropna=True)),
                "mlst_st_count": int(df["mlst_st"].replace("", np.nan).nunique(dropna=True)),
                "year_min": int(df["year"].min()) if df["year"].notna().any() else "",
                "year_max": int(df["year"].max()) if df["year"].notna().any() else "",
                "notes": "Computed over disrupted descendant tips in the primary ASR event-subtree package.",
            }
        )
    return pd.DataFrame(rows)


def event_hierarchy_label(event_row: pd.Series, tier: str) -> str:
    mechanism = clean_text(event_row.get("mechanism_call", ""))
    subcat = clean_text(event_row.get("event_subcategory", ""))
    gap = pd.to_numeric(pd.Series([event_row.get("insertion_subject_gap_bp", "")]), errors="coerce").iloc[0]
    orientation = clean_text(event_row.get("orientation", "")) or clean_text(event_row.get("hit_orientation", ""))
    basis = clean_text(event_row.get("breakpoint_coordinate_basis", ""))
    tsd_status = clean_text(event_row.get("tsd_or_flank_sequence_status", ""))
    event_id = clean_text(event_row.get("prn_event_id", ""))

    if tier == "exact_breakpoint_orientation_tsd":
        if basis == "read_reference" and orientation and tsd_status == "target_site_duplication_recovered":
            return f"read_exact_tsd|{event_id}"
        return f"not_exact_read_tsd|{mechanism}|{subcat}|{event_id}"
    if tier == "breakpoint_window_5bp":
        if pd.notna(gap):
            lower = int(math.floor((float(gap) - 1038) / 11) * 11 + 1038)
            upper = lower + 10
            if 1038 <= float(gap) <= 1048:
                lower, upper = 1038, 1048
            return f"{mechanism}|gap_{lower}_{upper}|orientation_{orientation or 'unknown'}"
        return f"{mechanism}|no_gap_window|orientation_{orientation or 'unknown'}"
    if tier == "mechanism_approximate_architecture":
        if pd.notna(gap):
            if 1038 <= float(gap) <= 1048:
                gap_group = "gap_1043_plusminus5"
            else:
                gap_group = f"gap_{int(gap)}"
        else:
            gap_group = clean_text(event_row.get("bp_category", "")) or "no_bp_category"
        return f"{mechanism}|{subcat}|{gap_group}"
    if tier == "broad_mechanism_class":
        return mechanism_group(mechanism, True)
    raise ValueError(f"Unknown tier: {tier}")


def concentration_metrics(counts: pd.Series) -> dict[str, object]:
    total = counts.sum()
    if total <= 0:
        return {}
    shares = counts / total
    return {
        "n_genomes": int(total),
        "n_event_definitions": int(counts.shape[0]),
        "dominant_definition": counts.index[0],
        "dominant_definition_count": int(counts.iloc[0]),
        "dominant_definition_share": float(counts.iloc[0] / total),
        "top3_definition_share": float(counts.head(3).sum() / total),
        "hhi": float((shares**2).sum()),
        "effective_number": float(1.0 / (shares**2).sum()) if (shares**2).sum() else np.nan,
    }


def build_event_definition_hierarchy(manifest: pd.DataFrame) -> pd.DataFrame:
    if not EVENT_MANIFEST_PATH.exists():
        return pd.DataFrame()
    event_manifest = read_tsv(EVENT_MANIFEST_PATH)
    sample_counts = (
        manifest.loc[manifest["prn_disrupted_bool"], "prn_event_id"]
        .replace("", np.nan)
        .dropna()
        .value_counts()
        .rename("manifest_disrupted_count")
        .reset_index()
        .rename(columns={"index": "prn_event_id"})
    )
    events = event_manifest.merge(sample_counts, on="prn_event_id", how="left")
    events["manifest_disrupted_count"] = pd.to_numeric(events["manifest_disrupted_count"], errors="coerce").fillna(
        pd.to_numeric(events.get("sample_count", 0), errors="coerce").fillna(0)
    )
    tiers = [
        ("exact_breakpoint_orientation_tsd", "Exact breakpoint + orientation + recovered TSD where read-reference evidence exists"),
        ("breakpoint_window_5bp", "Breakpoint/gap grouped in +/-5 bp windows with orientation"),
        ("mechanism_approximate_architecture", "Mechanism class plus approximate gap/architecture"),
        ("broad_mechanism_class", "Broad disruption mechanism only"),
    ]
    rows: list[dict[str, object]] = []
    for tier_id, tier_label in tiers:
        labels = events.apply(lambda row: event_hierarchy_label(row, tier_id), axis=1)
        counts = events.groupby(labels)["manifest_disrupted_count"].sum().sort_values(ascending=False)
        metrics = concentration_metrics(counts)
        if metrics:
            rows.append(
                {
                    "event_definition_tier": tier_id,
                    "event_definition_label": tier_label,
                    **metrics,
                    "notes": "Counts use disrupted genomes from the retained manifest mapped through the event evidence manifest.",
                }
            )
    return pd.DataFrame(rows)


def detectability_family_info(prn_event_id: str) -> dict[str, object]:
    event_id = clean_text(prn_event_id)
    if event_id == "prn_evt_coding_disrupted_is481__is481__gap1043":
        return {
            "family_key": "is481_1043",
            "family_label": "IS481-associated 1,043-bp architecture",
            "plot_group": "primary",
            "display_order": 1,
        }
    if event_id.startswith("prn_evt_rearrangement__within_contig__"):
        return {
            "family_key": "rearrangement_family",
            "family_label": "Rearrangement family",
            "plot_group": "primary",
            "display_order": 2,
        }
    if event_id.startswith("prn_evt_other_disruption__insertion_like__"):
        return {
            "family_key": "other_insertion_like",
            "family_label": "Other insertion-like disruptions",
            "plot_group": "primary",
            "display_order": 3,
        }
    if event_id == "prn_evt_coding_disrupted_is481__is481__gap1042":
        return {
            "family_key": "minor_is481",
            "family_label": "Minor IS481 architecture",
            "plot_group": "context",
            "display_order": 4,
        }
    if event_id == "prn_evt_intact":
        return {
            "family_key": "intact_control",
            "family_label": "Intact control",
            "plot_group": "context",
            "display_order": 5,
        }
    if event_id.startswith("prn_evt_insufficient__"):
        return {
            "family_key": "insufficient_data",
            "family_label": "Insufficient data",
            "plot_group": "context",
            "display_order": 6,
        }
    return {
        "family_key": "other_uncategorized",
        "family_label": event_id or "Uncategorized",
        "plot_group": "context",
        "display_order": 99,
    }


def detectability_status_bucket(read_validation_status: str) -> str:
    status = clean_text(read_validation_status)
    if status in {"supported", "supported_candidate", "supported_concordant"}:
        return "recovered"
    if status == "no_prn_is_signal_detected":
        return "true_nonrecovery"
    if status in {"tool_output_missing", "unresolved"}:
        return "compatibility_excluded"
    return "compatibility_excluded"


def detectability_compatibility_state(row: pd.Series) -> str:
    status = clean_text(row.get("read_validation_status", ""))
    plan_status = clean_text(row.get("recovery_plan_status", ""))
    if status in {"supported", "supported_candidate", "supported_concordant"}:
        return "observed_recovered"
    if status == "no_prn_is_signal_detected":
        return "observed_no_prn_signal"
    if plan_status:
        return plan_status
    if status == "tool_output_missing":
        return "tool_output_missing"
    if status == "unresolved":
        return "unresolved_no_followup_annotation"
    return "unresolved_other"


def confidence_interval(successes: int, total: int, z: float = 1.96) -> tuple[float, float]:
    if total <= 0:
        return (np.nan, np.nan)
    phat = successes / total
    denom = 1 + (z**2 / total)
    center = (phat + (z**2 / (2 * total))) / denom
    margin = z * math.sqrt((phat * (1 - phat) / total) + (z**2 / (4 * total**2))) / denom
    return (max(0.0, center - margin), min(1.0, center + margin))


def build_detectability_detail_rows(manifest: pd.DataFrame) -> pd.DataFrame:
    subset_path = STEP4_OUTPUTS / "bp_prn_validation_subset.tsv"
    read_validation_path = STEP4_OUTPUTS / "bp_prn_read_validation.tsv"
    followup_path = STEP4_OUTPUTS / "bp_prn_targeted_validation_followup_queue.tsv"

    if not subset_path.exists() or not read_validation_path.exists():
        return pd.DataFrame()

    subset = read_tsv(subset_path)
    read_validation = read_tsv(read_validation_path)
    followup = read_tsv(followup_path) if followup_path.exists() else pd.DataFrame()

    subset["analysis_layer"] = "empirical"
    if "notes" in subset.columns:
        subset = subset.rename(columns={"notes": "source_notes"})
    read_validation = read_validation.loc[
        :, ["sample_id_canonical", "prn_event_id", "read_validation_status", "read_support_class", "n_supporting_reads",
            "n_contradicting_reads", "targeted_locus_assembly_status", "validation_method", "validator_version", "notes"]
    ].copy()
    read_validation = read_validation.rename(columns={"notes": "validation_notes"})
    read_validation["sample_id_canonical"] = read_validation["sample_id_canonical"].map(clean_text)
    read_validation["prn_event_id"] = read_validation["prn_event_id"].map(clean_text)

    detail = subset.merge(
        read_validation,
        on=["sample_id_canonical", "prn_event_id"],
        how="left",
        suffixes=("", "_read"),
    )
    if not followup.empty:
        followup_subset = followup.loc[
            :,
            [
                "sample_id_canonical",
                "prn_event_id",
                "recovery_plan_status",
                "recovery_run_compatibility",
                "recovery_library_layout",
                "recovery_instrument_platform",
                "recovery_download_strategy",
                "followup_class",
                "recommended_action",
            ],
        ].copy()
        followup_subset["sample_id_canonical"] = followup_subset["sample_id_canonical"].map(clean_text)
        followup_subset["prn_event_id"] = followup_subset["prn_event_id"].map(clean_text)
        detail = detail.merge(followup_subset, on=["sample_id_canonical", "prn_event_id"], how="left")
    else:
        for col in [
            "recovery_plan_status",
            "recovery_run_compatibility",
            "recovery_library_layout",
            "recovery_instrument_platform",
            "recovery_download_strategy",
            "followup_class",
            "recommended_action",
        ]:
            detail[col] = ""

    detail["family_info"] = detail["prn_event_id"].map(detectability_family_info)
    detail["family_key"] = detail["family_info"].map(lambda x: clean_text(x.get("family_key", "")) if isinstance(x, dict) else "")
    detail["family_label"] = detail["family_info"].map(lambda x: clean_text(x.get("family_label", "")) if isinstance(x, dict) else "")
    detail["plot_group"] = detail["family_info"].map(lambda x: clean_text(x.get("plot_group", "")) if isinstance(x, dict) else "")
    detail["display_order"] = detail["family_info"].map(lambda x: int(x.get("display_order", 99)) if isinstance(x, dict) else 99)
    detail = detail.drop(columns=["family_info"])

    detail["read_validation_status"] = detail["read_validation_status"].map(clean_text)
    detail["read_support_class"] = detail["read_support_class"].map(clean_text)
    detail["targeted_locus_assembly_status"] = detail["targeted_locus_assembly_status"].map(clean_text)
    detail["recovery_plan_status"] = detail["recovery_plan_status"].map(clean_text)
    detail["recovery_run_compatibility"] = detail["recovery_run_compatibility"].map(clean_text)
    detail["recovery_library_layout"] = detail["recovery_library_layout"].map(clean_text)
    detail["recovery_instrument_platform"] = detail["recovery_instrument_platform"].map(clean_text)
    detail["recovery_download_strategy"] = detail["recovery_download_strategy"].map(clean_text)
    detail["followup_class"] = detail["followup_class"].map(clean_text)
    detail["recommended_action"] = detail["recommended_action"].map(clean_text)

    detail["analysis_layer"] = "empirical"
    detail["run_status"] = "empirical_validation"
    detail["status_bucket"] = detail["read_validation_status"].map(detectability_status_bucket)
    detail["compatibility_state"] = detail.apply(detectability_compatibility_state, axis=1)
    detail["count_in_total_denominator"] = 1
    detail["count_in_resolved_denominator"] = detail["status_bucket"].isin({"recovered", "true_nonrecovery"}).astype(int)
    detail["count_as_recovered"] = (detail["status_bucket"] == "recovered").astype(int)
    detail["count_as_true_nonrecovery"] = (detail["status_bucket"] == "true_nonrecovery").astype(int)
    detail["count_as_compatibility_excluded"] = (detail["status_bucket"] == "compatibility_excluded").astype(int)
    detail["downsample_fraction"] = ""
    detail["downsample_fraction_label"] = ""
    detail["downsample_replicate"] = ""
    detail["downsample_seed"] = ""
    detail["downsample_parent_sample_id"] = ""
    detail["downsample_n_read_pairs_total"] = ""
    detail["downsample_n_read_pairs_retained"] = ""
    detail["downsample_n_bam_records_retained"] = ""
    detail["downsample_reads_1_path"] = ""
    detail["downsample_reads_2_path"] = ""
    detail["downsample_bam_path"] = ""
    detail["parent_sample_id"] = ""
    detail["source_file"] = ";".join(
        [
            rel_source(STEP4_OUTPUTS / "bp_prn_validation_subset.tsv"),
            rel_source(STEP4_OUTPUTS / "bp_prn_read_validation.tsv"),
            rel_source(STEP4_OUTPUTS / "bp_prn_targeted_validation_followup_queue.tsv"),
        ]
    )
    detail["notes"] = detail.apply(
        lambda row: (
            f"legacy_status={clean_text(row.get('read_validation_status', ''))};"
            f"resolved_bucket={clean_text(row.get('status_bucket', ''))};"
            f"compatibility_state={clean_text(row.get('compatibility_state', ''))}"
        ),
        axis=1,
    )

    stress_path = DETECTABILITY_STRESS_RESULTS_PATH
    if stress_path.exists():
        stress = read_tsv(stress_path)
        stress["analysis_layer"] = "downsampling"
        stress["family_info"] = stress["prn_event_id"].map(detectability_family_info)
        stress["family_key"] = stress["family_info"].map(lambda x: clean_text(x.get("family_key", "")) if isinstance(x, dict) else "")
        stress["family_label"] = stress["family_info"].map(lambda x: clean_text(x.get("family_label", "")) if isinstance(x, dict) else "")
        stress["plot_group"] = stress["family_info"].map(lambda x: clean_text(x.get("plot_group", "")) if isinstance(x, dict) else "")
        stress["display_order"] = stress["family_info"].map(lambda x: int(x.get("display_order", 99)) if isinstance(x, dict) else 99)
        stress = stress.drop(columns=["family_info"], errors="ignore")
        for col in [
            "sample_id_canonical",
            "parent_sample_id",
            "prn_event_id",
            "prn_mechanism_call",
            "read_validation_status",
            "read_support_class",
            "targeted_locus_assembly_status",
            "downsample_fraction_label",
            "downsample_fraction",
            "downsample_replicate",
            "downsample_seed",
            "downsample_n_read_pairs_total",
            "downsample_n_read_pairs_retained",
            "run_status",
        ]:
            if col not in stress.columns:
                stress[col] = ""
        stress["read_validation_status"] = stress["read_validation_status"].map(clean_text)
        stress["read_support_class"] = stress["read_support_class"].map(clean_text)
        stress["targeted_locus_assembly_status"] = stress["targeted_locus_assembly_status"].map(clean_text)
        for col in ["validation_method", "validator_version", "validation_notes", "source_notes"]:
            if col in stress.columns:
                stress[col] = stress[col].map(clean_text)
            else:
                stress[col] = ""
        stress["status_bucket"] = stress["read_validation_status"].map(detectability_status_bucket)
        stress["compatibility_state"] = stress["status_bucket"].map(
            lambda bucket: {
                "recovered": "observed_recovered",
                "true_nonrecovery": "observed_no_prn_signal",
                "compatibility_excluded": "downsampled_compatibility_or_tool_failure",
            }.get(bucket, "downsampled_other")
        )
        stress["count_in_total_denominator"] = 1
        stress["count_in_resolved_denominator"] = stress["status_bucket"].isin({"recovered", "true_nonrecovery"}).astype(int)
        stress["count_as_recovered"] = (stress["status_bucket"] == "recovered").astype(int)
        stress["count_as_true_nonrecovery"] = (stress["status_bucket"] == "true_nonrecovery").astype(int)
        stress["count_as_compatibility_excluded"] = (stress["status_bucket"] == "compatibility_excluded").astype(int)
        stress["source_file"] = rel_source(stress_path)
        stress["notes"] = stress.apply(
            lambda row: (
                f"parent_sample_id={clean_text(row.get('parent_sample_id', ''))};"
                f"downsample_fraction={clean_text(row.get('downsample_fraction_label', ''))};"
                f"replicate={clean_text(row.get('downsample_replicate', ''))};"
                f"run_status={clean_text(row.get('run_status', ''))}"
            ),
            axis=1,
        )
        keep_cols = [
            "analysis_layer",
            "family_key",
            "family_label",
            "plot_group",
            "display_order",
            "sample_id_canonical",
            "parent_sample_id",
            "prn_event_id",
            "prn_mechanism_call",
            "read_validation_status",
            "read_support_class",
            "targeted_locus_assembly_status",
            "recovery_plan_status",
            "recovery_run_compatibility",
            "recovery_library_layout",
            "recovery_instrument_platform",
            "recovery_download_strategy",
            "followup_class",
            "recommended_action",
            "status_bucket",
            "compatibility_state",
            "count_in_total_denominator",
            "count_in_resolved_denominator",
            "count_as_recovered",
            "count_as_true_nonrecovery",
            "count_as_compatibility_excluded",
            "downsample_fraction_label",
            "downsample_fraction",
            "downsample_replicate",
            "downsample_seed",
            "downsample_n_read_pairs_total",
            "downsample_n_read_pairs_retained",
            "downsample_n_bam_records_retained",
            "downsample_reads_1_path",
            "downsample_reads_2_path",
            "downsample_bam_path",
            "validation_method",
            "validator_version",
            "source_notes",
            "validation_notes",
            "run_status",
            "source_file",
            "notes",
        ]
        for col in keep_cols:
            if col not in stress.columns:
                stress[col] = ""
        detail = pd.concat(
            [
                detail[keep_cols],
                stress[keep_cols],
            ],
            ignore_index=True,
            sort=False,
        )
    else:
        detail = detail[
            [
                "analysis_layer",
                "family_key",
                "family_label",
                "plot_group",
                "display_order",
                "sample_id_canonical",
                "parent_sample_id",
                "prn_event_id",
                "prn_mechanism_call",
                "read_validation_status",
                "read_support_class",
                "targeted_locus_assembly_status",
                "recovery_plan_status",
                "recovery_run_compatibility",
                "recovery_library_layout",
                "recovery_instrument_platform",
                "recovery_download_strategy",
                "followup_class",
                "recommended_action",
                "status_bucket",
                "compatibility_state",
                "count_in_total_denominator",
                "count_in_resolved_denominator",
                "count_as_recovered",
                "count_as_true_nonrecovery",
                "count_as_compatibility_excluded",
                "downsample_fraction_label",
                "downsample_fraction",
                "downsample_replicate",
                "downsample_seed",
                "downsample_n_read_pairs_total",
                "downsample_n_read_pairs_retained",
                "downsample_n_bam_records_retained",
                "downsample_reads_1_path",
                "downsample_reads_2_path",
                "downsample_bam_path",
                "validation_method",
                "validator_version",
                "source_notes",
                "validation_notes",
                "run_status",
                "source_file",
                "notes",
            ]
        ]
    return detail


def build_detectability_summary(detail: pd.DataFrame) -> pd.DataFrame:
    if detail.empty:
        return pd.DataFrame()

    empirical = detail.loc[detail["analysis_layer"] == "empirical"].copy()
    if empirical.empty:
        return pd.DataFrame()

    rows: list[dict[str, object]] = []
    grouped = empirical.groupby(["family_key", "family_label", "plot_group", "display_order"], dropna=False, sort=False)
    for (family_key, family_label, plot_group, display_order), frame in grouped:
        n_total = int(frame["count_in_total_denominator"].sum())
        n_resolved = int(frame["count_in_resolved_denominator"].sum())
        n_recovered = int(frame["count_as_recovered"].sum())
        n_true_nonrecovery = int(frame["count_as_true_nonrecovery"].sum())
        n_compatibility_excluded = int(frame["count_as_compatibility_excluded"].sum())
        n_incompatible = int((frame["recovery_plan_status"] == "linked_incompatible_run_current_short_read_validator").sum())
        n_no_download = int((frame["recovery_plan_status"] == "no_download_plan_match").sum())
        n_no_fastq = int((frame["recovery_plan_status"] == "linked_run_without_fastq_ftp").sum())
        legacy_rate = n_recovered / n_total if n_total else np.nan
        resolved_rate = n_recovered / n_resolved if n_resolved else np.nan
        ci_lower, ci_upper = confidence_interval(n_recovered, n_resolved)
        dominant_event = ""
        if "prn_event_id" in frame.columns and not frame.empty:
            dominant_event = frame["prn_event_id"].value_counts().index[0]
        rows.append(
            {
                "family_key": family_key,
                "family_label": family_label,
                "plot_group": plot_group,
                "display_order": display_order,
                "dominant_prn_event_id": dominant_event,
                "n_total": n_total,
                "n_resolved": n_resolved,
                "n_recovered": n_recovered,
                "n_true_nonrecovery": n_true_nonrecovery,
                "n_compatibility_excluded": n_compatibility_excluded,
                "n_incompatible": n_incompatible,
                "n_no_download_plan_match": n_no_download,
                "n_linked_run_without_fastq_ftp": n_no_fastq,
                "legacy_recovery_rate": legacy_rate,
                "resolved_recovery_rate": resolved_rate,
                "resolved_recovery_rate_ci_lower": ci_lower,
                "resolved_recovery_rate_ci_upper": ci_upper,
                "notes": (
                    "Resolved recovery excludes compatibility-excluded rows from the denominator; "
                    "legacy_recovery_rate reproduces the raw-count view for auditability."
                ),
            }
        )
    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary = summary.sort_values(["display_order", "family_key"]).reset_index(drop=True)
    return summary


def main() -> None:
    manifest = prepare_manifest()

    read_linked = build_read_linked_transportability(manifest)
    asr_representativeness = build_asr_representativeness_summary(manifest)
    local_support = build_local_origin_support_profile(manifest)
    within_origin = build_within_origin_structural_concentration(manifest)
    event_hierarchy = build_event_definition_hierarchy(manifest)
    detectability_detail = build_detectability_detail_rows(manifest)
    detectability_summary = build_detectability_summary(detectability_detail)
    missingness_visual = build_missingness_visual_summary()

    write_dual(read_linked, READ_LINKED_TRANSPORT_PATH, READ_LINKED_TRANSPORT_SUPP_PATH)
    write_dual(asr_representativeness, ASR_REPRESENTATIVENESS_PATH, ASR_REPRESENTATIVENESS_SUPP_PATH)
    write_dual(local_support, LOCAL_SUPPORT_PATH, LOCAL_SUPPORT_SUPP_PATH)
    write_dual(within_origin, WITHIN_ORIGIN_PATH, WITHIN_ORIGIN_SUPP_PATH)
    write_dual(event_hierarchy, EVENT_HIERARCHY_PATH, EVENT_HIERARCHY_SUPP_PATH)
    write_tsv(detectability_summary, DETECTABILITY_SUMMARY_PATH)
    write_tsv(detectability_detail, DETECTABILITY_DETAIL_PATH)
    write_tsv(missingness_visual, MISSINGNESS_VISUAL_PATH)

    print(f"Wrote read-linked transportability rows: {len(read_linked)}")
    print(f"Wrote ASR representativeness rows: {len(asr_representativeness)}")
    print(f"Wrote local package support rows: {len(local_support)}")
    print(f"Wrote within-origin concentration rows: {len(within_origin)}")
    print(f"Wrote event hierarchy rows: {len(event_hierarchy)}")
    print(f"Wrote detectability summary rows: {len(detectability_summary)}")
    print(f"Wrote detectability detail rows: {len(detectability_detail)}")
    print(f"Wrote missingness visual rows: {len(missingness_visual)}")


if __name__ == "__main__":
    main()
