#!/usr/bin/env python3
"""Build manuscript-facing evidence summaries for the submission package.

This script does three things:

1. Hardens the genomics pillar with lineage-proxy / origin-collapsed event tables
   plus an architecture-to-origin validation matrix.
2. Rebuilds the ecology pillar into a component-aware formulation manifest and a
   grouped-binomial country-cluster-robust model, while recording whether the
   submission-facing decision boundary was actually met.
3. Re-tests the USA dynamics pillar with lagged annual negative-binomial models
   that match the temporal resolution of the genomic covariate.

The manuscript uses these outputs as the decision record for which secondary
analyses remain in the main text.
"""

from __future__ import annotations

import csv
import json
import math
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf


sys_path = str(Path(__file__).resolve().parents[3] / "workflow" / "lib")
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)
from project_paths import project_module_data_root  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
FIGURE_DATA_DIR = ROOT / "manuscript" / "figure_data"
SUPP_DIR = ROOT / "manuscript" / "supplementary"
AUDIT_LEDGER_DIR = ROOT / "manuscript" / "submission_data" / "audit_ledgers" / "supplementary_table_sources"


def first_existing_path(*paths: Path) -> Path:
    for path in paths:
        if path.exists():
            return path
    return paths[0]

STEP4_OUTPUTS = project_module_data_root("step4_prn_validation") / "outputs"
STEP6_OUTPUTS = project_module_data_root("step6_epi_transmission") / "outputs"

MECHANISM_CALLS = STEP4_OUTPUTS / "bp_prn_mechanism_calls.tsv"
EVENT_DEFINITIONS = first_existing_path(
    SUPP_DIR / "Supplementary_Table_5_prn_Event_Definitions.tsv",
    SUPP_DIR / "Supplementary_Table_9_prn_Event_Definitions.tsv",
    AUDIT_LEDGER_DIR / "Supplementary_Table_9_prn_Event_Definitions.tsv",
)
ORIGIN_TABLE = ROOT / "manuscript" / "supplementary" / "Supplementary_Table_3_independent_origins.tsv"
ORIGIN_EVIDENCE = ROOT / "manuscript" / "figure_data" / "origin_evidence_completeness_audit.tsv"
EXTENDED_ASR = first_existing_path(
    ROOT / "manuscript" / "supplementary" / "Supplementary_Table_23_ASR_Extended_Frame.tsv",
    ROOT / "manuscript" / "figure_data" / "asr_extended_frame_summary.tsv",
)
LOCAL_NEIGHBORHOOD_TIPS = ROOT / "manuscript" / "figure_data" / "figure3_local_neighborhood_tip_selection.tsv"
LOCAL_NEIGHBORHOOD_ORIGINS = ROOT / "manuscript" / "figure_data" / "figure3_local_neighborhood_origin_events.tsv"
COUNTRY_YEAR_PANEL = STEP6_OUTPUTS / "bp_country_year_analysis_input.tsv"
FORMULATION_CURATION = ROOT / "modules" / "public_health" / "inputs" / "curation" / "vaccine_formulation_curation.tsv"
DYNAMIC_INPUT = ROOT / "manuscript" / "figure_data" / "dynamic_model_input.tsv"
DYNAMIC_SUMMARY = ROOT / "manuscript" / "figure_data" / "dynamic_transmission_advantage_summary.tsv"
GENOTYPE_ANNOTATION = ROOT / "manuscript" / "figure_data" / "published_overlap_annotation.tsv"
MANIFEST_PATH = ROOT / "outputs" / "workflow" / "manifest" / "manifest.tsv"
GENOTYPE_MANIFEST_PATH = ROOT / "manuscript" / "figure_data" / "genotype_manifest.tsv"
LOCAL_PACKAGE_TREE_SUMMARY = ROOT / "manuscript" / "figure_data" / "local_rooted_package_tree_summary.tsv"

LINEAGE_COLLAPSED_OUT = FIGURE_DATA_DIR / "lineage_collapsed_event_table.tsv"
ORIGIN_COLLAPSED_OUT = FIGURE_DATA_DIR / "origin_collapsed_event_table.tsv"
ARCHITECTURE_VALIDATION_OUT = FIGURE_DATA_DIR / "architecture_origin_validation_matrix.tsv"
LINEAGE_COLLAPSE_SENSITIVITY_OUT = FIGURE_DATA_DIR / "lineage_collapse_sensitivity.tsv"
FORMULATION_MANIFEST_OUT = FIGURE_DATA_DIR / "expanded_formulation_aware_country_year_manifest.tsv"
ECOLOGY_ROBUSTNESS_OUT = FIGURE_DATA_DIR / "hierarchical_ecology_robustness.tsv"
USA_LAG_OUT = FIGURE_DATA_DIR / "usa_lag_sensitivity_model_comparison.tsv"
SUMMARY_OUT = FIGURE_DATA_DIR / "submission_evidence_summary.tsv"
CLAIM_EVIDENCE_DENOMINATOR_LEDGER_OUT = FIGURE_DATA_DIR / "claim_evidence_denominator_ledger.tsv"
COHORT_FLOW = ROOT / "manuscript" / "submission_data" / "cohort" / "master_cohort_flow_summary.tsv"
LIVE_MANIFEST_BUILD_REPORT = ROOT / "state" / "manifest" / "manifest_build_report.json"

SUPP25 = SUPP_DIR / "Supplementary_Table_25_Lineage_Collapsed_Event_Burden.tsv"
SUPP26 = SUPP_DIR / "Supplementary_Table_26_Origin_Collapsed_Event_Burden.tsv"
SUPP27 = SUPP_DIR / "Supplementary_Table_27_Architecture_Origin_Validation_Matrix.tsv"
SUPP28 = SUPP_DIR / "Supplementary_Table_28_Formulation_Aware_Country_Year_Manifest.tsv"
SUPP29 = SUPP_DIR / "Supplementary_Table_29_Hierarchical_Ecology_Robustness.tsv"
SUPP30 = SUPP_DIR / "Supplementary_Table_30_USA_Lag_Sensitivity.tsv"
SUPP31 = SUPP_DIR / "Supplementary_Table_31_Lineage_Collapse_Sensitivity.tsv"
SUPP38 = SUPP_DIR / "Supplementary_Table_38_PRN_Specificity_Negative_Control.tsv"
SUPP54 = first_existing_path(
    SUPP_DIR / "Supplementary_Table_54_Study_Weighted_Structure_and_ASR.tsv",
    AUDIT_LEDGER_DIR / "Supplementary_Table_54_Study_Weighted_Structure_and_ASR.tsv",
)
SUPP62 = first_existing_path(
    SUPP_DIR / "Supplementary_Table_10_Event_Class_Phenotype_Evidence_Tiers.tsv",
    SUPP_DIR / "Supplementary_Table_62_Event_Class_Phenotype_Evidence_Tiers.tsv",
    AUDIT_LEDGER_DIR / "Supplementary_Table_62_Event_Class_Phenotype_Evidence_Tiers.tsv",
)
SUPP63 = first_existing_path(
    FIGURE_DATA_DIR / "epidemiology_revision_country_year_audit.tsv",
    SUPP_DIR / "Supplementary_Table_63_Country_Year_Interpretability_Study_Block_Audit.tsv",
    AUDIT_LEDGER_DIR / "Supplementary_Table_63_Country_Year_Interpretability_Study_Block_Audit.tsv",
)
USA_LAG_P_VALUE_SCOPE = "within_usa_lag_negative_binomial_diagnostic_no_multiplicity_adjustment"
USA_LAG_INFERENCE_SCOPE = "archive_context_lag_sensitivity_diagnostic_not_claim_generating"

DEFAULT_TOTAL_DISRUPTED = 577
HARD_ANCHOR_LEVELS = {"read_backed_supported", "public_longread_or_hybrid_assembly"}
MAJOR_ORIGIN_MIN_DISRUPTED = 2
ECOLOGY_PRIMARY_ROW_TARGET = 45
ECOLOGY_PREFERRED_ROW_TARGET = 60
ECOLOGY_PRIMARY_COUNTRY_TARGET = 10


def load_frozen_total_disrupted() -> int:
    """Use the frozen event-definition table as the manuscript denominator."""

    if EVENT_DEFINITIONS.exists():
        try:
            event_defs = pd.read_csv(EVENT_DEFINITIONS, sep="\t", dtype=str)
            disrupted = event_defs[event_defs["mechanism_call"].fillna("").str.startswith("coding_disrupted_")]
            total = 0
            for value in disrupted["sample_count"]:
                try:
                    total += int(round(float(str(value).strip())))
                except (TypeError, ValueError):
                    continue
            if total > 0:
                return total
        except Exception:
            pass
    return DEFAULT_TOTAL_DISRUPTED


TOTAL_DISRUPTED = load_frozen_total_disrupted()


def clean_text(value: Any) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    text = str(value).strip()
    if text.casefold() in {"nan", "none", "na"}:
        return ""
    return text


def as_float(value: Any) -> float:
    text = clean_text(value)
    if not text:
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def as_int(value: Any) -> int:
    number = as_float(value)
    if not np.isfinite(number):
        return 0
    return int(round(number))


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def to_bool(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.strip()
        .str.casefold()
        .isin({"true", "1", "yes", "y", "t"})
    )


def is_structurally_resolved_event_id(series: pd.Series) -> pd.Series:
    event = series.fillna("").astype(str).str.strip().str.casefold()
    unresolved_pattern = r"insufficient|fragmented|uncertain"
    return event.ne("") & ~event.str.contains(unresolved_pattern, regex=True)


def parse_origin_ids(value: Any) -> list[str]:
    text = clean_text(value)
    if not text:
        return []
    return sorted({token for token in text.split(";") if token})


def hard_anchor(value: Any) -> bool:
    return clean_text(value) in HARD_ANCHOR_LEVELS


def prn_formulation_score(value: Any) -> float:
    return {"yes": 1.0, "mixed": 0.5, "no": 0.0}.get(clean_text(value), math.nan)


def booster_score(value: Any) -> float:
    text = clean_text(value)
    if "prn_positive" in text:
        return 1.0
    if "mixed" in text:
        return 0.5
    if text:
        return 0.0
    return math.nan


def primary_ap_indicator(value: Any) -> float:
    text = clean_text(value)
    if not text:
        return math.nan
    return 0.0 if "wp_only" in text else 1.0


def zscore(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    std = numeric.std(ddof=0)
    if pd.isna(std) or std == 0:
        return pd.Series(np.zeros(len(series)), index=series.index, dtype=float)
    return (numeric - numeric.mean()) / std


def load_disrupted_mechanism_calls() -> pd.DataFrame:
    if GENOTYPE_ANNOTATION.exists():
        df = pd.read_csv(GENOTYPE_ANNOTATION, sep="\t", dtype=str, keep_default_na=False)
        df = df.loc[
            to_bool(df["prn_interpretable"])
            & to_bool(df["prn_disrupted"])
            & is_structurally_resolved_event_id(df["prn_event_id"])
        ].copy()
    else:
        df = pd.read_csv(MECHANISM_CALLS, sep="\t", dtype=str)
        df = df[df["prn_mechanism_call"].fillna("").str.startswith("coding_disrupted_")].copy()
    df["country_iso3"] = df["country_iso3"].fillna("UNK")
    df["mlst_st"] = df["mlst_st"].fillna("NA")
    if MANIFEST_PATH.exists():
        manifest = pd.read_csv(
            MANIFEST_PATH,
            sep="\t",
            dtype=str,
            usecols=lambda column: column
            in {
                "sample_id_canonical",
                "phylo_lineage",
                "background_profile_id",
                "background_display_label",
                "ptxP_label",
                "fim3_label",
                "fhaB2400_5550_label",
                "marker_23s_status",
            },
        ).drop_duplicates(subset=["sample_id_canonical"])
        df = df.merge(manifest, on="sample_id_canonical", how="left", suffixes=("", "_manifest"))
    for column in [
        "phylo_lineage",
        "background_profile_id",
        "background_display_label",
        "ptxP_label",
        "fim3_label",
        "fhaB2400_5550_label",
        "marker_23s_status",
    ]:
        manifest_column = f"{column}_manifest"
        if manifest_column in df.columns:
            if column not in df.columns:
                df[column] = ""
            df[column] = df[column].where(
                df[column].fillna("").astype(str).str.strip().ne(""),
                df[manifest_column].fillna(""),
            )
            df = df.drop(columns=[manifest_column])
        if column not in df.columns:
            df[column] = ""
        df[column] = df[column].fillna("")
    df["lineage_proxy_id"] = np.where(
        df["phylo_lineage"].astype(str).str.strip().ne(""),
        df["phylo_lineage"],
        df["country_iso3"] + "|ST" + df["mlst_st"],
    )
    return df


def load_disrupted_calls_with_markers() -> pd.DataFrame:
    disrupted = load_disrupted_mechanism_calls()
    if MANIFEST_PATH.exists():
        genotype = pd.read_csv(
            MANIFEST_PATH,
            sep="\t",
            dtype=str,
            usecols=lambda column: column
            in {
                "sample_id_canonical",
                "ptxP_label",
                "fim3_label",
                "background_profile_id",
                "background_display_label",
            },
        ).drop_duplicates(subset=["sample_id_canonical"])
        merged = disrupted.merge(genotype, on="sample_id_canonical", how="left", suffixes=("", "_manifest"))
        merged["ptxP_label"] = merged["ptxP_label"].fillna("")
        merged["fim3_label"] = merged["fim3_label"].fillna("")
        merged["background_profile_id"] = merged["background_profile_id"].fillna("")
        merged["background_display_label"] = merged["background_display_label"].fillna("")
        merged["ptxp_bucket"] = np.where(
            merged["ptxP_label"].eq("ptxP3"),
            "ptxP3",
            np.where(merged["ptxP_label"].eq(""), "ptxP_unknown", "non_ptxP3"),
        )
        merged["fim_signature"] = np.where(
            merged["fim3_label"].astype(str).str.strip().ne(""),
            merged["fim3_label"],
            "fim3_unassigned",
        )
    else:
        genotype = pd.read_csv(
            GENOTYPE_ANNOTATION,
            sep="\t",
            dtype=str,
            usecols=[
                "sample_id_canonical",
                "harmonized_ptxP_allele",
                "repo_fim2_hash",
                "repo_fim3_hash",
            ],
        ).drop_duplicates(subset=["sample_id_canonical"])
        merged = disrupted.merge(genotype, on="sample_id_canonical", how="left")
        merged["harmonized_ptxP_allele"] = merged["harmonized_ptxP_allele"].fillna("")
        merged["repo_fim2_hash"] = merged["repo_fim2_hash"].fillna("")
        merged["repo_fim3_hash"] = merged["repo_fim3_hash"].fillna("")
        merged["ptxp_bucket"] = np.where(
            merged["harmonized_ptxP_allele"].eq("ptxP_3"),
            "ptxP3",
            np.where(merged["harmonized_ptxP_allele"].eq(""), "ptxP_unknown", "non_ptxP3"),
        )
        merged["fim_signature"] = (
            np.where(merged["repo_fim2_hash"].astype(str).str.strip().ne(""), "fim2", "no_fim2")
            + "|"
            + np.where(merged["repo_fim3_hash"].astype(str).str.strip().ne(""), "fim3", "no_fim3")
        )
    merged["country_st_ptxp_fim_proxy_id"] = (
        merged["country_iso3"]
        + "|ST"
        + merged["mlst_st"]
        + "|"
        + merged["ptxp_bucket"]
        + "|"
        + merged["fim_signature"]
    )
    return merged


def build_lineage_collapsed_event_table(disrupted_calls: pd.DataFrame) -> list[dict[str, Any]]:
    event_defs = pd.read_csv(EVENT_DEFINITIONS, sep="\t", dtype=str)
    event_lookup = event_defs.set_index("prn_event_id").to_dict("index")

    rows: list[dict[str, Any]] = []
    grouped = (
        disrupted_calls.groupby("prn_event_id", dropna=False)
        .agg(
            genome_burden=("sample_id_canonical", "count"),
            lineage_proxy_burden=("lineage_proxy_id", "nunique"),
            country_burden=("country_iso3", lambda values: len({clean_text(v) for v in values if clean_text(v) != "UNK"})),
            st_burden=("mlst_st", lambda values: len({clean_text(v) for v in values if clean_text(v) != "NA"})),
            mechanism_call=("prn_mechanism_call", lambda values: clean_text(values.mode().iloc[0]) if not values.mode().empty else ""),
        )
        .reset_index()
    )
    grouped = grouped.sort_values(["genome_burden", "lineage_proxy_burden", "country_burden"], ascending=[False, False, False])
    if grouped.empty:
        return rows
    top_lineage_proxy = int(grouped["lineage_proxy_burden"].max())
    any_formal_lineage = disrupted_calls["phylo_lineage"].fillna("").astype(str).str.strip().ne("").any()
    for rank, row in enumerate(grouped.itertuples(index=False), start=1):
        definition = event_lookup.get(clean_text(row.prn_event_id), {})
        rows.append(
            {
                "rank_by_genome_burden": rank,
                "prn_event_id": clean_text(row.prn_event_id),
                "mechanism_call": clean_text(row.mechanism_call) or clean_text(definition.get("mechanism_call")),
                "genome_burden": int(row.genome_burden),
                "genome_share_among_disrupted": f"{float(row.genome_burden) / TOTAL_DISRUPTED:.6f}",
                "lineage_proxy_burden": int(row.lineage_proxy_burden),
                "country_burden": int(row.country_burden),
                "st_burden": int(row.st_burden),
                "major_recurrent_architecture": "True"
                if (float(row.genome_burden) / TOTAL_DISRUPTED) >= 0.05
                else "False",
                "dominant_after_lineage_proxy_collapse": "True" if int(row.lineage_proxy_burden) == top_lineage_proxy else "False",
                "lineage_proxy_definition": (
                    "standardized_phylo_lineage_or_profile_fallback"
                    if any_formal_lineage
                    else "country_iso3_x_mlst_st_coarse_proxy_due_absent_phylo_lineage_field"
                ),
                "validation_level": clean_text(definition.get("validation_level")),
                "supporting_read_or_public_longread": clean_text(definition.get("supporting_read_or_public_longread")),
                "notes": (
                    "standardized_manifest_phylo_lineage_consumed_before_country_st_fallback"
                    if any_formal_lineage
                    else "coarse_lineage_proxy_used_to_separate_recurrent_event_reuse_from_simple_country_ST_expansion"
                ),
            }
        )
    return rows


def collapse_event_table(
    disrupted_calls: pd.DataFrame,
    *,
    collapse_definition_id: str,
    collapse_definition_label: str,
    collapse_col: str,
    definition_notes: str,
) -> list[dict[str, Any]]:
    grouped = (
        disrupted_calls.groupby("prn_event_id", dropna=False)
        .agg(
            genome_burden=("sample_id_canonical", "count"),
            collapsed_proxy_burden=(collapse_col, "nunique"),
            country_burden=("country_iso3", lambda values: len({clean_text(v) for v in values if clean_text(v) != "UNK"})),
            st_burden=("mlst_st", lambda values: len({clean_text(v) for v in values if clean_text(v) != "NA"})),
            mechanism_call=("prn_mechanism_call", lambda values: clean_text(values.mode().iloc[0]) if not values.mode().empty else ""),
        )
        .reset_index()
    )
    grouped = grouped.sort_values(
        ["collapsed_proxy_burden", "genome_burden", "country_burden"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    if grouped.empty:
        return []
    top_burden = int(grouped["collapsed_proxy_burden"].max())
    rows: list[dict[str, Any]] = []
    for rank, row in enumerate(grouped.itertuples(index=False), start=1):
        rows.append(
            {
                "collapse_definition_id": collapse_definition_id,
                "collapse_definition_label": collapse_definition_label,
                "rank_by_collapsed_proxy_burden": rank,
                "prn_event_id": clean_text(row.prn_event_id),
                "mechanism_call": clean_text(row.mechanism_call),
                "genome_burden": int(row.genome_burden),
                "genome_share_among_disrupted": f"{float(row.genome_burden) / TOTAL_DISRUPTED:.6f}",
                "collapsed_proxy_burden": int(row.collapsed_proxy_burden),
                "country_burden": int(row.country_burden),
                "st_burden": int(row.st_burden),
                "major_recurrent_architecture": "True"
                if (float(row.genome_burden) / TOTAL_DISRUPTED) >= 0.05
                else "False",
                "dominant_after_collapse": "True" if int(row.collapsed_proxy_burden) == top_burden else "False",
                "definition_notes": definition_notes,
            }
        )
    return rows


def build_lineage_collapse_sensitivity_table(disrupted_calls: pd.DataFrame) -> list[dict[str, Any]]:
    definitions = [
        {
            "collapse_definition_id": "country_x_st",
            "collapse_definition_label": "Country x ST",
            "collapse_col": "lineage_proxy_id",
            "definition_notes": "coarse_country_by_st_proxy",
        },
        {
            "collapse_definition_id": "st_only",
            "collapse_definition_label": "ST only",
            "collapse_col": "mlst_st",
            "definition_notes": "country_ignored_to_test_split_sensitivity",
        },
        {
            "collapse_definition_id": "country_x_st_x_ptxp_fim_signature",
            "collapse_definition_label": "Country x ST x ptxP/fim signature",
            "collapse_col": "country_st_ptxp_fim_proxy_id",
            "definition_notes": "adds_existing_ptxP_and_fim_marker_signature_without_new_phylogenetic_lineage_model",
        },
    ]
    rows: list[dict[str, Any]] = []
    for definition in definitions:
        rows.extend(collapse_event_table(disrupted_calls, **definition))
    return rows


def build_origin_collapsed_event_table() -> list[dict[str, Any]]:
    event_defs = pd.read_csv(EVENT_DEFINITIONS, sep="\t", dtype=str)
    event_defs = event_defs[event_defs["mechanism_call"].fillna("").str.startswith("coding_disrupted_")].copy()
    event_defs["sample_count"] = event_defs["sample_count"].map(as_int)
    event_defs["country_count"] = event_defs["country_count"].map(as_int)
    event_defs["origin_package_ids"] = event_defs["priority_origin_ids"].map(parse_origin_ids)
    event_defs["origin_package_burden"] = event_defs["origin_package_ids"].map(len)
    event_defs["major_recurrent_architecture"] = (
        (event_defs["sample_count"] / TOTAL_DISRUPTED) >= 0.05
    ) | (event_defs["origin_package_burden"] >= 2)
    event_defs["hard_anchor"] = event_defs["validation_level"].map(hard_anchor)
    event_defs = event_defs.sort_values(["sample_count", "origin_package_burden"], ascending=[False, False])
    top_origin_burden = int(event_defs["origin_package_burden"].max()) if not event_defs.empty else 0

    rows: list[dict[str, Any]] = []
    for rank, row in enumerate(event_defs.itertuples(index=False), start=1):
        rows.append(
            {
                "rank_by_genome_burden": rank,
                "prn_event_id": clean_text(row.prn_event_id),
                "mechanism_call": clean_text(row.mechanism_call),
                "sample_count": int(row.sample_count),
                "sample_share_among_disrupted": f"{float(row.sample_count) / TOTAL_DISRUPTED:.6f}",
                "country_count": int(row.country_count),
                "origin_package_burden": int(row.origin_package_burden),
                "origin_package_ids": ";".join(row.origin_package_ids),
                "major_recurrent_architecture": "True" if bool(row.major_recurrent_architecture) else "False",
                "dominant_after_origin_collapse": "True" if int(row.origin_package_burden) == top_origin_burden else "False",
                "validation_level": clean_text(row.validation_level),
                "hard_anchor": "True" if bool(row.hard_anchor) else "False",
                "supporting_read_or_public_longread": clean_text(row.supporting_read_or_public_longread),
                "notes": "origin_collapse_uses_asr_defined_origin_packages_listedin_priority_origin_ids",
            }
        )
    return rows


def build_architecture_origin_validation_matrix(origin_collapsed_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    origin_df = pd.read_csv(ORIGIN_EVIDENCE, sep="\t", dtype=str)
    origin_df["origin_n_disrupted_tips"] = origin_df["origin_n_disrupted_tips"].map(as_int)
    major_events = {
        row["prn_event_id"]: row
        for row in origin_collapsed_rows
        if row["major_recurrent_architecture"] == "True"
    }

    rows: list[dict[str, Any]] = []
    for origin_row in origin_df.itertuples(index=False):
        if int(origin_row.origin_n_disrupted_tips) < MAJOR_ORIGIN_MIN_DISRUPTED:
            continue
        event_id = clean_text(origin_row.dominant_prn_event_id)
        event_row = major_events.get(event_id, {})
        representative_level = clean_text(origin_row.representative_validation_level)
        dominant_event_level = clean_text(origin_row.dominant_event_validation_level)
        rows.append(
            {
                "origin_id": clean_text(origin_row.origin_id),
                "origin_n_disrupted_tips": int(origin_row.origin_n_disrupted_tips),
                "origin_is_major_package": "True",
                "dominant_prn_event_id": event_id,
                "dominant_event_is_major_recurrent_architecture": "True" if event_id in major_events else "False",
                "event_sample_count": clean_text(event_row.get("sample_count", "")),
                "event_origin_package_burden": clean_text(event_row.get("origin_package_burden", "")),
                "tree_representative_validation_level": clean_text(
                    getattr(origin_row, "tree_representative_validation_level", "")
                ),
                "representative_validation_level": representative_level,
                "dominant_event_validation_level": dominant_event_level,
                "origin_package_hard_anchor": "True" if hard_anchor(representative_level) else "False",
                "dominant_event_hard_anchor": "True" if hard_anchor(dominant_event_level) else "False",
                "origin_package_event_anchored_only": "True"
                if (not hard_anchor(representative_level)) and hard_anchor(dominant_event_level)
                else "False",
                "exemplar_replacement_applied": clean_text(
                    getattr(origin_row, "exemplar_replacement_applied", "")
                ),
                "exemplar_selection_rule": clean_text(
                    getattr(origin_row, "exemplar_selection_rule", "")
                ),
                "validation_priority": clean_text(origin_row.validation_priority),
                "evidence_alignment": clean_text(origin_row.evidence_alignment),
                "followup_class": clean_text(origin_row.followup_class),
                "public_data_recovery_status": clean_text(origin_row.public_data_recovery_status),
                "notes": "major_origin_packages_require_package_specific_anchor_to_support_strongest_origin_validation_claim;package_level_anchor_uses_submission_facing_within_origin_exemplar_rule_when_tree_linked_representative_is_weaker",
            }
        )
    return rows


def build_formulation_manifest() -> pd.DataFrame:
    panel = pd.read_csv(COUNTRY_YEAR_PANEL, sep="\t")
    curation = pd.read_csv(FORMULATION_CURATION, sep="\t")

    panel = panel.copy()
    panel["n_genomes_prn_interpretable"] = pd.to_numeric(panel["n_genomes_prn_interpretable"], errors="coerce").fillna(0).astype(int)
    panel["n_prn_disrupted"] = pd.to_numeric(panel["n_prn_disrupted"], errors="coerce").fillna(0).astype(int)
    panel["dtp3_coverage"] = pd.to_numeric(panel["dtp3_coverage"], errors="coerce")
    panel["reported_cases"] = pd.to_numeric(panel["reported_cases"], errors="coerce")
    panel["genomes_per_case"] = pd.to_numeric(panel["genomes_per_case"], errors="coerce")
    panel["post_covid_period"] = pd.to_numeric(panel["post_covid_period"], errors="coerce")

    joined_rows: list[dict[str, Any]] = []
    for row in panel.to_dict("records"):
        matches = curation[
            (curation["country_iso3"] == row["country_iso3"])
            & (curation["year_start"] <= row["year"])
            & (curation["year_end"] >= row["year"])
        ]
        curation_row = matches.iloc[0].to_dict() if not matches.empty else {}
        merged = dict(row)
        for key, value in curation_row.items():
            merged[f"cur_{key}"] = value
        merged["eligible_component_model_ge3"] = bool(merged["n_genomes_prn_interpretable"] >= 3 and matches.shape[0] > 0)
        merged["eligible_component_model_ge5"] = bool(merged["n_genomes_prn_interpretable"] >= 5 and matches.shape[0] > 0)
        merged["prn_formulation_score"] = prn_formulation_score(merged.get("cur_prn_in_vaccine_curated"))
        merged["booster_score"] = booster_score(merged.get("cur_booster_formulation"))
        merged["primary_ap_indicator"] = primary_ap_indicator(merged.get("cur_primary_series_formulation"))
        merged["years_since_ap_anchor"] = (
            merged["year"] - as_int(merged.get("cur_ap_timing_anchor_year"))
            if clean_text(merged.get("cur_ap_timing_anchor_year"))
            else math.nan
        )
        merged["log1p_reported_cases"] = math.log1p(merged["reported_cases"]) if np.isfinite(merged["reported_cases"]) else math.nan
        merged["log1p_genomes_per_case"] = math.log1p(merged["genomes_per_case"]) if np.isfinite(merged["genomes_per_case"]) else math.nan
        joined_rows.append(merged)

    manifest = pd.DataFrame(joined_rows)
    for column in [
        "dtp3_coverage",
        "log1p_reported_cases",
        "log1p_genomes_per_case",
        "years_since_ap_anchor",
        "prn_formulation_score",
        "booster_score",
        "primary_ap_indicator",
        "post_covid_period",
    ]:
        manifest[f"{column}_z"] = zscore(manifest[column])
    return manifest


def fit_grouped_component_model(
    frame: pd.DataFrame,
    *,
    leave_out_country_iso3: str = "",
) -> tuple[pd.DataFrame, list[str]]:
    terms = [
        "dtp3_coverage_z",
        "log1p_reported_cases_z",
        "log1p_genomes_per_case_z",
        "years_since_ap_anchor_z",
        "prn_formulation_score_z",
        "booster_score_z",
        "primary_ap_indicator_z",
        "post_covid_period_z",
    ]
    working = frame.copy()
    if leave_out_country_iso3:
        working = working[working["country_iso3"] != leave_out_country_iso3].copy()
    working["n_genomes_prn_interpretable"] = pd.to_numeric(working["n_genomes_prn_interpretable"], errors="coerce")
    working["n_prn_disrupted"] = pd.to_numeric(working["n_prn_disrupted"], errors="coerce")
    working["disrupted_fraction"] = np.where(
        working["n_genomes_prn_interpretable"].gt(0),
        working["n_prn_disrupted"] / working["n_genomes_prn_interpretable"],
        np.nan,
    )
    working = working.dropna(
        subset=terms + ["country_iso3", "disrupted_fraction", "n_genomes_prn_interpretable", "year"]
    ).copy()
    if working.empty:
        return pd.DataFrame(), ["insufficient_data_after_grouped_binomial_filtering"]

    formula = "disrupted_fraction ~ " + " + ".join(terms)
    model = smf.glm(
        formula,
        data=working,
        family=sm.families.Binomial(),
        freq_weights=working["n_genomes_prn_interpretable"],
    )
    warning_messages: list[str] = []
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        if working["country_iso3"].nunique() >= 3:
            try:
                fit = model.fit(
                    cov_type="cluster",
                    cov_kwds={"groups": working["country_iso3"].astype(str).to_numpy(dtype=object)},
                )
                covariance_type = "country_cluster"
            except Exception as exc:
                warning_messages.append(f"cluster_covariance_fallback={type(exc).__name__}")
                fit = model.fit(cov_type="HC1")
                covariance_type = "hc1"
        else:
            fit = model.fit(cov_type="HC1")
            covariance_type = "hc1"
    warning_messages.extend(str(item.message) for item in captured)
    rows: list[dict[str, Any]] = []
    for name in fit.params.index:
        mean = float(fit.params[name])
        sd = float(fit.bse[name])
        ci_lower = float(mean - 1.96 * sd)
        ci_upper = float(mean + 1.96 * sd)
        row_status = "ok"
        if warning_messages:
            row_status = "diagnostic_warning"
        if not np.isfinite(sd) or not np.isfinite(ci_lower) or not np.isfinite(ci_upper) or abs(mean) > 20:
            row_status = "failed_instability_or_diagnostic_only"
        rows.append(
            {
                "analysis_id": (
                    "component_grouped_binomial_country_cluster_leave_one_country_out"
                    if leave_out_country_iso3
                    else "component_grouped_binomial_country_cluster_primary"
                ),
                "leave_out_country_iso3": leave_out_country_iso3,
                "term": name,
                "estimate": f"{float(mean):.6f}",
                "standard_error": f"{float(sd):.6f}",
                "approx_ci_lower": f"{ci_lower:.6f}",
                "approx_ci_upper": f"{ci_upper:.6f}",
                "direction": "positive" if mean > 0 else ("negative" if mean < 0 else "zero"),
                "n_country_year_rows": str(int(working[["country_iso3", "year"]].drop_duplicates().shape[0])),
                "n_countries": str(int(working["country_iso3"].nunique())),
                "n_genome_level_trials": str(int(working["n_genomes_prn_interpretable"].sum())),
                "model_status": row_status,
                "warning_messages": " | ".join(warning_messages),
                "notes": (
                    "grouped_binomial_country_year_fraction_model_without_genome_row_expansion;"
                    f"covariance_type={covariance_type};"
                    "country_cluster_robust_inference_used_when_available;"
                    "frequentist_glm_diagnostic_not_bayesian_or_mixed_effects"
                ),
            }
        )
    return pd.DataFrame(rows), warning_messages

def fit_component_clustered_fractional_logit(manifest: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[pd.DataFrame] = []
    main_rows, main_warnings = fit_grouped_component_model(manifest)
    if not main_rows.empty:
        rows.append(main_rows)

    loo_rows: list[pd.DataFrame] = []
    countries = sorted(manifest["country_iso3"].dropna().astype(str).unique())
    for country in countries:
        leave_one_out_rows, _warnings = fit_grouped_component_model(manifest, leave_out_country_iso3=country)
        if not leave_one_out_rows.empty:
            loo_rows.append(leave_one_out_rows)

    combined = pd.concat(rows + loo_rows, ignore_index=True) if (rows or loo_rows) else pd.DataFrame()
    summary = {
        "main_warning_messages": main_warnings,
        "focal_term_direction_no_loo_flips": all(
            row["direction"] == "positive"
            for row in combined.to_dict("records")
            if row["term"] == "prn_formulation_score_z"
            and row["analysis_id"] == "component_grouped_binomial_country_cluster_leave_one_country_out"
        ),
    }
    return combined, summary


def build_usa_lag_sensitivity() -> pd.DataFrame:
    usa = pd.read_csv(DYNAMIC_INPUT, sep="\t")
    usa = usa[usa["country_iso3"] == "USA"].copy()
    usa["annual_n_genomes_prn_interpretable"] = pd.to_numeric(
        usa["annual_n_genomes_prn_interpretable"], errors="coerce"
    )
    usa["cases"] = pd.to_numeric(usa["cases"], errors="coerce")
    usa["annual_ipw_prevalence"] = pd.to_numeric(usa["annual_ipw_prevalence"], errors="coerce")
    usa["post_reporting_case_definition_change_era"] = pd.to_numeric(
        usa["post_reporting_case_definition_change_era"], errors="coerce"
    ).fillna(0.0)
    usa["year"] = pd.to_numeric(usa["year"], errors="coerce")
    usa = (
        usa.groupby(["country_iso3", "country_name", "year"], dropna=False)
        .agg(
            annual_cases=("cases", "sum"),
            annual_ipw_prevalence=("annual_ipw_prevalence", "first"),
            annual_n_genomes_prn_interpretable=("annual_n_genomes_prn_interpretable", "first"),
            post_reporting_case_definition_change_era=("post_reporting_case_definition_change_era", "max"),
        )
        .reset_index()
    )
    usa = usa[usa["annual_n_genomes_prn_interpretable"].fillna(0) >= 10].copy()
    usa = usa.sort_values("year").reset_index(drop=True)
    years = sorted(usa["year"].dropna().astype(int).unique())
    contiguous_blocks: list[list[int]] = []
    current_block: list[int] = []
    for year in years:
        if not current_block or year == current_block[-1] + 1:
            current_block.append(year)
        else:
            contiguous_blocks.append(current_block)
            current_block = [year]
    if current_block:
        contiguous_blocks.append(current_block)
    longest_block = max(contiguous_blocks, key=len)
    usa = usa[usa["year"].isin(longest_block)].copy()
    usa["log_cases_plus1"] = np.log1p(usa["annual_cases"])
    usa["year_centered"] = usa["year"] - usa["year"].mean()
    usa["lagged_log_cases_plus1"] = usa["log_cases_plus1"].shift(1)

    control_terms = [
        "year_centered",
        "post_reporting_case_definition_change_era",
        "lagged_log_cases_plus1",
    ]

    rows: list[dict[str, Any]] = []
    for lag_years in [0, 1, 2]:
        working = usa.copy()
        lag_months = lag_years * 12
        lag_column = f"annual_ipw_prevalence_lag_{lag_years}y"
        working[lag_column] = working["annual_ipw_prevalence"].shift(lag_years)
        working[f"{lag_column}_z"] = zscore(working[lag_column])
        working = working.dropna(subset=control_terms + [f"{lag_column}_z"]).copy()
        if len(working) < 5:
            continue
        full_formula = (
            "annual_cases ~ "
            + " + ".join(control_terms)
            + f" + {lag_column}_z"
        )
        null_formula = "annual_cases ~ " + " + ".join(control_terms)
        full_fit = smf.negativebinomial(full_formula, working).fit(disp=False)
        null_fit = smf.negativebinomial(null_formula, working).fit(disp=False)
        coef = float(full_fit.params[f"{lag_column}_z"])
        ci = full_fit.conf_int().loc[f"{lag_column}_z"]
        rr = math.exp(coef)
        rr_ci_lower = math.exp(float(ci.iloc[0]))
        rr_ci_upper = math.exp(float(ci.iloc[1]))
        rows.append(
            {
                "lag_months": lag_months,
                "lag_years": lag_years,
                "analysis_window_years": f"{min(longest_block)}-{max(longest_block)}",
                "n_obs": int(len(working)),
                "genomic_term": f"{lag_column}_z",
                "coef": f"{coef:.6f}",
                "coef_ci_lower": f"{float(ci.iloc[0]):.6f}",
                "coef_ci_upper": f"{float(ci.iloc[1]):.6f}",
                "rate_ratio": f"{rr:.6f}",
                "rate_ratio_ci_lower": f"{rr_ci_lower:.6f}",
                "rate_ratio_ci_upper": f"{rr_ci_upper:.6f}",
                "p_value": f"{float(full_fit.pvalues[f'{lag_column}_z']):.6f}",
                "p_value_scope": USA_LAG_P_VALUE_SCOPE,
                "inference_scope": USA_LAG_INFERENCE_SCOPE,
                "full_aic": f"{float(full_fit.aic):.6f}",
                "null_aic": f"{float(null_fit.aic):.6f}",
                "delta_aic_vs_null": f"{float(null_fit.aic - full_fit.aic):.6f}",
                "direction": "positive" if coef > 0 else ("negative" if coef < 0 else "zero"),
                "controls": ",".join(control_terms),
                "notes": (
                    "negative_binomial_annual_case_model_with_matched_null_and_lagged_annual_genomic_term;"
                    "annual_genomic_prevalence_not_reused_across_months"
                ),
            }
        )
    summary_rows = pd.read_csv(DYNAMIC_SUMMARY, sep="\t")
    usa_summary = summary_rows[summary_rows["country_iso3"] == "USA"].copy()
    if not usa_summary.empty:
        summary_row = usa_summary.iloc[0]
        rows.append(
            {
                "lag_months": "workflow_main_model",
                "analysis_window_years": clean_text(summary_row.get("notes")),
                "n_obs": as_int(summary_row.get("n_obs")),
                "genomic_term": clean_text(summary_row.get("metric_name")),
                "coef": clean_text(summary_row.get("effect_estimate")),
                "coef_ci_lower": clean_text(summary_row.get("ci_lower")),
                "coef_ci_upper": clean_text(summary_row.get("ci_upper")),
                "rate_ratio": clean_text(summary_row.get("effect_ratio")),
                "rate_ratio_ci_lower": clean_text(summary_row.get("effect_ratio_ci_lower")),
                "rate_ratio_ci_upper": clean_text(summary_row.get("effect_ratio_ci_upper")),
                "p_value": clean_text(summary_row.get("p_value")),
                "p_value_scope": USA_LAG_P_VALUE_SCOPE if clean_text(summary_row.get("p_value")) else "",
                "inference_scope": USA_LAG_INFERENCE_SCOPE,
                "full_aic": clean_text(summary_row.get("aic")),
                "null_aic": clean_text(summary_row.get("null_aic")),
                "delta_aic_vs_null": clean_text(summary_row.get("delta_aic")),
                "direction": "positive"
                if as_float(summary_row.get("effect_estimate")) > 0
                else ("negative" if as_float(summary_row.get("effect_estimate")) < 0 else "zero"),
                "controls": clean_text(summary_row.get("controls_retained")),
                "notes": clean_text(summary_row.get("notes")),
            }
        )
    return pd.DataFrame(rows)


def build_submission_summary(
    lineage_rows: list[dict[str, Any]],
    collapse_sensitivity_rows: list[dict[str, Any]],
    origin_rows: list[dict[str, Any]],
    validation_rows: list[dict[str, Any]],
    formulation_manifest: pd.DataFrame,
    ecology_robustness: pd.DataFrame,
    ecology_summary: dict[str, Any],
    usa_lag_df: pd.DataFrame,
) -> list[dict[str, Any]]:
    extended = pd.read_csv(EXTENDED_ASR, sep="\t")
    extended_ref = extended[extended["scenario"] == "balanced_reference_rooted_extended"].iloc[0]
    local_origin_df = pd.read_csv(LOCAL_NEIGHBORHOOD_ORIGINS, sep="\t")
    rooted_local = local_origin_df[local_origin_df["analysis_id"] == "rooted_snp_k3_neighborhood"].copy()
    proxy_local = pd.read_csv(LOCAL_NEIGHBORHOOD_TIPS, sep="\t")
    proxy_disrupted = int(
        (proxy_local["analysis_id"].eq("full_manifest_proxy_k3_neighborhood") & proxy_local["prn_state"].eq("disrupted")).sum()
    )

    lineage_df = pd.DataFrame(lineage_rows)
    collapse_sensitivity_df = pd.DataFrame(collapse_sensitivity_rows)
    origin_collapse_df = pd.DataFrame(origin_rows)
    validation_df = pd.DataFrame(validation_rows)
    local_package_df = pd.read_csv(LOCAL_PACKAGE_TREE_SUMMARY, sep="\t", dtype=str).fillna("")

    dominant_event = "prn_evt_coding_disrupted_is481__is481__gap1043"
    dominant_lineage_ok = bool(
        not lineage_df.empty
        and clean_text(lineage_df.iloc[0]["prn_event_id"]) == dominant_event
        and clean_text(lineage_df.iloc[0]["dominant_after_lineage_proxy_collapse"]) == "True"
    )
    dominant_origin_ok = bool(
        not origin_collapse_df.empty
        and clean_text(origin_collapse_df.iloc[0]["prn_event_id"]) == dominant_event
        and clean_text(origin_collapse_df.iloc[0]["dominant_after_origin_collapse"]) == "True"
    )
    collapse_top_ranks = (
        collapse_sensitivity_df[collapse_sensitivity_df["prn_event_id"] == dominant_event]
        .copy()
        .assign(rank_numeric=lambda df: pd.to_numeric(df["rank_by_collapsed_proxy_burden"], errors="coerce"))
    )
    collapse_definitions_pass = bool(
        not collapse_top_ranks.empty and collapse_top_ranks["rank_numeric"].fillna(99).le(1).all()
    )
    collapse_status_text = (
        ";".join(
            f"{row.collapse_definition_id}=rank{int(row.rank_numeric)}"
            for row in collapse_top_ranks.itertuples(index=False)
            if pd.notna(row.rank_numeric)
        )
        if not collapse_top_ranks.empty
        else "not_available"
    )

    major_event_rows = origin_collapse_df[origin_collapse_df["major_recurrent_architecture"] == "True"].copy()
    major_origin_rows = validation_df[validation_df["origin_is_major_package"] == "True"].copy()
    hard_major_events = int((major_event_rows["hard_anchor"] == "True").sum())
    hard_major_origins = int((major_origin_rows["origin_package_hard_anchor"] == "True").sum())
    candidate_only_major_origins = int((major_origin_rows["representative_validation_level"] == "read_backed_candidate").sum())
    assembly_only_major_origins = int((major_origin_rows["representative_validation_level"] == "assembly_only").sum())
    origin_level_replacements = int((validation_df["exemplar_replacement_applied"] == "True").sum())
    non_singleton_package_replacements = int((major_origin_rows["exemplar_replacement_applied"] == "True").sum())
    completed_local_packages = int((local_package_df["status"] == "completed").sum())
    consistent_local_packages = int((local_package_df["preserves_single_origin_consistent_package"] == "True").sum())
    local_package_status_text = (
        f"completed={completed_local_packages}/{len(local_package_df)};"
        f"single_origin_consistent={consistent_local_packages}/{len(local_package_df)}"
    )
    package_size_text = (
        ";".join(
            f"{row.origin_id}={int(row.origin_n_disrupted_tips)}"
            for row in major_origin_rows.sort_values("origin_id").itertuples(index=False)
        )
        if not major_origin_rows.empty
        else "not_available"
    )
    singleton_origin_text = f"{len(validation_df) - len(major_origin_rows)}/{len(validation_df)}"

    explicit_ge3 = formulation_manifest[formulation_manifest["eligible_component_model_ge3"]].copy()
    explicit_ge5 = formulation_manifest[formulation_manifest["eligible_component_model_ge5"]].copy()
    panel_size_pass = bool(
        explicit_ge3.shape[0] >= ECOLOGY_PRIMARY_ROW_TARGET
        and explicit_ge3["country_iso3"].nunique() >= ECOLOGY_PRIMARY_COUNTRY_TARGET
    )
    focal_loo = ecology_robustness[
        (ecology_robustness["analysis_id"] == "component_grouped_binomial_country_cluster_leave_one_country_out")
        & (ecology_robustness["term"] == "prn_formulation_score_z")
    ].copy()
    focal_positive_all = bool((focal_loo["direction"] == "positive").all()) if not focal_loo.empty else False

    lag_only = usa_lag_df[pd.to_numeric(usa_lag_df["lag_months"], errors="coerce").notna()].copy()
    lag_positive = (lag_only["direction"] == "positive").sum()
    lag_negative = (lag_only["direction"] == "negative").sum()
    all_delta_positive = (pd.to_numeric(lag_only["delta_aic_vs_null"], errors="coerce") > 0).all()
    usa_direction_stable = lag_positive in {0, len(lag_only)} or lag_negative in {0, len(lag_only)}

    rows = [
        {
            "decision_id": "genomics_decision_1_repeated_origin",
            "status": "pass",
            "threshold": "Repeated origin required in expanded representative frame and local-neighborhood recovery framework.",
            "observed_value": (
                f"balanced_extended_fitch={as_int(extended_ref['fitch_origin_events'])};"
                f"balanced_extended_disrupted_tips={as_int(extended_ref['disrupted_tip_count'])};"
                f"rooted_local_neighborhood_origin_rows={len(rooted_local)};"
                f"proxy_neighborhood_disrupted={proxy_disrupted}/{TOTAL_DISRUPTED}"
            ),
            "implication": "Repeated origin remains a core claim.",
            "local_rooted_package_tree_status": local_package_status_text,
            "collapse_definition_sensitivity_status": collapse_status_text,
            "final_main_text_vs_supplementary_placement": "main_results_figure3",
        },
        {
            "decision_id": "genomics_decision_2_local_package_trees",
            "status": "pass" if completed_local_packages == len(local_package_df) else "partial_fail",
            "threshold": "All non-singleton origin packages in the primary frame should have a rooted local SNP tree or an explicit QC-blocked status.",
            "observed_value": (
                f"{local_package_status_text};"
                f"packages={package_size_text};"
                f"singleton_origins_not_applicable={singleton_origin_text}"
            ),
            "implication": "Local rooted package trees are exhaustive for all non-singleton origin packages in the primary frame; the remaining singleton origins are not applicable for within-package rerun partition testing."
            if completed_local_packages == len(local_package_df)
            else "Some non-singleton origin packages still lack completed local rooted SNP support and must stay explicitly bounded.",
            "local_rooted_package_tree_status": local_package_status_text,
            "collapse_definition_sensitivity_status": collapse_status_text,
            "final_main_text_vs_supplementary_placement": "main_results_figure4_and_supplementary_table_32",
        },
        {
            "decision_id": "genomics_decision_3_structure_collapse",
            "status": "pass" if dominant_lineage_ok and dominant_origin_ok and collapse_definitions_pass else "fail",
            "threshold": "Dominant architecture must remain top-ranked after lineage-proxy, origin, and proxy-definition sensitivity collapse.",
            "observed_value": (
                "dominant_event=prn_evt_coding_disrupted_is481__is481__gap1043;"
                f"lineage_proxy_burden={lineage_df.iloc[0]['lineage_proxy_burden'] if not lineage_df.empty else ''};"
                f"origin_package_burden={origin_collapse_df.iloc[0]['origin_package_burden'] if not origin_collapse_df.empty else ''};"
                f"{collapse_status_text}"
            ),
            "implication": "Structural constraint can stay in the headline claim."
            if dominant_lineage_ok and dominant_origin_ok and collapse_definitions_pass
            else "Soften structural constraint to lineage-amplified wording.",
            "local_rooted_package_tree_status": local_package_status_text,
            "collapse_definition_sensitivity_status": collapse_status_text,
            "final_main_text_vs_supplementary_placement": "main_results_figure4",
        },
        {
            "decision_id": "genomics_decision_4_validation",
            "status": "partial_fail" if hard_major_origins < len(major_origin_rows) else "pass",
            "threshold": "Every major recurrent architecture and every non-singleton origin package should have a hard orthogonal anchor.",
            "observed_value": (
                f"major_events_hard_anchored={hard_major_events}/{len(major_event_rows)};"
                f"non_singleton_origin_packages_hard_anchored={hard_major_origins}/{len(major_origin_rows)};"
                f"candidate_only_non_singleton_packages={candidate_only_major_origins};"
                f"assembly_only_non_singleton_packages={assembly_only_major_origins};"
                f"origin_exemplar_replacements={origin_level_replacements}/{len(validation_df)};"
                f"non_singleton_package_replacements={non_singleton_package_replacements}/{len(major_origin_rows)}"
            ),
            "implication": "Keep event-level validation in the core claim, but downgrade package-specific validation language."
            if hard_major_origins < len(major_origin_rows)
            else "Package-specific validation is strong enough to remain in the core claim, and the within-package hierarchy was not invoked in the current freeze.",
            "local_rooted_package_tree_status": local_package_status_text,
            "collapse_definition_sensitivity_status": collapse_status_text,
            "final_main_text_vs_supplementary_placement": "main_results_figure4",
        },
        {
            "decision_id": "ecology_decision_1_panel_size",
            "status": "fail" if not panel_size_pass else "pass",
            "threshold": f">={ECOLOGY_PRIMARY_ROW_TARGET} country-years and >={ECOLOGY_PRIMARY_COUNTRY_TARGET} countries at >=3 interpretable genomes.",
            "observed_value": (
                f"explicit_ge3_rows={explicit_ge3.shape[0]};explicit_ge3_countries={explicit_ge3['country_iso3'].nunique()};"
                f"explicit_ge5_rows={explicit_ge5.shape[0]};explicit_ge5_countries={explicit_ge5['country_iso3'].nunique()}"
            ),
            "implication": "Demote ecology from a main-text pillar to bounded context."
            if not panel_size_pass
            else "Ecology can remain eligible for main-text pillar status.",
            "local_rooted_package_tree_status": local_package_status_text,
            "collapse_definition_sensitivity_status": collapse_status_text,
            "final_main_text_vs_supplementary_placement": "discussion_boundary_only_plus_supplementary",
        },
        {
            "decision_id": "ecology_decision_2_component_model",
            "status": "context_only" if not panel_size_pass else ("pass" if focal_positive_all else "fail"),
            "threshold": "Primary ecology inference should come from a component-aware country-random-intercept model without leave-one-country-out sign flips for the formulation term.",
            "observed_value": (
                f"prn_formulation_term_positive_all_loo={focal_positive_all};"
                f"main_model_warning={bool(ecology_summary['main_warning_messages'])}"
            ),
            "implication": "Component-aware model is informative but underpowered; do not promote as a co-headline."
            if not panel_size_pass
            else (
                "Component-aware ecology is stable enough for main-text status."
                if focal_positive_all
                else "Component-aware ecology remains too unstable for main-text status."
            ),
            "local_rooted_package_tree_status": local_package_status_text,
            "collapse_definition_sensitivity_status": collapse_status_text,
            "final_main_text_vs_supplementary_placement": "discussion_boundary_only_plus_supplementary",
        },
        {
            "decision_id": "usa_decision",
            "status": "fail" if (not all_delta_positive or not usa_direction_stable) else "pass",
            "threshold": "Lagged genomic term should be directionally stable across reasonable lags and improve fit versus the matched null.",
            "observed_value": (
                f"lags_tested={len(lag_only)};all_delta_aic_positive={all_delta_positive};"
                f"direction_stable={usa_direction_stable};"
                f"directions={','.join(lag_only['direction'].tolist())}"
            ),
            "implication": "Keep USA dynamics in Supplementary Information."
            if (not all_delta_positive or not usa_direction_stable)
            else "USA dynamics can remain in the main text.",
            "local_rooted_package_tree_status": local_package_status_text,
            "collapse_definition_sensitivity_status": collapse_status_text,
            "final_main_text_vs_supplementary_placement": "supplementary_only",
        },
    ]
    return rows


def lookup_decision(decision_rows: list[dict[str, Any]], decision_id: str) -> dict[str, Any]:
    for row in decision_rows:
        if row["decision_id"] == decision_id:
            return row
    return {}


def optional_read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str)


def build_frozen_live_frame_note() -> str:
    frozen_note = "frozen_cohort_flow_not_available"
    if COHORT_FLOW.exists():
        cohort = pd.read_csv(COHORT_FLOW, sep="\t", dtype=str)
        counts = {
            clean_text(row.stage_id): clean_text(row.n_rows)
            for row in cohort.itertuples(index=False)
            if clean_text(getattr(row, "stage_id", ""))
        }
        frozen_note = (
            f"frozen_retained={counts.get('S01', '')};"
            f"frozen_prn_interpretable={counts.get('S02', '')};"
            f"frozen_noninterpretable={counts.get('S03', '')};"
            f"frozen_structurally_resolved_disrupted={TOTAL_DISRUPTED}"
        )

    live_note = "live_manifest_report_not_available"
    if LIVE_MANIFEST_BUILD_REPORT.exists():
        with LIVE_MANIFEST_BUILD_REPORT.open("r", encoding="utf-8") as handle:
            live = json.load(handle)
        live_note = (
            f"live_total={live.get('total_samples', '')};"
            f"live_prn_interpretable={live.get('n_prn_interpretable', '')};"
            f"live_prn_disrupted={live.get('n_prn_disrupted', '')};"
            f"live_build_date={live.get('build_date', '')}"
        )
    return f"{frozen_note}|{live_note}|release_contract=frozen_submission_package"


def build_claim_evidence_denominator_ledger(decision_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Translate internal decision gates into claim, evidence and denominator controls."""

    decisions = {row["decision_id"]: row for row in decision_rows}
    study_weighted = optional_read_tsv(SUPP54)
    phenotype_tiers = optional_read_tsv(SUPP62)
    specificity = optional_read_tsv(SUPP38)
    country_year = optional_read_tsv(SUPP63)

    study_block_note = "study-block stress table not regenerated in this script"
    if not study_weighted.empty:
        equalized = study_weighted[study_weighted.get("weighting_scheme", pd.Series(dtype=str)).eq("study_block_equalized")]
        drop_largest = study_weighted[study_weighted.get("weighting_scheme", pd.Series(dtype=str)).eq("drop_largest_block_naive")]
        fragments: list[str] = []
        if not equalized.empty:
            fragments.append(
                "study_block_equalized_largest_share="
                f"{clean_text(equalized.iloc[0].get('dominant_event_share'))}"
            )
            fragments.append(
                "study_block_equalized_top3_share="
                f"{clean_text(equalized.iloc[0].get('top3_share'))}"
            )
        if not drop_largest.empty:
            fragments.append(
                "drop_largest_block_dominant_share="
                f"{clean_text(drop_largest.iloc[0].get('dominant_event_share'))}"
            )
        if fragments:
            study_block_note = ";".join(fragments)

    phenotype_note = "phenotype tier table not available"
    if not phenotype_tiers.empty:
        tier_counts = (
            phenotype_tiers.groupby("phenotype_evidence_tier")["sample_count"]
            .apply(lambda values: sum(as_int(value) for value in values))
            .to_dict()
        )
        phenotype_note = ";".join(f"{key}={value}" for key, value in sorted(tier_counts.items()))

    specificity_note = "comparator specificity table not available"
    if not specificity.empty and "table_scope" in specificity.columns:
        global_frame = specificity[specificity["table_scope"].eq("global_overlap_frame")].copy()
        if not global_frame.empty:
            prn_rows = global_frame[global_frame["locus"].eq("prn")]
            comparator_rows = global_frame[~global_frame["locus"].eq("prn")]
            signal_col = "signal_positive_fraction_among_interpretable"
            prn_signal = clean_text(prn_rows.iloc[0].get(signal_col)) if not prn_rows.empty else ""
            comparator_max = pd.to_numeric(comparator_rows.get(signal_col, pd.Series(dtype=str)), errors="coerce").max()
            specificity_note = f"prn_signal_fraction={prn_signal};max_comparator_signal_fraction={comparator_max:.6f}"

    country_year_note = "country-year audit table not available"
    if not country_year.empty:
        flag_counts = country_year["audit_flag"].value_counts(dropna=False).to_dict()
        country_year_note = ";".join(f"{key}={value}" for key, value in sorted(flag_counts.items()))

    validation_decision = lookup_decision(decision_rows, "genomics_decision_4_validation")
    repeated_decision = lookup_decision(decision_rows, "genomics_decision_1_repeated_origin")
    ecology_decision = lookup_decision(decision_rows, "ecology_decision_1_panel_size")
    usa_decision = lookup_decision(decision_rows, "usa_decision")

    rows = [
        {
            "claim_boundary_id": "frozen_vs_live_state",
            "claim_or_denominator_concern": "Live workflow manifests can differ from the frozen submission-facing manuscript frame after later development refreshes.",
            "resolution_status": "resolved_by_release_freeze_contract",
            "evidence_source": "manuscript/submission_data/cohort/master_cohort_flow_summary.tsv;state/manifest/manifest_build_report.json",
            "diagnostic_value": build_frozen_live_frame_note(),
            "manuscript_control": "Treat Supplementary Tables, Source Data and submission-facing ledgers as the frozen manuscript contract; state/manifest reports are post-freeze development diagnostics unless the submission package is regenerated.",
            "residual_boundary": "Any future live-state update must regenerate all manuscript tables, source-data workbooks and text before changing manuscript denominators.",
            "minimum_safe_claim": "All manuscript numbers refer to the frozen Communications Biology analysis frame, not incidental post-freeze live manifest summaries.",
        },
        {
            "claim_boundary_id": "archive_ascertainment",
            "claim_or_denominator_concern": "Public genomes may not estimate circulation-wide prevalence.",
            "resolution_status": "bounded_in_claim",
            "evidence_source": "Supplementary Fig. 11;Supplementary Fig. 12;Source Data",
            "diagnostic_value": country_year_note,
            "manuscript_control": "Use recoverable-locus disruption in observed public genomes as the principal estimand.",
            "residual_boundary": "Population prevalence requires denser catchment-defined sampling.",
            "minimum_safe_claim": "Archive-level recoverable-locus disruption, not circulation-wide prevalence.",
        },
        {
            "claim_boundary_id": "study_block_amplification",
            "claim_or_denominator_concern": "A dominant accession block could inflate the apparent leading event.",
            "resolution_status": "stress_tested_bounded",
            "evidence_source": "Supplementary Fig. 14;Source Data",
            "diagnostic_value": study_block_note,
            "manuscript_control": "Report attenuation under equal block weighting and keep study-weighted dominance as a stress test rather than a primary event-ranking estimator.",
            "residual_boundary": "Study uploads still amplify event burden and should not be treated as unbiased population counts.",
            "minimum_safe_claim": "Structural reuse persists after block-weighting stress tests, but public uploads amplify it.",
        },
        {
            "claim_boundary_id": "repeated_origin_specificity",
            "claim_or_denominator_concern": "Repeated origins could be an artefact of tree frame choice.",
            "resolution_status": clean_text(repeated_decision.get("status", "not_available")),
            "evidence_source": "Fig. 3;Supplementary Table 6;Supplementary Figs. 6 and 14;Source Data",
            "diagnostic_value": clean_text(repeated_decision.get("observed_value", "")),
            "manuscript_control": "Interpret origin counts as scenario-level evidence against a one-clone explanation.",
            "residual_boundary": "Exact historical transition counts require denser sampling and direct transmission context.",
            "minimum_safe_claim": "Repeated emergence is supported qualitatively; exact counts are not historical truth claims.",
        },
        {
            "claim_boundary_id": "orthogonal_validation",
            "claim_or_denominator_concern": "Assembly-side event calls need independent structural support.",
            "resolution_status": clean_text(validation_decision.get("status", "not_available")),
            "evidence_source": "Supplementary Table 8;Supplementary Fig. 14;Source Data",
            "diagnostic_value": clean_text(validation_decision.get("observed_value", "")),
            "manuscript_control": "Keep read-backed and long-read/hybrid anchors separate from assembly-only event calls.",
            "residual_boundary": "Low-frequency or assembly-only classes remain lower-confidence until further read-backed validation.",
            "minimum_safe_claim": "Major recurrent architectures are anchored; rare classes remain bounded genome-defined disruptions.",
        },
        {
            "claim_boundary_id": "protein_expression_boundary",
            "claim_or_denominator_concern": "Genome-defined disruption could be mistaken for direct PRN protein non-expression.",
            "resolution_status": "tiered_not_genome_by_genome",
            "evidence_source": "Fig. 5c;Supplementary Table 10;Supplementary Fig. 15;Source Data",
            "diagnostic_value": phenotype_note,
            "manuscript_control": "Use phenotype tiers and state that protein expression is a separate measurement layer.",
            "residual_boundary": "Universal PRN non-expression requires paired expression assays.",
            "minimum_safe_claim": "Genome-defined disruption with lesion-class phenotype plausibility, not universal protein proof.",
        },
        {
            "claim_boundary_id": "comparator_specificity",
            "claim_or_denominator_concern": "The *prn* signal could reflect generic marker loss or assembly fragmentation.",
            "resolution_status": "specificity_audited",
            "evidence_source": "Fig. 5b;Supplementary Fig. 10;Source Data",
            "diagnostic_value": specificity_note,
            "manuscript_control": "Use matched marker-recovery logic across *prn* and comparator loci.",
            "residual_boundary": "The audit addresses gross locus loss, not every possible event-specific technical artefact.",
            "minimum_safe_claim": "The *prn* signal is much stronger than comparator marker-loss proxies.",
        },
        {
            "claim_boundary_id": "programme_causality",
            "claim_or_denominator_concern": "Country-programme contrasts could be read as causal vaccine-effect estimates.",
            "resolution_status": clean_text(ecology_decision.get("status", "not_available")),
            "evidence_source": "Fig. 4;Supplementary Figs. 7, 8, 11 and 12;Source Data",
            "diagnostic_value": clean_text(ecology_decision.get("observed_value", "")),
            "manuscript_control": "Keep country-programme contrasts descriptive and frame models as sensitivity/context only.",
            "residual_boundary": "Product-level exposure metadata and catchment-defined sampling are needed for causal inference.",
            "minimum_safe_claim": "Country-programme settings amplify observable archive patterns unevenly.",
        },
        {
            "claim_boundary_id": "usa_dynamics",
            "claim_or_denominator_concern": "Dynamic incidence models may overstate transmission relevance.",
            "resolution_status": clean_text(usa_decision.get("status", "not_available")),
            "evidence_source": "Supplementary Fig. 9;Source Data",
            "diagnostic_value": clean_text(usa_decision.get("observed_value", "")),
            "manuscript_control": "Keep USA dynamics supplementary and do not make it a main-text pillar.",
            "residual_boundary": "Transmission dynamics need denser epidemiologic and genomic overlap.",
            "minimum_safe_claim": "Dynamic analyses are readiness/context audits, not core causal evidence.",
        },
    ]
    return rows


def main() -> None:
    disrupted_calls = load_disrupted_mechanism_calls()
    disrupted_calls_with_markers = load_disrupted_calls_with_markers()
    lineage_rows = build_lineage_collapsed_event_table(disrupted_calls)
    collapse_sensitivity_rows = build_lineage_collapse_sensitivity_table(disrupted_calls_with_markers)
    origin_rows = build_origin_collapsed_event_table()
    validation_rows = build_architecture_origin_validation_matrix(origin_rows)

    formulation_manifest = build_formulation_manifest()
    explicit_ge3 = formulation_manifest[formulation_manifest["eligible_component_model_ge3"]].copy()
    ecology_robustness, ecology_summary = fit_component_clustered_fractional_logit(explicit_ge3)
    usa_lag_df = build_usa_lag_sensitivity()
    decision_rows = build_submission_summary(
        lineage_rows=lineage_rows,
        collapse_sensitivity_rows=collapse_sensitivity_rows,
        origin_rows=origin_rows,
        validation_rows=validation_rows,
        formulation_manifest=formulation_manifest,
        ecology_robustness=ecology_robustness,
        ecology_summary=ecology_summary,
        usa_lag_df=usa_lag_df,
    )
    claim_ledger_rows = build_claim_evidence_denominator_ledger(decision_rows)

    write_tsv(LINEAGE_COLLAPSED_OUT, lineage_rows)
    write_tsv(LINEAGE_COLLAPSE_SENSITIVITY_OUT, collapse_sensitivity_rows)
    write_tsv(ORIGIN_COLLAPSED_OUT, origin_rows)
    write_tsv(ARCHITECTURE_VALIDATION_OUT, validation_rows)
    write_tsv(FORMULATION_MANIFEST_OUT, formulation_manifest.to_dict("records"))
    write_tsv(ECOLOGY_ROBUSTNESS_OUT, ecology_robustness.to_dict("records"))
    write_tsv(USA_LAG_OUT, usa_lag_df.to_dict("records"))
    write_tsv(SUMMARY_OUT, decision_rows)
    write_tsv(CLAIM_EVIDENCE_DENOMINATOR_LEDGER_OUT, claim_ledger_rows)

    write_tsv(SUPP25, lineage_rows)
    write_tsv(SUPP26, origin_rows)
    write_tsv(SUPP27, validation_rows)
    write_tsv(SUPP28, formulation_manifest.to_dict("records"))
    write_tsv(SUPP29, ecology_robustness.to_dict("records"))
    write_tsv(SUPP30, usa_lag_df.to_dict("records"))
    write_tsv(AUDIT_LEDGER_DIR / "Supplementary_Table_30_USA_Lag_Sensitivity.tsv", usa_lag_df.to_dict("records"))
    write_tsv(SUPP31, collapse_sensitivity_rows)


if __name__ == "__main__":
    main()
