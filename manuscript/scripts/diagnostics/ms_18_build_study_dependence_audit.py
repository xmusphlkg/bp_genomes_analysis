#!/usr/bin/env python3
"""Build manuscript-facing study-dependence audit sidecars.

This audit addresses pseudo-replication concerns in the public-genome panel by:
1. Assigning BioProject/study/sample blocks to the retained manifest.
2. Re-running selected-country contrasts with block-aware summaries.
3. Re-weighting structural-event concentration across study blocks.
4. Staging study-block-balanced ASR resampling summaries.

The canonical workflow manifest is left unchanged. All outputs are manuscript-
facing sidecars derived from the retained manifest plus the archival inventory.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
FIGURE_DATA_DIR = REPO_ROOT / "manuscript" / "figure_data"
STUDY_DIR = FIGURE_DATA_DIR / "study_dependence"
SELECTED_DIR = FIGURE_DATA_DIR / "selected_country"
SUPP_DIR = REPO_ROOT / "manuscript" / "supplementary"
ARCHIVE_DIR = REPO_ROOT / "manuscript" / "supplementary" / "_archive"

MANIFEST_PATH = REPO_ROOT / "outputs" / "workflow" / "manifest" / "manifest.tsv"
INVENTORY_PATH = ARCHIVE_DIR / "supplementary_table_s1_data_inventory.tsv"
HISTORY_PATH = SELECTED_DIR / "country_program_history_manifest.tsv"
CONTRAST_PATH = SELECTED_DIR / "country_epoch_contrast_summary.tsv"
STRUCTURE_REUSE_PATH = SELECTED_DIR / "selected_country_structure_reuse.tsv"
STRUCTURAL_CONCENTRATION_PATH = FIGURE_DATA_DIR / "structural_event_concentration.tsv"
ASR_TIP_STATES_PATH = REPO_ROOT / "outputs" / "workflow" / "asr" / "tip_states.tsv"
ASR_RESAMPLING_REPLICATE_PATH = REPO_ROOT / "outputs" / "workflow" / "asr_resampling" / "resampling_replicates.tsv"
ASR_RESAMPLING_SUMMARY_PATH = REPO_ROOT / "outputs" / "workflow" / "asr_resampling" / "resampling_summary.tsv"

BLOCK_ASSIGNMENT_PATH = STUDY_DIR / "study_block_assignment.tsv"
LOO_PATH = STUDY_DIR / "selected_country_block_leave_one_out.tsv"
BOOTSTRAP_PATH = STUDY_DIR / "selected_country_block_bootstrap.tsv"
DOMINANCE_PATH = STUDY_DIR / "selected_country_block_dominance.tsv"
STRUCTURE_REWEIGHTED_PATH = STUDY_DIR / "structure_reuse_block_reweighted.tsv"
ASR_STUDY_BLOCK_RESAMPLING_PATH = STUDY_DIR / "asr_study_block_resampling.tsv"

SUPP52_PATH = SUPP_DIR / "Supplementary_Table_52_Study_Block_Assignment_and_Dominance.tsv"
SUPP53_PATH = SUPP_DIR / "Supplementary_Table_53_Selected_Country_Block_Audit.tsv"
SUPP54_PATH = SUPP_DIR / "Supplementary_Table_54_Study_Weighted_Structure_and_ASR.tsv"

N_BOOTSTRAP = 1000
BOOTSTRAP_SEED = 20260419
SELECTED_COUNTRIES = ("USA", "NZL", "AUS", "JPN")
PRIMARY_COUNTRY_PAIRS = {
    "USA": ("usa_wp_only", "usa_ap_prn_background"),
    "NZL": ("nzl_wp_only", "nzl_ap_with_prn"),
    "AUS": ("aus_wp_only", "aus_ap_with_prn"),
    "JPN": ("jpn_pre2012_mixed_ap", "jpn_ap_without_prn"),
}


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def clean_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = str(value).strip()
    return "" if text.casefold() in {"", "nan", "none", "na", "not available"} else text


def normalize_text_series(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).map(clean_text)


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


def safe_divide(numerator: object, denominator: object) -> float:
    numerator_value = pd.to_numeric(pd.Series([numerator]), errors="coerce").iloc[0]
    denominator_value = pd.to_numeric(pd.Series([denominator]), errors="coerce").iloc[0]
    if pd.isna(numerator_value) or pd.isna(denominator_value) or denominator_value == 0:
        return np.nan
    return float(numerator_value / denominator_value)


def sign_label(value: object) -> str:
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(numeric):
        return "non_estimable"
    if numeric > 0:
        return "increase"
    if numeric < 0:
        return "decrease"
    return "no_change"


def first_nonempty(values: Iterable[object]) -> str:
    for value in values:
        text = clean_text(value)
        if text:
            return text
    return ""


def year_text_from_numeric(value: object) -> str:
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(numeric):
        return ""
    return str(int(numeric))


def prepare_manifest() -> pd.DataFrame:
    manifest = read_tsv(MANIFEST_PATH)
    manifest["sample_id_canonical"] = normalize_text_series(manifest["sample_id_canonical"])
    manifest["country_iso3"] = normalize_text_series(manifest["country_iso3"])
    manifest["country_program_target"] = normalize_text_series(manifest["country_program_target"])
    manifest["bioproject_accession_manifest"] = normalize_text_series(manifest.get("bioproject_accession", pd.Series("", index=manifest.index)))
    manifest["year_num"] = to_numeric(manifest.get("year", pd.Series("", index=manifest.index)))
    manifest["prn_interpretable_bool"] = to_bool(manifest["prn_interpretable"])
    manifest["prn_disrupted_bool"] = to_bool(manifest["prn_disrupted"])
    manifest["has_reads_bool"] = to_bool(manifest["has_reads"])
    manifest["prn_event_id"] = normalize_text_series(manifest["prn_event_id"])
    manifest["prn_mechanism_call"] = normalize_text_series(manifest["prn_mechanism_call"])
    return manifest


def prepare_inventory() -> pd.DataFrame:
    inventory = read_tsv(INVENTORY_PATH)
    keep_columns = [
        "sample_id_canonical",
        "biosample_accession",
        "bioproject_accession",
        "study_accession",
        "country",
        "year",
        "month",
        "week_key",
        "collection_date_raw",
        "isolate",
        "strain",
    ]
    inventory = inventory.loc[:, [column for column in keep_columns if column in inventory.columns]].copy()
    for column in inventory.columns:
        inventory[column] = normalize_text_series(inventory[column])
    inventory = inventory.loc[inventory["sample_id_canonical"].ne("")].copy()
    aggregated = (
        inventory.groupby("sample_id_canonical", dropna=False)
        .agg({column: first_nonempty for column in inventory.columns if column != "sample_id_canonical"})
        .reset_index()
    )
    if "year" in aggregated.columns:
        aggregated = aggregated.rename(columns={"year": "year_inventory"})
    return aggregated


def build_epoch_lookup() -> pd.DataFrame:
    history = read_tsv(HISTORY_PATH)
    history["country_iso3"] = normalize_text_series(history["country_iso3"])
    history["epoch_id"] = normalize_text_series(history["epoch_id"])
    history["epoch_label"] = normalize_text_series(history["epoch_label"])
    history["start_year_num"] = to_numeric(history["start_year"])
    history["end_year_num"] = to_numeric(history["end_year"])
    return history.loc[history["country_iso3"].isin(SELECTED_COUNTRIES)].copy()


def first_row(frame: pd.DataFrame, columns: list[str]) -> dict[str, object]:
    if frame.empty:
        return {column: "" for column in columns}
    row = frame.iloc[0]
    return {column: row.get(column, "") for column in columns}


def choose_base_block(row: pd.Series) -> tuple[str, str]:
    bioproject = clean_text(row.get("bioproject_accession"))
    study = clean_text(row.get("study_accession"))
    sample_id = clean_text(row.get("sample_id_canonical"))
    if bioproject:
        return ("bioproject_accession", bioproject)
    if study:
        return ("study_accession", study)
    return ("sample_id_canonical_singleton", f"singleton:{sample_id}")


def choose_subblock(row: pd.Series) -> tuple[str, str]:
    base_block_id = clean_text(row.get("base_block_id"))
    week_key = clean_text(row.get("week_key"))
    month = clean_text(row.get("month"))
    year_text = clean_text(row.get("year"))
    if week_key:
        return ("base_block_plus_week_key", f"{base_block_id}::week={week_key}")
    if month:
        return ("base_block_plus_month", f"{base_block_id}::month={month}")
    if year_text:
        return ("base_block_plus_year", f"{base_block_id}::year={year_text}")
    return ("base_block_only", base_block_id)


def assign_selected_country_epochs(frame: pd.DataFrame, history: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    output["selected_country_epoch_id"] = ""
    output["selected_country_epoch_label"] = ""
    for item in history.itertuples(index=False):
        mask = (
            output["country_iso3"].eq(item.country_iso3)
            & output["year_num"].ge(item.start_year_num)
            & output["year_num"].le(item.end_year_num)
        )
        output.loc[mask, "selected_country_epoch_id"] = item.epoch_id
        output.loc[mask, "selected_country_epoch_label"] = item.epoch_label
    return output


def build_block_assignment(manifest: pd.DataFrame, inventory: pd.DataFrame, history: pd.DataFrame) -> pd.DataFrame:
    merged = manifest.merge(
        inventory,
        on="sample_id_canonical",
        how="left",
        suffixes=("", "_inventory"),
    )
    merged["bioproject_accession"] = normalize_text_series(
        merged.get("bioproject_accession", pd.Series("", index=merged.index))
    )
    merged["study_accession"] = normalize_text_series(
        merged.get("study_accession", pd.Series("", index=merged.index))
    )
    merged["year_inventory_num"] = to_numeric(merged.get("year_inventory", pd.Series("", index=merged.index)))
    merged["year_num"] = merged["year_num"].where(merged["year_num"].notna(), merged["year_inventory_num"])
    merged["year"] = merged["year_num"].map(year_text_from_numeric)
    merged["year"] = merged["year"].where(
        merged["year"].ne(""),
        normalize_text_series(merged.get("year_inventory", pd.Series("", index=merged.index))),
    )
    base_assignment = merged.apply(choose_base_block, axis=1, result_type="expand")
    merged["base_block_level"] = base_assignment[0]
    merged["base_block_id"] = base_assignment[1]
    subblock_assignment = merged.apply(choose_subblock, axis=1, result_type="expand")
    merged["subblock_level"] = subblock_assignment[0]
    merged["subblock_id"] = subblock_assignment[1]
    merged["base_block_size"] = merged.groupby("base_block_id", dropna=False)["sample_id_canonical"].transform("size")
    merged["subblock_size"] = merged.groupby("subblock_id", dropna=False)["sample_id_canonical"].transform("size")
    block_priority = {
        "bioproject_accession": 0,
        "study_accession": 1,
        "sample_id_canonical_singleton": 2,
    }
    merged["base_block_level_rank"] = merged["base_block_level"].map(block_priority).fillna(99)
    canonical_block_levels = (
        merged.sort_values(["base_block_id", "base_block_level_rank"])
        .drop_duplicates(subset=["base_block_id"])
        .loc[:, ["base_block_id", "base_block_level"]]
        .rename(columns={"base_block_level": "base_block_level_canonical"})
    )
    merged = merged.merge(canonical_block_levels, on="base_block_id", how="left")
    merged["base_block_level"] = merged["base_block_level_canonical"]
    merged = merged.drop(columns=["base_block_level_canonical", "base_block_level_rank"])
    merged["dominant_block_flag"] = False
    merged["country_iso3_key"] = merged["country_iso3"].replace("", "unknown_country")
    merged["target_key"] = merged["country_program_target"].replace("", "unknown_target")
    grouping = (
        merged.groupby(["country_iso3_key", "target_key", "base_block_id"], dropna=False)
        .size()
        .reset_index(name="block_n")
    )
    max_by_context = grouping.groupby(["country_iso3_key", "target_key"], dropna=False)["block_n"].transform("max")
    dominant_lookup = grouping.loc[grouping["block_n"].eq(max_by_context), ["country_iso3_key", "target_key", "base_block_id"]]
    dominant_lookup["dominant_block_flag"] = True
    merged = merged.merge(
        dominant_lookup,
        on=["country_iso3_key", "target_key", "base_block_id"],
        how="left",
        suffixes=("", "_lookup"),
    )
    merged["dominant_block_flag"] = merged["dominant_block_flag_lookup"].eq(True)
    merged = merged.drop(columns=["country_iso3_key", "target_key", "dominant_block_flag_lookup"])
    merged = assign_selected_country_epochs(merged, history)

    assignment = merged[
        [
            "sample_id_canonical",
            "country_iso3",
            "country_program_target",
            "selected_country_epoch_id",
            "selected_country_epoch_label",
            "base_block_id",
            "base_block_level",
            "subblock_id",
            "subblock_level",
            "base_block_size",
            "subblock_size",
            "dominant_block_flag",
            "bioproject_accession",
            "study_accession",
            "year",
            "month",
            "week_key",
            "collection_date_raw",
            "prn_interpretable",
            "prn_disrupted",
            "prn_event_id",
        ]
    ].copy()
    assignment = assignment.rename(
        columns={
            "selected_country_epoch_id": "epoch_id",
            "selected_country_epoch_label": "epoch_label",
        }
    )
    return assignment, merged


def summarize_epoch(frame: pd.DataFrame) -> dict[str, object]:
    total_n = int(len(frame))
    interpretable = frame.loc[frame["prn_interpretable_bool"]].copy()
    n_interpretable = int(len(interpretable))
    n_disrupted = int(interpretable["prn_disrupted_bool"].sum())
    naive_prevalence = safe_divide(n_disrupted, n_interpretable)
    if interpretable.empty:
        block_equalized = np.nan
        n_interpretable_blocks = 0
    else:
        block_prevalence = interpretable.groupby("base_block_id", dropna=False)["prn_disrupted_bool"].mean()
        block_equalized = float(block_prevalence.mean()) if len(block_prevalence) else np.nan
        n_interpretable_blocks = int(len(block_prevalence))
    return {
        "n_total_samples": total_n,
        "n_interpretable": n_interpretable,
        "n_disrupted": n_disrupted,
        "naive_prevalence": naive_prevalence,
        "block_equalized_prevalence": block_equalized,
        "n_total_blocks": int(frame["base_block_id"].nunique()),
        "n_interpretable_blocks": n_interpretable_blocks,
        "n_subblocks": int(frame["subblock_id"].nunique()),
    }


def compare_epochs(previous_frame: pd.DataFrame, next_frame: pd.DataFrame) -> dict[str, object]:
    previous = summarize_epoch(previous_frame)
    next_epoch = summarize_epoch(next_frame)
    return {
        **{f"previous_{key}": value for key, value in previous.items()},
        **{f"next_{key}": value for key, value in next_epoch.items()},
        "delta_naive_prevalence": (
            next_epoch["naive_prevalence"] - previous["naive_prevalence"]
            if not pd.isna(previous["naive_prevalence"]) and not pd.isna(next_epoch["naive_prevalence"])
            else np.nan
        ),
        "delta_block_equalized_prevalence": (
            next_epoch["block_equalized_prevalence"] - previous["block_equalized_prevalence"]
            if not pd.isna(previous["block_equalized_prevalence"]) and not pd.isna(next_epoch["block_equalized_prevalence"])
            else np.nan
        ),
    }


def resample_by_subblock(frame: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    subblocks = sorted(frame["subblock_id"].dropna().astype(str).unique().tolist())
    if not subblocks:
        return frame.iloc[0:0].copy()
    sampled = rng.choice(subblocks, size=len(subblocks), replace=True)
    chunks = [frame.loc[frame["subblock_id"].eq(subblock)].copy() for subblock in sampled]
    return pd.concat(chunks, ignore_index=True) if chunks else frame.iloc[0:0].copy()


def summarize_bootstrap_metric(values: list[float], observed: float) -> dict[str, object]:
    numeric = pd.Series(values, dtype=float)
    estimable = numeric.dropna()
    out: dict[str, object] = {
        "median": np.nan,
        "lower_95": np.nan,
        "upper_95": np.nan,
        "sign_flip_fraction": np.nan,
        "non_estimable_fraction": safe_divide(int(numeric.isna().sum()), int(len(numeric))),
    }
    if estimable.empty:
        return out
    out["median"] = float(estimable.median())
    out["lower_95"] = float(estimable.quantile(0.025))
    out["upper_95"] = float(estimable.quantile(0.975))
    observed_sign = sign_label(observed)
    if observed_sign in {"increase", "decrease"}:
        sign_flips = estimable.map(sign_label).ne(observed_sign).sum()
        out["sign_flip_fraction"] = safe_divide(int(sign_flips), int(len(numeric)))
    return out


def build_selected_country_block_audit(study_frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    contrasts = read_tsv(CONTRAST_PATH)
    eligible = contrasts.loc[
        normalize_text_series(contrasts["contrast_eligible"]).str.casefold().eq("true")
    ].copy()
    history = build_epoch_lookup()
    history_lookup = {
        (row.country_iso3, row.epoch_id): row
        for row in history.itertuples(index=False)
    }
    contrast_keys = set(zip(eligible["country_iso3"], eligible["previous_epoch_id"], eligible["next_epoch_id"]))
    supplemental_rows: list[dict[str, object]] = []
    for country_iso3, (previous_epoch_id, next_epoch_id) in PRIMARY_COUNTRY_PAIRS.items():
        key = (country_iso3, previous_epoch_id, next_epoch_id)
        if key in contrast_keys:
            continue
        previous_epoch = history_lookup.get((country_iso3, previous_epoch_id))
        next_epoch = history_lookup.get((country_iso3, next_epoch_id))
        if previous_epoch is None or next_epoch is None:
            continue
        contrast_rows = contrasts.loc[contrasts["country_iso3"].eq(country_iso3)]
        country_name = first_nonempty(contrast_rows["country_name"]) or clean_text(previous_epoch.country_name)
        supplemental_rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": country_name,
                "previous_epoch_id": previous_epoch_id,
                "next_epoch_id": next_epoch_id,
                "previous_epoch_label": clean_text(previous_epoch.epoch_label),
                "next_epoch_label": clean_text(next_epoch.epoch_label),
                "contrast_type": "primary_country_pair_nonadjacent",
                "contrast_eligible": "True",
                "notes": "Primary focal-country block audit assembled directly from the selected-country epoch manifest.",
            }
        )
    if supplemental_rows:
        eligible = pd.concat([eligible, pd.DataFrame(supplemental_rows)], ignore_index=True, sort=False)

    leave_one_out_rows: list[dict[str, object]] = []
    bootstrap_rows: list[dict[str, object]] = []
    dominance_rows: list[dict[str, object]] = []

    for contrast in eligible.itertuples(index=False):
        country_iso3 = contrast.country_iso3
        previous_epoch_id = contrast.previous_epoch_id
        next_epoch_id = contrast.next_epoch_id
        pair_label = f"{country_iso3}:{previous_epoch_id}->{next_epoch_id}"
        previous_frame = study_frame.loc[
            study_frame["country_iso3"].eq(country_iso3)
            & study_frame["selected_country_epoch_id"].eq(previous_epoch_id)
        ].copy()
        next_frame = study_frame.loc[
            study_frame["country_iso3"].eq(country_iso3)
            & study_frame["selected_country_epoch_id"].eq(next_epoch_id)
        ].copy()
        observed = compare_epochs(previous_frame, next_frame)
        observed_row = {
            "row_type": "observed",
            "comparison_id": pair_label,
            "country_iso3": country_iso3,
            "country_name": contrast.country_name,
            "previous_epoch_id": previous_epoch_id,
            "next_epoch_id": next_epoch_id,
            "previous_epoch_label": contrast.previous_epoch_label,
            "next_epoch_label": contrast.next_epoch_label,
            "dropped_block_id": "",
            "dropped_block_level": "",
            "dropped_block_size_combined": np.nan,
            "dropped_block_is_dominant_previous_epoch": False,
            "dropped_block_is_dominant_next_epoch": False,
            "dropped_block_is_dominant_any": False,
            **observed,
            "delta_naive_sign": sign_label(observed["delta_naive_prevalence"]),
            "delta_block_equalized_sign": sign_label(observed["delta_block_equalized_prevalence"]),
            "direction_matches_observed_naive": True,
            "direction_matches_observed_block_equalized": True,
        }
        leave_one_out_rows.append(observed_row)

        for epoch_slot, epoch_id, epoch_label, epoch_frame in [
            ("previous", previous_epoch_id, contrast.previous_epoch_label, previous_frame),
            ("next", next_epoch_id, contrast.next_epoch_label, next_frame),
        ]:
            block_rows = []
            for block_id, block_frame in epoch_frame.groupby("base_block_id", dropna=False):
                interpretable = block_frame.loc[block_frame["prn_interpretable_bool"]]
                block_rows.append(
                    {
                        "comparison_id": pair_label,
                        "country_iso3": country_iso3,
                        "country_name": contrast.country_name,
                        "epoch_slot": epoch_slot,
                        "epoch_id": epoch_id,
                        "epoch_label": epoch_label,
                        "block_id": block_id,
                        "block_level": first_nonempty(block_frame["base_block_level"]),
                        "block_total_samples": int(len(block_frame)),
                        "block_interpretable_samples": int(len(interpretable)),
                        "block_disrupted_samples": int(interpretable["prn_disrupted_bool"].sum()),
                        "block_total_share": safe_divide(int(len(block_frame)), int(len(epoch_frame))),
                        "block_interpretable_share": safe_divide(int(len(interpretable)), int(epoch_frame["prn_interpretable_bool"].sum())),
                        "block_disrupted_share": safe_divide(
                            int(interpretable["prn_disrupted_bool"].sum()),
                            int(epoch_frame.loc[epoch_frame["prn_interpretable_bool"], "prn_disrupted_bool"].sum()),
                        ),
                        "block_is_dominant": False,
                    }
                )
            if not block_rows:
                continue
            block_table = pd.DataFrame(block_rows).sort_values(
                ["block_total_samples", "block_interpretable_samples", "block_id"],
                ascending=[False, False, True],
            )
            if not block_table.empty:
                max_total = block_table["block_total_samples"].max()
                block_table["block_is_dominant"] = block_table["block_total_samples"].eq(max_total)
                dominance_rows.extend(block_table.to_dict(orient="records"))

        combined = pd.concat([previous_frame, next_frame], ignore_index=True)
        if combined.empty:
            continue
        for block_id, dropped_frame in combined.groupby("base_block_id", dropna=False):
            dropped_level = first_nonempty(dropped_frame["base_block_level"])
            previous_dominant = bool(
                pd.Series(dominance_rows)
                .map(lambda row: row.get("comparison_id") == pair_label and row.get("epoch_slot") == "previous" and row.get("block_id") == block_id and bool(row.get("block_is_dominant")))
                .any()
            )
            next_dominant = bool(
                pd.Series(dominance_rows)
                .map(lambda row: row.get("comparison_id") == pair_label and row.get("epoch_slot") == "next" and row.get("block_id") == block_id and bool(row.get("block_is_dominant")))
                .any()
            )
            leave_previous = previous_frame.loc[previous_frame["base_block_id"].ne(block_id)].copy()
            leave_next = next_frame.loc[next_frame["base_block_id"].ne(block_id)].copy()
            leave_summary = compare_epochs(leave_previous, leave_next)
            leave_one_out_rows.append(
                {
                    "row_type": "leave_one_out",
                    "comparison_id": pair_label,
                    "country_iso3": country_iso3,
                    "country_name": contrast.country_name,
                    "previous_epoch_id": previous_epoch_id,
                    "next_epoch_id": next_epoch_id,
                    "previous_epoch_label": contrast.previous_epoch_label,
                    "next_epoch_label": contrast.next_epoch_label,
                    "dropped_block_id": block_id,
                    "dropped_block_level": dropped_level,
                    "dropped_block_size_combined": int(len(dropped_frame)),
                    "dropped_block_is_dominant_previous_epoch": previous_dominant,
                    "dropped_block_is_dominant_next_epoch": next_dominant,
                    "dropped_block_is_dominant_any": previous_dominant or next_dominant,
                    **leave_summary,
                    "delta_naive_sign": sign_label(leave_summary["delta_naive_prevalence"]),
                    "delta_block_equalized_sign": sign_label(leave_summary["delta_block_equalized_prevalence"]),
                    "direction_matches_observed_naive": sign_label(leave_summary["delta_naive_prevalence"]) == observed_row["delta_naive_sign"],
                    "direction_matches_observed_block_equalized": (
                        sign_label(leave_summary["delta_block_equalized_prevalence"]) == observed_row["delta_block_equalized_sign"]
                    ),
                }
            )

        rng = np.random.default_rng(BOOTSTRAP_SEED + sum(ord(character) for character in pair_label))
        bootstrap_naive: list[float] = []
        bootstrap_block_equalized: list[float] = []
        previous_naive_values: list[float] = []
        next_naive_values: list[float] = []
        previous_block_values: list[float] = []
        next_block_values: list[float] = []
        for _ in range(N_BOOTSTRAP):
            previous_boot = resample_by_subblock(previous_frame, rng)
            next_boot = resample_by_subblock(next_frame, rng)
            boot_summary = compare_epochs(previous_boot, next_boot)
            bootstrap_naive.append(boot_summary["delta_naive_prevalence"])
            bootstrap_block_equalized.append(boot_summary["delta_block_equalized_prevalence"])
            previous_naive_values.append(boot_summary["previous_naive_prevalence"])
            next_naive_values.append(boot_summary["next_naive_prevalence"])
            previous_block_values.append(boot_summary["previous_block_equalized_prevalence"])
            next_block_values.append(boot_summary["next_block_equalized_prevalence"])

        bootstrap_naive_summary = summarize_bootstrap_metric(bootstrap_naive, observed_row["delta_naive_prevalence"])
        bootstrap_block_summary = summarize_bootstrap_metric(
            bootstrap_block_equalized, observed_row["delta_block_equalized_prevalence"]
        )
        dominant_candidates = [
            row
            for row in leave_one_out_rows
            if row["comparison_id"] == pair_label and row["row_type"] == "leave_one_out" and row["dropped_block_is_dominant_any"]
        ]
        dominant_row = max(dominant_candidates, key=lambda row: (row["dropped_block_size_combined"], row["dropped_block_id"])) if dominant_candidates else None
        dominant_naive = np.nan if dominant_row is None else dominant_row["delta_naive_prevalence"]
        dominant_block = np.nan if dominant_row is None else dominant_row["delta_block_equalized_prevalence"]
        dominant_non_estimable = dominant_row is None or sign_label(dominant_naive) == "non_estimable"
        dominant_reversal = dominant_row is not None and sign_label(dominant_naive) not in {
            observed_row["delta_naive_sign"],
            "no_change",
            "non_estimable",
        }
        bootstrap_cross_zero = (
            pd.notna(bootstrap_naive_summary["lower_95"])
            and pd.notna(bootstrap_naive_summary["upper_95"])
            and bootstrap_naive_summary["lower_95"] <= 0 <= bootstrap_naive_summary["upper_95"]
        )
        dominant_collapse = (
            dominant_row is not None
            and pd.notna(dominant_naive)
            and pd.notna(observed_row["delta_naive_prevalence"])
            and abs(dominant_naive) < 0.5 * abs(observed_row["delta_naive_prevalence"])
        )
        if dominant_non_estimable or dominant_reversal:
            study_call = "study_dominated_or_reversed"
        elif bootstrap_cross_zero or dominant_collapse:
            study_call = (
                "study_sensitive_positive_signal"
                if observed_row["delta_naive_sign"] == "increase"
                else "study_sensitive_directional_signal"
            )
        else:
            study_call = "study_dependence_audit_retained_direction"

        bootstrap_rows.append(
            {
                "comparison_id": pair_label,
                "country_iso3": country_iso3,
                "country_name": contrast.country_name,
                "previous_epoch_id": previous_epoch_id,
                "next_epoch_id": next_epoch_id,
                "previous_epoch_label": contrast.previous_epoch_label,
                "next_epoch_label": contrast.next_epoch_label,
                "n_bootstrap_replicates": N_BOOTSTRAP,
                "previous_subblock_count": int(previous_frame["subblock_id"].nunique()),
                "next_subblock_count": int(next_frame["subblock_id"].nunique()),
                "previous_base_block_count": int(previous_frame["base_block_id"].nunique()),
                "next_base_block_count": int(next_frame["base_block_id"].nunique()),
                "observed_previous_naive_prevalence": observed_row["previous_naive_prevalence"],
                "observed_next_naive_prevalence": observed_row["next_naive_prevalence"],
                "observed_delta_naive_prevalence": observed_row["delta_naive_prevalence"],
                "bootstrap_previous_naive_prevalence_median": summarize_bootstrap_metric(
                    previous_naive_values, observed_row["previous_naive_prevalence"]
                )["median"],
                "bootstrap_next_naive_prevalence_median": summarize_bootstrap_metric(
                    next_naive_values, observed_row["next_naive_prevalence"]
                )["median"],
                "bootstrap_delta_naive_prevalence_median": bootstrap_naive_summary["median"],
                "bootstrap_delta_naive_prevalence_lower_95": bootstrap_naive_summary["lower_95"],
                "bootstrap_delta_naive_prevalence_upper_95": bootstrap_naive_summary["upper_95"],
                "bootstrap_delta_naive_sign_flip_fraction": bootstrap_naive_summary["sign_flip_fraction"],
                "bootstrap_delta_naive_non_estimable_fraction": bootstrap_naive_summary["non_estimable_fraction"],
                "observed_previous_block_equalized_prevalence": observed_row["previous_block_equalized_prevalence"],
                "observed_next_block_equalized_prevalence": observed_row["next_block_equalized_prevalence"],
                "observed_delta_block_equalized_prevalence": observed_row["delta_block_equalized_prevalence"],
                "bootstrap_previous_block_equalized_prevalence_median": summarize_bootstrap_metric(
                    previous_block_values, observed_row["previous_block_equalized_prevalence"]
                )["median"],
                "bootstrap_next_block_equalized_prevalence_median": summarize_bootstrap_metric(
                    next_block_values, observed_row["next_block_equalized_prevalence"]
                )["median"],
                "bootstrap_delta_block_equalized_prevalence_median": bootstrap_block_summary["median"],
                "bootstrap_delta_block_equalized_prevalence_lower_95": bootstrap_block_summary["lower_95"],
                "bootstrap_delta_block_equalized_prevalence_upper_95": bootstrap_block_summary["upper_95"],
                "bootstrap_delta_block_equalized_sign_flip_fraction": bootstrap_block_summary["sign_flip_fraction"],
                "bootstrap_delta_block_equalized_non_estimable_fraction": bootstrap_block_summary["non_estimable_fraction"],
                "dominant_block_id": "" if dominant_row is None else dominant_row["dropped_block_id"],
                "dominant_block_level": "" if dominant_row is None else dominant_row["dropped_block_level"],
                "dominant_block_removed_delta_naive_prevalence": dominant_naive,
                "dominant_block_removed_delta_block_equalized_prevalence": dominant_block,
                "study_dependence_call": study_call,
                "notes": "Subblock bootstrap uses metadata-defined outbreak proxy blocks sampled with replacement within each epoch.",
            }
        )

    leave_one_out = pd.DataFrame(leave_one_out_rows).sort_values(
        ["country_iso3", "comparison_id", "row_type", "dropped_block_size_combined", "dropped_block_id"],
        ascending=[True, True, True, False, True],
    )
    bootstrap = pd.DataFrame(bootstrap_rows).sort_values(["country_iso3", "comparison_id"]).reset_index(drop=True)
    dominance = pd.DataFrame(dominance_rows).sort_values(
        ["country_iso3", "comparison_id", "epoch_slot", "block_total_samples", "block_id"],
        ascending=[True, True, True, False, True],
    )
    return leave_one_out, bootstrap, dominance


def concentration_metrics_from_weights(weights: pd.Series) -> dict[str, object]:
    weights = weights.loc[weights.gt(0)].sort_values(ascending=False)
    if weights.empty:
        return {
            "dominant_prn_event_id": "",
            "dominant_event_share": np.nan,
            "top3_share": np.nan,
            "hhi": np.nan,
            "effective_number": np.nan,
            "top3_prn_event_ids": "",
        }
    dominant_event_id = str(weights.index[0])
    dominant_share = float(weights.iloc[0])
    top3_share = float(weights.iloc[:3].sum())
    hhi = float((weights**2).sum())
    effective_number = safe_divide(1.0, hhi)
    top3_ids = ";".join(str(index) for index in weights.index[:3])
    return {
        "dominant_prn_event_id": dominant_event_id,
        "dominant_event_share": dominant_share,
        "top3_share": top3_share,
        "hhi": hhi,
        "effective_number": effective_number,
        "top3_prn_event_ids": top3_ids,
    }


def build_structure_reweighting(study_frame: pd.DataFrame) -> pd.DataFrame:
    disrupted = study_frame.loc[study_frame["prn_interpretable_bool"] & study_frame["prn_disrupted_bool"]].copy()
    disrupted = disrupted.loc[disrupted["prn_event_id"].ne("")]
    existing = read_tsv(STRUCTURAL_CONCENTRATION_PATH)
    existing = existing.loc[
        existing["scope"].eq("overall") & existing["mechanism_group"].eq("all")
    ].copy()
    existing["row_type"] = "current_naive_reference"
    existing["weighting_scheme"] = "none"

    rows: list[dict[str, object]] = existing.to_dict(orient="records")
    if disrupted.empty:
        return pd.DataFrame(rows)

    naive_weights = disrupted["prn_event_id"].value_counts(normalize=True).sort_values(ascending=False)
    block_event_weights = (
        disrupted.groupby(["base_block_id", "prn_event_id"], dropna=False)
        .size()
        .reset_index(name="event_count")
    )
    block_event_weights["within_block_event_share"] = (
        block_event_weights["event_count"]
        / block_event_weights.groupby("base_block_id", dropna=False)["event_count"].transform("sum")
    )
    weighted = (
        block_event_weights.groupby("prn_event_id", dropna=False)["within_block_event_share"].mean().sort_values(ascending=False)
    )
    weighted = weighted / weighted.sum()
    largest_block_share = safe_divide(
        int(disrupted["base_block_id"].value_counts().iloc[0]),
        int(len(disrupted)),
    )

    weighted_metrics = concentration_metrics_from_weights(weighted)
    naive_metrics = concentration_metrics_from_weights(naive_weights)
    top3_stable = naive_metrics["top3_prn_event_ids"] == weighted_metrics["top3_prn_event_ids"]
    null_reference = to_numeric(existing.get("null_dominant_event_share_mean", pd.Series(dtype=float))).max()
    interpretation = (
        "highly_concentrated_retained_under_study_weighting"
        if top3_stable and pd.notna(null_reference) and weighted_metrics["dominant_event_share"] > null_reference
        else "concentrated_but_partly_study_amplified"
    )
    rows.append(
        {
            "row_type": "study_block_equalized",
            "scope": "overall",
            "mechanism_group": "all",
            "n_genomes": int(len(disrupted)),
            "n_unique_events": int(disrupted["prn_event_id"].nunique()),
            "dominant_prn_event_id": weighted_metrics["dominant_prn_event_id"],
            "dominant_event_count": "",
            "dominant_event_share": weighted_metrics["dominant_event_share"],
            "top3_share": weighted_metrics["top3_share"],
            "top5_share": float(weighted.iloc[:5].sum()),
            "hhi": weighted_metrics["hhi"],
            "shannon_entropy": float(-(weighted * np.log(weighted)).sum()),
            "effective_number": weighted_metrics["effective_number"],
            "gini": np.nan,
            "null_draws": "",
            "null_model": "",
            "null_dominant_event_share_mean": "",
            "null_dominant_event_share_p_ge_observed": "",
            "null_top3_share_mean": "",
            "null_top3_share_p_ge_observed": "",
            "null_effective_number_mean": "",
            "null_effective_number_p_le_observed": "",
            "null_gini_mean": "",
            "null_gini_p_ge_observed": "",
            "stratum_definition": "study_block_equalized",
            "n_strata": int(disrupted["base_block_id"].nunique()),
            "dominant_event_stratum_count": "",
            "dominant_event_stratum_share": "",
            "top3_event_stratum_share": "",
            "null_dominant_event_stratum_share_mean": "",
            "null_dominant_event_stratum_share_p_ge_observed": "",
            "null_top3_event_stratum_share_mean": "",
            "null_top3_event_stratum_share_p_ge_observed": "",
            "notes": "Equal-weighted across study/BioProject blocks among interpretable disrupted genomes only.",
            "weighting_scheme": "study_block_equalized",
            "largest_block_share": largest_block_share,
            "top3_prn_event_ids": weighted_metrics["top3_prn_event_ids"],
            "top3_matches_current_naive": top3_stable,
            "interpretation_call": interpretation,
        }
    )

    if largest_block_share > 0.10:
        largest_block_id = disrupted["base_block_id"].value_counts().index[0]
        dropped = disrupted.loc[disrupted["base_block_id"].ne(largest_block_id)].copy()
        if not dropped.empty:
            dropped_naive = dropped["prn_event_id"].value_counts(normalize=True).sort_values(ascending=False)
            dropped_metrics = concentration_metrics_from_weights(dropped_naive)
            rows.append(
                {
                    "row_type": "drop_largest_block_naive",
                    "scope": "overall",
                    "mechanism_group": "all",
                    "n_genomes": int(len(dropped)),
                    "n_unique_events": int(dropped["prn_event_id"].nunique()),
                    "dominant_prn_event_id": dropped_metrics["dominant_prn_event_id"],
                    "dominant_event_count": int(dropped["prn_event_id"].value_counts().iloc[0]),
                    "dominant_event_share": dropped_metrics["dominant_event_share"],
                    "top3_share": dropped_metrics["top3_share"],
                    "top5_share": float(dropped_naive.iloc[:5].sum()),
                    "hhi": dropped_metrics["hhi"],
                    "shannon_entropy": float(-(dropped_naive * np.log(dropped_naive)).sum()),
                    "effective_number": dropped_metrics["effective_number"],
                    "gini": np.nan,
                    "null_draws": "",
                    "null_model": "",
                    "null_dominant_event_share_mean": "",
                    "null_dominant_event_share_p_ge_observed": "",
                    "null_top3_share_mean": "",
                    "null_top3_share_p_ge_observed": "",
                    "null_effective_number_mean": "",
                    "null_effective_number_p_le_observed": "",
                    "null_gini_mean": "",
                    "null_gini_p_ge_observed": "",
                    "stratum_definition": "drop_largest_block",
                    "n_strata": int(dropped["base_block_id"].nunique()),
                    "dominant_event_stratum_count": "",
                    "dominant_event_stratum_share": "",
                    "top3_event_stratum_share": "",
                    "null_dominant_event_stratum_share_mean": "",
                    "null_dominant_event_stratum_share_p_ge_observed": "",
                    "null_top3_event_stratum_share_mean": "",
                    "null_top3_event_stratum_share_p_ge_observed": "",
                    "notes": f"Observed event weights after removing largest block {largest_block_id}.",
                    "weighting_scheme": "drop_largest_block_naive",
                    "largest_block_share": largest_block_share,
                    "top3_prn_event_ids": dropped_metrics["top3_prn_event_ids"],
                    "top3_matches_current_naive": dropped_metrics["top3_prn_event_ids"] == naive_metrics["top3_prn_event_ids"],
                    "interpretation_call": interpretation,
                }
            )
    structure = pd.DataFrame(rows)
    preferred_columns = list(existing.columns) + [
        "weighting_scheme",
        "largest_block_share",
        "top3_prn_event_ids",
        "top3_matches_current_naive",
        "interpretation_call",
    ]
    ordered = [column for column in preferred_columns if column in structure.columns]
    ordered = list(dict.fromkeys(ordered))
    return structure.loc[:, ordered]


def build_asr_study_block_summary() -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    if ASR_RESAMPLING_REPLICATE_PATH.exists():
        replicate = read_tsv(ASR_RESAMPLING_REPLICATE_PATH)
        replicate = replicate.loc[replicate["scheme"].eq("study_block_balanced")].copy()
        if not replicate.empty:
            replicate.insert(0, "row_type", "replicate")
            frames.append(replicate)
    if ASR_RESAMPLING_SUMMARY_PATH.exists():
        summary = read_tsv(ASR_RESAMPLING_SUMMARY_PATH)
        summary = summary.loc[summary["scheme"].eq("study_block_balanced")].copy()
        if not summary.empty:
            summary.insert(0, "row_type", "summary")
            frames.append(summary)
    if not frames:
        return pd.DataFrame(
            [
                {
                    "row_type": "summary",
                    "scheme": "study_block_balanced",
                    "notes": "Study-block-balanced ASR resampling summary unavailable because no resampling outputs were found.",
                }
            ]
        )
    return pd.concat(frames, ignore_index=True, sort=False)


def build_combined_supplementary_tables(
    assignment: pd.DataFrame,
    dominance: pd.DataFrame,
    leave_one_out: pd.DataFrame,
    bootstrap: pd.DataFrame,
    structure: pd.DataFrame,
    asr_summary: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    assignment_supp = assignment.copy()
    assignment_supp.insert(0, "row_type", "sample_block_assignment")
    dominance_supp = dominance.copy()
    dominance_supp.insert(0, "row_type", "epoch_block_dominance")
    supp52 = pd.concat([assignment_supp, dominance_supp], ignore_index=True, sort=False)

    supp53 = pd.concat([leave_one_out, bootstrap], ignore_index=True, sort=False)
    structure_supp = structure.loc[:, ~structure.columns.duplicated()].copy()
    asr_supp = asr_summary.loc[:, ~asr_summary.columns.duplicated()].copy()
    supp54 = pd.concat([structure_supp, asr_supp], ignore_index=True, sort=False)
    return supp52, supp53, supp54


def main() -> None:
    manifest = prepare_manifest()
    inventory = prepare_inventory()
    history = build_epoch_lookup()
    assignment, study_frame = build_block_assignment(manifest, inventory, history)
    leave_one_out, bootstrap, dominance = build_selected_country_block_audit(study_frame)
    structure = build_structure_reweighting(study_frame)
    asr_summary = build_asr_study_block_summary()
    supp52, supp53, supp54 = build_combined_supplementary_tables(
        assignment, dominance, leave_one_out, bootstrap, structure, asr_summary
    )

    write_tsv(assignment, BLOCK_ASSIGNMENT_PATH)
    write_tsv(leave_one_out, LOO_PATH)
    write_tsv(bootstrap, BOOTSTRAP_PATH)
    write_tsv(dominance, DOMINANCE_PATH)
    write_tsv(structure, STRUCTURE_REWEIGHTED_PATH)
    write_tsv(asr_summary, ASR_STUDY_BLOCK_RESAMPLING_PATH)
    write_tsv(supp52, SUPP52_PATH)
    write_tsv(supp53, SUPP53_PATH)
    write_tsv(supp54, SUPP54_PATH)


if __name__ == "__main__":
    main()
