#!/usr/bin/env python3
"""Estimate country-year prn disruption prevalence with IPW correction."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_bool(value: object) -> bool:
    text = str(value).strip().lower()
    return text in {"true", "1", "yes"}


def coerce_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def load_prediction_map(missingness_model_path: str) -> tuple[dict[str, float], float]:
    with open(missingness_model_path, encoding="utf-8") as handle:
        payload = json.load(handle)
    prediction_map = {
        row["sample_id_canonical"]: float(row["prob_interpretable"])
        for row in payload.get("predictions", [])
        if row.get("sample_id_canonical")
    }
    fallback = float(payload.get("fallback_probability", payload.get("prediction_summary", {}).get("mean_probability", 0.5)))
    return prediction_map, fallback


def build_prevalence(
    manifest_path: str,
    missingness_model_path: str,
    ph_master_path: str,
    output_prevalence_path: str,
    output_figure_path: str,
    weight_truncation: float,
) -> pd.DataFrame:
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    ph_master = pd.read_csv(ph_master_path, sep="\t", dtype=str)
    prediction_map, fallback_probability = load_prediction_map(missingness_model_path)

    manifest["country_iso3"] = manifest.get("country_iso3", pd.Series(index=manifest.index, dtype=str)).fillna("")
    manifest["year"] = coerce_numeric(manifest.get("year", pd.Series(index=manifest.index, dtype=str)))
    manifest["prn_interpretable_flag"] = manifest.get("prn_interpretable", pd.Series(index=manifest.index, dtype=str)).map(parse_bool)
    manifest["prn_disrupted_flag"] = manifest.get("prn_disrupted", pd.Series(index=manifest.index, dtype=str)).map(parse_bool)
    manifest = manifest.loc[manifest["country_iso3"].ne("") & manifest["year"].notna()].copy()

    manifest["prob_interpretable"] = manifest["sample_id_canonical"].map(prediction_map).fillna(fallback_probability)
    manifest["prob_interpretable"] = manifest["prob_interpretable"].clip(lower=1e-6, upper=1.0)
    manifest["ipw_weight"] = (1.0 / manifest["prob_interpretable"]).clip(upper=float(weight_truncation))
    manifest["weighted_interpretable"] = np.where(manifest["prn_interpretable_flag"], manifest["ipw_weight"], 0.0)
    manifest["weighted_disrupted"] = np.where(
        manifest["prn_interpretable_flag"] & manifest["prn_disrupted_flag"],
        manifest["ipw_weight"],
        0.0,
    )

    aggregated = (
        manifest.groupby(["country_iso3", "year"], dropna=False)
        .agg(
            n_genomes_total=("sample_id_canonical", "count"),
            n_genomes_prn_interpretable=("prn_interpretable_flag", "sum"),
            n_prn_disrupted=("prn_disrupted_flag", "sum"),
            ipw_weight_total=("weighted_interpretable", "sum"),
            ipw_weighted_disrupted=("weighted_disrupted", "sum"),
            mean_probability=("prob_interpretable", "mean"),
            mean_ipw_weight=("ipw_weight", "mean"),
            max_ipw_weight=("ipw_weight", "max"),
        )
        .reset_index()
    )
    aggregated["year"] = aggregated["year"].astype(int)
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
    aggregated["boundary_lower_prevalence"] = aggregated["n_prn_disrupted"] / aggregated["n_genomes_total"]
    aggregated["boundary_upper_prevalence"] = (
        aggregated["n_prn_disrupted"]
        + (aggregated["n_genomes_total"] - aggregated["n_genomes_prn_interpretable"])
    ) / aggregated["n_genomes_total"]
    aggregated["n_missing_outcomes"] = aggregated["n_genomes_total"] - aggregated["n_genomes_prn_interpretable"]
    aggregated["weight_truncation"] = float(weight_truncation)
    aggregated["prediction_fallback_probability"] = fallback_probability

    ph_columns = [
        "country_iso3",
        "year",
        "country_name",
        "reported_cases",
        "dtp3_coverage",
        "genomes_count",
        "vaccine_program_type",
        "post_covid_period",
    ]
    ph_subset = ph_master[[column for column in ph_columns if column in ph_master.columns]].copy()
    ph_subset["year"] = coerce_numeric(ph_subset.get("year", pd.Series(dtype=str)))
    ph_subset = ph_subset.dropna(subset=["year"])
    ph_subset["year"] = ph_subset["year"].astype(int)
    aggregated = aggregated.merge(ph_subset, on=["country_iso3", "year"], how="left")

    output_path = Path(output_prevalence_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    aggregated.sort_values(["country_iso3", "year"]).to_csv(output_path, sep="\t", index=False)

    figure_path = Path(output_figure_path)
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    country_plot = (
        aggregated.groupby("country_iso3", dropna=False)
        .agg(
            naive_prevalence=("naive_prevalence", "mean"),
            ipw_prevalence=("ipw_prevalence", "mean"),
            boundary_lower_prevalence=("boundary_lower_prevalence", "mean"),
            boundary_upper_prevalence=("boundary_upper_prevalence", "mean"),
            n_country_years=("year", "count"),
        )
        .reset_index()
        .sort_values("ipw_prevalence", ascending=False)
        .head(20)
    )

    plt.rcParams.update({"figure.figsize": (10.5, 5.0)})
    figure, axis = plt.subplots()
    positions = np.arange(len(country_plot))
    axis.plot(positions, country_plot["naive_prevalence"], marker="o", label="Naive")
    axis.plot(positions, country_plot["ipw_prevalence"], marker="o", label="IPW")
    axis.vlines(
        positions,
        country_plot["boundary_lower_prevalence"],
        country_plot["boundary_upper_prevalence"],
        color="0.7",
        linewidth=2,
        label="Missing-outcome bounds",
    )
    axis.set_xticks(positions, country_plot["country_iso3"], rotation=60, ha="right")
    axis.set_ylabel("Disrupted prn prevalence")
    axis.set_title("Country-level prevalence: naive, IPW, and boundary bounds")
    axis.legend()
    figure.tight_layout()
    figure.savefig(figure_path)
    plt.close(figure)

    return aggregated


if "snakemake" in globals():
    build_prevalence(
        manifest_path=snakemake.input.manifest,
        missingness_model_path=snakemake.input.missingness_model,
        ph_master_path=snakemake.input.ph_master,
        output_prevalence_path=snakemake.output.prevalence,
        output_figure_path=snakemake.output.boundary_figure,
        weight_truncation=float(snakemake.params.weight_truncation),
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Build IPW-corrected prevalence table")
    parser.add_argument("--manifest", required=True, help="Manifest TSV")
    parser.add_argument("--missingness-model", required=True, help="Missingness model JSON")
    parser.add_argument("--ph-master", required=True, help="Public health master TSV")
    parser.add_argument("--prevalence-out", required=True, help="Output prevalence TSV")
    parser.add_argument("--boundary-figure-out", required=True, help="Output boundary sensitivity PDF")
    parser.add_argument("--weight-truncation", type=float, default=20.0, help="Maximum IP weight")
    arguments = parser.parse_args()

    build_prevalence(
        manifest_path=arguments.manifest,
        missingness_model_path=arguments.missingness_model,
        ph_master_path=arguments.ph_master,
        output_prevalence_path=arguments.prevalence_out,
        output_figure_path=arguments.boundary_figure_out,
        weight_truncation=arguments.weight_truncation,
    )