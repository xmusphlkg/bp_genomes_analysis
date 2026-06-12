#!/usr/bin/env python3
"""Assess interrupted time-series feasibility for routine aP introduction events."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def coerce_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def first_routine_ap_year(ph_master: pd.DataFrame) -> pd.DataFrame:
    routine_rows = ph_master.loc[ph_master["vaccine_program_type"].eq("ap_introduced_routine_or_mixed")].copy()
    return (
        routine_rows.groupby("country_iso3", dropna=False)
        .agg(first_routine_ap_year=("year", "min"), country_name=("country_name", "first"))
        .reset_index()
    )


def classify_feasibility(pre_rows: int, post_rows: int, observed_min: int | None, observed_max: int | None, event_year: int | None) -> tuple[str, str]:
    if event_year is None or observed_min is None or observed_max is None:
        return "not_feasible", "no_routine_ap_intro_detected"
    if event_year <= observed_min:
        return "not_feasible", "observed_panel_starts_at_or_after_intro"
    if event_year >= observed_max:
        return "not_feasible", "observed_panel_ends_before_post_intro_period"
    if pre_rows >= 3 and post_rows >= 3:
        return "feasible", "at_least_three_pre_and_post_country_years"
    if pre_rows >= 2 and post_rows >= 2:
        return "borderline", "only_two_pre_or_post_country_years"
    return "not_feasible", "insufficient_pre_post_country_years"


def build_feasibility_report(prevalence_path: str, ph_master_path: str, output_report_path: str) -> pd.DataFrame:
    prevalence = pd.read_csv(prevalence_path, sep="\t", dtype=str)
    ph_master = pd.read_csv(ph_master_path, sep="\t", dtype=str)

    prevalence["year"] = coerce_numeric(prevalence.get("year", pd.Series(dtype=str)))
    prevalence["n_genomes_prn_interpretable"] = coerce_numeric(
        prevalence.get("n_genomes_prn_interpretable", pd.Series(dtype=str))
    ).fillna(0)
    prevalence = prevalence.loc[prevalence["year"].notna() & prevalence["country_iso3"].notna()].copy()
    prevalence["year"] = prevalence["year"].astype(int)
    prevalence = prevalence.loc[prevalence["n_genomes_prn_interpretable"] >= 5].copy()

    ph_master["year"] = coerce_numeric(ph_master.get("year", pd.Series(dtype=str)))
    ph_master = ph_master.loc[ph_master["year"].notna() & ph_master["country_iso3"].notna()].copy()
    ph_master["year"] = ph_master["year"].astype(int)

    event_years = first_routine_ap_year(ph_master)

    records: list[dict[str, object]] = []
    for country_iso3, group in prevalence.groupby("country_iso3", dropna=False):
        event_row = event_years.loc[event_years["country_iso3"].eq(country_iso3)]
        event_year = int(event_row["first_routine_ap_year"].iloc[0]) if not event_row.empty else None
        country_name = str(event_row["country_name"].iloc[0]) if not event_row.empty else str(group.get("country_name", pd.Series()).dropna().iloc[0] if group.get("country_name") is not None and not group.get("country_name").dropna().empty else country_iso3)
        observed_min = int(group["year"].min()) if not group.empty else None
        observed_max = int(group["year"].max()) if not group.empty else None
        pre_rows = int((group["year"] < event_year).sum()) if event_year is not None else 0
        post_rows = int((group["year"] >= event_year).sum()) if event_year is not None else 0
        feasibility, reason = classify_feasibility(pre_rows, post_rows, observed_min, observed_max, event_year)
        records.append(
            {
                "country_iso3": country_iso3,
                "country_name": country_name,
                "policy_event": "first_routine_ap_intro",
                "first_routine_ap_year": event_year,
                "observed_year_min": observed_min,
                "observed_year_max": observed_max,
                "n_country_years": int(len(group)),
                "n_pre_intro_country_years": pre_rows,
                "n_post_intro_country_years": post_rows,
                "its_feasibility": feasibility,
                "reason": reason,
            }
        )

    report = pd.DataFrame.from_records(records).sort_values(["its_feasibility", "country_iso3"])
    output_path = Path(output_report_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(output_path, sep="\t", index=False)
    return report


if "snakemake" in globals():
    build_feasibility_report(
        prevalence_path=snakemake.input.prevalence,
        ph_master_path=snakemake.input.ph_master,
        output_report_path=snakemake.output.report,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Assess ITS feasibility")
    parser.add_argument("--prevalence", required=True, help="IPW prevalence TSV")
    parser.add_argument("--ph-master", required=True, help="Public health master TSV")
    parser.add_argument("--report-out", required=True, help="Output report TSV")
    arguments = parser.parse_args()

    build_feasibility_report(arguments.prevalence, arguments.ph_master, arguments.report_out)
