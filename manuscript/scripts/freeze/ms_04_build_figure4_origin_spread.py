#!/usr/bin/env python3
"""Build Figure 4 manuscript-facing tables for origin-clade expansion and event-centered spread."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
FIGURE_DATA_DIR = REPO_ROOT / "manuscript/figure_data"

ORIGIN_OUTPUT = FIGURE_DATA_DIR / "figure4_origin_clade_expansion.tsv"
EVENT_OUTPUT = FIGURE_DATA_DIR / "figure4_event_centered_country.tsv"
EVENT_POOLED_OUTPUT = FIGURE_DATA_DIR / "figure4_event_centered_pooled.tsv"

ORIGIN_COLUMNS = [
    "panel_id",
    "origin_id",
    "origin_rank",
    "mechanism_group",
    "dominant_prn_mechanism",
    "country_iso3",
    "country_name",
    "year",
    "relative_year",
    "yearly_disrupted_descendants",
    "cumulative_disrupted_descendants",
    "n_tips_disrupted",
    "n_tips_total",
    "n_countries",
    "first_year",
    "last_year",
    "duration_years",
    "branch_support",
    "origin_support_score",
    "notes",
]

EVENT_COLUMNS = [
    "event_type",
    "country_iso3",
    "country_name",
    "event_year",
    "year",
    "relative_year",
    "ipw_prevalence",
    "naive_prevalence",
    "n_genomes_prn_interpretable",
    "n_prn_disrupted",
    "n_origin_clade_descendants",
    "n_origin_clades_active",
    "n_new_origins_detected",
    "observed_relative_years",
    "notes",
]

EVENT_POOLED_COLUMNS = [
    "panel_id",
    "event_type",
    "metric_name",
    "relative_year",
    "n_countries",
    "mean_value",
    "median_value",
    "mean_difference",
    "median_difference",
    "n_pairs",
    "notes",
]


def mechanism_group(value: str) -> str:
    text = str(value or "")
    if "is481" in text:
        return "IS481"
    if "inversion" in text or "rearrangement" in text:
        return "Inversion/rearrangement"
    return "Other"


def mode_or_empty(values: pd.Series) -> str:
    non_missing = values.dropna().astype(str)
    if non_missing.empty:
        return ""
    return non_missing.mode().iloc[0]


def load_origin_descendants() -> pd.DataFrame:
    rows = []
    subtree_dir = REPO_ROOT / "outputs/workflow/asr/event_subtrees"
    for path in sorted(subtree_dir.glob("origin_*.descendant_tips.tsv")):
        df = pd.read_csv(path, sep="\t")
        df["origin_id"] = path.stem.split(".")[0]
        rows.append(df)

    descendants = pd.concat(rows, ignore_index=True)
    descendants["year"] = pd.to_numeric(descendants["year"], errors="coerce")
    return descendants


def build_origin_table() -> pd.DataFrame:
    origin_events = pd.read_csv(REPO_ROOT / "outputs/workflow/asr/origin_events.tsv", sep="\t")
    descendants = load_origin_descendants()
    disrupted = descendants[descendants["observed_prn_state"] == "disrupted"].copy()

    country_name_map = (
        pd.read_csv(REPO_ROOT / "outputs/workflow/epi/ipw_prevalence.tsv", sep="\t")[["country_iso3", "country_name"]]
        .dropna()
        .drop_duplicates()
        .set_index("country_iso3")["country_name"]
        .to_dict()
    )

    disrupted_summary = (
        disrupted.groupby("origin_id", dropna=False)
        .agg(
            disrupted_descendant_count=("tip_label", "count"),
            disrupted_country_count=("country_iso3", lambda values: values.dropna().nunique()),
            first_observed_year=("year", "min"),
            last_observed_year=("year", "max"),
            dominant_country_iso3=("country_iso3", mode_or_empty),
        )
        .reset_index()
    )

    summary = origin_events.merge(disrupted_summary, on="origin_id", how="left")
    summary["mechanism_group"] = summary["dominant_prn_mechanism"].map(mechanism_group)
    summary["first_year"] = pd.to_numeric(summary["first_year"], errors="coerce").fillna(summary["first_observed_year"])
    summary["last_year"] = pd.to_numeric(summary["last_year"], errors="coerce").fillna(summary["last_observed_year"])
    summary["duration_years"] = summary.apply(
        lambda row: int(row["last_year"] - row["first_year"] + 1)
        if pd.notna(row["first_year"]) and pd.notna(row["last_year"])
        else np.nan,
        axis=1,
    )
    summary["country_iso3"] = summary["dominant_country_iso3"].fillna("")
    summary["country_name"] = summary["country_iso3"].map(country_name_map).fillna("")
    summary["n_tips_disrupted"] = summary["n_tips_disrupted"].fillna(summary["disrupted_descendant_count"])
    rank_df = (
        summary.sort_values(["n_tips_disrupted", "duration_years", "origin_id"], ascending=[False, False, True])
        .reset_index(drop=True)
        .reset_index()
        .rename(columns={"index": "origin_rank"})[["origin_id", "origin_rank"]]
    )
    summary = summary.merge(rank_df, on="origin_id", how="left")

    summary_rows = []
    year_rows = []
    yearly_counts = (
        disrupted.groupby(["origin_id", "year"], dropna=False)
        .size()
        .rename("yearly_disrupted_descendants")
        .reset_index()
    )

    for row in summary.itertuples(index=False):
        summary_rows.append(
            {
                "panel_id": "origin_summary",
                "origin_id": row.origin_id,
                "origin_rank": int(row.origin_rank) + 1,
                "mechanism_group": row.mechanism_group,
                "dominant_prn_mechanism": row.dominant_prn_mechanism,
                "country_iso3": row.country_iso3,
                "country_name": row.country_name,
                "year": np.nan,
                "relative_year": np.nan,
                "yearly_disrupted_descendants": np.nan,
                "cumulative_disrupted_descendants": np.nan,
                "n_tips_disrupted": row.n_tips_disrupted,
                "n_tips_total": row.n_tips_total,
                "n_countries": row.n_countries,
                "first_year": row.first_year,
                "last_year": row.last_year,
                "duration_years": row.duration_years,
                "branch_support": row.branch_support,
                "origin_support_score": row.origin_support_score,
                "notes": (
                    "single_country_origin"
                    if pd.notna(row.n_countries) and row.n_countries <= 1
                    else "cross_country_origin"
                ),
            }
        )

        if pd.isna(row.first_year) or pd.isna(row.last_year):
            continue

        origin_counts = yearly_counts[yearly_counts["origin_id"] == row.origin_id].set_index("year")["yearly_disrupted_descendants"].to_dict()
        cumulative = 0
        for year in range(int(row.first_year), int(row.last_year) + 1):
            yearly_value = int(origin_counts.get(float(year), origin_counts.get(year, 0)))
            cumulative += yearly_value
            year_rows.append(
                {
                    "panel_id": "origin_year",
                    "origin_id": row.origin_id,
                    "origin_rank": int(row.origin_rank) + 1,
                    "mechanism_group": row.mechanism_group,
                    "dominant_prn_mechanism": row.dominant_prn_mechanism,
                    "country_iso3": row.country_iso3,
                    "country_name": row.country_name,
                    "year": year,
                    "relative_year": year - int(row.first_year),
                    "yearly_disrupted_descendants": yearly_value,
                    "cumulative_disrupted_descendants": cumulative,
                    "n_tips_disrupted": row.n_tips_disrupted,
                    "n_tips_total": row.n_tips_total,
                    "n_countries": row.n_countries,
                    "first_year": row.first_year,
                    "last_year": row.last_year,
                    "duration_years": row.duration_years,
                    "branch_support": row.branch_support,
                    "origin_support_score": row.origin_support_score,
                    "notes": "filled_origin_duration_window",
                }
            )

    output = pd.DataFrame(summary_rows + year_rows)
    return output[ORIGIN_COLUMNS].sort_values(["panel_id", "origin_rank", "year"])


def build_event_centered_table(origin_table: pd.DataFrame) -> pd.DataFrame:
    ipw = pd.read_csv(REPO_ROOT / "outputs/workflow/epi/ipw_prevalence.tsv", sep="\t")
    descendants = load_origin_descendants()
    disrupted = descendants[descendants["observed_prn_state"] == "disrupted"].copy()
    country_name_lookup = (
        ipw.groupby("country_iso3", dropna=False)["country_name"]
        .agg(mode_or_empty)
        .to_dict()
    )

    clade_country_year = (
        disrupted.groupby(["country_iso3", "year"], dropna=False)
        .agg(
            n_origin_clade_descendants=("tip_label", "count"),
            n_origin_clades_active=("origin_id", "nunique"),
        )
        .reset_index()
    )

    origin_first_years = (
        origin_table[origin_table["panel_id"] == "origin_summary"][["origin_id", "country_iso3", "first_year"]]
        .dropna(subset=["country_iso3", "first_year"])
        .copy()
    )
    new_origin_country_year = (
        origin_first_years.groupby(["country_iso3", "first_year"], dropna=False)
        .size()
        .rename("n_new_origins_detected")
        .reset_index()
        .rename(columns={"first_year": "year"})
    )

    detection_events = (
        ipw.loc[ipw["n_prn_disrupted"].fillna(0) > 0, ["country_iso3", "year"]]
        .groupby("country_iso3", dropna=False)["year"]
        .min()
        .reset_index()
        .rename(columns={"year": "event_year"})
        .assign(event_type="first_prn_detection")
    )
    detection_events["country_name"] = detection_events["country_iso3"].map(country_name_lookup).fillna("")

    local_origin_events = (
        origin_first_years.groupby("country_iso3", dropna=False)["first_year"]
        .min()
        .reset_index()
        .rename(columns={"first_year": "event_year"})
        .assign(event_type="first_local_origin")
    )
    local_origin_events["country_name"] = local_origin_events["country_iso3"].map(country_name_lookup).fillna("")

    prevalence = ipw[["country_iso3", "country_name", "year", "ipw_prevalence", "naive_prevalence", "n_genomes_prn_interpretable", "n_prn_disrupted"]].copy()

    rows = []
    for event_row in pd.concat([detection_events, local_origin_events], ignore_index=True).itertuples(index=False):
        for relative_year in range(-3, 4):
            year = int(event_row.event_year) + relative_year
            prevalence_match = prevalence[
                (prevalence["country_iso3"] == event_row.country_iso3) & (prevalence["year"] == year)
            ]
            clade_match = clade_country_year[
                (clade_country_year["country_iso3"] == event_row.country_iso3) & (clade_country_year["year"] == year)
            ]
            new_origin_match = new_origin_country_year[
                (new_origin_country_year["country_iso3"] == event_row.country_iso3) & (new_origin_country_year["year"] == year)
            ]

            prevalence_row = prevalence_match.iloc[0] if not prevalence_match.empty else None
            clade_row = clade_match.iloc[0] if not clade_match.empty else None
            new_origin_row = new_origin_match.iloc[0] if not new_origin_match.empty else None

            rows.append(
                {
                    "event_type": event_row.event_type,
                    "country_iso3": event_row.country_iso3,
                    "country_name": event_row.country_name if pd.notna(event_row.country_name) else "",
                    "event_year": int(event_row.event_year),
                    "year": year,
                    "relative_year": relative_year,
                    "ipw_prevalence": prevalence_row["ipw_prevalence"] if prevalence_row is not None else np.nan,
                    "naive_prevalence": prevalence_row["naive_prevalence"] if prevalence_row is not None else np.nan,
                    "n_genomes_prn_interpretable": prevalence_row["n_genomes_prn_interpretable"] if prevalence_row is not None else np.nan,
                    "n_prn_disrupted": prevalence_row["n_prn_disrupted"] if prevalence_row is not None else np.nan,
                    "n_origin_clade_descendants": clade_row["n_origin_clade_descendants"] if clade_row is not None else np.nan,
                    "n_origin_clades_active": clade_row["n_origin_clades_active"] if clade_row is not None else np.nan,
                    "n_new_origins_detected": new_origin_row["n_new_origins_detected"] if new_origin_row is not None else 0,
                    "observed_relative_years": np.nan,
                    "notes": "event_centered_window",
                }
            )

    event_df = pd.DataFrame(rows)
    observed_counts = (
        event_df.assign(has_signal=event_df["ipw_prevalence"].notna() | event_df["n_origin_clade_descendants"].notna())
        .groupby(["event_type", "country_iso3"], dropna=False)["has_signal"]
        .sum()
        .rename("observed_relative_years")
        .reset_index()
    )
    event_df = event_df.drop(columns=["observed_relative_years"]).merge(
        observed_counts, on=["event_type", "country_iso3"], how="left"
    )
    return event_df[EVENT_COLUMNS].sort_values(["event_type", "country_iso3", "relative_year"])


def build_event_pooled_table(event_table: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "ipw_prevalence",
        "naive_prevalence",
        "n_origin_clade_descendants",
        "n_origin_clades_active",
        "n_new_origins_detected",
    ]
    rows: list[dict[str, object]] = []

    for event_type, event_group in event_table.groupby("event_type", dropna=False):
        for metric_name in metrics:
            metric_group = event_group[["country_iso3", "relative_year", metric_name]].copy()
            metric_group[metric_name] = pd.to_numeric(metric_group[metric_name], errors="coerce")
            relative_summary = (
                metric_group.dropna(subset=[metric_name])
                .groupby("relative_year", dropna=False)[metric_name]
                .agg(["mean", "median", "count"])
                .reset_index()
            )
            for row in relative_summary.itertuples(index=False):
                rows.append(
                    {
                        "panel_id": "pooled_relative_year",
                        "event_type": event_type,
                        "metric_name": metric_name,
                        "relative_year": int(row.relative_year),
                        "n_countries": int(row.count),
                        "mean_value": row.mean,
                        "median_value": row.median,
                        "mean_difference": np.nan,
                        "median_difference": np.nan,
                        "n_pairs": np.nan,
                        "notes": "event_centered_country_mean_and_median",
                    }
                )

            pre = (
                metric_group[metric_group["relative_year"].between(-3, -1, inclusive="both")]
                .groupby("country_iso3", dropna=False)[metric_name]
                .mean()
            )
            post = (
                metric_group[metric_group["relative_year"].between(1, 3, inclusive="both")]
                .groupby("country_iso3", dropna=False)[metric_name]
                .mean()
            )
            paired = (
                pd.concat([pre.rename("pre"), post.rename("post")], axis=1)
                .dropna(subset=["pre", "post"])
                .assign(difference=lambda frame: frame["post"] - frame["pre"])
            )
            rows.append(
                {
                    "panel_id": "pre_post_difference",
                    "event_type": event_type,
                    "metric_name": metric_name,
                    "relative_year": np.nan,
                    "n_countries": int(len(paired)),
                    "mean_value": np.nan,
                    "median_value": np.nan,
                    "mean_difference": np.nan if paired.empty else paired["difference"].mean(),
                    "median_difference": np.nan if paired.empty else paired["difference"].median(),
                    "n_pairs": int(len(paired)),
                    "notes": "post_years_1_to_3_minus_pre_years_minus3_to_minus1",
                }
            )

    pooled = pd.DataFrame(rows)
    return pooled[EVENT_POOLED_COLUMNS].sort_values(["panel_id", "event_type", "metric_name", "relative_year"])


def main() -> int:
    FIGURE_DATA_DIR.mkdir(parents=True, exist_ok=True)

    origin_table = build_origin_table()
    origin_table.to_csv(ORIGIN_OUTPUT, sep="\t", index=False)

    event_table = build_event_centered_table(origin_table)
    event_table.to_csv(EVENT_OUTPUT, sep="\t", index=False)

    pooled_table = build_event_pooled_table(event_table)
    pooled_table.to_csv(EVENT_POOLED_OUTPUT, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
