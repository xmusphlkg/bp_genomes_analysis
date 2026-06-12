#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from workflow.lib.project_paths import project_module_data_root


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def norm(s: pd.Series) -> pd.Series:
    return s.astype(str).where(~s.isna(), "NA").replace({"nan": "NA", "None": "NA"})


def ensure_cols(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    df = df.copy()
    for c in cols:
        if c not in df.columns:
            df[c] = "NA"
        df[c] = norm(df[c])
    return df


def pick_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def frac_table(
    df: pd.DataFrame,
    group_cols: list[str],
    value_col: str,
    allowed_values: list[str] | None = None,
) -> pd.DataFrame:
    df = df.copy()
    df = ensure_cols(df, group_cols + [value_col])
    if allowed_values is not None:
        df = df[df[value_col].isin(allowed_values)].copy()

    grp = df.groupby(group_cols + [value_col], dropna=False).size().reset_index(name="n")
    totals = df.groupby(group_cols, dropna=False).size().reset_index(name="group_n")
    out = grp.merge(totals, on=group_cols, how="left")
    out["frac"] = out["n"] / out["group_n"]
    out = out.sort_values(group_cols + ["n"], ascending=[True] * len(group_cols) + [False])
    return out


def pivot_frac(
    df: pd.DataFrame,
    index_cols: list[str],
    value_col: str,
    frac_value: str,
    allowed_values: list[str] | None = None,
) -> pd.DataFrame:
    tab = frac_table(df, index_cols, value_col, allowed_values=allowed_values)
    wide = tab.pivot_table(index=index_cols, columns=value_col, values=frac_value, aggfunc="first")
    wide = wide.reset_index()
    wide.columns = [str(c) for c in wide.columns]
    return wide


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Step3D: turn prn_call into publication-ready trend tables (overall + by year×ST), and cross-tabs with AMR/ptxP"
    )
    step2_root = project_module_data_root("step2_typing")
    step3_root = project_module_data_root("step3_prn_scan")
    ap.add_argument(
        "--table",
        required=True,
        help="Merged table with prn columns (bp_step3 outputs/bp_qc_merged_mlst_markers_prn.tsv)",
    )
    ap.add_argument("--outdir", default=str(step3_root / "outputs"), help="Output directory")
    ap.add_argument("--prefix", default="bp", help="Output prefix")
    ap.add_argument(
        "--typing-manifest",
        default=str(step2_root / "outputs" / "bp_genotype_manifest.tsv"),
        help="Optional standardized Step2 genotype manifest TSV.",
    )
    ap.add_argument(
        "--min-group-n",
        type=int,
        default=20,
        help="Minimum group_n for year×ST table rows (filters noisy strata)",
    )
    ap.add_argument(
        "--top-st-n",
        type=int,
        default=20,
        help="Number of top STs to keep in ranking table",
    )
    ap.add_argument(
        "--early-max-year",
        type=int,
        default=2010,
        help="Early window: years <= this value",
    )
    ap.add_argument(
        "--late-min-year",
        type=int,
        default=2016,
        help="Late window: years >= this value",
    )
    args = ap.parse_args()

    table = Path(args.table)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(table, sep="\t", dtype=str)
    typing_path = Path(args.typing_manifest)
    if typing_path.exists() and typing_path.stat().st_size > 0:
        typing = pd.read_csv(typing_path, sep="\t", dtype=str)
        typing_cols = [
            "assembly_accession",
            "ptxP_label",
            "marker_ptxP_promoter_hash",
            "marker_23s_status",
            "23s_A2047G_call_raw",
        ]
        available_typing_cols = [col for col in typing_cols if col in typing.columns]
        if "assembly_accession" in available_typing_cols:
            acc_col = pick_col(df, ["Assembly Accession", "assembly_accession", "Current Accession"])
            if acc_col is not None:
                df = df.merge(
                    typing[available_typing_cols].drop_duplicates(subset=["assembly_accession"]),
                    left_on=acc_col,
                    right_on="assembly_accession",
                    how="left",
                )

    # Normalize column naming across step outputs
    year_col = pick_col(df, ["year", "year_x", "year_y"])
    country_col = pick_col(df, ["country", "country_x", "country_y"])
    st_col = pick_col(df, ["mlst_st", "mlst_st_x", "mlst_st_y"])
    ptxp_col = pick_col(df, ["ptxP_label", "marker_ptxP_promoter_allele", "marker_ptxP_promoter"])
    mr_col = pick_col(df, ["marker_23s_status", "23s_A2047G_call", "23s_A2047G_call_raw"])

    # Canonical working columns
    df["year"] = norm(df[year_col]) if year_col else "NA"
    df["country"] = norm(df[country_col]) if country_col else "NA"
    df["mlst_st"] = norm(df[st_col]) if st_col else "NA"
    df["ptxP"] = norm(df[ptxp_col]) if ptxp_col else "NA"
    df["23s_A2047G_call"] = norm(df[mr_col]) if mr_col else "NA"

    # Core columns used for stratification
    df = ensure_cols(df, ["year", "country", "mlst_st", "prn_call", "23s_A2047G_call", "ptxP"])

    # Keep only interpretable prn_call for most trend analyses
    prn_ok = ["intact", "disrupted_multi_hsp"]
    df_prn = df[df["prn_call"].isin(prn_ok)].copy()

    # Ensure year numeric for windows; non-numeric -> NA
    year_num = pd.to_numeric(df_prn["year"].where(df_prn["year"] != "NA"), errors="coerce")
    df_prn["year_num"] = year_num

    # 1) Overall by-year trend
    by_year = frac_table(df_prn, ["year"], "prn_call", allowed_values=prn_ok)
    by_year.to_csv(outdir / f"{args.prefix}_prn_trend_by_year_clean.tsv", sep="\t", index=False)

    # 2) By year×ST trend (filtered)
    by_year_st = frac_table(df_prn, ["year", "mlst_st"], "prn_call", allowed_values=prn_ok)
    by_year_st_filt = by_year_st[by_year_st["group_n"] >= int(args.min_group_n)].copy()
    by_year_st_filt.to_csv(outdir / f"{args.prefix}_prn_trend_by_year_mlst_st_min{int(args.min_group_n)}.tsv", sep="\t", index=False)

    # 3) ST ranking: overall disrupted fraction, plus early vs late windows
    st_overall = pivot_frac(df_prn, ["mlst_st"], "prn_call", "frac", allowed_values=prn_ok)
    st_counts = df_prn.groupby(["mlst_st"], dropna=False).size().reset_index(name="n_total")
    st_rank = st_overall.merge(st_counts, on=["mlst_st"], how="left")

    # Early/late windows
    early = df_prn[(df_prn["year_num"].notna()) & (df_prn["year_num"] <= int(args.early_max_year))].copy()
    late = df_prn[(df_prn["year_num"].notna()) & (df_prn["year_num"] >= int(args.late_min_year))].copy()

    early_w = pivot_frac(early, ["mlst_st"], "prn_call", "frac", allowed_values=prn_ok)
    early_n = early.groupby(["mlst_st"], dropna=False).size().reset_index(name="n_early")
    early_w = early_w.merge(early_n, on=["mlst_st"], how="left")

    late_w = pivot_frac(late, ["mlst_st"], "prn_call", "frac", allowed_values=prn_ok)
    late_n = late.groupby(["mlst_st"], dropna=False).size().reset_index(name="n_late")
    late_w = late_w.merge(late_n, on=["mlst_st"], how="left")

    st_rank = st_rank.merge(early_w, on=["mlst_st"], how="left", suffixes=("", "_early"))
    st_rank = st_rank.merge(late_w, on=["mlst_st"], how="left", suffixes=("", "_late"))

    # Normalize column names after merges
    # Expected columns: disrupted_multi_hsp, intact, n_total, disrupted_multi_hsp_early/intact_early, n_early, disrupted_multi_hsp_late/intact_late, n_late
    for c in ["disrupted_multi_hsp", "intact", "disrupted_multi_hsp_early", "intact_early", "disrupted_multi_hsp_late", "intact_late"]:
        if c not in st_rank.columns:
            st_rank[c] = pd.NA

    st_rank["delta_disrupted_late_minus_early"] = st_rank["disrupted_multi_hsp_late"].fillna(0) - st_rank["disrupted_multi_hsp_early"].fillna(0)

    st_rank = st_rank.sort_values(["n_total"], ascending=[False])
    st_rank.to_csv(outdir / f"{args.prefix}_prn_top_sts_rank.tsv", sep="\t", index=False)

    # Top ST list (for quick plotting / filtering)
    st_top = st_rank.sort_values(["n_total"], ascending=[False]).head(int(args.top_st_n)).copy()
    st_top[["mlst_st", "n_total", "disrupted_multi_hsp", "intact", "n_early", "n_late", "delta_disrupted_late_minus_early"]].to_csv(
        outdir / f"{args.prefix}_prn_top_sts_top{int(args.top_st_n)}.tsv", sep="\t", index=False
    )

    # 4) Cross-tabs: prn_call vs 23S / ptxP (overall and by year)
    ct_overall = frac_table(df_prn, ["prn_call"], "23s_A2047G_call")
    ct_overall.to_csv(outdir / f"{args.prefix}_prn_vs_23s_overall.tsv", sep="\t", index=False)

    ct_year = frac_table(df_prn, ["year", "prn_call"], "23s_A2047G_call")
    ct_year.to_csv(outdir / f"{args.prefix}_prn_vs_23s_by_year.tsv", sep="\t", index=False)

    ptxp_overall = frac_table(df_prn, ["prn_call"], "ptxP")
    ptxp_overall.to_csv(outdir / f"{args.prefix}_prn_vs_ptxP_overall.tsv", sep="\t", index=False)

    ptxp_year = frac_table(df_prn, ["year", "prn_call"], "ptxP")
    ptxp_year.to_csv(outdir / f"{args.prefix}_prn_vs_ptxP_by_year.tsv", sep="\t", index=False)

    # 5) A short text summary (for results section draft)
    overall_counts = df["prn_call"].value_counts(dropna=False)
    n_total_calls = int(len(df))
    n_prn_ok = int(len(df_prn))

    disrupted_frac = float((df_prn["prn_call"] == "disrupted_multi_hsp").mean()) if n_prn_ok else 0.0

    summary = []
    summary.append(f"Input rows: {n_total_calls}")
    summary.append(f"Rows with prn_call in {{intact, disrupted_multi_hsp}}: {n_prn_ok}")
    summary.append(f"Disrupted fraction among interpretable calls: {disrupted_frac:.3f}")
    summary.append("\nprn_call value counts (all rows):")
    for k, v in overall_counts.items():
        summary.append(f"- {k}: {int(v)}")

    (outdir / f"{args.prefix}_step3D_prn_trends_summary.txt").write_text("\n".join(summary) + "\n")

    print(f"Wrote Step3D tables into: {outdir}")


if __name__ == "__main__":
    main()
