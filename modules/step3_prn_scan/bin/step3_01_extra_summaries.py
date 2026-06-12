#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

import pandas as pd


MARKERS = ["prn", "ptxP_promoter", "fim2", "fim3"]


def norm(s: pd.Series) -> pd.Series:
    return s.astype(str).where(~s.isna(), "NA").replace({"nan": "NA", "None": "NA"})


def ensure_cols(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for c in cols:
        if c not in out.columns:
            out[c] = "NA"
        out[c] = norm(out[c])
    return out


def value_counts_table(df: pd.DataFrame, group_cols: list[str], value_col: str) -> pd.DataFrame:
    df = ensure_cols(df, group_cols + [value_col])

    if not group_cols:
        vc = df[value_col].value_counts(dropna=False).reset_index()
        vc.columns = [value_col, "n"]
        vc["group_n"] = int(len(df))
        vc["frac"] = vc["n"] / vc["group_n"]
        vc = vc.sort_values(["n"], ascending=[False])
        return vc

    grp = df.groupby(group_cols + [value_col], dropna=False).size().reset_index(name="n")
    totals = df.groupby(group_cols, dropna=False).size().reset_index(name="group_n")
    merged = grp.merge(totals, on=group_cols, how="left")
    merged["frac"] = merged["n"] / merged["group_n"]
    merged = merged.sort_values(group_cols + ["n"], ascending=[True] * len(group_cols) + [False])
    return merged


def simpson_diversity(counts: pd.Series) -> float:
    total = counts.sum()
    if total <= 0:
        return float("nan")
    p2 = ((counts / total) ** 2).sum()
    return float(1.0 - p2)


def shannon_diversity(counts: pd.Series) -> float:
    total = counts.sum()
    if total <= 0:
        return float("nan")
    ps = counts / total
    out = 0.0
    for p in ps:
        if p <= 0:
            continue
        out -= float(p) * math.log(float(p))
    return out


def marker_diversity(df: pd.DataFrame, marker: str, group_cols: list[str]) -> pd.DataFrame:
    status_col = f"marker_status_{marker}"
    hash_col = f"marker_{marker}"

    if status_col not in df.columns or hash_col not in df.columns:
        return pd.DataFrame()

    sub = df.copy()
    sub = ensure_cols(sub, group_cols + [status_col, hash_col])
    sub = sub[sub[status_col] == "ok"].copy()

    if sub.empty:
        cols = group_cols + [
            "n_total",
            "n_alleles",
            "top_allele",
            "top_n",
            "top_frac",
            "simpson",
            "shannon",
        ]
        return pd.DataFrame(columns=cols)

    rows = []
    if not group_cols:
        groups = [((), sub)]
    else:
        groups = list(sub.groupby(group_cols, dropna=False))

    for key, g in groups:
        vc = g[hash_col].value_counts(dropna=False)
        top_allele = str(vc.index[0]) if len(vc) else "NA"
        top_n = int(vc.iloc[0]) if len(vc) else 0
        total = int(vc.sum())
        rows.append(
            {
                **({} if not group_cols else {c: v for c, v in zip(group_cols, key)}),
                "n_total": total,
                "n_alleles": int(g[hash_col].nunique(dropna=False)),
                "top_allele": top_allele,
                "top_n": top_n,
                "top_frac": (top_n / total) if total else float("nan"),
                "simpson": simpson_diversity(vc),
                "shannon": shannon_diversity(vc),
            }
        )

    out = pd.DataFrame(rows)
    if group_cols:
        out = out.sort_values(group_cols)
    else:
        out = out.sort_values(["n_total"], ascending=[False])
    return out


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def write_txt(lines: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="Step3A: extra cross-tabs and diversity summaries")
    ap.add_argument("--table", required=True, help="Input merged table (bp_step2 output)")
    ap.add_argument("--outdir", default="outputs", help="Output directory")
    ap.add_argument("--prefix", default="bp", help="Output prefix")
    args = ap.parse_args()

    table = Path(args.table)
    outdir = Path(args.outdir)
    prefix = args.prefix

    df = pd.read_csv(table, sep="\t", dtype=str)

    # Basic cross-tabs
    write_tsv(value_counts_table(df, ["mlst_st"], "country"), outdir / f"{prefix}_country_by_mlst_st.tsv")
    write_tsv(value_counts_table(df, ["country"], "mlst_st"), outdir / f"{prefix}_mlst_st_by_country.tsv")
    write_tsv(value_counts_table(df, ["year"], "mlst_st"), outdir / f"{prefix}_mlst_st_by_year_step3.tsv")

    # 23S call by ST and by (country, year)
    if "23s_A2047G_call" in df.columns:
        write_tsv(value_counts_table(df, ["mlst_st"], "23s_A2047G_call"), outdir / f"{prefix}_23s_call_by_mlst_st.tsv")
        write_tsv(
            value_counts_table(df, ["country", "year"], "23s_A2047G_call"),
            outdir / f"{prefix}_23s_call_by_country_year.tsv",
        )

    # Marker status by year and by ST
    for m in MARKERS:
        status_col = f"marker_status_{m}"
        if status_col in df.columns:
            write_tsv(value_counts_table(df, ["year"], status_col), outdir / f"{prefix}_{status_col}_by_year_step3.tsv")
            write_tsv(value_counts_table(df, ["mlst_st"], status_col), outdir / f"{prefix}_{status_col}_by_mlst_st.tsv")

    # Marker allele counts by year (only ok)
    for m in MARKERS:
        status_col = f"marker_status_{m}"
        hash_col = f"marker_{m}"
        if status_col not in df.columns or hash_col not in df.columns:
            continue
        sub = df.copy()
        sub = ensure_cols(sub, ["year", status_col, hash_col])
        sub = sub[sub[status_col] == "ok"].copy()
        if sub.empty:
            continue
        tab = value_counts_table(sub, ["year"], hash_col)
        write_tsv(tab, outdir / f"{prefix}_marker_{m}_alleles_by_year_step3.tsv")

    # Diversity metrics for each marker (overall / by year / by country)
    for m in MARKERS:
        overall = marker_diversity(df, m, [])
        if not overall.empty:
            write_tsv(overall, outdir / f"{prefix}_marker_{m}_diversity_overall.tsv")
        by_year = marker_diversity(df, m, ["year"])
        if not by_year.empty:
            write_tsv(by_year, outdir / f"{prefix}_marker_{m}_diversity_by_year.tsv")
        by_country = marker_diversity(df, m, ["country"])
        if not by_country.empty:
            write_tsv(by_country, outdir / f"{prefix}_marker_{m}_diversity_by_country.tsv")

    # Simple text summary for quick reading
    lines: list[str] = []
    lines.append("Step3A summary")
    lines.append("")
    lines.append(f"Rows: {len(df)}")
    for c in ["country", "year", "month_key", "mlst_scheme", "mlst_st"]:
        if c in df.columns:
            lines.append(f"Unique {c}: {df[c].nunique(dropna=True)}")
    lines.append("")
    for m in MARKERS:
        status_col = f"marker_status_{m}"
        if status_col in df.columns:
            s = df[status_col].fillna("NA").astype(str)
            lines.append(f"{m} status counts: {dict(s.value_counts())}")
    write_txt(lines, outdir / f"{prefix}_step3A_summary.txt")


if __name__ == "__main__":
    main()
