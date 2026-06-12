#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def norm(s: pd.Series) -> pd.Series:
    return s.astype(str).where(~s.isna(), "NA").replace({"nan": "NA", "None": "NA"})


def value_counts(df: pd.DataFrame, group_cols: list[str], value_col: str) -> pd.DataFrame:
    df = df.copy()
    if value_col not in df.columns:
        df[value_col] = "NA"
    df[value_col] = norm(df[value_col])

    for c in group_cols:
        if c not in df.columns:
            df[c] = "NA"
        df[c] = norm(df[c])

    if not group_cols:
        vc = df[value_col].value_counts(dropna=False).reset_index()
        vc.columns = [value_col, "n"]
        vc["group_n"] = int(len(df))
        vc["frac"] = vc["n"] / vc["group_n"]
        return vc.sort_values(["n"], ascending=[False])

    grp = df.groupby(group_cols + [value_col], dropna=False).size().reset_index(name="n")
    totals = df.groupby(group_cols, dropna=False).size().reset_index(name="group_n")
    out = grp.merge(totals, on=group_cols, how="left")
    out["frac"] = out["n"] / out["group_n"]
    out = out.sort_values(group_cols + ["n"], ascending=[True] * len(group_cols) + [False])
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Summaries for Step3F breakpoint evidence categories")
    ap.add_argument("--evidence", required=True, help="Evidence TSV from step3_50")
    ap.add_argument("--outdir", default="outputs", help="Output directory")
    ap.add_argument("--prefix", default="bp", help="Prefix")
    args = ap.parse_args()

    ev = pd.read_csv(Path(args.evidence), sep="\t", dtype=str)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if "bp_category" not in ev.columns:
        raise SystemExit("ERROR: evidence table missing bp_category")

    value_counts(ev, [], "bp_category").to_csv(outdir / f"{args.prefix}_prn_bp_category_counts.tsv", sep="\t", index=False)

    for grp, name in [(["year"], "by_year"), (["mlst_st"], "by_mlst_st"), (["year", "mlst_st"], "by_year_mlst_st"), (["country"], "by_country")]:
        if all(c in ev.columns for c in grp):
            value_counts(ev, grp, "bp_category").to_csv(outdir / f"{args.prefix}_prn_bp_category_{name}.tsv", sep="\t", index=False)

    # Gap-length distribution for insertion-like calls (often indicates IS-sized insertions)
    if "bp_max_subject_gap" in ev.columns:
        sub = ev[ev["bp_category"] == "insertion_like"].copy()
        sub["bp_max_subject_gap"] = norm(sub["bp_max_subject_gap"])
        sub = sub[sub["bp_max_subject_gap"] != "NA"].copy()
        sub["gap"] = pd.to_numeric(sub["bp_max_subject_gap"], errors="coerce")
        sub = sub[sub["gap"].notna()].copy()
        if len(sub) > 0:
            vc = sub["gap"].value_counts().reset_index()
            vc.columns = ["bp_max_subject_gap", "n"]
            vc.to_csv(outdir / f"{args.prefix}_prn_bp_subject_gap_counts_insertion_like.tsv", sep="\t", index=False)

            # Write a tiny text summary for paper drafting
            gap_med = float(sub["gap"].median())
            gap_p25 = float(sub["gap"].quantile(0.25))
            gap_p75 = float(sub["gap"].quantile(0.75))
            gap_mode = float(vc.iloc[0]["bp_max_subject_gap"]) if len(vc) else float("nan")
            in_is_window = float(((sub["gap"] >= 900) & (sub["gap"] <= 1200)).mean())

            txt = []
            txt.append(f"insertion_like n: {len(sub)}")
            txt.append(f"bp_max_subject_gap median: {gap_med:.1f} (IQR {gap_p25:.1f}-{gap_p75:.1f})")
            txt.append(f"bp_max_subject_gap mode (most frequent): {gap_mode:.0f}")
            txt.append(f"fraction with gap in [900,1200] bp: {in_is_window:.3f}")
            (outdir / f"{args.prefix}_prn_bp_subject_gap_summary.txt").write_text("\n".join(txt) + "\n")


if __name__ == "__main__":
    main()
