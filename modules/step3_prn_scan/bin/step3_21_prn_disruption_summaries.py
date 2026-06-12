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
    ap = argparse.ArgumentParser(description="Summaries for prn disruption calls")
    ap.add_argument("--calls", required=True, help="Calls TSV from step3_20")
    ap.add_argument("--outdir", default="outputs", help="Output directory")
    ap.add_argument("--prefix", default="bp", help="Prefix")
    args = ap.parse_args()

    calls = pd.read_csv(Path(args.calls), sep="\t", dtype=str)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if "prn_call" not in calls.columns:
        raise SystemExit("ERROR: calls table missing prn_call")

    value_counts(calls, [], "prn_call").to_csv(outdir / f"{args.prefix}_prn_call_counts.tsv", sep="\t", index=False)

    for grp, name in [(["year"], "by_year"), (["country"], "by_country"), (["mlst_st"], "by_mlst_st"), (["year", "mlst_st"], "by_year_mlst_st")]:
        if all(c in calls.columns for c in grp):
            value_counts(calls, grp, "prn_call").to_csv(
                outdir / f"{args.prefix}_prn_call_{name}.tsv", sep="\t", index=False
            )


if __name__ == "__main__":
    main()
