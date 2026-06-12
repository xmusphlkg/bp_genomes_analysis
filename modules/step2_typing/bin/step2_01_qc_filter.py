#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def find_col(df: pd.DataFrame, patterns: list[str]) -> str | None:
    cols = list(df.columns)
    low = [c.lower() for c in cols]
    for p in patterns:
        p = p.lower()
        for c, cl in zip(cols, low):
            if p in cl:
                return c
    return None


def to_int_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").astype("Int64")


def to_float_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def main() -> None:
    ap = argparse.ArgumentParser(description="QC filter for B. pertussis assemblies using step1 metadata stats")
    ap.add_argument("--metadata", required=True, help="Input metadata CSV from step1 (bp_metadata_clean.csv)")
    ap.add_argument("--out-metadata", required=True, help="Output QC-passing metadata CSV")
    ap.add_argument("--out-accessions", required=True, help="Output accession list (one per line)")
    ap.add_argument("--out-rejects", default=None, help="Optional output rejects CSV")

    ap.add_argument("--min-len", type=int, default=3_500_000, help="Min total length (bp)")
    ap.add_argument("--max-len", type=int, default=4_500_000, help="Max total length (bp)")
    ap.add_argument("--max-contigs", type=int, default=200, help="Max number of contigs")
    ap.add_argument("--min-n50", type=int, default=20_000, help="Min contig N50")
    ap.add_argument(
        "--keep-missing-stats",
        action="store_true",
        help="If set, keep rows with missing assembly stats (they will be marked qc_reason=missing_stats)",
    )
    args = ap.parse_args()

    meta_path = Path(args.metadata)
    if not meta_path.exists() or meta_path.stat().st_size == 0:
        raise SystemExit(f"ERROR: metadata CSV missing or empty: {meta_path}")

    df = pd.read_csv(meta_path, dtype=str)

    col_acc = find_col(df, ["current accession", "assembly accession", "accession"])
    col_len = find_col(df, ["total sequence length", "stats total sequence length"])
    col_gc = find_col(df, ["gc percent", "stats gc"])
    col_contigs = find_col(df, ["number of contigs", "stats number of contigs"])
    col_n50 = find_col(df, ["contig n50", "stats contig n50", "n50"])

    if not col_acc:
        raise SystemExit("ERROR: could not detect accession column")
    if not col_len or not col_contigs or not col_n50:
        raise SystemExit("ERROR: could not detect length/contigs/n50 columns")

    df["_acc"] = df[col_acc].fillna("").astype(str).str.strip()
    df["_len"] = to_int_series(df[col_len])
    df["_contigs"] = to_int_series(df[col_contigs])
    df["_n50"] = to_int_series(df[col_n50])
    if col_gc:
        df["_gc"] = to_float_series(df[col_gc])

    missing_stats = df[["_len", "_contigs", "_n50"]].isna().any(axis=1)

    pass_len = (df["_len"] >= args.min_len) & (df["_len"] <= args.max_len)
    pass_contigs = df["_contigs"] <= args.max_contigs
    pass_n50 = df["_n50"] >= args.min_n50

    qc_pass = pass_len & pass_contigs & pass_n50

    df["qc_pass"] = qc_pass
    df["qc_reason"] = "pass"

    df.loc[missing_stats, "qc_reason"] = "missing_stats"
    df.loc[~missing_stats & ~pass_len, "qc_reason"] = "length_out_of_range"
    df.loc[~missing_stats & pass_len & ~pass_contigs, "qc_reason"] = "too_many_contigs"
    df.loc[~missing_stats & pass_len & pass_contigs & ~pass_n50, "qc_reason"] = "low_n50"

    if args.keep_missing_stats:
        df.loc[missing_stats, "qc_pass"] = True

    # Keep only rows with a real accession
    df = df[df["_acc"] != ""].copy()

    passed = df[df["qc_pass"]].copy()
    rejected = df[~df["qc_pass"]].copy()

    out_meta = Path(args.out_metadata)
    out_meta.parent.mkdir(parents=True, exist_ok=True)
    passed.drop(columns=[c for c in ["_acc", "_len", "_contigs", "_n50", "_gc"] if c in passed.columns]).to_csv(
        out_meta, index=False
    )

    out_acc = Path(args.out_accessions)
    out_acc.parent.mkdir(parents=True, exist_ok=True)
    passed[col_acc].dropna().astype(str).str.strip().drop_duplicates().to_csv(out_acc, index=False, header=False)

    if args.out_rejects:
        out_rej = Path(args.out_rejects)
        out_rej.parent.mkdir(parents=True, exist_ok=True)
        rejected.drop(columns=[c for c in ["_acc", "_len", "_contigs", "_n50", "_gc"] if c in rejected.columns]).to_csv(
            out_rej, index=False
        )

    print("[QC] input rows:", len(df))
    print("[QC] pass:", len(passed))
    print("[QC] reject:", len(rejected))
    print("[QC] reasons (rejects):")
    if len(rejected) == 0:
        print("  (none)")
    else:
        vc = rejected["qc_reason"].value_counts(dropna=False)
        for k, v in vc.items():
            print(f"  - {k}: {int(v)}")

    print("Wrote:")
    print(f"  - {out_meta}")
    print(f"  - {out_acc}")
    if args.out_rejects:
        print(f"  - {args.out_rejects}")


if __name__ == "__main__":
    main()
