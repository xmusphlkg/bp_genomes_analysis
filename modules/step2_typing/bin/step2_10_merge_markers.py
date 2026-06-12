#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser(description="Merge marker alleles + 23S calls into bp_qc_merged_mlst.tsv")
    ap.add_argument("--qc-mlst", required=True, help="QC+MLST table (outputs/bp_qc_merged_mlst.tsv)")
    ap.add_argument("--marker-alleles", required=True, help="Marker alleles long TSV from step2_08")
    ap.add_argument("--out", required=True, help="Output TSV")
    ap.add_argument("--23s-summary", default=None, help="Optional 23S per-genome summary TSV from step2_09")
    args = ap.parse_args()

    qc = Path(args.qc_mlst)
    if not qc.exists() or qc.stat().st_size == 0:
        raise SystemExit(f"ERROR: qc+mlst TSV missing or empty: {qc}")
    ma = Path(args.marker_alleles)
    if not ma.exists() or ma.stat().st_size == 0:
        raise SystemExit(f"ERROR: marker alleles TSV missing or empty: {ma}")

    df = pd.read_csv(qc, sep="\t", dtype=str)
    df_ma = pd.read_csv(ma, sep="\t", dtype=str)
    if "genome_fasta_path" not in df.columns:
        raise SystemExit("ERROR: qc table missing genome_fasta_path")
    if "genome_fasta_path" not in df_ma.columns or "marker" not in df_ma.columns:
        raise SystemExit("ERROR: marker allele table missing genome_fasta_path and/or marker columns")

    # Keep only the columns we want to merge, and only ok/below_threshold (ignore errors for now)
    keep_cols = [
        "genome_fasta_path",
        "marker",
        "status",
        "allele_hash",
        "allele_len",
        "pident",
        "qcov_pct",
    ]
    for c in keep_cols:
        if c not in df_ma.columns:
            df_ma[c] = ""
    df_ma = df_ma[keep_cols]

    # Pivot marker -> columns
    def pivot(col: str, prefix: str) -> pd.DataFrame:
        wide = df_ma.pivot_table(index="genome_fasta_path", columns="marker", values=col, aggfunc="first")
        wide.columns = [f"{prefix}{m}" for m in wide.columns]
        return wide.reset_index()

    df_hash = pivot("allele_hash", "marker_")
    df_len = pivot("allele_len", "marker_len_")
    df_pident = pivot("pident", "marker_pident_")
    df_qcov = pivot("qcov_pct", "marker_qcov_")
    df_status = pivot("status", "marker_status_")

    out = df.merge(df_hash, on="genome_fasta_path", how="left")
    out = out.merge(df_len, on="genome_fasta_path", how="left")
    out = out.merge(df_pident, on="genome_fasta_path", how="left")
    out = out.merge(df_qcov, on="genome_fasta_path", how="left")
    out = out.merge(df_status, on="genome_fasta_path", how="left")

    if args.__dict__.get("23s_summary"):
        s23 = Path(args.__dict__["23s_summary"])
        if s23.exists() and s23.stat().st_size > 0:
            df23 = pd.read_csv(s23, sep="\t", dtype=str)
            if "genome_fasta_path" in df23.columns and "call" in df23.columns:
                df23 = df23[["genome_fasta_path", "call"]].rename(columns={"call": "23s_A2047G_call"})
                out = out.merge(df23, on="genome_fasta_path", how="left")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    print("[MergeMarkers] rows:", len(out))
    print("Wrote:", out_path)


if __name__ == "__main__":
    main()
