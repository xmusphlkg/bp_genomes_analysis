#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def norm(s: pd.Series) -> pd.Series:
    out = s.astype(str)
    out = out.where(~s.isna(), None)
    out = out.str.strip()
    out = out.replace({"": None, "nan": None, "None": None})
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Summarize marker allele hashes from merged table")
    ap.add_argument("--merged", required=True, help="Merged table with marker_* columns (e.g. outputs/bp_qc_merged_mlst_markers.tsv)")
    ap.add_argument(
        "--genotype-manifest",
        default=None,
        help="Optional standardized genotype manifest from step2_14_harmonize_typing.py",
    )
    ap.add_argument("--outdir", default="outputs", help="Output directory")
    ap.add_argument("--prefix", default="bp", help="Output prefix")
    ap.add_argument(
        "--markers",
        default=None,
        help="Comma-separated marker names to summarize (default: auto-detect from marker_ columns)",
    )
    args = ap.parse_args()

    merged = Path(args.merged)
    if not merged.exists() or merged.stat().st_size == 0:
        raise SystemExit(f"ERROR: merged TSV missing or empty: {merged}")

    df = pd.read_csv(merged, sep="\t", dtype=str)
    for col in ["country", "year", "month_key"]:
        if col in df.columns:
            df[col] = norm(df[col])

    genotype_manifest = None
    genotype_manifest_path = Path(args.genotype_manifest) if args.genotype_manifest else merged.parent / f"{args.prefix}_genotype_manifest.tsv"
    if genotype_manifest_path.exists() and genotype_manifest_path.stat().st_size > 0:
        genotype_manifest = pd.read_csv(genotype_manifest_path, sep="\t", dtype=str)
        for col in ["assembly_accession", "country", "year", "month_key"]:
            if col in genotype_manifest.columns:
                genotype_manifest[col] = norm(genotype_manifest[col])

    # marker columns are marker_<name>
    marker_cols = [c for c in df.columns if c.startswith("marker_") and not c.startswith("marker_len_") and not c.startswith("marker_pident_") and not c.startswith("marker_qcov_") and not c.startswith("marker_status_")]
    if args.markers:
        wanted = [m.strip() for m in str(args.markers).split(",") if m.strip()]
        marker_cols = [f"marker_{m}" for m in wanted if f"marker_{m}" in df.columns]
    if not marker_cols:
        raise SystemExit("ERROR: no marker_<name> allele hash columns found")

    outdir = Path(args.outdir)
    prefix = str(args.prefix)

    for col in marker_cols:
        m = col.replace("marker_", "", 1)
        sub = df[[col]].copy()
        sub[col] = norm(sub[col])

        overall = (
            sub.dropna(subset=[col])
            .groupby(col, dropna=False)
            .size()
            .reset_index(name="n")
            .rename(columns={col: "allele_hash"})
            .sort_values(["n", "allele_hash"], ascending=[False, True])
        )
        write_tsv(overall, outdir / f"{prefix}_marker_{m}_allele_counts.tsv")

        if "country" in df.columns:
            by_cty = (
                df.dropna(subset=["country", col])
                .groupby(["country", col], dropna=False)
                .size()
                .reset_index(name="n")
                .rename(columns={col: "allele_hash"})
                .sort_values(["country", "n", "allele_hash"], ascending=[True, False, True])
            )
            write_tsv(by_cty, outdir / f"{prefix}_marker_{m}_by_country.tsv")

        if "year" in df.columns:
            by_year = (
                df.dropna(subset=["year", col])
                .groupby(["year", col], dropna=False)
                .size()
                .reset_index(name="n")
                .rename(columns={col: "allele_hash"})
                .sort_values(["year", "n", "allele_hash"], ascending=[True, False, True])
            )
            write_tsv(by_year, outdir / f"{prefix}_marker_{m}_by_year.tsv")

        if "month_key" in df.columns:
            by_mk = (
                df.dropna(subset=["month_key", col])
                .groupby(["month_key", col], dropna=False)
                .size()
                .reset_index(name="n")
                .rename(columns={col: "allele_hash"})
                .sort_values(["month_key", "n", "allele_hash"], ascending=[True, False, True])
            )
            write_tsv(by_mk, outdir / f"{prefix}_marker_{m}_by_month_key.tsv")

    if genotype_manifest is not None:
        label_specs = {
            "ptxP_promoter": "ptxP_label",
            "fim3": "fim3_label",
            "fhaB2400_5550": "fhaB2400_5550_label",
            "23s": "marker_23s_status",
        }
        for marker_name, label_col in label_specs.items():
            if label_col not in genotype_manifest.columns:
                continue
            sub = genotype_manifest[[label_col]].copy()
            sub[label_col] = norm(sub[label_col])
            overall = (
                sub.dropna(subset=[label_col])
                .groupby(label_col, dropna=False)
                .size()
                .reset_index(name="n")
                .rename(columns={label_col: "label"})
                .sort_values(["n", "label"], ascending=[False, True])
            )
            write_tsv(overall, outdir / f"{prefix}_marker_{marker_name}_label_counts.tsv")

            for strat_col, suffix in [("country", "by_country"), ("year", "by_year"), ("month_key", "by_month_key")]:
                if strat_col not in genotype_manifest.columns:
                    continue
                stratified = (
                    genotype_manifest.dropna(subset=[strat_col, label_col])
                    .groupby([strat_col, label_col], dropna=False)
                    .size()
                    .reset_index(name="n")
                    .rename(columns={label_col: "label"})
                    .sort_values([strat_col, "n", "label"], ascending=[True, False, True])
                )
                write_tsv(stratified, outdir / f"{prefix}_marker_{marker_name}_{suffix}_labels.tsv")

    print("[MarkerSummary] markers:", ",".join(c.replace("marker_", "", 1) for c in marker_cols))


if __name__ == "__main__":
    main()
