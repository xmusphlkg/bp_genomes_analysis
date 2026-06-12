#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def norm_str(s: pd.Series) -> pd.Series:
    out = s.astype(str)
    out = out.where(~s.isna(), None)
    out = out.str.strip()
    out = out.replace({"": None, "nan": None, "None": None})
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Summarize MLST results from bp_qc_merged_mlst.tsv")
    ap.add_argument(
        "--qc-mlst",
        required=True,
        help="Merged QC+MLST TSV (outputs/bp_qc_merged_mlst.tsv)",
    )
    ap.add_argument(
        "--outdir",
        default="outputs",
        help="Output directory (default: outputs)",
    )
    ap.add_argument(
        "--prefix",
        default="bp",
        help="Output file prefix (default: bp)",
    )
    args = ap.parse_args()

    qc_path = Path(args.qc_mlst)
    if not qc_path.exists() or qc_path.stat().st_size == 0:
        raise SystemExit(f"ERROR: qc+mlst TSV missing or empty: {qc_path}")

    outdir = Path(args.outdir)
    prefix = str(args.prefix)

    df = pd.read_csv(qc_path, sep="\t", dtype=str)
    for col in ["mlst_scheme", "mlst_st", "country", "year", "month", "iso_week", "date_resolution"]:
        if col in df.columns:
            df[col] = norm_str(df[col])

    if "mlst_st" not in df.columns:
        raise SystemExit("ERROR: input does not contain mlst_st column; did you run step2_06_merge_mlst.py?")

    # Overall ST counts
    df_st = (
        df.dropna(subset=["mlst_st"])
        .groupby(["mlst_st"], dropna=False)
        .size()
        .reset_index(name="n")
        .sort_values(["n", "mlst_st"], ascending=[False, True])
    )
    write_tsv(df_st, outdir / f"{prefix}_mlst_st_counts.tsv")

    # Scheme counts (if present)
    if "mlst_scheme" in df.columns:
        df_scheme = (
            df.dropna(subset=["mlst_scheme"])
            .groupby(["mlst_scheme"], dropna=False)
            .size()
            .reset_index(name="n")
            .sort_values(["n", "mlst_scheme"], ascending=[False, True])
        )
        write_tsv(df_scheme, outdir / f"{prefix}_mlst_scheme_counts.tsv")

    # By country
    if "country" in df.columns:
        df_cty = (
            df.dropna(subset=["country", "mlst_st"])
            .groupby(["country", "mlst_st"], dropna=False)
            .size()
            .reset_index(name="n")
            .sort_values(["country", "n", "mlst_st"], ascending=[True, False, True])
        )
        write_tsv(df_cty, outdir / f"{prefix}_mlst_st_by_country.tsv")

    # By year
    if "year" in df.columns:
        df_year = (
            df.dropna(subset=["year", "mlst_st"])
            .groupby(["year", "mlst_st"], dropna=False)
            .size()
            .reset_index(name="n")
            .sort_values(["year", "n", "mlst_st"], ascending=[True, False, True])
        )
        write_tsv(df_year, outdir / f"{prefix}_mlst_st_by_year.tsv")

    # By month_key (best for time series)
    if "month_key" in df.columns:
        df["month_key"] = norm_str(df["month_key"])
        df_mk = (
            df.dropna(subset=["month_key", "mlst_st"])
            .groupby(["month_key", "mlst_st"], dropna=False)
            .size()
            .reset_index(name="n")
            .sort_values(["month_key", "n", "mlst_st"], ascending=[True, False, True])
        )
        write_tsv(df_mk, outdir / f"{prefix}_mlst_st_by_month_key.tsv")

    # A tiny run summary text
    n_rows = len(df)
    n_typed = int(df["mlst_st"].notna().sum())
    top_st = df_st.head(10)
    summary_lines = [
        f"rows\t{n_rows}",
        f"typed_rows\t{n_typed}",
        f"typed_fraction\t{(n_typed / n_rows) if n_rows else 0.0:.6f}",
        "top10_st\t" + ";".join(f"{r.mlst_st}:{r.n}" for r in top_st.itertuples(index=False)),
    ]
    (outdir / f"{prefix}_mlst_summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    print("[MLSTSummary] wrote:")
    for name in [
        f"{prefix}_mlst_st_counts.tsv",
        f"{prefix}_mlst_scheme_counts.tsv",
        f"{prefix}_mlst_st_by_country.tsv",
        f"{prefix}_mlst_st_by_year.tsv",
        f"{prefix}_mlst_st_by_month_key.tsv",
        f"{prefix}_mlst_summary.txt",
    ]:
        p = outdir / name
        if p.exists():
            print(" -", p)


if __name__ == "__main__":
    main()
