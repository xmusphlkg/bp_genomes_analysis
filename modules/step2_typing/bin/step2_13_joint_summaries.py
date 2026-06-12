#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def _normalize_series(s: pd.Series) -> pd.Series:
    return s.astype(str).where(~s.isna(), "NA").replace({"nan": "NA", "None": "NA"})


def write_value_counts(
    df: pd.DataFrame,
    group_cols: list[str],
    value_col: str,
    out_path: Path,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    tmp = df.copy()
    tmp[value_col] = _normalize_series(tmp.get(value_col, pd.Series(["NA"] * len(tmp))))

    for c in group_cols:
        if c not in tmp.columns:
            tmp[c] = "NA"
        tmp[c] = _normalize_series(tmp[c])

    if not group_cols:
        counts = tmp[value_col].value_counts(dropna=False).reset_index()
        counts.columns = [value_col, "n"]
        counts["group_n"] = int(len(tmp))
        counts["frac"] = counts["n"] / counts["group_n"]
        counts = counts.sort_values(["n"], ascending=[False])
        counts.to_csv(out_path, sep="\t", index=False)
        return

    grp = tmp.groupby(group_cols + [value_col], dropna=False).size().reset_index(name="n")
    totals = tmp.groupby(group_cols, dropna=False).size().reset_index(name="group_n")
    merged = grp.merge(totals, on=group_cols, how="left")
    merged["frac"] = merged["n"] / merged["group_n"]

    merged = merged.sort_values(group_cols + ["n"], ascending=[True] * len(group_cols) + [False])
    merged.to_csv(out_path, sep="\t", index=False)


def top_genotype_combos(
    df: pd.DataFrame,
    out_path: Path,
    top_n: int,
    group_cols: list[str] | None = None,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    cols = {
        "mlst_st": "mlst_st",
        "prn": "marker_prn",
        "ptxP_promoter": "marker_ptxP_promoter",
        "fim2": "marker_fim2",
        "fim3": "marker_fim3",
        "23s_call": "23s_A2047G_call",
    }
    tmp = pd.DataFrame({k: _normalize_series(df.get(v, "NA")) for k, v in cols.items()})

    if group_cols:
        for c in group_cols:
            tmp[c] = _normalize_series(df.get(c, "NA"))
        keys = group_cols + list(cols.keys())
    else:
        keys = list(cols.keys())

    counts = tmp.groupby(keys, dropna=False).size().reset_index(name="n")

    if group_cols:
        totals = tmp.groupby(group_cols, dropna=False).size().reset_index(name="group_n")
        counts = counts.merge(totals, on=group_cols, how="left")
        counts["frac"] = counts["n"] / counts["group_n"]
        counts = counts.sort_values(group_cols + ["n"], ascending=[True] * len(group_cols) + [False])
        counts = counts.groupby(group_cols, dropna=False).head(top_n)
    else:
        counts["group_n"] = len(tmp)
        counts["frac"] = counts["n"] / counts["group_n"]
        counts = counts.sort_values(["n"], ascending=[False]).head(top_n)

    counts.to_csv(out_path, sep="\t", index=False)


def write_summary_txt(df: pd.DataFrame, out_path: Path, ref_23s_fa: Path | None, query_pos: int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    lines.append("Pertussis genomic summary (Step7)")
    lines.append("")
    lines.append(f"Rows (QC-pass assemblies): {len(df)}")

    if "Current Accession" in df.columns:
        lines.append(f"Unique current accessions: {df['Current Accession'].nunique(dropna=True)}")
    if "Assembly BioSample Accession" in df.columns:
        lines.append(f"Unique BioSamples: {df['Assembly BioSample Accession'].nunique(dropna=True)}")

    for c in ["country", "year", "month_key"]:
        if c in df.columns:
            lines.append(f"Unique {c}: {df[c].nunique(dropna=True)}")

    lines.append("")

    # 23S reference base sanity
    if ref_23s_fa and ref_23s_fa.exists():
        seq = ""
        for ln in ref_23s_fa.read_text(encoding="utf-8", errors="ignore").splitlines():
            if ln.startswith(">"):
                continue
            seq += ln.strip()
        seq = seq.upper().replace(" ", "")
        base = seq[query_pos - 1] if len(seq) >= query_pos else "?"
        lines.append(f"23S reference: {ref_23s_fa}")
        lines.append(f"23S query position used: {query_pos} (reference base={base})")
        lines.append("Note: the column name '23s_A2047G_call' reflects the intended mutation; interpret relative to the reference base above.")
        lines.append("")

    # Status rates
    def rate(col: str) -> str:
        if col not in df.columns:
            return "NA"
        s = df[col].fillna("NA").astype(str)
        denom = len(s)
        ok = (s == "ok").sum()
        na = (s == "NA").sum()
        below = (s == "below_threshold").sum()
        return f"ok={ok/denom:.3f} below_threshold={below/denom:.3f} NA={na/denom:.3f}"

    lines.append("Marker status rates (fraction of assemblies):")
    for m in ["prn", "ptxP_promoter", "fim2", "fim3"]:
        lines.append(f"- {m}: {rate(f'marker_status_{m}')}")

    lines.append("")
    if "23s_A2047G_call" in df.columns:
        s = df["23s_A2047G_call"].fillna("NA").astype(str)
        counts = s.value_counts(dropna=False)
        lines.append("23S call counts:")
        for k, v in counts.items():
            lines.append(f"- {k}: {v} ({v/len(s):.3f})")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="Step7: joint summaries across MLST, markers, and 23S calls")
    ap.add_argument(
        "--table",
        default="outputs/bp_qc_merged_mlst_markers.tsv",
        help="Input merged table (default: outputs/bp_qc_merged_mlst_markers.tsv)",
    )
    ap.add_argument("--outdir", default="outputs", help="Output directory (default: outputs)")
    ap.add_argument("--prefix", default="bp", help="Output prefix (default: bp)")
    ap.add_argument("--top-n", type=int, default=50, help="Top N genotype combos to report (default: 50)")
    ap.add_argument(
        "--ref-23s",
        default="references/23S_rRNA.fasta",
        help="23S reference FASTA used for calling (for report text; default: references/23S_rRNA.fasta)",
    )
    ap.add_argument("--query-pos", type=int, default=2047, help="Query position used for 23S call (default: 2047)")
    args = ap.parse_args()

    table = Path(args.table)
    outdir = Path(args.outdir)

    df = pd.read_csv(table, sep="\t", dtype=str)

    # 23S call summaries
    write_value_counts(df, [], "23s_A2047G_call", outdir / f"{args.prefix}_23s_call_counts.tsv")
    write_value_counts(df, ["country"], "23s_A2047G_call", outdir / f"{args.prefix}_23s_call_by_country.tsv")
    write_value_counts(df, ["year"], "23s_A2047G_call", outdir / f"{args.prefix}_23s_call_by_year.tsv")
    write_value_counts(df, ["month_key"], "23s_A2047G_call", outdir / f"{args.prefix}_23s_call_by_month_key.tsv")

    # Marker status summaries (assembly-level)
    for m in ["prn", "ptxP_promoter", "fim2", "fim3"]:
        col = f"marker_status_{m}"
        if col not in df.columns:
            continue
        write_value_counts(df, [], col, outdir / f"{args.prefix}_{col}_counts.tsv")
        write_value_counts(df, ["country"], col, outdir / f"{args.prefix}_{col}_by_country.tsv")
        write_value_counts(df, ["year"], col, outdir / f"{args.prefix}_{col}_by_year.tsv")
        write_value_counts(df, ["month_key"], col, outdir / f"{args.prefix}_{col}_by_month_key.tsv")

    # Top genotype combos (overall + by year)
    top_genotype_combos(df, outdir / f"{args.prefix}_genotype_combo_top.tsv", top_n=int(args.top_n))
    top_genotype_combos(
        df,
        outdir / f"{args.prefix}_genotype_combo_top_by_year.tsv",
        top_n=int(args.top_n),
        group_cols=["year"],
    )

    # Text report
    write_summary_txt(df, outdir / f"{args.prefix}_analysis_summary.txt", Path(args.ref_23s), int(args.query_pos))


if __name__ == "__main__":
    main()
