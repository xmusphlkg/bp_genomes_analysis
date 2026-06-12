#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def read_tsv(p: Path) -> pd.DataFrame:
    return pd.read_csv(p, sep="\t", dtype=str)


def fmt_pct(x: float) -> str:
    return f"{100.0 * x:.1f}%"


def safe_float(s: str) -> float | None:
    try:
        return float(s)
    except Exception:
        return None


def main() -> None:
    ap = argparse.ArgumentParser(description="Write a concise Markdown digest of Step3 outputs for plotting/writing")
    ap.add_argument("--outdir", default="outputs", help="Step3 outputs directory")
    ap.add_argument("--prefix", default="bp", help="Prefix")
    ap.add_argument("--out", default=None, help="Output markdown path (default: outputs/<prefix>_results_digest.md)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    prefix = args.prefix
    out_md = Path(args.out) if args.out else outdir / f"{prefix}_results_digest.md"

    # Inputs we expect
    p_summary = outdir / f"{prefix}_step3D_prn_trends_summary.txt"
    p_top_rank = outdir / f"{prefix}_prn_top_sts_rank.tsv"
    p_prn_23s = outdir / f"{prefix}_prn_vs_23s_overall.tsv"
    p_prn_ptxp = outdir / f"{prefix}_prn_vs_ptxP_overall.tsv"
    p_mashtree = outdir / f"{prefix}_mashtree.nwk"
    p_phylo_ann = outdir / f"{prefix}_phylo_annotations.tsv"

    p_bp_counts = outdir / f"{prefix}_prn_bp_category_counts.tsv"
    p_gap_summary = outdir / f"{prefix}_prn_bp_subject_gap_summary.txt"

    lines: list[str] = []
    lines.append(f"# Step3 Results Digest ({prefix})")
    lines.append("")
    lines.append("This file is auto-generated from Step3 outputs. Use it as a checklist for figures + results text.")
    lines.append("")

    # PRN summary
    if p_summary.exists():
        txt = p_summary.read_text().strip().splitlines()
        lines.append("## prn disruption (Step3C/3D)")
        lines.append("")
        for ln in txt[:30]:
            ln2 = ln.rstrip("\n")
            if not ln2.strip():
                continue
            if ln2.lstrip().startswith("-"):
                lines.append(ln2)
            else:
                lines.append(f"- {ln2}")
        lines.append("")
    else:
        lines.append("## prn disruption (Step3C/3D)")
        lines.append("")
        lines.append(f"- Missing summary file: {p_summary.name}")
        lines.append("")

    # Top ST table (take top 5 by n_total)
    if p_top_rank.exists():
        df = read_tsv(p_top_rank)
        # Ensure numeric ordering by n_total
        if "n_total" in df.columns:
            df["n_total_num"] = pd.to_numeric(df["n_total"], errors="coerce")
            df = df.sort_values(["n_total_num"], ascending=[False])
        top5 = df.head(5).copy()

        lines.append("## Top STs (where disrupted increases)")
        lines.append("")
        lines.append("Top STs by sample size with disrupted fractions (overall / early / late).")
        lines.append("")
        cols = [
            c
            for c in [
                "mlst_st",
                "n_total",
                "disrupted_multi_hsp",
                "intact",
                "n_early",
                "disrupted_multi_hsp_early",
                "n_late",
                "disrupted_multi_hsp_late",
                "delta_disrupted_late_minus_early",
            ]
            if c in top5.columns
        ]
        lines.append("```text")
        lines.append(top5[cols].to_string(index=False))
        lines.append("```")
        lines.append("")
        lines.append(f"Full table: `{p_top_rank.name}`")
        lines.append("")

    # prn vs 23S
    if p_prn_23s.exists():
        df = read_tsv(p_prn_23s)
        lines.append("## prn vs 23S A2047G (macrolide marker)")
        lines.append("")
        # show both prn_call rows
        show = df.head(20)
        lines.append("```text")
        lines.append(show.to_string(index=False))
        lines.append("```")
        lines.append("")
        # heuristic takeaway
        if {"prn_call", "23s_A2047G_call", "frac"}.issubset(df.columns):
            a2047g_rows = df[df["23s_A2047G_call"].str.contains("A2047G", na=False)]
            if len(a2047g_rows) == 0:
                lines.append("Takeaway: no obvious A2047G signal in current calls (mostly `other_base_T` / `no_call`).")
                lines.append("")

    # prn vs ptxP
    if p_prn_ptxp.exists():
        df = read_tsv(p_prn_ptxp)
        lines.append("## prn vs ptxP")
        lines.append("")
        # show top ptxP per prn_call
        if {"prn_call", "ptxP", "n"}.issubset(df.columns):
            df["n_num"] = pd.to_numeric(df["n"], errors="coerce")
            out_rows = []
            for prn_call in ["intact", "disrupted_multi_hsp"]:
                sub = df[df["prn_call"] == prn_call].sort_values(["n_num"], ascending=[False]).head(5)
                out_rows.append(sub)
            show = pd.concat(out_rows, ignore_index=True) if out_rows else df.head(20)
        else:
            show = df.head(20)
        lines.append("```text")
        lines.append(show.to_string(index=False))
        lines.append("```")
        lines.append("")
        lines.append("Note: `ptxP` here is the marker hash/allele ID from the Step6 marker extraction.")
        lines.append("")

    # Breakpoint evidence hardening
    if p_bp_counts.exists():
        lines.append("## Stronger evidence for prn disruption (Step3F)")
        lines.append("")
        df = read_tsv(p_bp_counts)
        lines.append("Breakpoint-evidence categories among disrupted calls that were re-scanned with detailed BLAST:")
        lines.append("")
        lines.append("```text")
        lines.append(df.to_string(index=False))
        lines.append("```")
        lines.append("")
        if p_gap_summary.exists():
            lines.append("Insertion-like gap length summary (subject gap; IS-sized insertions expected around ~1043bp):")
            lines.append("")
            lines.append("```text")
            lines.append(p_gap_summary.read_text().strip())
            lines.append("```")
            lines.append("")

    # Phylogeny: how to visualize
    lines.append("## Phylogeny visualization checklist (MashTree)")
    lines.append("")
    if p_mashtree.exists():
        lines.append(f"- Tree file: `{p_mashtree.name}`")
    else:
        lines.append(f"- Missing tree file: `{p_mashtree.name}`")

    if p_phylo_ann.exists():
        lines.append(f"- Annotation table aligned to tree sample IDs: `{p_phylo_ann.name}`")
        lines.append("  - Columns include: `sample_id`, `year`, `country`, `mlst_st`, `prn_call`, `23s_A2047G_call`, `ptxP`")
        lines.append("  - iTOL: upload the tree, then add datasets/metadata using this table to color strips by `prn_call` and label by `year`/`mlst_st`.")
    else:
        lines.append(f"- Missing phylogeny annotation table: `{p_phylo_ann.name}`")

    lines.append("")

    # Recommended next analysis
    lines.append("## What to do next (recommended)")
    lines.append("")
    lines.append("1) **Make 2 key figures**")
    lines.append("   - Figure A: disrupted fraction over time (overall + ST2). Use `bp_prn_trend_by_year_clean.tsv` and `bp_prn_trend_by_year_mlst_st_min20.tsv`.")
    lines.append("   - Figure B: MashTree colored by `prn_call` (and optionally `year`), using `bp_mashtree.nwk` + `bp_phylo_annotations.tsv`.")
    lines.append("")
    lines.append("2) **Upgrade from \"insertion-like\" to \"IS-confirmed\" (optional but makes it publication-grade)**")
    lines.append("   - Step3F shows most disrupted calls are within-contig with an IS-sized gap (~1043bp), arguing against simple assembly fragmentation.")
    lines.append("   - Gap sequences are already exported: `bp_prn_insertion_gap_plus_flanks.fasta` (+ metadata TSV).")
    lines.append("   - Next upgrade would BLAST these sequences against an IS reference set (e.g., IS481/IS1002/IS1663) and report top hits, or map reads (if available).")
    lines.append("")
    lines.append("3) **If you need a publication-grade phylogeny**")
    lines.append("   - Build a core-genome SNP tree (Snippy/Parsnp + IQ-TREE) and optionally date-scale it (TreeTime) when collection dates are precise enough.")
    lines.append("")

    out_md.write_text("\n".join(lines) + "\n")
    print(f"Wrote: {out_md}")


if __name__ == "__main__":
    main()
