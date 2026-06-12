#!/usr/bin/env python3

import argparse
import os
import random
import re
from pathlib import Path

import pandas as pd


def norm(s: pd.Series) -> pd.Series:
    return s.astype(str).where(~s.isna(), "NA").replace({"nan": "NA", "None": "NA"})


def sanitize(s: str) -> str:
    s = re.sub(r"\s+", "_", s.strip())
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s)
    s = re.sub(r"_+", "_", s)
    return s.strip("_") or "NA"


def pick_sample_id(row: pd.Series) -> str:
    for c in ["genome_resolved_accession", "Current Accession", "Assembly Accession"]:
        if c in row.index:
            v = str(row[c])
            if v and v != "NA" and v.lower() != "nan":
                return sanitize(v)
    return "sample"


def file_exists(path_str: str) -> bool:
    try:
        p = Path(path_str)
        return p.exists() and p.is_file() and p.stat().st_size > 0
    except Exception:
        return False


def stratified_sample(df: pd.DataFrame, group_cols: list[str], per_group: int, max_total: int, seed: int) -> pd.DataFrame:
    if max_total <= 0:
        return df

    rng = random.Random(seed)

    df = df.copy()
    for c in group_cols:
        if c not in df.columns:
            df[c] = "NA"
        df[c] = norm(df[c])

    chunks = []
    for _, g in df.groupby(group_cols, dropna=False):
        if len(g) <= per_group:
            chunks.append(g)
        else:
            idx = list(g.index)
            rng.shuffle(idx)
            chunks.append(g.loc[idx[:per_group]])

    out = pd.concat(chunks, axis=0)
    if len(out) <= max_total:
        return out

    idx = list(out.index)
    rng.shuffle(idx)
    return out.loc[idx[:max_total]]


def main() -> None:
    ap = argparse.ArgumentParser(description="Prepare phylogeny manifest + symlinked genome set")
    ap.add_argument("--table", required=True, help="Input merged table (bp_step2 output)")
    ap.add_argument("--outdir", default="outputs", help="Output directory")
    ap.add_argument("--prefix", default="bp", help="Output prefix")
    ap.add_argument("--max-genomes", type=int, default=0, help="Max genomes to include (0=all)")
    ap.add_argument("--per-group", type=int, default=5, help="Per-group cap used with --max-genomes")
    ap.add_argument(
        "--group-cols",
        default="mlst_st,year",
        help="Comma-separated columns for stratified sampling when --max-genomes > 0 (default: mlst_st,year)",
    )
    ap.add_argument("--seed", type=int, default=1, help="Random seed for sampling (default: 1)")
    args = ap.parse_args()

    table = Path(args.table)
    outdir = Path(args.outdir)
    prefix = args.prefix

    df = pd.read_csv(table, sep="\t", dtype=str)

    required = ["genome_status", "genome_fasta_path"]
    for c in required:
        if c not in df.columns:
            raise SystemExit(f"ERROR: missing required column in table: {c}")

    df["genome_status"] = norm(df["genome_status"])
    df["genome_fasta_path"] = norm(df["genome_fasta_path"])

    df = df[df["genome_status"] == "ok"].copy()
    df = df[df["genome_fasta_path"] != "NA"].copy()

    # Verify paths exist
    df["_path_exists"] = df["genome_fasta_path"].map(file_exists)
    missing = df[~df["_path_exists"]].copy()
    df = df[df["_path_exists"]].copy()

    group_cols = [c.strip() for c in str(args.group_cols).split(",") if c.strip()]
    if args.max_genomes and args.max_genomes > 0:
        df = stratified_sample(df, group_cols=group_cols, per_group=int(args.per_group), max_total=int(args.max_genomes), seed=int(args.seed))

    # Create symlink directory
    genomes_dir = outdir / f"{prefix}_phylo_genomes"
    genomes_dir.mkdir(parents=True, exist_ok=True)

    used: set[str] = set()
    sample_ids: list[str] = []
    genome_paths_out: list[str] = []

    for _, row in df.iterrows():
        sid = pick_sample_id(row)
        base = sid
        i = 1
        while sid in used:
            i += 1
            sid = f"{base}__{i}"
        used.add(sid)

        src = Path(str(row["genome_fasta_path"]))
        ext = src.suffix if src.suffix else ".fna"
        dest = genomes_dir / f"{sid}{ext}"
        if not dest.exists():
            os.symlink(src.resolve(), dest)

        sample_ids.append(sid)
        genome_paths_out.append(str(dest))

    out_manifest = df.copy()
    out_manifest.insert(0, "sample_id", sample_ids)
    out_manifest.insert(1, "genome_path", genome_paths_out)

    keep_cols = [
        "sample_id",
        "genome_path",
        "country",
        "year",
        "month_key",
        "mlst_scheme",
        "mlst_st",
        "mlst_profile",
        "marker_prn",
        "marker_ptxP_promoter",
        "marker_fim2",
        "marker_fim3",
        "marker_status_prn",
        "marker_status_ptxP_promoter",
        "marker_status_fim2",
        "marker_status_fim3",
        "23s_A2047G_call",
        "Current Accession",
        "Assembly BioSample Accession",
    ]

    for c in keep_cols:
        if c not in out_manifest.columns:
            out_manifest[c] = "NA"
        out_manifest[c] = norm(out_manifest[c])

    out_manifest = out_manifest[keep_cols].copy()

    outdir.mkdir(parents=True, exist_ok=True)

    manifest_path = outdir / f"{prefix}_phylo_manifest.tsv"
    out_manifest.to_csv(manifest_path, sep="\t", index=False)

    if not missing.empty:
        missing_path = outdir / f"{prefix}_phylo_missing_paths.tsv"
        missing[["genome_fasta_path", "genome_status"]].to_csv(missing_path, sep="\t", index=False)

    print(f"Wrote: {manifest_path}")
    print(f"Genomes dir: {genomes_dir}")
    print(f"N genomes: {len(out_manifest)}")


if __name__ == "__main__":
    main()
