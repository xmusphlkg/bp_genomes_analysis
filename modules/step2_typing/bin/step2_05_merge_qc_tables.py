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


def alt_accessions(acc: str) -> list[str]:
    acc = str(acc).strip()
    if not acc:
        return []
    out: list[str] = [acc]
    if acc.startswith("GCF_"):
        out.append("GCA_" + acc[len("GCF_") :])
    elif acc.startswith("GCA_"):
        out.append("GCF_" + acc[len("GCA_") :])
    # de-dupe, preserve order
    return list(dict.fromkeys(out))


def main() -> None:
    ap = argparse.ArgumentParser(description="Merge QC metadata with local genome FASTA paths")
    ap.add_argument("--metadata", required=True, help="QC metadata CSV from step2_01 (e.g. outputs/bp_metadata_qc.csv)")
    ap.add_argument("--genome-paths", required=True, help="Genome paths TSV from step2_02 (e.g. outputs/bp_genome_paths.tsv)")
    ap.add_argument("--out", required=True, help="Output merged TSV")
    ap.add_argument("--out-missing", default=None, help="Optional: write missing accessions (post-merge)")
    args = ap.parse_args()

    meta_path = Path(args.metadata)
    paths_path = Path(args.genome_paths)

    if not meta_path.exists() or meta_path.stat().st_size == 0:
        raise SystemExit(f"ERROR: metadata missing or empty: {meta_path}")
    if not paths_path.exists() or paths_path.stat().st_size == 0:
        raise SystemExit(f"ERROR: genome paths TSV missing or empty: {paths_path}")

    df_meta = pd.read_csv(meta_path, dtype=str)
    df_paths = pd.read_csv(paths_path, sep="\t", dtype=str)

    col_current = find_col(df_meta, ["current accession", "assembly accession"])
    col_asm = find_col(df_meta, ["assembly accession"])
    if not col_current:
        raise SystemExit("ERROR: could not find accession column in metadata")

    for required in ["input_accession", "status", "resolved_accession", "fasta_path"]:
        if required not in df_paths.columns:
            raise SystemExit(f"ERROR: genome paths TSV missing column: {required}")

    # Build lookup from genome paths
    lookup: dict[str, dict[str, str]] = {}
    for _, r in df_paths.iterrows():
        key = str(r.get("input_accession", "") or "").strip()
        if not key:
            continue
        lookup[key] = {
            "genome_status": str(r.get("status", "") or ""),
            "genome_resolved_accession": str(r.get("resolved_accession", "") or ""),
            "genome_fasta_path": str(r.get("fasta_path", "") or ""),
        }

    def resolve_row(row: pd.Series) -> dict[str, str]:
        candidates: list[str] = []
        candidates += alt_accessions(row.get(col_current, ""))
        if col_asm:
            candidates += alt_accessions(row.get(col_asm, ""))

        for c in candidates:
            hit = lookup.get(c)
            if hit and hit.get("genome_status") == "ok" and hit.get("genome_fasta_path"):
                return hit

        # If no ok hit, still return first non-empty status if present
        for c in candidates:
            hit = lookup.get(c)
            if hit:
                return hit

        return {"genome_status": "missing", "genome_resolved_accession": "", "genome_fasta_path": ""}

    resolved = df_meta.apply(resolve_row, axis=1, result_type="expand")
    df_out = pd.concat([df_meta, resolved], axis=1)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(out_path, sep="\t", index=False)

    missing_df = df_out[df_out["genome_status"] != "ok"]
    print("[Merge] rows:", len(df_out))
    print("[Merge] genome_status counts:")
    vc = df_out["genome_status"].value_counts(dropna=False)
    for k, v in vc.items():
        print(f"  - {k}: {int(v)}")

    if args.out_missing:
        miss_path = Path(args.out_missing)
        miss_path.parent.mkdir(parents=True, exist_ok=True)

        # Prefer current accession, fallback to assembly accession
        acc_series = df_out[col_current].fillna("").astype(str).str.strip()
        if col_asm:
            acc_series = acc_series.where(acc_series != "", df_out[col_asm].fillna("").astype(str).str.strip())

        missing_accs = acc_series[missing_df.index].dropna().astype(str).str.strip()
        missing_accs = missing_accs[missing_accs != ""].drop_duplicates()
        miss_path.write_text("\n".join(missing_accs) + ("\n" if len(missing_accs) else ""), encoding="utf-8")
        print(f"Wrote missing list: {miss_path} (n={len(missing_accs)})")

    print(f"Wrote merged table: {out_path}")


if __name__ == "__main__":
    main()
