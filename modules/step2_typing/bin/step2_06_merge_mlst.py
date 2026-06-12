#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def parse_mlst(path: Path) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            fasta_path, scheme, st = parts[0].strip(), parts[1].strip(), parts[2].strip()
            profile = ";".join(p.strip() for p in parts[3:] if p.strip())
            rows.append(
                {
                    "genome_fasta_path": fasta_path,
                    "mlst_scheme": scheme,
                    "mlst_st": st,
                    "mlst_profile": profile,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description="Merge mlst results into the QC merged table")
    ap.add_argument("--qc-merged", required=True, help="Merged QC TSV from step2_05 (e.g. outputs/bp_qc_merged.tsv)")
    ap.add_argument("--mlst", required=True, help="MLST output TSV from step2_03")
    ap.add_argument("--out", required=True, help="Output TSV")
    args = ap.parse_args()

    qc_path = Path(args.qc_merged)
    mlst_path = Path(args.mlst)

    if not qc_path.exists() or qc_path.stat().st_size == 0:
        raise SystemExit(f"ERROR: qc merged TSV missing or empty: {qc_path}")
    if not mlst_path.exists() or mlst_path.stat().st_size == 0:
        raise SystemExit(f"ERROR: mlst TSV missing or empty: {mlst_path}")

    df_qc = pd.read_csv(qc_path, sep="\t", dtype=str)
    df_mlst = parse_mlst(mlst_path)
    if df_mlst.empty:
        raise SystemExit(f"ERROR: mlst TSV has no parseable result lines: {mlst_path}")

    df_out = df_qc.merge(df_mlst, on="genome_fasta_path", how="left")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(out_path, sep="\t", index=False)

    n_with = int(df_out["mlst_st"].notna().sum())
    print("[MergeMLST] rows:", len(df_out))
    print("[MergeMLST] rows with mlst_st:", n_with)
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
