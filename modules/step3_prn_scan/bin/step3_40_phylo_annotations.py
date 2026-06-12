#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from workflow.lib.project_paths import project_module_data_root


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def norm(s: pd.Series) -> pd.Series:
    return s.astype(str).where(~s.isna(), "NA").replace({"nan": "NA", "None": "NA"})


def pick_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Step3E: build phylogeny annotation table by merging phylo manifest with prn disruption calls"
    )
    ap.add_argument("--manifest", required=True, help="step3_prn_scan/outputs/bp_phylo_manifest.tsv")
    ap.add_argument("--merged", required=True, help="step3_prn_scan/outputs/bp_qc_merged_mlst_markers_prn.tsv")
    step2_root = project_module_data_root("step2_typing")
    ap.add_argument(
        "--typing-manifest",
        default=str(step2_root / "outputs" / "bp_genotype_manifest.tsv"),
        help="Optional standardized Step2 genotype manifest TSV.",
    )
    ap.add_argument("--out", required=True, help="Output TSV for tree annotation")
    args = ap.parse_args()

    manifest = pd.read_csv(Path(args.manifest), sep="\t", dtype=str)
    merged = pd.read_csv(Path(args.merged), sep="\t", dtype=str)

    # Column mapping
    man_acc_col = pick_col(manifest, ["Current Accession", "Assembly Accession", "genome_resolved_accession", "sample_id"])
    if man_acc_col is None:
        raise SystemExit("ERROR: manifest missing an accession-like column")

    mer_acc_col = pick_col(merged, ["Current Accession", "Assembly Accession", "genome_resolved_accession"])
    if mer_acc_col is None:
        raise SystemExit("ERROR: merged table missing an accession-like column")

    for df, col in [(manifest, man_acc_col), (merged, mer_acc_col)]:
        df[col] = norm(df[col])

    typing_path = Path(args.typing_manifest)
    if typing_path.exists() and typing_path.stat().st_size > 0:
        typing = pd.read_csv(typing_path, sep="\t", dtype=str)
        typing_cols = [
            "assembly_accession",
            "ptxP_label",
            "marker_ptxP_promoter_hash",
            "fim3_label",
            "fhaB2400_5550_label",
            "marker_23s_status",
        ]
        available_typing_cols = [col for col in typing_cols if col in typing.columns]
        if "assembly_accession" in available_typing_cols:
            typing["assembly_accession"] = norm(typing["assembly_accession"])
            merged = merged.merge(
                typing[available_typing_cols].drop_duplicates(subset=["assembly_accession"]),
                left_on=mer_acc_col,
                right_on="assembly_accession",
                how="left",
            )

    # De-dup merged side to one row per accession
    merged_keyed = merged.drop_duplicates(subset=[mer_acc_col], keep="first").copy()

    ann = manifest.merge(merged_keyed, left_on=man_acc_col, right_on=mer_acc_col, how="left", suffixes=("", "_merged"))

    # Canonical fields to export
    year_col = pick_col(ann, ["year", "year_x", "year_y"])
    country_col = pick_col(ann, ["country", "country_x", "country_y"])
    st_col = pick_col(ann, ["mlst_st", "mlst_st_x", "mlst_st_y"])
    ptxp_label_col = pick_col(ann, ["ptxP_label"])
    ptxp_hash_col = pick_col(ann, ["marker_ptxP_promoter_hash", "marker_ptxP_promoter", "marker_ptxP_promoter_allele"])

    out = pd.DataFrame(
        {
            "sample_id": norm(ann["sample_id"]) if "sample_id" in ann.columns else norm(ann[man_acc_col]),
            "current_accession": norm(ann[man_acc_col]),
            "country": norm(ann[country_col]) if country_col else "NA",
            "year": norm(ann[year_col]) if year_col else "NA",
            "mlst_st": norm(ann[st_col]) if st_col else "NA",
            "prn_call": norm(ann["prn_call"]) if "prn_call" in ann.columns else "NA",
            "prn_qcov_union_pct": norm(ann["prn_qcov_union_pct"]) if "prn_qcov_union_pct" in ann.columns else "NA",
            "prn_best_single_qcov_pct": norm(ann["prn_best_single_qcov_pct"]) if "prn_best_single_qcov_pct" in ann.columns else "NA",
            "23s_A2047G_call": norm(ann["23s_A2047G_call"]) if "23s_A2047G_call" in ann.columns else "NA",
            "ptxP_label": norm(ann[ptxp_label_col]) if ptxp_label_col else "NA",
            "ptxP_hash": norm(ann[ptxp_hash_col]) if ptxp_hash_col else "NA",
        }
    )

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(Path(args.out), sep="\t", index=False)
    print(f"Wrote: {args.out} ({len(out)} rows)")


if __name__ == "__main__":
    main()
