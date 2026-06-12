#!/usr/bin/env python3
"""Assess whether the current candidate panel is ready for validation."""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.lib.project_paths import project_module_data_root


STEP4_DATA_ROOT = project_module_data_root("step4_prn_validation")
STEP5_DATA_ROOT = project_module_data_root("step5_phylogeny_asr")


def assess_validation_feasibility(
    mechanism_calls_path: str,
    origins_path: str = "",
    min_candidates: int = 20,
) -> dict:
    """Assess whether enough diverse candidates exist for validation."""

    mech = pd.read_csv(mechanism_calls_path, sep="\t", dtype=str)

    disrupted_mechs = {
        "coding_disrupted_is481",
        "coding_disrupted_other_is",
        "coding_disrupted_inversion_or_rearrangement",
        "coding_disrupted_deletion",
        "coding_disrupted_other",
    }
    disrupted = mech[mech["prn_mechanism_call"].isin(disrupted_mechs)].copy()
    intact = mech[mech["prn_mechanism_call"] == "intact"].copy()

    candidates: list[dict[str, object]] = []

    for mech_class in disrupted_mechs:
        subset = disrupted[disrupted["prn_mechanism_call"] == mech_class]
        if len(subset) == 0:
            continue

        has_reads = subset[subset.get("sra_run_accession", pd.Series(dtype=str)).fillna("") != ""]
        pool = has_reads if len(has_reads) >= 3 else subset
        if "country_iso3" in pool.columns:
            selected = pool.drop_duplicates(subset=["country_iso3"], keep="first")
            remaining = pool[~pool.index.isin(selected.index)]
            selected = pd.concat([selected, remaining]).head(5)
        else:
            selected = pool.head(5)

        for _, row in selected.iterrows():
            candidates.append(
                {
                    "sample_id_canonical": row["sample_id_canonical"],
                    "assembly_accession": row.get("assembly_accession", ""),
                    "mechanism": mech_class,
                    "country": row.get("country_iso3", ""),
                    "year": row.get("year", ""),
                    "confidence": row.get("prn_call_confidence", ""),
                    "has_reads": bool(str(row.get("sra_run_accession", "")).strip()),
                    "validation_role": "disrupted_target",
                }
            )

    if len(intact) > 0:
        intact_sample = intact.head(5)
        for _, row in intact_sample.iterrows():
            candidates.append(
                {
                    "sample_id_canonical": row["sample_id_canonical"],
                    "assembly_accession": row.get("assembly_accession", ""),
                    "mechanism": "intact",
                    "country": row.get("country_iso3", ""),
                    "year": row.get("year", ""),
                    "confidence": row.get("prn_call_confidence", ""),
                    "has_reads": bool(str(row.get("sra_run_accession", "")).strip()),
                    "validation_role": "intact_control",
                }
            )

    candidates_df = pd.DataFrame(candidates)
    n_candidates = len(candidates_df)
    n_with_reads = candidates_df["has_reads"].sum() if len(candidates_df) > 0 else 0
    n_mechanisms = candidates_df["mechanism"].nunique() if len(candidates_df) > 0 else 0
    n_countries = candidates_df["country"].nunique() if len(candidates_df) > 0 else 0

    candidate_panel_ready = n_candidates >= min_candidates

    return {
        "assessment": "validation_feasibility",
        "date": pd.Timestamp.now().isoformat(),
        "n_candidates": n_candidates,
        "n_with_reads": int(n_with_reads),
        "n_mechanism_classes": n_mechanisms,
        "n_countries": n_countries,
        "threshold_candidates": min_candidates,
        "candidate_panel_ready": candidate_panel_ready,
        "decision": "READY" if candidate_panel_ready else "NEEDS_DOWNGRADE",
        "candidates": candidates,
        "recommendation": (
            f"{n_candidates} validation candidates identified across {n_mechanisms} mechanism "
            f"classes and {n_countries} countries. Proceed with long-read validation."
            if candidate_panel_ready
            else f"Only {n_candidates}/{min_candidates} candidates found. Consider relaxing selection criteria or prioritizing public long-read data mining."
        ),
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Assess validation feasibility")
    parser.add_argument(
        "--mechanism-calls",
        default=str(STEP4_DATA_ROOT / "outputs" / "bp_prn_mechanism_calls.tsv"),
    )
    parser.add_argument(
        "--origins",
        default=str(STEP5_DATA_ROOT / "outputs" / "bp_prn_independent_origins.tsv"),
    )
    parser.add_argument("--min-candidates", type=int, default=20)
    parser.add_argument(
        "--out",
        default=str(STEP4_DATA_ROOT / "outputs" / "validation_feasibility_report.json"),
    )
    args = parser.parse_args()

    report = assess_validation_feasibility(args.mechanism_calls, args.origins, args.min_candidates)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, default=str)

    print(f"\nValidation feasibility — {report['decision']}")
    print(f"  Candidates: {report['n_candidates']}/{report['threshold_candidates']}")
    print(f"  With reads: {report['n_with_reads']}")
    print(f"  Mechanism classes: {report['n_mechanism_classes']}")
    print(f"  Countries: {report['n_countries']}")
