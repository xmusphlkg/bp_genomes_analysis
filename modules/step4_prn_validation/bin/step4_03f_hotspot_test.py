#!/usr/bin/env python3
"""Run a permutation-based clustering test for prn-local IS insertion positions."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np

from step4_02_scan_prn_mechanisms import (
    load_tsv_rows,
    normalize_text,
    parse_int,
    project_module_data_root,
    write_tsv,
)


STEP4_DATA_ROOT = project_module_data_root("step4_prn_validation")
STEP2_DATA_ROOT = project_module_data_root("step2_typing")


OUTPUT_COLUMNS = [
    "metric",
    "window_bp",
    "n_events",
    "prn_start",
    "prn_end",
    "observed_value",
    "expected_mean",
    "empirical_p_value",
    "top_window_start",
    "top_window_end",
    "top_window_event_count",
]


def load_prn_locus(gbff_path: Path) -> tuple[int, int]:
    current_feature: dict[str, str] | None = None

    def finalize_feature(feature: dict[str, str] | None) -> tuple[int, int] | None:
        if feature is None:
            return None
        qualifiers = "\n".join(feature["qualifiers"])
        if '/gene="prn"' not in qualifiers and "pertactin autotransporter" not in qualifiers.casefold():
            return None
        coordinates = [int(token) for token in re.findall(r"\d+", feature["location"])]
        if not coordinates:
            return None
        return min(coordinates), max(coordinates)

    with gbff_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if re.match(r"^     \S", line) and not line.startswith("                     /"):
                parsed = finalize_feature(current_feature)
                if parsed is not None:
                    return parsed
                current_feature = {
                    "type": line[5:21].strip(),
                    "location": line[21:].strip(),
                    "qualifiers": [],
                }
                continue
            if current_feature is not None and line.startswith("                     /"):
                current_feature["qualifiers"].append(line.strip())

    parsed = finalize_feature(current_feature)
    if parsed is not None:
        return parsed
    raise ValueError(f"prn locus not found in {gbff_path}")


def choose_sample_positions(evidence_rows: list[dict[str, str]]) -> list[int]:
    hits_by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in evidence_rows:
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        if sample_id:
            hits_by_sample[sample_id].append(row)

    positions: list[int] = []
    for sample_id in sorted(hits_by_sample):
        rows = hits_by_sample[sample_id]
        rows.sort(
            key=lambda row: (
                -(parse_int(row.get("prn_overlap_bp", "")) or 0),
                0 if normalize_text(row.get("tool", "")) == "panisa" else 1,
                -(parse_int(row.get("total_clipped_reads", "")) or 0),
                parse_int(row.get("locus_start", "")) or 0,
            )
        )
        chosen = rows[0]
        start = parse_int(chosen.get("locus_start", "")) or 0
        end = parse_int(chosen.get("locus_end", "")) or start
        positions.append((start + end) // 2)
    return positions


def max_window_count(positions: list[int], window_bp: int) -> tuple[int, int, int]:
    if not positions:
        return 0, 0, 0
    positions = sorted(positions)
    best_count = 0
    best_start = positions[0]
    left = 0
    for right, position in enumerate(positions):
        while position - positions[left] > window_bp:
            left += 1
        count = right - left + 1
        if count > best_count:
            best_count = count
            best_start = positions[left]
    return best_count, best_start, best_start + window_bp


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run a permutation-based clustering test on prn-local IS insertion evidence."
    )
    parser.add_argument(
        "--evidence",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_read_validation_is_calls.tsv",
        help="Per-hit evidence table from step4_03_validate_prn_with_reads.py.",
    )
    parser.add_argument(
        "--reference-gbff",
        type=Path,
        default=STEP2_DATA_ROOT
        / "outputs/_ref/GCF_000195715.1/ncbi_dataset/data/GCF_000195715.1/genomic.gbff",
        help="Reference GBFF used to locate the prn locus.",
    )
    parser.add_argument(
        "--window-bp",
        type=int,
        default=50,
        help="Sliding window width used to define clustering.",
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=10000,
        help="Number of uniform permutations across the prn locus.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260405,
        help="Random seed for reproducible permutation p-values.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_is_hotspot_results.tsv",
        help="Permutation summary TSV.",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_is_hotspot_density.pdf",
        help="Histogram-style PDF of prn-local insertion positions.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    evidence_rows = [
        row
        for row in load_tsv_rows(args.evidence)
        if normalize_text(row.get("tool", "")) in {"ismapper", "panisa"}
    ]
    prn_start, prn_end = load_prn_locus(args.reference_gbff)
    positions = choose_sample_positions(evidence_rows)
    observed_count, top_window_start, top_window_end = max_window_count(positions, args.window_bp)

    rng = np.random.default_rng(args.seed)
    simulated_counts: list[int] = []
    for _ in range(args.n_permutations):
        simulated = sorted(rng.integers(prn_start, prn_end + 1, size=len(positions)).tolist())
        simulated_counts.append(max_window_count(simulated, args.window_bp)[0])

    expected_mean = float(np.mean(simulated_counts)) if simulated_counts else 0.0
    empirical_p = (
        (1 + sum(count >= observed_count for count in simulated_counts)) / (1 + len(simulated_counts))
        if simulated_counts
        else 1.0
    )
    write_tsv(
        args.out,
        OUTPUT_COLUMNS,
        [
            {
                "metric": "max_events_in_window",
                "window_bp": str(args.window_bp),
                "n_events": str(len(positions)),
                "prn_start": str(prn_start),
                "prn_end": str(prn_end),
                "observed_value": str(observed_count),
                "expected_mean": f"{expected_mean:.6f}",
                "empirical_p_value": f"{empirical_p:.6f}",
                "top_window_start": str(top_window_start),
                "top_window_end": str(top_window_end),
                "top_window_event_count": str(observed_count),
            }
        ],
    )

    args.plot.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(8, 4.5))
    if positions:
        n_bins = max(5, min(20, max(1, (prn_end - prn_start + 1) // max(args.window_bp, 1))))
        axis.hist(positions, bins=n_bins, color="#2d6a4f", alpha=0.85)
        axis.vlines(positions, ymin=0, ymax=0.3, color="#1b4332", linewidth=1)
    axis.set_xlim(prn_start, prn_end)
    axis.set_xlabel("prn coordinate (Tohama I)")
    axis.set_ylabel("Validated insertion count")
    axis.set_title("prn-local IS insertion density")
    figure.tight_layout()
    figure.savefig(args.plot)
    plt.close(figure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
