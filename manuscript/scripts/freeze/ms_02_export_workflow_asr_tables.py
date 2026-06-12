#!/usr/bin/env python3
"""Export workflow-native ASR tables into manuscript staging outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def load_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    root = repo_root()
    parser = argparse.ArgumentParser(
        description="Stage workflow-native ASR outputs as manuscript supplementary tables and figure-data extracts."
    )
    parser.add_argument(
        "--origin-events",
        type=Path,
        default=root / "outputs/workflow/asr_sensitivity/composition_filtered/origin_events.tsv",
        help="Workflow Fitch origin-events TSV. Defaults to the composition-pruned primary ASR quality frame.",
    )
    parser.add_argument(
        "--sensitivity-summary",
        type=Path,
        default=root / "outputs/workflow/asr_sensitivity/sensitivity_summary.tsv",
        help="Workflow ASR sensitivity summary TSV. Raw workflow scenario labels are reframed for the revised manuscript.",
    )
    parser.add_argument(
        "--supp-origin-out",
        type=Path,
        default=root / "manuscript/supplementary/Supplementary_Table_3_independent_origins.tsv",
        help="Canonical submission-facing output path for the workflow origin-events supplementary table.",
    )
    parser.add_argument(
        "--legacy-supp-origin-out",
        type=Path,
        default=None,
        help="Optional legacy output path for the deprecated workflow_tree-suffixed origin table.",
    )
    parser.add_argument(
        "--supp-sensitivity-out",
        type=Path,
        default=root / "manuscript/supplementary/Supplementary_Table_6_ASR_Sensitivity.tsv",
        help="Canonical submission-facing output path for the workflow ASR sensitivity supplementary table.",
    )
    parser.add_argument(
        "--legacy-supp-sensitivity-out",
        type=Path,
        default=None,
        help="Optional legacy output path for the deprecated workflow_tree-suffixed ASR sensitivity table.",
    )
    parser.add_argument(
        "--figure-origin-out",
        type=Path,
        default=root / "manuscript/figure_data/figure3_workflow_origin_events.tsv",
        help="Output path for the staged workflow origin-events figure-data extract.",
    )
    parser.add_argument(
        "--figure-sensitivity-out",
        type=Path,
        default=root / "manuscript/figure_data/figure3_workflow_asr_sensitivity.tsv",
        help="Output path for the staged workflow ASR sensitivity figure-data extract.",
    )
    parser.add_argument(
        "--resampling-summary",
        type=Path,
        default=root / "outputs/workflow/asr_resampling/resampling_summary.tsv",
        help="Optional workflow ASR resampling summary TSV.",
    )
    parser.add_argument(
        "--figure-resampling-out",
        type=Path,
        default=root / "manuscript/figure_data/figure3_workflow_asr_resampling.tsv",
        help="Output path for the staged workflow ASR resampling summary extract.",
    )
    return parser


def require_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing required {label} file: {path}")


def reframe_sensitivity_rows(
    fieldnames: list[str], rows: list[dict[str, str]]
) -> tuple[list[str], list[dict[str, str]]]:
    """Promote the composition-pruned frame to primary while retaining unpruned comparability rows."""
    if "scenario" not in fieldnames:
        return fieldnames, rows
    by_scenario = {row.get("scenario", ""): row for row in rows}
    expected = {"primary", "support_70", "support_90", "composition_filtered"}
    if not expected.issubset(by_scenario):
        return fieldnames, rows

    renamed = [
        (
            "composition_filtered",
            "composition_pruned_primary_quality_frame",
            "Primary ASR quality frame; prunes the 33 nonreference IQ-TREE composition-failed tips",
        ),
        (
            "primary",
            "unpruned_reference_rooted_comparability",
            "Original Reference/Tohama I rooted unpruned tree retained as comparability sensitivity",
        ),
        (
            "support_70",
            "unpruned_support_ge_70",
            "Fitch restricted to branches with support >= 70 on the unpruned comparability frame",
        ),
        (
            "support_90",
            "unpruned_support_ge_90",
            "Fitch restricted to branches with support >= 90 on the unpruned comparability frame",
        ),
    ]
    reframed_rows: list[dict[str, str]] = []
    for source_name, output_name, notes in renamed:
        row = dict(by_scenario[source_name])
        row["scenario"] = output_name
        row["notes"] = notes
        reframed_rows.append(row)
    return fieldnames, reframed_rows


def main() -> int:
    args = build_arg_parser().parse_args()

    require_file(args.origin_events, "origin-events")
    require_file(args.sensitivity_summary, "sensitivity-summary")

    origin_fields, origin_rows = load_tsv(args.origin_events)
    sensitivity_fields, sensitivity_rows = load_tsv(args.sensitivity_summary)

    if not origin_fields:
        raise ValueError(f"No header found in {args.origin_events}")
    if not sensitivity_fields:
        raise ValueError(f"No header found in {args.sensitivity_summary}")

    sensitivity_fields, sensitivity_rows = reframe_sensitivity_rows(sensitivity_fields, sensitivity_rows)

    origin_outputs = [args.supp_origin_out, args.figure_origin_out]
    if args.legacy_supp_origin_out is not None:
        origin_outputs.append(args.legacy_supp_origin_out)

    for output in origin_outputs:
        write_tsv(output, origin_fields, origin_rows)

    sensitivity_outputs = [args.supp_sensitivity_out, args.figure_sensitivity_out]
    if args.legacy_supp_sensitivity_out is not None:
        sensitivity_outputs.append(args.legacy_supp_sensitivity_out)

    for output in sensitivity_outputs:
        write_tsv(output, sensitivity_fields, sensitivity_rows)

    if args.resampling_summary.exists():
        resampling_fields, resampling_rows = load_tsv(args.resampling_summary)
        if not resampling_fields:
            raise ValueError(f"No header found in {args.resampling_summary}")
        write_tsv(args.figure_resampling_out, resampling_fields, resampling_rows)

    print(f"Exported {len(origin_rows)} workflow origin rows")
    print(f"Exported {len(sensitivity_rows)} workflow sensitivity rows")
    print(f"Supplementary origin table: {args.supp_origin_out}")
    print(f"Supplementary sensitivity table: {args.supp_sensitivity_out}")
    if args.legacy_supp_origin_out is not None:
        print(f"Legacy supplementary origin table: {args.legacy_supp_origin_out}")
    if args.legacy_supp_sensitivity_out is not None:
        print(f"Legacy supplementary sensitivity table: {args.legacy_supp_sensitivity_out}")
    print(f"Figure origin extract: {args.figure_origin_out}")
    print(f"Figure sensitivity extract: {args.figure_sensitivity_out}")
    if args.resampling_summary.exists():
        print(f"Figure resampling extract: {args.figure_resampling_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
