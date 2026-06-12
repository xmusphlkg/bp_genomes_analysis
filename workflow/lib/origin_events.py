#!/usr/bin/env python3
"""Summarize intact-to-disrupted ASR transitions into origin-event tables."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


OUTPUT_COLUMNS = [
    "origin_id",
    "phylo_tree_id",
    "ancestral_state",
    "descendant_state",
    "clade_id",
    "sister_clade_id",
    "n_tips_total",
    "n_tips_disrupted",
    "n_countries",
    "first_year",
    "last_year",
    "major_lineage",
    "major_lineage_source",
    "major_mlst_st",
    "major_background_label",
    "major_ptxP_label",
    "major_fim3_label",
    "major_fhaB2400_5550_label",
    "major_23s_status",
    "dominant_prn_mechanism",
    "branch_support",
    "origin_support_score",
    "inference_method",
    "notes",
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def normalize_text(value: str) -> str:
    return (value or "").strip()


def load_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"no rows found in {path}")
    return rows


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_year(value: str) -> int | None:
    value = normalize_text(value)
    if not value:
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def parse_float(value: str) -> float | None:
    value = normalize_text(value)
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def majority(counter: Counter[str]) -> str:
    if not counter:
        return ""
    return counter.most_common(1)[0][0]


def build_child_map(rows: list[dict[str, str]]) -> dict[str, list[str]]:
    children: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        parent = normalize_text(row.get("parent_node_id", ""))
        node_id = normalize_text(row.get("node_id", ""))
        if parent:
            children[parent].append(node_id)
    return children


def collect_tip_rows(root_id: str, row_by_id: dict[str, dict[str, str]], children: dict[str, list[str]]) -> list[dict[str, str]]:
    stack = [root_id]
    tips: list[dict[str, str]] = []
    while stack:
        node_id = stack.pop()
        row = row_by_id[node_id]
        if normalize_text(row.get("node_type", "")) == "tip":
            tips.append(row)
            continue
        stack.extend(reversed(children.get(node_id, [])))
    return tips


def build_event_tip_rows(descendant_tip_rows: list[dict[str, str]], manifest_by_assembly: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    event_rows: list[dict[str, str]] = []
    for row in descendant_tip_rows:
        assembly = normalize_text(row.get("assembly_accession", ""))
        manifest_row = manifest_by_assembly.get(assembly, {})
        event_rows.append(
            {
                "tip_label": normalize_text(row.get("tip_label", "")),
                "sample_id_canonical": normalize_text(row.get("sample_id_canonical", "")) or normalize_text(manifest_row.get("sample_id_canonical", "")),
                "assembly_accession": assembly,
                "observed_prn_state": normalize_text(row.get("observed_prn_state", "")),
                "observed_prn_mechanism_call": normalize_text(row.get("observed_prn_mechanism_call", "")) or normalize_text(manifest_row.get("prn_mechanism_call", "")),
                "country_iso3": normalize_text(manifest_row.get("country_iso3", "")),
                "year": normalize_text(manifest_row.get("year", "")),
                "phylo_lineage": normalize_text(manifest_row.get("phylo_lineage", "")),
                "phylo_lineage_source": normalize_text(manifest_row.get("phylo_lineage_source", "")),
                "mlst_st": normalize_text(manifest_row.get("mlst_st", "")),
                "ptxP_label": normalize_text(manifest_row.get("ptxP_label", "")),
                "fim3_label": normalize_text(manifest_row.get("fim3_label", "")),
                "fhaB2400_5550_label": normalize_text(manifest_row.get("fhaB2400_5550_label", "")),
                "marker_23s_status": normalize_text(manifest_row.get("marker_23s_status", "")),
                "background_profile_id": normalize_text(manifest_row.get("background_profile_id", "")),
                "background_display_label": normalize_text(manifest_row.get("background_display_label", "")),
                "prn_call_confidence": normalize_text(manifest_row.get("prn_call_confidence", "")),
            }
        )
    return event_rows


def build_origin_rows(
    transition_rows: list[dict[str, str]],
    state_rows: list[dict[str, str]],
    manifest_rows: list[dict[str, str]],
    min_branch_support: float,
    event_dir: Path,
) -> list[dict[str, str]]:
    row_by_id = {normalize_text(row.get("node_id", "")): row for row in state_rows}
    children = build_child_map(state_rows)
    manifest_by_assembly = {
        normalize_text(row.get("assembly_accession", "")): row
        for row in manifest_rows
        if normalize_text(row.get("assembly_accession", ""))
    }

    output_rows: list[dict[str, str]] = []
    origin_index = 1
    event_dir.mkdir(parents=True, exist_ok=True)

    for row in transition_rows:
        if normalize_text(row.get("is_origin_candidate", "")) != "True":
            continue

        branch_support = parse_float(row.get("branch_support", ""))
        if branch_support is not None and branch_support < min_branch_support:
            continue

        node_id = normalize_text(row.get("node_id", ""))
        descendant_tip_rows = collect_tip_rows(node_id, row_by_id, children)
        event_tip_rows = build_event_tip_rows(descendant_tip_rows, manifest_by_assembly)
        disrupted_rows = [tip_row for tip_row in event_tip_rows if normalize_text(tip_row.get("observed_prn_state", "")) == "disrupted"]

        countries = {
            normalize_text(tip_row.get("country_iso3", ""))
            for tip_row in disrupted_rows
            if normalize_text(tip_row.get("country_iso3", ""))
        }
        years = [parse_year(tip_row.get("year", "")) for tip_row in disrupted_rows]
        years = [year for year in years if year is not None]
        lineage_counter = Counter(
            normalize_text(tip_row.get("phylo_lineage", ""))
            for tip_row in disrupted_rows
            if normalize_text(tip_row.get("phylo_lineage", ""))
        )
        lineage_source_counter = Counter(
            normalize_text(tip_row.get("phylo_lineage_source", ""))
            for tip_row in disrupted_rows
            if normalize_text(tip_row.get("phylo_lineage_source", ""))
        )
        mlst_counter = Counter(
            normalize_text(tip_row.get("mlst_st", ""))
            for tip_row in disrupted_rows
            if normalize_text(tip_row.get("mlst_st", ""))
        )
        background_counter = Counter(
            normalize_text(tip_row.get("background_display_label", ""))
            for tip_row in disrupted_rows
            if normalize_text(tip_row.get("background_display_label", ""))
        )
        background_profile_counter = Counter(
            normalize_text(tip_row.get("background_profile_id", ""))
            for tip_row in disrupted_rows
            if normalize_text(tip_row.get("background_profile_id", ""))
        )
        ptxp_counter = Counter(
            normalize_text(tip_row.get("ptxP_label", ""))
            for tip_row in disrupted_rows
            if normalize_text(tip_row.get("ptxP_label", ""))
        )
        fim3_counter = Counter(
            normalize_text(tip_row.get("fim3_label", ""))
            for tip_row in disrupted_rows
            if normalize_text(tip_row.get("fim3_label", ""))
        )
        fhab_counter = Counter(
            normalize_text(tip_row.get("fhaB2400_5550_label", ""))
            for tip_row in disrupted_rows
            if normalize_text(tip_row.get("fhaB2400_5550_label", ""))
        )
        status_23s_counter = Counter(
            normalize_text(tip_row.get("marker_23s_status", ""))
            for tip_row in disrupted_rows
            if normalize_text(tip_row.get("marker_23s_status", ""))
        )
        mechanism_counter = Counter(
            normalize_text(tip_row.get("observed_prn_mechanism_call", ""))
            for tip_row in disrupted_rows
            if normalize_text(tip_row.get("observed_prn_mechanism_call", ""))
        )

        parent_id = normalize_text(row.get("parent_node_id", ""))
        sister_clade_id = ""
        if parent_id:
            sibling_ids = [child_id for child_id in children.get(parent_id, []) if child_id != node_id]
            if sibling_ids:
                sister_clade_id = sibling_ids[0]

        n_tips_total = len(descendant_tip_rows)
        n_tips_disrupted = len(disrupted_rows)
        support_score = 0.0 if n_tips_total == 0 else n_tips_disrupted / n_tips_total
        origin_id = f"origin_{origin_index:04d}"
        major_lineage = majority(lineage_counter)
        if not major_lineage and background_profile_counter:
            major_lineage = f"profile::{majority(background_profile_counter)}"

        write_tsv(
            event_dir / f"{origin_id}.descendant_tips.tsv",
            [
                "tip_label",
                "sample_id_canonical",
                "assembly_accession",
                "observed_prn_state",
                "observed_prn_mechanism_call",
                "country_iso3",
                "year",
                "phylo_lineage",
                "phylo_lineage_source",
                "mlst_st",
                "ptxP_label",
                "fim3_label",
                "fhaB2400_5550_label",
                "marker_23s_status",
                "background_profile_id",
                "background_display_label",
                "prn_call_confidence",
            ],
            event_tip_rows,
        )

        notes = [
            f"transition={normalize_text(row.get('transition_from_parent', ''))}",
            f"descendant_tip_count={n_tips_total}",
            f"descendant_disrupted_tip_count={n_tips_disrupted}",
            f"tip_package={origin_id}.descendant_tips.tsv",
        ]
        if branch_support is None:
            notes.append("branch_support_not_available")

        output_rows.append(
            {
                "origin_id": origin_id,
                "phylo_tree_id": normalize_text(row.get("phylo_tree_id", "")),
                "ancestral_state": normalize_text(row.get("ancestral_state", "")),
                "descendant_state": normalize_text(row.get("descendant_state", "")),
                "clade_id": node_id,
                "sister_clade_id": sister_clade_id,
                "n_tips_total": str(n_tips_total),
                "n_tips_disrupted": str(n_tips_disrupted),
                "n_countries": str(len(countries)),
                "first_year": "" if not years else str(min(years)),
                "last_year": "" if not years else str(max(years)),
                "major_lineage": major_lineage,
                "major_lineage_source": majority(lineage_source_counter),
                "major_mlst_st": majority(mlst_counter),
                "major_background_label": majority(background_counter),
                "major_ptxP_label": majority(ptxp_counter),
                "major_fim3_label": majority(fim3_counter),
                "major_fhaB2400_5550_label": majority(fhab_counter),
                "major_23s_status": majority(status_23s_counter),
                "dominant_prn_mechanism": majority(mechanism_counter),
                "branch_support": normalize_text(row.get("branch_support", "")),
                "origin_support_score": f"{support_score:.6f}",
                "inference_method": "transition_scan_on_fitch_parsimony_ml_tree",
                "notes": ";".join(notes),
            }
        )
        origin_index += 1

    return output_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Summarize intact-to-disrupted ASR transitions as origin events.")
    parser.add_argument(
        "--states",
        type=Path,
        default=repo_root() / "outputs" / "workflow" / "asr" / "parsimony_states.tsv",
        help="Node-level ASR state table.",
    )
    parser.add_argument(
        "--transitions",
        type=Path,
        default=repo_root() / "outputs" / "workflow" / "asr" / "parsimony_transitions.tsv",
        help="Transition table from asr_parsimony.py.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=repo_root() / "outputs" / "workflow" / "manifest" / "manifest.tsv",
        help="Unified manifest TSV.",
    )
    parser.add_argument(
        "--min-branch-support",
        type=float,
        default=0.0,
        help="Optional branch-support threshold for internal-node origin events.",
    )
    parser.add_argument(
        "--event-dir",
        type=Path,
        required=True,
        help="Directory for per-origin descendant-tip packages.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output TSV for origin events.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    state_rows = load_tsv_rows(args.states)
    transition_rows = load_tsv_rows(args.transitions)
    manifest_rows = load_tsv_rows(args.manifest)
    output_rows = build_origin_rows(
        transition_rows=transition_rows,
        state_rows=state_rows,
        manifest_rows=manifest_rows,
        min_branch_support=args.min_branch_support,
        event_dir=args.event_dir,
    )
    write_tsv(args.out, OUTPUT_COLUMNS, output_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
