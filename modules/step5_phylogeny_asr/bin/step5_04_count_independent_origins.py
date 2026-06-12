#!/usr/bin/env python3
"""Count likely independent origins of PRN disruption from ancestral-state output."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.lib.project_paths import project_module_data_root


OUTPUT_COLUMNS = [
    "origin_id",
    "phylo_tree_id",
    "ancestral_state",
    "descendant_state",
    "clade_id",
    "sister_clade_id",
    "n_tips_disrupted",
    "n_countries",
    "first_year",
    "last_year",
    "major_lineage",
    "major_mlst_st",
    "dominant_prn_mechanism",
    "origin_support_score",
    "inference_method",
    "notes",
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


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


def parse_descendant_counts(value: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for token in normalize_text(value).split(";"):
        if not token or ":" not in token:
            continue
        key, raw = token.split(":", 1)
        try:
            counts[normalize_text(key)] = int(raw)
        except ValueError:
            continue
    return counts


def majority(counter: Counter) -> str:
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


def collect_tip_descendants(root_id: str, row_by_id: dict[str, dict[str, str]], children: dict[str, list[str]]) -> list[str]:
    stack = [root_id]
    tips: list[str] = []
    while stack:
        node_id = stack.pop()
        row = row_by_id[node_id]
        if normalize_text(row.get("node_type", "")) == "tip":
            tip = normalize_text(row.get("sample_id_canonical", "")) or normalize_text(row.get("tip_label", ""))
            if tip:
                tips.append(tip)
            continue
        stack.extend(reversed(children.get(node_id, [])))
    return tips


def build_origin_rows(
    ancestral_rows: list[dict[str, str]],
    mechanism_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    row_by_id = {normalize_text(row.get("node_id", "")): row for row in ancestral_rows}
    children = build_child_map(ancestral_rows)
    mechanism_by_sample = {
        normalize_text(row.get("sample_id_canonical", "")): row for row in mechanism_rows
    }

    output_rows: list[dict[str, str]] = []
    origin_index = 1

    for row in ancestral_rows:
        node_id = normalize_text(row.get("node_id", ""))
        inferred = normalize_text(row.get("inferred_prn_state", ""))
        parent_state = normalize_text(row.get("parent_inferred_prn_state", ""))
        if inferred != "disrupted":
            continue
        if parent_state == "disrupted":
            continue

        descendant_tips = collect_tip_descendants(node_id, row_by_id, children)
        tip_rows = [mechanism_by_sample[tip] for tip in descendant_tips if tip in mechanism_by_sample]
        disrupted_tip_rows = [tip_row for tip_row in tip_rows if normalize_text(tip_row.get("prn_mechanism_call", "")).startswith("coding_disrupted_")]

        countries = {
            normalize_text(tip_row.get("country_iso3", ""))
            for tip_row in disrupted_tip_rows
            if normalize_text(tip_row.get("country_iso3", ""))
        }
        years = [parse_year(tip_row.get("year", "")) for tip_row in disrupted_tip_rows]
        years = [year for year in years if year is not None]
        lineage_counter = Counter(
            normalize_text(tip_row.get("phylo_lineage", ""))
            for tip_row in disrupted_tip_rows
            if normalize_text(tip_row.get("phylo_lineage", ""))
        )
        mlst_counter = Counter(
            normalize_text(tip_row.get("mlst_st", ""))
            for tip_row in disrupted_tip_rows
            if normalize_text(tip_row.get("mlst_st", ""))
        )
        mechanism_counter = Counter(
            normalize_text(tip_row.get("prn_mechanism_call", ""))
            for tip_row in disrupted_tip_rows
            if normalize_text(tip_row.get("prn_mechanism_call", ""))
        )

        parent_id = normalize_text(row.get("parent_node_id", ""))
        sister_clade_id = ""
        if parent_id:
            sibling_ids = [child_id for child_id in children.get(parent_id, []) if child_id != node_id]
            if sibling_ids:
                sister_clade_id = sibling_ids[0]

        descendant_counts = parse_descendant_counts(row.get("descendant_state_counts", ""))
        disrupted_count = descendant_counts.get("disrupted", len(disrupted_tip_rows))
        total_descendant = int(normalize_text(row.get("descendant_tip_count", "0")) or "0")
        support_score = 0.0
        if total_descendant > 0:
            support_score = disrupted_count / total_descendant

        notes = [
            f"transition={normalize_text(row.get('transition_from_parent', '')) or 'root'}",
            f"descendant_tip_count={total_descendant}",
            f"mapped_tip_rows={len(tip_rows)}",
        ]
        if not disrupted_tip_rows:
            notes.append("no_disrupted_tip_metadata_recovered")

        output_rows.append(
            {
                "origin_id": f"origin_{origin_index:04d}",
                "phylo_tree_id": normalize_text(row.get("phylo_tree_id", "")),
                "ancestral_state": parent_state or "root",
                "descendant_state": inferred,
                "clade_id": node_id,
                "sister_clade_id": sister_clade_id,
                "n_tips_disrupted": str(disrupted_count),
                "n_countries": str(len(countries)),
                "first_year": "" if not years else str(min(years)),
                "last_year": "" if not years else str(max(years)),
                "major_lineage": majority(lineage_counter),
                "major_mlst_st": majority(mlst_counter),
                "dominant_prn_mechanism": majority(mechanism_counter),
                "origin_support_score": f"{support_score:.6f}",
                "inference_method": "transition_scan_on_fitch_parsimony_balanced_proxy_tree",
                "notes": ";".join(notes),
            }
        )
        origin_index += 1

    return output_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Count likely independent PRN disruption origins from ancestral-state transitions."
    )
    parser.add_argument(
        "--ancestral-states",
        type=Path,
        default=project_module_data_root("step5_phylogeny_asr")
        / "outputs"
        / "bp_prn_ancestral_states.tsv",
        help="Ancestral-state table from PHY-03.",
    )
    parser.add_argument(
        "--mechanism-calls",
        type=Path,
        default=project_module_data_root("step4_prn_validation") / "outputs" / "bp_prn_mechanism_calls.tsv",
        help="PRN mechanism call table.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=project_module_data_root("step5_phylogeny_asr")
        / "outputs"
        / "bp_prn_independent_origins.tsv",
        help="Independent-origin event table.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    ancestral_rows = load_tsv_rows(args.ancestral_states)
    mechanism_rows = load_tsv_rows(args.mechanism_calls)
    output_rows = build_origin_rows(ancestral_rows, mechanism_rows)
    write_tsv(args.out, OUTPUT_COLUMNS, output_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
