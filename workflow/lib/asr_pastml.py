#!/usr/bin/env python3
"""Normalize PastML output tables and summarize origin-event support."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path


STATE_OUTPUT_COLUMNS = [
    "phylo_tree_id",
    "node_id",
    "parent_node_id",
    "node_type",
    "tree_node_label",
    "tip_label",
    "sample_id_canonical",
    "assembly_accession",
    "observed_prn_mechanism_call",
    "observed_prn_state",
    "candidate_state_set",
    "prediction_class",
    "branch_support",
    "inference_method",
    "notes",
]
ORIGIN_OUTPUT_COLUMNS = [
    "origin_id",
    "phylo_tree_id",
    "origin_confidence",
    "ancestral_candidate_state_set",
    "descendant_candidate_state_set",
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
    "inference_method",
    "notes",
]
SUMMARY_COLUMNS = ["metric", "value", "notes"]
STATE_ORDER = ["intact", "disrupted", "insufficient_data", "uncertain"]


@dataclass
class Node:
    tree_label: str = ""
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None

    @property
    def is_tip(self) -> bool:
        return not self.children


def normalize_text(value: str) -> str:
    return (value or "").strip()


def load_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"no rows found in {path}")
    return rows


def count_tsv_rows(path: Path) -> int:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return sum(1 for _ in reader)


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


def mechanism_to_state(mechanism: str) -> str:
    mechanism = normalize_text(mechanism)
    if mechanism in {"intact", "reference_assumed_intact"}:
        return "intact"
    if mechanism.startswith("coding_disrupted_"):
        return "disrupted"
    if mechanism in {"insufficient_data", "uncertain_fragmented_assembly", "not_available_current_step3"}:
        return "insufficient_data"
    if not mechanism:
        return "insufficient_data"
    return "uncertain"


def sorted_states(states: set[str]) -> list[str]:
    return sorted(states, key=lambda item: STATE_ORDER.index(item) if item in STATE_ORDER else len(STATE_ORDER))


def majority(counter: Counter[str]) -> str:
    if not counter:
        return ""
    return counter.most_common(1)[0][0]


def parse_newick(text: str) -> Node:
    root = Node()
    current = root
    stack: list[Node] = []
    token = ""
    reading_length = False

    def flush_token() -> None:
        nonlocal token, current
        label = normalize_text(token)
        token = ""
        if label:
            current.tree_label = label

    for char in text.strip():
        if reading_length:
            if char in ",);":
                reading_length = False
            else:
                continue
        if char == "(":
            child = Node(parent=current)
            current.children.append(child)
            stack.append(current)
            current = child
        elif char == ",":
            flush_token()
            sibling = Node(parent=stack[-1])
            stack[-1].children.append(sibling)
            current = sibling
        elif char == ")":
            flush_token()
            current = stack.pop()
        elif char == ":":
            flush_token()
            reading_length = True
        elif char == ";":
            flush_token()
        else:
            token += char

    if len(root.children) == 1:
        return root.children[0]
    return root


def load_node_metadata(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None or not path.exists():
        return {}
    rows = load_tsv_rows(path)
    return {normalize_text(row.get("tree_node_label", "")): row for row in rows}


def group_states(raw_rows: list[dict[str, str]], column: str) -> dict[str, set[str]]:
    grouped: dict[str, set[str]] = defaultdict(set)
    for row in raw_rows:
        node_name = normalize_text(row.get("node", ""))
        value = normalize_text(row.get(column, ""))
        if node_name and value:
            grouped[node_name].add(value)
    return grouped


def collect_nodes(root: Node) -> list[Node]:
    nodes: list[Node] = []

    def walk(node: Node) -> None:
        nodes.append(node)
        for child in node.children:
            walk(child)

    walk(root)
    return nodes


def collect_tip_descendants(node: Node) -> list[Node]:
    tips: list[Node] = []
    stack = [node]
    while stack:
        current = stack.pop()
        if current.is_tip:
            tips.append(current)
            continue
        stack.extend(reversed(current.children))
    return tips


def build_state_rows(
    nodes: list[Node],
    grouped_states: dict[str, set[str]],
    manifest_rows: list[dict[str, str]],
    tree_id: str,
    reference_label: str,
    node_metadata: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    manifest_by_assembly = {
        normalize_text(row.get("assembly_accession", "")): row
        for row in manifest_rows
        if normalize_text(row.get("assembly_accession", ""))
    }
    output_rows: list[dict[str, str]] = []

    for node in nodes:
        tree_node_label = normalize_text(node.tree_label)
        state_set = grouped_states.get(tree_node_label, set())
        manifest_row = manifest_by_assembly.get(tree_node_label, {}) if node.is_tip and tree_node_label != reference_label else {}
        observed_mechanism = "reference_assumed_intact" if node.is_tip and tree_node_label == reference_label else normalize_text(manifest_row.get("prn_mechanism_call", ""))
        notes = []
        if not state_set:
            notes.append("pastml_state_missing")
        if node.is_tip and tree_node_label == reference_label:
            observed_state = "intact"
        elif node.is_tip:
            observed_state = mechanism_to_state(observed_mechanism)
            if not manifest_row:
                notes.append("tip_join_missing")
        else:
            observed_state = ""

        metadata_row = node_metadata.get(tree_node_label, {})
        if normalize_text(metadata_row.get("original_label", "")):
            notes.append(f"original_label={normalize_text(metadata_row.get('original_label', ''))}")

        output_rows.append(
            {
                "phylo_tree_id": tree_id,
                "node_id": tree_node_label,
                "parent_node_id": normalize_text(node.parent.tree_label) if node.parent else "",
                "node_type": "tip" if node.is_tip else "internal",
                "tree_node_label": tree_node_label,
                "tip_label": tree_node_label if node.is_tip else "",
                "sample_id_canonical": normalize_text(manifest_row.get("sample_id_canonical", "")),
                "assembly_accession": tree_node_label if node.is_tip and tree_node_label != reference_label else "",
                "observed_prn_mechanism_call": observed_mechanism,
                "observed_prn_state": observed_state,
                "candidate_state_set": ";".join(sorted_states(state_set)),
                "prediction_class": sorted_states(state_set)[0] if len(state_set) == 1 else ("ambiguous" if state_set else "missing"),
                "branch_support": normalize_text(metadata_row.get("branch_support", "")),
                "inference_method": "pastml_mppa_f81_on_reference_rooted_ml_tree",
                "notes": ";".join(notes),
            }
        )
    return output_rows


def build_origin_rows(
    root: Node,
    grouped_states: dict[str, set[str]],
    manifest_rows: list[dict[str, str]],
    tree_id: str,
    reference_label: str,
    node_metadata: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    manifest_by_assembly = {
        normalize_text(row.get("assembly_accession", "")): row
        for row in manifest_rows
        if normalize_text(row.get("assembly_accession", ""))
    }

    output_rows: list[dict[str, str]] = []
    origin_index = 1

    for node in collect_nodes(root):
        if node.parent is None:
            continue
        parent_set = grouped_states.get(normalize_text(node.parent.tree_label), set())
        child_set = grouped_states.get(normalize_text(node.tree_label), set())
        if "intact" not in parent_set or "disrupted" not in child_set or "disrupted" in parent_set:
            continue

        origin_confidence = "compatible"
        if "intact" not in child_set:
            origin_confidence = "strict"

        descendant_tips = collect_tip_descendants(node)
        disrupted_rows = []
        for tip in descendant_tips:
            tip_label = normalize_text(tip.tree_label)
            if tip_label == reference_label:
                continue
            manifest_row = manifest_by_assembly.get(tip_label, {})
            mechanism = normalize_text(manifest_row.get("prn_mechanism_call", ""))
            if mechanism_to_state(mechanism) == "disrupted":
                disrupted_rows.append(manifest_row)

        countries = {
            normalize_text(row.get("country_iso3", ""))
            for row in disrupted_rows
            if normalize_text(row.get("country_iso3", ""))
        }
        years = [parse_year(row.get("year", "")) for row in disrupted_rows]
        years = [year for year in years if year is not None]
        lineage_counter = Counter(
            normalize_text(row.get("phylo_lineage", ""))
            for row in disrupted_rows
            if normalize_text(row.get("phylo_lineage", ""))
        )
        lineage_source_counter = Counter(
            normalize_text(row.get("phylo_lineage_source", ""))
            for row in disrupted_rows
            if normalize_text(row.get("phylo_lineage_source", ""))
        )
        mlst_counter = Counter(
            normalize_text(row.get("mlst_st", ""))
            for row in disrupted_rows
            if normalize_text(row.get("mlst_st", ""))
        )
        background_counter = Counter(
            normalize_text(row.get("background_display_label", ""))
            for row in disrupted_rows
            if normalize_text(row.get("background_display_label", ""))
        )
        background_profile_counter = Counter(
            normalize_text(row.get("background_profile_id", ""))
            for row in disrupted_rows
            if normalize_text(row.get("background_profile_id", ""))
        )
        ptxp_counter = Counter(
            normalize_text(row.get("ptxP_label", ""))
            for row in disrupted_rows
            if normalize_text(row.get("ptxP_label", ""))
        )
        fim3_counter = Counter(
            normalize_text(row.get("fim3_label", ""))
            for row in disrupted_rows
            if normalize_text(row.get("fim3_label", ""))
        )
        fhab_counter = Counter(
            normalize_text(row.get("fhaB2400_5550_label", ""))
            for row in disrupted_rows
            if normalize_text(row.get("fhaB2400_5550_label", ""))
        )
        status_23s_counter = Counter(
            normalize_text(row.get("marker_23s_status", ""))
            for row in disrupted_rows
            if normalize_text(row.get("marker_23s_status", ""))
        )
        mechanism_counter = Counter(
            normalize_text(row.get("prn_mechanism_call", ""))
            for row in disrupted_rows
            if normalize_text(row.get("prn_mechanism_call", ""))
        )

        sibling_ids = [normalize_text(child.tree_label) for child in node.parent.children if child is not node]
        metadata_row = node_metadata.get(normalize_text(node.tree_label), {})
        major_lineage = majority(lineage_counter)
        if not major_lineage and background_profile_counter:
            major_lineage = f"profile::{majority(background_profile_counter)}"

        output_rows.append(
            {
                "origin_id": f"pastml_origin_{origin_index:04d}",
                "phylo_tree_id": tree_id,
                "origin_confidence": origin_confidence,
                "ancestral_candidate_state_set": ";".join(sorted_states(parent_set)),
                "descendant_candidate_state_set": ";".join(sorted_states(child_set)),
                "clade_id": normalize_text(node.tree_label),
                "sister_clade_id": sibling_ids[0] if sibling_ids else "",
                "n_tips_total": str(len(descendant_tips)),
                "n_tips_disrupted": str(len(disrupted_rows)),
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
                "branch_support": normalize_text(metadata_row.get("branch_support", "")),
                "inference_method": "transition_scan_on_pastml_mppa_f81_ml_tree",
                "notes": ";".join(
                    [
                        f"parent_state_set={';'.join(sorted_states(parent_set))}",
                        f"descendant_state_set={';'.join(sorted_states(child_set))}",
                        f"descendant_tip_count={len(descendant_tips)}",
                        f"descendant_disrupted_tip_count={len(disrupted_rows)}",
                    ]
                ),
            }
        )
        origin_index += 1

    return output_rows


def build_summary_rows(
    fitch_events_path: Path | None,
    pastml_state_rows: list[dict[str, str]],
    pastml_origin_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    fitch_event_count = ""
    if fitch_events_path is not None and fitch_events_path.exists():
        fitch_event_count = str(count_tsv_rows(fitch_events_path))

    strict_count = sum(1 for row in pastml_origin_rows if normalize_text(row.get("origin_confidence", "")) == "strict")
    compatible_count = sum(1 for row in pastml_origin_rows if normalize_text(row.get("origin_confidence", "")) == "compatible")
    ambiguous_count = sum(1 for row in pastml_state_rows if normalize_text(row.get("prediction_class", "")) == "ambiguous")

    return [
        {"metric": "fitch_origin_events", "value": fitch_event_count, "notes": ""},
            {"metric": "pastml_origin_events_strict", "value": str(strict_count), "notes": "parent excludes disrupted and child excludes intact"},
            {"metric": "pastml_origin_events_compatible", "value": str(compatible_count), "notes": "parent excludes disrupted and child contains both intact-compatible and disrupted-compatible support"},
        {"metric": "pastml_origin_events_total", "value": str(len(pastml_origin_rows)), "notes": "strict + compatible"},
        {"metric": "pastml_ambiguous_nodes", "value": str(ambiguous_count), "notes": "nodes with multi-state MPPA predictions"},
    ]


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Normalize PastML outputs and summarize transition support.")
    parser.add_argument("--tree", type=Path, required=True, help="Rooted tree used for PastML.")
    parser.add_argument("--raw-states", type=Path, required=True, help="PastML combined ancestral states table.")
    parser.add_argument("--manifest", type=Path, required=True, help="Unified manifest TSV.")
    parser.add_argument("--column", default="prn_state", help="Character column reconstructed by PastML.")
    parser.add_argument("--tree-id", default="workflow_ml_tree", help="Tree identifier written to output tables.")
    parser.add_argument("--reference-label", default="Reference", help="Reference tip label.")
    parser.add_argument("--node-metadata", type=Path, default=None, help="Optional rooted-tree node metadata TSV.")
    parser.add_argument("--fitch-events", type=Path, default=None, help="Optional Fitch origin-events TSV for comparison summary.")
    parser.add_argument("--out-states", type=Path, required=True, help="Normalized PastML state TSV.")
    parser.add_argument("--out-origin-events", type=Path, required=True, help="PastML origin-event TSV.")
    parser.add_argument("--out-summary", type=Path, required=True, help="Comparison summary TSV.")
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    root = parse_newick(args.tree.read_text(encoding="utf-8"))
    manifest_rows = load_tsv_rows(args.manifest)
    raw_rows = load_tsv_rows(args.raw_states)
    grouped = group_states(raw_rows, args.column)
    node_metadata = load_node_metadata(args.node_metadata)

    state_rows = build_state_rows(
        nodes=collect_nodes(root),
        grouped_states=grouped,
        manifest_rows=manifest_rows,
        tree_id=args.tree_id,
        reference_label=args.reference_label,
        node_metadata=node_metadata,
    )
    origin_rows = build_origin_rows(
        root=root,
        grouped_states=grouped,
        manifest_rows=manifest_rows,
        tree_id=args.tree_id,
        reference_label=args.reference_label,
        node_metadata=node_metadata,
    )
    summary_rows = build_summary_rows(args.fitch_events, state_rows, origin_rows)

    write_tsv(args.out_states, STATE_OUTPUT_COLUMNS, state_rows)
    write_tsv(args.out_origin_events, ORIGIN_OUTPUT_COLUMNS, origin_rows)
    write_tsv(args.out_summary, SUMMARY_COLUMNS, summary_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
