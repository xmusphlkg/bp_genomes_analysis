#!/usr/bin/env python3
"""Reconstruct PRN states on the balanced phylogeny using Fitch parsimony.

This task intentionally produces an auditable ancestral-state table rather than
counting independent origins. Rows include tip-level joins, node-level candidate
state sets, and parent/child state relationships for downstream PHY-04 logic.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.lib.project_paths import project_module_data_root


STATE_ORDER = ["intact", "disrupted", "insufficient_data", "uncertain"]
OUTPUT_COLUMNS = [
    "phylo_tree_id",
    "node_id",
    "parent_node_id",
    "node_type",
    "tip_label",
    "sample_id_canonical",
    "observed_prn_mechanism_call",
    "observed_prn_state",
    "candidate_state_set",
    "inferred_prn_state",
    "parent_inferred_prn_state",
    "transition_from_parent",
    "descendant_tip_count",
    "descendant_state_counts",
    "prn_call_confidence_mode",
    "inference_method",
    "notes",
]


@dataclass
class Node:
    name: str = ""
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None
    node_id: str = ""
    candidate_states: set[str] = field(default_factory=set)
    inferred_state: str = ""
    descendant_state_counts: Counter = field(default_factory=Counter)
    descendant_tip_count: int = 0

    @property
    def is_tip(self) -> bool:
        return not self.children


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


def mechanism_to_state(mechanism: str) -> str:
    mechanism = normalize_text(mechanism)
    if mechanism == "intact":
        return "intact"
    if mechanism.startswith("coding_disrupted_"):
        return "disrupted"
    if mechanism in {"insufficient_data", "uncertain_fragmented_assembly"}:
        return "insufficient_data"
    return "uncertain"


def preferred_state(states: set[str]) -> str:
    for state in STATE_ORDER:
        if state in states:
            return state
    return "uncertain"


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
        if not label:
            return
        if current.children:
            current.name = label
        else:
            current.name = label

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


def assign_node_ids(root: Node) -> list[Node]:
    nodes: list[Node] = []

    def walk(node: Node, counter: list[int]) -> None:
        node.node_id = f"node_{counter[0]:06d}"
        counter[0] += 1
        nodes.append(node)
        for child in node.children:
            walk(child, counter)

    walk(root, [1])
    return nodes


def fitch_downpass(node: Node, tip_state_map: dict[str, str]) -> set[str]:
    if node.is_tip:
        observed = tip_state_map.get(node.name, "uncertain")
        node.candidate_states = {observed}
        node.descendant_state_counts = Counter([observed])
        node.descendant_tip_count = 1
        return node.candidate_states

    child_sets = [fitch_downpass(child, tip_state_map) for child in node.children]
    intersection = set.intersection(*child_sets)
    node.candidate_states = intersection if intersection else set.union(*child_sets)
    node.descendant_state_counts = Counter()
    node.descendant_tip_count = 0
    for child in node.children:
        node.descendant_state_counts.update(child.descendant_state_counts)
        node.descendant_tip_count += child.descendant_tip_count
    return node.candidate_states


def fitch_uppass(node: Node, parent_state: str | None = None) -> None:
    if parent_state and parent_state in node.candidate_states:
        node.inferred_state = parent_state
    else:
        node.inferred_state = preferred_state(node.candidate_states)
    for child in node.children:
        fitch_uppass(child, node.inferred_state)


def serialize_counter(counter: Counter) -> str:
    parts = [f"{key}:{counter[key]}" for key in STATE_ORDER if counter.get(key)]
    return ";".join(parts)


def build_output_rows(nodes: list[Node], mechanism_rows: list[dict[str, str]], tree_id: str) -> list[dict[str, str]]:
    mechanism_by_sample = {
        normalize_text(row.get("sample_id_canonical", "")): row for row in mechanism_rows
    }
    output_rows: list[dict[str, str]] = []
    for node in nodes:
        sample_id = node.name if node.is_tip else ""
        mechanism_row = mechanism_by_sample.get(sample_id, {})
        observed_mechanism = normalize_text(mechanism_row.get("prn_mechanism_call", ""))
        observed_state = mechanism_to_state(observed_mechanism) if observed_mechanism else ""
        parent_state = node.parent.inferred_state if node.parent else ""
        transition = ""
        if parent_state:
            transition = f"{parent_state}->{node.inferred_state}"
        confidence_mode = normalize_text(mechanism_row.get("prn_call_confidence", "")) if node.is_tip else ""
        output_rows.append(
            {
                "phylo_tree_id": tree_id,
                "node_id": node.node_id,
                "parent_node_id": node.parent.node_id if node.parent else "",
                "node_type": "tip" if node.is_tip else "internal",
                "tip_label": node.name if node.is_tip else "",
                "sample_id_canonical": sample_id,
                "observed_prn_mechanism_call": observed_mechanism,
                "observed_prn_state": observed_state,
                "candidate_state_set": ";".join(sorted(node.candidate_states)),
                "inferred_prn_state": node.inferred_state,
                "parent_inferred_prn_state": parent_state,
                "transition_from_parent": transition,
                "descendant_tip_count": str(node.descendant_tip_count),
                "descendant_state_counts": serialize_counter(node.descendant_state_counts),
                "prn_call_confidence_mode": confidence_mode,
                "inference_method": "fitch_parsimony_on_balanced_proxy_tree",
                "notes": "tip_join_missing" if node.is_tip and not mechanism_row else "",
            }
        )
    return output_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Reconstruct PRN states on the balanced global phylogeny using Fitch parsimony."
    )
    parser.add_argument(
        "--tree",
        type=Path,
        default=project_module_data_root("step5_phylogeny_asr") / "outputs" / "bp_global_phylogeny.nwk",
        help="Balanced global phylogeny from PHY-02.",
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
        default=project_module_data_root("step5_phylogeny_asr") / "outputs" / "bp_prn_ancestral_states.tsv",
        help="Ancestral-state output TSV.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    tree_text = args.tree.read_text(encoding="utf-8")
    root = parse_newick(tree_text)
    nodes = assign_node_ids(root)
    mechanism_rows = load_tsv_rows(args.mechanism_calls)
    tip_state_map = {
        normalize_text(row.get("sample_id_canonical", "")): mechanism_to_state(row.get("prn_mechanism_call", ""))
        for row in mechanism_rows
    }
    fitch_downpass(root, tip_state_map)
    fitch_uppass(root)
    output_rows = build_output_rows(nodes, mechanism_rows, "balanced_main_tree")
    write_tsv(args.out, OUTPUT_COLUMNS, output_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
