#!/usr/bin/env python3
"""Run Fitch parsimony ASR on the current ML tree using the unified manifest."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path


STATE_ORDER = ["intact", "disrupted", "insufficient_data", "uncertain"]
TIP_STATE_COLUMNS = [
    "tree_tip_label",
    "sample_id_canonical",
    "assembly_accession",
    "pastml_id",
    "prn_state",
    "observed_prn_mechanism_call",
    "prn_call_confidence",
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
    "typing_source_tier",
    "is_reference",
    "notes",
]
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
    "inferred_prn_state",
    "parent_inferred_prn_state",
    "transition_from_parent",
    "descendant_tip_count",
    "descendant_state_counts",
    "prn_call_confidence_mode",
    "branch_support",
    "inference_method",
    "notes",
]
TRANSITION_OUTPUT_COLUMNS = [
    "phylo_tree_id",
    "node_id",
    "parent_node_id",
    "node_type",
    "tree_node_label",
    "tip_label",
    "sample_id_canonical",
    "assembly_accession",
    "ancestral_state",
    "descendant_state",
    "transition_from_parent",
    "is_origin_candidate",
    "descendant_tip_count",
    "descendant_disrupted_tip_count",
    "descendant_state_counts",
    "branch_support",
    "inference_method",
    "notes",
]


@dataclass
class Node:
    tree_label: str = ""
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None
    node_id: str = ""
    candidate_states: set[str] = field(default_factory=set)
    inferred_state: str = ""
    descendant_state_counts: Counter[str] = field(default_factory=Counter)
    descendant_tip_count: int = 0

    @property
    def is_tip(self) -> bool:
        return not self.children

    @property
    def branch_support(self) -> str:
        if self.is_tip:
            return ""
        return parse_support_value(self.tree_label)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def normalize_text(value: str) -> str:
    return (value or "").strip()


def parse_support_value(value: str) -> str:
    value = normalize_text(value)
    if not value:
        return ""
    try:
        return f"{float(value):g}"
    except ValueError:
        return ""


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


def load_node_metadata(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    rows = load_tsv_rows(path)
    return {normalize_text(row.get("tree_node_label", "")): row for row in rows}


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


def preferred_state(states: set[str]) -> str:
    for state in STATE_ORDER:
        if state in states:
            return state
    return "uncertain"


def sorted_states(states: set[str]) -> list[str]:
    return sorted(states, key=lambda item: STATE_ORDER.index(item) if item in STATE_ORDER else len(STATE_ORDER))


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
        observed = tip_state_map.get(node.tree_label, "insufficient_data")
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


def serialize_counter(counter: Counter[str]) -> str:
    parts = [f"{key}:{counter[key]}" for key in STATE_ORDER if counter.get(key)]
    return ";".join(parts)


def index_manifest_rows(rows: list[dict[str, str]], key_column: str) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        key = normalize_text(row.get(key_column, ""))
        if key:
            indexed[key] = row
    return indexed


def build_tip_rows(
    nodes: list[Node],
    manifest_rows: list[dict[str, str]],
    key_column: str,
    reference_label: str,
    reference_state: str,
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, str]]:
    manifest_by_key = index_manifest_rows(manifest_rows, key_column)
    tip_rows: list[dict[str, str]] = []
    pastml_rows: list[dict[str, str]] = []
    tip_state_map: dict[str, str] = {}

    for node in nodes:
        if not node.is_tip:
            continue
        tip_label = normalize_text(node.tree_label)
        manifest_row = manifest_by_key.get(tip_label, {})
        is_reference = tip_label == reference_label

        if is_reference:
            observed_mechanism = "reference_assumed_intact"
            observed_state = reference_state
            notes = "reference_tip_assumed_intact"
        else:
            observed_mechanism = normalize_text(manifest_row.get("prn_mechanism_call", ""))
            observed_state = mechanism_to_state(observed_mechanism)
            notes = ""
            if not manifest_row:
                notes = "tip_join_missing"

        tip_state_map[tip_label] = observed_state
        tip_rows.append(
            {
                "tree_tip_label": tip_label,
                "sample_id_canonical": normalize_text(manifest_row.get("sample_id_canonical", "")),
                "assembly_accession": tip_label if not is_reference else "",
                "pastml_id": tip_label,
                "prn_state": observed_state,
                "observed_prn_mechanism_call": observed_mechanism,
                "prn_call_confidence": normalize_text(manifest_row.get("prn_call_confidence", "")),
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
                "typing_source_tier": normalize_text(manifest_row.get("typing_source_tier", "")),
                "is_reference": "True" if is_reference else "False",
                "notes": notes,
            }
        )
        pastml_rows.append({"id": tip_label, "prn_state": observed_state})

    return tip_rows, pastml_rows, tip_state_map


def build_state_rows(
    nodes: list[Node],
    manifest_rows: list[dict[str, str]],
    key_column: str,
    tree_id: str,
    reference_label: str,
    reference_state: str,
    node_metadata: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    manifest_by_key = index_manifest_rows(manifest_rows, key_column)
    output_rows: list[dict[str, str]] = []

    for node in nodes:
        tip_label = normalize_text(node.tree_label) if node.is_tip else ""
        manifest_row = manifest_by_key.get(tip_label, {}) if tip_label and tip_label != reference_label else {}
        parent_state = node.parent.inferred_state if node.parent else ""
        transition = f"{parent_state}->{node.inferred_state}" if parent_state else ""

        if tip_label == reference_label:
            observed_mechanism = "reference_assumed_intact"
            observed_state = reference_state
            notes = "reference_tip_assumed_intact"
        else:
            observed_mechanism = normalize_text(manifest_row.get("prn_mechanism_call", ""))
            observed_state = mechanism_to_state(observed_mechanism) if node.is_tip else ""
            notes = "tip_join_missing" if node.is_tip and not manifest_row else ""

        metadata_row = node_metadata.get(normalize_text(node.tree_label), {})
        branch_support = normalize_text(metadata_row.get("branch_support", "")) or node.branch_support
        original_label = normalize_text(metadata_row.get("original_label", ""))
        if original_label:
            notes = ";".join(filter(None, [notes, f"original_label={original_label}"]))

        output_rows.append(
            {
                "phylo_tree_id": tree_id,
                "node_id": node.node_id,
                "parent_node_id": node.parent.node_id if node.parent else "",
                "node_type": "tip" if node.is_tip else "internal",
                "tree_node_label": normalize_text(node.tree_label),
                "tip_label": tip_label,
                "sample_id_canonical": normalize_text(manifest_row.get("sample_id_canonical", "")),
                "assembly_accession": tip_label if node.is_tip and tip_label != reference_label else "",
                "observed_prn_mechanism_call": observed_mechanism,
                "observed_prn_state": observed_state,
                "candidate_state_set": ";".join(sorted_states(node.candidate_states)),
                "inferred_prn_state": node.inferred_state,
                "parent_inferred_prn_state": parent_state,
                "transition_from_parent": transition,
                "descendant_tip_count": str(node.descendant_tip_count),
                "descendant_state_counts": serialize_counter(node.descendant_state_counts),
                "prn_call_confidence_mode": normalize_text(manifest_row.get("prn_call_confidence", "")) if node.is_tip else "",
                "branch_support": branch_support,
                "inference_method": "fitch_parsimony_on_ml_tree",
                "notes": notes,
            }
        )

    return output_rows


def build_transition_rows(state_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    output_rows: list[dict[str, str]] = []
    for row in state_rows:
        descendant_counts = Counter()
        for token in normalize_text(row.get("descendant_state_counts", "")).split(";"):
            if not token or ":" not in token:
                continue
            state, raw_count = token.split(":", 1)
            try:
                descendant_counts[normalize_text(state)] = int(raw_count)
            except ValueError:
                continue

        ancestral_state = normalize_text(row.get("parent_inferred_prn_state", ""))
        descendant_state = normalize_text(row.get("inferred_prn_state", ""))
        output_rows.append(
            {
                "phylo_tree_id": normalize_text(row.get("phylo_tree_id", "")),
                "node_id": normalize_text(row.get("node_id", "")),
                "parent_node_id": normalize_text(row.get("parent_node_id", "")),
                "node_type": normalize_text(row.get("node_type", "")),
                "tree_node_label": normalize_text(row.get("tree_node_label", "")),
                "tip_label": normalize_text(row.get("tip_label", "")),
                "sample_id_canonical": normalize_text(row.get("sample_id_canonical", "")),
                "assembly_accession": normalize_text(row.get("assembly_accession", "")),
                "ancestral_state": ancestral_state,
                "descendant_state": descendant_state,
                "transition_from_parent": normalize_text(row.get("transition_from_parent", "")),
                "is_origin_candidate": "True" if ancestral_state == "intact" and descendant_state == "disrupted" else "False",
                "descendant_tip_count": normalize_text(row.get("descendant_tip_count", "0")),
                "descendant_disrupted_tip_count": str(descendant_counts.get("disrupted", 0)),
                "descendant_state_counts": normalize_text(row.get("descendant_state_counts", "")),
                "branch_support": normalize_text(row.get("branch_support", "")),
                "inference_method": normalize_text(row.get("inference_method", "")),
                "notes": normalize_text(row.get("notes", "")),
            }
        )
    return output_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run Fitch parsimony ASR on the current ML tree.")
    parser.add_argument(
        "--tree",
        type=Path,
        default=repo_root() / "outputs" / "workflow" / "phylo" / "iqtree2" / "ml_tree.treefile",
        help="Rooted ML tree in Newick format.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=repo_root() / "outputs" / "workflow" / "manifest" / "manifest.tsv",
        help="Unified manifest TSV.",
    )
    parser.add_argument(
        "--tree-id",
        default="workflow_ml_tree",
        help="Identifier written to the output TSVs.",
    )
    parser.add_argument(
        "--tip-key-column",
        default="assembly_accession",
        help="Manifest column matching the tree tip labels.",
    )
    parser.add_argument(
        "--reference-label",
        default="Reference",
        help="Tree label used for the reference tip.",
    )
    parser.add_argument(
        "--reference-state",
        default="intact",
        choices=STATE_ORDER,
        help="State assigned to the reference tip.",
    )
    parser.add_argument(
        "--node-metadata",
        type=Path,
        default=None,
        help="Optional rooted-tree node metadata TSV with branch-support annotations.",
    )
    parser.add_argument(
        "--out-tip-states",
        type=Path,
        required=True,
        help="Tip-state table for auditing and downstream tools.",
    )
    parser.add_argument(
        "--out-pastml-input",
        type=Path,
        required=True,
        help="Minimal tip-state table formatted for later PastML runs.",
    )
    parser.add_argument(
        "--out-states",
        type=Path,
        required=True,
        help="Node-level parsimony state table.",
    )
    parser.add_argument(
        "--out-transitions",
        type=Path,
        required=True,
        help="Parent-child transition table.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    tree_text = args.tree.read_text(encoding="utf-8")
    root = parse_newick(tree_text)
    nodes = assign_node_ids(root)
    manifest_rows = load_tsv_rows(args.manifest)
    node_metadata = load_node_metadata(args.node_metadata)

    tip_rows, pastml_rows, tip_state_map = build_tip_rows(
        nodes=nodes,
        manifest_rows=manifest_rows,
        key_column=args.tip_key_column,
        reference_label=args.reference_label,
        reference_state=args.reference_state,
    )
    fitch_downpass(root, tip_state_map)
    fitch_uppass(root)

    state_rows = build_state_rows(
        nodes=nodes,
        manifest_rows=manifest_rows,
        key_column=args.tip_key_column,
        tree_id=args.tree_id,
        reference_label=args.reference_label,
        reference_state=args.reference_state,
        node_metadata=node_metadata,
    )
    transition_rows = build_transition_rows(state_rows)

    write_tsv(args.out_tip_states, TIP_STATE_COLUMNS, tip_rows)
    write_tsv(args.out_pastml_input, ["id", "prn_state"], pastml_rows)
    write_tsv(args.out_states, STATE_OUTPUT_COLUMNS, state_rows)
    write_tsv(args.out_transitions, TRANSITION_OUTPUT_COLUMNS, transition_rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
