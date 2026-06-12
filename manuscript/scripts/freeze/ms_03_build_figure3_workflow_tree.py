#!/usr/bin/env python3
"""Build workflow-native Figure 3 tree layout tables from the rooted 191-tip ML tree."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
from pathlib import Path


TREE_SEGMENT_COLUMNS = [
    "segment_id",
    "segment_type",
    "parent_node_id",
    "child_node_id",
    "x_cladogram",
    "xend_cladogram",
    "y_order",
    "yend_order",
]

TREE_NODE_COLUMNS = [
    "node_id",
    "parent_node_id",
    "node_type",
    "tip_label",
    "x_cladogram",
    "x_branch_length",
    "y_order",
    "branch_support",
    "is_reference",
    "is_fitch_origin",
    "origin_id",
    "n_tips_total",
    "n_tips_disrupted",
    "n_countries",
    "origin_first_year",
    "origin_last_year",
    "dominant_prn_mechanism",
    "origin_support_score",
    "observed_prn_state",
    "prn_mechanism_call",
    "country_iso3",
    "year",
    "mlst_st",
    "phylo_lineage",
]


@dataclass
class TreeNode:
    name: str = ""
    length: float = 0.0
    children: list["TreeNode"] = field(default_factory=list)
    parent: "TreeNode | None" = None
    x_cladogram: float = 0.0
    x_branch_length: float = 0.0
    y_order: float = 0.0

    @property
    def is_tip(self) -> bool:
        return not self.children


class NewickParser:
    def __init__(self, text: str) -> None:
        self.text = text
        self.index = 0

    def parse(self) -> TreeNode:
        node = self._parse_subtree()
        self._skip_whitespace()
        if self._peek() == ";":
            self.index += 1
        self._skip_whitespace()
        if self.index != len(self.text):
            raise ValueError(f"unexpected trailing content at position {self.index}")
        return node

    def _parse_subtree(self) -> TreeNode:
        self._skip_whitespace()
        if self._peek() == "(":
            self.index += 1
            children = [self._parse_subtree()]
            while True:
                self._skip_whitespace()
                token = self._peek()
                if token == ",":
                    self.index += 1
                    children.append(self._parse_subtree())
                    continue
                if token == ")":
                    self.index += 1
                    break
                raise ValueError(f"unexpected token {token!r} at position {self.index}")
            name = self._parse_label()
            length = self._parse_length()
            return TreeNode(name=name, length=length, children=children)

        name = self._parse_label()
        length = self._parse_length()
        return TreeNode(name=name, length=length)

    def _parse_label(self) -> str:
        self._skip_whitespace()
        start_index = self.index
        while self.index < len(self.text) and self.text[self.index] not in ":,();":
            self.index += 1
        return self.text[start_index:self.index].strip()

    def _parse_length(self) -> float:
        self._skip_whitespace()
        if self._peek() != ":":
            return 0.0
        self.index += 1
        start_index = self.index
        while self.index < len(self.text) and self.text[self.index] not in ",();":
            self.index += 1
        text = self.text[start_index:self.index].strip()
        return float(text) if text else 0.0

    def _skip_whitespace(self) -> None:
        while self.index < len(self.text) and self.text[self.index].isspace():
            self.index += 1

    def _peek(self) -> str:
        if self.index >= len(self.text):
            return ""
        return self.text[self.index]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def iter_nodes(node: TreeNode) -> list[TreeNode]:
    nodes = [node]
    for child in node.children:
        nodes.extend(iter_nodes(child))
    return nodes


def assign_root_label(root: TreeNode, node_metadata_rows: list[dict[str, str]]) -> None:
    if root.name:
        return
    used_labels = {node.name for node in iter_nodes(root) if node.name}
    for row in node_metadata_rows:
        candidate = row.get("tree_node_label", "")
        if row.get("node_type", "") != "internal":
            continue
        if candidate and candidate not in used_labels:
            root.name = candidate
            return
    root.name = "root"


def attach_parents(node: TreeNode, parent: TreeNode | None = None) -> None:
    node.parent = parent
    for child in node.children:
        attach_parents(child, node)


def assign_x_positions(node: TreeNode, cladogram_depth: int = 0, branch_depth: float = 0.0) -> None:
    node.x_cladogram = float(cladogram_depth)
    node.x_branch_length = branch_depth
    for child in node.children:
        assign_x_positions(child, cladogram_depth + 1, branch_depth + child.length)


def assign_y_positions(node: TreeNode, next_tip_order: int = 1) -> int:
    if node.is_tip:
        node.y_order = float(next_tip_order)
        return next_tip_order + 1

    current_tip_order = next_tip_order
    for child in node.children:
        current_tip_order = assign_y_positions(child, current_tip_order)
    node.y_order = sum(child.y_order for child in node.children) / len(node.children)
    return current_tip_order


def normalize_origin_node_id(raw_node_id: str) -> str:
    if raw_node_id.startswith("node_"):
        return raw_node_id.replace("node_", "m5node_", 1)
    return raw_node_id


def format_number(value: float) -> str:
    return f"{value:.6f}"


def build_arg_parser() -> argparse.ArgumentParser:
    root = repo_root()
    parser = argparse.ArgumentParser(description="Build workflow-native Figure 3 tree layout tables.")
    parser.add_argument(
        "--tree",
        type=Path,
        default=root / "outputs/workflow/asr_sensitivity/composition_filtered/rooted_ml_tree.reference_rooted.nwk",
        help="Rooted workflow ML tree in Newick format. Defaults to the composition-pruned primary ASR quality frame.",
    )
    parser.add_argument(
        "--node-metadata",
        type=Path,
        default=root / "outputs/workflow/asr_sensitivity/composition_filtered/rooted_tree_node_metadata.tsv",
        help="Per-node metadata emitted by the rooted-tree preprocessing step.",
    )
    parser.add_argument(
        "--tip-states",
        type=Path,
        default=root / "outputs/workflow/asr_sensitivity/composition_filtered/tip_states.tsv",
        help="Workflow-native ASR tip-state table.",
    )
    parser.add_argument(
        "--origin-events",
        type=Path,
        default=root / "outputs/workflow/asr_sensitivity/composition_filtered/origin_events.tsv",
        help="Workflow-native Fitch origin-event table.",
    )
    parser.add_argument(
        "--event-subtrees-dir",
        type=Path,
        default=root / "outputs/workflow/asr_sensitivity/composition_filtered/event_subtrees",
        help="Directory containing per-origin descendant tip packages.",
    )
    parser.add_argument(
        "--segments-out",
        type=Path,
        default=root / "manuscript/figure_data/figure3_workflow_tree_segments.tsv",
        help="Output path for the Figure 3 tree segment table.",
    )
    parser.add_argument(
        "--nodes-out",
        type=Path,
        default=root / "manuscript/figure_data/figure3_workflow_tree_nodes.tsv",
        help="Output path for the Figure 3 node-position table.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    tree_text = args.tree.read_text(encoding="utf-8").strip()
    root = NewickParser(tree_text).parse()

    node_metadata_rows = load_tsv(args.node_metadata)
    tip_state_rows = load_tsv(args.tip_states)
    origin_rows = load_tsv(args.origin_events)

    assign_root_label(root, node_metadata_rows)
    attach_parents(root)
    assign_x_positions(root)
    assign_y_positions(root)

    descendant_tip_sets: dict[str, frozenset[str]] = {}

    def compute_descendant_tip_sets(node: TreeNode) -> frozenset[str]:
        if node.is_tip:
            tip_set = frozenset([node.name])
        else:
            child_tip_sets = [compute_descendant_tip_sets(child) for child in node.children]
            merged: set[str] = set()
            for child_tip_set in child_tip_sets:
                merged.update(child_tip_set)
            tip_set = frozenset(merged)
        descendant_tip_sets[node.name] = tip_set
        return tip_set

    compute_descendant_tip_sets(root)

    node_metadata_by_id = {
        row.get("tree_node_label", ""): row
        for row in node_metadata_rows
        if row.get("tree_node_label", "")
    }
    tip_state_by_id = {
        row.get("tree_tip_label", ""): row
        for row in tip_state_rows
        if row.get("tree_tip_label", "")
    }
    def load_origin_tip_labels(origin_id: str) -> set[str]:
        path = args.event_subtrees_dir / f"{origin_id}.descendant_tips.tsv"
        if not path.exists():
            return set()
        rows = load_tsv(path)
        return {row.get("tip_label", "") for row in rows if row.get("tip_label", "")}

    def find_matching_tree_node(origin_row: dict[str, str]) -> str:
        normalized_id = normalize_origin_node_id(origin_row.get("clade_id", ""))
        if normalized_id in descendant_tip_sets:
            return normalized_id

        tip_labels = load_origin_tip_labels(origin_row.get("origin_id", ""))
        if not tip_labels:
            return normalized_id

        target_tip_set = frozenset(tip_labels)
        for node_id, tip_set in descendant_tip_sets.items():
            if tip_set == target_tip_set:
                return node_id
        return normalized_id

    origin_by_node_id = {}
    for row in origin_rows:
        mapped_node_id = find_matching_tree_node(row)
        if not mapped_node_id:
            continue
        row_with_mapping = dict(row)
        row_with_mapping["mapped_tree_node_id"] = mapped_node_id
        origin_by_node_id[mapped_node_id] = row_with_mapping

    segment_rows: list[dict[str, str]] = []
    node_rows: list[dict[str, str]] = []

    for node in iter_nodes(root):
        metadata_row = node_metadata_by_id.get(node.name, {})
        tip_row = tip_state_by_id.get(node.name, {})
        origin_row = origin_by_node_id.get(node.name, {})

        node_rows.append(
            {
                "node_id": node.name,
                "parent_node_id": node.parent.name if node.parent else "",
                "node_type": "tip" if node.is_tip else "internal",
                "tip_label": node.name if node.is_tip else "",
                "x_cladogram": format_number(node.x_cladogram),
                "x_branch_length": format_number(node.x_branch_length),
                "y_order": format_number(node.y_order),
                "branch_support": metadata_row.get("branch_support", ""),
                "is_reference": metadata_row.get("is_reference", "False"),
                "is_fitch_origin": "True" if origin_row else "False",
                "origin_id": origin_row.get("origin_id", ""),
                "n_tips_total": origin_row.get("n_tips_total", ""),
                "n_tips_disrupted": origin_row.get("n_tips_disrupted", ""),
                "n_countries": origin_row.get("n_countries", ""),
                "origin_first_year": origin_row.get("first_year", ""),
                "origin_last_year": origin_row.get("last_year", ""),
                "dominant_prn_mechanism": origin_row.get("dominant_prn_mechanism", ""),
                "origin_support_score": origin_row.get("origin_support_score", ""),
                "observed_prn_state": tip_row.get("prn_state", ""),
                "prn_mechanism_call": tip_row.get("prn_mechanism_call", "")
                or tip_row.get("observed_prn_mechanism_call", ""),
                "country_iso3": tip_row.get("country_iso3", ""),
                "year": tip_row.get("year", ""),
                "mlst_st": tip_row.get("mlst_st", ""),
                "phylo_lineage": tip_row.get("phylo_lineage", ""),
            }
        )

        if node.is_tip:
            continue

        child_positions = [child.y_order for child in node.children]
        segment_rows.append(
            {
                "segment_id": f"{node.name}.vertical",
                "segment_type": "vertical",
                "parent_node_id": node.name,
                "child_node_id": "",
                "x_cladogram": format_number(node.x_cladogram),
                "xend_cladogram": format_number(node.x_cladogram),
                "y_order": format_number(min(child_positions)),
                "yend_order": format_number(max(child_positions)),
            }
        )
        for child in node.children:
            segment_rows.append(
                {
                    "segment_id": f"{node.name}->{child.name}",
                    "segment_type": "horizontal",
                    "parent_node_id": node.name,
                    "child_node_id": child.name,
                    "x_cladogram": format_number(node.x_cladogram),
                    "xend_cladogram": format_number(child.x_cladogram),
                    "y_order": format_number(child.y_order),
                    "yend_order": format_number(child.y_order),
                }
            )

    write_tsv(args.segments_out, TREE_SEGMENT_COLUMNS, segment_rows)
    write_tsv(args.nodes_out, TREE_NODE_COLUMNS, node_rows)
    print(f"Wrote {len(segment_rows)} tree segments -> {args.segments_out}")
    print(f"Wrote {len(node_rows)} tree nodes -> {args.nodes_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
