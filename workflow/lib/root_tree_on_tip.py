#!/usr/bin/env python3
"""Re-root a Newick tree and assign stable internal node labels."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree


METADATA_COLUMNS = [
    "tree_node_label",
    "node_type",
    "original_label",
    "branch_support",
    "is_reference",
    "rooting_mode",
    "rooting_target",
]


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


def load_tree(path: Path) -> Tree:
    return Phylo.read(str(path), "newick")


def drop_empty_artifact_leaves(tree: Tree) -> int:
    """Remove unnamed leaves occasionally introduced after rerooting."""
    empty_leaves = [node for node in tree.get_terminals() if not normalize_text(node.name)]
    for node in empty_leaves:
        tree.prune(node)
    return len(empty_leaves)


def clade_original_label(clade: Clade) -> str:
    if clade.name:
        return normalize_text(clade.name)
    if clade.confidence is not None:
        return parse_support_value(str(clade.confidence))
    return ""


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=METADATA_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Re-root a tree and assign stable node labels.")
    parser.add_argument("--tree", type=Path, required=True, help="Input Newick tree path.")
    parser.add_argument("--out-tree", type=Path, required=True, help="Output rooted tree path.")
    parser.add_argument("--out-metadata", type=Path, required=True, help="Output node metadata TSV.")
    parser.add_argument(
        "--rooting-mode",
        choices=["reference", "midpoint"],
        default="reference",
        help="Rooting mode. reference roots on --outgroup; midpoint roots at the midpoint outgroup.",
    )
    parser.add_argument(
        "--outgroup",
        default="Reference",
        help="Tip label to use for reference rooting and reference-tip metadata.",
    )
    parser.add_argument("--prefix", default="m5node_", help="Prefix for generated internal node labels.")
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    tree = load_tree(args.tree)

    rooting_target = ""
    if args.rooting_mode == "reference":
        match = tree.find_any(name=args.outgroup)
        if match is None:
            raise ValueError(f"Outgroup tip not found in tree: {args.outgroup}")
        rooting_target = args.outgroup
        tree.root_with_outgroup(match)
    elif args.rooting_mode == "midpoint":
        tree.root_at_midpoint()
        rooting_target = "midpoint"
    else:  # pragma: no cover - argparse enforces choices
        raise ValueError(f"Unsupported rooting mode: {args.rooting_mode}")

    drop_empty_artifact_leaves(tree)

    metadata_rows: list[dict[str, str]] = []
    internal_index = 1
    for node in tree.find_clades(order="preorder"):
        original_label = clade_original_label(node)
        if node.is_terminal():
            tree_node_label = original_label
        else:
            tree_node_label = f"{args.prefix}{internal_index:06d}"
            internal_index += 1
            node.name = tree_node_label
            node.confidence = None

        metadata_rows.append(
            {
                "tree_node_label": tree_node_label,
                "node_type": "tip" if node.is_terminal() else "internal",
                "original_label": original_label,
                "branch_support": "" if node.is_terminal() else parse_support_value(original_label),
                "is_reference": "True" if node.is_terminal() and tree_node_label == args.outgroup else "False",
                "rooting_mode": args.rooting_mode,
                "rooting_target": rooting_target,
            }
        )

    args.out_tree.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(args.out_tree), "newick")
    write_tsv(args.out_metadata, metadata_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
