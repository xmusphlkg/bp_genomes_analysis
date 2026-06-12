#!/usr/bin/env python3
"""Prune specified tips from a Newick tree."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Node:
    children: list["Node"] = field(default_factory=list)
    name: str = ""
    length: str = ""

    @property
    def is_leaf(self) -> bool:
        return not self.children


class NewickParser:
    def __init__(self, text: str) -> None:
        self.text = text
        self.index = 0

    def parse(self) -> Node:
        node = self._parse_subtree()
        self._skip_whitespace()
        if self._peek() == ";":
            self.index += 1
        self._skip_whitespace()
        if self.index != len(self.text):
            raise ValueError(f"unexpected trailing content at position {self.index}")
        return node

    def _parse_subtree(self) -> Node:
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
            return Node(children=children, name=name, length=length)

        name = self._parse_label()
        length = self._parse_length()
        return Node(name=name, length=length)

    def _parse_label(self) -> str:
        self._skip_whitespace()
        start = self.index
        while self.index < len(self.text) and self.text[self.index] not in ":,();":
            self.index += 1
        return self.text[start:self.index].strip()

    def _parse_length(self) -> str:
        self._skip_whitespace()
        if self._peek() != ":":
            return ""
        self.index += 1
        start = self.index
        while self.index < len(self.text) and self.text[self.index] not in ",();":
            self.index += 1
        return self.text[start:self.index].strip()

    def _skip_whitespace(self) -> None:
        while self.index < len(self.text) and self.text[self.index].isspace():
            self.index += 1

    def _peek(self) -> str:
        if self.index >= len(self.text):
            return ""
        return self.text[self.index]


def load_tree(path: Path) -> Node:
    text = path.read_text(encoding="utf-8").strip()
    if not text:
        raise ValueError(f"empty tree file: {path}")
    return NewickParser(text).parse()


def combine_branch_lengths(parent_length: str, child_length: str) -> str:
    if not parent_length:
        return child_length
    if not child_length:
        return parent_length
    return f"{float(parent_length) + float(child_length):.16g}"


def prune_tree(node: Node, exclude_tips: set[str], *, is_root: bool = False) -> Node | None:
    if node.is_leaf:
        return None if node.name in exclude_tips else Node(name=node.name, length=node.length)

    kept_children: list[Node] = []
    for child in node.children:
        pruned_child = prune_tree(child, exclude_tips)
        if pruned_child is not None:
            kept_children.append(pruned_child)

    if not kept_children:
        return None

    if len(kept_children) == 1 and not is_root:
        collapsed = kept_children[0]
        collapsed.length = combine_branch_lengths(node.length, collapsed.length)
        if collapsed.children and not collapsed.name and node.name:
            collapsed.name = node.name
        return collapsed

    return Node(children=kept_children, name=node.name, length=node.length)


def count_leaves(node: Node) -> int:
    if node.is_leaf:
        return 1
    return sum(count_leaves(child) for child in node.children)


def to_newick(node: Node) -> str:
    children = ""
    if node.children:
        children = f"({','.join(to_newick(child) for child in node.children)})"
    name = node.name
    length = f":{node.length}" if node.length else ""
    return f"{children}{name}{length}"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prune one or more tips from a Newick tree.")
    parser.add_argument("--tree", type=Path, required=True, help="Input tree path.")
    parser.add_argument("--exclude-list", type=Path, required=True, help="One tip label per line to exclude.")
    parser.add_argument("--out-tree", type=Path, required=True, help="Output pruned tree path.")
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    tree = load_tree(args.tree)
    exclude_tips = {
        line.strip()
        for line in args.exclude_list.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }

    pruned = prune_tree(tree, exclude_tips, is_root=True)
    if pruned is None or count_leaves(pruned) < 2:
        raise ValueError("tree pruning would leave fewer than two tips")

    args.out_tree.parent.mkdir(parents=True, exist_ok=True)
    args.out_tree.write_text(f"{to_newick(pruned)};\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())