#!/usr/bin/env python3
"""Compare IQ-TREE2 and RAxML-NG outputs and emit a compact JSON report."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def get_tip_names(tree_path: Path) -> list[str]:
    newick = tree_path.read_text().strip()
    tips = re.findall(r"(?<=[(,])([^():;,\s]+)", newick)
    return sorted(set(tips))


def compare(iqtree_path: Path, raxml_path: Path) -> dict:
    iq_tips = get_tip_names(iqtree_path)
    rax_tips = get_tip_names(raxml_path)
    iq_set = set(iq_tips)
    rax_set = set(rax_tips)
    shared = sorted(iq_set & rax_set)

    report = {
        "iqtree_path": str(iqtree_path),
        "raxml_path": str(raxml_path),
        "iqtree_tip_count": len(iq_tips),
        "raxml_tip_count": len(rax_tips),
        "shared_tip_count": len(shared),
        "tips_match_exactly": iq_set == rax_set,
        "iqtree_only_tips": sorted(iq_set - rax_set),
        "raxml_only_tips": sorted(rax_set - iq_set),
    }

    try:
        from ete3 import Tree

        iq_tree = Tree(str(iqtree_path), format=1)
        rax_tree = Tree(str(raxml_path), format=1)
        rf, rf_max, *_ = iq_tree.robinson_foulds(raxml_tree := rax_tree, unrooted_trees=True)
        report["robinson_foulds"] = rf
        report["robinson_foulds_max"] = rf_max
        report["robinson_foulds_normalized"] = round(rf / rf_max, 6) if rf_max else 0.0
    except Exception as exc:  # pragma: no cover - optional dependency/runtime path
        report["robinson_foulds"] = None
        report["robinson_foulds_max"] = None
        report["robinson_foulds_normalized"] = None
        report["rf_error"] = str(exc)

    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare two Newick trees")
    parser.add_argument("--iqtree", required=True)
    parser.add_argument("--raxml", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    report = compare(Path(args.iqtree), Path(args.raxml))
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps(report, indent=2))


if "snakemake" in dir():
    report = compare(Path(str(snakemake.input.iq_tree)), Path(str(snakemake.input.rax_tree)))
    with Path(str(snakemake.output.comparison)).open("w") as handle:
        json.dump(report, handle, indent=2)


if __name__ == "__main__":
    main()