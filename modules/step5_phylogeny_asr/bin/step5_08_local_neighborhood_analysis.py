#!/usr/bin/env python3
"""Build local-neighborhood tree checks for disrupted prn tips.

The rooted SNP tree remains the primary ASR frame. This script adds two audit
layers requested for submission hardening:

1. A rooted-SNP local-neighborhood sensitivity containing every disrupted tip in
   the 191-tip ASR tree plus nearest intact neighbors.
2. A full-manifest proxy-tree neighborhood containing every disrupted tip
   present in the larger k-mer proxy tree plus nearest intact neighbors.

The proxy-tree layer is deliberately reported as a topology/neighborhood check,
not as a replacement for the rooted SNP ASR origin count.
"""

from __future__ import annotations

import argparse
import csv
import heapq
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.lib.project_paths import project_module_data_root, project_workflow_root


STATE_ORDER = ["intact", "disrupted", "insufficient_data", "uncertain"]
LONGREAD_RE = re.compile(r"nanopore|pacbio|hifi|rsii|sequel|minion|promethion", re.I)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def normalize(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() == "nan" else text


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: normalize(row.get(field, "")) for field in fieldnames})


def state_from_mechanism(mechanism: str) -> str:
    mechanism = normalize(mechanism)
    if mechanism in {"intact", "reference_assumed_intact"}:
        return "intact"
    if mechanism.startswith("coding_disrupted_"):
        return "disrupted"
    if mechanism in {"insufficient_data", "uncertain_fragmented_assembly", "not_available_current_step3"}:
        return "insufficient_data"
    if not mechanism:
        return "insufficient_data"
    return "uncertain"


def first_present(row: dict[str, str], columns: list[str]) -> str:
    values = [normalize(row.get(column, "")) for column in columns]
    return ";".join(value for value in values if value)


def load_rooted_snp_tip_metadata(tip_states_path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(tip_states_path)
    metadata: dict[str, dict[str, str]] = {}
    for row in rows:
        tip = normalize(row.get("tree_tip_label", ""))
        if not tip:
            continue
        metadata[tip] = {
            "tree_tip_label": tip,
            "sample_id_canonical": normalize(row.get("sample_id_canonical", "")),
            "assembly_accession": normalize(row.get("assembly_accession", "")),
            "prn_state": normalize(row.get("prn_state", "")),
            "prn_mechanism_call": normalize(row.get("observed_prn_mechanism_call", "")),
            "prn_event_id": "",
            "prn_call_confidence": normalize(row.get("prn_call_confidence", "")),
            "country_iso3": normalize(row.get("country_iso3", "")),
            "year": normalize(row.get("year", "")),
            "read_accessions": "",
            "sequencing_tech": "",
            "is_public_longread_or_hybrid": "False",
        }
    return metadata


def load_proxy_tip_metadata(manifest_path: Path, mechanism_path: Path, tree: Tree) -> dict[str, dict[str, str]]:
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    mechanism = pd.read_csv(mechanism_path, sep="\t", dtype=str)
    keep_cols = [
        "sample_id_canonical",
        "prn_mechanism_call",
        "prn_event_id",
        "prn_call_confidence",
        "country_iso3",
        "year",
        "assembly_accession",
        "sra_run_accession",
    ]
    mechanism = mechanism[[column for column in keep_cols if column in mechanism.columns]]
    merged = manifest.merge(mechanism, on="sample_id_canonical", how="left", suffixes=("", "_mechanism"))

    tip_labels = {terminal.name for terminal in tree.get_terminals()}
    metadata: dict[str, dict[str, str]] = {}
    for _, row in merged[merged["sample_id_canonical"].isin(tip_labels)].iterrows():
        tip = normalize(row.get("sample_id_canonical"))
        mechanism_call = normalize(row.get("prn_mechanism_call"))
        sequencing_tech = normalize(row.get("sequencing_tech"))
        assembly_accession = normalize(row.get("assembly_accession_mechanism")) or normalize(row.get("assembly_accession"))
        sra_run = normalize(row.get("sra_run_accession_mechanism")) or normalize(row.get("sra_run_accession"))
        ena_run = normalize(row.get("ena_run_accession"))
        metadata[tip] = {
            "tree_tip_label": tip,
            "sample_id_canonical": tip,
            "assembly_accession": assembly_accession,
            "prn_state": state_from_mechanism(mechanism_call),
            "prn_mechanism_call": mechanism_call,
            "prn_event_id": normalize(row.get("prn_event_id")),
            "prn_call_confidence": normalize(row.get("prn_call_confidence")),
            "country_iso3": normalize(row.get("country_iso3_mechanism")) or normalize(row.get("country_iso3")),
            "year": normalize(row.get("year_mechanism")) or normalize(row.get("year")),
            "read_accessions": ";".join(value for value in [sra_run, ena_run] if value),
            "sequencing_tech": sequencing_tech,
            "is_public_longread_or_hybrid": "True" if LONGREAD_RE.search(sequencing_tech) else "False",
        }
    return metadata


def build_graph(tree: Tree) -> tuple[dict[int, list[tuple[int, float]]], dict[int, str], dict[str, int]]:
    adjacency: dict[int, list[tuple[int, float]]] = defaultdict(list)
    label_by_id: dict[int, str] = {}
    id_by_label: dict[str, int] = {}
    for clade in tree.find_clades(order="level"):
        if clade.is_terminal():
            label_by_id[id(clade)] = clade.name
            id_by_label[clade.name] = id(clade)
        for child in clade.clades:
            weight = float(child.branch_length or 0.0)
            adjacency[id(clade)].append((id(child), weight))
            adjacency[id(child)].append((id(clade), weight))
    return adjacency, label_by_id, id_by_label


def nearest_intact_neighbors(
    tree: Tree,
    metadata: dict[str, dict[str, str]],
    neighbor_k: int,
) -> list[dict[str, Any]]:
    adjacency, label_by_id, id_by_label = build_graph(tree)
    intact = {tip for tip, row in metadata.items() if row.get("prn_state") == "intact"}
    disrupted = sorted(tip for tip, row in metadata.items() if row.get("prn_state") == "disrupted")
    rows: list[dict[str, Any]] = []

    for disrupted_tip in disrupted:
        start = id_by_label.get(disrupted_tip)
        if start is None:
            continue
        heap: list[tuple[float, int]] = [(0.0, start)]
        seen = {start: 0.0}
        found: list[tuple[str, float]] = []
        while heap and len(found) < neighbor_k:
            distance, node_id = heapq.heappop(heap)
            if distance != seen[node_id]:
                continue
            label = label_by_id.get(node_id)
            if label in intact:
                found.append((label, distance))
                continue
            for neighbor_id, weight in adjacency[node_id]:
                new_distance = distance + weight
                if new_distance < seen.get(neighbor_id, math.inf):
                    seen[neighbor_id] = new_distance
                    heapq.heappush(heap, (new_distance, neighbor_id))
        for rank, (neighbor_label, distance) in enumerate(found, start=1):
            rows.append(
                {
                    "source_disrupted_tip": disrupted_tip,
                    "neighbor_rank": rank,
                    "neighbor_tip": neighbor_label,
                    "patristic_distance": f"{distance:.8f}",
                }
            )
    return rows


def choose_oldest_intact(metadata: dict[str, dict[str, str]]) -> tuple[str, str]:
    candidates: list[tuple[float, str]] = []
    for tip, row in metadata.items():
        if row.get("prn_state") != "intact":
            continue
        year_text = normalize(row.get("year"))
        try:
            year = float(year_text)
        except ValueError:
            year = math.inf
        candidates.append((year, tip))
    if not candidates:
        raise ValueError("no intact tips available for proxy-tree root anchor")
    year, tip = min(candidates)
    basis = "oldest_intact_tip_with_year" if math.isfinite(year) else "lexicographic_intact_tip_without_year"
    return tip, basis


def prune_tree(tree_path: Path, keep_labels: set[str], root_anchor: str | None) -> Tree:
    tree = Phylo.read(tree_path, "newick")
    for tip in list(terminal.name for terminal in tree.get_terminals()):
        if tip not in keep_labels:
            tree.prune(tip)
    if root_anchor:
        anchor = next((terminal for terminal in tree.get_terminals() if terminal.name == root_anchor), None)
        if anchor is not None:
            tree.root_with_outgroup(anchor)
    return tree


def preferred_state(states: set[str]) -> str:
    for state in STATE_ORDER:
        if state in states:
            return state
    return sorted(states)[0] if states else "uncertain"


def run_fitch(tree: Tree, metadata: dict[str, dict[str, str]]) -> tuple[int, list[dict[str, Any]]]:
    score = 0
    node_ids = {id(clade): f"node_{idx:06d}" for idx, clade in enumerate(tree.find_clades(order="preorder"), start=1)}

    def downpass(clade: Clade) -> set[str]:
        nonlocal score
        if clade.is_terminal():
            state = metadata.get(clade.name, {}).get("prn_state", "insufficient_data")
            clade._candidate_states = {state}  # type: ignore[attr-defined]
            clade._descendant_counts = Counter([state])  # type: ignore[attr-defined]
            return clade._candidate_states  # type: ignore[attr-defined]
        child_sets: list[set[str]] = []
        counts: Counter[str] = Counter()
        for child in clade.clades:
            child_sets.append(downpass(child))
            counts.update(child._descendant_counts)  # type: ignore[attr-defined]
        intersection = set.intersection(*child_sets)
        if intersection:
            clade._candidate_states = intersection  # type: ignore[attr-defined]
        else:
            clade._candidate_states = set.union(*child_sets)  # type: ignore[attr-defined]
            score += 1
        clade._descendant_counts = counts  # type: ignore[attr-defined]
        return clade._candidate_states  # type: ignore[attr-defined]

    def uppass(clade: Clade, parent_state: str | None = None) -> None:
        candidates = clade._candidate_states  # type: ignore[attr-defined]
        clade._inferred_state = parent_state if parent_state and parent_state in candidates else preferred_state(candidates)  # type: ignore[attr-defined]
        for child in clade.clades:
            uppass(child, clade._inferred_state)  # type: ignore[attr-defined]

    downpass(tree.root)
    uppass(tree.root)

    transition_rows: list[dict[str, Any]] = []
    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            if parent._inferred_state == "intact" and child._inferred_state == "disrupted":  # type: ignore[attr-defined]
                descendant_tips = [terminal.name for terminal in child.get_terminals()]
                disrupted_descendants = [
                    tip for tip in descendant_tips if metadata.get(tip, {}).get("prn_state") == "disrupted"
                ]
                transition_rows.append(
                    {
                        "transition_node_id": node_ids[id(child)],
                        "parent_node_id": node_ids[id(parent)],
                        "descendant_tip_count": len(descendant_tips),
                        "descendant_disrupted_tip_count": len(disrupted_descendants),
                        "descendant_disrupted_tip_examples": ";".join(disrupted_descendants[:10]),
                    }
                )
    return score, transition_rows


def summarize_states(metadata: dict[str, dict[str, str]], labels: set[str]) -> Counter[str]:
    return Counter(metadata.get(label, {}).get("prn_state", "insufficient_data") for label in labels)


def build_selection_rows(
    analysis_id: str,
    metadata: dict[str, dict[str, str]],
    selected: set[str],
    disrupted_targets: set[str],
    neighbor_rows: list[dict[str, Any]],
    root_anchor: str,
) -> list[dict[str, Any]]:
    neighbor_summary: dict[str, dict[str, Any]] = defaultdict(lambda: {"count": 0, "min_rank": math.inf, "min_distance": math.inf, "examples": []})
    for row in neighbor_rows:
        neighbor_tip = normalize(row.get("neighbor_tip"))
        if neighbor_tip not in selected:
            continue
        summary = neighbor_summary[neighbor_tip]
        summary["count"] += 1
        summary["min_rank"] = min(summary["min_rank"], int(row["neighbor_rank"]))
        summary["min_distance"] = min(summary["min_distance"], float(row["patristic_distance"]))
        if len(summary["examples"]) < 10:
            summary["examples"].append(row["source_disrupted_tip"])

    rows: list[dict[str, Any]] = []
    for tip in sorted(selected):
        roles: list[str] = []
        if tip in disrupted_targets:
            roles.append("disrupted_target")
        if tip in neighbor_summary:
            roles.append("nearest_intact_neighbor")
        if tip == root_anchor:
            roles.append("root_anchor")
        meta = metadata.get(tip, {})
        summary = neighbor_summary.get(tip, {})
        rows.append(
            {
                "analysis_id": analysis_id,
                "tree_tip_label": tip,
                "sample_id_canonical": meta.get("sample_id_canonical", ""),
                "assembly_accession": meta.get("assembly_accession", ""),
                "prn_state": meta.get("prn_state", ""),
                "prn_mechanism_call": meta.get("prn_mechanism_call", ""),
                "prn_event_id": meta.get("prn_event_id", ""),
                "country_iso3": meta.get("country_iso3", ""),
                "year": meta.get("year", ""),
                "selection_roles": ";".join(roles),
                "nearest_to_disrupted_tip_count": summary.get("count", ""),
                "nearest_neighbor_min_rank": "" if summary.get("min_rank", math.inf) == math.inf else summary.get("min_rank"),
                "nearest_neighbor_min_distance": ""
                if summary.get("min_distance", math.inf) == math.inf
                else f"{summary.get('min_distance'):.8f}",
                "source_disrupted_tip_examples": ";".join(summary.get("examples", [])),
                "is_public_longread_or_hybrid": meta.get("is_public_longread_or_hybrid", "False"),
                "read_accessions": meta.get("read_accessions", ""),
            }
        )
    return rows


def run_analysis(
    *,
    analysis_id: str,
    tree_path: Path,
    metadata: dict[str, dict[str, str]],
    tree_method: str,
    root_anchor: str,
    root_anchor_basis: str,
    neighbor_k: int,
    output_tree_path: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    tree = Phylo.read(tree_path, "newick")
    source_labels = {terminal.name for terminal in tree.get_terminals()}
    disrupted_targets = {tip for tip in source_labels if metadata.get(tip, {}).get("prn_state") == "disrupted"}
    neighbor_rows = nearest_intact_neighbors(tree, metadata, neighbor_k)
    selected = set(disrupted_targets)
    selected.update(normalize(row.get("neighbor_tip")) for row in neighbor_rows)
    selected = {tip for tip in selected if tip}
    selected.add(root_anchor)

    pruned = prune_tree(tree_path, selected, root_anchor)
    output_tree_path.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(pruned, output_tree_path, "newick")
    fitch_score, transition_rows = run_fitch(pruned, metadata)

    source_counts = summarize_states(metadata, source_labels)
    selected_counts = summarize_states(metadata, {terminal.name for terminal in pruned.get_terminals()})
    total_disrupted_retained = sum(
        1
        for row in pd.read_csv(
            project_module_data_root("step4_prn_validation") / "outputs" / "bp_prn_mechanism_calls.tsv",
            sep="\t",
            dtype=str,
        )[
            "prn_mechanism_call"
        ].fillna("")
        if row.startswith("coding_disrupted_")
    )

    origin_interpretation = (
        "rooted_snp_local_neighborhood_sensitivity"
        if analysis_id == "rooted_snp_k3_neighborhood"
        else "proxy_tree_topology_check_not_absolute_origin_count"
    )
    summary = {
        "analysis_id": analysis_id,
        "tree_source": tree_path,
        "tree_method": tree_method,
        "neighbor_k": neighbor_k,
        "source_tree_tip_count": len(source_labels),
        "source_intact_tip_count": source_counts.get("intact", 0),
        "source_disrupted_tip_count": source_counts.get("disrupted", 0),
        "source_insufficient_tip_count": source_counts.get("insufficient_data", 0),
        "retained_cohort_disrupted_count": total_disrupted_retained,
        "source_tree_disrupted_coverage_fraction": f"{source_counts.get('disrupted', 0) / total_disrupted_retained:.4f}"
        if total_disrupted_retained
        else "",
        "selected_tip_count": len(pruned.get_terminals()),
        "selected_intact_tip_count": selected_counts.get("intact", 0),
        "selected_disrupted_tip_count": selected_counts.get("disrupted", 0),
        "selected_insufficient_tip_count": selected_counts.get("insufficient_data", 0),
        "root_anchor_tip": root_anchor,
        "root_anchor_basis": root_anchor_basis,
        "binary_fitch_minimum_change_count": fitch_score,
        "intact_to_disrupted_transition_count": len(transition_rows),
        "origin_count_interpretation": origin_interpretation,
        "output_tree": output_tree_path,
        "notes": (
            "All disrupted tips from the primary rooted SNP tree were retained with nearest intact neighbors."
            if analysis_id == "rooted_snp_k3_neighborhood"
            else (
                "All disrupted tips present in the full-manifest k-mer proxy tree were retained with nearest intact "
                "neighbors; counts are topology diagnostics and are not used as replacement ASR origin estimates."
            )
        ),
    }

    selection_rows = build_selection_rows(
        analysis_id=analysis_id,
        metadata=metadata,
        selected={terminal.name for terminal in pruned.get_terminals()},
        disrupted_targets=disrupted_targets,
        neighbor_rows=neighbor_rows,
        root_anchor=root_anchor,
    )
    for row in transition_rows:
        row["analysis_id"] = analysis_id
        row["origin_count_interpretation"] = origin_interpretation
    return summary, selection_rows, transition_rows


def build_arg_parser() -> argparse.ArgumentParser:
    workflow_root = project_workflow_root()
    step5_root = project_module_data_root("step5_phylogeny_asr")
    step4_root = project_module_data_root("step4_prn_validation")
    parser = argparse.ArgumentParser(description="Build local-neighborhood tree checks for prn-disrupted tips.")
    parser.add_argument("--neighbor-k", type=int, default=3, help="Nearest intact neighbors retained per disrupted tip.")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=workflow_root / "asr_local_neighborhood",
        help="Workflow output directory for neighborhood trees and ledgers.",
    )
    parser.add_argument(
        "--summary-out",
        type=Path,
        default=step5_root / "outputs" / "bp_prn_local_neighborhood_summary.tsv",
    )
    parser.add_argument(
        "--selection-out",
        type=Path,
        default=step5_root / "outputs" / "bp_prn_local_neighborhood_tip_selection.tsv",
    )
    parser.add_argument(
        "--transitions-out",
        type=Path,
        default=step5_root / "outputs" / "bp_prn_local_neighborhood_origin_events.tsv",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    workflow_root = project_workflow_root()
    step5_root = project_module_data_root("step5_phylogeny_asr")
    step4_root = project_module_data_root("step4_prn_validation")

    rooted_tree_path = workflow_root / "asr" / "rooted_ml_tree.reference_rooted.nwk"
    rooted_metadata = load_rooted_snp_tip_metadata(workflow_root / "asr" / "tip_states.tsv")
    proxy_tree_path = step5_root / "outputs" / "bp_global_phylogeny_full_sensitivity.nwk"
    proxy_tree = Phylo.read(proxy_tree_path, "newick")
    proxy_metadata = load_proxy_tip_metadata(
        step5_root / "outputs" / "bp_phylogeny_manifest_full.tsv",
        step4_root / "outputs" / "bp_prn_mechanism_calls.tsv",
        proxy_tree,
    )
    proxy_root_anchor, proxy_root_basis = choose_oldest_intact(proxy_metadata)

    analyses = [
        run_analysis(
            analysis_id="rooted_snp_k3_neighborhood",
            tree_path=rooted_tree_path,
            metadata=rooted_metadata,
            tree_method="reference_rooted_snippy_gubbins_iqtree2_ml_tree",
            root_anchor="Reference",
            root_anchor_basis="reference_tip",
            neighbor_k=args.neighbor_k,
            output_tree_path=args.out_dir / "rooted_snp_k3.neighborhood.nwk",
        ),
        run_analysis(
            analysis_id="full_manifest_proxy_k3_neighborhood",
            tree_path=proxy_tree_path,
            metadata=proxy_metadata,
            tree_method="sampled_canonical_kmer_composition_average_linkage_proxy_tree",
            root_anchor=proxy_root_anchor,
            root_anchor_basis=proxy_root_basis,
            neighbor_k=args.neighbor_k,
            output_tree_path=args.out_dir / "full_manifest_proxy_k3.neighborhood.nwk",
        ),
    ]

    summary_rows: list[dict[str, Any]] = []
    selection_rows: list[dict[str, Any]] = []
    transition_rows: list[dict[str, Any]] = []
    for summary, selection, transitions in analyses:
        summary_rows.append(summary)
        selection_rows.extend(selection)
        transition_rows.extend(transitions)

    summary_fields = [
        "analysis_id",
        "tree_source",
        "tree_method",
        "neighbor_k",
        "source_tree_tip_count",
        "source_intact_tip_count",
        "source_disrupted_tip_count",
        "source_insufficient_tip_count",
        "retained_cohort_disrupted_count",
        "source_tree_disrupted_coverage_fraction",
        "selected_tip_count",
        "selected_intact_tip_count",
        "selected_disrupted_tip_count",
        "selected_insufficient_tip_count",
        "root_anchor_tip",
        "root_anchor_basis",
        "binary_fitch_minimum_change_count",
        "intact_to_disrupted_transition_count",
        "origin_count_interpretation",
        "output_tree",
        "notes",
    ]
    selection_fields = [
        "analysis_id",
        "tree_tip_label",
        "sample_id_canonical",
        "assembly_accession",
        "prn_state",
        "prn_mechanism_call",
        "prn_event_id",
        "country_iso3",
        "year",
        "selection_roles",
        "nearest_to_disrupted_tip_count",
        "nearest_neighbor_min_rank",
        "nearest_neighbor_min_distance",
        "source_disrupted_tip_examples",
        "is_public_longread_or_hybrid",
        "read_accessions",
    ]
    transition_fields = [
        "analysis_id",
        "transition_node_id",
        "parent_node_id",
        "descendant_tip_count",
        "descendant_disrupted_tip_count",
        "descendant_disrupted_tip_examples",
        "origin_count_interpretation",
    ]
    write_tsv(args.summary_out, summary_fields, summary_rows)
    write_tsv(args.selection_out, selection_fields, selection_rows)
    write_tsv(args.transitions_out, transition_fields, transition_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
