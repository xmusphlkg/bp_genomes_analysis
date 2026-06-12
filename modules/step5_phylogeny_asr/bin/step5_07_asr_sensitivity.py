#!/usr/bin/env python3
"""
step5_07_asr_sensitivity.py — Sensitivity analysis for ancestral state reconstruction.

Tests robustness of the 34 independent prn- origin count to:
1. Tie-breaking rule (intact-biased vs disrupted-biased vs random)
2. State space (4-state vs binary collapse)
3. Random tie-breaking replicates for uncertainty quantification
"""
from __future__ import annotations

import csv
import json
import random
import sys
from collections import Counter
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR.parent / "outputs"
ANCESTRAL_STATES_PATH = OUTPUT_DIR / "bp_prn_ancestral_states.tsv"
ORIGINS_PATH = OUTPUT_DIR / "bp_prn_independent_origins.tsv"
TREE_PATH = OUTPUT_DIR / "bp_global_phylogeny.nwk"

INTACT = "intact"
DISRUPTED = "disrupted"
INSUFFICIENT = "insufficient_data"


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_newick_labels(tree_path: Path) -> dict[str, str]:
    """Extract tip label -> node mapping from Newick tree.
    Returns dict of {label: node_index} where node_index is position in tree.
    """
    s = tree_path.read_text().strip().rstrip(";").strip()

    labels = {}
    idx = [0]
    node_counter = [0]

    def parse():
        nid = node_counter[0]
        node_counter[0] += 1

        if idx[0] < len(s) and s[idx[0]] == '(':
            idx[0] += 1
            children = []
            while True:
                cid = parse()
                children.append(cid)
                if idx[0] < len(s) and s[idx[0]] == ',':
                    idx[0] += 1
                else:
                    break
            if idx[0] < len(s) and s[idx[0]] == ')':
                idx[0] += 1

        # Read label
        start = idx[0]
        while idx[0] < len(s) and s[idx[0]] not in ':,)':
            idx[0] += 1
        label = s[start:idx[0]].strip()
        if label:
            labels[label] = nid

        # Skip branch length
        if idx[0] < len(s) and s[idx[0]] == ':':
            idx[0] += 1
            while idx[0] < len(s) and s[idx[0]] not in ',)':
                idx[0] += 1

        return nid

    parse()
    return labels


def fitch_down(node_states, children_lists, n_nodes):
    """Bottom-up Fitch pass. Returns list of candidate sets per node."""
    candidates = [None] * n_nodes

    def recurse(nid):
        children = children_lists[nid]
        if not children:
            candidates[nid] = {node_states[nid]} if node_states[nid] else set()
            return

        for c in children:
            recurse(c)

        csets = [candidates[c] for c in children if candidates[c]]
        if not csets:
            candidates[nid] = set()
            return

        intersection = csets[0]
        for cs in csets[1:]:
            intersection = intersection & cs

        if intersection:
            candidates[nid] = intersection
        else:
            union = set()
            for cs in csets:
                union |= cs
            candidates[nid] = union

    recurse(0)
    return candidates


def fitch_up(candidates, children_lists, tie_break="intact", rng=None):
    """Top-down resolution pass."""
    resolved = [None] * len(candidates)

    def resolve(nid, parent_state=None):
        cands = candidates[nid]
        if not cands:
            resolved[nid] = INSUFFICIENT
        elif parent_state is not None and parent_state in cands:
            resolved[nid] = parent_state
        elif tie_break == "intact":
            if INTACT in cands:
                resolved[nid] = INTACT
            elif DISRUPTED in cands:
                resolved[nid] = DISRUPTED
            else:
                resolved[nid] = INSUFFICIENT
        elif tie_break == "disrupted":
            if DISRUPTED in cands:
                resolved[nid] = DISRUPTED
            elif INTACT in cands:
                resolved[nid] = INTACT
            else:
                resolved[nid] = INSUFFICIENT
        elif tie_break == "random":
            resolved[nid] = rng.choice(sorted(cands))
        else:
            resolved[nid] = sorted(cands)[0]

        for c in children_lists[nid]:
            resolve(c, resolved[nid])

    resolve(0)
    return resolved


def simple_parse_newick(tree_str):
    """Parse Newick into adjacency list and leaf labels."""
    s = tree_str.strip().rstrip(";").strip()
    idx = [0]
    node_counter = [0]
    children = [[]]
    labels = [None]

    def parse():
        nid = node_counter[0]
        node_counter[0] += 1
        children.append([])
        labels.append(None)

        if idx[0] < len(s) and s[idx[0]] == '(':
            idx[0] += 1
            while True:
                cid = parse()
                children[nid].append(cid)
                if idx[0] < len(s) and s[idx[0]] == ',':
                    idx[0] += 1
                else:
                    break
            if idx[0] < len(s) and s[idx[0]] == ')':
                idx[0] += 1

        start = idx[0]
        while idx[0] < len(s) and s[idx[0]] not in ':,)':
            idx[0] += 1
        lbl = s[start:idx[0]].strip()
        if lbl:
            labels[nid] = lbl

        if idx[0] < len(s) and s[idx[0]] == ':':
            idx[0] += 1
            while idx[0] < len(s) and s[idx[0]] not in ',)':
                idx[0] += 1

        return nid

    root = parse()
    return root, children, labels


def count_transitions(resolved, children, labels, tip_states_map):
    """Count intact->disrupted transitions."""
    origins = 0

    def scan(nid, parent_state=None):
        nonlocal origins
        state = resolved[nid]
        if parent_state == INTACT and state == DISRUPTED:
            origins += 1
        for c in children[nid]:
            scan(c, state)

    scan(0)
    return origins


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load tree
    tree_str = TREE_PATH.read_text().strip()
    root, children, labels = simple_parse_newick(tree_str)
    n_nodes = len(children)
    print(f"Tree: {n_nodes} nodes")

    # Load tip states from ASR output
    asr_rows = load_tsv(ANCESTRAL_STATES_PATH)
    tip_states_map = {}
    for row in asr_rows:
        if row.get("node_type") == "tip":
            lbl = row.get("tip_label", "").strip()
            st = row.get("observed_prn_state", INSUFFICIENT).strip()
            if lbl:
                tip_states_map[lbl] = st

    # Map labels to node indices
    label_to_node = {}
    for i, lbl in enumerate(labels):
        if lbl:
            label_to_node[lbl] = i

    print(f"Tip states: {len(tip_states_map)}")
    print(f"Label->node mappings: {len(label_to_node)}")

    # Assign observed states to nodes
    node_states = [None] * n_nodes
    for lbl, st in tip_states_map.items():
        if lbl in label_to_node:
            node_states[label_to_node[lbl]] = st

    # Count tips with assigned states
    n_assigned = sum(1 for s in node_states if s is not None)
    state_dist = Counter(s for s in node_states if s is not None)
    print(f"Assigned states: {n_assigned}, distribution: {dict(state_dist)}")

    results = []

    # 1. Intact-biased (default)
    print("\n1. Intact-biased (default)...")
    cands = fitch_down(node_states, children, n_nodes)
    resolved = fitch_up(cands, children, "intact")
    n_origins = count_transitions(resolved, children, labels, tip_states_map)
    results.append({"method": "intact_biased_default", "n_origins": n_origins})
    print(f"   Origins: {n_origins}")

    # 2. Disrupted-biased
    print("2. Disrupted-biased...")
    resolved = fitch_up(cands, children, "disrupted")
    n_origins = count_transitions(resolved, children, labels, tip_states_map)
    results.append({"method": "disrupted_biased", "n_origins": n_origins})
    print(f"   Origins: {n_origins}")

    # 3. Random tie-breaking (100 replicates)
    print("3. Random tie-breaking (100 replicates)...")
    random_counts = []
    rng = random.Random(42)
    for rep in range(100):
        resolved = fitch_up(cands, children, "random", rng)
        n_origins = count_transitions(resolved, children, labels, tip_states_map)
        random_counts.append(n_origins)

    results.append({
        "method": "random_100replicates",
        "n_origins_mean": round(float(np.mean(random_counts)), 1),
        "n_origins_sd": round(float(np.std(random_counts)), 1),
        "n_origins_min": int(np.min(random_counts)),
        "n_origins_max": int(np.max(random_counts)),
        "n_origins_median": round(float(np.median(random_counts)), 1),
    })
    print(f"   Mean={np.mean(random_counts):.1f}, SD={np.std(random_counts):.1f}, "
          f"range=[{np.min(random_counts)}, {np.max(random_counts)}]")

    # 4. Binary collapse (drop insufficient_data tips)
    print("4. Binary collapse (drop insufficient_data)...")
    node_states_binary = [None] * n_nodes
    for lbl, st in tip_states_map.items():
        if st in (INTACT, DISRUPTED) and lbl in label_to_node:
            node_states_binary[label_to_node[lbl]] = st

    n_binary = sum(1 for s in node_states_binary if s is not None)
    print(f"   Tips used: {n_binary} (dropped {n_assigned - n_binary} insufficient_data)")

    cands_b = fitch_down(node_states_binary, children, n_nodes)
    resolved = fitch_up(cands_b, children, "intact")
    n_origins = count_transitions(resolved, children, labels, tip_states_map)
    results.append({"method": "binary_intact_disrupted", "n_origins": n_origins,
                    "n_tips_used": n_binary})
    print(f"   Origins: {n_origins}")

    # 5. Random with binary
    print("5. Random tie-breaking, binary (100 replicates)...")
    rng2 = random.Random(99)
    random_counts_b = []
    for rep in range(100):
        resolved = fitch_up(cands_b, children, "random", rng2)
        n_origins = count_transitions(resolved, children, labels, tip_states_map)
        random_counts_b.append(n_origins)

    results.append({
        "method": "random_binary_100replicates",
        "n_origins_mean": round(float(np.mean(random_counts_b)), 1),
        "n_origins_sd": round(float(np.std(random_counts_b)), 1),
        "n_origins_min": int(np.min(random_counts_b)),
        "n_origins_max": int(np.max(random_counts_b)),
        "n_origins_median": round(float(np.median(random_counts_b)), 1),
    })
    print(f"   Mean={np.mean(random_counts_b):.1f}, SD={np.std(random_counts_b):.1f}, "
          f"range=[{np.min(random_counts_b)}, {np.max(random_counts_b)}]")

    # Save summary
    summary_path = OUTPUT_DIR / "bp_asr_sensitivity_summary.tsv"
    fieldnames = ["method", "n_origins", "n_origins_mean", "n_origins_sd",
                  "n_origins_min", "n_origins_max", "n_origins_median",
                  "n_tips_used"]
    with summary_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(results)
    print(f"\nSummary → {summary_path}")

    # Save JSON
    details = {
        "results": results,
        "tip_state_distribution": dict(Counter(tip_states_map.values())),
        "n_tips_total": len(tip_states_map),
        "n_nodes": n_nodes,
    }
    details_path = OUTPUT_DIR / "bp_asr_sensitivity_details.json"
    with details_path.open("w") as f:
        json.dump(details, f, indent=2)
    print(f"Details → {details_path}")

    print("\n=== SENSITIVITY SUMMARY ===")
    for r in results:
        method = r["method"]
        if "n_origins_mean" in r:
            print(f"  {method}: {r['n_origins_mean']:.1f} ± {r['n_origins_sd']:.1f} "
                  f"[{r['n_origins_min']}, {r['n_origins_max']}]")
        else:
            print(f"  {method}: {r['n_origins']} origins")


if __name__ == "__main__":
    main()
