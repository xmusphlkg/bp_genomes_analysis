#!/usr/bin/env python3
"""Build Mk-model uncertainty summaries for prn origin counts.

This script complements the Fitch + PastML sensitivity bundle with a binary
continuous-time Markov model on the existing rooted trees. Tip states are
collapsed to intact/disrupted, while insufficient or uncertain tips are
treated as missing observations. For each rooted analysis frame, the script:

1. fits an equal-rates two-state Mk model to the observed tip states;
2. computes the exact posterior expected number of intact->disrupted edges;
3. samples node-state histories conditional on the fitted model and observed
   tips to obtain a Monte Carlo interval for the number of origin edges.

The resulting table does not replace the primary Fitch/PastML readouts. It
adds a model-based interval so the manuscript can state more explicitly that
origin counts are analysis-frame summaries rather than fixed global truths.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

import numpy as np
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from scipy.optimize import minimize_scalar


ROOT = Path(__file__).resolve().parents[3]
SUPP_DIR = ROOT / "manuscript" / "supplementary"
FIGURE_DATA_DIR = ROOT / "manuscript" / "figure_data"

STATE_TO_INDEX = {"intact": 0, "disrupted": 1}
INDEX_TO_STATE = {0: "intact", 1: "disrupted"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.casefold() in {"", "nan", "none", "na"}:
        return ""
    return text


def fmt(value: Any, digits: int = 6) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    if not np.isfinite(number):
        return ""
    return f"{number:.{digits}f}"


def load_tip_states(path: Path) -> tuple[dict[str, str | None], int, int, int]:
    rows = read_tsv(path)
    tip_state_map: dict[str, str | None] = {}
    observed_intact = 0
    observed_disrupted = 0
    missing_tip_states = 0
    for row in rows:
        tip = clean_text(row.get("tree_tip_label"))
        state = clean_text(row.get("prn_state"))
        if not tip:
            continue
        if state in STATE_TO_INDEX:
            tip_state_map[tip] = state
            if state == "intact":
                observed_intact += 1
            else:
                observed_disrupted += 1
        else:
            tip_state_map[tip] = None
            missing_tip_states += 1
    return tip_state_map, observed_intact, observed_disrupted, missing_tip_states


def transition_matrix(branch_length: float, rate: float) -> np.ndarray:
    t = max(float(branch_length or 0.0), 0.0)
    q = max(float(rate), 1.0e-12)
    decay = math.exp(-2.0 * q * t)
    stay = 0.5 + 0.5 * decay
    switch = 0.5 - 0.5 * decay
    return np.asarray([[stay, switch], [switch, stay]], dtype=float)


def normalize_log_probs(log_probs: np.ndarray) -> np.ndarray:
    finite = np.isfinite(log_probs)
    if not finite.any():
        return np.asarray([0.5, 0.5], dtype=float)
    max_log = float(np.max(log_probs[finite]))
    probs = np.zeros_like(log_probs, dtype=float)
    probs[finite] = np.exp(log_probs[finite] - max_log)
    total = float(probs.sum())
    if total <= 0.0 or not np.isfinite(total):
        return np.asarray([0.5, 0.5], dtype=float)
    return probs / total


def build_transition_lookup(tree: Tree, rate: float) -> dict[Clade, np.ndarray]:
    lookup: dict[Clade, np.ndarray] = {}
    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            lookup[child] = np.clip(transition_matrix(child.branch_length or 0.0, rate), 1.0e-300, 1.0)
    return lookup


def compute_likelihoods(
    tree: Tree,
    tip_state_map: dict[str, str | None],
    rate: float,
) -> tuple[dict[Clade, np.ndarray], float]:
    log_likelihoods: dict[Clade, np.ndarray] = {}

    def recurse(clade: Clade) -> np.ndarray:
        if clade.is_terminal():
            observed = tip_state_map.get(clean_text(clade.name))
            if observed is None:
                log_raw = np.zeros(2, dtype=float)
            else:
                log_raw = np.full(2, -np.inf, dtype=float)
                log_raw[STATE_TO_INDEX[observed]] = 0.0
            log_likelihoods[clade] = log_raw
            return log_raw

        child_log_likelihoods = {child: recurse(child) for child in clade.clades}
        log_raw = np.zeros(2, dtype=float)
        for parent_state in range(2):
            total = 0.0
            for child in clade.clades:
                child_log_like = child_log_likelihoods[child]
                trans = np.clip(transition_matrix(child.branch_length or 0.0, rate), 1.0e-300, 1.0)
                child_terms = np.log(trans[parent_state, :]) + child_log_like
                child_finite = np.isfinite(child_terms)
                if not child_finite.any():
                    total = -np.inf
                    break
                child_max = float(np.max(child_terms[child_finite]))
                total += child_max + math.log(float(np.exp(child_terms[child_finite] - child_max).sum()))
            log_raw[parent_state] = total

        if not np.isfinite(log_raw).any():
            raise ValueError("Encountered non-finite node likelihoods while fitting Mk model.")
        log_likelihoods[clade] = log_raw
        return log_raw

    root_log_like = recurse(tree.root)
    root_prior = np.log(np.asarray([0.5, 0.5], dtype=float))
    root_terms = root_prior + root_log_like
    root_max = float(np.max(root_terms[np.isfinite(root_terms)]))
    tree_loglik = root_max + math.log(float(np.exp(root_terms - root_max).sum()))
    return log_likelihoods, tree_loglik


def fit_equal_rates_mk(tree: Tree, tip_state_map: dict[str, str | None]) -> tuple[float, float]:
    def objective(log_rate: float) -> float:
        rate = float(math.exp(log_rate))
        try:
            _, loglik = compute_likelihoods(tree, tip_state_map, rate)
        except ValueError:
            return math.inf
        return -loglik

    result = minimize_scalar(objective, bounds=(-18.0, 4.0), method="bounded", options={"xatol": 1.0e-2})
    if not result.success:
        raise RuntimeError(f"Mk equal-rates optimisation failed: {result.message}")
    rate = float(math.exp(result.x))
    return rate, float(-result.fun)


def sample_origin_distribution(
    tree: Tree,
    log_likelihoods: dict[Clade, np.ndarray],
    rate: float,
    *,
    n_draws: int,
    seed: int,
) -> tuple[list[int], float]:
    rng = np.random.default_rng(seed)
    transition_by_child = build_transition_lookup(tree, rate)
    root_log_posterior = np.log(np.asarray([0.5, 0.5], dtype=float)) + log_likelihoods[tree.root]
    root_posterior = normalize_log_probs(root_log_posterior)
    root_disrupted_prob = float(root_posterior[STATE_TO_INDEX["disrupted"]])

    origin_counts: list[int] = []
    for _ in range(n_draws):
        count = 0
        root_state = int(rng.choice([0, 1], p=root_posterior))
        stack: list[tuple[Clade, int]] = [(tree.root, root_state)]
        while stack:
            parent, parent_state = stack.pop()
            for child in parent.clades:
                trans = transition_by_child[child]
                child_log_probs = np.log(trans[parent_state, :]) + log_likelihoods[child]
                child_probs = normalize_log_probs(child_log_probs)
                child_state = int(rng.choice([0, 1], p=child_probs))
                if parent_state == STATE_TO_INDEX["intact"] and child_state == STATE_TO_INDEX["disrupted"]:
                    count += 1
                stack.append((child, child_state))
        origin_counts.append(count)
    return origin_counts, root_disrupted_prob


def expected_origin_count(
    tree: Tree,
    log_likelihoods: dict[Clade, np.ndarray],
    rate: float,
) -> tuple[float, dict[Clade, np.ndarray]]:
    transition_by_child = build_transition_lookup(tree, rate)
    root_log_posterior = np.log(np.asarray([0.5, 0.5], dtype=float)) + log_likelihoods[tree.root]
    root_posterior = normalize_log_probs(root_log_posterior)
    posterior_by_clade: dict[Clade, np.ndarray] = {tree.root: root_posterior}
    expected = 0.0

    stack = [tree.root]
    while stack:
        parent = stack.pop()
        parent_post = posterior_by_clade[parent]
        for child in parent.clades:
            trans = transition_by_child[child]
            conditional = np.zeros((2, 2), dtype=float)
            child_post = np.zeros(2, dtype=float)
            for parent_state in range(2):
                child_log_probs = np.log(trans[parent_state, :]) + log_likelihoods[child]
                probs = normalize_log_probs(child_log_probs)
                conditional[parent_state, :] = probs
                child_post += parent_post[parent_state] * probs
            posterior_by_clade[child] = child_post
            expected += float(parent_post[STATE_TO_INDEX["intact"]] * conditional[STATE_TO_INDEX["intact"], STATE_TO_INDEX["disrupted"]])
            stack.append(child)
    return expected, posterior_by_clade


def load_origin_counts(frame_dir: Path) -> tuple[str, str]:
    fitch_path = frame_dir / "origin_events.tsv"
    pastml_path = frame_dir / "pastml_origin_events.tsv"

    def count_rows(path: Path) -> str:
        if not path.exists():
            return ""
        with path.open(encoding="utf-8") as handle:
            return str(max(0, sum(1 for _ in handle) - 1))

    return count_rows(fitch_path), count_rows(pastml_path)


def build_row(
    *,
    scenario: str,
    analysis_frame: str,
    rooting_mode: str,
    tree_path: Path,
    tip_path: Path,
    notes: str,
    n_draws: int,
    seed: int,
) -> dict[str, Any]:
    if not tree_path.exists() or not tip_path.exists():
        return {
            "scenario": scenario,
            "analysis_frame": analysis_frame,
            "rooting_mode": rooting_mode,
            "status": "missing_tree_or_tip_states",
            "n_mk_draws": n_draws,
            "seed": seed,
            "source_tree": str(tree_path.relative_to(ROOT)) if tree_path.exists() else str(tree_path),
            "tip_state_file": str(tip_path.relative_to(ROOT)) if tip_path.exists() else str(tip_path),
            "notes": notes,
        }

    tree = Phylo.read(str(tree_path), "newick")
    tip_state_map, observed_intact, observed_disrupted, missing_tip_states = load_tip_states(tip_path)
    total_tips = len(list(tree.get_terminals()))
    rate, loglik = fit_equal_rates_mk(tree, tip_state_map)
    log_likelihoods, _ = compute_likelihoods(tree, tip_state_map, rate)
    expected_origins, _ = expected_origin_count(tree, log_likelihoods, rate)
    draws, root_disrupted_prob = sample_origin_distribution(
        tree,
        log_likelihoods,
        rate,
        n_draws=n_draws,
        seed=seed,
    )
    fitch_origin_events, pastml_origin_events = load_origin_counts(tree_path.parent)
    lower, upper = np.quantile(draws, [0.025, 0.975], method="nearest")

    return {
        "scenario": scenario,
        "analysis_frame": analysis_frame,
        "rooting_mode": rooting_mode,
        "status": "ok",
        "tip_count": total_tips,
        "observed_intact_tips": observed_intact,
        "observed_disrupted_tips": observed_disrupted,
        "missing_tip_states": missing_tip_states,
        "fitch_origin_events": fitch_origin_events,
        "pastml_origin_events": pastml_origin_events,
        "mk_equal_rates_log_likelihood": fmt(loglik, digits=4),
        "mk_equal_rates_q": f"{rate:.8g}",
        "root_disrupted_posterior_probability": fmt(root_disrupted_prob),
        "mk_expected_origin_events": fmt(expected_origins, digits=4),
        "mk_origin_count_mean": fmt(float(np.mean(draws)), digits=4),
        "mk_origin_count_median": int(np.median(draws)),
        "mk_origin_count_lower_95": int(lower),
        "mk_origin_count_upper_95": int(upper),
        "mk_origin_count_min": int(min(draws)),
        "mk_origin_count_max": int(max(draws)),
        "n_mk_draws": n_draws,
        "seed": seed,
        "source_tree": str(tree_path.relative_to(ROOT)),
        "tip_state_file": str(tip_path.relative_to(ROOT)),
        "notes": notes,
    }


def main() -> None:
    n_draws = 250
    seed = 20260412
    scenarios = [
        {
            "scenario": "composition_filtered_reference_rooted_primary",
            "analysis_frame": "composition_pruned_primary_quality_frame",
            "rooting_mode": "reference",
            "tree_path": ROOT / "outputs" / "workflow" / "asr_sensitivity" / "composition_filtered" / "rooted_ml_tree.reference_rooted.nwk",
            "tip_path": ROOT / "outputs" / "workflow" / "asr_sensitivity" / "composition_filtered" / "tip_states.tsv",
            "notes": "Binary equal-rates Mk summary on the composition-pruned primary ASR frame.",
        },
        {
            "scenario": "composition_filtered_midpoint_rooted",
            "analysis_frame": "composition_pruned_primary_quality_frame",
            "rooting_mode": "midpoint",
            "tree_path": ROOT / "outputs" / "workflow" / "asr_rooting_sensitivity" / "composition_filtered_midpoint" / "rooted_ml_tree.midpoint_rooted.nwk",
            "tip_path": ROOT / "outputs" / "workflow" / "asr_rooting_sensitivity" / "composition_filtered_midpoint" / "tip_states.tsv",
            "notes": "Binary equal-rates Mk summary on the composition-pruned midpoint rerun.",
        },
        {
            "scenario": "unpruned_reference_rooted_comparability",
            "analysis_frame": "unpruned_comparability_frame",
            "rooting_mode": "reference",
            "tree_path": ROOT / "outputs" / "workflow" / "asr" / "rooted_ml_tree.reference_rooted.nwk",
            "tip_path": ROOT / "outputs" / "workflow" / "asr" / "tip_states.tsv",
            "notes": "Binary equal-rates Mk summary on the original unpruned comparability tree.",
        },
        {
            "scenario": "unpruned_midpoint_rooted",
            "analysis_frame": "unpruned_comparability_frame",
            "rooting_mode": "midpoint",
            "tree_path": ROOT / "outputs" / "workflow" / "asr_rooting_sensitivity" / "unpruned_midpoint" / "rooted_ml_tree.midpoint_rooted.nwk",
            "tip_path": ROOT / "outputs" / "workflow" / "asr_rooting_sensitivity" / "unpruned_midpoint" / "tip_states.tsv",
            "notes": "Binary equal-rates Mk summary on the unpruned midpoint rerun.",
        },
    ]

    rows = [
        build_row(
            scenario=scenario["scenario"],
            analysis_frame=scenario["analysis_frame"],
            rooting_mode=scenario["rooting_mode"],
            tree_path=scenario["tree_path"],
            tip_path=scenario["tip_path"],
            notes=scenario["notes"],
            n_draws=n_draws,
            seed=seed,
        )
        for scenario in scenarios
    ]

    fieldnames = [
        "scenario",
        "analysis_frame",
        "rooting_mode",
        "status",
        "tip_count",
        "observed_intact_tips",
        "observed_disrupted_tips",
        "missing_tip_states",
        "fitch_origin_events",
        "pastml_origin_events",
        "mk_equal_rates_log_likelihood",
        "mk_equal_rates_q",
        "root_disrupted_posterior_probability",
        "mk_expected_origin_events",
        "mk_origin_count_mean",
        "mk_origin_count_median",
        "mk_origin_count_lower_95",
        "mk_origin_count_upper_95",
        "mk_origin_count_min",
        "mk_origin_count_max",
        "n_mk_draws",
        "seed",
        "source_tree",
        "tip_state_file",
        "notes",
    ]
    write_tsv(SUPP_DIR / "Supplementary_Table_20_Mk_Origin_Uncertainty.tsv", fieldnames, rows)
    write_tsv(FIGURE_DATA_DIR / "asr_mk_origin_uncertainty.tsv", fieldnames, rows)


if __name__ == "__main__":
    main()
