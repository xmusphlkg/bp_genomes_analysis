#!/usr/bin/env python3
"""Build manuscript-facing ASR, validation, and context ledgers.

This module bundles the submission diagnostics that sit alongside the main
manuscript outputs. It covers three functional groups:

1. ASR diagnostics and rooted-origin ledgers;
2. validation-evidence manifests and targeted follow-up audits; and
3. context diagnostics such as missingness, genotype background, and clock
   signal checks.

Rooting sensitivity is kept on the same Fitch + PastML footing by rerunning
midpoint scenarios through the full M5 wrapper instead of relying on
Fitch-only shortcuts.
"""

from __future__ import annotations

import csv
import math
import random
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path
from statistics import median
from typing import Any

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.stats import fisher_exact


ROOT = Path(__file__).resolve().parents[3]
SUPP_DIR = ROOT / "manuscript" / "supplementary"
FIGURE_DATA_DIR = ROOT / "manuscript" / "figure_data"
ASR_ROOTING_DIR = ROOT / "outputs" / "workflow" / "asr_rooting_sensitivity"
PRN_EVENT_EVIDENCE_MANIFEST = FIGURE_DATA_DIR / "prn_event_evidence_manifest.tsv"

sys.path.insert(0, str(ROOT / "workflow" / "lib"))
from project_paths import project_module_data_root  # noqa: E402
from asr_parsimony import Node, assign_node_ids, parse_newick  # noqa: E402

STEP4_OUTPUTS = project_module_data_root("step4_prn_validation") / "outputs"


STATE_ORDER = ["intact", "disrupted", "insufficient_data", "uncertain"]
LONGREAD_RE = re.compile(r"nanopore|pacbio|hifi|rsii|sequel|minion|promethion", re.I)


def clean_text(value: Any) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    text = str(value).strip()
    if text.casefold() in {"nan", "none", "na"}:
        return ""
    return text


def as_bool(value: Any) -> bool:
    return clean_text(value).casefold() in {"true", "1", "yes", "y"}


def as_float(value: Any) -> float:
    text = clean_text(value)
    if not text:
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def as_int(value: Any) -> int:
    number = as_float(value)
    if not np.isfinite(number):
        return 0
    return int(round(number))


def fmt(value: Any, digits: int = 4) -> str:
    number = as_float(value)
    if not np.isfinite(number):
        return ""
    return f"{number:.{digits}f}"


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


def copy_tsv(source: Path, destination: Path) -> None:
    rows = read_tsv(source)
    if not rows:
        raise ValueError(f"No rows in {source}")
    write_tsv(destination, list(rows[0]), rows)


def run_command(command: list[str]) -> None:
    subprocess.run(command, cwd=ROOT, check=True)


def count_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open(encoding="utf-8") as handle:
        return max(0, sum(1 for _ in handle) - 1)


def mode_text(values: Any) -> str:
    texts = [clean_text(value) for value in values if clean_text(value)]
    if not texts:
        return ""
    counts = Counter(texts)
    max_count = max(counts.values())
    return sorted(text for text, count in counts.items() if count == max_count)[0]


def summarize_asr_dir(outdir: Path) -> dict[str, Any]:
    tip_path = outdir / "tip_states.tsv"
    fitch_path = outdir / "origin_events.tsv"
    pastml_path = outdir / "pastml_origin_events.tsv"
    tip_rows = read_tsv(tip_path)
    pastml_rows = read_tsv(pastml_path) if pastml_path.exists() else []
    return {
        "tip_count": len(tip_rows),
        "disrupted_tip_count": sum(1 for row in tip_rows if row.get("prn_state") == "disrupted"),
        "fitch_origin_events": count_rows(fitch_path),
        "pastml_origin_events": len(pastml_rows) if pastml_path.exists() else "",
        "pastml_strict_origin_events": sum(1 for row in pastml_rows if row.get("origin_confidence") == "strict")
        if pastml_path.exists()
        else "",
        "pastml_compatible_origin_events": sum(
            1 for row in pastml_rows if row.get("origin_confidence") == "compatible"
        )
        if pastml_path.exists()
        else "",
        "pastml_status": "available" if pastml_path.exists() else "not_available",
    }


def run_full_rooting_scenario(
    *,
    tree: Path,
    outdir: Path,
    tree_id: str,
    rooting_mode: str,
    manifest: Path,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            "bash",
            str(ROOT / "workflow" / "bin" / "m5_asr.sh"),
            "--tree",
            str(tree),
            "--manifest",
            str(manifest),
            "--outdir",
            str(outdir),
            "--tree-id",
            tree_id,
            "--reference-label",
            "Reference",
            "--reference-state",
            "intact",
            "--rooting-mode",
            rooting_mode,
            "--pastml-threads",
            "2",
        ]
    )


def build_asr_rooting_sensitivity() -> None:
    manifest = ROOT / "outputs" / "workflow" / "manifest" / "manifest.tsv"
    scenarios = [
        {
            "scenario": "composition_filtered_reference_rooted_primary",
            "analysis_frame": "composition_pruned_quality_frame",
            "rooting_mode": "reference",
            "excluded_tip_count": 33,
            "source_tree": ROOT / "outputs" / "workflow" / "asr_sensitivity" / "composition_filtered.treefile",
            "outdir": ROOT / "outputs" / "workflow" / "asr_sensitivity" / "composition_filtered",
            "status": "existing_full_asr",
            "notes": "Primary ASR quality frame; 33 nonreference IQ-TREE composition-failed tips pruned.",
        },
        {
            "scenario": "composition_filtered_midpoint_rooted",
            "analysis_frame": "composition_pruned_quality_frame",
            "rooting_mode": "midpoint",
            "excluded_tip_count": 33,
            "source_tree": ROOT / "outputs" / "workflow" / "asr_sensitivity" / "composition_filtered.treefile",
            "outdir": ASR_ROOTING_DIR / "composition_filtered_midpoint",
            "status": "full_midpoint_rerun",
            "notes": "Midpoint-rooted full rerun under the same M5 Fitch + PastML wrapper used for the primary ASR frame.",
        },
        {
            "scenario": "unpruned_reference_rooted_comparability",
            "analysis_frame": "unpruned_comparability_frame",
            "rooting_mode": "reference",
            "excluded_tip_count": 0,
            "source_tree": ROOT / "outputs" / "workflow" / "phylo" / "iqtree2" / "ml_tree.treefile",
            "outdir": ROOT / "outputs" / "workflow" / "asr",
            "status": "existing_full_asr",
            "notes": "Original Reference/Tohama I rooted analysis retained as comparability sensitivity.",
        },
        {
            "scenario": "unpruned_midpoint_rooted",
            "analysis_frame": "unpruned_comparability_frame",
            "rooting_mode": "midpoint",
            "excluded_tip_count": 0,
            "source_tree": ROOT / "outputs" / "workflow" / "phylo" / "iqtree2" / "ml_tree.treefile",
            "outdir": ASR_ROOTING_DIR / "unpruned_midpoint",
            "status": "full_midpoint_rerun",
            "notes": "Midpoint-rooted full rerun under the same M5 Fitch + PastML wrapper used for the unpruned comparability frame.",
        },
    ]

    rows: list[dict[str, Any]] = []
    for scenario in scenarios:
        if scenario["status"] == "full_midpoint_rerun":
            run_full_rooting_scenario(
                tree=scenario["source_tree"],
                outdir=scenario["outdir"],
                tree_id=scenario["scenario"],
                rooting_mode=scenario["rooting_mode"],
                manifest=manifest,
            )
        summary = summarize_asr_dir(scenario["outdir"])
        rows.append(
            {
                "scenario": scenario["scenario"],
                "analysis_frame": scenario["analysis_frame"],
                "rooting_mode": scenario["rooting_mode"],
                "excluded_tip_count": scenario["excluded_tip_count"],
                "tip_count": summary["tip_count"],
                "disrupted_tip_count": summary["disrupted_tip_count"],
                "fitch_origin_events": summary["fitch_origin_events"],
                "pastml_origin_events": summary["pastml_origin_events"],
                "pastml_strict_origin_events": summary["pastml_strict_origin_events"],
                "pastml_compatible_origin_events": summary["pastml_compatible_origin_events"],
                "pastml_status": summary["pastml_status"],
                "source_tree": str(scenario["source_tree"].relative_to(ROOT)),
                "output_dir": str(scenario["outdir"].relative_to(ROOT)),
                "notes": scenario["notes"],
            }
        )

    mad_status = "mad_cli_available" if shutil.which("mad") else "mad_cli_not_available"
    rows.append(
        {
            "scenario": "mad_rooting_feasibility",
            "analysis_frame": "composition_pruned_quality_frame",
            "rooting_mode": "MAD",
            "excluded_tip_count": 33,
            "tip_count": "",
            "disrupted_tip_count": "",
            "fitch_origin_events": "",
            "pastml_origin_events": "",
            "pastml_strict_origin_events": "",
            "pastml_compatible_origin_events": "",
            "pastml_status": "not_run",
            "source_tree": "outputs/workflow/asr_sensitivity/composition_filtered.treefile",
            "output_dir": "",
            "notes": f"MAD rooting recorded as feasibility check; status={mad_status}.",
        }
    )

    fieldnames = [
        "scenario",
        "analysis_frame",
        "rooting_mode",
        "excluded_tip_count",
        "tip_count",
        "disrupted_tip_count",
        "fitch_origin_events",
        "pastml_origin_events",
        "pastml_strict_origin_events",
        "pastml_compatible_origin_events",
        "pastml_status",
        "source_tree",
        "output_dir",
        "notes",
    ]
    write_tsv(SUPP_DIR / "Supplementary_Table_13_ASR_Rooting_Sensitivity.tsv", fieldnames, rows)
    write_tsv(FIGURE_DATA_DIR / "asr_rooting_sensitivity.tsv", fieldnames, rows)


def build_reframed_asr_sensitivity_table() -> None:
    summary_rows = {row["scenario"]: row for row in read_tsv(ROOT / "outputs" / "workflow" / "asr_sensitivity" / "sensitivity_summary.tsv")}
    order = [
        (
            "composition_filtered",
            "composition_pruned_primary_quality_frame",
            "Primary ASR quality frame; prunes the 33 nonreference IQ-TREE composition-failed tips.",
        ),
        (
            "primary",
            "unpruned_reference_rooted_comparability",
            "Original Reference/Tohama I rooted unpruned tree retained as comparability sensitivity.",
        ),
        (
            "support_70",
            "unpruned_support_ge_70",
            "Unpruned comparability tree with Fitch origins restricted to branches with support >=70.",
        ),
        (
            "support_90",
            "unpruned_support_ge_90",
            "Unpruned comparability tree with Fitch origins restricted to branches with support >=90.",
        ),
    ]
    rows: list[dict[str, Any]] = []
    for source_key, scenario, notes in order:
        source = summary_rows[source_key]
        row = dict(source)
        row["scenario"] = scenario
        row["notes"] = notes
        rows.append(row)
    fieldnames = [
        "scenario",
        "excluded_tip_count",
        "tip_count",
        "disrupted_tip_count",
        "fitch_origin_events",
        "pastml_origin_events",
        "pastml_strict_origin_events",
        "pastml_compatible_origin_events",
        "notes",
    ]
    write_tsv(SUPP_DIR / "Supplementary_Table_6_ASR_Sensitivity.tsv", fieldnames, rows)
    write_tsv(FIGURE_DATA_DIR / "figure3_workflow_asr_sensitivity.tsv", fieldnames, rows)


def stage_primary_origin_events() -> None:
    source = ROOT / "outputs" / "workflow" / "asr_sensitivity" / "composition_filtered" / "origin_events.tsv"
    copy_tsv(source, SUPP_DIR / "Supplementary_Table_3_independent_origins.tsv")
    copy_tsv(source, FIGURE_DATA_DIR / "figure3_workflow_origin_events.tsv")


def fitch_downpass_sets(node: Node, tip_state_map: dict[str, str]) -> set[str]:
    if node.is_tip:
        node.candidate_states = {tip_state_map.get(node.tree_label, "insufficient_data")}
        return node.candidate_states
    child_sets = [fitch_downpass_sets(child, tip_state_map) for child in node.children]
    intersection = set.intersection(*child_sets)
    node.candidate_states = intersection if intersection else set.union(*child_sets)
    return node.candidate_states


def choose_state(states: set[str], preferred_order: list[str], rng: random.Random | None = None) -> str:
    if rng is not None:
        return rng.choice(sorted(states))
    for state in preferred_order:
        if state in states:
            return state
    return sorted(states)[0]


def fitch_uppass_tiebreak(
    node: Node,
    preferred_order: list[str],
    *,
    parent_state: str | None = None,
    rng: random.Random | None = None,
) -> None:
    if parent_state and parent_state in node.candidate_states:
        node.inferred_state = parent_state
    else:
        node.inferred_state = choose_state(node.candidate_states, preferred_order, rng)
    for child in node.children:
        fitch_uppass_tiebreak(child, preferred_order, parent_state=node.inferred_state, rng=rng)


def count_fitch_origins_with_tiebreak(
    tree_text: str,
    tip_state_map: dict[str, str],
    preferred_order: list[str],
    *,
    rng: random.Random | None = None,
) -> int:
    root = parse_newick(tree_text)
    nodes = assign_node_ids(root)
    fitch_downpass_sets(root, tip_state_map)
    fitch_uppass_tiebreak(root, preferred_order, rng=rng)
    return sum(
        1
        for node in nodes
        if node.parent is not None and node.parent.inferred_state == "intact" and node.inferred_state == "disrupted"
    )


def build_tiebreak_sensitivity() -> None:
    tree_path = ROOT / "outputs" / "workflow" / "asr_sensitivity" / "composition_filtered" / "rooted_ml_tree.reference_rooted.nwk"
    tip_path = ROOT / "outputs" / "workflow" / "asr_sensitivity" / "composition_filtered" / "tip_states.tsv"
    tree_text = tree_path.read_text(encoding="utf-8")
    tip_state_map = {row["tree_tip_label"]: row["prn_state"] for row in read_tsv(tip_path)}
    deterministic = {
        "intact_biased_primary": ["intact", "disrupted", "insufficient_data", "uncertain"],
        "disrupted_biased_pessimistic": ["disrupted", "intact", "insufficient_data", "uncertain"],
        "insufficient_biased_diagnostic": ["insufficient_data", "intact", "disrupted", "uncertain"],
    }
    rows: list[dict[str, Any]] = []
    for scenario, order in deterministic.items():
        count = count_fitch_origins_with_tiebreak(tree_text, tip_state_map, order)
        rows.append(
            {
                "tree_frame": "composition_pruned_reference_rooted_primary",
                "scenario": scenario,
                "iterations": 1,
                "fitch_origin_events": count,
                "origin_count_min": count,
                "origin_count_median": count,
                "origin_count_max": count,
                "tip_count": len(tip_state_map),
                "disrupted_tip_count": sum(1 for state in tip_state_map.values() if state == "disrupted"),
                "notes": "Deterministic Fitch uppass tie rule.",
            }
        )

    random_counts = [
        count_fitch_origins_with_tiebreak(
            tree_text,
            tip_state_map,
            STATE_ORDER,
            rng=random.Random(42 + index),
        )
        for index in range(1000)
    ]
    rows.append(
        {
            "tree_frame": "composition_pruned_reference_rooted_primary",
            "scenario": "random_unbiased_tie_break_seed42_1000",
            "iterations": 1000,
            "fitch_origin_events": fmt(float(np.mean(random_counts)), 2),
            "origin_count_min": min(random_counts),
            "origin_count_median": median(random_counts),
            "origin_count_max": max(random_counts),
            "tip_count": len(tip_state_map),
            "disrupted_tip_count": sum(1 for state in tip_state_map.values() if state == "disrupted"),
            "notes": "Randomly selects among candidate states when parent state is not admissible; fixed seed series 42..1041.",
        }
    )
    fieldnames = [
        "tree_frame",
        "scenario",
        "iterations",
        "fitch_origin_events",
        "origin_count_min",
        "origin_count_median",
        "origin_count_max",
        "tip_count",
        "disrupted_tip_count",
        "notes",
    ]
    write_tsv(SUPP_DIR / "Supplementary_Table_14_ASR_TieBreak_Sensitivity.tsv", fieldnames, rows)


def build_composition_failed_tip_audit() -> None:
    failed_path = ROOT / "outputs" / "workflow" / "asr_sensitivity" / "composition_failed_nonreference.txt"
    failed = [line.strip() for line in failed_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    manifest = pd.read_csv(ROOT / "outputs" / "workflow" / "manifest" / "manifest.tsv", sep="\t", dtype=str)
    tip_states = pd.read_csv(ROOT / "outputs" / "workflow" / "asr" / "tip_states.tsv", sep="\t", dtype=str)
    tip_states = tip_states[
        ["tree_tip_label", "prn_state", "observed_prn_mechanism_call", "prn_call_confidence"]
    ].rename(columns={"tree_tip_label": "assembly_accession", "prn_state": "asr_tip_state"})
    merged = pd.DataFrame({"assembly_accession": failed}).merge(manifest, on="assembly_accession", how="left")
    merged = merged.merge(
        tip_states[["assembly_accession", "asr_tip_state", "observed_prn_mechanism_call", "prn_call_confidence"]],
        on="assembly_accession",
        how="left",
    )
    rows = []
    for row in merged.fillna("").to_dict("records"):
        rows.append(
            {
                "assembly_accession": row.get("assembly_accession", ""),
                "sample_id_canonical": row.get("sample_id_canonical", ""),
                "country_iso3": row.get("country_iso3", ""),
                "year": row.get("year", ""),
                "asr_tip_state": row.get("asr_tip_state", ""),
                "prn_mechanism_call": row.get("prn_mechanism_call") or row.get("observed_prn_mechanism_call", ""),
                "prn_call_confidence": row.get("prn_call_confidence", ""),
                "sequencing_tech": row.get("sequencing_tech", ""),
                "raw_read_link_status": row.get("raw_read_link_status", ""),
                "included_in_unpruned_asr": "True" if row.get("asr_tip_state") else "False",
                "included_in_composition_pruned_primary_asr": "False",
                "notes": "Nonreference IQ-TREE composition chi-square failure pruned from the primary composition-quality ASR frame.",
            }
        )
    fieldnames = [
        "assembly_accession",
        "sample_id_canonical",
        "country_iso3",
        "year",
        "asr_tip_state",
        "prn_mechanism_call",
        "prn_call_confidence",
        "sequencing_tech",
        "raw_read_link_status",
        "included_in_unpruned_asr",
        "included_in_composition_pruned_primary_asr",
        "notes",
    ]
    write_tsv(SUPP_DIR / "composition_failed_tip_audit.tsv", fieldnames, rows)
    write_tsv(SUPP_DIR / "Supplementary_Table_15_Composition_Failed_Tip_Audit.tsv", fieldnames, rows)


def aggregate_read_validation_by_event() -> dict[str, dict[str, Any]]:
    path = STEP4_OUTPUTS / "bp_prn_read_validation.tsv"
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    out: dict[str, dict[str, Any]] = {}
    for event_id, group in df.groupby("prn_event_id", dropna=False):
        if not clean_text(event_id):
            continue
        supporting = [as_int(value) for value in group.get("n_supporting_reads", pd.Series(dtype=str))]
        out[event_id] = {
            "supporting_read_count": max(supporting) if supporting else 0,
            "supporting_validation_rows": int(
                group["read_validation_status"].isin(["supported", "supported_candidate", "supported_concordant"]).sum()
            )
            if "read_validation_status" in group
            else 0,
            "read_validation_statuses": ";".join(sorted({clean_text(value) for value in group["read_validation_status"] if clean_text(value)})),
            "read_support_classes": ";".join(sorted({clean_text(value) for value in group["read_support_class"] if clean_text(value)})),
        }
    return out


def aggregate_tsd_by_event() -> dict[str, dict[str, Any]]:
    path = STEP4_OUTPUTS / "bp_prn_read_validation_tsd.tsv"
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    out: dict[str, dict[str, Any]] = {}
    for event_id, group in df.groupby("prn_event_id", dropna=False):
        if not clean_text(event_id):
            continue
        clipped = [as_int(value) for value in group.get("total_clipped_reads", pd.Series(dtype=str))]
        out[event_id] = {
            "tsd_supported_validation_rows": len(group),
            "representative_tsd_direct_repeats": mode_text(group.get("direct_repeats", pd.Series(dtype=str))),
            "representative_left_flank_sequence": mode_text(group.get("left_sequence", pd.Series(dtype=str))),
            "representative_right_flank_sequence": mode_text(group.get("right_sequence", pd.Series(dtype=str))),
            "representative_inverted_repeat_signature": mode_text(
                group.get("inverted_repeats", pd.Series(dtype=str))
            ),
            "max_total_clipped_reads": max(clipped) if clipped else 0,
        }
    return out


def derive_event_evidence_type(row: dict[str, Any]) -> str:
    validation_level = clean_text(row.get("validation_level"))
    supporting_reads = as_int(row.get("supporting_read_count"))
    if supporting_reads > 0:
        return "read_backed_targeted_validation"
    if "longread" in validation_level or "hybrid" in validation_level or clean_text(row.get("longread_exemplar")):
        return "public_longread_or_hybrid_assembly"
    if validation_level == "assembly_only":
        return "assembly_only"
    if validation_level:
        return validation_level
    return "not_recovered_current_public_data"


def tsd_or_flank_status(row: dict[str, Any]) -> str:
    if clean_text(row.get("tsd_direct_repeats")):
        return "target_site_duplication_recovered"
    if clean_text(row.get("read_locus_start")) or clean_text(row.get("read_locus_end")):
        return "read_interval_without_tsd"
    if clean_text(row.get("example_gap_start")) or clean_text(row.get("example_gap_end")):
        return "assembly_coordinate_only"
    return "not_recovered_current_public_data"


def extend_event_evidence_manifest() -> None:
    path = SUPP_DIR / "Supplementary_Table_9_prn_Event_Definitions.tsv"
    rows = read_tsv(path)
    read_by_event = aggregate_read_validation_by_event()
    tsd = pd.read_csv(STEP4_OUTPUTS / "bp_prn_read_validation_tsd.tsv", sep="\t", dtype=str).fillna("")
    tsd_by_event = {
        event_id: group.iloc[0].to_dict() for event_id, group in tsd.groupby("prn_event_id", dropna=False) if clean_text(event_id)
    }
    followup = pd.read_csv(
        STEP4_OUTPUTS / "bp_prn_targeted_validation_followup_queue.tsv",
        sep="\t",
        dtype=str,
    ).fillna("")
    origins_by_event = {
        event_id: ";".join(sorted({clean_text(value) for value in group["origin_id"] if clean_text(value)}))
        for event_id, group in followup.groupby("prn_event_id", dropna=False)
        if clean_text(event_id)
    }

    extended: list[dict[str, Any]] = []
    for row in rows:
        event_id = row["prn_event_id"]
        read_summary = read_by_event.get(event_id, {})
        tsd_row = tsd_by_event.get(event_id, {})
        longread_exemplar = row.get("supporting_read_or_public_longread", "")
        next_row = dict(row)
        next_row["rooting_frame"] = "event_catalog_not_rooting_specific;primary_asr=composition_pruned"
        next_row["breakpoint_left"] = row.get("read_locus_start") or row.get("example_gap_start", "")
        next_row["breakpoint_right"] = row.get("read_locus_end") or row.get("example_gap_end", "")
        next_row["breakpoint_coordinate_basis"] = "read_reference" if row.get("read_locus_start") else "example_assembly"
        next_row["orientation"] = row.get("hit_orientation", "")
        next_row["supporting_read_count"] = read_summary.get("supporting_read_count", "")
        next_row["supporting_validation_rows"] = read_summary.get("supporting_validation_rows", "")
        next_row["read_validation_statuses"] = read_summary.get("read_validation_statuses", "")
        next_row["read_support_classes"] = read_summary.get("read_support_classes", "")
        if not next_row.get("tsd_direct_repeats"):
            next_row["tsd_direct_repeats"] = tsd_row.get("direct_repeats", "")
        next_row["tsd_or_flank_sequence_status"] = tsd_or_flank_status(next_row)
        next_row["longread_exemplar"] = longread_exemplar
        next_row["validation_status"] = read_summary.get("read_validation_statuses") or row.get("validation_level", "")
        next_row["evidence_type"] = derive_event_evidence_type(next_row)
        next_row["priority_origin_ids"] = origins_by_event.get(event_id, "")
        next_row["validation_priority"] = (
            "dominant_1043bp_architecture"
            if event_id == "prn_evt_coding_disrupted_is481__is481__gap1043"
            else ("origin_linked_event" if next_row["priority_origin_ids"] else "")
        )
        extended.append(next_row)

    fieldnames = list(rows[0].keys())
    for column in [
        "rooting_frame",
        "evidence_type",
        "breakpoint_left",
        "breakpoint_right",
        "breakpoint_coordinate_basis",
        "orientation",
        "tsd_or_flank_sequence_status",
        "supporting_read_count",
        "supporting_validation_rows",
        "read_validation_statuses",
        "read_support_classes",
        "longread_exemplar",
        "validation_status",
        "priority_origin_ids",
        "validation_priority",
    ]:
        if column not in fieldnames:
            fieldnames.append(column)
    write_tsv(path, fieldnames, extended)
    write_tsv(PRN_EVENT_EVIDENCE_MANIFEST, fieldnames, extended)


def build_is481_event_evidence_audit() -> None:
    manifest_path = PRN_EVENT_EVIDENCE_MANIFEST
    if not manifest_path.exists():
        return
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    tsd_by_event = aggregate_tsd_by_event()
    mask = (
        manifest.get("is_element_name", pd.Series(dtype=str)).eq("IS481")
        | manifest.get("event_definition_rule", pd.Series(dtype=str)).str.contains("IS481", case=False, na=False)
    )
    subset = manifest.loc[mask].copy()
    if subset.empty:
        return

    subset["sample_count_numeric"] = subset.get("sample_count", "").map(as_int)
    subset["gap_numeric"] = subset.get("insertion_subject_gap_bp", "").map(as_int)
    subset = subset.sort_values(
        by=["sample_count_numeric", "gap_numeric", "prn_event_id"],
        ascending=[False, False, True],
        kind="stable",
    )

    rows: list[dict[str, Any]] = []
    for row in subset.to_dict("records"):
        event_id = clean_text(row.get("prn_event_id"))
        tsd_row = tsd_by_event.get(event_id, {})
        representative_tsd = clean_text(row.get("tsd_direct_repeats")) or clean_text(
            tsd_row.get("representative_tsd_direct_repeats")
        )
        notes: list[str] = []
        tsd_status = clean_text(row.get("tsd_or_flank_sequence_status"))
        if representative_tsd:
            notes.append("read_level_tsd_recovered")
        elif tsd_status == "read_interval_without_tsd":
            notes.append("read_interval_recovered_without_tsd")
        elif "longread" in clean_text(row.get("evidence_type")) or "hybrid" in clean_text(row.get("evidence_type")):
            notes.append("public_longread_or_hybrid_anchor_without_read_junction")
        elif clean_text(row.get("evidence_type")) == "assembly_only":
            notes.append("assembly_only_without_public_read_junction")
        else:
            notes.append("targeted_validation_pending_or_unresolved")
        if event_id == "prn_evt_coding_disrupted_is481__is481__gap1043":
            notes.append("dominant_1043bp_architecture")
        if clean_text(row.get("validation_priority")):
            notes.append(clean_text(row.get("validation_priority")))

        rows.append(
            {
                "prn_event_id": event_id,
                "event_subcategory": row.get("event_subcategory", ""),
                "sample_count": row.get("sample_count", ""),
                "country_count": row.get("country_count", ""),
                "year_min": row.get("year_min", ""),
                "year_max": row.get("year_max", ""),
                "insertion_subject_gap_bp": row.get("insertion_subject_gap_bp", ""),
                "orientation": row.get("orientation", ""),
                "evidence_type": row.get("evidence_type", ""),
                "validation_level": row.get("validation_level", ""),
                "breakpoint_coordinate_basis": row.get("breakpoint_coordinate_basis", ""),
                "breakpoint_left": row.get("breakpoint_left", ""),
                "breakpoint_right": row.get("breakpoint_right", ""),
                "tsd_or_flank_sequence_status": tsd_status,
                "representative_tsd_direct_repeats": representative_tsd,
                "representative_left_flank_sequence": clean_text(
                    tsd_row.get("representative_left_flank_sequence")
                ),
                "representative_right_flank_sequence": clean_text(
                    tsd_row.get("representative_right_flank_sequence")
                ),
                "representative_inverted_repeat_signature": clean_text(
                    tsd_row.get("representative_inverted_repeat_signature")
                ),
                "max_total_clipped_reads": tsd_row.get("max_total_clipped_reads", ""),
                "tsd_supported_validation_rows": tsd_row.get("tsd_supported_validation_rows", ""),
                "supporting_read_count": row.get("supporting_read_count", ""),
                "supporting_validation_rows": row.get("supporting_validation_rows", ""),
                "read_validation_statuses": row.get("read_validation_statuses", ""),
                "longread_exemplar": row.get("longread_exemplar", ""),
                "supporting_read_or_public_longread": row.get("supporting_read_or_public_longread", ""),
                "priority_origin_ids": row.get("priority_origin_ids", ""),
                "validation_priority": row.get("validation_priority", ""),
                "notes": ";".join(dict.fromkeys(notes)),
            }
        )

    fieldnames = [
        "prn_event_id",
        "event_subcategory",
        "sample_count",
        "country_count",
        "year_min",
        "year_max",
        "insertion_subject_gap_bp",
        "orientation",
        "evidence_type",
        "validation_level",
        "breakpoint_coordinate_basis",
        "breakpoint_left",
        "breakpoint_right",
        "tsd_or_flank_sequence_status",
        "representative_tsd_direct_repeats",
        "representative_left_flank_sequence",
        "representative_right_flank_sequence",
        "representative_inverted_repeat_signature",
        "max_total_clipped_reads",
        "tsd_supported_validation_rows",
        "supporting_read_count",
        "supporting_validation_rows",
        "read_validation_statuses",
        "longread_exemplar",
        "supporting_read_or_public_longread",
        "priority_origin_ids",
        "validation_priority",
        "notes",
    ]
    write_tsv(SUPP_DIR / "Supplementary_Table_21_IS481_Event_Evidence_Audit.tsv", fieldnames, rows)
    write_tsv(FIGURE_DATA_DIR / "is481_event_evidence_audit.tsv", fieldnames, rows)


def enrich_validation_followup_queue() -> None:
    path = SUPP_DIR / "Supplementary_Table_12_Targeted_Validation_Followup.tsv"
    rows = read_tsv(path)
    enriched: list[dict[str, Any]] = []
    for row in rows:
        followup_class = clean_text(row.get("followup_class"))
        recovery_status = clean_text(row.get("recovery_plan_status"))
        validation_status = clean_text(row.get("read_validation_status"))
        raw_reads = as_bool(row.get("raw_reads_available"))
        read_accession = clean_text(row.get("read_accession_primary"))
        origin_id = clean_text(row.get("origin_id"))
        mechanism = clean_text(row.get("prn_mechanism_call"))
        if validation_status in {"supported", "supported_concordant"}:
            recovery = "targeted_read_validation_completed_read_backed"
        elif validation_status == "supported_candidate":
            recovery = "targeted_read_validation_completed_candidate_signal"
        elif validation_status == "no_prn_is_signal_detected":
            recovery = "targeted_read_validation_completed_no_local_signal"
        elif validation_status == "tool_output_missing":
            recovery = "targeted_read_validation_attempted_tool_output_missing"
        elif followup_class == "public_longread_or_hybrid_exemplar_present":
            recovery = "public_longread_or_hybrid_anchor_present"
        elif recovery_status == "recoverable_paired_illumina":
            recovery = "public_reads_recoverable_via_paired_illumina_fastq"
        elif recovery_status == "linked_incompatible_run_current_short_read_validator":
            recovery = "public_run_present_but_requires_non_illumina_or_nonpaired_validator"
        elif recovery_status == "linked_run_without_fastq_ftp":
            recovery = "public_run_metadata_present_fastq_not_indexed"
        elif raw_reads and read_accession:
            recovery = "public_reads_linked_but_recovery_mode_not_curated"
        elif origin_id == "origin_0006":
            recovery = "not_recoverable_from_current_public_data;highest_priority_for_new_assay_or_read_recovery"
        else:
            recovery = "not_recoverable_from_current_public_data"
        priority = []
        if as_bool(row.get("is_priority_origin_exemplar")):
            priority.append("origin_exemplar")
        if origin_id == "origin_0006":
            priority.append("assembly_only_origin_0006")
        if "inversion" in mechanism or "rearrangement" in mechanism:
            priority.append("rearrangement_breakpoint")
        if as_bool(row.get("is_unresolved_validation_row")):
            priority.append("unresolved_read_validation")
        if recovery_status == "linked_incompatible_run_current_short_read_validator":
            priority.append("incompatible_validator_input")
        if recovery_status == "linked_run_without_fastq_ftp":
            priority.append("no_fastq_ftp")
        next_row = dict(row)
        next_row["public_data_recovery_status"] = recovery
        next_row["validation_priority"] = ";".join(priority)
        next_row["explicit_unresolved_note"] = (
            "Unresolved status is not treated as negative evidence; current public-data action is recorded separately."
            if "unresolved" in priority
            else ""
        )
        enriched.append(next_row)
    fieldnames = list(rows[0].keys())
    for column in ["public_data_recovery_status", "validation_priority", "explicit_unresolved_note"]:
        if column not in fieldnames:
            fieldnames.append(column)
    write_tsv(path, fieldnames, enriched)
    write_tsv(FIGURE_DATA_DIR / "figure6_targeted_validation_followup.tsv", fieldnames, enriched)


def origin_evidence_alignment_note(
    representative_level: str,
    event_level: str,
    followup_class: str,
) -> str:
    representative = clean_text(representative_level)
    event = clean_text(event_level)
    followup = clean_text(followup_class)
    if representative == event and representative:
        return "origin_exemplar_matches_dominant_event_evidence_tier"
    if representative == "assembly_only" and event:
        return "dominant_event_has_external_support_but_origin_exemplar_remains_assembly_only"
    if representative == "public_longread_or_hybrid_assembly" and event.startswith("read_backed"):
        return "origin_exemplar_is_longread_anchored_while_same_event_is_read_backed_elsewhere"
    if representative.startswith("read_backed") and event == "public_longread_or_hybrid_assembly":
        return "origin_exemplar_is_read_backed_while_event_summary_is_longread_anchored"
    if followup == "assembly_only":
        return "assembly_only_origin_exemplar_requires_new_data_or_same_event_anchor"
    if followup == "public_longread_or_hybrid_exemplar_present":
        return "origin_exemplar_is_public_longread_or_hybrid_anchor"
    if followup == "read_backed_or_candidate_available":
        return "origin_exemplar_has_read_backed_or_candidate_support"
    return "origin_exemplar_requires_manual_audit"


def followup_row_has_longread_signal(row: dict[str, str]) -> bool:
    return bool(LONGREAD_RE.search(clean_text(row.get("sequencing_tech"))))


def followup_validation_level(row: dict[str, str]) -> str:
    status = clean_text(row.get("read_validation_status"))
    if status in {"supported", "supported_concordant"}:
        return "read_backed_supported"
    if followup_row_has_longread_signal(row):
        return "public_longread_or_hybrid_assembly"
    if status == "supported_candidate":
        return "read_backed_candidate"
    if status == "no_prn_is_signal_detected":
        return "read_backed_no_local_signal"
    if status == "unresolved":
        return "read_validation_unresolved"
    return "assembly_only"


def validation_level_rank(level: str) -> int:
    ordering = {
        "read_backed_supported": 0,
        "public_longread_or_hybrid_assembly": 1,
        "read_backed_candidate": 2,
        "read_backed_no_local_signal": 3,
        "read_validation_unresolved": 4,
        "assembly_only": 5,
    }
    return ordering.get(clean_text(level), 99)


def followup_supporting_hook(row: dict[str, str], validation_level: str) -> str:
    if validation_level.startswith("read_backed"):
        return (
            clean_text(row.get("read_accession_primary"))
            or clean_text(row.get("sra_run_accession"))
            or clean_text(row.get("ena_run_accession"))
            or clean_text(row.get("sample_id_canonical"))
        )
    if validation_level == "public_longread_or_hybrid_assembly":
        assembly = clean_text(row.get("assembly_accession"))
        sequencing_tech = clean_text(row.get("sequencing_tech"))
        if assembly and sequencing_tech:
            return f"{assembly}::{sequencing_tech}"
        return assembly or sequencing_tech
    return clean_text(row.get("assembly_accession"))


def choose_package_validation_exemplar(origin_rows: list[dict[str, str]]) -> dict[str, str]:
    if not origin_rows:
        return {}

    priority_rows = [row for row in origin_rows if as_bool(row.get("is_priority_origin_exemplar"))]
    candidates = priority_rows or origin_rows

    def rank(row: dict[str, str]) -> tuple[Any, ...]:
        validation_level = followup_validation_level(row)
        if validation_level in {
            "read_backed_supported",
            "public_longread_or_hybrid_assembly",
            "read_backed_candidate",
        }:
            evidence_rank = validation_level_rank(validation_level)
        elif clean_text(row.get("recovery_plan_status")) == "recoverable_paired_illumina":
            evidence_rank = 3
        elif clean_text(row.get("prn_call_confidence")) in {"assembly_high", "assembly_moderate"}:
            evidence_rank = 4
        else:
            evidence_rank = 5
        year = as_int(row.get("year")) if clean_text(row.get("year")) else 9999
        return (
            evidence_rank,
            -as_int(row.get("origin_n_disrupted_descendants")),
            year,
            clean_text(row.get("sample_id_canonical")),
        )

    best = dict(sorted(candidates, key=rank)[0])
    validation_level = followup_validation_level(best)
    best["selected_validation_level"] = validation_level
    best["selected_supporting_read_or_public_longread"] = followup_supporting_hook(best, validation_level)
    return best


def build_origin_evidence_completeness_audit() -> None:
    origin_rows = read_tsv(FIGURE_DATA_DIR / "fig03_independent_origins.tsv")
    followup_rows = read_tsv(SUPP_DIR / "Supplementary_Table_12_Targeted_Validation_Followup.tsv")
    event_rows = read_tsv(PRN_EVENT_EVIDENCE_MANIFEST)

    followup_by_origin: dict[str, list[dict[str, str]]] = {}
    for row in followup_rows:
        origin_id = clean_text(row.get("origin_id"))
        if not origin_id or not as_bool(row.get("is_origin_defining_tip")):
            continue
        followup_by_origin.setdefault(origin_id, []).append(row)
    event_by_id = {clean_text(row.get("prn_event_id")): row for row in event_rows}
    exemplar_selection_rule = (
        "within_origin_package_origin_defining_tip_ranked_as_"
        "read_backed_supported>public_longread_or_hybrid_assembly>read_backed_candidate>"
        "recoverable_paired_illumina>assembly_high_or_moderate>other"
    )

    rows: list[dict[str, Any]] = []
    for origin_row in origin_rows:
        origin_id = clean_text(origin_row.get("origin_id"))
        event_id = clean_text(origin_row.get("dominant_prn_event_id"))
        tree_representative_sample = clean_text(origin_row.get("representative_sample_id_canonical"))
        followup_row = choose_package_validation_exemplar(followup_by_origin.get(origin_id, []))
        event_row = event_by_id.get(event_id, {})
        tree_representative_level = clean_text(origin_row.get("validation_level"))
        representative_level = clean_text(followup_row.get("selected_validation_level")) or tree_representative_level
        representative_support = (
            clean_text(followup_row.get("selected_supporting_read_or_public_longread"))
            or clean_text(origin_row.get("supporting_read_or_public_longread"))
        )
        exemplar_replacement_applied = (
            "True"
            if clean_text(followup_row.get("sample_id_canonical"))
            and clean_text(followup_row.get("sample_id_canonical")) != tree_representative_sample
            else "False"
        )
        event_level = clean_text(event_row.get("validation_level"))
        rows.append(
            {
                "origin_id": origin_id,
                "phylo_tree_id": clean_text(origin_row.get("phylo_tree_id")),
                "dominant_prn_mechanism": clean_text(origin_row.get("dominant_prn_mechanism")),
                "dominant_prn_event_id": event_id,
                "origin_n_disrupted_tips": clean_text(origin_row.get("n_tips_disrupted")),
                "tree_representative_sample_id_canonical": tree_representative_sample,
                "tree_representative_assembly_accession": clean_text(
                    origin_row.get("representative_assembly_accession")
                ),
                "tree_representative_country_iso3": clean_text(origin_row.get("representative_country_iso3")),
                "tree_representative_year": clean_text(origin_row.get("representative_year")),
                "tree_representative_validation_level": tree_representative_level,
                "tree_representative_supporting_read_or_public_longread": clean_text(
                    origin_row.get("supporting_read_or_public_longread")
                ),
                "representative_sample_id_canonical": clean_text(
                    followup_row.get("sample_id_canonical")
                )
                or tree_representative_sample,
                "representative_assembly_accession": clean_text(
                    followup_row.get("assembly_accession")
                )
                or clean_text(origin_row.get("representative_assembly_accession")),
                "representative_country_iso3": clean_text(followup_row.get("country_iso3"))
                or clean_text(origin_row.get("representative_country_iso3")),
                "representative_year": clean_text(followup_row.get("year"))
                or clean_text(origin_row.get("representative_year")),
                "representative_validation_level": representative_level,
                "representative_supporting_read_or_public_longread": representative_support,
                "exemplar_replacement_applied": exemplar_replacement_applied,
                "exemplar_selection_rule": exemplar_selection_rule,
                "dominant_event_validation_level": event_level,
                "dominant_event_supporting_validation_rows": clean_text(
                    event_row.get("supporting_validation_rows")
                ),
                "dominant_event_supporting_read_or_public_longread": clean_text(
                    event_row.get("supporting_read_or_public_longread")
                ),
                "dominant_event_longread_exemplar": clean_text(event_row.get("longread_exemplar")),
                "followup_class": clean_text(followup_row.get("followup_class")),
                "public_data_recovery_status": clean_text(
                    followup_row.get("public_data_recovery_status")
                ),
                "recovery_plan_status": clean_text(followup_row.get("recovery_plan_status")),
                "read_validation_status": clean_text(followup_row.get("read_validation_status")),
                "read_support_class": clean_text(followup_row.get("read_support_class")),
                "sequencing_tech": clean_text(followup_row.get("sequencing_tech")),
                "evidence_alignment": origin_evidence_alignment_note(
                    representative_level,
                    event_level,
                    clean_text(followup_row.get("followup_class")),
                ),
                "validation_priority": clean_text(
                    followup_row.get("validation_priority")
                ),
                "notes": clean_text(origin_row.get("notes")),
            }
        )

    fieldnames = [
        "origin_id",
        "phylo_tree_id",
        "dominant_prn_mechanism",
        "dominant_prn_event_id",
        "origin_n_disrupted_tips",
        "tree_representative_sample_id_canonical",
        "tree_representative_assembly_accession",
        "tree_representative_country_iso3",
        "tree_representative_year",
        "tree_representative_validation_level",
        "tree_representative_supporting_read_or_public_longread",
        "representative_sample_id_canonical",
        "representative_assembly_accession",
        "representative_country_iso3",
        "representative_year",
        "representative_validation_level",
        "representative_supporting_read_or_public_longread",
        "exemplar_replacement_applied",
        "exemplar_selection_rule",
        "dominant_event_validation_level",
        "dominant_event_supporting_validation_rows",
        "dominant_event_supporting_read_or_public_longread",
        "dominant_event_longread_exemplar",
        "followup_class",
        "public_data_recovery_status",
        "recovery_plan_status",
        "read_validation_status",
        "read_support_class",
        "sequencing_tech",
        "evidence_alignment",
        "validation_priority",
        "notes",
    ]
    write_tsv(SUPP_DIR / "Supplementary_Table_24_Origin_Evidence_Completeness.tsv", fieldnames, rows)
    write_tsv(FIGURE_DATA_DIR / "origin_evidence_completeness_audit.tsv", fieldnames, rows)


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    keep = values.notna() & weights.notna() & (weights > 0)
    if not keep.any():
        return math.nan
    return float(np.average(values.loc[keep], weights=weights.loc[keep]))


def build_ipw_diagnostics() -> None:
    manifest = pd.read_csv(ROOT / "outputs" / "workflow" / "manifest" / "manifest.tsv", sep="\t", dtype=str)
    predictions = pd.read_csv(
        ROOT / "outputs" / "workflow" / "missingness_model" / "missingness_model_predictions.tsv",
        sep="\t",
        dtype=str,
    )
    ipw = pd.read_csv(ROOT / "outputs" / "workflow" / "epi" / "ipw_prevalence.tsv", sep="\t", dtype=str)
    rows: list[dict[str, Any]] = []

    perf = predictions.copy()
    perf["prob_interpretable"] = pd.to_numeric(perf["prob_interpretable"], errors="coerce")
    perf["y_actual"] = pd.to_numeric(perf["y_actual"], errors="coerce")
    perf["in_model_flag"] = perf.get("in_model", pd.Series(index=perf.index, dtype=str)).map(as_bool)
    oof_column = next(
        (
            column
            for column in ["prob_interpretable_oof", "oof_prob_interpretable", "prob_interpretable_out_of_fold"]
            if column in perf.columns
        ),
        "",
    )
    metric_scope = "training_only_predictions"
    metric_probability_column = "prob_interpretable"
    if oof_column:
        perf[oof_column] = pd.to_numeric(perf[oof_column], errors="coerce")
        oof_perf = perf.dropna(subset=[oof_column, "y_actual"]).copy()
        if not oof_perf.empty:
            perf = oof_perf
            metric_scope = "out_of_fold_predictions"
            metric_probability_column = oof_column
    else:
        perf = perf.loc[perf["in_model_flag"]].dropna(subset=["prob_interpretable", "y_actual"]).copy()

    if not perf.empty:
        brier = float(np.mean((perf[metric_probability_column] - perf["y_actual"]) ** 2))
        pred = (perf[metric_probability_column] >= 0.5).astype(int)
        accuracy = float((pred == perf["y_actual"]).mean())
        rows.append(
            {
                "diagnostic_scope": "missingness_model_performance",
                "metric": "out_of_fold_accuracy" if metric_scope == "out_of_fold_predictions" else "training_only_accuracy",
                "value": fmt(accuracy),
                "n": len(perf),
                "notes": (
                    "Computed from staged missingness_model_predictions.tsv rows with out-of-fold probabilities."
                    if metric_scope == "out_of_fold_predictions"
                    else "Computed from training-scope fitted probabilities only; not a held-out generalization metric."
                ),
            }
        )
        rows.append(
            {
                "diagnostic_scope": "missingness_model_performance",
                "metric": "out_of_fold_brier_score" if metric_scope == "out_of_fold_predictions" else "training_only_brier_score",
                "value": fmt(brier),
                "n": len(perf),
                "notes": (
                    "Lower is better; complements out-of-fold discrimination if available."
                    if metric_scope == "out_of_fold_predictions"
                    else "Lower is better, but this row is training-only unless an out-of-fold probability column is staged."
                ),
            }
        )

    merged = manifest.merge(predictions[["sample_id_canonical", "prob_interpretable"]], on="sample_id_canonical", how="left")
    merged["interpretable"] = merged["prn_interpretable"].map(as_bool)
    merged["prob_interpretable"] = pd.to_numeric(merged["prob_interpretable"], errors="coerce")
    merged["ipw_weight"] = (1.0 / merged["prob_interpretable"].clip(lower=0.05)).clip(upper=20.0)
    merged["year_numeric"] = pd.to_numeric(merged["year"], errors="coerce")
    merged["has_reads_numeric"] = merged["has_reads"].map(as_bool).astype(float)
    merged["log_total_length"] = np.log1p(pd.to_numeric(merged["total_sequence_length"], errors="coerce"))
    merged["log_n_contigs"] = np.log1p(pd.to_numeric(merged["n_contigs"], errors="coerce"))
    for covariate in ["year_numeric", "has_reads_numeric", "log_total_length", "log_n_contigs"]:
        all_mean = float(merged[covariate].mean())
        interp = merged.loc[merged["interpretable"]].copy()
        interp_mean = float(interp[covariate].mean())
        weighted = weighted_mean(interp[covariate], interp["ipw_weight"])
        rows.append(
            {
                "diagnostic_scope": "covariate_balance",
                "metric": covariate,
                "value": fmt(weighted),
                "n": int(interp[covariate].notna().sum()),
                "unweighted_interpretable_mean": fmt(interp_mean),
                "target_retained_cohort_mean": fmt(all_mean),
                "weighted_interpretable_mean": fmt(weighted),
                "notes": "IPW-weighted interpretable mean is reported against the retained-cohort target mean as a diagnostic, not a pass/fail rule.",
            }
        )

    for country in ["CHN", "JPN", "USA"]:
        subset = ipw.loc[ipw["country_iso3"] == country].copy()
        if subset.empty:
            continue
        for col in ["n_genomes_total", "n_genomes_prn_interpretable", "n_missing_outcomes", "max_ipw_weight"]:
            subset[col] = pd.to_numeric(subset[col], errors="coerce")
        rows.append(
            {
                "diagnostic_scope": "country_year_overlap",
                "country_iso3": country,
                "metric": "country_year_rows_with_genomes",
                "value": int((subset["n_genomes_total"] > 0).sum()),
                "n": int(subset["n_genomes_total"].sum()),
                "n_interpretable": int(subset["n_genomes_prn_interpretable"].sum()),
                "n_missing_outcomes": int(subset["n_missing_outcomes"].sum()),
                "max_ipw_weight": fmt(subset["max_ipw_weight"].max()),
                "notes": "Country-level overlap diagnostic for high-leverage focal countries.",
            }
        )
    for col in ["mean_ipw_weight", "max_ipw_weight", "n_missing_outcomes"]:
        ipw[col] = pd.to_numeric(ipw[col], errors="coerce")
    rows.append(
        {
            "diagnostic_scope": "weight_distribution",
            "metric": "max_ipw_weight_global",
            "value": fmt(ipw["max_ipw_weight"].max()),
            "n": int(len(ipw)),
            "notes": "Weights are truncated at the workflow-specified upper bound.",
        }
    )
    fieldnames = [
        "diagnostic_scope",
        "country_iso3",
        "metric",
        "value",
        "n",
        "n_interpretable",
        "n_missing_outcomes",
        "unweighted_interpretable_mean",
        "target_retained_cohort_mean",
        "weighted_interpretable_mean",
        "max_ipw_weight",
        "notes",
    ]
    write_tsv(SUPP_DIR / "Supplementary_Table_16_IPW_Diagnostics.tsv", fieldnames, rows)


def feature_flags(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["ptxP3"] = out.get("repo_ptxP_allele", pd.Series(index=out.index, dtype=str)).fillna("").eq("ptxP_3")
    out["fim2_present"] = out.get("repo_fim2_hash", pd.Series(index=out.index, dtype=str)).fillna("").str.strip().ne("")
    out["fim3_present"] = out.get("repo_fim3_hash", pd.Series(index=out.index, dtype=str)).fillna("").str.strip().ne("")
    out["ptxP3_fim2_fim3"] = out["ptxP3"] & out["fim2_present"] & out["fim3_present"]
    return out


def fisher_or_not_testable(origin_pos: int, origin_total: int, baseline_pos: int, baseline_total: int) -> tuple[str, str, str]:
    if origin_total == 0 or baseline_total == 0:
        return "", "", "not_testable_empty_group"
    table = [[origin_pos, origin_total - origin_pos], [baseline_pos, baseline_total - baseline_pos]]
    if sum(sum(row) for row in table) == 0 or (baseline_pos == 0 and origin_pos == 0):
        return "", "", "not_testable_no_positive_counts"
    if origin_pos == origin_total and baseline_pos == baseline_total:
        return "", "", "not_testable_no_negative_counts"
    odds_ratio, p_value = fisher_exact(table)
    return fmt(odds_ratio), f"{p_value:.4g}", "fisher_exact"


def build_genotype_background_enrichment() -> None:
    annotation = pd.read_csv(FIGURE_DATA_DIR / "published_overlap_annotation.tsv", sep="\t", dtype=str).fillna("")
    annotation = feature_flags(annotation)
    followup = pd.read_csv(SUPP_DIR / "Supplementary_Table_12_Targeted_Validation_Followup.tsv", sep="\t", dtype=str).fillna("")
    origin_accessions = set(
        followup.loc[followup["is_priority_origin_exemplar"].map(as_bool), "assembly_accession"].astype(str)
    )
    origins = annotation.loc[annotation["assembly_accession"].isin(origin_accessions)].copy()
    composition_tips = pd.read_csv(
        ROOT / "outputs" / "workflow" / "asr_sensitivity" / "composition_filtered" / "tip_states.tsv",
        sep="\t",
        dtype=str,
    ).fillna("")
    tree_accessions = set(composition_tips.loc[composition_tips["is_reference"] != "True", "assembly_accession"])
    baselines = {
        "composition_pruned_tree_nonreference": annotation.loc[annotation["assembly_accession"].isin(tree_accessions)].copy(),
        "full_interpretable_retained_cohort": annotation.loc[annotation["prn_interpretable"].map(as_bool)].copy(),
    }
    rows: list[dict[str, Any]] = []
    for baseline_name, baseline in baselines.items():
        for feature in ["ptxP3", "fim2_present", "fim3_present", "ptxP3_fim2_fim3"]:
            origin_pos = int(origins[feature].sum())
            origin_total = int(origins[feature].notna().sum())
            baseline_pos = int(baseline[feature].sum())
            baseline_total = int(baseline[feature].notna().sum())
            odds_ratio, p_value, status = fisher_or_not_testable(origin_pos, origin_total, baseline_pos, baseline_total)
            rows.append(
                {
                    "baseline_frame": baseline_name,
                    "feature": feature,
                    "origin_positive": origin_pos,
                    "origin_total": origin_total,
                    "baseline_positive": baseline_pos,
                    "baseline_total": baseline_total,
                    "odds_ratio": odds_ratio,
                    "p_value": p_value,
                    "test_status": status,
                    "notes": "Exploratory enrichment among priority origin exemplars; not interpreted as causal mechanism.",
                }
            )
    fieldnames = [
        "baseline_frame",
        "feature",
        "origin_positive",
        "origin_total",
        "baseline_positive",
        "baseline_total",
        "odds_ratio",
        "p_value",
        "test_status",
        "notes",
    ]
    write_tsv(SUPP_DIR / "Supplementary_Table_17_Genotype_Background_Enrichment.tsv", fieldnames, rows)


def root_to_tip_clock_row(tree_path: Path, tip_path: Path, frame: str) -> dict[str, Any]:
    tree = Phylo.read(str(tree_path), "newick")
    depths = tree.depths()
    tip_years = {
        row["tree_tip_label"]: as_float(row.get("year"))
        for row in read_tsv(tip_path)
        if row.get("tree_tip_label") != "Reference"
    }
    xs: list[float] = []
    ys: list[float] = []
    for terminal in tree.get_terminals():
        year = tip_years.get(terminal.name, math.nan)
        if np.isfinite(year):
            xs.append(year)
            ys.append(float(depths[terminal]))
    if len(xs) < 3:
        return {
            "tree_frame": frame,
            "n_dated_tips": len(xs),
            "slope_substitutions_per_year": "",
            "intercept": "",
            "r_squared": "",
            "timing_recommendation": "not_testable_too_few_dated_tips",
            "notes": "Root-to-tip diagnostic requires at least three dated tips.",
        }
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    predicted = slope * x + intercept
    ss_res = float(np.sum((y - predicted) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot else math.nan
    recommendation = (
        "exploratory_tip_dating_feasible"
        if slope > 0 and np.isfinite(r_squared) and r_squared >= 0.1
        else "do_not_force_tip_dated_primary;report_observed_tip_year_ranges"
    )
    return {
        "tree_frame": frame,
        "n_dated_tips": len(xs),
        "slope_substitutions_per_year": f"{slope:.8g}",
        "intercept": f"{intercept:.8g}",
        "r_squared": fmt(r_squared),
        "timing_recommendation": recommendation,
        "notes": "Simple root-to-tip linear diagnostic; not a formal dated phylogeny.",
    }


def build_clock_signal_diagnostics() -> None:
    rows = [
        root_to_tip_clock_row(
            ROOT / "outputs" / "workflow" / "asr_sensitivity" / "composition_filtered" / "rooted_ml_tree.reference_rooted.nwk",
            ROOT / "outputs" / "workflow" / "asr_sensitivity" / "composition_filtered" / "tip_states.tsv",
            "composition_pruned_reference_rooted_primary",
        ),
        root_to_tip_clock_row(
            ROOT / "outputs" / "workflow" / "asr" / "rooted_ml_tree.reference_rooted.nwk",
            ROOT / "outputs" / "workflow" / "asr" / "tip_states.tsv",
            "unpruned_reference_rooted_comparability",
        ),
    ]
    fieldnames = [
        "tree_frame",
        "n_dated_tips",
        "slope_substitutions_per_year",
        "intercept",
        "r_squared",
        "timing_recommendation",
        "notes",
    ]
    write_tsv(SUPP_DIR / "Supplementary_Table_18_Clock_Signal_Diagnostics.tsv", fieldnames, rows)


def disrupted_mechanism_bucket(mechanism_call: str) -> str:
    mechanism_call = clean_text(mechanism_call)
    if mechanism_call == "coding_disrupted_is481":
        return "IS481 insertion"
    if mechanism_call == "coding_disrupted_inversion_or_rearrangement":
        return "Inversion / rearrangement"
    if mechanism_call == "coding_disrupted_other":
        return "Other disruptions"
    return ""


def build_country_mechanism_interpretability_audit() -> None:
    manifest = pd.read_csv(ROOT / "outputs" / "workflow" / "manifest" / "manifest.tsv", sep="\t", dtype=str).fillna("")
    manifest["country_iso3"] = manifest["country_iso3"].map(clean_text).replace("", "unknown")
    manifest["country_name"] = manifest["country"].map(clean_text)
    manifest["country_name"] = manifest["country_name"].where(manifest["country_name"].ne(""), manifest["country_iso3"])
    manifest["prn_interpretable_flag"] = manifest["prn_interpretable"].map(as_bool)
    manifest["prn_disrupted_flag"] = manifest["prn_disrupted"].map(as_bool)
    manifest["mechanism_bucket"] = manifest["prn_mechanism_call"].map(disrupted_mechanism_bucket)

    rows: list[dict[str, Any]] = []
    for country_iso3, frame in sorted(manifest.groupby("country_iso3"), key=lambda item: item[0]):
        country_name_counts = Counter(clean_text(value) for value in frame["country_name"] if clean_text(value))
        country_name = country_name_counts.most_common(1)[0][0] if country_name_counts else country_iso3
        n_total = int(len(frame))
        interpretable = frame.loc[frame["prn_interpretable_flag"]].copy()
        disrupted = interpretable.loc[interpretable["prn_disrupted_flag"]].copy()
        n_interpretable = int(len(interpretable))
        n_disrupted = int(len(disrupted))
        n_intact = int(n_interpretable - n_disrupted)
        n_uninterpretable = int(n_total - n_interpretable)
        is481_n = int(disrupted["mechanism_bucket"].eq("IS481 insertion").sum())
        rearrangement_n = int(disrupted["mechanism_bucket"].eq("Inversion / rearrangement").sum())
        other_n = int(disrupted["mechanism_bucket"].eq("Other disruptions").sum())
        dominant_mechanism = (
            Counter(bucket for bucket in disrupted["mechanism_bucket"] if clean_text(bucket)).most_common(1)[0][0]
            if n_disrupted > 0
            else "none_observed"
        )

        notes: list[str] = []
        if n_interpretable < 5:
            notes.append("low_interpretable_n")
        if n_disrupted == 0:
            notes.append("no_disrupted_genomes")
        elif n_disrupted < 5:
            notes.append("low_disrupted_n")
        if n_total > 0 and n_interpretable / n_total < 0.5:
            notes.append("interpretability_below_half_of_country_records")

        rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": country_name,
                "n_retained_genomes": n_total,
                "n_prn_interpretable": n_interpretable,
                "interpretability_fraction": fmt(n_interpretable / n_total if n_total else math.nan, digits=6),
                "n_prn_disrupted": n_disrupted,
                "disrupted_fraction_within_interpretable": fmt(
                    n_disrupted / n_interpretable if n_interpretable else math.nan,
                    digits=6,
                ),
                "n_prn_intact": n_intact,
                "n_prn_uninterpretable_or_uncertain": n_uninterpretable,
                "n_disrupted_is481": is481_n,
                "n_disrupted_rearrangement": rearrangement_n,
                "n_disrupted_other": other_n,
                "is481_fraction_within_disrupted": fmt(is481_n / n_disrupted if n_disrupted else math.nan, digits=6),
                "rearrangement_fraction_within_disrupted": fmt(
                    rearrangement_n / n_disrupted if n_disrupted else math.nan,
                    digits=6,
                ),
                "other_fraction_within_disrupted": fmt(other_n / n_disrupted if n_disrupted else math.nan, digits=6),
                "dominant_disrupted_mechanism": dominant_mechanism,
                "mechanism_composition_basis": "interpretable_disrupted_genomes_only",
                "source_file": "outputs/workflow/manifest/manifest.tsv",
                "notes": ";".join(notes),
            }
        )

    fieldnames = [
        "country_iso3",
        "country_name",
        "n_retained_genomes",
        "n_prn_interpretable",
        "interpretability_fraction",
        "n_prn_disrupted",
        "disrupted_fraction_within_interpretable",
        "n_prn_intact",
        "n_prn_uninterpretable_or_uncertain",
        "n_disrupted_is481",
        "n_disrupted_rearrangement",
        "n_disrupted_other",
        "is481_fraction_within_disrupted",
        "rearrangement_fraction_within_disrupted",
        "other_fraction_within_disrupted",
        "dominant_disrupted_mechanism",
        "mechanism_composition_basis",
        "source_file",
        "notes",
    ]
    write_tsv(SUPP_DIR / "Supplementary_Table_19_Country_Interpretability_Mechanism_Audit.tsv", fieldnames, rows)
    write_tsv(FIGURE_DATA_DIR / "country_mechanism_interpretability_audit.tsv", fieldnames, rows)


def build_extended_frame_asr_summary() -> None:
    current_primary = ROOT / "outputs" / "workflow" / "asr"
    current_quality = ROOT / "outputs" / "workflow" / "asr_sensitivity" / "composition_filtered"
    balanced_primary = ROOT / "outputs" / "workflow" / "asr_balanced_ml"
    balanced_support_70 = ROOT / "outputs" / "workflow" / "asr_balanced_ml_sensitivity" / "support_70"
    balanced_support_90 = ROOT / "outputs" / "workflow" / "asr_balanced_ml_sensitivity" / "support_90"
    balanced_quality = ROOT / "outputs" / "workflow" / "asr_balanced_ml_sensitivity" / "composition_filtered"

    required = [current_primary, current_quality, balanced_primary, balanced_support_70, balanced_support_90, balanced_quality]
    if any(not path.exists() for path in required):
        return

    baseline = summarize_asr_dir(current_primary)
    scenarios = [
        {
            "scenario": "current_reference_rooted_comparability",
            "tree_scope": "current_rooted_ml_tree",
            "analysis_frame": "unpruned_reference_rooted_comparability",
            "minimum_branch_support": "0",
            "outdir": current_primary,
            "notes": "Original manuscript-scale rooted ML tree kept as the comparability frame.",
        },
        {
            "scenario": "current_composition_pruned_primary",
            "tree_scope": "current_rooted_ml_tree",
            "analysis_frame": "composition_pruned_primary_quality_frame",
            "minimum_branch_support": "0",
            "outdir": current_quality,
            "notes": "Primary quality frame used in the main manuscript after pruning 33 nonreference composition-failed tips.",
        },
        {
            "scenario": "balanced_reference_rooted_extended",
            "tree_scope": "balanced_extended_ml_tree",
            "analysis_frame": "balanced_reference_rooted_extended_frame",
            "minimum_branch_support": "0",
            "outdir": balanced_primary,
            "notes": "Expanded balanced ML tree used to test whether repeated origin survives a broader intact-tip context.",
        },
        {
            "scenario": "balanced_support_ge_70",
            "tree_scope": "balanced_extended_ml_tree",
            "analysis_frame": "balanced_reference_rooted_support_ge_70",
            "minimum_branch_support": "70",
            "outdir": balanced_support_70,
            "notes": "Balanced extended frame with Fitch origins restricted to branches with support >=70.",
        },
        {
            "scenario": "balanced_support_ge_90",
            "tree_scope": "balanced_extended_ml_tree",
            "analysis_frame": "balanced_reference_rooted_support_ge_90",
            "minimum_branch_support": "90",
            "outdir": balanced_support_90,
            "notes": "Balanced extended frame with Fitch origins restricted to branches with support >=90.",
        },
        {
            "scenario": "balanced_composition_pruned",
            "tree_scope": "balanced_extended_ml_tree",
            "analysis_frame": "balanced_composition_pruned_quality_frame",
            "minimum_branch_support": "0",
            "outdir": balanced_quality,
            "notes": "Balanced extended frame after pruning the nonreference composition-failed tips identified by IQ-TREE.",
        },
    ]

    rows: list[dict[str, Any]] = []
    for scenario in scenarios:
        summary = summarize_asr_dir(scenario["outdir"])
        rows.append(
            {
                "scenario": scenario["scenario"],
                "tree_scope": scenario["tree_scope"],
                "analysis_frame": scenario["analysis_frame"],
                "minimum_branch_support": scenario["minimum_branch_support"],
                "tip_count": summary["tip_count"],
                "disrupted_tip_count": summary["disrupted_tip_count"],
                "delta_tip_count_vs_current_reference_rooted": summary["tip_count"] - baseline["tip_count"],
                "delta_disrupted_tip_count_vs_current_reference_rooted": (
                    summary["disrupted_tip_count"] - baseline["disrupted_tip_count"]
                ),
                "fitch_origin_events": summary["fitch_origin_events"],
                "pastml_origin_events": summary["pastml_origin_events"],
                "pastml_strict_origin_events": summary["pastml_strict_origin_events"],
                "pastml_compatible_origin_events": summary["pastml_compatible_origin_events"],
                "source_dir": str(scenario["outdir"].relative_to(ROOT)),
                "notes": scenario["notes"],
            }
        )

    fieldnames = [
        "scenario",
        "tree_scope",
        "analysis_frame",
        "minimum_branch_support",
        "tip_count",
        "disrupted_tip_count",
        "delta_tip_count_vs_current_reference_rooted",
        "delta_disrupted_tip_count_vs_current_reference_rooted",
        "fitch_origin_events",
        "pastml_origin_events",
        "pastml_strict_origin_events",
        "pastml_compatible_origin_events",
        "source_dir",
        "notes",
    ]
    write_tsv(SUPP_DIR / "Supplementary_Table_23_ASR_Extended_Frame.tsv", fieldnames, rows)
    write_tsv(FIGURE_DATA_DIR / "asr_extended_frame_summary.tsv", fieldnames, rows)


def build_asr_diagnostics() -> None:
    build_asr_rooting_sensitivity()
    build_reframed_asr_sensitivity_table()
    stage_primary_origin_events()
    build_tiebreak_sensitivity()
    build_composition_failed_tip_audit()
    build_extended_frame_asr_summary()


def build_validation_diagnostics() -> None:
    extend_event_evidence_manifest()
    build_is481_event_evidence_audit()
    enrich_validation_followup_queue()
    build_origin_evidence_completeness_audit()


def build_context_diagnostics() -> None:
    build_ipw_diagnostics()
    build_genotype_background_enrichment()
    build_clock_signal_diagnostics()
    build_country_mechanism_interpretability_audit()


def main() -> None:
    build_asr_diagnostics()
    build_validation_diagnostics()
    build_context_diagnostics()


if __name__ == "__main__":
    main()
