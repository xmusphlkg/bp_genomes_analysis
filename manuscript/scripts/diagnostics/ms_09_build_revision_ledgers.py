#!/usr/bin/env python3
"""Build revision-support ledgers for missingness bounds and event definitions."""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.api as sm


sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "workflow" / "lib"))
from project_paths import project_module_data_root  # noqa: E402
from asr_parsimony import Node, assign_node_ids, parse_newick  # noqa: E402


ROOT = Path(__file__).resolve().parents[3]
SUPP_DIR = ROOT / "manuscript" / "supplementary"
FIGURE_DATA_DIR = ROOT / "manuscript" / "figure_data"
ECOLOGY_COUNTRY_YEAR_DATA = FIGURE_DATA_DIR / "ecology_country_year_observations.tsv"
PRN_COUNTRY_YEAR_DATA = FIGURE_DATA_DIR / "fig01_prn_country_year_summary.tsv"
STEP3_OUTPUTS = project_module_data_root("step3_prn_scan") / "outputs"
STEP4_OUTPUTS = project_module_data_root("step4_prn_validation") / "outputs"
STEP5_OUTPUTS = project_module_data_root("step5_phylogeny_asr") / "outputs"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def as_float(value: Any) -> float:
    if value is None:
        return math.nan
    text = str(value).strip()
    if text in {"", "---", "NA", "nan", "None"}:
        return math.nan
    return float(text)


def as_int(value: Any) -> int | None:
    number = as_float(value)
    if math.isnan(number):
        return None
    return int(round(number))


def fmt(value: Any, digits: int = 4) -> str:
    number = as_float(value)
    if math.isnan(number):
        return ""
    return f"{number:.{digits}f}"


def fraction(numerator: int | None, denominator: int | None) -> str:
    if numerator is None or denominator in {None, 0}:
        return ""
    return fmt(numerator / denominator)


def build_country_missingness_rows() -> list[dict[str, Any]]:
    country_year = pd.read_csv(PRN_COUNTRY_YEAR_DATA, sep="\t", dtype=str).fillna("")
    ecology = pd.read_csv(ECOLOGY_COUNTRY_YEAR_DATA, sep="\t", dtype=str).fillna("")
    name_lookup = (
        ecology.loc[:, ["country_iso3", "country_name"]]
        .drop_duplicates()
        .set_index("country_iso3")["country_name"]
        .to_dict()
    )

    numeric_cols = [
        "n_genomes_total",
        "n_prn_intact",
        "n_prn_disrupted",
        "n_prn_uncertain_fragmented",
        "n_prn_insufficient",
    ]
    for col in numeric_cols:
        country_year[col] = pd.to_numeric(country_year[col], errors="coerce").fillna(0).astype(int)

    grouped = (
        country_year.groupby("country_iso3", dropna=False)[numeric_cols]
        .sum()
        .reset_index()
        .sort_values(["country_iso3"])
    )
    rows: list[dict[str, Any]] = []
    for _, row in grouped.iterrows():
        country_iso3 = str(row["country_iso3"]).strip()
        if not country_iso3:
            country_iso3 = "unknown"
        n_intact = int(row["n_prn_intact"])
        n_disrupted = int(row["n_prn_disrupted"])
        n_uncertain = int(row["n_prn_uncertain_fragmented"])
        n_insufficient = int(row["n_prn_insufficient"])
        rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": name_lookup.get(country_iso3, country_iso3),
                "n_total": int(row["n_genomes_total"]),
                "n_interpretable": n_intact + n_disrupted,
                "n_uninterpretable": n_uncertain + n_insufficient,
                "n_prn_disrupted": n_disrupted,
            }
        )
    return rows


def build_missingness_bounds() -> None:
    mechanism_rows = read_tsv(SUPP_DIR / "Supplementary_Table_2_prn_mechanism_classification.tsv")
    country_rows = build_country_missingness_rows()
    country_year = pd.read_csv(ECOLOGY_COUNTRY_YEAR_DATA, sep="\t")

    total = sum(as_int(row["sample_count"]) or 0 for row in mechanism_rows)
    interpretable = sum(
        as_int(row["sample_count"]) or 0 for row in mechanism_rows if str(row["is_interpretable"]).lower() == "true"
    )
    disrupted = sum(
        as_int(row["sample_count"]) or 0
        for row in mechanism_rows
        if str(row["is_definitive_disrupted"]).lower() == "true"
    )
    uninterpretable = total - interpretable

    rows: list[dict[str, Any]] = []
    rows.append(
        {
            "analysis_scope": "global_prn_prevalence_bounds",
            "unit_id": "retained_cohort",
            "n_total": total,
            "n_interpretable": interpretable,
            "n_uninterpretable_or_uncertain": uninterpretable,
            "n_disrupted_observed": disrupted,
            "observed_interpretable_prevalence": fraction(disrupted, interpretable),
            "lower_if_all_uninterpretable_intact": fraction(disrupted, total),
            "upper_if_all_uninterpretable_disrupted": fraction(disrupted + uninterpretable, total),
            "interpretation": "Retained-cohort disruption prevalence is highly sensitive to how uninterpretable loci are assigned.",
        }
    )

    for row in mechanism_rows:
        if str(row["is_definitive_disrupted"]).lower() != "true":
            continue
        count = as_int(row["sample_count"]) or 0
        rows.append(
            {
                "analysis_scope": "mechanism_class_bounds",
                "unit_id": row["prn_mechanism_call"],
                "n_total": total,
                "n_interpretable": interpretable,
                "n_uninterpretable_or_uncertain": uninterpretable,
                "n_disrupted_observed": count,
                "observed_interpretable_prevalence": fraction(count, disrupted),
                "lower_if_all_uninterpretable_intact": fraction(count, total),
                "upper_if_all_uninterpretable_disrupted": fraction(count + uninterpretable, total),
                "interpretation": (
                    "Known mechanism composition is defined for interpretable disrupted genomes; "
                    "upper bounds assign all uninterpretable loci to this class one class at a time."
                ),
            }
        )

    for row in country_rows:
        n_total = as_int(row["n_total"])
        n_interpretable = as_int(row["n_interpretable"])
        n_uninterpretable = as_int(row["n_uninterpretable"])
        n_disrupted = as_int(row["n_prn_disrupted"]) or 0
        rows.append(
            {
                "analysis_scope": "country_prn_prevalence_bounds",
                "unit_id": row["country_iso3"],
                "country_iso3": row["country_iso3"],
                "country_name": row["country_name"],
                "n_total": n_total,
                "n_interpretable": n_interpretable,
                "n_uninterpretable_or_uncertain": n_uninterpretable,
                "n_disrupted_observed": n_disrupted,
                "observed_interpretable_prevalence": fraction(n_disrupted, n_interpretable),
                "lower_if_all_uninterpretable_intact": fraction(n_disrupted, n_total),
                "upper_if_all_uninterpretable_disrupted": fraction(
                    n_disrupted + (n_uninterpretable or 0), n_total
                ),
                "interpretation": "Country-level prevalence bound from extreme assignment of uninterpretable loci.",
            }
        )

    for _, record in country_year.iterrows():
        n_total = as_int(record.get("response_n_genomes_total", record.get("genomes_count")))
        n_interpretable = as_int(record.get("response_n_genomes_prn_interpretable", record.get("n_genomes_prn_interpretable")))
        n_disrupted = as_int(record.get("n_prn_disrupted"))
        if n_disrupted is None and n_interpretable:
            n_disrupted = int(round(as_float(record.get("response_naive_prevalence")) * n_interpretable))
        n_disrupted = n_disrupted or 0
        n_uninterpretable = None
        if n_total is not None and n_interpretable is not None:
            n_uninterpretable = max(0, n_total - n_interpretable)
        rows.append(
            {
                "analysis_scope": "country_year_prn_prevalence_bounds",
                "unit_id": f"{record.get('country_iso3')}:{as_int(record.get('year'))}",
                "country_iso3": record.get("country_iso3"),
                "country_name": record.get("country_name"),
                "year": as_int(record.get("year")),
                "n_total": n_total,
                "n_interpretable": n_interpretable,
                "n_uninterpretable_or_uncertain": n_uninterpretable,
                "n_disrupted_observed": n_disrupted,
                "observed_interpretable_prevalence": fraction(n_disrupted, n_interpretable),
                "lower_if_all_uninterpretable_intact": fraction(n_disrupted, n_total),
                "upper_if_all_uninterpretable_disrupted": fraction(n_disrupted + (n_uninterpretable or 0), n_total),
                "interpretation": "Primary ecology-panel country-year response bound under extreme missing-locus assignment.",
            }
        )

    rows.extend(build_ecology_model_bounds(country_year))

    fieldnames = [
        "analysis_scope",
        "unit_id",
        "scenario",
        "country_iso3",
        "country_name",
        "year",
        "n_total",
        "n_interpretable",
        "n_uninterpretable_or_uncertain",
        "n_disrupted_observed",
        "observed_interpretable_prevalence",
        "lower_if_all_uninterpretable_intact",
        "upper_if_all_uninterpretable_disrupted",
        "ap_exposure_v2_log_or",
        "ap_exposure_v2_ci_lower",
        "ap_exposure_v2_ci_upper",
        "ap_exposure_v2_p_value",
        "interpretation",
    ]
    write_tsv(SUPP_DIR / "Supplementary_Table_8_Missingness_Bounds.tsv", fieldnames, rows)


def zscore(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    std = numeric.std(ddof=0)
    if not np.isfinite(std) or std == 0:
        return numeric * 0.0
    return (numeric - numeric.mean()) / std


def build_ecology_model_bounds(country_year: pd.DataFrame) -> list[dict[str, Any]]:
    df = country_year.copy()
    for col in [
        "response_n_genomes_total",
        "response_n_genomes_prn_interpretable",
        "n_prn_disrupted",
        "response_naive_prevalence",
        "ap_exposure_v2_score",
        "reported_cases",
        "post_covid_period",
        "genomes_per_case_effective",
        "workflow_genomes_per_case",
    ]:
        df[col] = pd.to_numeric(df.get(col), errors="coerce")

    df["n_total"] = df["response_n_genomes_total"].fillna(df.get("genomes_count"))
    df["n_interpretable"] = df["response_n_genomes_prn_interpretable"].fillna(df.get("n_genomes_prn_interpretable"))
    df["n_disrupted"] = df["n_prn_disrupted"]
    missing_disrupted = df["n_disrupted"].isna() & df["response_naive_prevalence"].notna() & df["n_interpretable"].notna()
    df.loc[missing_disrupted, "n_disrupted"] = (
        df.loc[missing_disrupted, "response_naive_prevalence"] * df.loc[missing_disrupted, "n_interpretable"]
    ).round()
    df["n_disrupted"] = df["n_disrupted"].fillna(0)
    df["n_uninterpretable"] = (df["n_total"] - df["n_interpretable"]).clip(lower=0)
    df["log1p_reported_cases_z"] = zscore(np.log1p(df["reported_cases"]))
    df["ap_exposure_v2_score_z"] = zscore(df["ap_exposure_v2_score"])
    df["genomes_per_case_z"] = zscore(df["genomes_per_case_effective"].fillna(df["workflow_genomes_per_case"]))

    covariates = ["ap_exposure_v2_score_z", "log1p_reported_cases_z", "post_covid_period", "genomes_per_case_z"]
    keep = df.dropna(subset=["n_total", "n_interpretable", "n_disrupted", *covariates]).copy()
    scenarios = {
        "observed_interpretable_only": (keep["n_disrupted"], keep["n_interpretable"]),
        "all_uninterpretable_intact": (keep["n_disrupted"], keep["n_total"]),
        "all_uninterpretable_disrupted": (keep["n_disrupted"] + keep["n_uninterpretable"], keep["n_total"]),
    }

    rows: list[dict[str, Any]] = []
    for scenario, (successes, trials) in scenarios.items():
        valid = trials > 0
        x = sm.add_constant(keep.loc[valid, covariates], has_constant="add")
        y = successes.loc[valid] / trials.loc[valid]
        try:
            fit = sm.GLM(y, x, family=sm.families.Binomial(), freq_weights=trials.loc[valid]).fit()
            estimate = fit.params["ap_exposure_v2_score_z"]
            ci_lower, ci_upper = fit.conf_int().loc["ap_exposure_v2_score_z"]
            p_value = fit.pvalues["ap_exposure_v2_score_z"]
            interpretation = (
                "Sensitivity refit of the same four-covariate ecology structure using an extreme "
                "uninterpretable-locus assignment; not a replacement for the IPW primary model."
            )
        except Exception as exc:  # pragma: no cover - diagnostic table should still be written
            estimate = ci_lower = ci_upper = p_value = math.nan
            interpretation = f"Ecology sensitivity model failed to converge: {exc}"

        rows.append(
            {
                "analysis_scope": "ecology_coefficient_missingness_sensitivity",
                "unit_id": "ap_exposure_v2_score_z",
                "scenario": scenario,
                "n_total": int(trials.loc[valid].sum()),
                "n_interpretable": int(keep.loc[valid, "n_interpretable"].sum()),
                "n_uninterpretable_or_uncertain": int(keep.loc[valid, "n_uninterpretable"].sum()),
                "n_disrupted_observed": int(successes.loc[valid].sum()),
                "ap_exposure_v2_log_or": fmt(estimate),
                "ap_exposure_v2_ci_lower": fmt(ci_lower),
                "ap_exposure_v2_ci_upper": fmt(ci_upper),
                "ap_exposure_v2_p_value": f"{p_value:.3g}" if np.isfinite(p_value) else "",
                "interpretation": interpretation,
            }
        )
    return rows


def first_by_key(rows: list[dict[str, str]], key: str) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        value = row.get(key, "")
        if value and value not in out:
            out[value] = row
    return out


def build_event_definitions() -> None:
    figure_rows = [
        row
        for row in read_tsv(FIGURE_DATA_DIR / "figure2_prn_structural_landscape.tsv")
        if row.get("panel_id") == "event_catalog"
    ]
    event_rows = first_by_key(read_tsv(STEP4_OUTPUTS / "bp_prn_event_catalog.tsv"), "prn_event_id")
    evidence_by_accession = first_by_key(
        read_tsv(STEP3_OUTPUTS / "bp_prn_breakpoint_evidence.tsv"),
        "genome_resolved_accession",
    )
    read_is_by_event = first_by_key(
        read_tsv(STEP4_OUTPUTS / "bp_prn_read_validation_is_calls.tsv"),
        "prn_event_id",
    )
    tsd_by_event = first_by_key(
        read_tsv(STEP4_OUTPUTS / "bp_prn_read_validation_tsd.tsv"),
        "prn_event_id",
    )

    rows: list[dict[str, Any]] = []
    for row in figure_rows:
        event_id = row["prn_event_id"]
        raw_event = event_rows.get(event_id, {})
        bp_evidence = evidence_by_accession.get(row["example_assembly_accession"], {})
        read_is = read_is_by_event.get(event_id, {})
        tsd = tsd_by_event.get(event_id, {})
        rule_parts = [
            row["prn_mechanism_call"],
            row.get("bp_category", ""),
            row.get("is_element_name", ""),
            f"gap={row.get('insertion_subject_gap_bp', '')}",
            f"orientation={row.get('hit_orientation', '')}",
        ]
        rows.append(
            {
                "prn_event_id": event_id,
                "event_definition_rule": "|".join(part for part in rule_parts if part and part != "gap=" and part != "orientation="),
                "mechanism_call": row["prn_mechanism_call"],
                "event_subcategory": row.get("subcategory", ""),
                "sample_count": row.get("sample_count", ""),
                "country_count": row.get("country_count", ""),
                "year_min": row.get("year_min", ""),
                "year_max": row.get("year_max", ""),
                "bp_category": row.get("bp_category", ""),
                "insertion_subject_gap_bp": row.get("insertion_subject_gap_bp", ""),
                "is_element_name": row.get("is_element_name", ""),
                "hit_support_tier": row.get("hit_support_tier", ""),
                "hit_orientation": row.get("hit_orientation", ""),
                "call_confidence": row.get("prn_call_confidence", ""),
                "evidence_flags": raw_event.get("evidence_flags_signature", ""),
                "example_sample_id_canonical": row.get("example_sample_id_canonical", ""),
                "example_assembly_accession": row.get("example_assembly_accession", ""),
                "example_sequencing_tech": row.get("representative_sequencing_tech", ""),
                "example_contig_id": bp_evidence.get("bp_contig_id", ""),
                "example_gap_start": bp_evidence.get("bp_gap_start", ""),
                "example_gap_end": bp_evidence.get("bp_gap_end", ""),
                "read_reference_record": read_is.get("reference_record", tsd.get("chromosome", "")),
                "read_locus_start": read_is.get("locus_start", tsd.get("insertion_start", "")),
                "read_locus_end": read_is.get("locus_end", tsd.get("insertion_end", "")),
                "tsd_direct_repeats": read_is.get("direct_repeats", tsd.get("direct_repeats", "")),
                "validation_level": row.get("validation_level", ""),
                "supporting_read_or_public_longread": row.get("supporting_read_or_public_longread", ""),
                "definition_limitations": (
                    "Reference-coordinate fields are populated only when read-level or example-assembly coordinates "
                    "are present in the current ledgers; blank fields should not be interpreted as absence of a junction."
                ),
            }
        )

    fieldnames = [
        "prn_event_id",
        "event_definition_rule",
        "mechanism_call",
        "event_subcategory",
        "sample_count",
        "country_count",
        "year_min",
        "year_max",
        "bp_category",
        "insertion_subject_gap_bp",
        "is_element_name",
        "hit_support_tier",
        "hit_orientation",
        "call_confidence",
        "evidence_flags",
        "example_sample_id_canonical",
        "example_assembly_accession",
        "example_sequencing_tech",
        "example_contig_id",
        "example_gap_start",
        "example_gap_end",
        "read_reference_record",
        "read_locus_start",
        "read_locus_end",
        "tsd_direct_repeats",
        "validation_level",
        "supporting_read_or_public_longread",
        "definition_limitations",
    ]
    write_tsv(SUPP_DIR / "Supplementary_Table_9_prn_Event_Definitions.tsv", fieldnames, rows)


def preferred_binary_state(states: set[str]) -> str:
    for state in ["intact", "disrupted", "insufficient_data", "uncertain"]:
        if state in states:
            return state
    return sorted(states)[0] if states else "uncertain"


def fitch_downpass_sets(node: Node, tip_state_sets: dict[str, set[str]]) -> set[str]:
    if node.is_tip:
        node.candidate_states = set(tip_state_sets.get(node.tree_label, {"intact", "disrupted"}))
        return node.candidate_states
    child_sets = [fitch_downpass_sets(child, tip_state_sets) for child in node.children]
    intersection = set.intersection(*child_sets)
    node.candidate_states = intersection if intersection else set.union(*child_sets)
    return node.candidate_states


def fitch_uppass_sets(node: Node, parent_state: str | None = None) -> None:
    if parent_state and parent_state in node.candidate_states:
        node.inferred_state = parent_state
    else:
        node.inferred_state = preferred_binary_state(node.candidate_states)
    for child in node.children:
        fitch_uppass_sets(child, node.inferred_state)


def prune_tree(node: Node, drop_labels: set[str]) -> Node | None:
    if node.is_tip:
        return None if node.tree_label in drop_labels else node
    kept_children: list[Node] = []
    for child in node.children:
        kept = prune_tree(child, drop_labels)
        if kept is None:
            continue
        kept.parent = node
        kept_children.append(kept)
    node.children = kept_children
    if not kept_children:
        return None
    return node


def count_fitch_origins(tree_text: str, tip_state_sets: dict[str, set[str]], drop_labels: set[str] | None = None) -> int:
    root = parse_newick(tree_text)
    if drop_labels:
        root = prune_tree(root, drop_labels) or root
        root.parent = None
    nodes = assign_node_ids(root)
    fitch_downpass_sets(root, tip_state_sets)
    fitch_uppass_sets(root)
    return sum(
        1
        for node in nodes
        if node.parent is not None and node.parent.inferred_state == "intact" and node.inferred_state == "disrupted"
    )


def build_asr_trait_coding_sensitivity() -> None:
    tree_text = (ROOT / "outputs" / "workflow" / "asr" / "rooted_ml_tree.reference_rooted.nwk").read_text()
    tip_rows = read_tsv(ROOT / "outputs" / "workflow" / "asr" / "tip_states.tsv")
    base_states = {row["tree_tip_label"]: row["prn_state"] for row in tip_rows}
    insufficient_labels = {label for label, state in base_states.items() if state == "insufficient_data"}

    scenario_maps: dict[str, tuple[dict[str, set[str]], set[str], str]] = {
        "third_state_primary": (
            {label: {state} for label, state in base_states.items()},
            set(),
            "Primary Fitch coding with insufficient_data retained as an explicit third state.",
        ),
        "insufficient_as_missing_ambiguous": (
            {
                label: ({"intact", "disrupted"} if state == "insufficient_data" else {state})
                for label, state in base_states.items()
            },
            set(),
            "Insufficient-data tips treated as ambiguous between intact and disrupted.",
        ),
        "insufficient_as_intact": (
            {label: ({"intact"} if state == "insufficient_data" else {state}) for label, state in base_states.items()},
            set(),
            "Insufficient-data tips assigned to intact.",
        ),
        "insufficient_as_disrupted": (
            {label: ({"disrupted"} if state == "insufficient_data" else {state}) for label, state in base_states.items()},
            set(),
            "Insufficient-data tips assigned to disrupted.",
        ),
        "insufficient_pruned": (
            {label: {state} for label, state in base_states.items() if state != "insufficient_data"},
            insufficient_labels,
            "Insufficient-data tips pruned from the rooted tree before Fitch counting.",
        ),
    }

    rows: list[dict[str, Any]] = []
    for scenario, (state_sets, drop_labels, notes) in scenario_maps.items():
        origins = count_fitch_origins(tree_text, state_sets, drop_labels)
        retained_states = {
            label: next(iter(states)) if len(states) == 1 else "ambiguous"
            for label, states in state_sets.items()
            if label not in drop_labels
        }
        rows.append(
            {
                "scenario": scenario,
                "tip_count": len(retained_states),
                "disrupted_tip_count": sum(1 for state in retained_states.values() if state == "disrupted"),
                "insufficient_or_ambiguous_tip_count": sum(
                    1 for state in retained_states.values() if state in {"insufficient_data", "ambiguous"}
                ),
                "fitch_origin_events": origins,
                "notes": notes,
            }
        )

    write_tsv(
        SUPP_DIR / "Supplementary_Table_10_ASR_Trait_Coding_Sensitivity.tsv",
        [
            "scenario",
            "tip_count",
            "disrupted_tip_count",
            "insufficient_or_ambiguous_tip_count",
            "fitch_origin_events",
            "notes",
        ],
        rows,
    )


def export_tsv(source: Path, destination: Path) -> None:
    rows = read_tsv(source)
    if not rows:
        raise ValueError(f"no rows found in {source}")
    write_tsv(destination, list(rows[0].keys()), rows)


def export_local_neighborhood_outputs() -> None:
    local_dir = STEP5_OUTPUTS
    export_tsv(
        local_dir / "bp_prn_local_neighborhood_summary.tsv",
        SUPP_DIR / "Supplementary_Table_11_Local_Neighborhood_Tree_Check.tsv",
    )
    export_tsv(
        local_dir / "bp_prn_local_neighborhood_tip_selection.tsv",
        FIGURE_DATA_DIR / "figure3_local_neighborhood_tip_selection.tsv",
    )
    export_tsv(
        local_dir / "bp_prn_local_neighborhood_origin_events.tsv",
        FIGURE_DATA_DIR / "figure3_local_neighborhood_origin_events.tsv",
    )


def export_validation_followup_queue() -> None:
    source = STEP4_OUTPUTS / "bp_prn_targeted_validation_followup_queue.tsv"
    export_tsv(source, SUPP_DIR / "Supplementary_Table_12_Targeted_Validation_Followup.tsv")
    export_tsv(source, FIGURE_DATA_DIR / "figure6_targeted_validation_followup.tsv")


def main() -> None:
    build_missingness_bounds()
    build_event_definitions()
    build_asr_trait_coding_sensitivity()
    export_local_neighborhood_outputs()
    export_validation_followup_queue()


if __name__ == "__main__":
    main()
