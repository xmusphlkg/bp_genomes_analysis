#!/usr/bin/env python3
"""Assign deterministic PRN confidence tiers from PRN-02/03 outputs."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path

from step4_02_scan_prn_mechanisms import (
    EVENT_COLUMNS,
    MECHANISM_COLUMNS,
    STEP4_DATA_ROOT,
    load_tsv_rows,
    normalize_text,
    parse_float,
    parse_int,
    write_tsv,
)


def load_best_hit_by_sample(path: Path) -> dict[str, dict[str, str]]:
    best_hit_by_sample: dict[str, dict[str, str]] = {}
    for row in load_tsv_rows(path):
        if normalize_text(row.get("is_best_hit", "")).casefold() != "true":
            continue
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        if sample_id:
            best_hit_by_sample[sample_id] = row
    return best_hit_by_sample


def parse_tokens(value: str) -> list[str]:
    return [token for token in normalize_text(value).split(";") if token]


def build_rule_note(notes: str, rule_id: str) -> str:
    tokens: list[str] = []
    for token in parse_tokens(notes):
        if token == "confidence_is_pre_prn04_rule_tier":
            continue
        if token.startswith("prn04_rule="):
            continue
        tokens.append(token)
    tokens.append(f"prn04_rule={rule_id}")
    return ";".join(dict.fromkeys(tokens))


def bool_field(value: str) -> bool:
    return normalize_text(value).casefold() == "true"


def score_prn_call(
    mechanism_row: dict[str, str],
    best_hit_row: dict[str, str] | None,
) -> tuple[str, str]:
    mechanism_call = normalize_text(mechanism_row.get("prn_mechanism_call", ""))
    prn_call_initial = normalize_text(mechanism_row.get("prn_call_initial", ""))
    bp_category = normalize_text(mechanism_row.get("bp_category", ""))
    query_cov = parse_float(mechanism_row.get("prn_query_cov_pct", "")) or 0.0
    best_single_cov = parse_float(mechanism_row.get("prn_best_single_cov_pct", "")) or 0.0
    prn_hsp_n = parse_int(mechanism_row.get("prn_hsp_n", "")) or 0
    evidence_flags = set(parse_tokens(mechanism_row.get("evidence_flags", "")))

    hit_support_tier = normalize_text((best_hit_row or {}).get("hit_support_tier", ""))
    passes_prn02_support_rule = bool_field((best_hit_row or {}).get("passes_prn02_support_rule", ""))
    supports_assigned_mechanism = bool_field((best_hit_row or {}).get("supports_assigned_mechanism", ""))

    if mechanism_call == "insufficient_data":
        if prn_call_initial == "not_available_current_step3" or "no_step3_prn_input" in evidence_flags:
            return "insufficient_evidence", "no_current_step3_prn_input"
        if prn_call_initial == "missing_fasta" or "missing_fasta" in evidence_flags:
            return "insufficient_evidence", "assembly_sequence_unavailable"
        if prn_call_initial == "partial" or "partial_prn_call" in evidence_flags:
            return "insufficient_evidence", "partial_prn_alignment_no_structural_upgrade"
        return "insufficient_evidence", "other_insufficient_signal"

    if mechanism_call == "intact":
        if prn_call_initial == "intact" and prn_hsp_n == 1 and query_cov >= 95.0 and best_single_cov >= 95.0:
            return "assembly_high", "intact_single_hsp_ge95cov"
        return "assembly_moderate", "intact_but_suboptimal_alignment"

    if mechanism_call == "coding_disrupted_is481":
        if (
            bp_category == "insertion_like"
            and prn_hsp_n >= 2
            and query_cov >= 95.0
            and passes_prn02_support_rule
            and supports_assigned_mechanism
            and hit_support_tier == "strong"
        ):
            return "assembly_high", "is481_supported_strong"
        if (
            bp_category == "insertion_like"
            and prn_hsp_n >= 2
            and query_cov >= 95.0
            and passes_prn02_support_rule
            and supports_assigned_mechanism
            and hit_support_tier == "moderate"
        ):
            return "assembly_moderate", "is481_supported_moderate"
        return "assembly_low", "is481_label_without_supported_best_hit"

    if mechanism_call == "coding_disrupted_other_is":
        if (
            bp_category == "insertion_like"
            and prn_hsp_n >= 2
            and query_cov >= 95.0
            and passes_prn02_support_rule
            and supports_assigned_mechanism
            and hit_support_tier == "strong"
        ):
            return "assembly_moderate", "other_is_supported_strong"
        if (
            bp_category == "insertion_like"
            and prn_hsp_n >= 2
            and query_cov >= 95.0
            and passes_prn02_support_rule
            and supports_assigned_mechanism
            and hit_support_tier == "moderate"
        ):
            return "assembly_low", "other_is_supported_moderate"
        return "assembly_low", "other_is_incomplete_support_profile"

    if mechanism_call == "coding_disrupted_inversion_or_rearrangement":
        if bp_category == "within_contig" and prn_hsp_n >= 2 and query_cov >= 95.0:
            return "assembly_moderate", "within_contig_split_alignment"
        return "assembly_low", "rearrangement_but_incomplete_structural_support"

    if mechanism_call == "coding_disrupted_other":
        if bp_category == "insertion_like" and prn_hsp_n >= 2 and query_cov >= 95.0:
            return "assembly_low", "insertion_like_without_supported_is"
        return "assembly_low", "other_disruption_limited_support"

    if mechanism_call == "uncertain_fragmented_assembly":
        return "assembly_low", "fragmented_or_near_contig_end"

    return "insufficient_evidence", "fallback_unclassified"


def build_event_catalog(scored_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in scored_rows:
        grouped[row["prn_event_id"]].append(row)

    catalog_rows: list[dict[str, str]] = []
    for event_id, members in sorted(grouped.items(), key=lambda item: (-len(item[1]), item[0])):
        years = [parse_int(row.get("year", "")) for row in members if parse_int(row.get("year", "")) is not None]
        countries = sorted({row.get("country_iso3", "") for row in members if normalize_text(row.get("country_iso3", ""))})
        assemblies = sorted({row.get("assembly_accession", "") for row in members if normalize_text(row.get("assembly_accession", ""))})
        confidence_counts = Counter(row.get("prn_call_confidence", "") for row in members)
        confidence_mode = confidence_counts.most_common(1)[0][0]
        confidence_values = ",".join(sorted(confidence_counts))
        notes = [f"prn04_confidence_values={confidence_values}"]
        base_rule_notes = sorted(
            {
                token
                for row in members
                for token in parse_tokens(row.get("notes", ""))
                if token.startswith("prn04_rule=")
            }
        )
        notes.extend(base_rule_notes)

        catalog_rows.append(
            {
                "prn_event_id": event_id,
                "prn_mechanism_call": members[0]["prn_mechanism_call"],
                "prn_call_confidence_mode": confidence_mode,
                "bp_category": members[0]["bp_category"],
                "is_element_best_hit": members[0]["is_element_best_hit"],
                "insertion_subject_gap_bp": members[0]["insertion_subject_gap_bp"],
                "evidence_flags_signature": members[0]["evidence_flags"],
                "sample_count": str(len(members)),
                "assembly_count": str(len(assemblies)),
                "country_iso3_values": ";".join(countries),
                "year_min": "" if not years else str(min(years)),
                "year_max": "" if not years else str(max(years)),
                "example_sample_id_canonical": members[0]["sample_id_canonical"],
                "example_assembly_accession": members[0]["assembly_accession"],
                "notes": ";".join(notes),
            }
        )
    return catalog_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Assign deterministic PRN confidence tiers using PRN-02 mechanism calls and "
            "PRN-03 best-hit support summaries."
        )
    )
    parser.add_argument(
        "--mechanism-calls",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_mechanism_calls.tsv",
        help="PRN-02 mechanism calls TSV to rescore.",
    )
    parser.add_argument(
        "--is-hits",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_is_hits.tsv",
        help="PRN-03 insertion-sequence hit summary TSV.",
    )
    parser.add_argument(
        "--mechanism-out",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_mechanism_calls.tsv",
        help="Rescored mechanism output path.",
    )
    parser.add_argument(
        "--event-catalog-out",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_event_catalog.tsv",
        help="Updated event catalog output path.",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    mechanism_rows = load_tsv_rows(args.mechanism_calls)
    best_hit_by_sample = load_best_hit_by_sample(args.is_hits)

    scored_rows: list[dict[str, str]] = []
    for row in mechanism_rows:
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        confidence_label, rule_id = score_prn_call(row, best_hit_by_sample.get(sample_id))
        rescored = dict(row)
        rescored["prn_call_confidence"] = confidence_label
        rescored["notes"] = build_rule_note(row.get("notes", ""), rule_id)
        scored_rows.append(rescored)

    write_tsv(args.mechanism_out, MECHANISM_COLUMNS, scored_rows)
    write_tsv(args.event_catalog_out, EVENT_COLUMNS, build_event_catalog(scored_rows))

    confidence_counts = Counter(row["prn_call_confidence"] for row in scored_rows)
    for label, count in sorted(confidence_counts.items()):
        print(f"{label}\t{count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
