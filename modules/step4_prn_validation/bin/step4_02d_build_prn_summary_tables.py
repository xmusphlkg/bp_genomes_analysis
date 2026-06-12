#!/usr/bin/env python3
"""Build figure-ready and integration-ready PRN summary tables."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from step4_02_scan_prn_mechanisms import (
    STEP4_DATA_ROOT,
    load_tsv_rows,
    normalize_text,
    parse_int,
    write_tsv,
)


MECHANISM_SUMMARY_COLUMNS = [
    "prn_mechanism_call",
    "sample_count",
    "event_count",
    "country_count",
    "year_min",
    "year_max",
    "is_interpretable",
    "is_definitive_disrupted",
    "is_uncertain_fragmented",
    "sample_fraction_all",
    "sample_fraction_interpretable",
]

CONFIDENCE_SUMMARY_COLUMNS = [
    "prn_mechanism_call",
    "prn_call_confidence",
    "sample_count",
    "event_count",
    "country_count",
    "year_min",
    "year_max",
    "sample_fraction_all",
    "sample_fraction_within_mechanism",
]

BREAKPOINT_SUMMARY_COLUMNS = [
    "prn_mechanism_call",
    "prn_call_confidence",
    "bp_category",
    "insertion_size_bin",
    "insertion_subject_gap_bp",
    "is_element_best_hit",
    "best_hit_support_tier",
    "sample_count",
    "event_count",
    "country_count",
    "year_min",
    "year_max",
]

UNRESOLVED_SUMMARY_COLUMNS = [
    "prn_mechanism_call",
    "prn_call_initial",
    "prn_call_confidence",
    "prn04_rule",
    "bp_category",
    "sample_count",
    "event_count",
    "country_count",
    "year_min",
    "year_max",
    "example_prn_event_id",
    "example_assembly_accession",
]

COUNTRY_YEAR_SUMMARY_COLUMNS = [
    "country_iso3",
    "year",
    "country_year_key",
    "n_genomes_total",
    "n_genomes_prn_interpretable",
    "n_prn_intact",
    "n_prn_disrupted",
    "n_prn_uncertain_fragmented",
    "n_prn_insufficient",
    "n_prn_disrupted_assembly_high",
    "n_prn_disrupted_assembly_moderate",
    "n_prn_disrupted_assembly_low",
    "frac_prn_disrupted",
]

DEFINITIVE_DISRUPTED_CALLS = {
    "coding_disrupted_is481",
    "coding_disrupted_other_is",
    "coding_disrupted_inversion_or_rearrangement",
    "coding_disrupted_other",
}


def parse_tokens(value: str) -> list[str]:
    return [token for token in normalize_text(value).split(";") if token]


def parse_rule_id(notes: str) -> str:
    for token in parse_tokens(notes):
        if token.startswith("prn04_rule="):
            return token.split("=", 1)[1]
    return ""


def country_label(value: str) -> str:
    return normalize_text(value) or "unknown"


def year_label(value: str) -> str:
    return normalize_text(value) or "unknown"


def collect_year_bounds(rows: list[dict[str, str]]) -> tuple[str, str]:
    years = [parse_int(row.get("year", "")) for row in rows if parse_int(row.get("year", "")) is not None]
    if not years:
        return "", ""
    return str(min(years)), str(max(years))


def format_fraction(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return ""
    return f"{numerator / denominator:.6f}"


def is_interpretable(row: dict[str, str]) -> bool:
    mechanism = normalize_text(row.get("prn_mechanism_call", ""))
    return mechanism in {"intact", *DEFINITIVE_DISRUPTED_CALLS}


def is_definitive_disrupted(row: dict[str, str]) -> bool:
    return normalize_text(row.get("prn_mechanism_call", "")) in DEFINITIVE_DISRUPTED_CALLS


def is_uncertain_fragmented(row: dict[str, str]) -> bool:
    return normalize_text(row.get("prn_mechanism_call", "")) == "uncertain_fragmented_assembly"


def insertion_size_bin(value: str) -> str:
    gap_bp = parse_int(value)
    if gap_bp is None:
        return "not_applicable"
    if gap_bp < 50:
        return "lt50"
    if gap_bp < 200:
        return "50_199"
    if gap_bp < 500:
        return "200_499"
    if gap_bp < 900:
        return "500_899"
    if gap_bp < 1100:
        return "900_1099"
    return "ge1100"


def best_hit_by_sample(is_hit_rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in is_hit_rows:
        if normalize_text(row.get("is_best_hit", "")).casefold() != "true":
            continue
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        if sample_id:
            lookup[sample_id] = row
    return lookup


def build_mechanism_summary(mechanism_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in mechanism_rows:
        grouped[row["prn_mechanism_call"]].append(row)

    total_rows = len(mechanism_rows)
    interpretable_total = sum(1 for row in mechanism_rows if is_interpretable(row))
    summary_rows: list[dict[str, str]] = []

    for mechanism_call in sorted(grouped):
        members = grouped[mechanism_call]
        year_min, year_max = collect_year_bounds(members)
        summary_rows.append(
            {
                "prn_mechanism_call": mechanism_call,
                "sample_count": str(len(members)),
                "event_count": str(len({row["prn_event_id"] for row in members if normalize_text(row.get("prn_event_id", ""))})),
                "country_count": str(len({row["country_iso3"] for row in members if normalize_text(row.get("country_iso3", ""))})),
                "year_min": year_min,
                "year_max": year_max,
                "is_interpretable": "true" if is_interpretable(members[0]) else "false",
                "is_definitive_disrupted": "true" if mechanism_call in DEFINITIVE_DISRUPTED_CALLS else "false",
                "is_uncertain_fragmented": "true" if mechanism_call == "uncertain_fragmented_assembly" else "false",
                "sample_fraction_all": format_fraction(len(members), total_rows),
                "sample_fraction_interpretable": format_fraction(len(members), interpretable_total)
                if is_interpretable(members[0])
                else "",
            }
        )
    return summary_rows


def build_confidence_summary(mechanism_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    totals_by_mechanism: dict[str, int] = defaultdict(int)
    total_rows = len(mechanism_rows)

    for row in mechanism_rows:
        key = (row["prn_mechanism_call"], row["prn_call_confidence"])
        grouped[key].append(row)
        totals_by_mechanism[row["prn_mechanism_call"]] += 1

    summary_rows: list[dict[str, str]] = []
    for mechanism_call, confidence in sorted(grouped):
        members = grouped[(mechanism_call, confidence)]
        year_min, year_max = collect_year_bounds(members)
        summary_rows.append(
            {
                "prn_mechanism_call": mechanism_call,
                "prn_call_confidence": confidence,
                "sample_count": str(len(members)),
                "event_count": str(len({row["prn_event_id"] for row in members if normalize_text(row.get("prn_event_id", ""))})),
                "country_count": str(len({row["country_iso3"] for row in members if normalize_text(row.get("country_iso3", ""))})),
                "year_min": year_min,
                "year_max": year_max,
                "sample_fraction_all": format_fraction(len(members), total_rows),
                "sample_fraction_within_mechanism": format_fraction(len(members), totals_by_mechanism[mechanism_call]),
            }
        )
    return summary_rows


def build_breakpoint_summary(
    mechanism_rows: list[dict[str, str]],
    hit_lookup: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str, str, str, str, str, str], list[dict[str, str]]] = defaultdict(list)

    for row in mechanism_rows:
        mechanism_call = normalize_text(row.get("prn_mechanism_call", ""))
        if mechanism_call in {"intact", "insufficient_data"}:
            continue
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        hit_row = hit_lookup.get(sample_id)
        grouped[
            (
                mechanism_call,
                row["prn_call_confidence"],
                normalize_text(row.get("bp_category", "")) or "not_applicable",
                insertion_size_bin(row.get("insertion_subject_gap_bp", "")),
                normalize_text(row.get("insertion_subject_gap_bp", "")) or "not_applicable",
                normalize_text(row.get("is_element_best_hit", "")) or "none",
                normalize_text((hit_row or {}).get("hit_support_tier", "")) or "not_applicable",
            )
        ].append(row)

    summary_rows: list[dict[str, str]] = []
    for key in sorted(grouped):
        members = grouped[key]
        year_min, year_max = collect_year_bounds(members)
        summary_rows.append(
            {
                "prn_mechanism_call": key[0],
                "prn_call_confidence": key[1],
                "bp_category": key[2],
                "insertion_size_bin": key[3],
                "insertion_subject_gap_bp": key[4],
                "is_element_best_hit": key[5],
                "best_hit_support_tier": key[6],
                "sample_count": str(len(members)),
                "event_count": str(len({row["prn_event_id"] for row in members if normalize_text(row.get("prn_event_id", ""))})),
                "country_count": str(len({row["country_iso3"] for row in members if normalize_text(row.get("country_iso3", ""))})),
                "year_min": year_min,
                "year_max": year_max,
            }
        )
    return summary_rows


def build_unresolved_summary(mechanism_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str, str, str, str], list[dict[str, str]]] = defaultdict(list)

    for row in mechanism_rows:
        mechanism_call = normalize_text(row.get("prn_mechanism_call", ""))
        if mechanism_call not in {"insufficient_data", "uncertain_fragmented_assembly"}:
            continue
        grouped[
            (
                mechanism_call,
                normalize_text(row.get("prn_call_initial", "")) or "not_applicable",
                normalize_text(row.get("prn_call_confidence", "")),
                parse_rule_id(row.get("notes", "")) or "unknown_rule",
                normalize_text(row.get("bp_category", "")) or "not_applicable",
            )
        ].append(row)

    summary_rows: list[dict[str, str]] = []
    for key in sorted(grouped):
        members = grouped[key]
        year_min, year_max = collect_year_bounds(members)
        summary_rows.append(
            {
                "prn_mechanism_call": key[0],
                "prn_call_initial": key[1],
                "prn_call_confidence": key[2],
                "prn04_rule": key[3],
                "bp_category": key[4],
                "sample_count": str(len(members)),
                "event_count": str(len({row["prn_event_id"] for row in members if normalize_text(row.get("prn_event_id", ""))})),
                "country_count": str(len({row["country_iso3"] for row in members if normalize_text(row.get("country_iso3", ""))})),
                "year_min": year_min,
                "year_max": year_max,
                "example_prn_event_id": members[0]["prn_event_id"],
                "example_assembly_accession": members[0]["assembly_accession"],
            }
        )
    return summary_rows


def build_country_year_summary(mechanism_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)

    for row in mechanism_rows:
        country = country_label(row.get("country_iso3", ""))
        year = year_label(row.get("year", ""))
        grouped[(country, year)].append(row)

    summary_rows: list[dict[str, str]] = []
    for country, year in sorted(
        grouped,
        key=lambda item: (
            item[0],
            int(item[1]) if item[1].isdigit() else 999999,
            item[1],
        ),
    ):
        members = grouped[(country, year)]
        n_total = len(members)
        n_interpretable = sum(1 for row in members if is_interpretable(row))
        n_intact = sum(1 for row in members if normalize_text(row.get("prn_mechanism_call", "")) == "intact")
        n_disrupted = sum(1 for row in members if is_definitive_disrupted(row))
        n_uncertain = sum(1 for row in members if is_uncertain_fragmented(row))
        n_insufficient = sum(1 for row in members if normalize_text(row.get("prn_mechanism_call", "")) == "insufficient_data")
        n_disrupted_high = sum(
            1
            for row in members
            if is_definitive_disrupted(row) and normalize_text(row.get("prn_call_confidence", "")) == "assembly_high"
        )
        n_disrupted_moderate = sum(
            1
            for row in members
            if is_definitive_disrupted(row) and normalize_text(row.get("prn_call_confidence", "")) == "assembly_moderate"
        )
        n_disrupted_low = sum(
            1
            for row in members
            if is_definitive_disrupted(row) and normalize_text(row.get("prn_call_confidence", "")) == "assembly_low"
        )

        summary_rows.append(
            {
                "country_iso3": country,
                "year": year,
                "country_year_key": f"{country}|{year}",
                "n_genomes_total": str(n_total),
                "n_genomes_prn_interpretable": str(n_interpretable),
                "n_prn_intact": str(n_intact),
                "n_prn_disrupted": str(n_disrupted),
                "n_prn_uncertain_fragmented": str(n_uncertain),
                "n_prn_insufficient": str(n_insufficient),
                "n_prn_disrupted_assembly_high": str(n_disrupted_high),
                "n_prn_disrupted_assembly_moderate": str(n_disrupted_moderate),
                "n_prn_disrupted_assembly_low": str(n_disrupted_low),
                "frac_prn_disrupted": format_fraction(n_disrupted, n_interpretable),
            }
        )
    return summary_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Materialize PRN summary tables for mechanism composition, breakpoint patterns, "
            "unresolved cases, and country-year aggregation."
        )
    )
    parser.add_argument(
        "--mechanism-calls",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_mechanism_calls.tsv",
        help="PRN mechanism calls table from PRN-02/04.",
    )
    parser.add_argument(
        "--is-hits",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_is_hits.tsv",
        help="PRN IS-hit summary table from PRN-03.",
    )
    parser.add_argument(
        "--mechanism-summary-out",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_mechanism_summary.tsv",
        help="Overall mechanism summary TSV.",
    )
    parser.add_argument(
        "--confidence-summary-out",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_confidence_summary.tsv",
        help="Mechanism-by-confidence summary TSV.",
    )
    parser.add_argument(
        "--breakpoint-summary-out",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_breakpoint_summary.tsv",
        help="Breakpoint and insertion-pattern summary TSV.",
    )
    parser.add_argument(
        "--unresolved-summary-out",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_unresolved_summary.tsv",
        help="Unresolved and uncertain-call summary TSV.",
    )
    parser.add_argument(
        "--country-year-summary-out",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_country_year_summary.tsv",
        help="Country-year genomic summary TSV for later integration tasks.",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    mechanism_rows = load_tsv_rows(args.mechanism_calls)
    is_hit_rows = load_tsv_rows(args.is_hits)
    hit_lookup = best_hit_by_sample(is_hit_rows)

    mechanism_summary_rows = build_mechanism_summary(mechanism_rows)
    confidence_summary_rows = build_confidence_summary(mechanism_rows)
    breakpoint_summary_rows = build_breakpoint_summary(mechanism_rows, hit_lookup)
    unresolved_summary_rows = build_unresolved_summary(mechanism_rows)
    country_year_summary_rows = build_country_year_summary(mechanism_rows)

    write_tsv(args.mechanism_summary_out, MECHANISM_SUMMARY_COLUMNS, mechanism_summary_rows)
    write_tsv(args.confidence_summary_out, CONFIDENCE_SUMMARY_COLUMNS, confidence_summary_rows)
    write_tsv(args.breakpoint_summary_out, BREAKPOINT_SUMMARY_COLUMNS, breakpoint_summary_rows)
    write_tsv(args.unresolved_summary_out, UNRESOLVED_SUMMARY_COLUMNS, unresolved_summary_rows)
    write_tsv(args.country_year_summary_out, COUNTRY_YEAR_SUMMARY_COLUMNS, country_year_summary_rows)

    print(f"mechanism_summary_rows\t{len(mechanism_summary_rows)}")
    print(f"confidence_summary_rows\t{len(confidence_summary_rows)}")
    print(f"breakpoint_summary_rows\t{len(breakpoint_summary_rows)}")
    print(f"unresolved_summary_rows\t{len(unresolved_summary_rows)}")
    print(f"country_year_summary_rows\t{len(country_year_summary_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
