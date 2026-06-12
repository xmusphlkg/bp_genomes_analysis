#!/usr/bin/env python3
"""Build figure-ready and integration-ready PRN summary tables with confidence intervals.

IMPROVED VERSION: Adds Wilson score confidence intervals for disruption rates
to quantify uncertainty, especially important for small sample sizes.
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path
from pathlib import PurePath
from typing import List, Dict, Tuple

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
    "sample_fraction_within_mechanism",
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
    "frac_prn_disrupted_ci_lower",  # NEW: 95% CI lower bound
    "frac_prn_disrupted_ci_upper",  # NEW: 95% CI upper bound
    "frac_prn_disrupted_n",         # NEW: Sample size for transparency
]

DEFINITIVE_DISRUPTED_CALLS = {
    "coding_disrupted_is481",
    "coding_disrupted_other_is",
    "coding_disrupted_inversion_or_rearrangement",
    "coding_disrupted_other",
}


def wilson_score_interval(successes: int, n: int, z: float = 1.96) -> Tuple[float, float]:
    """
    Calculate Wilson score confidence interval for a proportion.
    
    The Wilson score interval provides better coverage than the Wald interval,
    especially for small sample sizes or proportions near 0 or 1.
    
    Args:
        successes: Number of successes (e.g., disrupted PRN calls)
        n: Total sample size (e.g., interpretable genomes)
        z: Z-score for desired confidence level (1.96 for 95% CI)
    
    Returns:
        Tuple of (lower_bound, upper_bound)
    """
    if n <= 0:
        return (float('nan'), float('nan'))
    
    p = successes / n
    
    # Wilson score interval formula
    denominator = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denominator
    margin = (z / denominator) * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))
    
    lower = max(0.0, centre - margin)
    upper = min(1.0, centre + margin)
    
    return (lower, upper)


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
        return "unknown"
    if gap_bp <= 0:
        return "0"
    if gap_bp <= 100:
        return "1-100"
    if gap_bp <= 500:
        return "101-500"
    if gap_bp <= 1000:
        return "501-1000"
    return ">1000"


def build_country_year_summaries_with_ci(
    mechanism_rows: List[dict[str, str]],
    is_hit_rows: List[dict[str, str]]
) -> List[dict[str, str]]:
    """Build country-year summaries WITH Wilson confidence intervals."""
    
    is_hits_by_sample = {
        row["assembly_accession"]: row
        for row in is_hit_rows
    }
    
    grouped_data: Dict[Tuple[str, int], List[dict]] = defaultdict(list)
    for row in mechanism_rows:
        country = row.get("country_iso3", "")
        year_str = row.get("year", "")
        year = parse_int(year_str)
        if not country or year is None:
            continue
        grouped_data[(country, year)].append(row)
    
    summary_rows = []
    for (country, year), rows in sorted(grouped_data.items()):
        n_total = len(rows)
        n_interpretable = sum(1 for r in rows if is_interpretable(r))
        n_intact = sum(1 for r in rows if normalize_text(r.get("prn_mechanism_call", "")) == "intact")
        n_disrupted = sum(1 for r in rows if is_definitive_disrupted(r))
        n_uncertain = sum(1 for r in rows if is_uncertain_fragmented(r))
        n_insufficient = n_total - n_interpretable - n_uncertain
        
        assembly_quality_counts = {"high": 0, "moderate": 0, "low": 0}
        for row in rows:
            if is_definitive_disrupted(row):
                acc = row.get("assembly_accession", "")
                is_row = is_hits_by_sample.get(acc, {})
                quality = normalize_text(is_row.get("overall_assembly_quality", ""))
                if quality in assembly_quality_counts:
                    assembly_quality_counts[quality] += 1
        
        n_disrupted_high = assembly_quality_counts["high"]
        n_disrupted_moderate = assembly_quality_counts["moderate"]
        n_disrupted_low = assembly_quality_counts["low"]
        
        # Calculate Wilson confidence interval for disruption rate
        frac_disrupted = n_disrupted / n_interpretable if n_interpretable > 0 else float('nan')
        ci_lower, ci_upper = wilson_score_interval(n_disrupted, n_interpretable)
        
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
                "frac_prn_disrupted_ci_lower": f"{ci_lower:.6f}" if not math.isnan(ci_lower) else "",
                "frac_prn_disrupted_ci_upper": f"{ci_upper:.6f}" if not math.isnan(ci_upper) else "",
                "frac_prn_disrupted_n": str(n_interpretable),  # Explicit sample size
            }
        )
    return summary_rows


# Export the improved function
build_country_year_summaries = build_country_year_summaries_with_ci


def build_mechanism_summaries(
    mechanism_rows: list[dict[str, str]],
    total_rows: int,
    interpretable_total: int,
) -> list[dict[str, str]]:
    all_countries = set()
    all_years = []
    for row in mechanism_rows:
        country = row.get("country_iso3", "")
        year_str = row.get("year", "")
        year = parse_int(year_str)
        if country:
            all_countries.add(country)
        if year is not None:
            all_years.append(year)
    
    return [
        {
            "prn_mechanism_call": normalize_text(mechanism_rows[0].get("prn_mechanism_call", "")),
            "sample_count": str(len(mechanism_rows)),
            "event_count": str(len(set(row.get("prn_event_id", "") for row in mechanism_rows))),
            "country_count": str(len(all_countries)),
            "year_min": str(min(all_years)) if all_years else "",
            "year_max": str(max(all_years)) if all_years else "",
            "is_interpretable": str(is_interpretable(mechanism_rows[0])),
            "is_definitive_disrupted": str(is_definitive_disrupted(mechanism_rows[0])),
            "is_uncertain_fragmented": str(is_uncertain_fragmented(mechanism_rows[0])),
            "sample_fraction_all": format_fraction(len(mechanism_rows), total_rows),
            "sample_fraction_within_mechanism": (
                format_fraction(len(mechanism_rows), interpretable_total)
                if is_interpretable(mechanism_rows[0])
                else ""
            ),
        }
    ]


def build_confidence_summaries(mechanism_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    totals_by_mechanism: Dict[str, int] = defaultdict(int)
    total_rows = len(mechanism_rows)
    for row in mechanism_rows:
        mechanism = normalize_text(row.get("prn_mechanism_call", ""))
        confidence = normalize_text(row.get("prn_call_confidence", ""))
        grouped[(mechanism, confidence)].append(row)
        totals_by_mechanism[mechanism] += 1
    
    summary_rows = []
    for (mechanism, confidence), rows in sorted(grouped.items()):
        all_countries = set()
        all_years = []
        for row in rows:
            country = row.get("country_iso3", "")
            year_str = row.get("year", "")
            year = parse_int(year_str)
            if country:
                all_countries.add(country)
            if year is not None:
                all_years.append(year)
        
        summary_rows.append(
            {
                "prn_mechanism_call": mechanism,
                "prn_call_confidence": confidence,
                "sample_count": str(len(rows)),
                "event_count": str(len(set(row.get("prn_event_id", "") for row in rows))),
                "country_count": str(len(all_countries)),
                "year_min": str(min(all_years)) if all_years else "",
                "year_max": str(max(all_years)) if all_years else "",
                "sample_fraction_all": format_fraction(len(rows), total_rows),
                "sample_fraction_within_mechanism": format_fraction(len(rows), totals_by_mechanism[mechanism]),
            }
        )
    return summary_rows


def build_breakpoint_summaries(mechanism_rows: list[dict[str, str]], is_hit_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    is_hits_by_sample = {row["assembly_accession"]: row for row in is_hit_rows}
    
    grouped: Dict[Tuple[str, str, str, str, str, str, str], List[dict]] = defaultdict(list)
    for row in mechanism_rows:
        mechanism = normalize_text(row.get("prn_mechanism_call", ""))
        confidence = normalize_text(row.get("prn_call_confidence", ""))
        bp_cat = normalize_text(row.get("bp_category", ""))
        ins_size = insertion_size_bin(row.get("insertion_subject_gap_bp", ""))
        ins_gap = row.get("insertion_subject_gap_bp", "")
        acc = row.get("assembly_accession", "")
        is_row = is_hits_by_sample.get(acc, {})
        best_hit = normalize_text(is_row.get("is_element_best_hit", ""))
        support = normalize_text(is_row.get("best_hit_support_tier", ""))
        grouped[(mechanism, confidence, bp_cat, ins_size, ins_gap, best_hit, support)].append(row)
    
    summary_rows = []
    for key, rows in sorted(grouped.items()):
        (mechanism, confidence, bp_cat, ins_size, ins_gap, best_hit, support) = key
        all_countries = set()
        all_years = []
        for row in rows:
            country = row.get("country_iso3", "")
            year_str = row.get("year", "")
            year = parse_int(year_str)
            if country:
                all_countries.add(country)
            if year is not None:
                all_years.append(year)
        
        summary_rows.append(
            {
                "prn_mechanism_call": mechanism,
                "prn_call_confidence": confidence,
                "bp_category": bp_cat,
                "insertion_size_bin": ins_size,
                "insertion_subject_gap_bp": ins_gap,
                "is_element_best_hit": best_hit,
                "best_hit_support_tier": support,
                "sample_count": str(len(rows)),
                "event_count": str(len(set(row.get("prn_event_id", "") for row in rows))),
                "country_count": str(len(all_countries)),
                "year_min": str(min(all_years)) if all_years else "",
                "year_max": str(max(all_years)) if all_years else "",
            }
        )
    return summary_rows


def build_unresolved_summaries(mechanism_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    uncertain_rows = [row for row in mechanism_rows if is_uncertain_fragmented(row)]
    
    grouped: Dict[Tuple[str, str, str, str], List[dict]] = defaultdict(list)
    for row in uncertain_rows:
        mechanism = normalize_text(row.get("prn_mechanism_call", ""))
        initial = normalize_text(row.get("prn_call_initial", ""))
        confidence = normalize_text(row.get("prn_call_confidence", ""))
        prn04 = normalize_text(row.get("prn04_rule_applied", ""))
        bp_cat = normalize_text(row.get("bp_category", ""))
        grouped[(mechanism, initial, confidence, prn04, bp_cat)].append(row)
    
    summary_rows = []
    for key, rows in sorted(grouped.items()):
        (mechanism, initial, confidence, prn04, bp_cat) = key
        all_countries = set()
        all_years = []
        example_event = ""
        example_acc = ""
        for row in rows:
            country = row.get("country_iso3", "")
            year_str = row.get("year", "")
            year = parse_int(year_str)
            if country:
                all_countries.add(country)
            if year is not None:
                all_years.append(year)
            if not example_event:
                example_event = row.get("prn_event_id", "")
                example_acc = row.get("assembly_accession", "")
        
        summary_rows.append(
            {
                "prn_mechanism_call": mechanism,
                "prn_call_initial": initial,
                "prn_call_confidence": confidence,
                "prn04_rule": prn04,
                "bp_category": bp_cat,
                "sample_count": str(len(rows)),
                "event_count": str(len(set(row.get("prn_event_id", "") for row in rows))),
                "country_count": str(len(all_countries)),
                "year_min": str(min(all_years)) if all_years else "",
                "year_max": str(max(all_years)) if all_years else "",
                "example_prn_event_id": example_event,
                "example_assembly_accession": example_acc,
            }
        )
    return summary_rows


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    
    mechanism_rows = load_tsv_rows(args.mechanism_calls)
    is_hit_rows = load_tsv_rows(args.is_hits)
    
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    mechanism_groups: Dict[str, List[dict]] = defaultdict(list)
    for row in mechanism_rows:
        mechanism = normalize_text(row.get("prn_mechanism_call", ""))
        mechanism_groups[mechanism].append(row)
    interpretable_total = sum(1 for row in mechanism_rows if is_interpretable(row))
    
    mechanism_summaries = []
    for mechanism, rows in sorted(mechanism_groups.items()):
        mechanism_summaries.extend(build_mechanism_summaries(rows, len(mechanism_rows), interpretable_total))
    write_tsv(output_dir / "bp_prn_mechanism_summary.tsv", MECHANISM_SUMMARY_COLUMNS, mechanism_summaries)
    
    confidence_summaries = build_confidence_summaries(mechanism_rows)
    write_tsv(output_dir / "bp_prn_confidence_summary.tsv", CONFIDENCE_SUMMARY_COLUMNS, confidence_summaries)
    
    breakpoint_summaries = build_breakpoint_summaries(mechanism_rows, is_hit_rows)
    write_tsv(output_dir / "bp_prn_breakpoint_summary.tsv", BREAKPOINT_SUMMARY_COLUMNS, breakpoint_summaries)
    
    unresolved_summaries = build_unresolved_summaries(mechanism_rows)
    write_tsv(output_dir / "bp_prn_unresolved_summary.tsv", UNRESOLVED_SUMMARY_COLUMNS, unresolved_summaries)
    
    # Use IMPROVED function with confidence intervals
    country_year_summaries = build_country_year_summaries_with_ci(mechanism_rows, is_hit_rows)
    write_tsv(output_dir / "bp_prn_country_year_summary.tsv", COUNTRY_YEAR_SUMMARY_COLUMNS, country_year_summaries)
    
    return 0


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Materialize PRN summary tables for mechanism composition, breakpoint patterns, "
            "unresolved cases, and country-year aggregation. "
            "IMPROVED: Now includes Wilson confidence intervals for disruption rates."
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
        "--output-dir",
        "-o",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs",
        help="Output directory for summary tables.",
    )
    return parser


if __name__ == "__main__":
    raise SystemExit(main())
