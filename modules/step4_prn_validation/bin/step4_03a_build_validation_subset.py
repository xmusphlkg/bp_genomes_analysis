#!/usr/bin/env python3
"""Build a balanced raw-read validation subset manifest for PRN validation."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from step4_02_scan_prn_mechanisms import (
    STEP1_DATA_ROOT,
    STEP4_DATA_ROOT,
    load_tsv_rows,
    normalize_text,
    parse_int,
    write_tsv,
)


OUTPUT_COLUMNS = [
    "sample_id_canonical",
    "biosample_accession",
    "assembly_accession",
    "current_accession",
    "country_name",
    "country_iso3",
    "year",
    "country_year_key",
    "analysis_cohort_id",
    "analysis_cohort_name",
    "cohort_priority_flag",
    "country_year_cell_genome_n",
    "raw_reads_available",
    "raw_read_run_count",
    "raw_read_link_status",
    "raw_read_link_source",
    "raw_read_lookup_date",
    "sra_run_accession",
    "ena_run_accession",
    "sra_sample_accession",
    "ena_sample_accession",
    "read_accession_primary",
    "read_accession_source",
    "prn_call_initial",
    "prn_mechanism_call",
    "prn_call_confidence",
    "prn_event_id",
    "prn_query_cov_pct",
    "prn_best_single_cov_pct",
    "prn_hsp_n",
    "bp_category",
    "is_element_best_hit",
    "is_element_best_hit_pident",
    "is_element_best_hit_qcov",
    "mlst_st",
    "mlst_present",
    "selection_stratum",
    "selection_reason",
    "country_balance_key",
    "country_stratum_rank",
    "country_stratum_total",
    "selection_round",
    "selection_order",
    "stratum_available_n",
    "stratum_selected_n",
    "selection_notes",
    "evidence_flags",
    "notes",
]

STRATUM_QUOTAS = {
    "insufficient_partial": 1,
    "insufficient_missing_fasta": 5,
    "insufficient_no_step3": 8,
    "other_disruption_ambiguous": 19,
    "is481_moderate_boundary": 1,
    "rearrangement_structural": 6,
    "is481_strong_control": 8,
    "intact_control": 9,
}

STRATUM_REASONS = {
    "insufficient_partial": "partial_prn_call_unresolved",
    "insufficient_missing_fasta": "assembly_sequence_missing",
    "insufficient_no_step3": "no_step3_prn_input_unresolved",
    "other_disruption_ambiguous": "low_confidence_disrupted_other",
    "is481_moderate_boundary": "boundary_is481_case",
    "rearrangement_structural": "within_contig_structural_rearrangement",
    "is481_strong_control": "strong_is481_positive_control",
    "intact_control": "intact_control_for_validation",
}

COUNTRY_SORT_SENTINEL = "unknown"
YEAR_SORT_SENTINEL = 999999
MLST_MISSING = {"", "-", "na", "n/a", "none", "missing"}
COUNTRY_MISSING = {"", "-", "na", "n/a", "none", "missing", "unknown"}


def accession_source(row: dict[str, str]) -> tuple[str, str, str]:
    sra_run = normalize_text(row.get("sra_run_accession", ""))
    ena_run = normalize_text(row.get("ena_run_accession", ""))
    sra_sample = normalize_text(row.get("sra_sample_accession", ""))
    ena_sample = normalize_text(row.get("ena_sample_accession", ""))
    if sra_run:
        return sra_run, "SRA", sra_sample or sra_run
    if ena_run:
        return ena_run, "ENA", ena_sample or ena_run
    if sra_sample:
        return sra_sample, "SRA_SAMPLE_ONLY", sra_sample
    if ena_sample:
        return ena_sample, "ENA_SAMPLE_ONLY", ena_sample
    return "", "unavailable", ""


def is_present(value: str) -> bool:
    return normalize_text(value).casefold() not in MLST_MISSING


def country_balance_key(row: dict[str, str]) -> str:
    country_iso3 = normalize_text(row.get("country_iso3", ""))
    country_name = normalize_text(row.get("country_name", ""))
    if country_iso3 and country_iso3.casefold() not in COUNTRY_MISSING:
        return country_iso3
    if country_name and country_name.casefold() not in COUNTRY_MISSING:
        return country_name
    return COUNTRY_SORT_SENTINEL


def year_sort_key(value: str) -> int:
    parsed = parse_int(value)
    if parsed is None:
        return YEAR_SORT_SENTINEL
    return parsed


def mlst_sort_key(value: str) -> str:
    value = normalize_text(value)
    if not is_present(value):
        return ""
    return value


def build_records(cohort_rows: list[dict[str, str]], mechanism_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    mechanism_by_sample = {row["sample_id_canonical"]: row for row in mechanism_rows}
    records: list[dict[str, str]] = []

    for cohort_row in cohort_rows:
        sample_id = normalize_text(cohort_row.get("sample_id_canonical", ""))
        mechanism_row = mechanism_by_sample.get(sample_id)
        if mechanism_row is None:
            continue

        read_primary, read_source, sample_primary = accession_source(cohort_row)
        country_name = normalize_text(cohort_row.get("country", ""))
        country_iso3 = normalize_text(mechanism_row.get("country_iso3", ""))
        if not country_iso3 and country_name:
            country_iso3 = country_name

        mlst = normalize_text(mechanism_row.get("mlst_st", ""))
        mlst_clean = "" if not is_present(mlst) else mlst
        selection_stratum = build_selection_stratum(mechanism_row)

        records.append(
            {
                "sample_id_canonical": sample_id,
                "biosample_accession": normalize_text(cohort_row.get("biosample_accession", "")),
                "assembly_accession": normalize_text(cohort_row.get("assembly_accession", "")),
                "current_accession": normalize_text(cohort_row.get("current_accession", "")),
                "country_name": country_name,
                "country_iso3": normalize_text(mechanism_row.get("country_iso3", "")),
                "year": normalize_text(cohort_row.get("year", "")),
                "country_year_key": normalize_text(cohort_row.get("country_year_key", "")),
                "analysis_cohort_id": normalize_text(cohort_row.get("analysis_cohort_id", "")),
                "analysis_cohort_name": normalize_text(cohort_row.get("analysis_cohort_name", "")),
                "cohort_priority_flag": normalize_text(cohort_row.get("cohort_priority_flag", "")),
                "country_year_cell_genome_n": normalize_text(cohort_row.get("country_year_cell_genome_n", "")),
                "raw_reads_available": normalize_text(cohort_row.get("raw_reads_available", "")),
                "raw_read_run_count": normalize_text(cohort_row.get("raw_read_run_count", "")),
                "raw_read_link_status": normalize_text(cohort_row.get("raw_read_link_status", "")),
                "raw_read_link_source": normalize_text(cohort_row.get("raw_read_link_source", "")),
                "raw_read_lookup_date": normalize_text(cohort_row.get("raw_read_lookup_date", "")),
                "sra_run_accession": normalize_text(cohort_row.get("sra_run_accession", "")),
                "ena_run_accession": normalize_text(cohort_row.get("ena_run_accession", "")),
                "sra_sample_accession": normalize_text(cohort_row.get("sra_sample_accession", "")),
                "ena_sample_accession": normalize_text(cohort_row.get("ena_sample_accession", "")),
                "read_accession_primary": read_primary,
                "read_accession_source": read_source,
                "prn_call_initial": normalize_text(mechanism_row.get("prn_call_initial", "")),
                "prn_mechanism_call": normalize_text(mechanism_row.get("prn_mechanism_call", "")),
                "prn_call_confidence": normalize_text(mechanism_row.get("prn_call_confidence", "")),
                "prn_event_id": normalize_text(mechanism_row.get("prn_event_id", "")),
                "prn_query_cov_pct": normalize_text(mechanism_row.get("prn_query_cov_pct", "")),
                "prn_best_single_cov_pct": normalize_text(mechanism_row.get("prn_best_single_cov_pct", "")),
                "prn_hsp_n": normalize_text(mechanism_row.get("prn_hsp_n", "")),
                "bp_category": normalize_text(mechanism_row.get("bp_category", "")),
                "is_element_best_hit": normalize_text(mechanism_row.get("is_element_best_hit", "")),
                "is_element_best_hit_pident": normalize_text(mechanism_row.get("is_element_best_hit_pident", "")),
                "is_element_best_hit_qcov": normalize_text(mechanism_row.get("is_element_best_hit_qcov", "")),
                "mlst_st": mlst_clean,
                "mlst_present": "true" if mlst_clean else "false",
                "selection_stratum": selection_stratum,
                "selection_reason": STRATUM_REASONS[selection_stratum],
                "country_balance_key": country_balance_key({**cohort_row, **mechanism_row, "country_name": country_name}),
                "country_stratum_rank": "",
                "country_stratum_total": "",
                "selection_round": "",
                "selection_order": "",
                "stratum_available_n": "",
                "stratum_selected_n": "",
                "selection_notes": "country_round_robin;mlst_rare_first_if_present;year_ascending;run_provenance_explicit",
                "evidence_flags": normalize_text(mechanism_row.get("evidence_flags", "")),
                "notes": normalize_text(mechanism_row.get("notes", "")),
            }
        )

    return records


def build_selection_stratum(row: dict[str, str]) -> str:
    mechanism = normalize_text(row.get("prn_mechanism_call", ""))
    confidence = normalize_text(row.get("prn_call_confidence", ""))
    initial = normalize_text(row.get("prn_call_initial", ""))
    if mechanism == "intact":
        return "intact_control"
    if mechanism == "coding_disrupted_is481":
        if confidence == "assembly_moderate":
            return "is481_moderate_boundary"
        return "is481_strong_control"
    if mechanism == "coding_disrupted_inversion_or_rearrangement":
        return "rearrangement_structural"
    if mechanism == "coding_disrupted_other":
        return "other_disruption_ambiguous"
    if mechanism == "insufficient_data":
        if initial == "missing_fasta":
            return "insufficient_missing_fasta"
        if initial == "partial":
            return "insufficient_partial"
        return "insufficient_no_step3"
    return "ignored"


def sort_key_for_country_group(row: dict[str, str]) -> tuple:
    mlst = normalize_text(row.get("mlst_st", ""))
    mlst_present = is_present(mlst)
    mlst_key = mlst_sort_key(mlst)
    year_key = year_sort_key(row.get("year", ""))
    sample_id = normalize_text(row.get("sample_id_canonical", ""))
    run_source = normalize_text(row.get("read_accession_source", ""))
    return (
        0 if mlst_present else 1,
        mlst_key if mlst_present else "zzzzzz",
        year_key,
        0 if run_source == "ENA" else 1,
        sample_id,
    )


def round_robin_select(records: list[dict[str, str]], quota: int) -> list[dict[str, str]]:
    country_groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in records:
        country_groups[row["country_balance_key"]].append(row)

    country_order = sorted(country_groups, key=lambda c: (c == COUNTRY_SORT_SENTINEL, -len(country_groups[c]), c))
    for country in country_groups:
        country_groups[country].sort(key=sort_key_for_country_group)

    country_offsets = {country: 0 for country in country_groups}
    selected: list[dict[str, str]] = []
    country_selected_counts: dict[str, int] = defaultdict(int)
    country_available = {country: len(rows) for country, rows in country_groups.items()}
    round_index = 0

    while len(selected) < quota:
        progressed = False
        for country in country_order:
            idx = country_offsets[country]
            rows = country_groups[country]
            if idx >= len(rows):
                continue
            row = dict(rows[idx])
            country_offsets[country] += 1
            country_selected_counts[country] += 1
            row["country_stratum_rank"] = str(country_selected_counts[country])
            row["country_stratum_total"] = str(country_available[country])
            row["selection_order"] = str(len(selected) + 1)
            row["selection_round"] = str(round_index + 1)
            selected.append(row)
            progressed = True
            if len(selected) >= quota:
                break
        if not progressed:
            break
        round_index += 1

    return selected


def select_subset(records: list[dict[str, str]]) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    global_order = 1
    for stratum in [
        "insufficient_partial",
        "insufficient_missing_fasta",
        "insufficient_no_step3",
        "other_disruption_ambiguous",
        "is481_moderate_boundary",
        "rearrangement_structural",
        "is481_strong_control",
        "intact_control",
    ]:
        stratum_records = [row for row in records if row["selection_stratum"] == stratum]
        quota = STRATUM_QUOTAS[stratum]
        if not stratum_records:
            continue
        chosen = round_robin_select(stratum_records, min(quota, len(stratum_records)))
        for row in chosen:
            row["stratum_available_n"] = str(len(stratum_records))
            row["stratum_selected_n"] = str(len(chosen))
            row["selection_order"] = str(global_order)
            global_order += 1
        selected.extend(chosen)
    return selected


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build a balanced raw-read validation subset spanning major PRN mechanism classes, "
            "countries, confidence levels, and accessible read provenance."
        )
    )
    parser.add_argument(
        "--cohort-d",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_cohort_D_validation.tsv",
        help="Raw-read validation cohort manifest from GC-04.",
    )
    parser.add_argument(
        "--mechanism-calls",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_mechanism_calls.tsv",
        help="PRN mechanism calls from PRN-04.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_validation_subset.tsv",
        help="Validation subset TSV to write.",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    cohort_rows = load_tsv_rows(args.cohort_d)
    mechanism_rows = load_tsv_rows(args.mechanism_calls)

    records = build_records(cohort_rows, mechanism_rows)
    usable_records = [row for row in records if row["selection_stratum"] != "ignored" and row["raw_reads_available"] == "true"]
    selected = select_subset(usable_records)

    write_tsv(args.out, OUTPUT_COLUMNS, selected)

    stratum_counts = Counter(row["selection_stratum"] for row in selected)
    country_counts = Counter(row["country_balance_key"] for row in selected)
    print(f"selected_rows\t{len(selected)}")
    for stratum, count in sorted(stratum_counts.items()):
        print(f"stratum\t{stratum}\t{count}")
    print(f"country_diversity\t{len(country_counts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
