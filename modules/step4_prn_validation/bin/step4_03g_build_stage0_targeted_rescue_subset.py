#!/usr/bin/env python3
"""Build the Stage 0 targeted-rescue subset for retained legacy Step3 PRN gaps.

This is the plan-specific rescue manifest for the current manuscript revision.
It identifies retained `AUS`, `GBR`, and `JPN` samples that still carry the
legacy `no_current_step3_prn_input` flag, excludes rows already rescued through
curated overrides, and emits a subset that can be passed directly to the
existing Step4 batch builder and read-validation runner.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from step4_02_scan_prn_mechanisms import (
    STEP4_DATA_ROOT,
    WORKFLOW_DATA_ROOT,
    load_tsv_rows,
    normalize_text,
    write_tsv,
)


LEGACY_GAP_RULE = "prn04_rule=no_current_step3_prn_input"
DEFAULT_TARGET_COUNTRIES = ("AUS", "GBR", "JPN")

SUBSET_COLUMNS = [
    "sample_id_canonical",
    "biosample_accession",
    "assembly_accession",
    "country_name",
    "country_iso3",
    "year",
    "sra_run_accession",
    "ena_run_accession",
    "sra_sample_accession",
    "ena_sample_accession",
    "read_accession_primary",
    "read_accession_source",
    "raw_reads_available",
    "raw_read_run_count",
    "raw_read_link_status",
    "raw_read_link_source",
    "prn_call_initial",
    "prn_mechanism_call",
    "prn_call_confidence",
    "prn_event_id",
    "prn_interpretable",
    "prn_disrupted",
    "manifest_read_validation_status",
    "existing_validation_status",
    "existing_validation_method",
    "existing_validation_note",
    "evidence_flags",
    "notes",
    "record_decision",
    "phylogeny_selected_for_tree",
    "prn_rescue_status",
    "prn_rescue_source",
    "rescued_prn_call",
    "data_origin",
    "country_program_target",
    "culture_status",
    "specimen_type",
    "ct_or_dna_input",
    "stage0_track",
    "stage0_priority",
    "stage0_reason",
    "selection_note",
]

SUMMARY_COLUMNS = [
    "country_iso3",
    "stage0_track",
    "manifest_legacy_gap_rows",
    "retained_legacy_gap_rows",
    "linked_retained_rows",
    "rescued_override_rows",
    "existing_validation_rows",
    "pending_stage0_subset_rows",
    "nonlinked_retained_rows",
    "selection_note",
]


def parse_bool(value: str) -> bool:
    return normalize_text(value).lower() in {"true", "1", "yes", "y", "t"}


def parse_year(value: str) -> int | None:
    text = normalize_text(value)
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def split_countries(value: str) -> tuple[str, ...]:
    countries = []
    for token in value.split(","):
        cleaned = normalize_text(token).upper()
        if cleaned:
            countries.append(cleaned)
    return tuple(countries)


def read_accession_primary(row: dict[str, str]) -> tuple[str, str]:
    for field, source in [
        ("sra_run_accession", "SRA"),
        ("ena_run_accession", "ENA"),
        ("sra_sample_accession", "SRA_SAMPLE_ONLY"),
        ("ena_sample_accession", "ENA_SAMPLE_ONLY"),
    ]:
        token = normalize_text(row.get(field, ""))
        if token:
            return token, source
    return "", "unavailable"


def is_retained(row: dict[str, str]) -> bool:
    decision = normalize_text(row.get("record_decision", ""))
    return decision.startswith("retain_") or parse_bool(row.get("phylogeny_selected_for_tree", ""))


def is_legacy_gap(row: dict[str, str]) -> bool:
    notes = normalize_text(row.get("notes", ""))
    evidence_flags = normalize_text(row.get("evidence_flags", ""))
    return LEGACY_GAP_RULE in notes or "no_step3_prn_input" in evidence_flags


def stage0_track(country_iso3: str, year: int | None) -> str:
    if country_iso3 == "GBR":
        return "gbr_fast_wp_anchor_rescue"
    if country_iso3 == "JPN":
        if year is not None and year >= 2012:
            return "jpn_post2012_prn_free_anchor"
        return "jpn_pre2012_prn_free_anchor_rescue"
    if country_iso3 == "AUS":
        if year is not None and 1997 <= year <= 2014:
            return "aus_transition_window_densification"
        return "aus_legacy_gap_cleanup"
    return "non_target"


def stage0_priority(country_iso3: str, year: int | None) -> str:
    if country_iso3 == "GBR":
        return "high"
    if country_iso3 == "JPN":
        return "high" if year is not None and year < 2012 else "medium"
    if country_iso3 == "AUS":
        return "high" if year is not None and year >= 1997 else "medium"
    return "low"


def stage0_reason(country_iso3: str, year: int | None) -> str:
    if country_iso3 == "GBR":
        return "recover_historical_wp_anchor_before_new_uk_collection"
    if country_iso3 == "JPN":
        if year is not None and year >= 2012:
            return "triage_post2012_prn_free_context_for_japan_anchor"
        return "rescue_pre2012_comparator_before_new_post2012_acquisition"
    if country_iso3 == "AUS":
        if year is not None and 1997 <= year <= 2014:
            return "densify_transition_and_early_ap_window_before_interpretation"
        return "clear_legacy_gap_outside_primary_transition_window"
    return "outside_current_stage0_targets"


def load_validation_index(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}

    index: dict[str, dict[str, str]] = {}
    for row in load_tsv_rows(path):
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        if sample_id and sample_id not in index:
            index[sample_id] = row
    return index


def build_subset_row(row: dict[str, str], validation_row: dict[str, str] | None) -> dict[str, str]:
    year = parse_year(row.get("year", ""))
    country_iso3 = normalize_text(row.get("country_iso3", "")).upper()
    primary_accession, accession_source = read_accession_primary(row)
    validation_status = normalize_text((validation_row or {}).get("read_validation_status", ""))
    validation_method = normalize_text((validation_row or {}).get("validation_method", ""))
    validation_note = normalize_text((validation_row or {}).get("notes", ""))

    return {
        "sample_id_canonical": normalize_text(row.get("sample_id_canonical", "")),
        "biosample_accession": normalize_text(row.get("biosample_accession", "")),
        "assembly_accession": normalize_text(row.get("assembly_accession", "")),
        "country_name": normalize_text(row.get("country", "")) or country_iso3,
        "country_iso3": country_iso3,
        "year": normalize_text(row.get("year", "")),
        "sra_run_accession": normalize_text(row.get("sra_run_accession", "")),
        "ena_run_accession": normalize_text(row.get("ena_run_accession", "")),
        "sra_sample_accession": normalize_text(row.get("sra_sample_accession", "")),
        "ena_sample_accession": normalize_text(row.get("ena_sample_accession", "")),
        "read_accession_primary": primary_accession,
        "read_accession_source": accession_source,
        "raw_reads_available": normalize_text(row.get("raw_reads_available", "")),
        "raw_read_run_count": normalize_text(row.get("raw_read_run_count", "")),
        "raw_read_link_status": normalize_text(row.get("raw_read_link_status", "")),
        "raw_read_link_source": normalize_text(row.get("raw_read_link_source", "")),
        "prn_call_initial": normalize_text(row.get("prn_call_initial", "")),
        "prn_mechanism_call": normalize_text(row.get("prn_mechanism_call", "")),
        "prn_call_confidence": normalize_text(row.get("prn_call_confidence", "")),
        "prn_event_id": normalize_text(row.get("prn_event_id", "")),
        "prn_interpretable": normalize_text(row.get("prn_interpretable", "")),
        "prn_disrupted": normalize_text(row.get("prn_disrupted", "")),
        "manifest_read_validation_status": normalize_text(row.get("read_validation_status", "")),
        "existing_validation_status": validation_status,
        "existing_validation_method": validation_method,
        "existing_validation_note": validation_note,
        "evidence_flags": normalize_text(row.get("evidence_flags", "")),
        "notes": normalize_text(row.get("notes", "")),
        "record_decision": normalize_text(row.get("record_decision", "")),
        "phylogeny_selected_for_tree": normalize_text(row.get("phylogeny_selected_for_tree", "")),
        "prn_rescue_status": normalize_text(row.get("prn_rescue_status", "")),
        "prn_rescue_source": normalize_text(row.get("prn_rescue_source", "")),
        "rescued_prn_call": normalize_text(row.get("rescued_prn_call", "")),
        "data_origin": normalize_text(row.get("data_origin", "")),
        "country_program_target": normalize_text(row.get("country_program_target", "")),
        "culture_status": normalize_text(row.get("culture_status", "")),
        "specimen_type": normalize_text(row.get("specimen_type", "")),
        "ct_or_dna_input": normalize_text(row.get("ct_or_dna_input", "")),
        "stage0_track": stage0_track(country_iso3, year),
        "stage0_priority": stage0_priority(country_iso3, year),
        "stage0_reason": stage0_reason(country_iso3, year),
        "selection_note": "retained_legacy_gap_linked_without_promoted_rescue_override",
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build the Stage 0 targeted rescue subset for retained AUS/GBR/JPN legacy Step3 PRN gaps."
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=WORKFLOW_DATA_ROOT / "manifest" / "manifest.tsv",
        help="Canonical manifest TSV.",
    )
    parser.add_argument(
        "--read-validation",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_read_validation.tsv",
        help="Existing Step4 read-validation TSV used for annotation only.",
    )
    parser.add_argument(
        "--target-countries",
        default=",".join(DEFAULT_TARGET_COUNTRIES),
        help="Comma-separated ISO3 list for Stage 0 rescue targeting.",
    )
    parser.add_argument(
        "--include-rescued",
        action="store_true",
        help="Include already rescued override rows in the subset output.",
    )
    parser.add_argument(
        "--out-subset",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_stage0_targeted_rescue_subset.tsv",
        help="Output subset TSV for Stage 0 targeted rescue.",
    )
    parser.add_argument(
        "--out-summary",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_stage0_targeted_rescue_summary.tsv",
        help="Output per-country rescue summary TSV.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    target_countries = set(split_countries(args.target_countries))
    manifest_rows = load_tsv_rows(args.manifest)
    validation_index = load_validation_index(args.read_validation)

    counters: dict[str, Counter] = {}
    subset_rows: list[dict[str, str]] = []

    for row in manifest_rows:
        country_iso3 = normalize_text(row.get("country_iso3", "")).upper()
        if country_iso3 not in target_countries:
            continue
        if not is_legacy_gap(row):
            continue

        country_counter = counters.setdefault(country_iso3, Counter())
        country_counter["manifest_legacy_gap_rows"] += 1

        retained = is_retained(row)
        linked = normalize_text(row.get("raw_read_link_status", "")) == "linked"
        rescued = normalize_text(row.get("prn_rescue_status", "")) == "rescued_override" or normalize_text(
            row.get("rescued_prn_call", "")
        ) != ""
        has_existing_validation = normalize_text(
            (validation_index.get(normalize_text(row.get("sample_id_canonical", "")), {}) or {}).get(
                "read_validation_status", ""
            )
        ) != ""

        if retained:
            country_counter["retained_legacy_gap_rows"] += 1
        if retained and linked:
            country_counter["linked_retained_rows"] += 1
        if retained and not linked:
            country_counter["nonlinked_retained_rows"] += 1
        if rescued:
            country_counter["rescued_override_rows"] += 1
        if has_existing_validation:
            country_counter["existing_validation_rows"] += 1

        if not retained or not linked:
            continue
        if rescued and not args.include_rescued:
            continue

        subset_rows.append(
            build_subset_row(
                row=row,
                validation_row=validation_index.get(normalize_text(row.get("sample_id_canonical", ""))),
            )
        )
        country_counter["pending_stage0_subset_rows"] += 1

    subset_rows.sort(
        key=lambda row: (
            row["country_iso3"],
            {"high": 0, "medium": 1, "low": 2}.get(row["stage0_priority"], 9),
            parse_year(row["year"]) if parse_year(row["year"]) is not None else 999999,
            row["assembly_accession"],
        )
    )

    summary_rows: list[dict[str, str]] = []
    overall = Counter()
    for country_iso3 in sorted(counters):
        counter = counters[country_iso3]
        overall.update(counter)
        summary_rows.append(
            {
                "country_iso3": country_iso3,
                "stage0_track": stage0_track(country_iso3, None),
                "manifest_legacy_gap_rows": str(counter.get("manifest_legacy_gap_rows", 0)),
                "retained_legacy_gap_rows": str(counter.get("retained_legacy_gap_rows", 0)),
                "linked_retained_rows": str(counter.get("linked_retained_rows", 0)),
                "rescued_override_rows": str(counter.get("rescued_override_rows", 0)),
                "existing_validation_rows": str(counter.get("existing_validation_rows", 0)),
                "pending_stage0_subset_rows": str(counter.get("pending_stage0_subset_rows", 0)),
                "nonlinked_retained_rows": str(counter.get("nonlinked_retained_rows", 0)),
                "selection_note": "summary_excludes_non_target_countries_and_non_legacy_gap_rows",
            }
        )

    if counters:
        summary_rows.append(
            {
                "country_iso3": "ALL",
                "stage0_track": "combined_stage0_targeted_rescue",
                "manifest_legacy_gap_rows": str(overall.get("manifest_legacy_gap_rows", 0)),
                "retained_legacy_gap_rows": str(overall.get("retained_legacy_gap_rows", 0)),
                "linked_retained_rows": str(overall.get("linked_retained_rows", 0)),
                "rescued_override_rows": str(overall.get("rescued_override_rows", 0)),
                "existing_validation_rows": str(overall.get("existing_validation_rows", 0)),
                "pending_stage0_subset_rows": str(overall.get("pending_stage0_subset_rows", 0)),
                "nonlinked_retained_rows": str(overall.get("nonlinked_retained_rows", 0)),
                "selection_note": "combined_across_target_countries",
            }
        )

    write_tsv(args.out_subset, SUBSET_COLUMNS, subset_rows)
    write_tsv(args.out_summary, SUMMARY_COLUMNS, summary_rows)

    print(f"Stage 0 targeted rescue subset rows: {len(subset_rows)}")
    print(f"Subset written to: {args.out_subset}")
    print(f"Summary written to: {args.out_summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
