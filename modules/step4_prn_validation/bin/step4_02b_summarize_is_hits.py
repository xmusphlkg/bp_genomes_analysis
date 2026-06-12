#!/usr/bin/env python3
"""Summarize PRN-02 insertion-sequence scan evidence into a reviewable hit table."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from step4_02_scan_prn_mechanisms import (
    STEP3_DATA_ROOT,
    STEP4_DATA_ROOT,
    accession_root,
    choose_candidate,
    load_gap_flank_sequences,
    load_reference_fasta,
    load_tsv_rows,
    normalize_text,
    parse_float,
    recover_gap_sequence,
    scan_is_gap_sequence_all,
    write_tsv,
)


IS_HIT_COLUMNS = [
    "prn_event_id",
    "sample_id_canonical",
    "biosample_accession",
    "assembly_accession",
    "country_iso3",
    "year",
    "mlst_st",
    "prn_mechanism_call",
    "bp_category",
    "insertion_subject_gap_bp",
    "gap_sequence_source",
    "is_scan_status",
    "reference_id",
    "is_element_name",
    "reference_source_accession",
    "hit_rank",
    "is_best_hit",
    "passes_prn02_support_rule",
    "supports_assigned_mechanism",
    "hit_support_tier",
    "hit_orientation",
    "hit_pident",
    "hit_qcov",
    "hit_score",
    "event_hit_sample_n",
    "event_best_hit_sample_n",
    "notes",
]


def support_tier(hit_pident: float | None, hit_qcov: float | None) -> str:
    if hit_pident is None or hit_qcov is None:
        return "no_signal"
    if hit_qcov >= 80.0 and hit_pident >= 90.0:
        return "strong"
    if hit_qcov >= 60.0 and hit_pident >= 85.0:
        return "moderate"
    if hit_qcov > 0.0 and hit_pident > 0.0:
        return "weak"
    return "no_signal"


def hit_supports_assigned_mechanism(mechanism_call: str, is_element_name: str, passes_rule: bool) -> str:
    if not passes_rule:
        return "false"
    if mechanism_call == "coding_disrupted_is481":
        return "true" if is_element_name == "IS481" else "false"
    if mechanism_call == "coding_disrupted_other_is":
        return "true" if is_element_name and is_element_name != "IS481" else "false"
    if mechanism_call in {"coding_disrupted_other", "coding_disrupted_inversion_or_rearrangement"}:
        return "false"
    return ""


def build_scan_status(mechanism_row: dict[str, str], *, gap_source: str, hit_rank: int | None, passes_rule: bool) -> str:
    if normalize_text(mechanism_row.get("bp_category", "")) != "insertion_like":
        return "not_applicable_non_insertion_like"
    if gap_source in {"missing_gap_metadata", "gap_sequence_source_missing"}:
        return "not_scannable_missing_gap_sequence"
    if hit_rank is None:
        return "scanned_no_hit_rows_emitted"
    if hit_rank == 1 and passes_rule:
        return "scanned_supported_best_hit"
    if hit_rank == 1:
        return "scanned_candidate_best_hit"
    return "scanned_lower_rank_candidate"


def build_base_row(mechanism_row: dict[str, str], *, gap_source: str) -> dict[str, str]:
    return {
        "prn_event_id": mechanism_row.get("prn_event_id", ""),
        "sample_id_canonical": mechanism_row.get("sample_id_canonical", ""),
        "biosample_accession": mechanism_row.get("biosample_accession", ""),
        "assembly_accession": mechanism_row.get("assembly_accession", ""),
        "country_iso3": mechanism_row.get("country_iso3", ""),
        "year": mechanism_row.get("year", ""),
        "mlst_st": mechanism_row.get("mlst_st", ""),
        "prn_mechanism_call": mechanism_row.get("prn_mechanism_call", ""),
        "bp_category": mechanism_row.get("bp_category", ""),
        "insertion_subject_gap_bp": mechanism_row.get("insertion_subject_gap_bp", ""),
        "gap_sequence_source": gap_source,
        "is_scan_status": "",
        "reference_id": "",
        "is_element_name": "",
        "reference_source_accession": "",
        "hit_rank": "",
        "is_best_hit": "",
        "passes_prn02_support_rule": "",
        "supports_assigned_mechanism": "",
        "hit_support_tier": "",
        "hit_orientation": "",
        "hit_pident": "",
        "hit_qcov": "",
        "hit_score": "",
        "event_hit_sample_n": "",
        "event_best_hit_sample_n": "",
        "notes": "",
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize PRN-02 insertion-sequence evidence into a sample-level event-to-hit table "
            "that preserves hit ranks and best-hit status."
        )
    )
    parser.add_argument(
        "--mechanism-calls",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_mechanism_calls.tsv",
        help="PRN-02 mechanism calls TSV.",
    )
    parser.add_argument(
        "--gap-metadata",
        type=Path,
        default=STEP3_DATA_ROOT / "outputs" / "bp_prn_insertion_gap_plus_flanks.tsv",
        help="Step3 gap extraction metadata TSV.",
    )
    parser.add_argument(
        "--gap-flank-fasta",
        type=Path,
        default=STEP3_DATA_ROOT / "outputs" / "bp_prn_insertion_gap_plus_flanks.fasta",
        help="Step3 extracted gap+flank FASTA.",
    )
    parser.add_argument(
        "--is-reference-fasta",
        type=Path,
        default=STEP4_DATA_ROOT / "references" / "is_elements" / "bp_is_reference.fasta",
        help="PRN-01 reference FASTA.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_is_hits.tsv",
        help="Output event-to-IS-hit TSV.",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    mechanism_rows = load_tsv_rows(args.mechanism_calls)
    gap_rows = load_tsv_rows(args.gap_metadata)
    gap_flank_sequences = load_gap_flank_sequences(args.gap_flank_fasta)
    references = load_reference_fasta(args.is_reference_fasta)

    gap_by_root: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in gap_rows:
        gap_by_root[accession_root(row.get("genome_resolved_accession", ""))].append(row)

    hit_rows: list[dict[str, str]] = []
    event_hit_samples: dict[tuple[str, str], set[str]] = defaultdict(set)
    event_best_hit_samples: dict[tuple[str, str], set[str]] = defaultdict(set)

    for mechanism_row in mechanism_rows:
        assembly_accession = mechanism_row.get("assembly_accession", "")
        root = accession_root(assembly_accession)
        gap_row = choose_candidate(gap_by_root.get(root, []), assembly_accession, "genome_resolved_accession")
        gap_sequence, gap_source = recover_gap_sequence(gap_row, gap_flank_sequences)
        base = build_base_row(mechanism_row, gap_source=gap_source)

        if normalize_text(mechanism_row.get("bp_category", "")) != "insertion_like":
            row = dict(base)
            row["is_scan_status"] = build_scan_status(mechanism_row, gap_source=gap_source, hit_rank=None, passes_rule=False)
            row["notes"] = "is_screen_not_attempted_for_non_insertion_like_pattern"
            hit_rows.append(row)
            continue

        if not gap_sequence:
            row = dict(base)
            row["is_scan_status"] = build_scan_status(mechanism_row, gap_source=gap_source, hit_rank=None, passes_rule=False)
            row["notes"] = "unable_to_recover_gap_sequence_for_is_screen"
            hit_rows.append(row)
            continue

        hits = scan_is_gap_sequence_all(gap_sequence, references)
        if not hits:
            row = dict(base)
            row["is_scan_status"] = "scanned_no_hit_rows_emitted"
            row["notes"] = "gap_sequence_scanned_but_no_reference_hits_returned"
            hit_rows.append(row)
            continue

        for rank, hit in enumerate(hits, start=1):
            hit_pid = parse_float(hit["hit_pident"])
            hit_qcov = parse_float(hit["hit_qcov"])
            tier = support_tier(hit_pid, hit_qcov)
            passes_rule = rank == 1 and tier in {"strong", "moderate"}
            row = dict(base)
            row.update(
                {
                    "is_scan_status": build_scan_status(
                        mechanism_row,
                        gap_source=gap_source,
                        hit_rank=rank,
                        passes_rule=passes_rule,
                    ),
                    "reference_id": hit["reference_id"],
                    "is_element_name": hit["is_element_name"],
                    "reference_source_accession": hit["reference_source_accession"],
                    "hit_rank": str(rank),
                    "is_best_hit": "true" if rank == 1 else "false",
                    "passes_prn02_support_rule": "true" if passes_rule else "false",
                    "supports_assigned_mechanism": hit_supports_assigned_mechanism(
                        mechanism_row.get("prn_mechanism_call", ""),
                        hit["is_element_name"],
                        passes_rule,
                    ),
                    "hit_support_tier": tier,
                    "hit_orientation": hit["hit_orientation"],
                    "hit_pident": hit["hit_pident"],
                    "hit_qcov": hit["hit_qcov"],
                    "hit_score": hit["hit_score"],
                    "notes": "gap_sequence_scanned_against_prn01_reference_set",
                }
            )
            hit_rows.append(row)
            if tier in {"strong", "moderate"}:
                event_hit_samples[(row["prn_event_id"], row["is_element_name"])].add(row["sample_id_canonical"])
            if passes_rule:
                event_best_hit_samples[(row["prn_event_id"], row["is_element_name"])].add(row["sample_id_canonical"])

    for row in hit_rows:
        key = (row["prn_event_id"], row["is_element_name"])
        row["event_hit_sample_n"] = str(len(event_hit_samples.get(key, set()))) if row["is_element_name"] else "0"
        row["event_best_hit_sample_n"] = (
            str(len(event_best_hit_samples.get(key, set()))) if row["is_element_name"] else "0"
        )

    write_tsv(args.out, IS_HIT_COLUMNS, hit_rows)

    status_counts = Counter(row["is_scan_status"] for row in hit_rows)
    print(f"Wrote IS-hit summary: {args.out}")
    print(f"Rows written: {len(hit_rows)}")
    print("Scan status counts: " + ", ".join(f"{k}={v}" for k, v in sorted(status_counts.items())))
    print(
        "Supported best-hit rows: "
        f"{sum(1 for row in hit_rows if row['passes_prn02_support_rule'] == 'true')}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
