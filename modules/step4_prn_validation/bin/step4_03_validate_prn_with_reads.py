#!/usr/bin/env python3
"""Build a machine-readable PRN read-validation table for the validation subset."""

from __future__ import annotations

import argparse
import csv
import os
import platform
import re
from pathlib import Path

from step4_02_scan_prn_mechanisms import (
    load_tsv_rows,
    normalize_text,
    parse_float,
    parse_int,
    repo_root,
    write_tsv,
)


PROJECT_ROOT = Path(
    os.environ.get(
        "PERTUSSIS_PROJECT_DATA_ROOT",
        str(repo_root() / "pertussis_data" / "pertussis_gene"),
    )
)


OUTPUT_COLUMNS = [
    "sample_id_canonical",
    "sra_run_accession",
    "prn_event_id",
    "prn_mechanism_call",
    "read_validation_status",
    "read_support_class",
    "n_supporting_reads",
    "n_contradicting_reads",
    "junction_supported",
    "targeted_locus_assembly_status",
    "validation_method",
    "validator_version",
    "notes",
]

EVIDENCE_COLUMNS = [
    "sample_id_canonical",
    "prn_event_id",
    "prn_mechanism_call",
    "tool",
    "is_element_name",
    "query_reference_id",
    "reference_record",
    "locus_start",
    "locus_end",
    "prn_overlap_bp",
    "call_label",
    "orientation",
    "gap_bp",
    "percent_id",
    "percent_cov",
    "gene_interruption",
    "left_gene",
    "left_description",
    "right_gene",
    "right_description",
    "start_clipped_reads",
    "end_clipped_reads",
    "total_clipped_reads",
    "direct_repeats",
    "inverted_repeats",
    "left_sequence",
    "right_sequence",
    "evidence_note",
]

TSD_COLUMNS = [
    "sample_id_canonical",
    "prn_event_id",
    "prn_mechanism_call",
    "inferred_is_element_name",
    "chromosome",
    "insertion_start",
    "insertion_end",
    "total_clipped_reads",
    "direct_repeats",
    "inverted_repeats",
    "left_sequence",
    "right_sequence",
]

AUDIT_VALIDATION_METHOD = "accession_linked_read_validation_manifest_audit_only"
READ_VALIDATION_METHOD = "ismapper_panisa_stage4_prn_validation"
VALIDATOR_VERSION = "step4_03_validate_prn_with_reads.py::2026-04-05"


def overlaps(start: int, end: int, locus_start: int, locus_end: int) -> int:
    return max(0, min(end, locus_end) - max(start, locus_start) + 1)


def infer_is_element(value: str) -> str:
    value = normalize_text(value)
    upper = value.upper()
    for token in ("IS481", "IS1002", "IS1663"):
        if token in upper:
            return token
    return value.split("_", 1)[0] if value else "unknown"


def infer_panisa_element(row: dict[str, str]) -> str:
    signature = " ".join(
        [
            normalize_text(row.get("Inverted repeats", "")),
            normalize_text(row.get("Left sequence", "")),
            normalize_text(row.get("Right sequence", "")),
        ]
    ).upper()
    if "TGTGAA" in signature or "TTCACA" in signature:
        return "IS481"
    return "untyped"


def targeted_locus_assembly_status(row: dict[str, str]) -> str:
    initial = normalize_text(row.get("prn_call_initial", ""))
    mechanism = normalize_text(row.get("prn_mechanism_call", ""))
    if initial == "missing_fasta":
        return "assembly_sequence_missing"
    if initial == "partial":
        return "partial_prn_alignment"
    if initial == "not_available_current_step3":
        return "no_current_step3_prn_input"
    if mechanism == "intact":
        return "assembly_locus_intact"
    if mechanism == "insufficient_data":
        return "assembly_locus_unresolved"
    return "assembly_event_called"


def base_notes(row: dict[str, str]) -> list[str]:
    parts = [
        f"raw_read_link_status={normalize_text(row.get('raw_read_link_status', '')) or 'unknown'}",
        f"read_accession_source={normalize_text(row.get('read_accession_source', '')) or 'unknown'}",
        f"prn_call_confidence={normalize_text(row.get('prn_call_confidence', '')) or 'unknown'}",
    ]
    evidence_flags = normalize_text(row.get("evidence_flags", ""))
    if evidence_flags:
        parts.append(f"evidence_flags={evidence_flags}")
    upstream_notes = normalize_text(row.get("notes", ""))
    if upstream_notes:
        parts.append(f"upstream_notes={upstream_notes}")
    return parts


def load_prn_locus(gbff_path: Path) -> dict[str, str | int]:
    record_id = ""
    current_feature: dict[str, str] | None = None

    def finalize_feature(feature: dict[str, str] | None) -> dict[str, str | int] | None:
        if feature is None:
            return None
        qualifiers = "\n".join(feature["qualifiers"])
        if '/gene="prn"' not in qualifiers and "pertactin autotransporter" not in qualifiers.casefold():
            return None
        coordinates = [int(token) for token in re.findall(r"\d+", feature["location"])]
        if not coordinates:
            return None
        return {
            "record_id": record_id,
            "start": min(coordinates),
            "end": max(coordinates),
        }

    with gbff_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("VERSION"):
                parts = line.split()
                if len(parts) > 1:
                    record_id = parts[1]
                continue

            if re.match(r"^     \S", line) and not line.startswith("                     /"):
                parsed = finalize_feature(current_feature)
                if parsed is not None:
                    return parsed
                feature_type = line[5:21].strip()
                location = line[21:].strip()
                current_feature = {"type": feature_type, "location": location, "qualifiers": []}
                continue

            if current_feature is not None and line.startswith("                     /"):
                current_feature["qualifiers"].append(line.strip())

    parsed = finalize_feature(current_feature)
    if parsed is not None:
        return parsed
    raise ValueError(f"prn locus not found in {gbff_path}")


def load_batch_samples(batch_path: Path | None, work_root: Path | None) -> set[str]:
    if batch_path is not None and batch_path.exists():
        return {
            normalize_text(row.get("sample_id_canonical", ""))
            for row in load_tsv_rows(batch_path)
            if normalize_text(row.get("batch_status", "")) == "selected"
        }
    if work_root is None or not work_root.exists():
        return set()

    batch_samples: set[str] = set()
    ismapper_root = work_root / "ismapper"
    if ismapper_root.exists():
        for sample_dir in ismapper_root.iterdir():
            if sample_dir.is_dir():
                batch_samples.add(sample_dir.name)
    panisa_root = work_root / "panisa"
    if panisa_root.exists():
        for panisa_file in panisa_root.glob("*.panisa.tsv"):
            batch_samples.add(panisa_file.name.replace(".panisa.tsv", ""))
    return batch_samples


def merge_unique_sample_rows(
    existing_rows: list[dict[str, str]],
    refreshed_rows: list[dict[str, str]],
    *,
    sample_field: str = "sample_id_canonical",
) -> list[dict[str, str]]:
    refreshed_by_sample = {
        normalize_text(row.get(sample_field, "")): row
        for row in refreshed_rows
        if normalize_text(row.get(sample_field, ""))
    }
    merged: list[dict[str, str]] = []
    consumed: set[str] = set()

    for row in existing_rows:
        sample_id = normalize_text(row.get(sample_field, ""))
        replacement = refreshed_by_sample.get(sample_id)
        if replacement is None:
            merged.append(row)
            continue
        merged.append(replacement)
        consumed.add(sample_id)

    for row in refreshed_rows:
        sample_id = normalize_text(row.get(sample_field, ""))
        if sample_id in consumed:
            continue
        merged.append(row)

    return merged


def merge_multirow_sample_rows(
    existing_rows: list[dict[str, str]],
    refreshed_rows: list[dict[str, str]],
    *,
    sample_field: str = "sample_id_canonical",
) -> list[dict[str, str]]:
    refreshed_samples = {
        normalize_text(row.get(sample_field, ""))
        for row in refreshed_rows
        if normalize_text(row.get(sample_field, ""))
    }
    preserved = [
        row
        for row in existing_rows
        if normalize_text(row.get(sample_field, "")) not in refreshed_samples
    ]
    return preserved + refreshed_rows


def parse_ismapper_hits(
    sample_id: str,
    sample_row: dict[str, str],
    work_root: Path,
    prn_locus: dict[str, str | int],
) -> tuple[list[dict[str, str]], dict[str, int]]:
    evidence_rows: list[dict[str, str]] = []
    table_files = sorted((work_root / "ismapper" / sample_id).glob("**/*_table.txt"))
    metrics = {"ismapper_tables_found": len(table_files), "panisa_file_found": 0}
    prn_start = int(prn_locus["start"])
    prn_end = int(prn_locus["end"])
    prn_record = normalize_text(str(prn_locus["record_id"]))

    for table_path in table_files:
        query_reference_id = table_path.parent.name
        is_element_name = infer_is_element(query_reference_id)
        with table_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                x = parse_int(row.get("x", ""))
                y = parse_int(row.get("y", ""))
                if x is None or y is None:
                    continue
                locus_start = min(x, y)
                locus_end = max(x, y)
                left_description = normalize_text(row.get("left_description", ""))
                right_description = normalize_text(row.get("right_description", ""))
                prn_overlap_bp = overlaps(locus_start, locus_end, prn_start, prn_end)
                prn_named = "pertactin" in left_description.casefold() or "pertactin" in right_description.casefold()
                if prn_overlap_bp <= 0 and not prn_named:
                    continue
                reference_record = normalize_text(table_path.name.split("__", 1)[-1].removesuffix("_table.txt"))
                if reference_record and prn_record and reference_record != prn_record:
                    continue
                evidence_rows.append(
                    {
                        "sample_id_canonical": sample_id,
                        "prn_event_id": normalize_text(sample_row.get("prn_event_id", "")),
                        "prn_mechanism_call": normalize_text(sample_row.get("prn_mechanism_call", "")),
                        "tool": "ismapper",
                        "is_element_name": is_element_name,
                        "query_reference_id": query_reference_id,
                        "reference_record": reference_record,
                        "locus_start": str(locus_start),
                        "locus_end": str(locus_end),
                        "prn_overlap_bp": str(prn_overlap_bp),
                        "call_label": normalize_text(row.get("call", "")),
                        "orientation": normalize_text(row.get("orientation", "")),
                        "gap_bp": normalize_text(row.get("gap", "")),
                        "percent_id": normalize_text(row.get("percent_ID", "")),
                        "percent_cov": normalize_text(row.get("percent_cov", "")),
                        "gene_interruption": normalize_text(row.get("gene_interruption", "")),
                        "left_gene": normalize_text(row.get("left_gene", "")),
                        "left_description": left_description,
                        "right_gene": normalize_text(row.get("right_gene", "")),
                        "right_description": right_description,
                        "start_clipped_reads": "",
                        "end_clipped_reads": "",
                        "total_clipped_reads": "",
                        "direct_repeats": "",
                        "inverted_repeats": "",
                        "left_sequence": "",
                        "right_sequence": "",
                        "evidence_note": "prn_overlap_detected_by_ismapper",
                    }
                )
    return evidence_rows, metrics


def parse_panisa_hits(
    sample_id: str,
    sample_row: dict[str, str],
    work_root: Path,
    prn_locus: dict[str, str | int],
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, int]]:
    evidence_rows: list[dict[str, str]] = []
    tsd_rows: list[dict[str, str]] = []
    panisa_path = work_root / "panisa" / f"{sample_id}.panisa.tsv"
    metrics = {"ismapper_tables_found": 0, "panisa_file_found": 1 if panisa_path.exists() else 0}
    if not panisa_path.exists():
        return evidence_rows, tsd_rows, metrics

    prn_start = int(prn_locus["start"])
    prn_end = int(prn_locus["end"])
    prn_record = normalize_text(str(prn_locus["record_id"]))

    with panisa_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chromosome = normalize_text(row.get("Chromosome", ""))
            if chromosome and prn_record and chromosome != prn_record:
                continue
            start_pos = parse_int(row.get("Start position", ""))
            end_pos = parse_int(row.get("End position", ""))
            if start_pos is None or end_pos is None:
                continue
            locus_start = min(start_pos, end_pos)
            locus_end = max(start_pos, end_pos)
            prn_overlap_bp = overlaps(locus_start, locus_end, prn_start, prn_end)
            if prn_overlap_bp <= 0:
                continue

            start_clipped = parse_int(row.get("Start clipped reads", "")) or 0
            end_clipped = parse_int(row.get("End clipped reads", "")) or 0
            total_clipped = start_clipped + end_clipped
            inferred_is_element_name = infer_panisa_element(row)
            evidence_rows.append(
                {
                    "sample_id_canonical": sample_id,
                    "prn_event_id": normalize_text(sample_row.get("prn_event_id", "")),
                    "prn_mechanism_call": normalize_text(sample_row.get("prn_mechanism_call", "")),
                    "tool": "panisa",
                    "is_element_name": inferred_is_element_name,
                    "query_reference_id": "",
                    "reference_record": chromosome,
                    "locus_start": str(locus_start),
                    "locus_end": str(locus_end),
                    "prn_overlap_bp": str(prn_overlap_bp),
                    "call_label": "panisa_insertion_site",
                    "orientation": "",
                    "gap_bp": str(locus_end - locus_start + 1),
                    "percent_id": "",
                    "percent_cov": "",
                    "gene_interruption": "",
                    "left_gene": "",
                    "left_description": "",
                    "right_gene": "",
                    "right_description": "",
                    "start_clipped_reads": str(start_clipped),
                    "end_clipped_reads": str(end_clipped),
                    "total_clipped_reads": str(total_clipped),
                    "direct_repeats": normalize_text(row.get("Direct repeats", "")),
                    "inverted_repeats": normalize_text(row.get("Inverted repeats", "")),
                    "left_sequence": normalize_text(row.get("Left sequence", "")),
                    "right_sequence": normalize_text(row.get("Right sequence", "")),
                    "evidence_note": "prn_overlap_detected_by_panisa",
                }
            )
            tsd_rows.append(
                {
                    "sample_id_canonical": sample_id,
                    "prn_event_id": normalize_text(sample_row.get("prn_event_id", "")),
                    "prn_mechanism_call": normalize_text(sample_row.get("prn_mechanism_call", "")),
                    "inferred_is_element_name": inferred_is_element_name,
                    "chromosome": chromosome,
                    "insertion_start": str(locus_start),
                    "insertion_end": str(locus_end),
                    "total_clipped_reads": str(total_clipped),
                    "direct_repeats": normalize_text(row.get("Direct repeats", "")),
                    "inverted_repeats": normalize_text(row.get("Inverted repeats", "")),
                    "left_sequence": normalize_text(row.get("Left sequence", "")),
                    "right_sequence": normalize_text(row.get("Right sequence", "")),
                }
            )
    return evidence_rows, tsd_rows, metrics


def choose_primary_hit(evidence_rows: list[dict[str, str]]) -> dict[str, str]:
    def sort_key(row: dict[str, str]) -> tuple[int, int, float, int]:
        is481_rank = 0 if normalize_text(row.get("is_element_name", "")).upper() == "IS481" else 1
        overlap_rank = -(parse_int(row.get("prn_overlap_bp", "")) or 0)
        tool_rank = 0 if normalize_text(row.get("tool", "")) == "panisa" else 1
        clipped = parse_int(row.get("total_clipped_reads", "")) or 0
        percent_cov = parse_float(row.get("percent_cov", "")) or 0.0
        locus_start = parse_int(row.get("locus_start", "")) or 0
        return (is481_rank, overlap_rank, tool_rank, -max(clipped, percent_cov), locus_start)

    return sorted(evidence_rows, key=sort_key)[0]


def build_status_fields(
    sample_id: str,
    evidence_rows: list[dict[str, str]],
    batch_samples: set[str],
    metrics: dict[str, int],
) -> dict[str, str]:
    has_ismapper = any(row.get("tool") == "ismapper" for row in evidence_rows)
    has_panisa = any(row.get("tool") == "panisa" for row in evidence_rows)
    total_supporting_reads = sum(parse_int(row.get("total_clipped_reads", "")) or 0 for row in evidence_rows)

    if evidence_rows:
        primary_hit = choose_primary_hit(evidence_rows)
        element = normalize_text(primary_hit.get("is_element_name", "")).casefold() or "untyped"
        if has_ismapper and has_panisa:
            return {
                "read_validation_status": "supported_concordant",
                "read_support_class": f"{element}_ismapper_panisa",
                "n_supporting_reads": str(total_supporting_reads) if total_supporting_reads else "",
                "n_contradicting_reads": "",
                "junction_supported": "yes",
            }
        if has_panisa:
            return {
                "read_validation_status": "supported",
                "read_support_class": f"{element}_panisa_only",
                "n_supporting_reads": str(total_supporting_reads) if total_supporting_reads else "",
                "n_contradicting_reads": "",
                "junction_supported": "yes",
            }
        return {
            "read_validation_status": "supported_candidate",
            "read_support_class": f"{element}_ismapper_only",
            "n_supporting_reads": "",
            "n_contradicting_reads": "",
            "junction_supported": "candidate_only",
        }

    if sample_id in batch_samples:
        if metrics.get("ismapper_tables_found", 0) or metrics.get("panisa_file_found", 0):
            return {
                "read_validation_status": "no_prn_is_signal_detected",
                "read_support_class": "no_prn_local_is_signal",
                "n_supporting_reads": "",
                "n_contradicting_reads": "",
                "junction_supported": "no",
            }
        return {
            "read_validation_status": "tool_output_missing",
            "read_support_class": "run_incomplete",
            "n_supporting_reads": "",
            "n_contradicting_reads": "",
            "junction_supported": "not_evaluated",
        }

    return {
        "read_validation_status": "unresolved",
        "read_support_class": "not_run",
        "n_supporting_reads": "",
        "n_contradicting_reads": "",
        "junction_supported": "not_evaluated",
    }


def build_notes(
    row: dict[str, str],
    *,
    batch_label: str,
    evidence_rows: list[dict[str, str]],
    metrics: dict[str, int],
    used_actual_validation: bool,
) -> str:
    parts = base_notes(row)
    if not used_actual_validation:
        parts.insert(0, "read_alignment_not_run_in_this_iteration")
        return ";".join(parts)

    parts.append(f"batch_label={batch_label or 'unspecified'}")
    parts.append(f"ismapper_tables_found={metrics.get('ismapper_tables_found', 0)}")
    parts.append(f"panisa_file_found={metrics.get('panisa_file_found', 0)}")
    parts.append(f"prn_hit_count={len(evidence_rows)}")
    if evidence_rows:
        primary_hit = choose_primary_hit(evidence_rows)
        parts.append(
            "primary_hit="
            + ",".join(
                [
                    normalize_text(primary_hit.get("tool", "")) or "unknown_tool",
                    normalize_text(primary_hit.get("is_element_name", "")) or "unknown_is",
                    normalize_text(primary_hit.get("locus_start", "")) or "",
                    normalize_text(primary_hit.get("locus_end", "")) or "",
                ]
            )
        )
    return ";".join(parts)


def build_output_rows(
    validation_rows: list[dict[str, str]],
    *,
    work_root: Path | None,
    batch_samples: set[str],
    batch_label: str,
    prn_locus: dict[str, str | int] | None,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    output_rows: list[dict[str, str]] = []
    evidence_rows: list[dict[str, str]] = []
    tsd_rows: list[dict[str, str]] = []
    actual_validation_enabled = work_root is not None and prn_locus is not None

    for row in validation_rows:
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        sample_evidence: list[dict[str, str]] = []
        sample_tsd_rows: list[dict[str, str]] = []
        sample_metrics = {"ismapper_tables_found": 0, "panisa_file_found": 0}

        if actual_validation_enabled and sample_id:
            ismapper_rows, ismapper_metrics = parse_ismapper_hits(sample_id, row, work_root, prn_locus)
            panisa_rows, panisa_tsd_rows, panisa_metrics = parse_panisa_hits(sample_id, row, work_root, prn_locus)
            sample_evidence.extend(ismapper_rows)
            sample_evidence.extend(panisa_rows)
            sample_tsd_rows.extend(panisa_tsd_rows)
            sample_metrics["ismapper_tables_found"] = ismapper_metrics["ismapper_tables_found"]
            sample_metrics["panisa_file_found"] = panisa_metrics["panisa_file_found"]

        status_fields = build_status_fields(sample_id, sample_evidence, batch_samples, sample_metrics)
        output_rows.append(
            {
                "sample_id_canonical": sample_id,
                "sra_run_accession": normalize_text(row.get("sra_run_accession", "")),
                "prn_event_id": normalize_text(row.get("prn_event_id", "")),
                "prn_mechanism_call": normalize_text(row.get("prn_mechanism_call", "")),
                "read_validation_status": status_fields["read_validation_status"],
                "read_support_class": status_fields["read_support_class"],
                "n_supporting_reads": status_fields["n_supporting_reads"],
                "n_contradicting_reads": status_fields["n_contradicting_reads"],
                "junction_supported": status_fields["junction_supported"],
                "targeted_locus_assembly_status": targeted_locus_assembly_status(row),
                "validation_method": READ_VALIDATION_METHOD if actual_validation_enabled else AUDIT_VALIDATION_METHOD,
                "validator_version": f"{VALIDATOR_VERSION};python={platform.python_version()}",
                "notes": build_notes(
                    row,
                    batch_label=batch_label,
                    evidence_rows=sample_evidence,
                    metrics=sample_metrics,
                    used_actual_validation=actual_validation_enabled,
                ),
            }
        )
        evidence_rows.extend(sample_evidence)
        tsd_rows.extend(sample_tsd_rows)

    return output_rows, evidence_rows, tsd_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a machine-readable PRN read-validation table for the selected validation subset."
    )
    parser.add_argument(
        "--subset",
        type=Path,
        default=PROJECT_ROOT / "step4_prn_validation" / "outputs" / "bp_prn_validation_subset.tsv",
        help="Validation subset manifest from VAL-01.",
    )
    parser.add_argument(
        "--is-work-root",
        type=Path,
        default=None,
        help="Stage-local read-validation work root containing ismapper/ and panisa/.",
    )
    parser.add_argument(
        "--batch",
        type=Path,
        default=None,
        help="Batch manifest produced by step4_03d_build_read_validation_batch.py.",
    )
    parser.add_argument(
        "--reference-gbff",
        type=Path,
        default=PROJECT_ROOT
        / "step2_typing"
        / "_ref/GCF_000195715.1/ncbi_dataset/data/GCF_000195715.1/genomic.gbff",
        help="Reference GBFF used to locate the prn locus in Tohama I.",
    )
    parser.add_argument(
        "--batch-label",
        default="",
        help="Human-readable batch label appended into notes when actual read validation is parsed.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=PROJECT_ROOT / "step4_prn_validation" / "outputs" / "bp_prn_read_validation.tsv",
        help="Read-validation output table.",
    )
    parser.add_argument(
        "--evidence-out",
        type=Path,
        default=PROJECT_ROOT / "step4_prn_validation" / "outputs" / "bp_prn_read_validation_is_calls.tsv",
        help="Per-hit evidence table joining ISMapper and panISa prn-local signals.",
    )
    parser.add_argument(
        "--tsd-out",
        type=Path,
        default=PROJECT_ROOT / "step4_prn_validation" / "outputs" / "bp_prn_read_validation_tsd.tsv",
        help="panISa-derived direct-repeat and inverted-repeat calls overlapping prn.",
    )
    parser.add_argument(
        "--merge-base",
        type=Path,
        default=None,
        help="Existing read-validation table to merge untouched samples from.",
    )
    parser.add_argument(
        "--evidence-merge-base",
        type=Path,
        default=None,
        help="Existing per-hit evidence table to merge untouched samples from.",
    )
    parser.add_argument(
        "--tsd-merge-base",
        type=Path,
        default=None,
        help="Existing TSD evidence table to merge untouched samples from.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    validation_rows = load_tsv_rows(args.subset)
    sample_ids = [normalize_text(row.get("sample_id_canonical", "")) for row in validation_rows]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("duplicate sample_id_canonical values found in validation subset")

    batch_samples = load_batch_samples(args.batch, args.is_work_root)
    prn_locus = load_prn_locus(args.reference_gbff) if args.is_work_root is not None else None
    output_rows, evidence_rows, tsd_rows = build_output_rows(
        validation_rows,
        work_root=args.is_work_root,
        batch_samples=batch_samples,
        batch_label=args.batch_label,
        prn_locus=prn_locus,
    )

    if args.merge_base is not None and args.merge_base.exists():
        output_rows = merge_unique_sample_rows(load_tsv_rows(args.merge_base), output_rows)
    if args.evidence_merge_base is not None and args.evidence_merge_base.exists():
        evidence_rows = merge_multirow_sample_rows(load_tsv_rows(args.evidence_merge_base), evidence_rows)
    if args.tsd_merge_base is not None and args.tsd_merge_base.exists():
        tsd_rows = merge_multirow_sample_rows(load_tsv_rows(args.tsd_merge_base), tsd_rows)

    write_tsv(args.out, OUTPUT_COLUMNS, output_rows)
    write_tsv(args.evidence_out, EVIDENCE_COLUMNS, evidence_rows)
    write_tsv(args.tsd_out, TSD_COLUMNS, tsd_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
