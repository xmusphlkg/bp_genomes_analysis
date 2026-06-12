#!/usr/bin/env python3
"""Build a recovery plan for blocked Step4 read-validation rows.

The blocked manifest only tells us that local FASTQ/BAM inputs are missing.
This planner bridges those blocked rows back to the broader raw-read download
catalog, picks the best public run per sample, and explicitly distinguishes:

1. paired-Illumina runs that can feed the current ISMapper/panISa validator;
2. linked but incompatible runs (for example long-read or non-paired inputs);
3. linked runs without indexed FASTQ FTP endpoints; and
4. samples with no matching download-plan row.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from step4_02_scan_prn_mechanisms import (
    load_tsv_rows,
    normalize_text,
    project_module_data_root,
    project_workflow_root,
    write_tsv,
)


STEP1_DATA_ROOT = project_module_data_root("step1_ingest")
STEP4_DATA_ROOT = project_module_data_root("step4_prn_validation")
WORKFLOW_DATA_ROOT = project_workflow_root()


OUTPUT_COLUMNS = [
    "sample_id_canonical",
    "assembly_accession",
    "prn_event_id",
    "prn_mechanism_call",
    "prn_call_confidence",
    "raw_read_link_status",
    "missing_inputs",
    "preferred_run_accessions",
    "candidate_run_count",
    "candidate_run_accessions",
    "selected_run_accession",
    "selected_run_source",
    "selected_run_compatibility",
    "selected_library_layout",
    "selected_instrument_platform",
    "selected_download_strategy",
    "ena_fastq_ftp",
    "ena_fastq_md5",
    "estimated_total_bytes",
    "fallback_run_accession",
    "fallback_run_source",
    "fallback_run_compatibility",
    "fallback_library_layout",
    "fallback_instrument_platform",
    "fallback_download_strategy",
    "fallback_ena_fastq_ftp",
    "fallback_ena_fastq_md5",
    "fallback_estimated_total_bytes",
    "recovery_plan_status",
    "recovery_priority_tier",
    "recovery_plan_note",
]


def split_tokens(value: str) -> list[str]:
    return [token for token in (normalize_text(part) for part in normalize_text(value).split(";")) if token]


def parse_int(value: str) -> int:
    text = normalize_text(value)
    if not text:
        return 10**18
    try:
        return int(float(text))
    except ValueError:
        return 10**18


def ftp_pair_count(row: dict[str, str]) -> int:
    return len(split_tokens(row.get("ena_fastq_ftp", "")))


def recovery_status_for_plan_row(row: dict[str, str] | None) -> str:
    if row is None:
        return "no_download_plan_match"

    compatibility = normalize_text(row.get("run_compatibility", ""))
    if compatibility == "paired_illumina_fastq" and ftp_pair_count(row) == 2:
        return "recoverable_paired_illumina"
    if compatibility == "no_fastq_ftp" or not split_tokens(row.get("ena_fastq_ftp", "")):
        return "linked_run_without_fastq_ftp"
    return "linked_incompatible_run_current_short_read_validator"


def recovery_priority_tier(status: str) -> str:
    if status == "recoverable_paired_illumina":
        return "1"
    if status == "linked_incompatible_run_current_short_read_validator":
        return "2"
    if status == "linked_run_without_fastq_ftp":
        return "3"
    return "4"


def status_rank(status: str) -> int:
    return {
        "recoverable_paired_illumina": 0,
        "linked_incompatible_run_current_short_read_validator": 1,
        "linked_run_without_fastq_ftp": 2,
        "no_download_plan_match": 3,
    }.get(status, 9)


def build_plan_index(plan_rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in plan_rows:
        sample_id = normalize_text(row.get("sample_id_canonical", ""))
        if sample_id:
            by_sample[sample_id].append(row)
    return by_sample


def plan_sort_key(row: dict[str, str], preferred_set: set[str]) -> tuple[int, int, int, str]:
    run = normalize_text(row.get("run_accession", ""))
    status = recovery_status_for_plan_row(row)
    return (
        status_rank(status),
        0 if run in preferred_set else 1,
        parse_int(row.get("estimated_total_bytes", "")),
        run,
    )


def ranked_plan_rows(
    sample_id: str,
    preferred_runs: list[str],
    plan_index: dict[str, list[dict[str, str]]],
) -> list[dict[str, str]]:
    candidates = plan_index.get(sample_id, [])
    if not candidates:
        return []

    preferred_set = set(preferred_runs)
    return sorted(candidates, key=lambda row: plan_sort_key(row, preferred_set))


def choose_best_plan_row(sample_id: str, preferred_runs: list[str], plan_index: dict[str, list[dict[str, str]]]) -> dict[str, str] | None:
    ranked = ranked_plan_rows(sample_id, preferred_runs, plan_index)
    return ranked[0] if ranked else None


def choose_fallback_plan_row(
    selected: dict[str, str] | None,
    ranked_candidates: list[dict[str, str]],
) -> dict[str, str] | None:
    if selected is None:
        return None

    selected_run = normalize_text(selected.get("run_accession", ""))
    for candidate in ranked_candidates:
        run = normalize_text(candidate.get("run_accession", ""))
        if not run or run == selected_run:
            continue
        if recovery_status_for_plan_row(candidate) != "recoverable_paired_illumina":
            continue
        return candidate
    return None


def build_recovery_row(row: dict[str, str], plan_index: dict[str, list[dict[str, str]]]) -> dict[str, str]:
    sample_id = normalize_text(row.get("sample_id_canonical", ""))
    preferred_runs = []
    for field in ("sra_run_accession", "ena_run_accession", "read_accession_primary"):
        preferred_runs.extend(split_tokens(row.get(field, "")))
    preferred_runs = list(dict.fromkeys(preferred_runs))

    candidates = plan_index.get(sample_id, [])
    ranked_candidates = ranked_plan_rows(sample_id, preferred_runs, plan_index)
    selected = ranked_candidates[0] if ranked_candidates else None
    fallback = choose_fallback_plan_row(selected, ranked_candidates)
    status = recovery_status_for_plan_row(selected)

    note_parts = [f"candidate_runs={len(candidates)}"]
    if preferred_runs:
        note_parts.append(f"preferred_runs={';'.join(preferred_runs)}")
    if selected is not None:
        note_parts.extend(
            [
                f"selected_run={normalize_text(selected.get('run_accession', ''))}",
                f"selected_compatibility={normalize_text(selected.get('run_compatibility', '')) or 'unknown'}",
            ]
        )
    else:
        note_parts.append("selected_run=none")
    if fallback is not None:
        note_parts.extend(
            [
                f"fallback_run={normalize_text(fallback.get('run_accession', ''))}",
                f"fallback_compatibility={normalize_text(fallback.get('run_compatibility', '')) or 'unknown'}",
            ]
        )

    return {
        "sample_id_canonical": sample_id,
        "assembly_accession": normalize_text(row.get("assembly_accession", "")),
        "prn_event_id": normalize_text(row.get("prn_event_id", "")),
        "prn_mechanism_call": normalize_text(row.get("prn_mechanism_call", "")),
        "prn_call_confidence": normalize_text(row.get("prn_call_confidence", "")),
        "raw_read_link_status": normalize_text(row.get("raw_read_link_status", "")),
        "missing_inputs": normalize_text(row.get("missing_inputs", "")),
        "preferred_run_accessions": ";".join(preferred_runs),
        "candidate_run_count": str(len(candidates)),
        "candidate_run_accessions": ";".join(
            normalize_text(candidate.get("run_accession", "")) for candidate in candidates if normalize_text(candidate.get("run_accession", ""))
        ),
        "selected_run_accession": normalize_text(selected.get("run_accession", "")) if selected else "",
        "selected_run_source": normalize_text(selected.get("run_source", "")) if selected else "",
        "selected_run_compatibility": normalize_text(selected.get("run_compatibility", "")) if selected else "",
        "selected_library_layout": normalize_text(selected.get("ena_library_layout", "")) if selected else "",
        "selected_instrument_platform": normalize_text(selected.get("ena_instrument_platform", "")) if selected else "",
        "selected_download_strategy": normalize_text(selected.get("download_strategy", "")) if selected else "",
        "ena_fastq_ftp": normalize_text(selected.get("ena_fastq_ftp", "")) if selected else "",
        "ena_fastq_md5": normalize_text(selected.get("ena_fastq_md5", "")) if selected else "",
        "estimated_total_bytes": normalize_text(selected.get("estimated_total_bytes", "")) if selected else "",
        "fallback_run_accession": normalize_text(fallback.get("run_accession", "")) if fallback else "",
        "fallback_run_source": normalize_text(fallback.get("run_source", "")) if fallback else "",
        "fallback_run_compatibility": normalize_text(fallback.get("run_compatibility", "")) if fallback else "",
        "fallback_library_layout": normalize_text(fallback.get("ena_library_layout", "")) if fallback else "",
        "fallback_instrument_platform": normalize_text(fallback.get("ena_instrument_platform", "")) if fallback else "",
        "fallback_download_strategy": normalize_text(fallback.get("download_strategy", "")) if fallback else "",
        "fallback_ena_fastq_ftp": normalize_text(fallback.get("ena_fastq_ftp", "")) if fallback else "",
        "fallback_ena_fastq_md5": normalize_text(fallback.get("ena_fastq_md5", "")) if fallback else "",
        "fallback_estimated_total_bytes": normalize_text(fallback.get("estimated_total_bytes", "")) if fallback else "",
        "recovery_plan_status": status,
        "recovery_priority_tier": recovery_priority_tier(status),
        "recovery_plan_note": ";".join(note_parts),
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a recovery plan for blocked Step4 read-validation rows."
    )
    parser.add_argument(
        "--blocked",
        type=Path,
        default=STEP4_DATA_ROOT / "work" / "read_validation" / "current" / "bp_prn_read_validation_missing_inputs.tsv",
        help="Blocked-manifest TSV from step4_03d_build_read_validation_batch.py.",
    )
    parser.add_argument(
        "--download-plan",
        type=Path,
        default=STEP4_DATA_ROOT / "inputs" / "bp_raw_reads_download_plan.tsv",
        help="Global raw-read download plan used to select recoverable runs.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=STEP4_DATA_ROOT / "work" / "read_validation" / "current" / "bp_prn_blocked_recovery_plan.tsv",
        help="Recovery-plan output TSV.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    blocked_rows = [
        row
        for row in load_tsv_rows(args.blocked)
        if normalize_text(row.get("batch_status", "")) == "blocked_missing_local_inputs"
    ]
    plan_rows = load_tsv_rows(args.download_plan)
    plan_index = build_plan_index(plan_rows)
    output_rows = [build_recovery_row(row, plan_index) for row in blocked_rows]
    write_tsv(args.out, OUTPUT_COLUMNS, output_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
