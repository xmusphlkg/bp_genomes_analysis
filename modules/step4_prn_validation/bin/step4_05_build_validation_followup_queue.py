#!/usr/bin/env python3
"""Build a targeted validation follow-up queue.

The queue prioritizes origin-defining disrupted tips and unresolved rows from
the current 57-row raw-read validation subset. It classifies rows by practical
next action: recoverable reads, public long-read/hybrid exemplar already
available, assembly-only, or likely unrecoverable with current public inputs.
"""

from __future__ import annotations

import argparse
import glob
import math
import re
from pathlib import Path
from typing import Any

import pandas as pd

from step4_00_build_blocked_recovery_plan import (
    build_plan_index,
    choose_best_plan_row,
    recovery_status_for_plan_row,
)
from step4_02_scan_prn_mechanisms import STEP4_DATA_ROOT, WORKFLOW_DATA_ROOT


LONGREAD_RE = re.compile(r"nanopore|pacbio|hifi|rsii|sequel|minion|promethion", re.I)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def norm(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() == "nan" else text


def join_values(values: list[Any]) -> str:
    out: list[str] = []
    for value in values:
        text = norm(value)
        if text and text not in out:
            out.append(text)
    return ";".join(out)


def has_reads(row: pd.Series) -> bool:
    if norm(row.get("raw_reads_available")).lower() == "true":
        return True
    return bool(
        join_values(
            [
                row.get("sra_run_accession"),
                row.get("ena_run_accession"),
                row.get("read_accession_primary"),
                row.get("sra_run_accession_validation"),
                row.get("ena_run_accession_validation"),
                row.get("read_accession_primary_validation"),
            ]
        )
    )


def split_accessions(value: Any) -> list[str]:
    return [token for token in (norm(part) for part in norm(value).split(";")) if token]


def has_longread(row: pd.Series) -> bool:
    return bool(LONGREAD_RE.search(norm(row.get("sequencing_tech"))))


def has_longread_signal(row: pd.Series) -> bool:
    if has_longread(row):
        return True
    platform = norm(row.get("recovery_instrument_platform"))
    tech = norm(row.get("recovery_selected_run_accession"))
    return bool(LONGREAD_RE.search(platform)) or bool(LONGREAD_RE.search(tech))


def validation_level_for_queue_row(row: pd.Series) -> str:
    status = norm(row.get("read_validation_status"))
    if status in {"supported", "supported_concordant"}:
        return "read_backed_supported"
    if has_longread_signal(row):
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
    return ordering.get(norm(level), 99)


def load_download_plan_index(root: Path) -> dict[str, list[dict[str, str]]]:
    path = STEP4_DATA_ROOT / "inputs" / "bp_raw_reads_download_plan.tsv"
    if not path.exists():
        return {}
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return build_plan_index(frame.to_dict("records"))


def load_origin_tip_rows(root: Path) -> pd.DataFrame:
    origin = pd.read_csv(WORKFLOW_DATA_ROOT / "asr" / "origin_events.tsv", sep="\t", dtype=str)
    frames: list[pd.DataFrame] = []
    for path in sorted(glob.glob(str(WORKFLOW_DATA_ROOT / "asr" / "event_subtrees" / "origin_*.descendant_tips.tsv"))):
        origin_id = Path(path).name.split(".", 1)[0]
        frame = pd.read_csv(path, sep="\t", dtype=str)
        frame["origin_id"] = origin_id
        frames.append(frame[frame["observed_prn_state"].eq("disrupted")].copy())
    if not frames:
        return pd.DataFrame()
    tips = pd.concat(frames, ignore_index=True)
    origin_keep = origin[
        [
            "origin_id",
            "n_tips_disrupted",
            "n_tips_total",
            "dominant_prn_mechanism",
            "branch_support",
            "origin_support_score",
        ]
    ].copy()
    return tips.merge(origin_keep, on="origin_id", how="left")


def load_validation_rows(root: Path) -> pd.DataFrame:
    path = root / "manuscript" / "figure_data" / "figure6_targeted_validation_followup.tsv"
    if not path.exists():
        return pd.DataFrame()
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    for column in [
        "sample_id_canonical",
        "read_validation_status",
        "read_support_class",
        "n_supporting_reads",
        "n_contradicting_reads",
        "junction_supported",
        "targeted_locus_assembly_status",
        "read_accession_source",
        "evidence_flags",
        "notes",
    ]:
        if column not in frame.columns:
            frame[column] = ""
    frame["read_accession_source"] = frame["read_accession_source"].replace("", "targeted_followup_queue")
    return frame.drop_duplicates("sample_id_canonical", keep="first").reset_index(drop=True)


def load_sample_metadata(root: Path) -> pd.DataFrame:
    qc_manifest = pd.read_csv(WORKFLOW_DATA_ROOT / "manifest" / "genome_catalog.tsv", sep="\t", dtype=str)
    mechanism = pd.read_csv(root / "manuscript" / "figure_data" / "fig02_prn_mechanism_calls.tsv", sep="\t", dtype=str)
    qc_columns = [
        "sample_id_canonical",
        "assembly_accession",
        "assembly_level",
        "sra_run_accession",
        "ena_run_accession",
        "country_iso3",
        "country",
        "year",
        "sequencing_tech",
        "raw_reads_available",
        "raw_read_link_status",
    ]
    mechanism_columns = [
        "sample_id_canonical",
        "prn_mechanism_call",
        "prn_event_id",
        "prn_call_confidence",
        "country_iso3",
        "year",
        "evidence_flags",
        "notes",
    ]
    qc_subset = qc_manifest[[column for column in qc_columns if column in qc_manifest.columns]]
    mechanism_subset = mechanism[[column for column in mechanism_columns if column in mechanism.columns]]
    merged = mechanism_subset.merge(qc_subset, on="sample_id_canonical", how="left", suffixes=("", "_qc"))
    return merged.drop_duplicates("sample_id_canonical")


def attach_recovery_plan(queue: pd.DataFrame, root: Path) -> pd.DataFrame:
    plan_index = load_download_plan_index(root)
    recovery_rows: list[dict[str, str]] = []
    for _, row in queue.iterrows():
        preferred_runs: list[str] = []
        for field in ("sra_run_accession", "ena_run_accession", "read_accession_primary"):
            preferred_runs.extend(split_accessions(row.get(field, "")))
        sample_id = norm(row.get("sample_id_canonical"))
        selected = choose_best_plan_row(sample_id, preferred_runs, plan_index)
        recovery_rows.append(
            {
                "sample_id_canonical": sample_id,
                "recovery_selected_run_accession": norm(selected.get("run_accession", "")) if selected else "",
                "recovery_run_compatibility": norm(selected.get("run_compatibility", "")) if selected else "",
                "recovery_library_layout": norm(selected.get("ena_library_layout", "")) if selected else "",
                "recovery_instrument_platform": norm(selected.get("ena_instrument_platform", "")) if selected else "",
                "recovery_download_strategy": norm(selected.get("download_strategy", "")) if selected else "",
                "recovery_ena_fastq_ftp": norm(selected.get("ena_fastq_ftp", "")) if selected else "",
                "recovery_estimated_total_bytes": norm(selected.get("estimated_total_bytes", "")) if selected else "",
                "recovery_plan_status": recovery_status_for_plan_row(selected),
            }
        )
    recovery = pd.DataFrame(recovery_rows)
    return queue.merge(recovery, on="sample_id_canonical", how="left")


def choose_priority_origin_exemplars(queue: pd.DataFrame) -> set[str]:
    origin_rows = queue[queue["is_origin_defining_tip"].eq("True")].copy()
    if origin_rows.empty:
        return set()

    def rank(row: pd.Series) -> tuple:
        # Prespecified package-exemplar rule: prefer the strongest within-package
        # orthogonal anchor before falling back to candidate, recoverable, or
        # assembly-only rows.
        validation_level = validation_level_for_queue_row(row)
        if validation_level in {
            "read_backed_supported",
            "public_longread_or_hybrid_assembly",
            "read_backed_candidate",
        }:
            evidence_rank = validation_level_rank(validation_level)
        elif norm(row.get("recovery_plan_status")) == "recoverable_paired_illumina":
            evidence_rank = 3
        elif norm(row.get("prn_call_confidence")) in {"assembly_high", "assembly_moderate"}:
            evidence_rank = 4
        else:
            evidence_rank = 5
        try:
            disrupted = -int(float(norm(row.get("n_tips_disrupted")) or 0))
        except ValueError:
            disrupted = 0
        try:
            year = int(float(norm(row.get("year"))))
        except ValueError:
            year = 9999
        return (norm(row.get("origin_id")), evidence_rank, disrupted, year, norm(row.get("sample_id_canonical")))

    selected: set[str] = set()
    for origin_id, frame in origin_rows.groupby("origin_id", dropna=False):
        ranked = sorted(frame.to_dict("records"), key=lambda row: rank(pd.Series(row)))
        if ranked:
            selected.add(norm(ranked[0].get("sample_id_canonical")))
    return selected


def classify_followup(row: pd.Series) -> tuple[str, str, str, str]:
    status = norm(row.get("read_validation_status"))
    support = norm(row.get("read_support_class"))
    mechanism = norm(row.get("prn_mechanism_call"))
    confidence = norm(row.get("prn_call_confidence"))
    target_status = norm(row.get("targeted_locus_assembly_status"))
    recovery_status = norm(row.get("recovery_plan_status"))
    compatibility = norm(row.get("recovery_run_compatibility"))
    platform = norm(row.get("recovery_instrument_platform"))

    if status in {"supported", "supported_concordant", "supported_candidate"}:
        return (
            "read_backed_or_candidate_available",
            "retain_as_read_supported_or_candidate; optional junction polish only if this sample anchors a major origin package",
            "resolved_or_candidate_read_signal_present",
            "1",
        )
    if has_longread(row):
        return (
            "public_longread_or_hybrid_exemplar_present",
            "use public long-read or hybrid assembly as the event anchor; extract junction evidence if a hard anchor is needed",
            "sequencing_tech_contains_long_read_keyword",
            "1",
        )
    if recovery_status == "recoverable_paired_illumina":
        if status == "unresolved" or status == "not_run" or not status:
            return (
                "can_recover_reads",
                "rerun targeted ISMapper/panISa and junction extraction on listed read accessions",
                "download plan provides paired-Illumina FASTQ inputs for the current short-read validator",
                "2",
            )
        return (
            "can_recheck_reads",
            "inspect read-level negative or ambiguous call if this row becomes an origin-package anchor",
            "download plan provides paired-Illumina FASTQ inputs for the current short-read validator",
            "3",
        )
    if recovery_status == "linked_incompatible_run_current_short_read_validator":
        return (
            "linked_incompatible_run_current_short_read_validator",
            "retain as a separate long-read/nonpaired follow-up track; current ISMapper/panISa workflow assumes paired Illumina and should not be forced onto this run",
            f"download plan selects {compatibility or 'incompatible'} input on platform {platform or 'unknown'}",
            "2" if has_longread_signal(row) else "3",
        )
    if recovery_status == "linked_run_without_fastq_ftp":
        return (
            "linked_run_without_fastq_ftp",
            "retain as public run metadata only; recovery requires submitted files or a different accession because indexed FASTQ FTP endpoints are absent",
            "download plan resolves a linked run but no paired FASTQ FTP endpoints are indexed",
            "3",
        )
    if mechanism.startswith("coding_disrupted_"):
        return (
            "assembly_only",
            "retain as assembly-defined; seek alternate raw-read links or a same-event long-read exemplar",
            f"{confidence or 'assembly'} disrupted assembly call but no public read accession in current manifests",
            "3",
        )
    if mechanism in {"insufficient_data", "uncertain_fragmented_assembly"} or target_status in {
        "assembly_sequence_missing",
        "no_current_step3_prn_input",
        "partial_prn_alignment",
    }:
        return (
            "likely_unrecoverable_current_inputs",
            "do not interpret as negative evidence; recovery requires new input assembly, alternate accession, or original reads",
            "locus/input unavailable or fragmented and no current recoverable read accession",
            "4",
        )
    return (
        "manual_review",
        "manual review before assigning validation effort",
        "row does not meet automatic follow-up rules",
        "4",
    )


def build_queue(root: Path) -> pd.DataFrame:
    origin_tips = load_origin_tip_rows(root)
    validation_rows = load_validation_rows(root)
    sample_meta = load_sample_metadata(root)

    origin_samples = set(origin_tips["sample_id_canonical"].dropna()) if not origin_tips.empty else set()
    unresolved_samples = set(validation_rows.loc[validation_rows["read_validation_status"].eq("unresolved"), "sample_id_canonical"].dropna())
    queue_samples = sorted(origin_samples | unresolved_samples)

    base = pd.DataFrame({"sample_id_canonical": queue_samples})
    if not origin_tips.empty:
        origin_agg = (
            origin_tips.sort_values(["origin_id", "sample_id_canonical"])
            .groupby("sample_id_canonical", dropna=False)
            .agg(
                origin_id=("origin_id", lambda x: ";".join(sorted(set(norm(v) for v in x if norm(v))))),
                n_tips_disrupted=("n_tips_disrupted", lambda x: join_values(list(x))),
                n_tips_total=("n_tips_total", lambda x: join_values(list(x))),
                dominant_prn_mechanism=("dominant_prn_mechanism", lambda x: join_values(list(x))),
                branch_support=("branch_support", lambda x: join_values(list(x))),
                origin_support_score=("origin_support_score", lambda x: join_values(list(x))),
                origin_tip_label=("tip_label", lambda x: join_values(list(x))),
                origin_observed_prn_mechanism_call=("observed_prn_mechanism_call", lambda x: join_values(list(x))),
            )
            .reset_index()
        )
        base = base.merge(origin_agg, on="sample_id_canonical", how="left")
    else:
        base["origin_id"] = ""

    validation_keep = validation_rows.add_suffix("_validation").rename(columns={"sample_id_canonical_validation": "sample_id_canonical"})
    base = base.merge(validation_keep, on="sample_id_canonical", how="left")
    base = base.merge(sample_meta.add_suffix("_manifest").rename(columns={"sample_id_canonical_manifest": "sample_id_canonical"}), on="sample_id_canonical", how="left")

    for column in [
        "read_validation_status",
        "read_support_class",
        "n_supporting_reads",
        "n_contradicting_reads",
        "junction_supported",
        "targeted_locus_assembly_status",
        "read_accession_source",
        "evidence_flags",
        "notes",
    ]:
        validation_column = f"{column}_validation"
        if validation_column in base.columns:
            base[column] = base[validation_column].where(base[validation_column].fillna("").astype(str).str.strip().ne(""), base.get(column, ""))

    def pick(row: pd.Series, field: str) -> str:
        return norm(row.get(f"{field}_validation")) or norm(row.get(f"{field}_manifest")) or norm(row.get(field))

    rows: list[dict[str, Any]] = []
    for _, row in base.iterrows():
        out: dict[str, Any] = {
            "sample_id_canonical": norm(row.get("sample_id_canonical")),
            "is_origin_defining_tip": "True" if norm(row.get("origin_id")) else "False",
            "is_unresolved_validation_row": "True"
            if norm(row.get("read_validation_status_validation")) == "unresolved"
            else "False",
            "origin_id": norm(row.get("origin_id")),
            "origin_tip_label": norm(row.get("origin_tip_label")),
            "origin_n_disrupted_descendants": norm(row.get("n_tips_disrupted")),
            "origin_n_total_descendants": norm(row.get("n_tips_total")),
            "origin_dominant_prn_mechanism": norm(row.get("dominant_prn_mechanism")),
            "origin_branch_support": norm(row.get("branch_support")),
            "origin_support_score": norm(row.get("origin_support_score")),
            "assembly_accession": pick(row, "assembly_accession"),
            "assembly_level": norm(row.get("assembly_level_manifest")),
            "country_iso3": pick(row, "country_iso3"),
            "year": pick(row, "year"),
            "prn_mechanism_call": pick(row, "prn_mechanism_call")
            or norm(row.get("origin_observed_prn_mechanism_call")),
            "prn_event_id": pick(row, "prn_event_id"),
            "prn_call_confidence": pick(row, "prn_call_confidence"),
            "sra_run_accession": pick(row, "sra_run_accession"),
            "ena_run_accession": pick(row, "ena_run_accession"),
            "read_accession_primary": norm(row.get("read_accession_primary_validation")),
            "raw_reads_available": pick(row, "raw_reads_available"),
            "raw_read_link_status": pick(row, "raw_read_link_status"),
            "sequencing_tech": norm(row.get("sequencing_tech_manifest")),
            "read_validation_status": norm(row.get("read_validation_status_validation")),
            "read_support_class": norm(row.get("read_support_class_validation")),
            "n_supporting_reads": norm(row.get("n_supporting_reads_validation")),
            "n_contradicting_reads": norm(row.get("n_contradicting_reads_validation")),
            "junction_supported": norm(row.get("junction_supported_validation")),
            "targeted_locus_assembly_status": norm(row.get("targeted_locus_assembly_status_validation")),
            "selection_stratum": norm(row.get("selection_stratum_validation")),
            "selection_reason": norm(row.get("selection_reason_validation")),
            "evidence_flags": pick(row, "evidence_flags"),
            "notes": pick(row, "notes"),
        }
        rows.append(out)

    queue = pd.DataFrame(rows)
    queue = attach_recovery_plan(queue, root)
    priority_exemplars = choose_priority_origin_exemplars(queue)
    followup = queue.apply(lambda row: classify_followup(row), axis=1, result_type="expand")
    followup.columns = ["followup_class", "recommended_action", "followup_rationale", "priority_tier"]
    queue = pd.concat([queue, followup], axis=1)
    queue["is_priority_origin_exemplar"] = queue["sample_id_canonical"].isin(priority_exemplars).map(lambda value: "True" if value else "False")
    queue["queue_scope"] = queue.apply(
        lambda row: ";".join(
            scope
            for scope, present in [
                ("origin_defining_disrupted_tip", row["is_origin_defining_tip"] == "True"),
                ("unresolved_validation_row", row["is_unresolved_validation_row"] == "True"),
            ]
            if present
        ),
        axis=1,
    )
    queue["queue_id"] = [f"val_followup_{idx:04d}" for idx in range(1, len(queue) + 1)]
    sort_cols = ["priority_tier", "is_priority_origin_exemplar", "origin_id", "followup_class", "sample_id_canonical"]
    queue = queue.sort_values(sort_cols, ascending=[True, False, True, True, True]).reset_index(drop=True)
    queue["queue_id"] = [f"val_followup_{idx:04d}" for idx in range(1, len(queue) + 1)]
    return queue


def build_arg_parser() -> argparse.ArgumentParser:
    root = repo_root()
    parser = argparse.ArgumentParser(description="Build targeted prn validation follow-up queue.")
    parser.add_argument(
        "--out",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_targeted_validation_followup_queue.tsv",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    queue = build_queue(repo_root())
    columns = [
        "queue_id",
        "queue_scope",
        "priority_tier",
        "is_priority_origin_exemplar",
        "is_origin_defining_tip",
        "is_unresolved_validation_row",
        "origin_id",
        "origin_tip_label",
        "origin_n_disrupted_descendants",
        "origin_n_total_descendants",
        "origin_dominant_prn_mechanism",
        "origin_branch_support",
        "origin_support_score",
        "sample_id_canonical",
        "assembly_accession",
        "assembly_level",
        "country_iso3",
        "year",
        "prn_mechanism_call",
        "prn_event_id",
        "prn_call_confidence",
        "sra_run_accession",
        "ena_run_accession",
        "read_accession_primary",
        "raw_reads_available",
        "raw_read_link_status",
        "sequencing_tech",
        "recovery_plan_status",
        "recovery_selected_run_accession",
        "recovery_run_compatibility",
        "recovery_library_layout",
        "recovery_instrument_platform",
        "recovery_download_strategy",
        "recovery_estimated_total_bytes",
        "read_validation_status",
        "read_support_class",
        "n_supporting_reads",
        "n_contradicting_reads",
        "junction_supported",
        "targeted_locus_assembly_status",
        "recovery_ena_fastq_ftp",
        "selection_stratum",
        "selection_reason",
        "followup_class",
        "recommended_action",
        "followup_rationale",
        "evidence_flags",
        "notes",
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    queue.to_csv(args.out, sep="\t", index=False, columns=columns, lineterminator="\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
