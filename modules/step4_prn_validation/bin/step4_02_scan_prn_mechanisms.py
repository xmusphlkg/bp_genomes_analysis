#!/usr/bin/env python3
"""Classify prn disruption mechanisms from legacy step3 evidence and IS references."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.lib.project_paths import project_module_data_root, project_workflow_root


MISSING_TOKENS = {"", "na", "n/a", "none", "missing", "unknown", "not applicable"}
MECHANISM_COLUMNS = [
    "sample_id_canonical",
    "biosample_accession",
    "assembly_accession",
    "sra_run_accession",
    "country_iso3",
    "year",
    "mlst_st",
    "phylo_lineage",
    "prn_call_initial",
    "prn_mechanism_call",
    "prn_call_confidence",
    "prn_event_id",
    "prn_query_cov_pct",
    "prn_best_single_cov_pct",
    "prn_hsp_n",
    "bp_category",
    "insertion_subject_gap_bp",
    "is_element_best_hit",
    "is_element_best_hit_pident",
    "is_element_best_hit_qcov",
    "read_validation_status",
    "read_validation_support",
    "evidence_flags",
    "notes",
]

EVENT_COLUMNS = [
    "prn_event_id",
    "prn_mechanism_call",
    "prn_call_confidence_mode",
    "bp_category",
    "is_element_best_hit",
    "insertion_subject_gap_bp",
    "evidence_flags_signature",
    "sample_count",
    "assembly_count",
    "country_iso3_values",
    "year_min",
    "year_max",
    "example_sample_id_canonical",
    "example_assembly_accession",
    "notes",
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


STEP1_DATA_ROOT = project_module_data_root("step1_ingest")
STEP3_DATA_ROOT = project_module_data_root("step3_prn_scan")
STEP4_DATA_ROOT = project_module_data_root("step4_prn_validation")
PUBLIC_HEALTH_DATA_ROOT = project_module_data_root("public_health")
WORKFLOW_DATA_ROOT = project_workflow_root()


def normalize_text(value: str) -> str:
    return (value or "").strip()


def is_missing(value: str) -> bool:
    return normalize_text(value).casefold() in MISSING_TOKENS


def parse_int(value: str) -> int | None:
    value = normalize_text(value)
    if not value or is_missing(value):
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def parse_float(value: str) -> float | None:
    value = normalize_text(value)
    if not value or is_missing(value):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def accession_root(accession: str) -> str:
    accession = normalize_text(accession)
    if not accession:
        return ""
    if "." in accession:
        accession = accession.split(".", 1)[0]
    if accession.startswith("GCA_") or accession.startswith("GCF_"):
        return accession[4:]
    return accession


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1]


def fasta_iter(path: Path):
    name = None
    seq_chunks: list[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_chunks).upper()
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if name is not None:
            yield name, "".join(seq_chunks).upper()


def extract_contig_subseq(fasta_path: Path, contig_id: str, start_1: int, end_1: int) -> str | None:
    for name, sequence in fasta_iter(fasta_path):
        if name != contig_id:
            continue
        start_1 = max(1, start_1)
        end_1 = min(len(sequence), end_1)
        if end_1 < start_1:
            return None
        return sequence[start_1 - 1 : end_1]
    return None


def load_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"no rows found in {path}")
    return rows


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_country_iso3_map(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            iso3 = normalize_text(row.get("country_iso3", ""))
            if not iso3:
                continue
            for key in (
                row.get("raw_country_string", ""),
                row.get("normalized_lookup_key", ""),
                row.get("normalized_country_name", ""),
            ):
                normalized = normalize_text(key).casefold()
                if normalized:
                    mapping[normalized] = iso3
    return mapping


def choose_candidate(candidates: list[dict[str, str]], preferred_accession: str, accession_field: str) -> dict[str, str] | None:
    if not candidates:
        return None
    preferred_accession = normalize_text(preferred_accession)
    ranked = sorted(
        candidates,
        key=lambda row: (
            0 if normalize_text(row.get(accession_field, "")) == preferred_accession else 1,
            normalize_text(row.get(accession_field, "")),
        ),
    )
    return ranked[0]


def choose_prn_candidate(candidates: list[dict[str, str]], preferred_accession: str) -> dict[str, str] | None:
    if not candidates:
        return None
    preferred_accession = normalize_text(preferred_accession)
    informative_calls = {"intact", "disrupted_multi_hsp", "partial"}
    ranked = sorted(
        candidates,
        key=lambda row: (
            0 if normalize_text(row.get("prn_call", "")) in informative_calls else 1,
            0 if normalize_text(row.get("genome_resolved_accession", "")) == preferred_accession else 1,
            normalize_text(row.get("genome_resolved_accession", "")),
        ),
    )
    return ranked[0]


def load_reference_fasta(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    for header, sequence in fasta_iter(path):
        parts = header.split("|")
        reference_id = parts[0]
        is_element_name = parts[1] if len(parts) > 1 else reference_id
        source_accession = parts[2] if len(parts) > 2 else ""
        records.append(
            {
                "reference_id": reference_id,
                "is_element_name": is_element_name,
                "source_accession": source_accession,
                "sequence": sequence,
            }
        )
    if not records:
        raise ValueError(f"no reference sequences found in {path}")
    return records


def parse_gap_flank_header(header: str) -> dict[str, str]:
    parts = header.split("|")
    parsed = {
        "genome_resolved_accession": parts[0] if parts else "",
        "contig_id": "",
        "extract_start": "",
        "extract_end": "",
        "gap_start": "",
        "gap_end": "",
        "gap_len": "",
    }
    for part in parts[1:]:
        if ":" not in part:
            continue
        key, value = part.split(":", 1)
        if key in {"gap", "gaplen", "year", "country", "st"}:
            if key == "gap":
                if "-" in value:
                    start, end = value.split("-", 1)
                    parsed["gap_start"] = start
                    parsed["gap_end"] = end
            elif key == "gaplen":
                parsed["gap_len"] = value
        else:
            if "-" in value:
                start, end = value.split("-", 1)
                parsed["contig_id"] = key
                parsed["extract_start"] = start
                parsed["extract_end"] = end
    return parsed


def load_gap_flank_sequences(path: Path) -> dict[tuple[str, str, str, str], str]:
    sequences: dict[tuple[str, str, str, str], str] = {}
    for header, sequence in fasta_iter(path):
        parsed = parse_gap_flank_header(header)
        key = (
            normalize_text(parsed["genome_resolved_accession"]),
            normalize_text(parsed["contig_id"]),
            normalize_text(parsed["extract_start"]),
            normalize_text(parsed["extract_end"]),
        )
        if all(key):
            sequences[key] = sequence
    return sequences


def recover_gap_sequence(
    gap_row: dict[str, str] | None,
    gap_flank_sequences: dict[tuple[str, str, str, str], str],
) -> tuple[str, str]:
    if gap_row is None:
        return "", "missing_gap_metadata"

    contig_id = normalize_text(gap_row.get("contig_id", ""))
    extract_start = normalize_text(gap_row.get("extract_start", ""))
    extract_end = normalize_text(gap_row.get("extract_end", ""))
    extract_start_int = parse_int(extract_start)
    gap_start = parse_int(gap_row.get("gap_start", ""))
    gap_end = parse_int(gap_row.get("gap_end", ""))
    gap_flank_key = (
        normalize_text(gap_row.get("genome_resolved_accession", "")),
        contig_id,
        extract_start,
        extract_end,
    )
    flank_sequence = gap_flank_sequences.get(gap_flank_key, "")
    if flank_sequence and extract_start_int is not None and gap_start is not None and gap_end is not None:
        offset_start = gap_start - extract_start_int
        offset_end = gap_end - extract_start_int + 1
        if offset_start >= 0:
            gap_sequence = flank_sequence[offset_start:offset_end]
            if gap_sequence:
                return gap_sequence, "step3_gap_flank_fasta"

    fasta_path = Path(normalize_text(gap_row.get("genome_fasta_path", "")))
    if fasta_path.exists() and contig_id and gap_start is not None and gap_end is not None:
        extracted = extract_contig_subseq(fasta_path, contig_id, gap_start, gap_end)
        if extracted:
            return extracted, "legacy_genome_fasta"

    return "", "gap_sequence_source_missing"


def build_kmer_index(sequence: str, k: int) -> dict[str, list[int]]:
    index: dict[str, list[int]] = defaultdict(list)
    for i in range(0, len(sequence) - k + 1):
        index[sequence[i : i + k]].append(i)
    return index


def best_offset_alignment(query_seq: str, reference_seq: str, k: int) -> tuple[float, float]:
    if not query_seq or not reference_seq:
        return 0.0, 0.0
    if len(query_seq) < k or len(reference_seq) < k:
        overlap = min(len(query_seq), len(reference_seq))
        if overlap == 0:
            return 0.0, 0.0
        matches = sum(1 for a, b in zip(query_seq[:overlap], reference_seq[:overlap]) if a == b)
        return 100.0 * matches / overlap, 100.0 * overlap / len(query_seq)

    ref_index = build_kmer_index(reference_seq, k)
    offset_counter: Counter[int] = Counter()
    for q_pos in range(0, len(query_seq) - k + 1):
        kmer = query_seq[q_pos : q_pos + k]
        for r_pos in ref_index.get(kmer, []):
            offset_counter[r_pos - q_pos] += 1

    if not offset_counter:
        return 0.0, 0.0

    best_pid = 0.0
    best_qcov = 0.0
    for offset, _count in offset_counter.most_common(8):
        q_start = max(0, -offset)
        q_end = min(len(query_seq), len(reference_seq) - offset)
        overlap = q_end - q_start
        if overlap <= 0:
            continue
        matches = 0
        for q_index in range(q_start, q_end):
            r_index = q_index + offset
            if query_seq[q_index] == reference_seq[r_index]:
                matches += 1
        pid = 100.0 * matches / overlap
        qcov = 100.0 * overlap / len(query_seq)
        if (pid * qcov, pid, qcov) > (best_pid * best_qcov, best_pid, best_qcov):
            best_pid = pid
            best_qcov = qcov
    return best_pid, best_qcov


def scan_is_gap_sequence(query_seq: str, references: list[dict[str, str]]) -> dict[str, str]:
    hits = scan_is_gap_sequence_all(query_seq, references)
    if not hits:
        return {"is_element_best_hit": "", "is_element_best_hit_pident": "", "is_element_best_hit_qcov": ""}
    top_hit = hits[0]
    return {
        "is_element_best_hit": top_hit["is_element_name"],
        "is_element_best_hit_pident": top_hit["hit_pident"],
        "is_element_best_hit_qcov": top_hit["hit_qcov"],
    }


def scan_is_gap_sequence_all(query_seq: str, references: list[dict[str, str]]) -> list[dict[str, str]]:
    query_seq = normalize_text(query_seq).upper()
    if not query_seq:
        return []

    rc_query = reverse_complement(query_seq)
    k = max(7, min(15, min(len(query_seq), min(len(ref["sequence"]) for ref in references)) // 6 or 7))
    hits: list[dict[str, str]] = []

    for ref in references:
        best_pid = 0.0
        best_qcov = 0.0
        best_orientation = "+"
        for orientation, oriented_query in (("+", query_seq), ("-", rc_query)):
            pid, qcov = best_offset_alignment(oriented_query, ref["sequence"], k)
            if (pid * qcov, pid, qcov) > (best_pid * best_qcov, best_pid, best_qcov):
                best_pid = pid
                best_qcov = qcov
                best_orientation = orientation
        hits.append(
            {
                "reference_id": ref["reference_id"],
                "is_element_name": ref["is_element_name"],
                "reference_source_accession": ref["source_accession"],
                "hit_orientation": best_orientation,
                "hit_pident": f"{best_pid:.2f}",
                "hit_qcov": f"{best_qcov:.2f}",
                "hit_score": f"{best_pid * best_qcov:.2f}",
            }
        )

    hits.sort(
        key=lambda hit: (
            float(hit["hit_score"]),
            float(hit["hit_pident"]),
            float(hit["hit_qcov"]),
            hit["is_element_name"],
        ),
        reverse=True,
    )
    return hits


def confidence_from_call(mechanism_call: str, is_hit_qcov: float | None, is_hit_pident: float | None) -> str:
    if mechanism_call == "intact":
        return "rule_high"
    if mechanism_call in {"coding_disrupted_is481", "coding_disrupted_other_is"}:
        if (is_hit_qcov or 0.0) >= 90.0 and (is_hit_pident or 0.0) >= 95.0:
            return "rule_high"
        return "rule_medium"
    if mechanism_call == "coding_disrupted_inversion_or_rearrangement":
        return "rule_medium"
    if mechanism_call == "coding_disrupted_other":
        return "rule_medium"
    if mechanism_call == "uncertain_fragmented_assembly":
        return "rule_low"
    return "rule_low"


def build_event_id(row: dict[str, str]) -> str:
    mechanism = row["prn_mechanism_call"]
    if mechanism == "intact":
        return "prn_evt_intact"
    if mechanism == "insufficient_data":
        initial = normalize_text(row["prn_call_initial"]).casefold().replace(" ", "_")
        return f"prn_evt_insufficient__{initial or 'unspecified'}"
    if mechanism == "uncertain_fragmented_assembly":
        return "prn_evt_fragmented_contigs"
    if mechanism in {"coding_disrupted_is481", "coding_disrupted_other_is"}:
        hit = normalize_text(row["is_element_best_hit"]).lower() or "unknown_is"
        gap = normalize_text(row["insertion_subject_gap_bp"]) or "na"
        return f"prn_evt_{mechanism}__{hit}__gap{gap}"
    if mechanism == "coding_disrupted_inversion_or_rearrangement":
        cov = normalize_text(row["prn_best_single_cov_pct"]).split(".")[0] or "na"
        return f"prn_evt_rearrangement__{normalize_text(row['bp_category']) or 'na'}__cov{cov}"
    if mechanism == "coding_disrupted_other":
        gap = normalize_text(row["insertion_subject_gap_bp"]) or "na"
        return f"prn_evt_other_disruption__{normalize_text(row['bp_category']) or 'na'}__gap{gap}"
    return f"prn_evt_{mechanism}"


def classify_row(
    qc_row: dict[str, str],
    prn_row: dict[str, str] | None,
    bp_row: dict[str, str] | None,
    gap_row: dict[str, str] | None,
    *,
    country_iso3_map: dict[str, str],
    references: list[dict[str, str]],
    gap_flank_sequences: dict[tuple[str, str, str, str], str],
) -> dict[str, str]:
    country = normalize_text(qc_row.get("country", ""))
    country_iso3 = country_iso3_map.get(country.casefold(), "")
    evidence_flags: list[str] = []
    notes: list[str] = ["confidence_is_pre_prn04_rule_tier"]

    if prn_row is None:
        evidence_flags.extend(["no_step3_prn_input", "legacy_step3_coverage_gap"])
        mechanism_call = "insufficient_data"
        classified = {
            "sample_id_canonical": qc_row.get("sample_id_canonical", ""),
            "biosample_accession": qc_row.get("biosample_accession", ""),
            "assembly_accession": qc_row.get("current_accession", ""),
            "sra_run_accession": qc_row.get("sra_run_accession", ""),
            "country_iso3": country_iso3,
            "year": qc_row.get("year", ""),
            "mlst_st": "",
            "phylo_lineage": "",
            "prn_call_initial": "not_available_current_step3",
            "prn_mechanism_call": mechanism_call,
            "prn_call_confidence": "rule_low",
            "prn_event_id": "",
            "prn_query_cov_pct": "",
            "prn_best_single_cov_pct": "",
            "prn_hsp_n": "",
            "bp_category": "",
            "insertion_subject_gap_bp": "",
            "is_element_best_hit": "",
            "is_element_best_hit_pident": "",
            "is_element_best_hit_qcov": "",
            "read_validation_status": "not_run",
            "read_validation_support": "not_evaluated",
            "evidence_flags": ";".join(sorted(set(evidence_flags))),
            "notes": ";".join(notes),
        }
        classified["prn_event_id"] = build_event_id(classified)
        return classified

    if normalize_text(prn_row.get("genome_resolved_accession", "")) != normalize_text(qc_row.get("current_accession", "")):
        evidence_flags.append("mapped_from_legacy_mirror_or_alt_accession")
        notes.append(f"step3_source_accession={normalize_text(prn_row.get('genome_resolved_accession', ''))}")

    prn_call_initial = normalize_text(prn_row.get("prn_call", ""))
    mechanism_call = "insufficient_data"
    bp_category = normalize_text(bp_row.get("bp_category", "")) if bp_row else ""
    is_hit = {"is_element_best_hit": "", "is_element_best_hit_pident": "", "is_element_best_hit_qcov": ""}

    if prn_call_initial == "intact":
        mechanism_call = "intact"
        evidence_flags.append("prn_intact_call")
    elif prn_call_initial == "missing_fasta":
        mechanism_call = "insufficient_data"
        evidence_flags.extend(["missing_fasta", "assembly_sequence_unavailable"])
    elif prn_call_initial == "partial":
        mechanism_call = "insufficient_data"
        evidence_flags.extend(["partial_prn_call", "no_step4_structural_upgrade"])
    elif prn_call_initial == "disrupted_multi_hsp":
        evidence_flags.append("disrupted_multi_hsp")
        if bp_row is None:
            mechanism_call = "insufficient_data"
            evidence_flags.append("missing_breakpoint_evidence")
        else:
            if bp_category == "fragmented_contigs" or bp_category == "near_contig_end":
                mechanism_call = "uncertain_fragmented_assembly"
                evidence_flags.append("bp_fragmentation_signal")
            elif bp_category == "within_contig":
                mechanism_call = "coding_disrupted_inversion_or_rearrangement"
                evidence_flags.extend(["bp_within_contig", "bp_opposite_strand_or_complex"])
            elif bp_category == "insertion_like":
                evidence_flags.append("bp_insertion_like")
                gap_sequence = ""
                if gap_row is not None:
                    gap_len = parse_int(gap_row.get("gap_len", ""))
                    gap_sequence, gap_source = recover_gap_sequence(gap_row, gap_flank_sequences)
                    if gap_sequence:
                        evidence_flags.append("gap_sequence_extracted")
                    else:
                        evidence_flags.append(gap_source)
                    if gap_len is not None and gap_len >= 50:
                        evidence_flags.append("insertion_gap_ge_50bp")
                else:
                    evidence_flags.append("missing_gap_metadata")

                if gap_sequence:
                    is_hit = scan_is_gap_sequence(gap_sequence, references)
                    hit_qcov = parse_float(is_hit["is_element_best_hit_qcov"])
                    hit_pid = parse_float(is_hit["is_element_best_hit_pident"])
                    if hit_qcov is not None and hit_pid is not None and hit_qcov >= 80.0 and hit_pid >= 90.0:
                        evidence_flags.append("strong_is_hit")
                        if normalize_text(is_hit["is_element_best_hit"]) == "IS481":
                            mechanism_call = "coding_disrupted_is481"
                            evidence_flags.append("best_hit_is481")
                        else:
                            mechanism_call = "coding_disrupted_other_is"
                            evidence_flags.append("best_hit_other_is")
                    elif hit_qcov is not None and hit_pid is not None and hit_qcov >= 60.0 and hit_pid >= 85.0:
                        evidence_flags.append("moderate_is_hit")
                        if normalize_text(is_hit["is_element_best_hit"]) == "IS481":
                            mechanism_call = "coding_disrupted_is481"
                            evidence_flags.append("best_hit_is481")
                        else:
                            mechanism_call = "coding_disrupted_other_is"
                            evidence_flags.append("best_hit_other_is")
                    else:
                        mechanism_call = "coding_disrupted_other"
                        evidence_flags.append("no_supported_is_hit")
                else:
                    mechanism_call = "coding_disrupted_other"
                    evidence_flags.append("no_gap_sequence_for_is_screen")
            else:
                mechanism_call = "coding_disrupted_other"
                evidence_flags.append("unclassified_breakpoint_pattern")
    else:
        mechanism_call = "insufficient_data"
        evidence_flags.append("unknown_prn_initial_call")

    hit_qcov = parse_float(is_hit["is_element_best_hit_qcov"])
    hit_pid = parse_float(is_hit["is_element_best_hit_pident"])
    classified = {
        "sample_id_canonical": qc_row.get("sample_id_canonical", ""),
        "biosample_accession": qc_row.get("biosample_accession", ""),
        "assembly_accession": qc_row.get("current_accession", ""),
        "sra_run_accession": qc_row.get("sra_run_accession", ""),
        "country_iso3": country_iso3,
        "year": qc_row.get("year", ""),
        "mlst_st": normalize_text(prn_row.get("mlst_st", "")) or qc_row.get("mlst_st", ""),
        "phylo_lineage": "",
        "prn_call_initial": prn_call_initial,
        "prn_mechanism_call": mechanism_call,
        "prn_call_confidence": confidence_from_call(mechanism_call, hit_qcov, hit_pid),
        "prn_event_id": "",
        "prn_query_cov_pct": normalize_text(prn_row.get("prn_qcov_union_pct", "")),
        "prn_best_single_cov_pct": normalize_text(prn_row.get("prn_best_single_qcov_pct", "")),
        "prn_hsp_n": normalize_text(prn_row.get("prn_hsp_n_contrib", "")) or normalize_text(prn_row.get("prn_hsp_n", "")),
        "bp_category": bp_category,
        "insertion_subject_gap_bp": normalize_text(bp_row.get("bp_max_subject_gap", "")) if bp_row else "",
        "is_element_best_hit": is_hit["is_element_best_hit"],
        "is_element_best_hit_pident": is_hit["is_element_best_hit_pident"],
        "is_element_best_hit_qcov": is_hit["is_element_best_hit_qcov"],
        "read_validation_status": "not_run",
        "read_validation_support": "not_evaluated",
        "evidence_flags": ";".join(sorted(set(flag for flag in evidence_flags if flag))),
        "notes": ";".join(dict.fromkeys(note for note in notes if note)),
    }
    classified["prn_event_id"] = build_event_id(classified)
    return classified


def build_event_catalog(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["prn_event_id"]].append(row)

    catalog_rows: list[dict[str, str]] = []
    for event_id, members in sorted(grouped.items(), key=lambda item: (-len(item[1]), item[0])):
        years = [parse_int(row["year"]) for row in members if parse_int(row["year"]) is not None]
        countries = sorted({row["country_iso3"] for row in members if normalize_text(row["country_iso3"])})
        assemblies = sorted({row["assembly_accession"] for row in members if normalize_text(row["assembly_accession"])})
        catalog_rows.append(
            {
                "prn_event_id": event_id,
                "prn_mechanism_call": members[0]["prn_mechanism_call"],
                "prn_call_confidence_mode": Counter(row["prn_call_confidence"] for row in members).most_common(1)[0][0],
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
                "notes": members[0]["notes"],
            }
        )
    return catalog_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Upgrade legacy step3 prn calls into canonical sample-level mechanism calls using "
            "breakpoint heuristics and the PRN-01 IS reference set."
        )
    )
    parser.add_argument(
        "--qc-manifest",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_public_genome_qc_manifest.tsv",
        help="Retained canonical sample manifest from GC-03.",
    )
    parser.add_argument(
        "--prn-calls",
        type=Path,
        default=STEP3_DATA_ROOT / "outputs" / "bp_prn_disruption_calls.tsv",
        help="Legacy step3 coarse prn calls TSV.",
    )
    parser.add_argument(
        "--breakpoint-evidence",
        type=Path,
        default=STEP3_DATA_ROOT / "outputs" / "bp_prn_breakpoint_evidence.tsv",
        help="Legacy step3 breakpoint evidence TSV.",
    )
    parser.add_argument(
        "--gap-metadata",
        type=Path,
        default=STEP3_DATA_ROOT / "outputs" / "bp_prn_insertion_gap_plus_flanks.tsv",
        help="Step3 gap extraction metadata TSV used to recover insertion-only sequences.",
    )
    parser.add_argument(
        "--gap-flank-fasta",
        type=Path,
        default=STEP3_DATA_ROOT / "outputs" / "bp_prn_insertion_gap_plus_flanks.fasta",
        help="Step3 FASTA of extracted gap+flank sequences.",
    )
    parser.add_argument(
        "--is-reference-fasta",
        type=Path,
        default=STEP4_DATA_ROOT / "references" / "is_elements" / "bp_is_reference.fasta",
        help="PRN-01 reference FASTA.",
    )
    parser.add_argument(
        "--country-map",
        type=Path,
        default=PUBLIC_HEALTH_DATA_ROOT / "outputs" / "ph_country_name_map.tsv",
        help="Country normalization table used to map country names to ISO3.",
    )
    parser.add_argument(
        "--mechanism-out",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_mechanism_calls.tsv",
        help="Main mechanism calls output TSV.",
    )
    parser.add_argument(
        "--event-out",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_event_catalog.tsv",
        help="Event catalog output TSV.",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    qc_rows = load_tsv_rows(args.qc_manifest)
    prn_rows = load_tsv_rows(args.prn_calls)
    bp_rows = load_tsv_rows(args.breakpoint_evidence)
    gap_rows = load_tsv_rows(args.gap_metadata)
    gap_flank_sequences = load_gap_flank_sequences(args.gap_flank_fasta)
    references = load_reference_fasta(args.is_reference_fasta)
    country_iso3_map = load_country_iso3_map(args.country_map)

    prn_by_root: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in prn_rows:
        prn_by_root[accession_root(row.get("genome_resolved_accession", ""))].append(row)

    bp_by_root: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in bp_rows:
        bp_by_root[accession_root(row.get("genome_resolved_accession", ""))].append(row)

    gap_by_root: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in gap_rows:
        gap_by_root[accession_root(row.get("genome_resolved_accession", ""))].append(row)

    mechanism_rows: list[dict[str, str]] = []
    for qc_row in qc_rows:
        preferred_accession = qc_row.get("current_accession", "")
        root = qc_row.get("assembly_accession_root", "") or accession_root(preferred_accession)
        prn_row = choose_prn_candidate(prn_by_root.get(root, []), preferred_accession)
        bp_row = choose_candidate(bp_by_root.get(root, []), preferred_accession, "genome_resolved_accession")
        gap_row = choose_candidate(gap_by_root.get(root, []), preferred_accession, "genome_resolved_accession")
        mechanism_rows.append(
            classify_row(
                qc_row,
                prn_row,
                bp_row,
                gap_row,
                country_iso3_map=country_iso3_map,
                references=references,
                gap_flank_sequences=gap_flank_sequences,
            )
        )

    event_rows = build_event_catalog(mechanism_rows)
    write_tsv(args.mechanism_out, MECHANISM_COLUMNS, mechanism_rows)
    write_tsv(args.event_out, EVENT_COLUMNS, event_rows)

    counts = Counter(row["prn_mechanism_call"] for row in mechanism_rows)
    print(f"Wrote mechanism calls: {args.mechanism_out}")
    print(f"Wrote event catalog: {args.event_out}")
    print(f"Retained samples processed: {len(mechanism_rows)}")
    print("Mechanism counts: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    print(f"Unique prn events: {len(event_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
