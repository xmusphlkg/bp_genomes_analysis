#!/usr/bin/env python3
"""Run a full prn caller threshold grid from cached broad BLAST HSP evidence."""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import os
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile


ROOT = Path(__file__).resolve().parents[3]
DATA_HOME = Path(os.environ.get("PERTUSSIS_PROJECT_DATA_ROOT", ROOT / "pertussis_data" / "pertussis_gene"))
FIGURE_DATA = ROOT / "manuscript" / "figure_data"
SUPPLEMENTARY = ROOT / "manuscript" / "supplementary"
AUDIT_LEDGER_DIR = ROOT / "manuscript" / "submission_data" / "audit_ledgers" / "supplementary_table_sources"

QC_MANIFEST = DATA_HOME / "step1_ingest" / "outputs" / "bp_public_genome_qc_manifest.tsv"
GENOME_PATHS = DATA_HOME / "step2_typing" / "outputs" / "bp_genome_paths_qc.tsv"
RAW_GENOME_PATHS = DATA_HOME / "step1_ingest" / "outputs" / "bp_raw_read_step3_genome_paths.tsv"
PRN_QUERY = ROOT / "modules" / "step2_typing" / "refs" / "markers" / "prn_maker.fasta"
IS_REFERENCE = ROOT / "modules" / "step4_prn_validation" / "refs" / "is_elements" / "bp_is_reference.fasta"
MANUSCRIPT_ANNOTATION = FIGURE_DATA / "published_overlap_annotation.tsv"

GRID_DIR = DATA_HOME / "step4_prn_validation" / "threshold_grid"
HSP_CACHE = GRID_DIR / "prn_blast_hsps_min80.tsv"
SUMMARY_OUT = FIGURE_DATA / "prn_threshold_grid_full.tsv"
EVENT_OUT = FIGURE_DATA / "prn_threshold_grid_top_events.tsv"
SUPPLEMENTARY_OUT = AUDIT_LEDGER_DIR / "Supplementary_Table_67_prn_Threshold_Grid_Sensitivity.tsv"

HSP_PIDENT_GRID = [85.0, 88.0, 90.0, 92.0, 95.0]
LOCUS_QCOV_GRID = [90.0, 92.5, 95.0, 97.5, 99.0]
IS_SUPPORT_PROFILES = [
    {
        "is_support_profile": "relaxed",
        "strong_is_qcov": 70.0,
        "strong_is_pident": 85.0,
        "moderate_is_qcov": 50.0,
        "moderate_is_pident": 80.0,
    },
    {
        "is_support_profile": "default",
        "strong_is_qcov": 80.0,
        "strong_is_pident": 90.0,
        "moderate_is_qcov": 60.0,
        "moderate_is_pident": 85.0,
    },
    {
        "is_support_profile": "strict",
        "strong_is_qcov": 90.0,
        "strong_is_pident": 95.0,
        "moderate_is_qcov": 70.0,
        "moderate_is_pident": 90.0,
    },
]

DISRUPTED_MECHANISMS = {
    "coding_disrupted_is481",
    "coding_disrupted_other_is",
    "coding_disrupted_inversion_or_rearrangement",
    "coding_disrupted_other",
}
INTERPRETABLE_MECHANISMS = DISRUPTED_MECHANISMS | {"intact"}


@dataclass(frozen=True)
class Hsp:
    sample_id: str
    accession: str
    sseqid: str
    pident: float
    length: int
    qlen: int
    qstart: int
    qend: int
    slen: int
    sstart: int
    send: int
    bitscore: float

    @property
    def q_a(self) -> int:
        return min(self.qstart, self.qend)

    @property
    def q_b(self) -> int:
        return max(self.qstart, self.qend)

    @property
    def s_a(self) -> int:
        return min(self.sstart, self.send)

    @property
    def s_b(self) -> int:
        return max(self.sstart, self.send)

    @property
    def same_strand(self) -> bool:
        return (self.qend - self.qstart) * (self.send - self.sstart) > 0


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", newline="", encoding="utf-8", delete=False, dir=path.parent) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
        temp_path = Path(handle.name)
    temp_path.replace(path)


def normalize(value: str | None) -> str:
    return (value or "").strip()


def parse_float(value: str | None) -> float | None:
    value = normalize(value)
    if value == "" or value.upper() == "NA":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_int(value: str | None) -> int | None:
    parsed = parse_float(value)
    return None if parsed is None else int(parsed)


def fasta_iter(path: Path):
    header = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks).upper()


def load_reference_fasta(path: Path) -> list[dict[str, str]]:
    refs: list[dict[str, str]] = []
    for header, sequence in fasta_iter(path):
        parts = header.split("|")
        refs.append(
            {
                "reference_id": parts[0],
                "is_element_name": parts[1] if len(parts) > 1 else parts[0],
                "sequence": sequence,
            }
        )
    if not refs:
        raise ValueError(f"No IS reference records in {path}")
    return refs


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1]


def build_kmer_index(sequence: str, k: int) -> dict[str, list[int]]:
    index: dict[str, list[int]] = defaultdict(list)
    for idx in range(0, len(sequence) - k + 1):
        index[sequence[idx : idx + k]].append(idx)
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
    offsets: Counter[int] = Counter()
    for q_pos in range(0, len(query_seq) - k + 1):
        kmer = query_seq[q_pos : q_pos + k]
        for r_pos in ref_index.get(kmer, []):
            offsets[r_pos - q_pos] += 1

    best_pid = 0.0
    best_qcov = 0.0
    for offset, _count in offsets.most_common(8):
        q_start = max(0, -offset)
        q_end = min(len(query_seq), len(reference_seq) - offset)
        overlap = q_end - q_start
        if overlap <= 0:
            continue
        matches = 0
        for q_index in range(q_start, q_end):
            if query_seq[q_index] == reference_seq[q_index + offset]:
                matches += 1
        pid = 100.0 * matches / overlap
        qcov = 100.0 * overlap / len(query_seq)
        if (pid * qcov, pid, qcov) > (best_pid * best_qcov, best_pid, best_qcov):
            best_pid = pid
            best_qcov = qcov
    return best_pid, best_qcov


def scan_is_gap_sequence(query_seq: str, references: list[dict[str, str]]) -> dict[str, str]:
    query_seq = normalize(query_seq).upper()
    if not query_seq:
        return {"is_element_best_hit": "", "is_element_best_hit_pident": "", "is_element_best_hit_qcov": ""}
    rc_query = reverse_complement(query_seq)
    k = max(7, min(15, min(len(query_seq), min(len(ref["sequence"]) for ref in references)) // 6 or 7))
    hits: list[dict[str, str]] = []
    for ref in references:
        best_pid = 0.0
        best_qcov = 0.0
        for oriented_query in (query_seq, rc_query):
            pid, qcov = best_offset_alignment(oriented_query, ref["sequence"], k)
            if (pid * qcov, pid, qcov) > (best_pid * best_qcov, best_pid, best_qcov):
                best_pid = pid
                best_qcov = qcov
        hits.append(
            {
                "is_element_best_hit": ref["is_element_name"],
                "is_element_best_hit_pident": f"{best_pid:.2f}",
                "is_element_best_hit_qcov": f"{best_qcov:.2f}",
                "score": best_pid * best_qcov,
            }
        )
    hits.sort(key=lambda hit: (hit["score"], hit["is_element_best_hit"]), reverse=True)
    top = hits[0]
    return {
        "is_element_best_hit": top["is_element_best_hit"],
        "is_element_best_hit_pident": top["is_element_best_hit_pident"],
        "is_element_best_hit_qcov": top["is_element_best_hit_qcov"],
    }


def extract_contig_subseq(fasta_path: Path, contig_id: str, start_1: int, end_1: int) -> str:
    for name, sequence in fasta_iter(fasta_path):
        if name != contig_id:
            continue
        start_1 = max(1, start_1)
        end_1 = min(len(sequence), end_1)
        if end_1 < start_1:
            return ""
        return sequence[start_1 - 1 : end_1]
    return ""


def load_genome_paths(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    out: dict[str, str] = {}
    for row in rows:
        if normalize(row.get("status")) != "ok":
            continue
        accession = normalize(row.get("resolved_accession"))
        fasta = normalize(row.get("fasta_path"))
        if accession and fasta:
            out[accession] = resolve_path_text(fasta)
    return out


def resolve_path_text(path_text: str) -> str:
    path = Path(path_text)
    if path.is_absolute():
        return str(path)
    for candidate in (ROOT / path, DATA_HOME / path, Path.cwd() / path):
        if candidate.exists():
            return str(candidate)
    return str(ROOT / path)


def load_raw_manuscript_records(path: Path, raw_genome_paths: dict[str, str]) -> list[dict[str, str]]:
    if not path.exists():
        return []
    records: list[dict[str, str]] = []
    for row in read_tsv(path):
        accession = normalize(row.get("assembly_accession"))
        if not accession.startswith("RRASM_"):
            continue
        fasta = raw_genome_paths.get(accession, "")
        if not fasta or not Path(fasta).exists():
            continue
        records.append(
            {
                "sample_id": normalize(row.get("sample_id_canonical")) or normalize(row.get("biosample_accession")) or accession,
                "biosample_accession": normalize(row.get("biosample_accession")),
                "assembly_accession": accession,
                "fasta_path": fasta,
                "country": normalize(row.get("country_iso3")),
                "year": normalize(row.get("year")),
            }
        )
    return records


def load_qc_records(path: Path, genome_paths: dict[str, str], max_genomes: int | None) -> list[dict[str, str]]:
    rows = read_tsv(path)
    records: list[dict[str, str]] = []
    for row in rows:
        accession = normalize(row.get("current_accession"))
        if not accession:
            continue
        fasta = genome_paths.get(accession, "")
        records.append(
            {
                "sample_id": normalize(row.get("sample_id_canonical")) or normalize(row.get("biosample_accession")) or accession,
                "biosample_accession": normalize(row.get("biosample_accession")),
                "assembly_accession": accession,
                "fasta_path": fasta,
                "country": normalize(row.get("country")),
                "year": normalize(row.get("year")),
            }
        )
    return records


def blast_one(record: dict[str, str], query: Path, blastn: str, min_cache_pident: float, max_targets: int) -> list[dict[str, object]]:
    fasta = Path(record["fasta_path"])
    if not fasta.exists() or fasta.stat().st_size == 0:
        return []
    outfmt = "6 pident length qlen qstart qend sseqid slen sstart send bitscore"
    cmd = [
        blastn,
        "-query",
        str(query),
        "-subject",
        str(fasta),
        "-outfmt",
        outfmt,
        "-task",
        "blastn",
        "-perc_identity",
        str(min_cache_pident),
        "-max_target_seqs",
        str(max_targets),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError((proc.stderr or "").strip() or "blastn failed")
    rows: list[dict[str, object]] = []
    for line in (proc.stdout or "").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 10:
            continue
        rows.append(
            {
                "sample_id": record["sample_id"],
                "assembly_accession": record["assembly_accession"],
                "sseqid": parts[5],
                "pident": parts[0],
                "length": parts[1],
                "qlen": parts[2],
                "qstart": parts[3],
                "qend": parts[4],
                "slen": parts[6],
                "sstart": parts[7],
                "send": parts[8],
                "bitscore": parts[9],
            }
        )
    return rows


def build_hsp_cache(
    records: list[dict[str, str]],
    cache_path: Path,
    query: Path,
    blastn: str,
    min_cache_pident: float,
    jobs: int,
    max_targets: int,
    force: bool,
) -> None:
    if cache_path.exists() and cache_path.stat().st_size > 0 and not force:
        print(f"[threshold-grid] reusing HSP cache: {cache_path}", flush=True)
        return
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample_id",
        "assembly_accession",
        "sseqid",
        "pident",
        "length",
        "qlen",
        "qstart",
        "qend",
        "slen",
        "sstart",
        "send",
        "bitscore",
    ]
    rows_written = 0
    with NamedTemporaryFile("w", newline="", encoding="utf-8", delete=False, dir=cache_path.parent) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        with cf.ThreadPoolExecutor(max_workers=max(1, jobs)) as executor:
            futures = [
                executor.submit(blast_one, record, query, blastn, min_cache_pident, max_targets)
                for record in records
            ]
            for idx, future in enumerate(cf.as_completed(futures), start=1):
                for row in future.result():
                    writer.writerow(row)
                    rows_written += 1
                if idx % 100 == 0:
                    print(f"[threshold-grid] BLAST finished {idx}/{len(records)}", flush=True)
        temp_path = Path(handle.name)
    temp_path.replace(cache_path)
    print(f"[threshold-grid] wrote {rows_written} HSP rows: {cache_path}", flush=True)


def load_hsps(path: Path) -> dict[str, list[Hsp]]:
    grouped: dict[str, list[Hsp]] = defaultdict(list)
    for row in read_tsv(path):
        try:
            hsp = Hsp(
                sample_id=row["sample_id"],
                accession=row["assembly_accession"],
                sseqid=row["sseqid"],
                pident=float(row["pident"]),
                length=int(float(row["length"])),
                qlen=int(float(row["qlen"])),
                qstart=int(float(row["qstart"])),
                qend=int(float(row["qend"])),
                slen=int(float(row["slen"])),
                sstart=int(float(row["sstart"])),
                send=int(float(row["send"])),
                bitscore=float(row["bitscore"]),
            )
        except (KeyError, ValueError):
            continue
        grouped[hsp.sample_id].append(hsp)
    return grouped


def interval_len(a: int, b: int) -> int:
    a2, b2 = sorted((a, b))
    return max(0, b2 - a2 + 1)


def greedy_contributing_hsps(hsps: list[Hsp], qcov_target: float) -> tuple[list[Hsp], float]:
    if not hsps:
        return [], 0.0
    qlen = hsps[0].qlen
    if qlen <= 0:
        return [], 0.0
    intervals: list[tuple[int, int]] = []
    contributing: list[Hsp] = []

    def merge(intervals_in: list[tuple[int, int]]) -> list[tuple[int, int]]:
        merged: list[tuple[int, int]] = []
        for a, b in sorted(intervals_in):
            if not merged or a > merged[-1][1] + 1:
                merged.append((a, b))
            else:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        return merged

    def cov_len(intervals_in: list[tuple[int, int]]) -> int:
        return sum(b - a + 1 for a, b in intervals_in)

    for hsp in sorted(hsps, key=lambda h: (h.bitscore, h.pident, h.length), reverse=True):
        before = cov_len(intervals)
        intervals2 = merge(intervals + [(hsp.q_a, hsp.q_b)])
        after = cov_len(intervals2)
        if after > before:
            intervals = intervals2
            contributing.append(hsp)
        if 100.0 * after / qlen >= qcov_target:
            break
    return contributing, 100.0 * cov_len(intervals) / qlen


def classify_breakpoint(contrib: list[Hsp]) -> dict[str, object]:
    if not contrib:
        return {"bp_category": "no_contrib"}
    best = max(contrib, key=lambda h: (h.bitscore, h.pident, h.length))
    best_single_qcov = 100.0 * interval_len(best.qstart, best.qend) / best.qlen if best.qlen > 0 else 0.0
    contigs = sorted({h.sseqid for h in contrib})
    end_dists = [min(h.s_a - 1, h.slen - h.s_b) for h in contrib if h.slen > 0]
    min_end_dist = min(end_dists) if end_dists else None
    ordered = sorted(contrib, key=lambda h: (h.q_a, h.q_b, -h.bitscore))
    max_subject_gap = None
    max_query_gap = None
    max_gap_pair = None
    same_contig_pairs = 0
    same_strand_pairs = 0
    for left, right in zip(ordered, ordered[1:]):
        q_gap = max(0, right.q_a - left.q_b - 1)
        if left.sseqid == right.sseqid:
            same_contig_pairs += 1
            s_gap = max(0, right.s_a - left.s_b - 1)
            if left.same_strand and right.same_strand:
                same_strand_pairs += 1
            if max_subject_gap is None or s_gap > max_subject_gap:
                max_subject_gap = s_gap
                max_gap_pair = (left, right)
            max_query_gap = q_gap if max_query_gap is None else max(max_query_gap, q_gap)
        else:
            max_query_gap = q_gap if max_query_gap is None else max(max_query_gap, q_gap)

    if len(contigs) >= 2:
        category = "fragmented_contigs"
    elif min_end_dist is not None and min_end_dist <= 200:
        category = "near_contig_end"
    else:
        category = "within_contig"
    if max_subject_gap is not None:
        qg = max_query_gap or 0
        if max_subject_gap >= 50 and (qg <= 10 or max_subject_gap >= 5 * max(1, qg)):
            category = "insertion_like"

    gap_start = None
    gap_end = None
    if category == "insertion_like" and max_gap_pair is not None:
        left, right = max_gap_pair
        gap_start = left.s_b + 1
        gap_end = right.s_a - 1
        if gap_end < gap_start:
            gap_start = None
            gap_end = None
    return {
        "bp_category": category,
        "bp_best_single_qcov_pct": best_single_qcov,
        "bp_contrib_hsp_n": len(contrib),
        "bp_contig_id": contigs[0] if contigs else "",
        "bp_gap_start": gap_start,
        "bp_gap_end": gap_end,
        "bp_max_subject_gap": max_subject_gap,
        "bp_max_query_gap": max_query_gap,
        "bp_same_contig_pairs": same_contig_pairs,
        "bp_same_strand_pairs": same_strand_pairs,
    }


def build_event_id(mechanism: str, bp: dict[str, object], is_hit: dict[str, str], best_single_qcov: float) -> str:
    if mechanism == "intact":
        return "prn_evt_intact"
    if mechanism == "insufficient_data":
        return "prn_evt_insufficient__threshold_grid"
    if mechanism == "uncertain_fragmented_assembly":
        return "prn_evt_fragmented_contigs"
    if mechanism in {"coding_disrupted_is481", "coding_disrupted_other_is"}:
        hit = normalize(is_hit.get("is_element_best_hit")).lower() or "unknown_is"
        gap = str(bp.get("bp_max_subject_gap") or "na")
        return f"prn_evt_{mechanism}__{hit}__gap{gap}"
    if mechanism == "coding_disrupted_inversion_or_rearrangement":
        return f"prn_evt_rearrangement__{bp.get('bp_category') or 'na'}__cov{int(best_single_qcov)}"
    if mechanism == "coding_disrupted_other":
        gap = str(bp.get("bp_max_subject_gap") or "na")
        return f"prn_evt_other_disruption__{bp.get('bp_category') or 'na'}__gap{gap}"
    return f"prn_evt_{mechanism}"


def classify_sample(
    record: dict[str, str],
    hsps_all: list[Hsp],
    *,
    hsp_min_pident: float,
    locus_qcov_threshold: float,
    is_profile: dict[str, float | str],
    references: list[dict[str, str]],
    gap_cache: dict[tuple[str, str, int, int], str],
    is_cache: dict[tuple[str, str, int, int], dict[str, str]],
) -> dict[str, str]:
    hsps = [hsp for hsp in hsps_all if hsp.pident >= hsp_min_pident]
    if not hsps:
        return {
            "sample_id": record["sample_id"],
            "assembly_accession": record["assembly_accession"],
            "prn_call_initial": "no_hit",
            "prn_mechanism_call": "insufficient_data",
            "prn_event_id": "prn_evt_insufficient__no_hit",
            "is_support_tier": "none",
        }
    best = max(hsps, key=lambda h: (h.bitscore, h.pident, h.length))
    best_single_qcov = 100.0 * interval_len(best.qstart, best.qend) / best.qlen if best.qlen > 0 else 0.0
    contrib, union_qcov = greedy_contributing_hsps(hsps, locus_qcov_threshold)

    if best_single_qcov >= locus_qcov_threshold:
        mechanism = "intact"
        initial = "intact"
        bp = {"bp_category": ""}
        is_hit = {"is_element_best_hit": "", "is_element_best_hit_pident": "", "is_element_best_hit_qcov": ""}
        is_support = "none"
    elif union_qcov >= locus_qcov_threshold and len(contrib) > 1:
        initial = "disrupted_multi_hsp"
        bp = classify_breakpoint(contrib)
        category = str(bp.get("bp_category") or "")
        is_hit = {"is_element_best_hit": "", "is_element_best_hit_pident": "", "is_element_best_hit_qcov": ""}
        is_support = "none"
        if category in {"fragmented_contigs", "near_contig_end"}:
            mechanism = "uncertain_fragmented_assembly"
        elif category == "within_contig":
            mechanism = "coding_disrupted_inversion_or_rearrangement"
        elif category == "insertion_like":
            contig_id = str(bp.get("bp_contig_id") or "")
            gap_start = bp.get("bp_gap_start")
            gap_end = bp.get("bp_gap_end")
            gap_seq = ""
            if isinstance(gap_start, int) and isinstance(gap_end, int) and contig_id:
                gap_key = (record["assembly_accession"], contig_id, gap_start, gap_end)
                if gap_key not in gap_cache:
                    gap_cache[gap_key] = extract_contig_subseq(Path(record["fasta_path"]), contig_id, gap_start, gap_end)
                gap_seq = gap_cache[gap_key]
                if gap_seq:
                    if gap_key not in is_cache:
                        is_cache[gap_key] = scan_is_gap_sequence(gap_seq, references)
                    is_hit = is_cache[gap_key]
            hit_qcov = parse_float(is_hit.get("is_element_best_hit_qcov"))
            hit_pid = parse_float(is_hit.get("is_element_best_hit_pident"))
            if hit_qcov is not None and hit_pid is not None:
                if (
                    hit_qcov >= float(is_profile["strong_is_qcov"])
                    and hit_pid >= float(is_profile["strong_is_pident"])
                ):
                    is_support = "strong"
                elif (
                    hit_qcov >= float(is_profile["moderate_is_qcov"])
                    and hit_pid >= float(is_profile["moderate_is_pident"])
                ):
                    is_support = "moderate"
            if is_support in {"strong", "moderate"}:
                if normalize(is_hit.get("is_element_best_hit")) == "IS481":
                    mechanism = "coding_disrupted_is481"
                else:
                    mechanism = "coding_disrupted_other_is"
            else:
                mechanism = "coding_disrupted_other"
        else:
            mechanism = "coding_disrupted_other"
    elif union_qcov >= 50.0:
        initial = "partial"
        mechanism = "insufficient_data"
        bp = {"bp_category": ""}
        is_hit = {"is_element_best_hit": "", "is_element_best_hit_pident": "", "is_element_best_hit_qcov": ""}
        is_support = "none"
    else:
        initial = "low_coverage"
        mechanism = "insufficient_data"
        bp = {"bp_category": ""}
        is_hit = {"is_element_best_hit": "", "is_element_best_hit_pident": "", "is_element_best_hit_qcov": ""}
        is_support = "none"

    return {
        "sample_id": record["sample_id"],
        "assembly_accession": record["assembly_accession"],
        "prn_call_initial": initial,
        "prn_mechanism_call": mechanism,
        "prn_event_id": build_event_id(mechanism, bp, is_hit, best_single_qcov),
        "prn_query_cov_pct": f"{union_qcov:.2f}",
        "prn_best_single_cov_pct": f"{best_single_qcov:.2f}",
        "prn_hsp_n": str(len(contrib)),
        "bp_category": str(bp.get("bp_category") or ""),
        "insertion_subject_gap_bp": "" if bp.get("bp_max_subject_gap") is None else str(bp.get("bp_max_subject_gap")),
        "is_element_best_hit": normalize(is_hit.get("is_element_best_hit")),
        "is_element_best_hit_pident": normalize(is_hit.get("is_element_best_hit_pident")),
        "is_element_best_hit_qcov": normalize(is_hit.get("is_element_best_hit_qcov")),
        "is_support_tier": is_support,
    }


def prn_status(row: dict[str, str]) -> str:
    mechanism = normalize(row.get("prn_mechanism_call"))
    if mechanism == "intact":
        return "PRN+"
    if mechanism in DISRUPTED_MECHANISMS:
        return "PRN-"
    return ""


def load_manuscript_default(path: Path) -> dict[str, dict[str, str]]:
    defaults: dict[str, dict[str, str]] = {}
    if not path.exists():
        return defaults
    for row in read_tsv(path):
        sample_id = normalize(row.get("sample_id_canonical"))
        if sample_id:
            defaults[sample_id] = row
    return defaults


def summarize_grid(
    calls: list[dict[str, str]],
    manuscript_default: dict[str, dict[str, str]],
    *,
    grid_id: str,
    hsp_min_pident: float,
    locus_qcov_threshold: float,
    is_profile: dict[str, float | str],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    mechanisms = Counter(row["prn_mechanism_call"] for row in calls)
    events = Counter(
        row["prn_event_id"]
        for row in calls
        if row["prn_mechanism_call"] in DISRUPTED_MECHANISMS and row["prn_event_id"]
    )
    top_events = events.most_common(5)
    top3_count = sum(count for _event, count in events.most_common(3))
    resolved = sum(mechanisms[mechanism] for mechanism in DISRUPTED_MECHANISMS)
    interpretable = resolved + mechanisms["intact"]
    status_changed = 0
    mechanism_changed = 0
    event_changed = 0
    comparable_status = 0
    comparable_mechanism = 0
    for row in calls:
        default = manuscript_default.get(row["sample_id"])
        if not default:
            continue
        grid_status = prn_status(row)
        default_status = normalize(default.get("repo_prn_status"))
        if grid_status and default_status:
            comparable_status += 1
            if grid_status != default_status:
                status_changed += 1
        default_mechanism = normalize(default.get("prn_mechanism_call"))
        if row["prn_mechanism_call"] and default_mechanism:
            comparable_mechanism += 1
            if row["prn_mechanism_call"] != default_mechanism:
                mechanism_changed += 1
        default_event = normalize(default.get("prn_event_id"))
        if row["prn_event_id"] and default_event:
            if row["prn_event_id"] != default_event:
                event_changed += 1

    dominant_event, dominant_event_count = ("", 0)
    if events:
        dominant_event, dominant_event_count = events.most_common(1)[0]
    summary = {
        "grid_id": grid_id,
        "hsp_min_pident": f"{hsp_min_pident:.1f}",
        "locus_qcov_threshold": f"{locus_qcov_threshold:.1f}",
        "is_support_profile": is_profile["is_support_profile"],
        "strong_is_qcov": f"{float(is_profile['strong_is_qcov']):.1f}",
        "strong_is_pident": f"{float(is_profile['strong_is_pident']):.1f}",
        "moderate_is_qcov": f"{float(is_profile['moderate_is_qcov']):.1f}",
        "moderate_is_pident": f"{float(is_profile['moderate_is_pident']):.1f}",
        "n_retained": len(calls),
        "n_interpretable": interpretable,
        "n_intact": mechanisms["intact"],
        "n_structural_disrupted": resolved,
        "n_is481": mechanisms["coding_disrupted_is481"],
        "n_rearrangement": mechanisms["coding_disrupted_inversion_or_rearrangement"],
        "n_other_disruption": mechanisms["coding_disrupted_other"],
        "n_other_is": mechanisms["coding_disrupted_other_is"],
        "n_uncertain_fragmented": mechanisms["uncertain_fragmented_assembly"],
        "n_insufficient": mechanisms["insufficient_data"],
        "n_strong_is_support": sum(1 for row in calls if row["is_support_tier"] == "strong"),
        "n_moderate_is_support": sum(1 for row in calls if row["is_support_tier"] == "moderate"),
        "dominant_event_id": dominant_event,
        "dominant_event_count": dominant_event_count,
        "dominant_event_share": "" if resolved == 0 else f"{dominant_event_count / resolved:.3f}",
        "top3_event_count": top3_count,
        "top3_event_share": "" if resolved == 0 else f"{top3_count / resolved:.3f}",
        "gap1043_count": events["prn_evt_coding_disrupted_is481__is481__gap1043"],
        "cov58_count": events["prn_evt_rearrangement__within_contig__cov58"],
        "cov91_count": events["prn_evt_rearrangement__within_contig__cov91"],
        "status_compared_to_manuscript_n": comparable_status,
        "status_changed_vs_manuscript_n": status_changed,
        "status_changed_vs_manuscript_fraction": "" if comparable_status == 0 else f"{status_changed / comparable_status:.3f}",
        "mechanism_compared_to_manuscript_n": comparable_mechanism,
        "mechanism_changed_vs_manuscript_n": mechanism_changed,
        "mechanism_changed_vs_manuscript_fraction": "" if comparable_mechanism == 0 else f"{mechanism_changed / comparable_mechanism:.3f}",
        "event_changed_vs_manuscript_n": event_changed,
    }
    event_rows = [
        {
            "grid_id": grid_id,
            "event_rank": rank,
            "prn_event_id": event_id,
            "sample_count": count,
            "sample_share_structural_disrupted": "" if resolved == 0 else f"{count / resolved:.3f}",
        }
        for rank, (event_id, count) in enumerate(top_events, start=1)
    ]
    return summary, event_rows


def add_default_deltas(rows: list[dict[str, object]]) -> None:
    defaults = [
        row
        for row in rows
        if row["hsp_min_pident"] == "90.0"
        and row["locus_qcov_threshold"] == "95.0"
        and row["is_support_profile"] == "default"
    ]
    if not defaults:
        return
    default = defaults[0]
    for row in rows:
        for field in [
            "n_interpretable",
            "n_intact",
            "n_structural_disrupted",
            "n_is481",
            "n_rearrangement",
            "n_other_disruption",
            "n_uncertain_fragmented",
            "n_insufficient",
            "gap1043_count",
            "cov58_count",
            "cov91_count",
        ]:
            row[f"delta_{field}_vs_default_grid"] = int(row[field]) - int(default[field])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=48)
    parser.add_argument("--blastn", default=os.environ.get("BLASTN", "blastn"))
    parser.add_argument("--min-cache-pident", type=float, default=80.0)
    parser.add_argument("--max-targets", type=int, default=250)
    parser.add_argument("--force-blast", action="store_true")
    parser.add_argument("--max-genomes", type=int, default=0, help="Development subset size; 0 means all retained genomes.")
    args = parser.parse_args()

    genome_paths = load_genome_paths(GENOME_PATHS)
    raw_genome_paths = load_genome_paths(RAW_GENOME_PATHS)
    records = load_qc_records(QC_MANIFEST, genome_paths, None)
    existing_samples = {record["sample_id"] for record in records}
    for record in load_raw_manuscript_records(MANUSCRIPT_ANNOTATION, raw_genome_paths):
        if record["sample_id"] not in existing_samples:
            records.append(record)
            existing_samples.add(record["sample_id"])
    if args.max_genomes:
        records = records[: args.max_genomes]
    missing = [record["assembly_accession"] for record in records if not record["fasta_path"]]
    if missing:
        raise SystemExit(f"Missing FASTA paths for {len(missing)} retained records; examples: {', '.join(missing[:5])}")
    cache_path = HSP_CACHE if args.max_genomes == 0 else GRID_DIR / f"prn_blast_hsps_min80_dev{args.max_genomes}.tsv"
    build_hsp_cache(
        records,
        cache_path,
        PRN_QUERY,
        args.blastn,
        args.min_cache_pident,
        args.jobs,
        args.max_targets,
        args.force_blast,
    )
    hsps_by_sample = load_hsps(cache_path)
    manuscript_default = load_manuscript_default(MANUSCRIPT_ANNOTATION)
    references = load_reference_fasta(IS_REFERENCE)

    summary_rows: list[dict[str, object]] = []
    event_rows: list[dict[str, object]] = []
    gap_cache: dict[tuple[str, str, int, int], str] = {}
    is_cache: dict[tuple[str, str, int, int], dict[str, str]] = {}
    combo_index = 0
    for hsp_min_pident in HSP_PIDENT_GRID:
        for locus_qcov_threshold in LOCUS_QCOV_GRID:
            for is_profile in IS_SUPPORT_PROFILES:
                combo_index += 1
                grid_id = f"grid_{combo_index:03d}"
                calls = [
                    classify_sample(
                        record,
                        hsps_by_sample.get(record["sample_id"], []),
                        hsp_min_pident=hsp_min_pident,
                        locus_qcov_threshold=locus_qcov_threshold,
                        is_profile=is_profile,
                        references=references,
                        gap_cache=gap_cache,
                        is_cache=is_cache,
                    )
                    for record in records
                ]
                summary, top_events = summarize_grid(
                    calls,
                    manuscript_default,
                    grid_id=grid_id,
                    hsp_min_pident=hsp_min_pident,
                    locus_qcov_threshold=locus_qcov_threshold,
                    is_profile=is_profile,
                )
                summary_rows.append(summary)
                event_rows.extend(top_events)
    add_default_deltas(summary_rows)

    summary_fields = [
        "grid_id",
        "hsp_min_pident",
        "locus_qcov_threshold",
        "is_support_profile",
        "strong_is_qcov",
        "strong_is_pident",
        "moderate_is_qcov",
        "moderate_is_pident",
        "n_retained",
        "n_interpretable",
        "n_intact",
        "n_structural_disrupted",
        "n_is481",
        "n_rearrangement",
        "n_other_disruption",
        "n_other_is",
        "n_uncertain_fragmented",
        "n_insufficient",
        "n_strong_is_support",
        "n_moderate_is_support",
        "dominant_event_id",
        "dominant_event_count",
        "dominant_event_share",
        "top3_event_count",
        "top3_event_share",
        "gap1043_count",
        "cov58_count",
        "cov91_count",
        "status_compared_to_manuscript_n",
        "status_changed_vs_manuscript_n",
        "status_changed_vs_manuscript_fraction",
        "mechanism_compared_to_manuscript_n",
        "mechanism_changed_vs_manuscript_n",
        "mechanism_changed_vs_manuscript_fraction",
        "event_changed_vs_manuscript_n",
        "delta_n_interpretable_vs_default_grid",
        "delta_n_intact_vs_default_grid",
        "delta_n_structural_disrupted_vs_default_grid",
        "delta_n_is481_vs_default_grid",
        "delta_n_rearrangement_vs_default_grid",
        "delta_n_other_disruption_vs_default_grid",
        "delta_n_uncertain_fragmented_vs_default_grid",
        "delta_n_insufficient_vs_default_grid",
        "delta_gap1043_count_vs_default_grid",
        "delta_cov58_count_vs_default_grid",
        "delta_cov91_count_vs_default_grid",
    ]
    event_fields = ["grid_id", "event_rank", "prn_event_id", "sample_count", "sample_share_structural_disrupted"]
    write_tsv(SUMMARY_OUT, summary_fields, summary_rows)
    write_tsv(EVENT_OUT, event_fields, event_rows)
    write_tsv(SUPPLEMENTARY_OUT, summary_fields, summary_rows)
    print(f"[threshold-grid] wrote {SUMMARY_OUT.relative_to(ROOT)}")
    print(f"[threshold-grid] wrote {EVENT_OUT.relative_to(ROOT)}")
    print(f"[threshold-grid] wrote {SUPPLEMENTARY_OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
