#!/usr/bin/env python3

import argparse
import concurrent.futures as cf
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


def norm(s: pd.Series) -> pd.Series:
    return s.astype(str).where(~s.isna(), "NA").replace({"nan": "NA", "None": "NA"})


@dataclass(frozen=True)
class Hsp:
    pident: float
    length: int
    qlen: int
    qstart: int
    qend: int
    sseqid: str
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


def parse_hsps(raw: str, min_pident: float) -> list[Hsp]:
    hsps: list[Hsp] = []
    for line in raw.splitlines():
        parts = line.rstrip("\n").split("\t")
        # outfmt: pident length qlen qstart qend sseqid slen sstart send bitscore
        if len(parts) != 10:
            continue
        try:
            pident = float(parts[0])
            if pident < min_pident:
                continue
            hsps.append(
                Hsp(
                    pident=pident,
                    length=int(float(parts[1])),
                    qlen=int(float(parts[2])),
                    qstart=int(float(parts[3])),
                    qend=int(float(parts[4])),
                    sseqid=str(parts[5]),
                    slen=int(float(parts[6])),
                    sstart=int(float(parts[7])),
                    send=int(float(parts[8])),
                    bitscore=float(parts[9]),
                )
            )
        except Exception:
            continue
    return hsps


def greedy_contributing_hsps(hsps: list[Hsp], qcov_target: float = 95.0) -> list[Hsp]:
    if not hsps:
        return []
    qlen = hsps[0].qlen
    if qlen <= 0:
        return []

    hsps_sorted = sorted(hsps, key=lambda h: (h.bitscore, h.pident, h.length), reverse=True)

    intervals: list[tuple[int, int]] = []

    def merge_intervals(intervals_in: list[tuple[int, int]]) -> list[tuple[int, int]]:
        if not intervals_in:
            return []
        intervals2 = sorted(intervals_in)
        merged: list[tuple[int, int]] = []
        cur_a, cur_b = intervals2[0]
        for a, b in intervals2[1:]:
            if a <= cur_b + 1:
                cur_b = max(cur_b, b)
            else:
                merged.append((cur_a, cur_b))
                cur_a, cur_b = a, b
        merged.append((cur_a, cur_b))
        return merged

    def cov_len(intervals_in: list[tuple[int, int]]) -> int:
        return sum(b - a + 1 for a, b in intervals_in)

    contributing: list[Hsp] = []
    for h in hsps_sorted:
        before = cov_len(intervals)
        intervals2 = merge_intervals(intervals + [(h.q_a, h.q_b)])
        after = cov_len(intervals2)
        if after > before:
            contributing.append(h)
            intervals = intervals2
        if 100.0 * (after / qlen) >= qcov_target:
            break

    return contributing


def blast_prn_detailed(query_fa: Path, subject_fa: Path, blast_threads: int, max_targets: int) -> str:
    outfmt = "6 pident length qlen qstart qend sseqid slen sstart send bitscore"
    cmd = [
        "blastn",
        "-query",
        str(query_fa),
        "-subject",
        str(subject_fa),
        "-outfmt",
        outfmt,
        "-num_threads",
        str(max(1, int(blast_threads))),
        "-max_target_seqs",
        str(max_targets),
        "-task",
        "blastn",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError((proc.stderr or "").strip() or "blastn failed")
    return proc.stdout or ""


def classify_breakpoint(contrib: list[Hsp]) -> tuple[str, dict[str, str]]:
    if not contrib:
        return "no_contrib", {}

    # Best single HSP query coverage
    best = max(contrib, key=lambda h: (h.bitscore, h.pident, h.length))
    best_single_qcov = 0.0
    if best.qlen > 0:
        best_single_qcov = 100.0 * ((best.q_b - best.q_a + 1) / best.qlen)

    contigs = sorted({h.sseqid for h in contrib})
    n_contigs = len(contigs)

    # Distance to contig ends for each contributing HSP
    end_dists = []
    for h in contrib:
        if h.slen > 0:
            end_dists.append(min(h.s_a - 1, h.slen - h.s_b))
    min_end_dist = min(end_dists) if end_dists else None

    # Compute subject gaps between query-ordered HSPs when on same contig and same strand
    ordered = sorted(contrib, key=lambda h: (h.q_a, h.q_b, -h.bitscore))
    max_subject_gap = None
    max_query_gap = None
    max_gap_pair = None  # (left_hsp, right_hsp)
    same_contig_pairs = 0
    same_strand_pairs = 0

    for a, b in zip(ordered, ordered[1:]):
        q_gap = max(0, b.q_a - a.q_b - 1)
        if a.sseqid == b.sseqid:
            same_contig_pairs += 1
            # subject gap on forward coordinate
            s_gap = max(0, b.s_a - a.s_b - 1)
            if a.same_strand and b.same_strand:
                same_strand_pairs += 1
            if max_subject_gap is None or s_gap > max_subject_gap:
                max_subject_gap = s_gap
                max_gap_pair = (a, b)
            max_query_gap = q_gap if max_query_gap is None else max(max_query_gap, q_gap)
        else:
            # across contigs we don't compute a meaningful subject gap
            max_query_gap = q_gap if max_query_gap is None else max(max_query_gap, q_gap)

    # Heuristic categories
    # - fragmented_contigs: contributing HSPs span multiple contigs
    # - contig_end_supported: breakpoints appear near contig ends (assembly fragmentation plausible)
    # - insertion_like: same contig, union coverage needs >1 HSP, and subject gap substantially exceeds query gap
    # - complex: otherwise multi-HSP but not clearly insertion-like

    if n_contigs >= 2:
        cat = "fragmented_contigs"
    else:
        # single contig
        if min_end_dist is not None and min_end_dist <= 200:
            cat = "near_contig_end"
        else:
            cat = "within_contig"

        # insertion-like signal
        if max_subject_gap is not None:
            qg = max_query_gap or 0
            sg = max_subject_gap
            # if query has tiny/zero gap but subject has a sizeable gap -> insertion
            if sg >= 50 and (qg <= 10 or sg >= 5 * max(1, qg)):
                cat = "insertion_like"

    gap_start = None
    gap_end = None
    if cat == "insertion_like" and max_gap_pair is not None:
        left, right = max_gap_pair
        # 1-based inclusive coordinates for the gap on the subject contig
        gap_start = left.s_b + 1
        gap_end = right.s_a - 1
        if gap_end < gap_start:
            gap_start, gap_end = None, None

    meta = {
        "bp_category": cat,
        "bp_contrib_hsp_n": str(len(contrib)),
        "bp_contig_n": str(n_contigs),
        "bp_best_single_qcov_pct": f"{best_single_qcov:.2f}",
        "bp_min_dist_to_contig_end": "NA" if min_end_dist is None else str(int(min_end_dist)),
        "bp_max_subject_gap": "NA" if max_subject_gap is None else str(int(max_subject_gap)),
        "bp_max_query_gap": "NA" if max_query_gap is None else str(int(max_query_gap)),
        "bp_same_contig_pairs": str(int(same_contig_pairs)),
        "bp_same_strand_pairs": str(int(same_strand_pairs)),
        "bp_contig_id": contigs[0] if contigs else "NA",
        "bp_gap_start": "NA" if gap_start is None else str(int(gap_start)),
        "bp_gap_end": "NA" if gap_end is None else str(int(gap_end)),
    }
    return cat, meta


def scan_one(
    acc: str,
    fasta_path: str,
    year: str,
    country: str,
    mlst_st: str,
    query_fa: str,
    blast_threads: int,
    max_targets: int,
    min_pident: float,
) -> dict[str, str]:
    p = Path(fasta_path)
    if not p.exists() or p.stat().st_size == 0:
        return {
            "genome_resolved_accession": acc,
            "genome_fasta_path": fasta_path,
            "year": year,
            "country": country,
            "mlst_st": mlst_st,
            "bp_category": "missing_fasta",
        }

    try:
        raw = blast_prn_detailed(Path(query_fa), p, blast_threads=blast_threads, max_targets=max_targets)
        hsps = parse_hsps(raw, min_pident=min_pident)
        contrib = greedy_contributing_hsps(hsps, qcov_target=95.0)
        _, meta = classify_breakpoint(contrib)
        return {
            "genome_resolved_accession": acc,
            "genome_fasta_path": fasta_path,
            "year": year,
            "country": country,
            "mlst_st": mlst_st,
            **meta,
        }
    except Exception as e:
        return {
            "genome_resolved_accession": acc,
            "genome_fasta_path": fasta_path,
            "year": year,
            "country": country,
            "mlst_st": mlst_st,
            "bp_category": "error",
            "error": str(e),
        }


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Step3F: stronger evidence for prn disruption: contig fragmentation vs insertion-like within-contig gaps"
    )
    ap.add_argument("--calls", required=True, help="Step3C calls TSV (bp_prn_disruption_calls.tsv)")
    ap.add_argument("--prn-query", required=True, help="prn query FASTA")
    ap.add_argument("--out", required=True, help="Output evidence TSV")
    ap.add_argument("--jobs", type=int, default=40, help="Parallelism across genomes")
    ap.add_argument("--executor", choices=["thread", "process"], default="thread")
    ap.add_argument("--blast-threads", type=int, default=1, help="blastn -num_threads per task")
    ap.add_argument("--max-targets", type=int, default=200, help="blastn max_target_seqs")
    ap.add_argument("--min-pident", type=float, default=90.0, help="Minimum pident for HSP inclusion")
    args = ap.parse_args()

    if subprocess.run(["bash", "-lc", "command -v blastn"], capture_output=True).returncode != 0:
        raise SystemExit("ERROR: blastn not found in PATH")

    query = Path(args.prn_query)
    if not query.exists() or query.stat().st_size == 0:
        raise SystemExit(f"ERROR: prn query missing/empty: {query}")

    calls = pd.read_csv(Path(args.calls), sep="\t", dtype=str)
    for c in ["genome_resolved_accession", "genome_fasta_path", "prn_call"]:
        if c not in calls.columns:
            raise SystemExit(f"ERROR: calls TSV missing required column: {c}")

    # Only analyze disrupted calls (where the interpretation matters)
    calls["prn_call"] = norm(calls["prn_call"])
    subset = calls[calls["prn_call"] == "disrupted_multi_hsp"].copy()

    for c in ["genome_resolved_accession", "genome_fasta_path", "year", "country", "mlst_st"]:
        if c not in subset.columns:
            subset[c] = "NA"
        subset[c] = norm(subset[c])

    subset = subset.drop_duplicates(subset=["genome_resolved_accession"], keep="first")

    jobs = max(1, int(args.jobs))
    jobs = min(jobs, len(subset)) if len(subset) > 0 else 1

    Executor = cf.ThreadPoolExecutor if args.executor == "thread" else cf.ProcessPoolExecutor

    rows: list[dict[str, str]] = []
    with Executor(max_workers=jobs) as ex:
        futs = [
            ex.submit(
                scan_one,
                str(r["genome_resolved_accession"]),
                str(r["genome_fasta_path"]),
                str(r.get("year", "NA")),
                str(r.get("country", "NA")),
                str(r.get("mlst_st", "NA")),
                str(query),
                int(args.blast_threads),
                int(args.max_targets),
                float(args.min_pident),
            )
            for _, r in subset.iterrows()
        ]

        for i, fut in enumerate(cf.as_completed(futs), start=1):
            rows.append(fut.result())
            if i % 100 == 0:
                print(f"[bp] finished {i}/{len(futs)}", flush=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"Wrote: {out}")


if __name__ == "__main__":
    main()
