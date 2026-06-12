#!/usr/bin/env python3

import argparse
import concurrent.futures as cf
import subprocess
from dataclasses import dataclass
from pathlib import Path


FASTA_SUFFIXES = {".fa", ".fasta", ".fna"}


def read_genome_fasta_paths(genome_paths_tsv: Path) -> list[str]:
    paths: list[str] = []
    seen: set[str] = set()
    with genome_paths_tsv.open("r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx_status = header.index("status")
        idx_path = header.index("fasta_path")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx_status, idx_path):
                continue
            if parts[idx_status] != "ok":
                continue
            p = parts[idx_path].strip()
            if not p or p in seen:
                continue
            seen.add(p)
            paths.append(p)
    return paths


@dataclass(frozen=True)
class Hit:
    genome_fasta_path: str
    qseqid: str
    sseqid: str
    pident: float
    length: int
    qlen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    bitscore: float
    qseq: str
    sseq: str

    @property
    def qcov_pct(self) -> float:
        return 100.0 * (self.length / self.qlen) if self.qlen else 0.0


def blast_23s(query_fa: Path, subject_fa: Path, max_targets: int, blast_threads: int) -> str:
    outfmt = "6 qseqid sseqid pident length qlen qstart qend sstart send bitscore qseq sseq"
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
        raise RuntimeError(proc.stderr.strip() or "blastn failed")
    return proc.stdout or ""


def parse_hits(genome_fa: str, raw: str) -> list[Hit]:
    hits: list[Hit] = []
    for line in raw.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 12:
            continue
        try:
            hits.append(
                Hit(
                    genome_fasta_path=genome_fa,
                    qseqid=parts[0],
                    sseqid=parts[1],
                    pident=float(parts[2]),
                    length=int(float(parts[3])),
                    qlen=int(float(parts[4])),
                    qstart=int(float(parts[5])),
                    qend=int(float(parts[6])),
                    sstart=int(float(parts[7])),
                    send=int(float(parts[8])),
                    bitscore=float(parts[9]),
                    qseq=parts[10],
                    sseq=parts[11],
                )
            )
        except Exception:
            continue
    hits.sort(key=lambda h: (h.bitscore, h.pident, h.length), reverse=True)
    return hits


def base_at_query_pos(hit: Hit, query_pos_1based: int) -> str | None:
    # Only valid if the hit covers the query position.
    q_lo = min(hit.qstart, hit.qend)
    q_hi = max(hit.qstart, hit.qend)
    if not (q_lo <= query_pos_1based <= q_hi):
        return None

    qpos = hit.qstart - 1
    for qc, sc in zip(hit.qseq, hit.sseq):
        if qc != "-":
            qpos += 1
        if qpos == query_pos_1based:
            if sc == "-":
                return None
            return sc.upper()
    return None


def scan_one(
    genome_fa: str,
    query_fa: str,
    query_pos: int,
    max_targets: int,
    min_pident: float,
    min_qcov: float,
    blast_threads: int,
) -> tuple[list[dict[str, str]], dict[str, str]]:
    subject = Path(genome_fa)
    per_hit: list[dict[str, str]] = []
    per_genome: dict[str, str] = {"genome_fasta_path": genome_fa, "call": "unknown", "n_hits": "0"}
    if not subject.exists() or subject.stat().st_size == 0:
        per_genome["call"] = "missing_fasta"
        return per_hit, per_genome

    try:
        raw = blast_23s(Path(query_fa), subject, max_targets=max_targets, blast_threads=blast_threads)
    except Exception as e:
        per_genome["call"] = "error"
        per_genome["error"] = str(e)
        return per_hit, per_genome

    hits = parse_hits(genome_fa, raw)
    per_genome["n_hits"] = str(len(hits))
    if not hits:
        per_genome["call"] = "no_hit"
        return per_hit, per_genome

    bases: list[str] = []
    for idx, h in enumerate(hits, start=1):
        qcov = h.qcov_pct
        passed = (h.pident >= min_pident) and (qcov >= min_qcov)
        b = base_at_query_pos(h, query_pos)
        if passed and b in {"A", "C", "G", "T"}:
            bases.append(b)
        per_hit.append(
            {
                "genome_fasta_path": genome_fa,
                "hit_index": str(idx),
                "qseqid": h.qseqid,
                "sseqid": h.sseqid,
                "pident": f"{h.pident:.2f}",
                "qcov_pct": f"{qcov:.2f}",
                "qstart": str(h.qstart),
                "qend": str(h.qend),
                "sstart": str(h.sstart),
                "send": str(h.send),
                "bitscore": f"{h.bitscore:.2f}",
                "pass_threshold": "1" if passed else "0",
                "base_at_pos": b or "",
            }
        )

    # Genome-level call: any G among good hits -> A2047G
    if not bases:
        per_genome["call"] = "no_call"
        return per_hit, per_genome

    uniq = sorted(set(bases))
    if "G" in uniq:
        per_genome["call"] = "A2047G" if uniq == ["G"] else "mixed_includes_A2047G"
    elif uniq == ["A"]:
        per_genome["call"] = "WT_A2047"
    else:
        per_genome["call"] = "other_base_" + "".join(uniq)
    per_genome["bases"] = ",".join(bases)
    return per_hit, per_genome


def main() -> None:
    ap = argparse.ArgumentParser(description="Call 23S rRNA A2047G from assemblies using BLAST alignment mapping")
    ap.add_argument("--genome-paths", required=True, help="TSV from step2_02_index_genomes.py")
    ap.add_argument("--query-23s", required=True, help="23S rRNA reference FASTA (must be >= query_pos)")
    ap.add_argument("--query-pos", type=int, default=2047, help="Position in the 23S reference (1-based, default: 2047)")
    ap.add_argument("--out-hits", required=True, help="Output TSV: per-hit details")
    ap.add_argument("--out-summary", required=True, help="Output TSV: per-genome call")
    ap.add_argument("--jobs", type=int, default=8, help="Parallel genomes (default: 8)")
    ap.add_argument(
        "--executor",
        choices=["thread", "process"],
        default="thread",
        help="Concurrency backend (default: thread). Thread is usually faster for spawning blastn subprocesses.",
    )
    ap.add_argument("--max-targets", type=int, default=10, help="blastn max_target_seqs (default: 10)")
    ap.add_argument("--blast-threads", type=int, default=1, help="blastn -num_threads per task (default: 1)")
    ap.add_argument("--min-pident", type=float, default=95.0, help="Minimum percent identity (default: 95)")
    ap.add_argument("--min-qcov", type=float, default=80.0, help="Minimum query coverage (default: 80)")
    args = ap.parse_args()

    genome_paths = Path(args.genome_paths)
    if not genome_paths.exists() or genome_paths.stat().st_size == 0:
        raise SystemExit(f"ERROR: genome paths TSV missing or empty: {genome_paths}")

    query = Path(args.query_23s)
    if not query.exists() or query.stat().st_size == 0:
        raise SystemExit(f"ERROR: query FASTA missing or empty: {query}")

    if subprocess.run(["bash", "-lc", "command -v blastn"], capture_output=True).returncode != 0:
        raise SystemExit("ERROR: blastn not found. Install ncbi-blast+ (conda: conda install -c bioconda ncbi-blast)")

    genomes = read_genome_fasta_paths(genome_paths)
    if not genomes:
        raise SystemExit("ERROR: no genomes with status=ok found in genome paths TSV")

    jobs = max(1, int(args.jobs))
    jobs = min(jobs, len(genomes))

    out_hits = Path(args.out_hits)
    out_hits.parent.mkdir(parents=True, exist_ok=True)
    out_sum = Path(args.out_summary)
    out_sum.parent.mkdir(parents=True, exist_ok=True)

    hit_header = [
        "genome_fasta_path",
        "hit_index",
        "qseqid",
        "sseqid",
        "pident",
        "qcov_pct",
        "qstart",
        "qend",
        "sstart",
        "send",
        "bitscore",
        "pass_threshold",
        "base_at_pos",
    ]
    sum_header = ["genome_fasta_path", "call", "bases", "n_hits", "error"]

    with out_hits.open("w", encoding="utf-8") as fh_hits, out_sum.open("w", encoding="utf-8") as fh_sum:
        fh_hits.write("\t".join(hit_header) + "\n")
        fh_sum.write("\t".join(sum_header) + "\n")

        Executor = cf.ThreadPoolExecutor if args.executor == "thread" else cf.ProcessPoolExecutor
        with Executor(max_workers=jobs) as ex:
            futs = [
                ex.submit(
                    scan_one,
                    g,
                    str(query),
                    int(args.query_pos),
                    int(args.max_targets),
                    float(args.min_pident),
                    float(args.min_qcov),
                    int(args.blast_threads),
                )
                for g in genomes
            ]
            for i, fut in enumerate(cf.as_completed(futs), start=1):
                per_hit, per_genome = fut.result()
                for r in per_hit:
                    fh_hits.write("\t".join(r.get(c, "") for c in hit_header) + "\n")
                fh_sum.write("\t".join(per_genome.get(c, "") for c in sum_header) + "\n")
                if i % 100 == 0:
                    print(f"[23S] finished {i}/{len(genomes)} genomes", flush=True)

    print(f"Wrote: {out_hits}")
    print(f"Wrote: {out_sum}")
    print(f"Genomes (unique FASTA): {len(genomes)}")


if __name__ == "__main__":
    main()
