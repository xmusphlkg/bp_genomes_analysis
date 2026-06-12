#!/usr/bin/env python3

import argparse
import concurrent.futures as cf
import hashlib
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


FASTA_SUFFIXES = {".fa", ".fasta", ".fna"}


def md5_hex(s: str) -> str:
    return hashlib.md5(s.encode("utf-8")).hexdigest()


def read_genome_fasta_paths(genome_paths_tsv: Path) -> list[str]:
    paths: list[str] = []
    seen: set[str] = set()
    with genome_paths_tsv.open("r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            idx_status = header.index("status")
            idx_path = header.index("fasta_path")
        except ValueError:
            raise SystemExit("ERROR: genome_paths TSV missing required columns: status, fasta_path")

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


def read_done_genomes(out_path: Path, markers: set[str]) -> set[str]:
    """Return genomes that already have rows for all markers.

    We treat status=error as NOT done so it can be retried.
    """
    if not out_path.exists() or out_path.stat().st_size == 0:
        return set()
    done: dict[str, set[str]] = {}
    with out_path.open("r", encoding="utf-8", errors="ignore") as f:
        header = f.readline().rstrip("\n").split("\t")
        if not header or "genome_fasta_path" not in header or "marker" not in header:
            return set()
        idx_g = header.index("genome_fasta_path")
        idx_m = header.index("marker")
        idx_s = header.index("status") if "status" in header else None
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx_g, idx_m):
                continue
            g = parts[idx_g].strip()
            m = parts[idx_m].strip()
            if not g or not m:
                continue
            if idx_s is not None and idx_s < len(parts) and parts[idx_s].strip() == "error":
                continue
            done.setdefault(g, set()).add(m)
    return {g for g, ms in done.items() if markers.issubset(ms)}


def list_marker_fastas(ref_dir: Path) -> list[Path]:
    if not ref_dir.exists() or not ref_dir.is_dir():
        return []
    out: list[Path] = []
    for p in sorted(ref_dir.iterdir()):
        if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES and p.stat().st_size > 0:
            out.append(p)
    return out


def read_first_fasta_seq(path: Path) -> str:
    """Read the first FASTA record's sequence (uppercase, no whitespace)."""
    seq: list[str] = []
    started = False
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if started:
                    break
                started = True
                continue
            if started:
                seq.append(line)
    return "".join(seq).replace(" ", "").replace("\t", "").upper()


def write_fasta_records(path: Path, records: Iterable[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            seq = seq.strip().replace(" ", "").replace("\t", "").replace("\r", "").replace("\n", "").upper()
            for i in range(0, len(seq), 80):
                out.write(seq[i : i + 80] + "\n")


@dataclass(frozen=True)
class BlastHit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    qlen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: str
    bitscore: float
    qseq: str
    sseq: str

    @property
    def qcov_pct(self) -> float:
        return 100.0 * (self.length / self.qlen) if self.qlen else 0.0


def parse_hits(raw: str) -> list[BlastHit]:
    hits: list[BlastHit] = []
    for line in raw.splitlines():
        parts = line.rstrip("\n").split("\t")
        # outfmt: qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qseq sseq
        if len(parts) != 13:
            continue
        try:
            hits.append(
                BlastHit(
                    qseqid=parts[0],
                    sseqid=parts[1],
                    pident=float(parts[2]),
                    length=int(float(parts[3])),
                    qlen=int(float(parts[4])),
                    qstart=int(float(parts[5])),
                    qend=int(float(parts[6])),
                    sstart=int(float(parts[7])),
                    send=int(float(parts[8])),
                    evalue=parts[9],
                    bitscore=float(parts[10]),
                    qseq=parts[11],
                    sseq=parts[12],
                )
            )
        except Exception:
            continue
    return hits


def best_hit(hits: list[BlastHit]) -> BlastHit | None:
    if not hits:
        return None
    return max(hits, key=lambda h: (h.bitscore, h.pident, h.length))


def blast_marker_against_subject(query_fa: Path, subject_fa: Path, task: str, max_targets: int) -> str:
    outfmt = "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qseq sseq"
    cmd = [
        "blastn",
        "-query",
        str(query_fa),
        "-subject",
        str(subject_fa),
        "-outfmt",
        outfmt,
        "-num_threads",
        "1",
        "-max_target_seqs",
        str(max_targets),
        "-task",
        task,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"blastn failed: query={query_fa} subject={subject_fa} stderr={proc.stderr.strip()}")
    return proc.stdout or ""


def ungap(seq: str) -> str:
    return seq.replace("-", "").replace(" ", "").upper()


def marker_name_from_path(p: Path) -> str:
    # e.g. prn_maker.fasta -> prn
    name = p.name
    for suf in FASTA_SUFFIXES:
        if name.lower().endswith(suf):
            return name[: -len(suf)]
    return p.stem


def scan_one_genome(
    genome_fa: str,
    marker_names: list[str],
    combined_query_fa: str,
    task: str,
    max_targets: int,
    min_pident: float,
    min_qcov: float,
    emit_fasta_dir: str | None,
    blast_threads: int,
) -> list[dict[str, str]]:
    subject = Path(genome_fa)
    rows: list[dict[str, str]] = []
    if not subject.exists() or subject.stat().st_size == 0:
        return rows

    if emit_fasta_dir:
        Path(emit_fasta_dir).mkdir(parents=True, exist_ok=True)
    try:
        outfmt = "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qseq sseq"
        cmd = [
            "blastn",
            "-query",
            str(combined_query_fa),
            "-subject",
            str(subject),
            "-outfmt",
            outfmt,
            "-num_threads",
            str(max(1, int(blast_threads))),
            "-max_target_seqs",
            str(max_targets),
            "-task",
            task,
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                f"blastn failed: query={combined_query_fa} subject={subject} stderr={(proc.stderr or '').strip()}"
            )
        raw = proc.stdout or ""
    except Exception as e:
        for marker in marker_names:
            rows.append(
                {
                    "genome_fasta_path": genome_fa,
                    "marker": marker,
                    "status": "error",
                    "error": str(e),
                }
            )
        return rows

    hits = parse_hits(raw)
    hits_by_marker: dict[str, list[BlastHit]] = {}
    for h in hits:
        hits_by_marker.setdefault(h.qseqid, []).append(h)

    for marker in marker_names:
        hit = best_hit(hits_by_marker.get(marker, []))
        if hit is None:
            rows.append(
                {
                    "genome_fasta_path": genome_fa,
                    "marker": marker,
                    "status": "no_hit",
                }
            )
            continue

        qcov = hit.qcov_pct
        passed = (hit.pident >= min_pident) and (qcov >= min_qcov)
        seq = ungap(hit.sseq)
        allele_hash = md5_hex(seq) if seq else ""

        if emit_fasta_dir and seq and passed:
            out_path = Path(emit_fasta_dir) / f"{marker}.fasta"
            genome_label = subject.name
            header = f">{genome_label}|{allele_hash}|pident={hit.pident:.2f}|qcov={qcov:.2f}"
            lines = [header]
            for i in range(0, len(seq), 80):
                lines.append(seq[i : i + 80])
            with out_path.open("a", encoding="utf-8") as fh:
                fh.write("\n".join(lines) + "\n")

        rows.append(
            {
                "genome_fasta_path": genome_fa,
                "marker": marker,
                "status": "ok" if passed else "below_threshold",
                "pident": f"{hit.pident:.2f}",
                "qcov_pct": f"{qcov:.2f}",
                "length": str(hit.length),
                "qlen": str(hit.qlen),
                "qstart": str(hit.qstart),
                "qend": str(hit.qend),
                "sstart": str(hit.sstart),
                "send": str(hit.send),
                "bitscore": f"{hit.bitscore:.2f}",
                "allele_len": str(len(seq)),
                "allele_hash": allele_hash,
            }
        )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Extract marker allele sequences by BLASTing marker FASTAs against each assembly (subject mode)."
    )
    ap.add_argument("--genome-paths", required=True, help="TSV from step2_02_index_genomes.py")
    ap.add_argument("--markers-dir", required=True, help="Directory with marker query FASTAs (e.g. references/markers)")
    ap.add_argument("--out", required=True, help="Output TSV (long format: genome_fasta_path x marker)")
    ap.add_argument("--jobs", type=int, default=8, help="Parallel genomes (default: 8)")
    ap.add_argument("--resume", action="store_true", help="Resume from existing --out by skipping genomes that already have all markers")
    ap.add_argument(
        "--executor",
        choices=["thread", "process"],
        default="thread",
        help="Concurrency backend (default: thread). Thread is usually faster for spawning blastn subprocesses.",
    )
    ap.add_argument("--blast-task", default="blastn", help="blastn task (default: blastn)")
    ap.add_argument("--blast-threads", type=int, default=1, help="blastn -num_threads per task (default: 1)")
    ap.add_argument("--max-targets", type=int, default=5, help="blastn max_target_seqs (default: 5)")
    ap.add_argument("--min-pident", type=float, default=90.0, help="Minimum percent identity")
    ap.add_argument("--min-qcov", type=float, default=80.0, help="Minimum query coverage percent")
    ap.add_argument(
        "--emit-fasta-dir",
        default=None,
        help="Optional output directory to write extracted allele sequences as FASTA (one file per marker; only threshold-passing hits)",
    )
    args = ap.parse_args()

    genome_paths = Path(args.genome_paths)
    if not genome_paths.exists() or genome_paths.stat().st_size == 0:
        raise SystemExit(f"ERROR: genome paths TSV missing or empty: {genome_paths}")

    markers_dir = Path(args.markers_dir)
    marker_fastas = list_marker_fastas(markers_dir)
    if not marker_fastas:
        raise SystemExit(
            f"ERROR: no marker FASTA files found in {markers_dir}. Put e.g. prn_maker.fasta, ptxP_promoter.fasta, fim2.fasta, fim3.fasta"
        )

    genomes = read_genome_fasta_paths(genome_paths)
    if not genomes:
        raise SystemExit("ERROR: no genomes with status=ok found in genome paths TSV")

    jobs = max(1, int(args.jobs))
    jobs = min(jobs, len(genomes))

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Ensure blastn exists early
    if subprocess.run(["bash", "-lc", "command -v blastn"], capture_output=True).returncode != 0:
        raise SystemExit("ERROR: blastn not found. Install ncbi-blast+ (conda: conda install -c bioconda ncbi-blast)")

    marker_fa_strs = [str(p) for p in marker_fastas]
    marker_names_set = {marker_name_from_path(p) for p in marker_fastas}
    marker_names = [marker_name_from_path(p) for p in marker_fastas]

    # Build a combined multi-FASTA query once to reduce blastn invocations
    combined_query = out_path.parent / "_markers_combined_query.fasta"
    records: list[tuple[str, str]] = []
    for p in marker_fastas:
        name = marker_name_from_path(p)
        seq = read_first_fasta_seq(p)
        if not seq:
            raise SystemExit(f"ERROR: marker FASTA has no sequence: {p}")
        records.append((name, seq))
    write_fasta_records(combined_query, records)

    done_genomes = read_done_genomes(out_path, marker_names_set) if args.resume else set()
    genomes_todo = [g for g in genomes if g not in done_genomes]
    if args.resume:
        print(f"[Markers] completed genomes: {len(done_genomes)}", flush=True)
        print(f"[Markers] remaining genomes: {len(genomes_todo)}", flush=True)

    if not genomes_todo:
        print("[Markers] nothing to do", flush=True)
        return

    header = [
        "genome_fasta_path",
        "marker",
        "status",
        "pident",
        "qcov_pct",
        "length",
        "qlen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "bitscore",
        "allele_len",
        "allele_hash",
        "error",
    ]

    wrote = 0
    mode = "a" if args.resume else "w"
    with out_path.open(mode, encoding="utf-8") as out:
        if mode == "w" or out_path.stat().st_size == 0:
            out.write("\t".join(header) + "\n")

        Executor = cf.ThreadPoolExecutor if args.executor == "thread" else cf.ProcessPoolExecutor
        with Executor(max_workers=jobs) as ex:
            futs = [
                ex.submit(
                    scan_one_genome,
                    g,
                    marker_names,
                    str(combined_query),
                    args.blast_task,
                    int(args.max_targets),
                    float(args.min_pident),
                    float(args.min_qcov),
                    args.emit_fasta_dir,
                    int(args.blast_threads),
                )
                for g in genomes_todo
            ]
            for i, fut in enumerate(cf.as_completed(futs), start=1):
                rows = fut.result()
                for r in rows:
                    out.write("\t".join(r.get(c, "") for c in header) + "\n")
                    wrote += 1
                out.flush()
                if i % 50 == 0:
                    print(f"[Markers] finished {i}/{len(genomes_todo)} genomes", flush=True)

    print(f"Wrote: {out_path}")
    print(f"Genomes (unique FASTA): {len(genomes)}")
    print(f"Markers: {len(marker_fastas)}")
    print(f"Rows written: {wrote}")


if __name__ == "__main__":
    # Make multiprocessing more stable on some Linux setups
    os.environ.setdefault("PYTHONUNBUFFERED", "1")
    main()
