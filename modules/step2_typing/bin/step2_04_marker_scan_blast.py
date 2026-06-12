#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path


def list_query_fastas(ref_dir: Path) -> list[Path]:
    if not ref_dir.exists() or not ref_dir.is_dir():
        return []
    out: list[Path] = []
    for p in sorted(ref_dir.iterdir()):
        if p.is_file() and p.suffix.lower() in {".fa", ".fasta", ".fna"}:
            if p.stat().st_size > 0:
                out.append(p)
    return out


def read_ok_genomes(tsv_path: Path) -> list[tuple[str, Path]]:
    genomes: list[tuple[str, Path]] = []
    with tsv_path.open("r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx_in = header.index("input_accession") if "input_accession" in header else None
        idx_res = header.index("resolved_accession") if "resolved_accession" in header else None
        idx_status = header.index("status")
        idx_path = header.index("fasta_path")

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx_status, idx_path):
                continue
            if parts[idx_status] != "ok":
                continue
            fasta = parts[idx_path].strip()
            if not fasta:
                continue
            name = None
            if idx_res is not None and idx_res < len(parts) and parts[idx_res].strip():
                name = parts[idx_res].strip()
            elif idx_in is not None and idx_in < len(parts) and parts[idx_in].strip():
                name = parts[idx_in].strip()
            else:
                name = Path(fasta).name
            genomes.append((name, Path(fasta)))

    return genomes


def ensure_blast_db(fasta: Path, db_prefix: Path) -> None:
    # Create db only if all expected files exist.
    expected = [db_prefix.with_suffix(".nhr"), db_prefix.with_suffix(".nin"), db_prefix.with_suffix(".nsq")]
    if all(p.exists() for p in expected):
        return

    cmd = [
        "makeblastdb",
        "-dbtype",
        "nucl",
        "-in",
        str(fasta),
        "-out",
        str(db_prefix),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"ERROR: makeblastdb failed for {fasta}\nSTDERR:\n{proc.stderr}")


def blast_one(query_fa: Path, db_prefix: Path, max_targets: int) -> str:
    outfmt = "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"
    cmd = [
        "blastn",
        "-query",
        str(query_fa),
        "-db",
        str(db_prefix),
        "-outfmt",
        outfmt,
        "-max_target_seqs",
        str(max_targets),
        "-task",
        "blastn",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"ERROR: blastn failed for query {query_fa}\nSTDERR:\n{proc.stderr}")
    return proc.stdout


def main() -> None:
    ap = argparse.ArgumentParser(description="Scan assemblies for marker sequences using BLAST")
    ap.add_argument("--genome-paths", required=True, help="TSV from step2_02_index_genomes.py")
    ap.add_argument("--references-dir", required=True, help="Directory containing query FASTAs")
    ap.add_argument("--out", required=True, help="Output TSV of BLAST hits")

    ap.add_argument("--min-pident", type=float, default=90.0, help="Minimum percent identity")
    ap.add_argument("--min-qcov", type=float, default=80.0, help="Minimum query coverage percent")
    ap.add_argument("--max-targets", type=int, default=5, help="blastn max_target_seqs")
    ap.add_argument("--db-dir", default="outputs/blastdb", help="Where to store per-genome BLAST databases")
    args = ap.parse_args()

    genome_paths = Path(args.genome_paths)
    if not genome_paths.exists() or genome_paths.stat().st_size == 0:
        raise SystemExit(f"ERROR: genome paths TSV missing or empty: {genome_paths}")

    ref_dir = Path(args.references_dir)
    queries = list_query_fastas(ref_dir)
    if not queries:
        raise SystemExit(f"ERROR: no query FASTA files found in: {ref_dir}")

    genomes = read_ok_genomes(genome_paths)
    if not genomes:
        raise SystemExit("ERROR: no genomes with status=ok found in genome paths TSV")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    db_dir = Path(args.db_dir)
    db_dir.mkdir(parents=True, exist_ok=True)

    header = [
        "genome",
        "query_file",
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "qlen",
        "qcov_pct",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "pass_threshold",
    ]

    with out_path.open("w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")

        for genome_name, fasta in genomes:
            if not fasta.exists() or fasta.stat().st_size == 0:
                continue

            db_prefix = db_dir / genome_name
            ensure_blast_db(fasta, db_prefix)

            for q in queries:
                raw = blast_one(q, db_prefix, args.max_targets)
                if not raw.strip():
                    continue

                for line in raw.splitlines():
                    parts = line.split("\t")
                    if len(parts) != 12:
                        continue

                    qseqid, sseqid = parts[0], parts[1]
                    pident = float(parts[2])
                    length = int(float(parts[3]))
                    qlen = int(float(parts[4]))
                    qstart, qend = int(float(parts[6])), int(float(parts[7]))
                    sstart, send = int(float(parts[8])), int(float(parts[9]))
                    evalue, bitscore = parts[10], parts[11]

                    qcov = 100.0 * (length / qlen) if qlen else 0.0
                    passed = (pident >= args.min_pident) and (qcov >= args.min_qcov)

                    out.write(
                        "\t".join(
                            [
                                genome_name,
                                q.name,
                                qseqid,
                                sseqid,
                                f"{pident:.2f}",
                                str(length),
                                str(qlen),
                                f"{qcov:.2f}",
                                str(qstart),
                                str(qend),
                                str(sstart),
                                str(send),
                                evalue,
                                bitscore,
                                "1" if passed else "0",
                            ]
                        )
                        + "\n"
                    )

    print(f"Wrote: {out_path}")
    print(f"Genomes scanned: {len(genomes)}")
    print(f"Query FASTAs: {len(queries)}")


if __name__ == "__main__":
    main()
