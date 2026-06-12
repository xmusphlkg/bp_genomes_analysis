#!/usr/bin/env python3
"""T08 — Assembly QC: compute basic genome statistics for all available assemblies.

Reads the assembly_map.tsv and computes per-genome stats without external tools:
  - total_length, n_contigs, gc_pct, longest_contig, n50, l50

For deeper QC (QUAST), use the Snakemake reads_qc/assembly_qc rules.

Outputs:
    workflow/assembly_qc/assembly_qc_stats.tsv
    workflow/assembly_qc/assembly_qc_pass_fail.tsv  (with pass/fail flags)

Usage:
    python modules/step1_ingest/bin/raw_reads/26_run_assembly_qc.py [--threads 4]
"""

import csv
import sys
import math
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import argparse

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from workflow.lib.project_paths import project_workflow_root

MANIFEST = ROOT / "state" / "manifest" / "manifest.tsv"
ASSEMBLY_DIR = ROOT / "pertussis_data" / "bp_genomes_qc" / "assemblies"
OUTPUT_DIR = project_workflow_root() / "assembly_qc"

# QC thresholds (from config.yaml)
QC = {
    "min_assembly_length": 3_800_000,
    "max_assembly_length": 4_600_000,
    "max_contigs": 400,
    "min_n50": 10_000,
}


def parse_fasta(fasta_path: str) -> list[tuple[str, str]]:
    """Parse FASTA file, return list of (header, sequence) tuples."""
    sequences = []
    header = ""
    seq_parts = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header:
                    sequences.append((header, "".join(seq_parts)))
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
    if header:
        sequences.append((header, "".join(seq_parts)))
    return sequences


def compute_stats(fasta_path: str) -> dict:
    """Compute assembly statistics for a single FASTA file."""
    try:
        contigs = parse_fasta(fasta_path)
    except Exception as e:
        return {"error": str(e)}

    if not contigs:
        return {"error": "empty_fasta"}

    lengths = sorted([len(seq) for _, seq in contigs], reverse=True)
    total_length = sum(lengths)
    n_contigs = len(lengths)
    longest = lengths[0]

    # GC content
    total_gc = sum(seq.count("G") + seq.count("C") for _, seq in contigs)
    total_bases = sum(seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C") for _, seq in contigs)
    gc_pct = (total_gc / total_bases * 100) if total_bases > 0 else 0.0

    # N50 / L50
    cumsum = 0
    n50 = 0
    l50 = 0
    half = total_length / 2
    for i, length in enumerate(lengths):
        cumsum += length
        if cumsum >= half:
            n50 = length
            l50 = i + 1
            break

    return {
        "total_length": total_length,
        "n_contigs": n_contigs,
        "gc_pct": round(gc_pct, 2),
        "longest_contig": longest,
        "n50": n50,
        "l50": l50,
    }


def qc_verdict(stats: dict) -> tuple[str, list[str]]:
    """Apply QC thresholds, return (PASS|FAIL, [reasons])."""
    if "error" in stats:
        return "FAIL", [f"error:{stats['error']}"]
    reasons = []
    if stats["total_length"] < QC["min_assembly_length"]:
        reasons.append(f"short:{stats['total_length']}")
    if stats["total_length"] > QC["max_assembly_length"]:
        reasons.append(f"long:{stats['total_length']}")
    if stats["n_contigs"] > QC["max_contigs"]:
        reasons.append(f"fragmented:{stats['n_contigs']}")
    if stats["n50"] < QC["min_n50"]:
        reasons.append(f"low_n50:{stats['n50']}")
    return ("FAIL" if reasons else "PASS"), reasons


def process_one(row: dict) -> dict:
    """Process one manifest entry."""
    sample = row["sample_id_canonical"]
    acc = row["assembly_accession"]
    fasta = ASSEMBLY_DIR / f"{acc}.fasta"

    if not fasta.is_file() or fasta.stat().st_size == 0:
        return {
            "sample_id_canonical": sample,
            "assembly_accession": acc,
            "qc_status": "NO_FASTA",
            "qc_reasons": "missing_fasta",
            "total_length": "",
            "n_contigs": "",
            "gc_pct": "",
            "longest_contig": "",
            "n50": "",
            "l50": "",
        }

    stats = compute_stats(str(fasta))
    verdict, reasons = qc_verdict(stats)

    return {
        "sample_id_canonical": sample,
        "assembly_accession": acc,
        "qc_status": verdict,
        "qc_reasons": ";".join(reasons) if reasons else "",
        "total_length": stats.get("total_length", ""),
        "n_contigs": stats.get("n_contigs", ""),
        "gc_pct": stats.get("gc_pct", ""),
        "longest_contig": stats.get("longest_contig", ""),
        "n50": stats.get("n50", ""),
        "l50": stats.get("l50", ""),
    }


def main():
    parser = argparse.ArgumentParser(description="Assembly QC")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    if not MANIFEST.exists():
        print(f"ERROR: Manifest not found: {MANIFEST}", file=sys.stderr)
        sys.exit(1)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load manifest
    with open(MANIFEST, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    print(f"Processing {len(rows)} samples ({args.threads} threads)...")

    # Process in parallel
    with ProcessPoolExecutor(max_workers=args.threads) as pool:
        results = list(pool.map(process_one, rows))

    # Write results
    out_cols = [
        "sample_id_canonical", "assembly_accession", "qc_status", "qc_reasons",
        "total_length", "n_contigs", "gc_pct", "longest_contig", "n50", "l50",
    ]
    out_path = OUTPUT_DIR / "assembly_qc_stats.tsv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    # Summary
    from collections import Counter
    status_counts = Counter(r["qc_status"] for r in results)
    print(f"\nQC Results:")
    for status, count in sorted(status_counts.items()):
        print(f"  {status}: {count}")
    print(f"\nOutput: {out_path}")


if __name__ == "__main__":
    main()
