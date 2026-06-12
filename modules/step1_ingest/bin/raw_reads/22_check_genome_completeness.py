#!/usr/bin/env python3
"""22_check_genome_completeness.py — Verify genome database integrity.

Checks that every sample in the manifest has a valid assembly file in
pertussis_data/bp_genomes_qc/assemblies/ and reports:
  - Total vs available vs missing counts
  - File integrity (non-empty, valid FASTA header)
  - Duplicate or orphan files
  - Size statistics

Usage:
    python modules/step1_ingest/bin/raw_reads/22_check_genome_completeness.py
    python modules/step1_ingest/bin/raw_reads/22_check_genome_completeness.py --strict   # exit 1 if missing > 0
"""

import csv
import os
import sys
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parents[4]
MANIFEST = ROOT / "outputs" / "workflow" / "manifest" / "manifest.tsv"
ASSEMBLY_DIR = ROOT / "pertussis_data" / "bp_genomes_qc" / "assemblies"


def is_valid_fasta(path: Path) -> tuple[bool, str]:
    """Quick validation: file exists, non-empty, starts with '>'."""
    if not path.is_file():
        return False, "not_found"
    sz = path.stat().st_size
    if sz == 0:
        return False, "empty"
    try:
        with open(path, "r") as f:
            first = f.readline()
            if not first.startswith(">"):
                return False, "no_fasta_header"
    except Exception:
        return False, "read_error"
    return True, "ok"


def main():
    strict = "--strict" in sys.argv

    # Load manifest accessions
    manifest_accs = set()
    with open(MANIFEST, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            manifest_accs.add(row["assembly_accession"])
    print(f"Manifest accessions: {len(manifest_accs)}")

    # Check each accession
    status = Counter()
    problems = []
    sizes = []

    for acc in sorted(manifest_accs):
        fasta = ASSEMBLY_DIR / f"{acc}.fasta"
        valid, reason = is_valid_fasta(fasta)
        if valid:
            status["ok"] += 1
            sizes.append(fasta.stat().st_size)
        else:
            status[reason] += 1
            problems.append((acc, reason))

    # Check for orphan files (in directory but not in manifest)
    dir_files = set()
    if ASSEMBLY_DIR.is_dir():
        for f in ASSEMBLY_DIR.iterdir():
            if f.suffix == ".fasta" and f.is_file():
                acc = f.stem
                dir_files.add(acc)

    orphans = dir_files - manifest_accs
    in_manifest_only = manifest_accs - dir_files

    # Report
    print(f"\n{'='*55}")
    print(f"  Genome Database Integrity Report")
    print(f"{'='*55}")
    print(f"\n  Directory: {ASSEMBLY_DIR}/")
    print(f"  Manifest accessions: {len(manifest_accs)}")
    print(f"  Files in directory:  {len(dir_files)}")
    print(f"  Orphan files:        {len(orphans)}")
    print(f"\n  Status:")
    for k, v in sorted(status.items()):
        marker = "  " if k == "ok" else "⚠ "
        print(f"    {marker}{k:20s}: {v}")

    if sizes:
        import statistics
        print(f"\n  Size statistics (valid files):")
        print(f"    min:    {min(sizes)/1e6:.1f} MB")
        print(f"    median: {statistics.median(sizes)/1e6:.1f} MB")
        print(f"    max:    {max(sizes)/1e6:.1f} MB")
        print(f"    total:  {sum(sizes)/1e9:.1f} GB")

    pct = status["ok"] / len(manifest_accs) * 100 if manifest_accs else 0
    print(f"\n  Completeness: {status['ok']}/{len(manifest_accs)} ({pct:.1f}%)")

    if problems:
        print(f"\n  Missing/invalid ({len(problems)}):")
        for acc, reason in problems[:20]:
            print(f"    {acc}: {reason}")
        if len(problems) > 20:
            print(f"    ... and {len(problems)-20} more")

        # Write missing list for retry
        missing_path = ASSEMBLY_DIR / "_metadata" / "missing.txt"
        missing_path.parent.mkdir(parents=True, exist_ok=True)
        with open(missing_path, "w") as f:
            for acc, _ in problems:
                f.write(acc + "\n")
        print(f"\n  Missing list: {missing_path}")

    if orphans:
        print(f"\n  Orphan files (not in manifest): {len(orphans)}")
        for o in sorted(orphans)[:10]:
            print(f"    {o}.fasta")

    print()

    if strict and status.get("not_found", 0) + status.get("empty", 0) > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
