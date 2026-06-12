#!/usr/bin/env python3
"""Consolidate all genome assemblies into pertussis_data/bp_genomes_qc/assemblies/.

Merges all data sources into a single flat directory:
    pertussis_data/bp_genomes_qc/assemblies/<GCA_accession>.fasta

Sources (checked in priority order):
  1. Previously consolidated genomes/ dir (from earlier run)
  2. NCBI local downloads:  pertussis_data/bp_genomes_qc/ncbi_dataset/data/<GCA_*>/*.fna
  3. NAS remote downloads:  pertussis_data/pertussis_gene/step1_ingest/outputs/ncbi_downloads/<GCA_*>/*.fna
  4. NAS shovill assemblies: pertussis_data/pertussis_gene/step1_ingest/outputs/assemblies/<batch>/<run>/contigs.fa

Outputs:
    pertussis_data/bp_genomes_qc/assemblies/<GCA_accession>.fasta      — one file per sample
    pertussis_data/bp_genomes_qc/assemblies/_metadata/source_map.tsv   — tracks origin of each file
    pertussis_data/bp_genomes_qc/assemblies/_metadata/missing.txt      — accessions still not available

Usage:
    python modules/step1_ingest/bin/raw_reads/20_build_genome_database.py [--dry-run]
"""

import csv
import os
import shutil
import sys
from pathlib import Path
from collections import Counter

from raw_read_utils import project_module_data_root, project_workflow_root

ROOT = Path(__file__).resolve().parents[4]
STEP1_DATA_ROOT = project_module_data_root("step1_ingest")
WORKFLOW_DATA_ROOT = project_workflow_root()
MANIFEST = WORKFLOW_DATA_ROOT / "manifest" / "manifest.tsv"

# Output
OUT_DIR = ROOT / "pertussis_data" / "bp_genomes_qc" / "assemblies"
META_DIR = OUT_DIR / "_metadata"

# ── Data sources ──
GENOMES_DIR = ROOT / "genomes"                       # previous consolidation (symlinks)
NCBI_LOCAL = STEP1_DATA_ROOT / "bp_genomes_month_ready" / "ncbi_dataset" / "data"
NAS_DL = STEP1_DATA_ROOT / "outputs" / "ncbi_downloads"
NAS_SHOVILL = STEP1_DATA_ROOT / "outputs" / "assemblies"


def _resolve_fasta(path: str) -> str | None:
    """Resolve a path (possibly a symlink) to a real readable file."""
    p = Path(path)
    try:
        real = p.resolve(strict=True)
        if real.is_file() and real.stat().st_size > 0:
            return str(real)
    except (OSError, RuntimeError):
        pass
    return None


def find_in_genomes_dir(acc: str) -> str | None:
    """Check the previously-consolidated genomes/ directory."""
    f = GENOMES_DIR / f"{acc}.fasta"
    return _resolve_fasta(str(f))


def find_ncbi_local(acc: str) -> str | None:
    d = NCBI_LOCAL / acc
    if d.is_dir():
        for f in d.glob("*.fna"):
            r = _resolve_fasta(str(f))
            if r:
                return r
    return None


def find_nas_download(acc: str) -> str | None:
    d = NAS_DL / acc
    if not d.is_dir():
        return None
    for f in d.glob("*.fna"):
        r = _resolve_fasta(str(f))
        if r:
            return r
    return None


def build_shovill_index() -> dict[str, str]:
    idx = {}
    if not NAS_SHOVILL.is_dir():
        return idx
    for path in NAS_SHOVILL.rglob("contigs.fa"):
        try:
            idx[path.parent.name] = str(path)
        except OSError:
            pass
    return idx


def main():
    dry_run = "--dry-run" in sys.argv

    if not MANIFEST.exists():
        print(f"ERROR: {MANIFEST} not found", file=sys.stderr)
        sys.exit(1)

    rows = []
    with open(MANIFEST, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append(row)
    print(f"Manifest: {len(rows)} samples")

    print("Indexing NAS shovill assemblies …")
    shovill = build_shovill_index()
    print(f"  {len(shovill)} run accessions")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    META_DIR.mkdir(parents=True, exist_ok=True)

    stats = Counter()
    report = []

    for row in rows:
        acc = row["assembly_accession"]
        sra = row.get("sra_run_accession", "").strip()
        ena = row.get("ena_run_accession", "").strip()
        dst = OUT_DIR / f"{acc}.fasta"

        # Already done?
        if dst.is_file() and dst.stat().st_size > 0:
            stats["already_present"] += 1
            report.append((acc, "already_present", str(dst)))
            continue

        # Try sources in order
        src, source = None, "missing"

        # 1. genomes/ dir (previous consolidation, may be symlinks)
        src = find_in_genomes_dir(acc)
        if src:
            source = "genomes_dir"
        else:
            # 2. NCBI local
            src = find_ncbi_local(acc)
            if src:
                source = "ncbi_local"
            else:
                # 3. NAS downloaded
                src = find_nas_download(acc)
                if src:
                    source = "nas_download"
                else:
                    # 4. NAS shovill
                    for rid in [sra, ena]:
                        if rid and rid in shovill:
                            r = _resolve_fasta(shovill[rid])
                            if r:
                                src, source = r, "nas_shovill"
                                break

        if src:
            if not dry_run:
                try:
                    shutil.copy2(src, dst)
                except PermissionError:
                    # NAS files owned by another user; skip, will be re-downloaded
                    stats["permission_denied"] += 1
                    report.append((acc, "permission_denied", src))
                    continue
            stats[source] += 1
            report.append((acc, source, src))
        else:
            stats["missing"] += 1
            report.append((acc, "missing", ""))

    # Write source map
    map_path = META_DIR / "source_map.tsv"
    with open(map_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["assembly_accession", "source", "original_path"])
        w.writerows(report)

    # Write missing list
    missing = sorted(set(r[0] for r in report if r[1] == "missing"))
    missing_path = META_DIR / "missing.txt"
    with open(missing_path, "w") as f:
        for a in missing:
            f.write(a + "\n")

    avail = sum(v for k, v in stats.items() if k != "missing")
    mode = "[DRY RUN]" if dry_run else ""
    print(f"\n{'='*55}")
    print(f"  Genome Database Summary {mode}")
    print(f"{'='*55}")
    print(f"  Directory: pertussis_data/bp_genomes_qc/assemblies/")
    print(f"  Format:    <GCA_accession>.fasta\n")
    for k in ["already_present", "genomes_dir", "ncbi_local", "nas_download", "nas_shovill", "permission_denied", "missing"]:
        if stats[k]:
            print(f"  {k:20s}: {stats[k]:>5d}")
    print(f"\n  Total available: {avail}")
    print(f"  Still missing:   {stats['missing']}")
    print(f"\n  {map_path}")
    print(f"  {missing_path}")


if __name__ == "__main__":
    main()
