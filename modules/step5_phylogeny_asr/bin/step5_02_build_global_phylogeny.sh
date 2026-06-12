#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$REPO_ROOT"

STEP1_DATA_ROOT="$(project_module_data_root step1_ingest)"
STEP2_DATA_ROOT="$(project_module_data_root step2_typing)"
STEP5_DATA_ROOT="$(project_module_data_root step5_phylogeny_asr)"

BALANCED_MANIFEST="${BALANCED_MANIFEST:-$STEP5_DATA_ROOT/outputs/bp_phylogeny_manifest_balanced.tsv}"
FULL_MANIFEST="${FULL_MANIFEST:-$STEP5_DATA_ROOT/outputs/bp_phylogeny_manifest_full.tsv}"
GENOME_PATHS_QC="${GENOME_PATHS_QC:-$STEP2_DATA_ROOT/outputs/bp_genome_paths_qc.tsv}"
GENOME_PATHS_ALL="${GENOME_PATHS_ALL:-$STEP2_DATA_ROOT/outputs/bp_genome_paths.tsv}"
RAW_GENOME_PATHS="${RAW_GENOME_PATHS:-$STEP1_DATA_ROOT/outputs/bp_raw_read_step3_genome_paths.tsv}"
ASSEMBLY_ROOT="${ASSEMBLY_ROOT:-$REPO_ROOT/pertussis_data/bp_genomes_qc/assemblies}"
BALANCED_OUT="${BALANCED_OUT:-$STEP5_DATA_ROOT/outputs/bp_global_phylogeny.nwk}"
FULL_OUT="${FULL_OUT:-$STEP5_DATA_ROOT/outputs/bp_global_phylogeny_full_sensitivity.nwk}"
METADATA_OUT="${METADATA_OUT:-$STEP5_DATA_ROOT/outputs/bp_global_phylogeny_run_metadata.tsv}"
KMER_SIZE="${KMER_SIZE:-21}"
SAMPLE_STRIDE="${SAMPLE_STRIDE:-100}"
FEATURE_BINS="${FEATURE_BINS:-512}"

mkdir -p "$(dirname "$BALANCED_OUT")"

project_env_python bio_tools - "$BALANCED_MANIFEST" "$FULL_MANIFEST" "$GENOME_PATHS_QC" "$GENOME_PATHS_ALL" "$RAW_GENOME_PATHS" "$ASSEMBLY_ROOT" "$BALANCED_OUT" "$FULL_OUT" "$METADATA_OUT" "$KMER_SIZE" "$SAMPLE_STRIDE" "$FEATURE_BINS" <<'PY'
from __future__ import annotations

import csv
import hashlib
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import pdist


def normalize_text(value: str) -> str:
    return (value or "").strip()


def load_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"no rows found in {path}")
    return rows


def write_metadata(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "tree_name",
        "manifest_path",
        "output_tree_path",
        "method",
        "kmer_size",
        "sample_stride",
        "feature_bins",
        "manifest_row_count",
        "included_tip_count",
        "excluded_missing_fasta_count",
        "software_versions",
        "run_utc",
        "notes",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def resolve_genome_map(paths_files: list[Path]) -> dict[str, Path]:
    genome_map: dict[str, Path] = {}
    for path in paths_files:
        if not path.exists():
            continue
        for row in load_tsv_rows(path):
            if normalize_text(row.get("status", "")) != "ok":
                continue
            fasta_path = Path(normalize_text(row.get("fasta_path", "")))
            if not fasta_path.exists():
                continue
            for key in (row.get("resolved_accession", ""), row.get("input_accession", "")):
                key = normalize_text(key)
                if key and key not in genome_map:
                    genome_map[key] = fasta_path
    if not genome_map:
        raise ValueError("no usable genome fasta paths resolved from genome path tables")
    return genome_map


def fasta_sequence(path: Path) -> str:
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                chunks.append(line.strip().upper())
    return "".join(chunks)


TRANSLATION = str.maketrans("ACGT", "TGCA")


def canonical_kmer(kmer: str) -> str:
    rc = kmer.translate(TRANSLATION)[::-1]
    return kmer if kmer <= rc else rc


def hashed_feature_vector(sequence: str, kmer_size: int, sample_stride: int, feature_bins: int) -> np.ndarray:
    vector = np.zeros(feature_bins, dtype=np.float32)
    limit = len(sequence) - kmer_size + 1
    if limit <= 0:
        return vector
    for start in range(0, limit, sample_stride):
        kmer = sequence[start : start + kmer_size]
        if "N" in kmer:
            continue
        canonical = canonical_kmer(kmer)
        digest = hashlib.blake2b(canonical.encode("ascii"), digest_size=8).digest()
        vector[int.from_bytes(digest, "little") % feature_bins] += 1.0
    norm = np.linalg.norm(vector)
    if norm > 0:
        vector /= norm
    else:
        # Keep downstream cosine distances finite even for empty or highly degraded inputs.
        vector[0] = 1.0
    return vector


def linkage_to_newick(node, labels: list[str], parent_distance: float, is_root: bool = False) -> str:
    if node.is_leaf():
        branch_length = max(parent_distance - node.dist, 0.0)
        return f"{labels[node.id]}:{branch_length:.8f}"
    left = linkage_to_newick(node.left, labels, node.dist, False)
    right = linkage_to_newick(node.right, labels, node.dist, False)
    if is_root:
        return f"({left},{right});"
    branch_length = max(parent_distance - node.dist, 0.0)
    return f"({left},{right}):{branch_length:.8f}"


def software_versions() -> str:
    versions = [
        f"python={platform.python_version()}",
        f"numpy={np.__version__}",
    ]
    try:
        import scipy

        versions.append(f"scipy={scipy.__version__}")
    except Exception:
        pass
    try:
        iqtree = subprocess.run(["iqtree2", "--version"], check=False, capture_output=True, text=True)
        version_text = " ".join((iqtree.stdout or iqtree.stderr).split())
        if version_text:
            versions.append(f"iqtree2={version_text}")
    except FileNotFoundError:
        versions.append("iqtree2=not_found")
    return ";".join(versions)


def manifest_label(row: dict[str, str]) -> str:
    return normalize_text(row.get("sample_id_canonical", "")) or normalize_text(row.get("current_accession", ""))


def manifest_accession(row: dict[str, str]) -> str:
    return normalize_text(row.get("current_accession", "")) or normalize_text(row.get("assembly_accession", ""))


def build_tree(
    manifest_rows: list[dict[str, str]],
    genome_map: dict[str, Path],
    *,
    assembly_root: Path,
    kmer_size: int,
    sample_stride: int,
    feature_bins: int,
) -> tuple[str, int]:
    labels: list[str] = []
    vectors: list[np.ndarray] = []
    missing = 0
    for row in manifest_rows:
        accession = manifest_accession(row)
        fasta_path = genome_map.get(accession)
        if fasta_path is None:
            fallback = assembly_root / f"{accession}.fasta"
            if fallback.exists():
                fasta_path = fallback
        if fasta_path is None:
            missing += 1
            continue
        labels.append(manifest_label(row))
        sequence = fasta_sequence(fasta_path)
        vectors.append(hashed_feature_vector(sequence, kmer_size, sample_stride, feature_bins))

    if len(labels) < 2:
        raise ValueError("need at least two included genomes to build a tree")

    matrix = np.vstack(vectors)
    condensed = pdist(matrix, metric="cosine")
    condensed = np.nan_to_num(condensed, nan=0.0, posinf=1.0, neginf=1.0)
    linked = linkage(condensed, method="average")
    tree = to_tree(linked, rd=False)
    return linkage_to_newick(tree, labels, tree.dist, True), missing


def main(argv: list[str]) -> int:
    balanced_manifest = Path(argv[1])
    full_manifest = Path(argv[2])
    genome_paths_qc = Path(argv[3])
    genome_paths_all = Path(argv[4])
    raw_genome_paths = Path(argv[5])
    assembly_root = Path(argv[6])
    balanced_out = Path(argv[7])
    full_out = Path(argv[8])
    metadata_out = Path(argv[9])
    kmer_size = int(argv[10])
    sample_stride = int(argv[11])
    feature_bins = int(argv[12])

    genome_map = resolve_genome_map([genome_paths_qc, genome_paths_all, raw_genome_paths])
    version_string = software_versions()
    run_utc = datetime.now(timezone.utc).replace(microsecond=0).isoformat()

    metadata_rows: list[dict[str, str]] = []
    for tree_name, manifest_path, output_path in (
        ("balanced_main_tree", balanced_manifest, balanced_out),
        ("full_sensitivity_tree", full_manifest, full_out),
    ):
        manifest_rows = load_tsv_rows(manifest_path)
        newick, missing = build_tree(
            manifest_rows,
            genome_map,
            assembly_root=assembly_root,
            kmer_size=kmer_size,
            sample_stride=sample_stride,
            feature_bins=feature_bins,
        )
        output_path.write_text(newick + "\n", encoding="utf-8")
        metadata_rows.append(
            {
                "tree_name": tree_name,
                "manifest_path": str(manifest_path),
                "output_tree_path": str(output_path),
                "method": "sampled_canonical_kmer_composition_tree_average_linkage",
                "kmer_size": str(kmer_size),
                "sample_stride": str(sample_stride),
                "feature_bins": str(feature_bins),
                "manifest_row_count": str(len(manifest_rows)),
                "included_tip_count": str(len(manifest_rows) - missing),
                "excluded_missing_fasta_count": str(missing),
                "software_versions": version_string,
                "run_utc": run_utc,
                "notes": "distance=cosine_on_l2_normalized_kmer_bin_vectors;documented_proxy_tree_for_reproducible_PHY-02_bootstrap",
            }
        )

    write_metadata(metadata_out, metadata_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
PY
