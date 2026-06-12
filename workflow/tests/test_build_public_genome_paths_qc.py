from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "workflow" / "lib" / "build_public_genome_paths_qc.py"


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_build_public_genome_paths_qc_from_public_manifest(tmp_path: Path) -> None:
    manifest = tmp_path / "bp_public_genome_qc_manifest.tsv"
    assembly_root = tmp_path / "bp_genomes_qc" / "assemblies"
    output = tmp_path / "bp_genome_paths_qc.tsv"
    summary = tmp_path / "bp_genome_paths_qc_summary.json"

    assembly_root.mkdir(parents=True)
    fasta_1 = assembly_root / "GCA_000001.1.fasta"
    fasta_2 = assembly_root / "GCA_000002.1.fasta"
    fasta_1.write_text(">public_1\nACGT\n", encoding="utf-8")
    fasta_2.write_text(">public_2\nTGCA\n", encoding="utf-8")

    write_tsv(
        manifest,
        ["sample_id_canonical", "assembly_accession", "current_accession"],
        [
            {
                "sample_id_canonical": "sample_1",
                "assembly_accession": "GCA_000001.1",
                "current_accession": "GCA_000001.1",
            },
            {
                "sample_id_canonical": "sample_2",
                "assembly_accession": "GCA_000002.1",
                "current_accession": "GCA_000002.1",
            },
        ],
    )

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--manifest",
            str(manifest),
            "--assembly-root",
            str(assembly_root),
            "--output",
            str(output),
            "--summary-output",
            str(summary),
        ],
        check=True,
    )

    rows = list(csv.DictReader(output.open(encoding="utf-8"), delimiter="\t"))
    assert len(rows) == 2
    assert rows[0]["status"] == "ok"
    assert rows[1]["status"] == "ok"
    assert rows[0]["fasta_path"] == str(fasta_1)
    assert rows[1]["fasta_path"] == str(fasta_2)

    payload = json.loads(summary.read_text(encoding="utf-8"))
    assert payload["manifest_rows"] == 2
    assert payload["public_path_rows"] == 2
    assert payload["missing_path_rows"] == 0
