#!/usr/bin/env python3
"""Build a curated insertion-sequence reference layer for prn mechanism work."""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from step4_02_scan_prn_mechanisms import project_module_data_root


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
REFERENCE_BUILD_ID = "bp_is_reference_v1"
REFERENCE_VERSION = "2026-03-21"


@dataclass(frozen=True)
class ReferenceSpec:
    reference_id: str
    is_element_name: str
    source_accession: str
    source_db: str
    source_definition: str
    source_organism: str
    sequence_scope: str
    scan_role: str
    provenance_note: str
    source_start: int | None = None
    source_end: int | None = None
    source_strand: str = "+"
    feature_product: str = ""

    @property
    def source_url(self) -> str:
        if self.source_start is None or self.source_end is None:
            return f"https://www.ncbi.nlm.nih.gov/nuccore/{self.source_accession}"
        return (
            f"https://www.ncbi.nlm.nih.gov/nuccore/{self.source_accession}"
            f"?report=fasta&from={self.source_start}&to={self.source_end}"
        )


REFERENCE_SPECS = [
    ReferenceSpec(
        reference_id="IS481_M28220_full_element",
        is_element_name="IS481",
        source_accession="M28220.1",
        source_db="NCBI nuccore",
        source_definition="Bordetella pertussis insertion sequence IS481 homolog.",
        source_organism="Bordetella pertussis",
        sequence_scope="full_element_homolog",
        scan_role="primary_blast_subject",
        provenance_note="Classic standalone IS481 homolog accession used as a publication-grade seed reference.",
    ),
    ReferenceSpec(
        reference_id="IS1002_Z54268_tnpA",
        is_element_name="IS1002",
        source_accession="Z54268.1",
        source_db="NCBI nuccore",
        source_definition="B.parapertussis tnpA gene (insertion sequence IS1002).",
        source_organism="Bordetella parapertussis",
        sequence_scope="tnpA_family_marker",
        scan_role="primary_blast_subject",
        provenance_note=(
            "Standalone IS1002 tnpA accession retained as a curated family marker; "
            "organism differs from B. pertussis and is recorded explicitly."
        ),
    ),
    ReferenceSpec(
        reference_id="IS1663_CP085969_marker_cds",
        is_element_name="IS1663",
        source_accession="CP085969.1",
        source_db="NCBI nuccore",
        source_definition="Bordetella pertussis strain FDAARGOS_1543 chromosome, complete genome.",
        source_organism="Bordetella pertussis",
        source_start=169190,
        source_end=170206,
        source_strand="-",
        feature_product="IS110-like element IS1663 family transposase",
        sequence_scope="family_marker_cds",
        scan_role="primary_blast_subject",
        provenance_note=(
            "Marker region extracted from an annotated complete B. pertussis genome feature table "
            "because a clean standalone IS1663 nuccore seed was not available in the planning inputs."
        ),
    ),
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


STEP4_DATA_ROOT = project_module_data_root("step4_prn_validation")


def fetch_sequence(spec: ReferenceSpec, *, tool_name: str) -> str:
    params = {
        "db": "nuccore",
        "id": spec.source_accession,
        "rettype": "fasta",
        "retmode": "text",
        "tool": tool_name,
    }
    if spec.source_start is not None and spec.source_end is not None:
        params["seq_start"] = str(spec.source_start)
        params["seq_stop"] = str(spec.source_end)
        params["strand"] = "2" if spec.source_strand == "-" else "1"

    url = f"{EUTILS_BASE}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=60) as response:
        text = response.read().decode("utf-8")

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"unexpected FASTA response for {spec.reference_id}: {text[:200]!r}")
    sequence = "".join(lines[1:]).upper()
    if not sequence:
        raise ValueError(f"empty sequence returned for {spec.reference_id}")
    return sequence


def wrap_fasta(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def md5_hex(text: str) -> str:
    return hashlib.md5(text.encode("utf-8")).hexdigest()


def sha256_hex(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            handle.write(f"{wrap_fasta(sequence)}\n")


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build the bp_step4 insertion-sequence reference FASTA and provenance table from curated NCBI accessions."
    )
    parser.add_argument(
        "--out-fasta",
        type=Path,
        default=STEP4_DATA_ROOT / "references" / "is_elements" / "bp_is_reference.fasta",
        help="Output FASTA for downstream BLAST-style scanning.",
    )
    parser.add_argument(
        "--out-metadata",
        type=Path,
        default=STEP4_DATA_ROOT / "references" / "is_elements" / "bp_is_reference_metadata.tsv",
        help="Output machine-readable reference provenance table.",
    )
    parser.add_argument(
        "--tool-name",
        default="codex_prn01_builder",
        help="NCBI E-utilities tool name recorded in requests.",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    metadata_rows: list[dict[str, str]] = []
    fasta_records: list[tuple[str, str]] = []

    for spec in REFERENCE_SPECS:
        sequence = fetch_sequence(spec, tool_name=args.tool_name)
        fasta_header = f"{spec.reference_id}|{spec.is_element_name}|{spec.source_accession}"
        fasta_records.append((fasta_header, sequence))
        metadata_rows.append(
            {
                "reference_build_id": REFERENCE_BUILD_ID,
                "reference_version": REFERENCE_VERSION,
                "reference_id": spec.reference_id,
                "is_element_name": spec.is_element_name,
                "source_db": spec.source_db,
                "source_accession": spec.source_accession,
                "source_start": "" if spec.source_start is None else str(spec.source_start),
                "source_end": "" if spec.source_end is None else str(spec.source_end),
                "source_strand": spec.source_strand,
                "source_organism": spec.source_organism,
                "source_definition": spec.source_definition,
                "feature_product": spec.feature_product,
                "sequence_scope": spec.sequence_scope,
                "scan_role": spec.scan_role,
                "sequence_length_bp": str(len(sequence)),
                "sequence_md5": md5_hex(sequence),
                "sequence_sha256": sha256_hex(sequence),
                "source_url": spec.source_url,
                "provenance_note": spec.provenance_note,
            }
        )

    fieldnames = [
        "reference_build_id",
        "reference_version",
        "reference_id",
        "is_element_name",
        "source_db",
        "source_accession",
        "source_start",
        "source_end",
        "source_strand",
        "source_organism",
        "source_definition",
        "feature_product",
        "sequence_scope",
        "scan_role",
        "sequence_length_bp",
        "sequence_md5",
        "sequence_sha256",
        "source_url",
        "provenance_note",
    ]

    write_fasta(args.out_fasta, fasta_records)
    write_tsv(args.out_metadata, metadata_rows, fieldnames)

    print(f"Wrote FASTA: {args.out_fasta}")
    print(f"Wrote metadata: {args.out_metadata}")
    print(f"Reference build id: {REFERENCE_BUILD_ID}")
    print(f"Reference version: {REFERENCE_VERSION}")
    print(f"Reference count: {len(metadata_rows)}")
    print("Reference ids: " + ", ".join(row["reference_id"] for row in metadata_rows))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted", file=sys.stderr)
        raise SystemExit(130)
