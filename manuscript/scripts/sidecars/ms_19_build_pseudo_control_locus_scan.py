#!/usr/bin/env python3
"""Build structure-matched pseudo-control marker calls for the PRN specificity audit.

This manuscript-facing sidecar intentionally leaves the core Step2 marker tables
untouched. It preferentially downloads the public NCBI RefSeq Tohama-I GBFF
through a checksum-pinned cache, falls back to the existing local Step2 copy if
the online source is unavailable, scans the same Step2 assembly set with the
existing Step2 BLAST allele extractor, and writes a compact per-BioSample
marker-status table that ms_16 can merge into Supplementary Table 38.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from urllib.request import urlopen

import pandas as pd
from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[3]
DATA_HOME = Path(
    os.environ.get(
        "PERTUSSIS_PROJECT_DATA_ROOT",
        str(ROOT / "pertussis_data" / "pertussis_gene"),
    )
)
STEP2_DIR = DATA_HOME / "step2_typing"
STEP2_CODE_DIR = ROOT / "modules" / "step2_typing"
GBFF_PATH = (
    STEP2_DIR
    / "outputs"
    / "_ref"
    / "GCF_000195715.1"
    / "ncbi_dataset"
    / "data"
    / "GCF_000195715.1"
    / "genomic.gbff"
)
PUBLIC_REFSEQ_ASSEMBLY_NAME = "GCF_000195715.1_ASM19571v1"
PUBLIC_REFSEQ_GBFF_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/715/"
    f"{PUBLIC_REFSEQ_ASSEMBLY_NAME}/{PUBLIC_REFSEQ_ASSEMBLY_NAME}_genomic.gbff.gz"
)
PUBLIC_REFSEQ_GBFF_MD5 = "d7587bcace69be598a07d1cef6835183"
GENOME_PATHS_PATH = STEP2_DIR / "outputs" / "bp_genome_paths_qc.tsv"
STEP2_MERGED_PATH = STEP2_DIR / "outputs" / "bp_qc_merged_mlst_markers.tsv"
STEP2_EXTRACTOR = STEP2_CODE_DIR / "bin" / "step2_08_extract_marker_alleles.py"
ASSEMBLY_FASTA_DIR = ROOT / "pertussis_data" / "bp_genomes_qc" / "assemblies"
OUTDIR = ROOT / "manuscript" / "figure_data" / "pseudo_control_loci"
PUBLIC_REFSEQ_CACHE_DIR = OUTDIR / "reference_cache"
PUBLIC_REFSEQ_CACHE_GBFF = PUBLIC_REFSEQ_CACHE_DIR / f"{PUBLIC_REFSEQ_ASSEMBLY_NAME}_genomic.gbff"
MARKER_DIR = OUTDIR / "marker_references"
ALLELES_PATH = OUTDIR / "pseudo_control_marker_alleles.tsv"
STATUS_PATH = OUTDIR / "pseudo_control_marker_status.tsv"
CANDIDATE_PATH = OUTDIR / "pseudo_control_candidate_loci.tsv"
RESOLVED_GENOME_PATHS_PATH = OUTDIR / "pseudo_control_genome_paths_resolved.tsv"
REFERENCE_AUDIT_PATH = OUTDIR / "pseudo_control_marker_reference_audit.tsv"


@dataclass(frozen=True)
class CandidateLocus:
    locus: str
    label: str
    product_contains: str
    structural_match_class: str
    analysis_role: str
    literature_rationale: str
    source_url: str


@dataclass(frozen=True)
class ReferenceContext:
    gbff_path: Path
    downloaded_by_ms19: bool
    extraction_source: str
    extraction_method: str
    notes: str


CANDIDATES = [
    CandidateLocus(
        locus="brkA",
        label="BrkA",
        product_contains="BrkA",
        structural_match_class="type_V_autotransporter_surface_virulence_factor",
        analysis_role="primary_structure_matched_pseudo_control",
        literature_rationale=(
            "Autotransporter virulence factor included in published virulence-marker panels; "
            "same secretion-family class as PRN but not a licensed aP component."
        ),
        source_url="https://pubmed.ncbi.nlm.nih.gov/15096543/",
    ),
    CandidateLocus(
        locus="tcfA",
        label="TcfA",
        product_contains="TcfA",
        structural_match_class="type_V_autotransporter_colonization_factor",
        analysis_role="primary_structure_matched_pseudo_control",
        literature_rationale=(
            "Autotransporter tracheal colonization factor included in virulence-marker panels; "
            "a closer structural pseudo-control than fimbrial subunits."
        ),
        source_url="https://pubmed.ncbi.nlm.nih.gov/15096543/",
    ),
    CandidateLocus(
        locus="vag8",
        label="Vag8",
        product_contains="Vag8",
        structural_match_class="type_V_autotransporter_surface_virulence_factor",
        analysis_role="primary_structure_matched_pseudo_control",
        literature_rationale=(
            "Autotransporter virulence-associated gene used in pertussis marker-surveillance panels; "
            "tests whether PRN-like loss signals generalize across related surface autotransporters."
        ),
        source_url="https://pubmed.ncbi.nlm.nih.gov/15096543/",
    ),
    CandidateLocus(
        locus="sphB1",
        label="SphB1",
        product_contains="SphB1",
        structural_match_class="type_V_autotransporter_serine_protease",
        analysis_role="primary_structure_matched_pseudo_control",
        literature_rationale=(
            "Serine protease autotransporter with assigned host-interaction/virulence function; "
            "a PRN-sized type-V secretion comparator that is not the primary PRN antigen locus."
        ),
        source_url="https://pubmed.ncbi.nlm.nih.gov/21554944/",
    ),
    CandidateLocus(
        locus="phg",
        label="Phg",
        product_contains="Phg",
        structural_match_class="pertactin_homologous_autotransporter",
        analysis_role="primary_structure_matched_pseudo_control",
        literature_rationale=(
            "Pertactin-homologous autotransporter, providing the most direct family-level pseudo-control "
            "available from the local reference annotation."
        ),
        source_url="https://www.sciencedirect.com/science/article/pii/S0944501305000236",
    ),
    CandidateLocus(
        locus="bapC",
        label="BapC",
        product_contains="BapC",
        structural_match_class="type_V_autotransporter_with_reference_pseudogene_caveat",
        analysis_role="secondary_structure_matched_pseudo_control",
        literature_rationale=(
            "Bordetella autotransporter C is biologically close to PRN/BrkA, but the Tohama-I reference "
            "annotation carries a pseudogene caveat, so this is secondary rather than headline control evidence."
        ),
        source_url="https://pubmed.ncbi.nlm.nih.gov/21554944/",
    ),
    CandidateLocus(
        locus="fhaB",
        label="FHA/FhaB",
        product_contains="filamentous hemagglutinin",
        structural_match_class="large_surface_or_secreted_acellular_vaccine_adhesin",
        analysis_role="secondary_vaccine_antigen_pseudo_control",
        literature_rationale=(
            "FHA is an acellular-vaccine antigen and large adhesin; it is not an autotransporter, but it tests "
            "whether another large vaccine adhesin shows a PRN-like loss pattern on the same assembly frame."
        ),
        source_url="https://wwwnc.cdc.gov/eid/article/22/2/15-1332_article",
    ),
]


def file_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_binary(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(url, timeout=120) as response, destination.open("wb") as handle:
        shutil.copyfileobj(response, handle)


def ensure_public_refseq_gbff() -> Path:
    PUBLIC_REFSEQ_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if PUBLIC_REFSEQ_CACHE_GBFF.exists() and file_md5(PUBLIC_REFSEQ_CACHE_GBFF) == PUBLIC_REFSEQ_GBFF_MD5:
        return PUBLIC_REFSEQ_CACHE_GBFF

    tmp_gz = None
    tmp_gbff = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=PUBLIC_REFSEQ_CACHE_DIR, suffix=".gbff.gz", delete=False
        ) as gz_handle:
            tmp_gz = Path(gz_handle.name)
        download_binary(PUBLIC_REFSEQ_GBFF_URL, tmp_gz)

        tmp_gbff = PUBLIC_REFSEQ_CACHE_GBFF.with_suffix(".gbff.tmp")
        with gzip.open(tmp_gz, "rb") as src, tmp_gbff.open("wb") as dst:
            shutil.copyfileobj(src, dst)

        if file_md5(tmp_gbff) != PUBLIC_REFSEQ_GBFF_MD5:
            raise RuntimeError(
                "downloaded RefSeq GBFF checksum did not match the pinned md5 "
                f"({PUBLIC_REFSEQ_GBFF_MD5})"
            )

        tmp_gbff.replace(PUBLIC_REFSEQ_CACHE_GBFF)
        return PUBLIC_REFSEQ_CACHE_GBFF
    finally:
        if tmp_gz is not None:
            tmp_gz.unlink(missing_ok=True)
        if tmp_gbff is not None and tmp_gbff.exists() and tmp_gbff != PUBLIC_REFSEQ_CACHE_GBFF:
            tmp_gbff.unlink(missing_ok=True)


def resolve_reference_context(reference_source: str) -> ReferenceContext:
    if reference_source in {"auto", "online"}:
        try:
            gbff_path = ensure_public_refseq_gbff()
            return ReferenceContext(
                gbff_path=gbff_path,
                downloaded_by_ms19=True,
                extraction_source="ncbi_refseq_ftp_checksum_pinned",
                extraction_method="Bio.SeqIO GenBank feature extraction from checksum-pinned NCBI RefSeq FTP cache",
                notes=(
                    "Marker FASTA was extracted by ms_19 from the public NCBI RefSeq "
                    "Tohama-I GBFF (checksum-pinned FTP cache); the local Step2 copy is "
                    "retained only as a fallback."
                ),
            )
        except Exception as exc:
            if reference_source == "online":
                raise SystemExit(f"ERROR: public RefSeq download failed: {exc}") from exc
            print(
                f"WARNING: public RefSeq download failed ({exc}); falling back to the local Step2 GBFF.",
                file=sys.stderr,
            )

    if not GBFF_PATH.exists():
        raise SystemExit(f"ERROR: local Step2 GBFF not found: {GBFF_PATH}")

    return ReferenceContext(
        gbff_path=GBFF_PATH,
        downloaded_by_ms19=False,
        extraction_source="local_step2_tohama_i_gbff",
        extraction_method="Bio.SeqIO GenBank feature extraction from local Step2 Tohama-I GBFF",
        notes=(
            "Marker FASTA was extracted by ms_19 from an existing local Tohama-I GBFF "
            "reference; public RefSeq download was not used."
        ),
    )


def clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "na"} else text


def marker_hash(sequence: str) -> str:
    return hashlib.md5(sequence.encode("utf-8")).hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def accession_swap(accession: str) -> str:
    text = clean_text(accession)
    if text.startswith("GCF_"):
        return "GCA_" + text[4:]
    if text.startswith("GCA_"):
        return "GCF_" + text[4:]
    return text


def accession_base(accession: str) -> str:
    return clean_text(accession).split(".", maxsplit=1)[0]


def accession_version(accession: str) -> int:
    text = clean_text(accession)
    if "." not in text:
        return -1
    try:
        return int(text.rsplit(".", maxsplit=1)[1])
    except ValueError:
        return -1


def build_assembly_index() -> tuple[dict[str, Path], dict[str, list[Path]]]:
    exact: dict[str, Path] = {}
    by_base: dict[str, list[Path]] = {}
    if not ASSEMBLY_FASTA_DIR.exists():
        return exact, by_base

    for path in sorted(ASSEMBLY_FASTA_DIR.glob("*.fasta")):
        stem = path.stem
        names = {stem, accession_swap(stem)}
        for name in names:
            exact.setdefault(name, path)
            by_base.setdefault(accession_base(name), []).append(path)

    for base, paths in by_base.items():
        by_base[base] = sorted(
            set(paths),
            key=lambda p: (accession_version(p.stem), p.stem),
            reverse=True,
        )
    return exact, by_base


def resolve_fasta_path(row: pd.Series, exact: dict[str, Path], by_base: dict[str, list[Path]]) -> tuple[str, str, str]:
    original = clean_text(row.get("fasta_path", ""))
    resolved_accession = clean_text(row.get("resolved_accession", ""))
    input_accession = clean_text(row.get("input_accession", ""))

    if original:
        original_path = Path(original)
        for candidate_path in [original_path, STEP2_DIR / original_path, ROOT / original_path]:
            if candidate_path.exists() and candidate_path.stat().st_size > 0:
                return str(candidate_path), "ok", "direct_existing_path"

    accessions = [resolved_accession, accession_swap(resolved_accession), input_accession, accession_swap(input_accession)]
    for accession in accessions:
        if accession in exact and exact[accession].exists():
            return str(exact[accession]), "ok", f"assembly_index_exact:{accession}"

    for accession in accessions:
        base = accession_base(accession)
        if base in by_base and by_base[base]:
            return (
                str(by_base[base][0]),
                "ok",
                f"assembly_index_version_fallback:{accession}->{by_base[base][0].stem}",
            )

    return original, "missing", "unresolved_local_fasta"


def write_resolved_genome_paths() -> pd.DataFrame:
    genome_paths = pd.read_csv(GENOME_PATHS_PATH, sep="\t", dtype=str)
    exact, by_base = build_assembly_index()
    rows: list[dict[str, object]] = []
    for _, row in genome_paths.iterrows():
        original_status = clean_text(row.get("status", ""))
        resolved_path, resolved_status, note = resolve_fasta_path(row, exact, by_base)
        status = "ok" if original_status == "ok" and resolved_status == "ok" else "missing"
        rows.append(
            {
                "input_accession": clean_text(row.get("input_accession", "")),
                "resolved_accession": clean_text(row.get("resolved_accession", "")),
                "status": status,
                "fasta_path": resolved_path if status == "ok" else "",
                "note": note,
                "original_status": original_status,
                "original_fasta_path": clean_text(row.get("fasta_path", "")),
            }
        )
    resolved = pd.DataFrame(rows)
    RESOLVED_GENOME_PATHS_PATH.parent.mkdir(parents=True, exist_ok=True)
    resolved.to_csv(RESOLVED_GENOME_PATHS_PATH, sep="\t", index=False)
    return resolved


def qualifier_matches(feature, key: str, needle: str) -> bool:
    return any(str(value).strip().lower() == needle.lower() for value in feature.qualifiers.get(key, []))


def qualifier_contains(feature, key: str, needle: str) -> bool:
    return any(needle.lower() in str(value).lower() for value in feature.qualifiers.get(key, []))


def find_candidate_feature(records, candidate: CandidateLocus):
    best = None
    best_record_id = ""
    for record in records:
        for feature_type in ("CDS", "gene"):
            for feature in record.features:
                if feature.type != feature_type:
                    continue
                if (
                    qualifier_matches(feature, "gene", candidate.locus)
                    or qualifier_matches(feature, "locus_tag", candidate.locus)
                    or qualifier_matches(feature, "old_locus_tag", candidate.locus)
                    or qualifier_contains(feature, "product", candidate.product_contains)
                ):
                    sequence = str(feature.extract(record.seq)).upper()
                    if best is None or len(sequence) > len(best["sequence"]):
                        best = {
                            "sequence": sequence,
                            "feature_type": feature_type,
                            "gene": ";".join(feature.qualifiers.get("gene", [])),
                            "locus_tag": ";".join(feature.qualifiers.get("locus_tag", [])),
                            "old_locus_tag": ";".join(feature.qualifiers.get("old_locus_tag", [])),
                            "product": ";".join(feature.qualifiers.get("product", [])),
                            "feature_location": str(feature.location),
                            "is_pseudo": "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers,
                            "note": ";".join(feature.qualifiers.get("note", [])),
                        }
                        best_record_id = record.id
    if best is None:
        raise SystemExit(f"ERROR: could not extract pseudo-control locus from GBFF: {candidate.locus}")
    best["record_id"] = best_record_id
    return best


def write_fasta(path: Path, locus: str, sequence: str, description: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f">{locus} {description}".strip()]
    for i in range(0, len(sequence), 80):
        lines.append(sequence[i : i + 80])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def extract_marker_references(reference_ctx: ReferenceContext) -> pd.DataFrame:
    records = list(SeqIO.parse(reference_ctx.gbff_path, "genbank"))
    rows: list[dict[str, object]] = []
    for candidate in CANDIDATES:
        feature = find_candidate_feature(records, candidate)
        sequence = feature["sequence"]
        write_fasta(
            MARKER_DIR / f"{candidate.locus}.fasta",
            candidate.locus,
            sequence,
            f"ref=GCF_000195715.1 {feature['record_id']} {feature['feature_type']} product={feature['product']}",
        )
        rows.append(
            {
                "locus": candidate.locus,
                "locus_label": candidate.label,
                "locus_length_bp": len(sequence),
                "reference_record": feature["record_id"],
                "feature_type": feature["feature_type"],
                "gene": feature["gene"],
                "locus_tag": feature["locus_tag"],
                "old_locus_tag": feature["old_locus_tag"],
                "product": feature["product"],
                "feature_location": feature["feature_location"],
                "is_pseudo": feature["is_pseudo"],
                "note": feature["note"],
                "reference_sequence_md5": marker_hash(sequence),
                "marker_reference_path": display_path(MARKER_DIR / f"{candidate.locus}.fasta"),
                "extraction_source": reference_ctx.extraction_source,
                "structural_match_class": candidate.structural_match_class,
                "analysis_role": candidate.analysis_role,
                "literature_rationale": candidate.literature_rationale,
                "source_url": candidate.source_url,
            }
        )
    out = pd.DataFrame(rows)
    CANDIDATE_PATH.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(CANDIDATE_PATH, sep="\t", index=False)
    return out


def write_reference_audit(candidate_table: pd.DataFrame, reference_ctx: ReferenceContext) -> pd.DataFrame:
    """Record how marker FASTA references were generated and verify sequence identity."""

    candidate_by_locus = {candidate.locus: candidate for candidate in CANDIDATES}
    gbff_records = list(SeqIO.parse(reference_ctx.gbff_path, "genbank"))
    rows: list[dict[str, object]] = []
    for _, candidate_row in candidate_table.iterrows():
        locus = clean_text(candidate_row.get("locus", ""))
        fasta_path = MARKER_DIR / f"{locus}.fasta"
        fasta_records = list(SeqIO.parse(fasta_path, "fasta")) if fasta_path.exists() else []
        fasta_sequence = str(fasta_records[0].seq).upper() if len(fasta_records) == 1 else ""
        fasta_md5 = marker_hash(fasta_sequence) if fasta_sequence else ""

        gbff_feature = find_candidate_feature(gbff_records, candidate_by_locus[locus]) if locus in candidate_by_locus else {}
        gbff_sequence = clean_text(gbff_feature.get("sequence", ""))
        gbff_md5 = marker_hash(gbff_sequence) if gbff_sequence else ""
        candidate_md5 = clean_text(candidate_row.get("reference_sequence_md5", ""))
        matches_candidate = fasta_md5 != "" and fasta_md5 == candidate_md5
        matches_gbff = fasta_md5 != "" and fasta_md5 == gbff_md5
        audit_status = "pass" if len(fasta_records) == 1 and matches_candidate and matches_gbff else "fail"

        rows.append(
            {
                "locus": locus,
                "locus_label": clean_text(candidate_row.get("locus_label", "")),
                "fasta_path": display_path(fasta_path),
                "fasta_exists": fasta_path.exists(),
                "fasta_record_count": len(fasta_records),
                "fasta_length_bp": len(fasta_sequence) if fasta_sequence else "",
                "fasta_md5": fasta_md5,
                "candidate_length_bp": clean_text(candidate_row.get("locus_length_bp", "")),
                "candidate_md5": candidate_md5,
                "gbff_path": display_path(GBFF_PATH),
                "gbff_reference_record": clean_text(gbff_feature.get("record_id", "")),
                "gbff_feature_type": clean_text(gbff_feature.get("feature_type", "")),
                "gbff_gene": clean_text(gbff_feature.get("gene", "")),
                "gbff_locus_tag": clean_text(gbff_feature.get("locus_tag", "")),
                "gbff_old_locus_tag": clean_text(gbff_feature.get("old_locus_tag", "")),
                "gbff_product": clean_text(gbff_feature.get("product", "")),
                "gbff_feature_location": clean_text(gbff_feature.get("feature_location", "")),
                "gbff_is_pseudo": gbff_feature.get("is_pseudo", ""),
                "gbff_note": clean_text(gbff_feature.get("note", "")),
                "gbff_md5": gbff_md5,
                "sequence_matches_candidate_md5": matches_candidate,
                "sequence_matches_gbff_feature": matches_gbff,
                "source_reference_path": display_path(reference_ctx.gbff_path),
                "source_downloaded_by_ms19": "yes" if reference_ctx.downloaded_by_ms19 else "no",
                "extraction_method": reference_ctx.extraction_method,
                "audit_status": audit_status,
                "notes": reference_ctx.notes,
            }
        )
    out = pd.DataFrame(rows)
    REFERENCE_AUDIT_PATH.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(REFERENCE_AUDIT_PATH, sep="\t", index=False)
    return out


def run_marker_scan(jobs: int, min_pident: float, min_qcov: float) -> None:
    cmd = [
        sys.executable,
        str(STEP2_EXTRACTOR),
        "--genome-paths",
        str(RESOLVED_GENOME_PATHS_PATH),
        "--markers-dir",
        str(MARKER_DIR),
        "--out",
        str(ALLELES_PATH),
        "--jobs",
        str(jobs),
        "--blast-threads",
        "1",
        "--min-pident",
        str(min_pident),
        "--min-qcov",
        str(min_qcov),
    ]
    proc = subprocess.run(cmd, cwd=STEP2_DIR, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"ERROR: pseudo-control marker scan failed with exit code {proc.returncode}")


def build_status_table(candidate_table: pd.DataFrame) -> pd.DataFrame:
    metadata = pd.read_csv(STEP2_MERGED_PATH, sep="\t", dtype=str)
    alleles = pd.read_csv(ALLELES_PATH, sep="\t", dtype=str)
    resolved_paths = pd.read_csv(RESOLVED_GENOME_PATHS_PATH, sep="\t", dtype=str)
    alleles["genome_fasta_path"] = alleles["genome_fasta_path"].map(clean_text)
    alleles["marker"] = alleles["marker"].map(clean_text)
    metadata["genome_fasta_path"] = metadata["genome_fasta_path"].map(clean_text)
    metadata["Current Accession"] = metadata["Current Accession"].map(clean_text)
    accession_to_path = (
        resolved_paths.loc[resolved_paths["status"].eq("ok"), ["resolved_accession", "fasta_path"]]
        .drop_duplicates(subset=["resolved_accession"])
        .set_index("resolved_accession")["fasta_path"]
        .to_dict()
    )
    metadata["resolved_genome_fasta_path"] = metadata["Current Accession"].map(accession_to_path).fillna("")

    keep = [
        "Assembly BioSample Accession",
        "Current Accession",
        "genome_fasta_path",
        "resolved_genome_fasta_path",
        "country",
        "year",
        "mlst_st",
    ]
    status = metadata[[column for column in keep if column in metadata.columns]].copy()
    status = status.rename(
        columns={
            "Assembly BioSample Accession": "biosample_accession",
            "Current Accession": "current_accession",
        }
    )
    status["biosample_accession"] = status["biosample_accession"].map(clean_text)
    status = status.loc[status["biosample_accession"] != ""].drop_duplicates(subset=["biosample_accession"])

    for value_col, prefix in [
        ("status", "marker_status_"),
        ("pident", "marker_pident_"),
        ("qcov_pct", "marker_qcov_"),
        ("allele_len", "marker_len_"),
        ("allele_hash", "marker_"),
    ]:
        wide = alleles.pivot_table(
            index="genome_fasta_path",
            columns="marker",
            values=value_col,
            aggfunc="first",
        )
        wide.columns = [f"{prefix}{column}" for column in wide.columns]
        status = status.merge(
            wide.reset_index(),
            left_on="resolved_genome_fasta_path",
            right_on="genome_fasta_path",
            how="left",
            suffixes=("", "_marker_scan"),
        )
        if "genome_fasta_path_marker_scan" in status.columns:
            status = status.drop(columns=["genome_fasta_path_marker_scan"])

    status["pseudo_control_loci_scanned"] = ",".join(candidate_table["locus"].tolist())
    STATUS_PATH.parent.mkdir(parents=True, exist_ok=True)
    status.to_csv(STATUS_PATH, sep="\t", index=False)
    return status


def main() -> None:
    parser = argparse.ArgumentParser(description="Build structure-matched pseudo-control marker scan sidecars.")
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--min-pident", type=float, default=90.0)
    parser.add_argument("--min-qcov", type=float, default=80.0)
    parser.add_argument(
        "--reference-source",
        choices=("auto", "online", "local"),
        default="auto",
        help="Source for the Tohama-I GBFF reference (default: auto, try public NCBI RefSeq then fall back local).",
    )
    parser.add_argument("--skip-scan", action="store_true")
    args = parser.parse_args()

    reference_ctx = resolve_reference_context(args.reference_source)
    candidate_table = extract_marker_references(reference_ctx)
    reference_audit = write_reference_audit(candidate_table, reference_ctx)
    resolved_paths = write_resolved_genome_paths()
    if not args.skip_scan:
        run_marker_scan(jobs=max(1, int(args.jobs)), min_pident=float(args.min_pident), min_qcov=float(args.min_qcov))
    if not ALLELES_PATH.exists() or ALLELES_PATH.stat().st_size == 0:
        raise SystemExit(f"ERROR: marker allele scan is missing or empty: {ALLELES_PATH}")
    status_table = build_status_table(candidate_table)

    print("Built pseudo-control locus sidecars:")
    for path in [CANDIDATE_PATH, REFERENCE_AUDIT_PATH, ALLELES_PATH, STATUS_PATH]:
        print(f" - {display_path(path)}")
    print(f"Reference source: {reference_ctx.extraction_source} ({display_path(reference_ctx.gbff_path)})")
    print(f"Scanned loci: {','.join(candidate_table['locus'].tolist())}")
    print(f"Reference audit pass rows: {int(reference_audit['audit_status'].eq('pass').sum())}/{len(reference_audit)}")
    print(f"Resolved genomes for scan: {int(resolved_paths['status'].eq('ok').sum())}/{len(resolved_paths)}")
    print(f"Status rows: {len(status_table)}")


if __name__ == "__main__":
    main()
