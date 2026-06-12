#!/usr/bin/env python3

import argparse
import gzip
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ModuleNotFoundError as exc:
    SeqIO = None
    Seq = Any
    BIOPYTHON_IMPORT_ERROR = exc
else:
    BIOPYTHON_IMPORT_ERROR = None


@dataclass(frozen=True)
class Extracted:
    name: str
    seq: Seq
    description: str


def require_biopython() -> None:
    if BIOPYTHON_IMPORT_ERROR is not None:
        raise SystemExit(
            "ERROR: Biopython is not installed. Recommended install: "
            "`conda install -c conda-forge biopython` (or create the repo env via "
            "`conda env create -f environment.yml`)."
        )


def run(cmd: list[str]) -> None:
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"ERROR: command failed: {' '.join(cmd)}\nSTDERR:\n{proc.stderr}")


def find_gbff(dataset_dir: Path) -> Path:
    # Prefer uncompressed .gbff, else .gbff.gz
    for p in dataset_dir.rglob("*.gbff"):
        if p.is_file() and p.stat().st_size > 0:
            return p
    for p in dataset_dir.rglob("*.gbff.gz"):
        if p.is_file() and p.stat().st_size > 0:
            return p
    raise SystemExit(f"ERROR: no .gbff/.gbff.gz found under: {dataset_dir}")


def iter_records(gbff_path: Path):
    require_biopython()
    if gbff_path.suffix == ".gz":
        with gzip.open(gbff_path, "rt") as fh:
            yield from SeqIO.parse(fh, "genbank")
    else:
        with gbff_path.open("r", encoding="utf-8", errors="ignore") as fh:
            yield from SeqIO.parse(fh, "genbank")


def qualifier_contains(feature, key: str, needle: str) -> bool:
    vals = feature.qualifiers.get(key, [])
    for v in vals:
        if needle.lower() in str(v).lower():
            return True
    return False


def qualifier_equals(feature, key: str, value: str) -> bool:
    vals = feature.qualifiers.get(key, [])
    for v in vals:
        if str(v).strip().lower() == value.strip().lower():
            return True
    return False


def extract_gene(record, gene_name: str, allow_product_contains: str | None = None) -> Extracted:
    best = None
    # Prefer CDS, then gene features
    for ftype in ("CDS", "gene", "rRNA"):
        for feat in record.features:
            if feat.type != ftype:
                continue
            if qualifier_equals(feat, "gene", gene_name) or qualifier_equals(feat, "locus_tag", gene_name):
                seq = feat.extract(record.seq)
                if best is None or len(seq) > len(best.seq):
                    best = Extracted(gene_name, seq, f"{record.id} {ftype} gene={gene_name}")
            elif allow_product_contains and qualifier_contains(feat, "product", allow_product_contains):
                seq = feat.extract(record.seq)
                if best is None or len(seq) > len(best.seq):
                    best = Extracted(gene_name, seq, f"{record.id} {ftype} product~{allow_product_contains}")
    if best is None:
        raise SystemExit(f"ERROR: could not find gene in GBFF: {gene_name}")
    return best


def extract_23s_rrna(record) -> Extracted:
    best = None
    for feat in record.features:
        if feat.type != "rRNA":
            continue
        if qualifier_contains(feat, "product", "23S"):
            seq = feat.extract(record.seq)
            if best is None or len(seq) > len(best.seq):
                best = Extracted("23S_rRNA", seq, f"{record.id} rRNA product~23S")
    if best is None:
        raise SystemExit("ERROR: could not find 23S rRNA feature in GBFF")
    return best


def extract_ptxP_promoter(record, upstream_len: int) -> Extracted:
    # Approximation: take upstream_len bases upstream of ptxA start, in transcription direction.
    ptxA_feat = None
    for feat in record.features:
        if feat.type not in {"CDS", "gene"}:
            continue
        if qualifier_equals(feat, "gene", "ptxA") or qualifier_contains(feat, "product", "pertussis toxin subunit 1"):
            ptxA_feat = feat
            break
    if ptxA_feat is None:
        raise SystemExit("ERROR: could not find ptxA in GBFF; cannot derive ptxP promoter")

    strand = int(getattr(ptxA_feat.location, "strand", 1) or 1)
    start = int(ptxA_feat.location.start)
    end = int(ptxA_feat.location.end)
    seq = record.seq
    if strand >= 0:
        lo = max(0, start - upstream_len)
        hi = start
        prom = seq[lo:hi]
    else:
        lo = end
        hi = min(len(seq), end + upstream_len)
        prom = seq[lo:hi].reverse_complement()
    return Extracted("ptxP_promoter", prom, f"{record.id} upstream_len={upstream_len} from ptxA")


def write_fasta(path: Path, name: str, seq: Seq, desc: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    s = str(seq).upper()
    if not s:
        raise SystemExit(f"ERROR: empty sequence for {name}")
    header = f">{name} {desc}".strip()
    # wrap 80
    lines = [header]
    for i in range(0, len(s), 80):
        lines.append(s[i : i + 80])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="Build marker reference FASTAs from an annotated reference assembly (GBFF)")
    ap.add_argument(
        "--ref-accession",
        default="GCF_000195715.1",
        help="Reference assembly accession to download (default: GCF_000195715.1; adjust if download fails)",
    )
    ap.add_argument(
        "--workdir",
        default="outputs/_ref",
        help="Work directory for reference download/extraction (default: outputs/_ref)",
    )
    ap.add_argument(
        "--markers-outdir",
        default="references/markers",
        help="Where to write marker FASTAs (default: references/markers)",
    )
    ap.add_argument(
        "--out-23s",
        default="references/23S_rRNA.fasta",
        help="Where to write 23S FASTA (default: references/23S_rRNA.fasta)",
    )
    ap.add_argument(
        "--ptxp-upstream-len",
        type=int,
        default=500,
        help="Length upstream of ptxA to use as ptxP promoter query (default: 500)",
    )
    args = ap.parse_args()
    require_biopython()

    # Ensure datasets exists
    if subprocess.run(["bash", "-lc", "command -v datasets"], capture_output=True).returncode != 0:
        raise SystemExit("ERROR: datasets CLI not found")

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    zip_path = workdir / f"{args.ref_accession}.zip"
    out_dir = workdir / args.ref_accession
    out_dir.mkdir(parents=True, exist_ok=True)

    if not zip_path.exists() or zip_path.stat().st_size == 0:
        print(f"[Ref] downloading: {args.ref_accession}")
        run(
            [
                "datasets",
                "download",
                "genome",
                "accession",
                str(args.ref_accession),
                "--include",
                "genome,gbff",
                "--filename",
                str(zip_path),
            ]
        )

    # Extract
    if not any(out_dir.iterdir()):
        print(f"[Ref] extracting: {zip_path}")
        run(["unzip", "-q", "-o", str(zip_path), "-d", str(out_dir)])

    data_dir = out_dir / "ncbi_dataset" / "data"
    if not data_dir.exists():
        raise SystemExit(f"ERROR: unexpected datasets layout, missing: {data_dir}")

    gbff = find_gbff(data_dir)
    print(f"[Ref] GBFF: {gbff}")

    # Parse first record as main chromosome; if multiple, we still search all records.
    records = list(iter_records(gbff))
    if not records:
        raise SystemExit(f"ERROR: no records parsed from GBFF: {gbff}")

    # Search across records and pick best match per marker
    best_prn = None
    best_fim2 = None
    best_fim3 = None
    best_ptxp = None
    best_23s = None
    for rec in records:
        if best_prn is None:
            try:
                best_prn = extract_gene(rec, "prn", allow_product_contains="pertactin")
            except SystemExit:
                pass
        if best_fim2 is None:
            try:
                best_fim2 = extract_gene(rec, "fim2", allow_product_contains="fimbrial")
            except SystemExit:
                pass
        if best_fim3 is None:
            try:
                best_fim3 = extract_gene(rec, "fim3", allow_product_contains="fimbrial")
            except SystemExit:
                pass
        if best_ptxp is None:
            try:
                best_ptxp = extract_ptxP_promoter(rec, upstream_len=int(args.ptxp_upstream_len))
            except SystemExit:
                pass
        if best_23s is None:
            try:
                best_23s = extract_23s_rrna(rec)
            except SystemExit:
                pass

    if best_prn is None or best_fim2 is None or best_fim3 is None or best_ptxp is None or best_23s is None:
        missing = [
            n
            for n, v in [
                ("prn", best_prn),
                ("fim2", best_fim2),
                ("fim3", best_fim3),
                ("ptxP_promoter", best_ptxp),
                ("23S_rRNA", best_23s),
            ]
            if v is None
        ]
        raise SystemExit(
            "ERROR: failed to extract some markers from reference GBFF: " + ",".join(missing) +
            "\nTry a different --ref-accession. You can list candidates with: datasets summary genome taxon 'Bordetella pertussis' | head"
        )

    markers_out = Path(args.markers_outdir)
    write_fasta(markers_out / "prn_maker.fasta", "prn", best_prn.seq, f"ref={args.ref_accession} {best_prn.description}")
    write_fasta(markers_out / "fim2.fasta", "fim2", best_fim2.seq, f"ref={args.ref_accession} {best_fim2.description}")
    write_fasta(markers_out / "fim3.fasta", "fim3", best_fim3.seq, f"ref={args.ref_accession} {best_fim3.description}")
    write_fasta(
        markers_out / "ptxP_promoter.fasta",
        "ptxP_promoter",
        best_ptxp.seq,
        f"ref={args.ref_accession} {best_ptxp.description}",
    )

    write_fasta(Path(args.out_23s), "23S_rRNA", best_23s.seq, f"ref={args.ref_accession} {best_23s.description}")

    print("[Ref] wrote marker FASTAs:")
    for p in [
        markers_out / "prn_maker.fasta",
        markers_out / "fim2.fasta",
        markers_out / "fim3.fasta",
        markers_out / "ptxP_promoter.fasta",
        Path(args.out_23s),
    ]:
        print(" -", p)


if __name__ == "__main__":
    main()
