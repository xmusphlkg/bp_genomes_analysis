import argparse
import subprocess
from pathlib import Path

def run_to_file(cmd: list[str], outfile: Path) -> None:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed:\n{' '.join(cmd)}\n\nSTDERR:\n{p.stderr.strip()}")
    outfile.write_text(p.stdout, encoding="utf-8")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True, help="JSONL report from datasets summary")
    ap.add_argument("--prefix", default="bp", help="Output prefix (default bp)")
    args = ap.parse_args()

    report = Path(args.report)
    if not report.exists():
        raise SystemExit(f"ERROR: report not found: {report}")

    min_tsv = Path(f"{args.prefix}_min_metadata.tsv")
    ext_tsv = Path(f"{args.prefix}_extended_metadata.tsv")

    # Min TSV (quick inspection)
    cmd_min = [
        "dataformat", "tsv", "genome", "--inputfile", str(report), "--force",
        "--fields",
        "accession,organism-name,assminfo-name,assminfo-submitter,assminfo-level,assminfo-release-date"
    ]
    run_to_file(cmd_min, min_tsv)
    print(f"Wrote: {min_tsv}")

    # Extended TSV (date/geo + useful QC fields)
    cmd_ext = [
        "dataformat", "tsv", "genome", "--inputfile", str(report), "--force",
        "--fields",
        ",".join([
            "accession","current-accession","source_database","organism-name",
            "assminfo-name","assminfo-level","assminfo-status","assminfo-refseq-category",
            "assminfo-release-date","assminfo-bioproject",
            "assminfo-biosample-accession","assminfo-biosample-collection-date","assminfo-biosample-geo-loc-name",
            "assminfo-biosample-lat-lon","assminfo-biosample-host","assminfo-biosample-host-disease",
            "assminfo-biosample-isolation-source","assminfo-biosample-strain","assminfo-biosample-isolate",
            "assminfo-sequencing-tech",
            "assmstats-total-sequence-len","assmstats-gc-percent","assmstats-number-of-contigs","assmstats-contig-n50"
        ])
    ]
    run_to_file(cmd_ext, ext_tsv)
    print(f"Wrote: {ext_tsv}")

if __name__ == "__main__":
    main()
