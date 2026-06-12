#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def norm(s: pd.Series) -> pd.Series:
    return s.astype(str).where(~s.isna(), "NA").replace({"nan": "NA", "None": "NA"})


def fasta_iter(path: Path):
    name = None
    seq_chunks = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_chunks)
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if name is not None:
            yield name, "".join(seq_chunks)


def extract_contig_subseq(fasta_path: Path, contig_id: str, start_1: int, end_1: int) -> str | None:
    """Extract 1-based inclusive subsequence from contig_id."""
    for name, seq in fasta_iter(fasta_path):
        if name == contig_id:
            if start_1 < 1:
                start_1 = 1
            if end_1 > len(seq):
                end_1 = len(seq)
            if end_1 < start_1:
                return None
            return seq[start_1 - 1 : end_1]
    return None


def main() -> None:
    ap = argparse.ArgumentParser(description="Extract insertion-like prn gap sequences (with flanks) from Step3F evidence")
    ap.add_argument("--evidence", required=True, help="bp_prn_breakpoint_evidence.tsv")
    ap.add_argument("--out-fasta", required=True, help="Output FASTA of extracted gap+flanks")
    ap.add_argument("--out-tsv", required=True, help="Output TSV metadata for extracted sequences")
    ap.add_argument("--flank", type=int, default=200, help="Flank bp to include on each side")
    ap.add_argument("--min-gap", type=int, default=50, help="Minimum bp_max_subject_gap to extract")
    args = ap.parse_args()

    ev = pd.read_csv(Path(args.evidence), sep="\t", dtype=str)
    for c in [
        "genome_resolved_accession",
        "genome_fasta_path",
        "bp_category",
        "bp_contig_id",
        "bp_gap_start",
        "bp_gap_end",
        "bp_max_subject_gap",
    ]:
        if c not in ev.columns:
            raise SystemExit(f"ERROR: evidence missing column: {c}")

    for c in ["year", "country", "mlst_st"]:
        if c not in ev.columns:
            ev[c] = "NA"

    for c in ev.columns:
        ev[c] = norm(ev[c])

    sub = ev[ev["bp_category"] == "insertion_like"].copy()
    sub = sub[sub["bp_gap_start"] != "NA"]
    sub = sub[sub["bp_gap_end"] != "NA"]
    sub["gap_len"] = pd.to_numeric(sub["bp_max_subject_gap"], errors="coerce")
    sub = sub[sub["gap_len"].notna()].copy()
    sub = sub[sub["gap_len"] >= int(args.min_gap)].copy()

    out_fa = Path(args.out_fasta)
    out_tsv = Path(args.out_tsv)
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    records = []
    with out_fa.open("w", encoding="utf-8") as fa:
        for _, r in sub.iterrows():
            acc = str(r["genome_resolved_accession"])
            fasta_path = Path(str(r["genome_fasta_path"]))
            contig = str(r["bp_contig_id"])
            gap_start = int(float(r["bp_gap_start"]))
            gap_end = int(float(r["bp_gap_end"]))
            flank = int(args.flank)
            start = max(1, gap_start - flank)
            end = gap_end + flank

            if not fasta_path.exists() or fasta_path.stat().st_size == 0:
                continue

            seq = extract_contig_subseq(fasta_path, contig, start, end)
            if not seq:
                continue

            header = f"{acc}|{contig}:{start}-{end}|gap:{gap_start}-{gap_end}|gaplen:{int(r['gap_len'])}|year:{r['year']}|country:{r['country']}|st:{r['mlst_st']}"
            fa.write(f">{header}\n")
            # Wrap at 80 columns
            for i in range(0, len(seq), 80):
                fa.write(seq[i : i + 80] + "\n")

            records.append(
                {
                    "genome_resolved_accession": acc,
                    "genome_fasta_path": str(fasta_path),
                    "contig_id": contig,
                    "extract_start": start,
                    "extract_end": end,
                    "gap_start": gap_start,
                    "gap_end": gap_end,
                    "gap_len": int(r["gap_len"]),
                    "year": r["year"],
                    "country": r["country"],
                    "mlst_st": r["mlst_st"],
                }
            )

    pd.DataFrame(records).to_csv(out_tsv, sep="\t", index=False)
    print(f"Wrote: {out_fa} ({len(records)} sequences)")
    print(f"Wrote: {out_tsv}")


if __name__ == "__main__":
    main()
