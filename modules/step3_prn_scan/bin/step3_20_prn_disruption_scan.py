#!/usr/bin/env python3

import argparse
import concurrent.futures as cf
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class Hsp:
    pident: float
    length: int
    qlen: int
    qstart: int
    qend: int
    bitscore: float


def norm(s: pd.Series) -> pd.Series:
    return s.astype(str).where(~s.isna(), "NA").replace({"nan": "NA", "None": "NA"})


def norm_nonempty(s: pd.Series) -> pd.Series:
    out = norm(s).str.strip()
    return out.where(out != "", "NA")


def normalize_input_table_schema(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize legacy Step2 tables and current QC manifests to scanner columns."""
    df = df.copy()
    resolved = pd.Series(["NA"] * len(df), index=df.index, dtype="object")
    for c in ["genome_resolved_accession", "current_accession", "assembly_accession"]:
        if c not in df.columns:
            continue
        values = norm_nonempty(df[c])
        unresolved = resolved == "NA"
        resolved = resolved.where(~unresolved, values)

    if (resolved == "NA").all():
        raise SystemExit(
            "ERROR: missing required accession column in table: expected genome_resolved_accession, "
            "current_accession, or assembly_accession"
        )
    df["genome_resolved_accession"] = resolved

    if "genome_status" in df.columns:
        df["genome_status"] = norm_nonempty(df["genome_status"])
    else:
        df["genome_status"] = "ok"

    return df


def resolve_fasta_path(raw_path: str, genome_paths_tsv: Path) -> str:
    """Resolve a FASTA path from a genome-paths table.

    Current QC path tables use absolute FASTA paths. This fallback keeps older
    relative tables interpretable while still allowing validation to catch stale
    pre-migration paths.
    """
    path_text = str(raw_path).strip()
    path = Path(path_text)
    if path.is_absolute():
        return str(path)

    candidates = [
        genome_paths_tsv.parent / path,
        Path.cwd() / path,
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate.resolve())
    return str(candidates[0].resolve())


def read_genome_fasta_paths(genome_paths_tsv: Path) -> dict[str, str]:
    """Return mapping: resolved_accession -> fasta_path for status=ok."""
    if not genome_paths_tsv.exists():
        raise SystemExit(f"ERROR: genome paths TSV not found: {genome_paths_tsv}")
    df = pd.read_csv(genome_paths_tsv, sep="\t", dtype=str)
    for c in ["resolved_accession", "fasta_path", "status"]:
        if c not in df.columns:
            raise SystemExit(f"ERROR: genome paths TSV missing column: {c}")
    df["status"] = norm(df["status"])
    df["resolved_accession"] = norm(df["resolved_accession"])
    df["fasta_path"] = norm(df["fasta_path"])
    df = df[df["status"] == "ok"].copy()
    df = df[df["resolved_accession"] != "NA"].copy()
    df = df[df["fasta_path"] != "NA"].copy()
    df["fasta_path"] = df["fasta_path"].map(lambda p: resolve_fasta_path(str(p), genome_paths_tsv))
    # de-dup keeping first
    df = df.drop_duplicates(subset=["resolved_accession"], keep="first")
    return dict(zip(df["resolved_accession"], df["fasta_path"]))


def validate_requested_genome_paths(df: pd.DataFrame, acc_to_fa: dict[str, str], genome_paths_tsv: Path) -> None:
    requested = sorted({str(v) for v in norm(df["genome_resolved_accession"]).tolist() if str(v) != "NA"})
    missing_from_table = [acc for acc in requested if acc not in acc_to_fa]
    if missing_from_table:
        examples = ", ".join(missing_from_table[:10])
        raise SystemExit(
            f"ERROR: genome paths TSV lacks FASTA paths for {len(missing_from_table)} input genomes "
            f"with genome_status=ok: {genome_paths_tsv}. Examples: {examples}. "
            "This usually indicates a stale pre-migration path table; use bp_genome_paths_qc.tsv."
        )

    missing_fastas: list[str] = []
    for acc in requested:
        fasta = Path(acc_to_fa[acc])
        if not fasta.exists() or fasta.stat().st_size == 0:
            missing_fastas.append(f"{acc}={fasta}")
    if missing_fastas:
        examples = "; ".join(missing_fastas[:10])
        raise SystemExit(
            f"ERROR: genome paths TSV points to {len(missing_fastas)} missing/empty FASTA files: "
            f"{genome_paths_tsv}. Examples: {examples}"
        )


def parse_hsps(raw: str, min_pident: float) -> list[Hsp]:
    hsps: list[Hsp] = []
    for line in raw.splitlines():
        parts = line.rstrip("\n").split("\t")
        # outfmt: pident length qlen qstart qend bitscore
        if len(parts) != 6:
            continue
        try:
            pident = float(parts[0])
            if pident < min_pident:
                continue
            hsps.append(
                Hsp(
                    pident=pident,
                    length=int(float(parts[1])),
                    qlen=int(float(parts[2])),
                    qstart=int(float(parts[3])),
                    qend=int(float(parts[4])),
                    bitscore=float(parts[5]),
                )
            )
        except Exception:
            continue
    return hsps


def union_coverage_pct(hsps: list[Hsp]) -> float:
    if not hsps:
        return 0.0
    qlen = hsps[0].qlen
    if qlen <= 0:
        return 0.0
    intervals = []
    for h in hsps:
        a, b = sorted((h.qstart, h.qend))
        intervals.append((a, b))
    intervals.sort()
    merged = []
    cur_a, cur_b = intervals[0]
    for a, b in intervals[1:]:
        if a <= cur_b + 1:
            cur_b = max(cur_b, b)
        else:
            merged.append((cur_a, cur_b))
            cur_a, cur_b = a, b
    merged.append((cur_a, cur_b))
    covered = sum(b - a + 1 for a, b in merged)
    return 100.0 * (covered / qlen)


def _interval_len(a: int, b: int) -> int:
    a2, b2 = sorted((a, b))
    return max(0, b2 - a2 + 1)


def greedy_union_coverage(hsps: list[Hsp]) -> tuple[float, int]:
    """Return (union_qcov_pct, contributing_hsp_n) using a greedy bitscore ordering.

    Counts only HSPs that add at least 1 new bp of query coverage.
    """
    if not hsps:
        return 0.0, 0
    qlen = hsps[0].qlen
    if qlen <= 0:
        return 0.0, 0

    # Sort by alignment quality; prefer long/high bitscore HSPs first
    hsps_sorted = sorted(hsps, key=lambda h: (h.bitscore, h.pident, h.length), reverse=True)

    intervals: list[tuple[int, int]] = []
    contributing = 0

    def add_interval(intervals_in: list[tuple[int, int]], new_a: int, new_b: int) -> tuple[list[tuple[int, int]], int]:
        a, b = sorted((new_a, new_b))
        if b < a:
            return intervals_in, 0
        before_cov = sum(bb - aa + 1 for aa, bb in intervals_in)
        intervals2 = intervals_in + [(a, b)]
        intervals2.sort()
        merged: list[tuple[int, int]] = []
        cur_a, cur_b = intervals2[0]
        for aa, bb in intervals2[1:]:
            if aa <= cur_b + 1:
                cur_b = max(cur_b, bb)
            else:
                merged.append((cur_a, cur_b))
                cur_a, cur_b = aa, bb
        merged.append((cur_a, cur_b))
        after_cov = sum(bb - aa + 1 for aa, bb in merged)
        return merged, max(0, after_cov - before_cov)

    for h in hsps_sorted:
        intervals, added = add_interval(intervals, h.qstart, h.qend)
        if added > 0:
            contributing += 1
        covered = sum(bb - aa + 1 for aa, bb in intervals)
        if 100.0 * (covered / qlen) >= 95.0:
            break

    covered = sum(bb - aa + 1 for aa, bb in intervals)
    return 100.0 * (covered / qlen), contributing


def call_from_hsps(hsps: list[Hsp]) -> tuple[str, dict[str, str]]:
    if not hsps:
        return "no_hit", {}

    best = max(hsps, key=lambda h: (h.bitscore, h.pident, h.length))
    best_single_qcov = 0.0
    if best.qlen > 0:
        best_single_qcov = 100.0 * (_interval_len(best.qstart, best.qend) / best.qlen)

    # Use greedy union coverage so small spurious HSPs don't force a multi-HSP call.
    qcov, contributing_hsp_n = greedy_union_coverage(hsps)

    # Simple interpretation:
    # - intact: essentially full-length in a single HSP
    # - disrupted: full-length but split into multiple HSPs (suggests insertion/deletion/rearrangement)
    # - partial: incomplete coverage
    if best_single_qcov >= 95.0:
        call = "intact"
    elif qcov >= 95.0 and contributing_hsp_n > 1:
        call = "disrupted_multi_hsp"
    elif qcov >= 50.0:
        call = "partial"
    else:
        call = "low_coverage"

    meta = {
        "prn_call": call,
        "prn_hsp_n": str(len(hsps)),
        "prn_hsp_n_contrib": str(contributing_hsp_n),
        "prn_qcov_union_pct": f"{qcov:.2f}",
        "prn_best_single_qcov_pct": f"{best_single_qcov:.2f}",
        "prn_best_pident": f"{best.pident:.2f}",
        "prn_best_bitscore": f"{best.bitscore:.2f}",
    }
    return call, meta


def blast_prn(query_fa: Path, subject_fa: Path, blast_threads: int, max_targets: int) -> str:
    outfmt = "6 pident length qlen qstart qend bitscore"
    cmd = [
        "blastn",
        "-query",
        str(query_fa),
        "-subject",
        str(subject_fa),
        "-outfmt",
        outfmt,
        "-num_threads",
        str(max(1, int(blast_threads))),
        "-max_target_seqs",
        str(max_targets),
        "-task",
        "blastn",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError((proc.stderr or "").strip() or "blastn failed")
    return proc.stdout or ""


def scan_one(
    resolved_acc: str,
    fasta_path: str,
    year: str,
    country: str,
    mlst_st: str,
    query_fa: str,
    blast_threads: int,
    max_targets: int,
    min_pident: float,
) -> dict[str, str]:
    p = Path(fasta_path)
    if not p.exists() or p.stat().st_size == 0:
        return {
            "genome_resolved_accession": resolved_acc,
            "genome_fasta_path": fasta_path,
            "year": year,
            "country": country,
            "mlst_st": mlst_st,
            "prn_call": "missing_fasta",
        }

    try:
        raw = blast_prn(Path(query_fa), p, blast_threads=blast_threads, max_targets=max_targets)
        hsps = parse_hsps(raw, min_pident=min_pident)
        _, meta = call_from_hsps(hsps)
        return {
            "genome_resolved_accession": resolved_acc,
            "genome_fasta_path": fasta_path,
            "year": year,
            "country": country,
            "mlst_st": mlst_st,
            **(meta or {"prn_call": "no_hit"}),
        }
    except Exception as e:
        return {
            "genome_resolved_accession": resolved_acc,
            "genome_fasta_path": fasta_path,
            "year": year,
            "country": country,
            "mlst_st": mlst_st,
            "prn_call": "error",
            "error": str(e),
        }


def main() -> None:
    ap = argparse.ArgumentParser(description="Step3C: scan prn disruption by BLAST HSP fragmentation and union query coverage")
    ap.add_argument("--table", required=True, help="QC manifest or merged Step2 table used for metadata + accession")
    ap.add_argument("--genome-paths", required=True, help="step2_typing/outputs/bp_genome_paths_qc.tsv")
    ap.add_argument("--prn-query", required=True, help="prn query FASTA")
    ap.add_argument("--out", required=True, help="Output calls TSV")
    ap.add_argument("--out-merged", required=True, help="Output merged table with prn_call columns")
    ap.add_argument("--jobs", type=int, default=40, help="Parallelism across genomes")
    ap.add_argument("--executor", choices=["thread", "process"], default="thread")
    ap.add_argument("--blast-threads", type=int, default=1, help="blastn -num_threads per task")
    ap.add_argument("--max-targets", type=int, default=50, help="blastn max_target_seqs")
    ap.add_argument("--min-pident", type=float, default=90.0, help="Minimum pident for HSP inclusion")
    args = ap.parse_args()

    table = Path(args.table)
    genome_paths = Path(args.genome_paths)
    query = Path(args.prn_query)
    out = Path(args.out)
    out_merged = Path(args.out_merged)

    if subprocess.run(["bash", "-lc", "command -v blastn"], capture_output=True).returncode != 0:
        raise SystemExit("ERROR: blastn not found in PATH")
    if not query.exists() or query.stat().st_size == 0:
        raise SystemExit(f"ERROR: prn query missing/empty: {query}")

    df = normalize_input_table_schema(pd.read_csv(table, sep="\t", dtype=str))

    for c in ["year", "country", "mlst_st"]:
        if c in df.columns:
            df[c] = norm(df[c])

    df = df[df["genome_status"] == "ok"].copy()
    df = df[df["genome_resolved_accession"] != "NA"].copy()

    acc_to_fa = read_genome_fasta_paths(genome_paths)
    validate_requested_genome_paths(df, acc_to_fa, genome_paths)

    # Resolve fasta path
    df["genome_fasta_path"] = df["genome_resolved_accession"].map(lambda a: acc_to_fa.get(str(a), "NA"))

    # Build one record per resolved accession (carry along metadata for stratified summaries)
    meta_cols = [c for c in ["year", "country", "mlst_st"] if c in df.columns]

    def first_non_na(s: pd.Series) -> str:
        s2 = norm(s)
        for v in s2.tolist():
            if v != "NA":
                return str(v)
        return "NA"

    agg_spec = {"genome_fasta_path": "first"}
    for c in meta_cols:
        agg_spec[c] = first_non_na

    uniq = df.groupby("genome_resolved_accession", as_index=False).agg(agg_spec)
    if (uniq["genome_fasta_path"] == "NA").any():
        missing_n = int((uniq["genome_fasta_path"] == "NA").sum())
        raise SystemExit(f"ERROR: internal genome path resolution failure for {missing_n} genomes")
    if uniq.empty:
        raise SystemExit("ERROR: no genomes selected for prn disruption scan after filtering genome_status=ok")

    jobs = max(1, int(args.jobs))
    jobs = min(jobs, len(uniq))

    Executor = cf.ThreadPoolExecutor if args.executor == "thread" else cf.ProcessPoolExecutor

    rows: list[dict[str, str]] = []
    with Executor(max_workers=jobs) as ex:
        futs = [
            ex.submit(
                scan_one,
                str(r["genome_resolved_accession"]),
                str(r["genome_fasta_path"]),
                str(r.get("year", "NA")),
                str(r.get("country", "NA")),
                str(r.get("mlst_st", "NA")),
                str(query),
                int(args.blast_threads),
                int(args.max_targets),
                float(args.min_pident),
            )
            for _, r in uniq.iterrows()
        ]
        for i, fut in enumerate(cf.as_completed(futs), start=1):
            rows.append(fut.result())
            if i % 100 == 0:
                print(f"[prn] finished {i}/{len(futs)}", flush=True)

    calls = pd.DataFrame(rows)
    out.parent.mkdir(parents=True, exist_ok=True)
    calls.to_csv(out, sep="\t", index=False)

    # Merge back into the input table after normalizing legacy/current schemas.
    merged = normalize_input_table_schema(pd.read_csv(table, sep="\t", dtype=str))
    calls_keyed = calls.drop_duplicates(subset=["genome_resolved_accession"], keep="first")
    # Avoid clobbering/overlap with existing columns in the main table
    calls_keyed = calls_keyed.drop(columns=["genome_fasta_path"], errors="ignore")
    merged = merged.merge(calls_keyed, on=["genome_resolved_accession"], how="left")
    out_merged.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_merged, sep="\t", index=False)

    print(f"Wrote: {out}")
    print(f"Wrote: {out_merged}")


if __name__ == "__main__":
    main()
