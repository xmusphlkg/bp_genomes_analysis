#!/usr/bin/env python3
"""
Runtime environment:
  PROJECT_ENV_KEY: bio_tools
  PROJECT_ENV_NAME: pertussis-bio-tools

T04 / Reads availability: trace reads (SRR/ERR/DRR) for each sample.

Strategy:
1. Use existing linkage from Step5 phylogeny manifest (already has ena_run_accession, sra_run_accession)
2. For samples missing run accessions, query ENA programmatic API by BioSample
3. Output: runs.tsv with all traced runs + reads availability report

Readiness threshold: >=30% of 2,247 genomes must have traceable reads.
"""
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd

try:
    import requests
except ImportError:
    requests = None


ENA_FILEREPORT_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_SEARCH_URL = "https://www.ebi.ac.uk/ena/portal/api/search"


def query_ena_runs_by_biosamples(
    biosamples: list[str],
    max_retries: int = 3,
    timeout_seconds: int = 30,
) -> dict[str, list[dict]]:
    """Query ENA for run accessions linked to one or more BioSamples."""
    if requests is None:
        return {}

    biosamples = [value.strip() for value in biosamples if str(value).strip()]
    if not biosamples:
        return {}

    params = {
        "result": "read_run",
        "query": " OR ".join(f'sample_accession="{biosample}"' for biosample in biosamples),
        "fields": "run_accession,experiment_accession,instrument_platform,"
                  "library_layout,fastq_ftp,fastq_md5,fastq_bytes,"
                  "sample_accession,study_accession",
        "format": "tsv",
        "limit": max(100, len(biosamples) * 20),
    }

    for attempt in range(max_retries):
        try:
            resp = requests.get(ENA_SEARCH_URL, params=params, timeout=timeout_seconds)
            if resp.status_code == 200 and resp.text.strip():
                lines = resp.text.strip().split("\n")
                if len(lines) > 1:
                    header = lines[0].split("\t")
                    grouped: dict[str, list[dict]] = defaultdict(list)
                    for line in lines[1:]:
                        row = dict(zip(header, line.split("\t")))
                        biosample = str(row.get("sample_accession", "")).strip()
                        if biosample:
                            grouped[biosample].append(row)
                    return dict(grouped)
            return {}
        except (requests.RequestException, ConnectionError):
            if attempt < max_retries - 1:
                time.sleep(2 ** attempt)
    return {}


def chunked(values: list[tuple[str, str]], size: int) -> list[list[tuple[str, str]]]:
    """Yield successive fixed-size chunks from a list."""
    if size < 1:
        size = 1
    return [values[i:i + size] for i in range(0, len(values), size)]


def trace_reads(
    manifest_path: str,
    min_reads_pct: float = 30.0,
    ena_batch_size: int = 25,
    progress_every: int = 100,
) -> tuple[pd.DataFrame, dict]:
    """Trace reads availability for all samples in the manifest."""

    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    n_total = len(manifest)

    # ── Phase 1: Use existing linkage ────────────────────────────────────
    runs_records = []

    for _, row in manifest.iterrows():
        sample_id = row["sample_id_canonical"]
        # Check all known run accession fields
        sra_run = str(row.get("sra_run_accession", "")).strip()
        ena_run = str(row.get("ena_run_accession", "")).strip()

        if sra_run and sra_run != "nan" and sra_run != "":
            runs_records.append({
                "sample_id_canonical": sample_id,
                "run_accession": sra_run,
                "source": "step5_sra_linkage",
                "platform": str(row.get("sequencing_tech", "")),
                "reads_available": str(row.get("raw_reads_available", "")),
            })
        if ena_run and ena_run != "nan" and ena_run != "" and ena_run != sra_run:
            runs_records.append({
                "sample_id_canonical": sample_id,
                "run_accession": ena_run,
                "source": "step5_ena_linkage",
                "platform": str(row.get("sequencing_tech", "")),
                "reads_available": str(row.get("raw_reads_available", "")),
            })

    linked_samples = {r["sample_id_canonical"] for r in runs_records}
    n_linked_phase1 = len(linked_samples)
    print(f"Phase 1 (existing linkage): {n_linked_phase1}/{n_total} "
          f"({100*n_linked_phase1/n_total:.1f}%) samples have run accessions")

    # ── Phase 2: Query ENA for unlinked samples ─────────────────────────
    unlinked = manifest[~manifest["sample_id_canonical"].isin(linked_samples)]
    n_to_query = len(unlinked)
    n_found_phase2 = 0

    if n_to_query > 0 and requests is not None:
        biosample_queries: list[tuple[str, str]] = []
        for _, row in unlinked.iterrows():
            biosample = str(row.get("biosample_accession", "")).strip()
            if not biosample or biosample == "nan":
                continue
            biosample_queries.append((biosample, row["sample_id_canonical"]))

        batches = chunked(biosample_queries, ena_batch_size)
        total_batches = len(batches)
        print(
            f"Phase 2: Querying ENA for {len(biosample_queries)} unlinked samples "
            f"in {total_batches} batches (batch_size={ena_batch_size})..."
        )

        queried = 0
        for batch_index, batch in enumerate(batches, start=1):
            results_by_biosample = query_ena_runs_by_biosamples(
                [biosample for biosample, _ in batch]
            )
            for biosample, sample_id in batch:
                for r in results_by_biosample.get(biosample, []):
                    runs_records.append({
                        "sample_id_canonical": sample_id,
                        "run_accession": r.get("run_accession", ""),
                        "source": "ena_api_query",
                        "platform": r.get("instrument_platform", ""),
                        "library_layout": r.get("library_layout", ""),
                        "fastq_ftp": r.get("fastq_ftp", ""),
                        "fastq_md5": r.get("fastq_md5", ""),
                        "fastq_bytes": r.get("fastq_bytes", ""),
                        "study_accession": r.get("study_accession", ""),
                        "reads_available": "true",
                    })
                    n_found_phase2 += 1

            queried += len(batch)
            if queried % progress_every == 0 or queried == len(biosample_queries):
                print(
                    f"  Queried {queried}/{len(biosample_queries)} samples "
                    f"across {batch_index}/{total_batches} batches..."
                )
                time.sleep(0.2)
    elif requests is None:
        print("Phase 2: Skipped (requests library not available)")

    # ── Build output ─────────────────────────────────────────────────────
    runs_df = pd.DataFrame(runs_records)

    # Deduplicate: one row per sample × run
    if len(runs_df) > 0:
        runs_df = runs_df.drop_duplicates(subset=["sample_id_canonical", "run_accession"])

    # Count samples with at least one run
    samples_with_runs = runs_df["sample_id_canonical"].nunique() if len(runs_df) > 0 else 0
    pct_with_runs = 100 * samples_with_runs / n_total if n_total > 0 else 0

    # ── Readiness decision ───────────────────────────────────────────────
    reads_available = pct_with_runs >= min_reads_pct

    reads_report = {
        "assessment": "reads_availability",
        "date": pd.Timestamp.now().isoformat(),
        "total_samples": n_total,
        "samples_with_runs": samples_with_runs,
        "pct_with_runs": round(pct_with_runs, 1),
        "total_runs": len(runs_df),
        "phase1_linked": n_linked_phase1,
        "phase2_found": n_found_phase2,
        "threshold_pct": min_reads_pct,
        "reads_available": reads_available,
        "decision": "READY" if reads_available else "NEEDS_DOWNGRADE",
        "recommendation": (
            "Reads availability meets threshold. Proceed with download and QC pipeline."
            if reads_available else
            f"Only {pct_with_runs:.1f}% have reads (threshold: {min_reads_pct}%). "
            "Consider: (1) strengthen public long-read data mining, "
            "(2) conservative frequency/association claims, "
            "(3) focus mechanism analysis on assembly-confident subset only."
        ),
    }

    if len(runs_df) > 0:
        reads_report["platform_distribution"] = runs_df["platform"].value_counts().to_dict()
        if "library_layout" in runs_df.columns:
            reads_report["layout_distribution"] = (
                runs_df["library_layout"].fillna("unknown").value_counts().to_dict()
            )

    return runs_df, reads_report


# ── Snakemake entry point ────────────────────────────────────────────────────
if "snakemake" in dir():
    runs_df, reads_report = trace_reads(
        manifest_path=snakemake.input.manifest,
        min_reads_pct=snakemake.params.min_reads_pct,
    )
    runs_df.to_csv(snakemake.output.runs_table, sep="\t", index=False)
    with open(snakemake.output.reads_report, "w") as f:
        json.dump(reads_report, f, indent=2, default=str)


# ── CLI entry point ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Trace reads availability")
    parser.add_argument("--manifest", required=True, help="Unified manifest TSV")
    parser.add_argument("--min-pct", type=float, default=30.0, help="Minimum reads %% for readiness")
    parser.add_argument("--ena-batch-size", type=int, default=25, help="Number of BioSamples per ENA query batch")
    parser.add_argument("--progress-every", type=int, default=100, help="Emit progress after this many queried samples")
    parser.add_argument("--out-runs", required=True, help="Output runs TSV")
    parser.add_argument("--out-report", dest="out_report", required=True, help="Output reads report JSON")
    args = parser.parse_args()

    runs_df, reads_report = trace_reads(
        args.manifest,
        args.min_pct,
        ena_batch_size=args.ena_batch_size,
        progress_every=args.progress_every,
    )
    runs_df.to_csv(args.out_runs, sep="\t", index=False)
    with open(args.out_report, "w") as f:
        json.dump(reads_report, f, indent=2, default=str)

    status = "READY" if reads_report["reads_available"] else "NEEDS_DOWNGRADE"
    print(f"\n{'='*60}")
    print(f"Reads availability — {status}")
    print(f"  Samples with runs: {reads_report['samples_with_runs']}/{reads_report['total_samples']} "
          f"({reads_report['pct_with_runs']}%)")
    print(f"  Threshold: {reads_report['threshold_pct']}%")
    print(f"  Decision: {reads_report['decision']}")
    print(f"{'='*60}")
