#!/usr/bin/env python3
"""Run matched downsampling stress tests for PRN read validation.

This helper keeps the same Step4 read-validation path, but applies it to
paired-Illumina exemplars at 100%, 50%, 25%, and 10% read-pair fractions with
three replicates per depth. The resulting rows are later folded into the
manuscript-facing detectability detail table.
"""

from __future__ import annotations

import argparse
import gzip
import os
import random
import re
import shutil
import subprocess
import sys
from itertools import zip_longest
from pathlib import Path

import pandas as pd
import pysam
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from step4_02_scan_prn_mechanisms import (
    STEP4_DATA_ROOT,
    WORKFLOW_DATA_ROOT,
    load_tsv_rows,
    normalize_text,
    write_tsv,
)

SCRIPT_DIR = Path(__file__).resolve().parent


TARGET_EXEMPLARS = [
    {
        "family_key": "is481_1043",
        "family_label": "IS481-associated 1,043-bp architecture",
        "candidate_sample_ids": [
            "SAMN03455350",
            "SAMN03249376",
            "SAMN03216671",
            "SAMN03455360",
            "SAMD01593374",
        ],
    },
    {
        "family_key": "rearrangement_family",
        "family_label": "Rearrangement family",
        "candidate_sample_ids": [
            "SAMN03854490",
            "SAMN03249363",
            "SAMN03249368",
            "SAMN03455361",
            "SAMN03216670",
        ],
    },
    {
        "family_key": "other_insertion_like",
        "family_label": "Other insertion-like disruptions",
        "candidate_sample_ids": [
            "SAMN08136991",
            "SAMN06007543",
            "SAMN12585531",
            "SAMN12108514",
            "SAMN12385786",
        ],
    },
]

DEFAULT_FRACTIONS = (1.0, 0.5, 0.25, 0.10)
DEFAULT_REPLICATES = 3
DEFAULT_SEED = 20260419


def normalize_read_name(text: str) -> str:
    token = normalize_text(text).split()[0]
    return re.sub(r"(/[12])$", "", token)


def fraction_token(fraction: float) -> str:
    if fraction >= 0.999999:
        return "100"
    return str(int(round(fraction * 100))).zfill(2)


def selected_indices(total_pairs: int, fraction: float, seed: int) -> set[int]:
    if total_pairs <= 0:
        return set()
    if fraction >= 0.999999:
        return set(range(total_pairs))
    rng = random.Random(seed)
    indices = {index for index in range(total_pairs) if rng.random() < fraction}
    if not indices:
        indices = {rng.randrange(total_pairs)}
    return indices


def count_fastq_pairs(path_r1: Path, path_r2: Path) -> int:
    total = 0
    with gzip.open(path_r1, "rt") as handle_r1, gzip.open(path_r2, "rt") as handle_r2:
        for record_r1, record_r2 in zip_longest(FastqGeneralIterator(handle_r1), FastqGeneralIterator(handle_r2)):
            if record_r1 is None or record_r2 is None:
                raise ValueError(f"FASTQ pair length mismatch: {path_r1} vs {path_r2}")
            total += 1
    return total


def copy_or_symlink(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        dst.symlink_to(src)
    except OSError:
        shutil.copy2(src, dst)


def downsample_fastq_pair(
    path_r1: Path,
    path_r2: Path,
    out_r1: Path,
    out_r2: Path,
    *,
    fraction: float,
    seed: int,
) -> tuple[int, int, set[str]]:
    total_pairs = count_fastq_pairs(path_r1, path_r2)
    keep_indices = selected_indices(total_pairs, fraction, seed)

    selected_names: set[str] = set()
    retained_pairs = 0
    with gzip.open(path_r1, "rt") as handle_r1, gzip.open(path_r2, "rt") as handle_r2:
        with gzip.open(out_r1, "wt", compresslevel=6) as out_handle_r1, gzip.open(
            out_r2, "wt", compresslevel=6
        ) as out_handle_r2:
            for index, (record_r1, record_r2) in enumerate(
                zip_longest(FastqGeneralIterator(handle_r1), FastqGeneralIterator(handle_r2))
            ):
                if record_r1 is None or record_r2 is None:
                    raise ValueError(f"FASTQ pair length mismatch: {path_r1} vs {path_r2}")
                if index not in keep_indices:
                    continue
                title_r1, seq_r1, qual_r1 = record_r1
                title_r2, seq_r2, qual_r2 = record_r2
                name_r1 = normalize_read_name(title_r1)
                name_r2 = normalize_read_name(title_r2)
                if name_r1 != name_r2:
                    raise ValueError(f"FASTQ mate name mismatch at pair {index}: {name_r1} != {name_r2}")
                selected_names.add(name_r1)
                retained_pairs += 1
                out_handle_r1.write(f"@{title_r1}\n{seq_r1}\n+\n{qual_r1}\n")
                out_handle_r2.write(f"@{title_r2}\n{seq_r2}\n+\n{qual_r2}\n")
    return total_pairs, retained_pairs, selected_names


def downsample_bam(path_bam: Path, out_bam: Path, selected_names: set[str]) -> int:
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    if out_bam.exists():
        out_bam.unlink()
    out_bai = out_bam.with_suffix(out_bam.suffix + ".bai")
    if out_bai.exists():
        out_bai.unlink()

    selected_records = 0
    with pysam.AlignmentFile(path_bam, "rb") as in_bam, pysam.AlignmentFile(out_bam, "wb", template=in_bam) as out_handle:
        for record in in_bam.fetch(until_eof=True):
            if normalize_read_name(record.query_name) not in selected_names:
                continue
            out_handle.write(record)
            selected_records += 1

    try:
        pysam.index(str(out_bam))
    except Exception as exc:  # pragma: no cover - only exercised if BAM sorting assumptions break
        raise RuntimeError(f"Failed to index downsampled BAM {out_bam}") from exc
    return selected_records


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate matched downsampling stress tests for Step4 read validation."
    )
    parser.add_argument(
        "--subset",
        type=Path,
        default=STEP4_DATA_ROOT / "outputs" / "bp_prn_validation_subset.tsv",
        help="Validation subset manifest used as the source of exemplar metadata.",
    )
    parser.add_argument(
        "--reads-root",
        type=Path,
        default=WORKFLOW_DATA_ROOT / "reads_clean",
        help="Root containing the original paired FASTQ files.",
    )
    parser.add_argument(
        "--snippy-root",
        type=Path,
        default=WORKFLOW_DATA_ROOT / "snippy",
        help="Root containing the original read-mode Snippy BAMs.",
    )
    parser.add_argument(
        "--stress-root",
        type=Path,
        default=STEP4_DATA_ROOT / "work" / "read_validation" / "detectability_stress",
        help="Working directory for the downsampling stress test.",
    )
    parser.add_argument(
        "--fractions",
        nargs="+",
        type=float,
        default=list(DEFAULT_FRACTIONS),
        help="Read-pair fractions to test.",
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=DEFAULT_REPLICATES,
        help="Replicates per read fraction.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help="Base random seed used to generate reproducible replicate seeds.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=3,
        help="Concurrent sample jobs passed to the Step4 read-validation runner.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=2,
        help="Threads per sample passed to the Step4 read-validation runner.",
    )
    parser.add_argument(
        "--min-support",
        type=int,
        default=3,
        help="Minimum panISa clipped-read support passed to the Step4 runner.",
    )
    parser.add_argument(
        "--batch-label",
        default="detectability_stress",
        help="Work-root label used by Step4 runner and parser.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Remove any existing stress-root outputs before rerunning.",
    )
    parser.add_argument(
        "--skip-tool-check",
        action="store_true",
        help="Skip the Step4 validation-environment checks in the shell runner.",
    )
    parser.add_argument(
        "--step4-runner",
        type=Path,
        default=SCRIPT_DIR / "step4_03e_run_is_read_validation.sh",
        help="Shell runner that executes ISMapper and panISa on the downsampled files.",
    )
    parser.add_argument(
        "--step4-parser",
        type=Path,
        default=SCRIPT_DIR / "step4_03_validate_prn_with_reads.py",
        help="Parser that converts the runner outputs into the read-validation tables.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Merged detectability stress results TSV. Defaults inside the stress root.",
    )
    return parser


def choose_exemplar_rows(subset_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_sample = {
        normalize_text(row.get("sample_id_canonical", "")): row
        for row in subset_rows
        if normalize_text(row.get("sample_id_canonical", ""))
    }
    chosen_rows: list[dict[str, str]] = []
    for target in TARGET_EXEMPLARS:
        selected = None
        for sample_id in target["candidate_sample_ids"]:
            row = by_sample.get(sample_id)
            if row is None:
                continue
            if normalize_text(row.get("raw_read_link_status", "")) != "linked":
                continue
            selected = dict(row)
            selected["target_family_key"] = target["family_key"]
            selected["target_family_label"] = target["family_label"]
            selected["parent_sample_id"] = sample_id
            break
        if selected is None:
            raise ValueError(f"No paired-Illumina exemplar found for {target['family_label']}")
        chosen_rows.append(selected)
    return chosen_rows


def make_fraction_label(fraction: float) -> str:
    if fraction >= 0.999999:
        return "100%"
    return f"{int(round(fraction * 100))}%"


def build_replicate_rows(
    exemplar_rows: list[dict[str, str]],
    *,
    reads_root: Path,
    snippy_root: Path,
    stress_root: Path,
    fractions: list[float],
    replicates: int,
    base_seed: int,
) -> list[dict[str, str]]:
    reads_out = stress_root / "reads_clean"
    snippy_out = stress_root / "snippy"
    reads_out.mkdir(parents=True, exist_ok=True)
    snippy_out.mkdir(parents=True, exist_ok=True)

    replicate_rows: list[dict[str, str]] = []
    for family_index, exemplar in enumerate(exemplar_rows, start=1):
        parent_sample_id = normalize_text(exemplar.get("sample_id_canonical", ""))
        parent_reads_r1 = reads_root / f"{parent_sample_id}_1.fastq.gz"
        parent_reads_r2 = reads_root / f"{parent_sample_id}_2.fastq.gz"
        parent_bam = snippy_root / parent_sample_id / "snps.bam"
        if not parent_reads_r1.exists() or not parent_reads_r2.exists():
            raise FileNotFoundError(f"Paired FASTQ not found for exemplar {parent_sample_id}")
        if not parent_bam.exists():
            raise FileNotFoundError(f"Snippy BAM not found for exemplar {parent_sample_id}")

        for fraction_index, fraction in enumerate(fractions, start=1):
            fraction_label = make_fraction_label(fraction)
            for replicate in range(1, replicates + 1):
                replicate_seed = base_seed + family_index * 10_000 + fraction_index * 100 + replicate
                replicate_sample_id = (
                    f"{parent_sample_id}__{exemplar['target_family_key']}__f{fraction_token(fraction)}__r{replicate:02d}"
                )
                out_reads_r1 = reads_out / f"{replicate_sample_id}_1.fastq.gz"
                out_reads_r2 = reads_out / f"{replicate_sample_id}_2.fastq.gz"
                out_bam = snippy_out / replicate_sample_id / "snps.bam"
                out_bam.parent.mkdir(parents=True, exist_ok=True)

                if fraction >= 0.999999:
                    copy_or_symlink(parent_reads_r1, out_reads_r1)
                    copy_or_symlink(parent_reads_r2, out_reads_r2)
                    copy_or_symlink(parent_bam, out_bam)
                    parent_bai = parent_bam.with_suffix(parent_bam.suffix + ".bai")
                    if parent_bai.exists():
                        copy_or_symlink(parent_bai, out_bam.with_suffix(out_bam.suffix + ".bai"))
                    total_pairs = count_fastq_pairs(parent_reads_r1, parent_reads_r2)
                    retained_pairs = total_pairs
                    with gzip.open(parent_reads_r1, "rt") as parent_fastq:
                        selected_names = {
                            normalize_read_name(title) for title, _, _ in FastqGeneralIterator(parent_fastq)
                        }
                    selected_records = 0
                    with pysam.AlignmentFile(parent_bam, "rb") as bam_handle:
                        for record in bam_handle.fetch(until_eof=True):
                            if normalize_read_name(record.query_name) in selected_names:
                                selected_records += 1
                else:
                    total_pairs, retained_pairs, selected_names = downsample_fastq_pair(
                        parent_reads_r1,
                        parent_reads_r2,
                        out_reads_r1,
                        out_reads_r2,
                        fraction=fraction,
                        seed=replicate_seed,
                    )
                    selected_records = downsample_bam(parent_bam, out_bam, selected_names)

                replicate_rows.append(
                    {
                        **{key: clean_value for key, clean_value in exemplar.items()},
                        "sample_id_canonical": replicate_sample_id,
                        "parent_sample_id": parent_sample_id,
                        "target_family_key": exemplar["target_family_key"],
                        "target_family_label": exemplar["target_family_label"],
                        "downsample_fraction": f"{fraction:.2f}",
                        "downsample_fraction_label": fraction_label,
                        "downsample_replicate": str(replicate),
                        "downsample_seed": str(replicate_seed),
                        "downsample_n_read_pairs_total": str(total_pairs),
                        "downsample_n_read_pairs_retained": str(retained_pairs),
                        "downsample_n_bam_records_retained": str(selected_records),
                        "downsample_reads_1_path": str(out_reads_r1),
                        "downsample_reads_2_path": str(out_reads_r2),
                        "downsample_bam_path": str(out_bam),
                    }
                )
    return replicate_rows


def run_step4_validation(
    *,
    runner: Path,
    parser: Path,
    subset: Path,
    stress_root: Path,
    batch_label: str,
    jobs: int,
    threads: int,
    min_support: int,
    force: bool,
    skip_tool_check: bool,
    out_prefix: Path,
) -> None:
    runner_cmd = [
        str(runner),
        "--batch-label",
        batch_label,
        "--subset",
        str(subset),
        "--reads-root",
        str(stress_root / "reads_clean"),
        "--snippy-root",
        str(stress_root / "snippy"),
        "--jobs",
        str(jobs),
        "--threads",
        str(threads),
        "--min-support",
        str(min_support),
    ]
    if force:
        runner_cmd.append("--force")
    if skip_tool_check:
        runner_cmd.append("--skip-tool-check")

    env = os.environ.copy()
    print("[1/3] Running Step4 read-validation tools on the downsampled exemplars")
    subprocess.run(runner_cmd, check=True, env=env)

    batch_path = stress_root / "bp_prn_read_validation_batch.tsv"
    work_root = stress_root
    parse_cmd = [
        sys.executable,
        str(parser),
        "--subset",
        str(subset),
        "--is-work-root",
        str(work_root),
        "--batch",
        str(batch_path),
        "--batch-label",
        batch_label,
        "--out",
        str(out_prefix.with_name(f"{out_prefix.stem}_read_validation.tsv")),
        "--evidence-out",
        str(out_prefix.with_name(f"{out_prefix.stem}_read_validation_is_calls.tsv")),
        "--tsd-out",
        str(out_prefix.with_name(f"{out_prefix.stem}_read_validation_tsd.tsv")),
    ]
    print("[2/3] Parsing the Step4 outputs into read-validation tables")
    subprocess.run(parse_cmd, check=True, env=env)


def build_results_table(
    manifest_rows: list[dict[str, str]],
    validation_rows: list[dict[str, str]],
    *,
    source_path: Path,
) -> pd.DataFrame:
    manifest = pd.DataFrame(manifest_rows).fillna("")
    validation = pd.DataFrame(validation_rows).fillna("")
    if manifest.empty or validation.empty:
        return pd.DataFrame()

    if "notes" in manifest.columns:
        manifest = manifest.rename(columns={"notes": "source_notes"})

    validation = validation.loc[
        :,
        [
            "sample_id_canonical",
            "read_validation_status",
            "read_support_class",
            "n_supporting_reads",
            "n_contradicting_reads",
            "targeted_locus_assembly_status",
            "validation_method",
            "validator_version",
            "notes",
        ],
    ].copy()
    validation = validation.rename(columns={"notes": "validation_notes"})
    merged = manifest.merge(validation, on="sample_id_canonical", how="left")
    merged["analysis_layer"] = "downsampling"
    status_key = merged["read_validation_status"].map(normalize_text)
    merged["status_bucket"] = status_key.map(
        lambda value: {
            "supported": "recovered",
            "supported_candidate": "recovered",
            "supported_concordant": "recovered",
            "no_prn_is_signal_detected": "true_nonrecovery",
            "tool_output_missing": "compatibility_excluded",
            "unresolved": "compatibility_excluded",
        }.get(value, "compatibility_excluded")
    )
    merged["compatibility_state"] = status_key.map(
        lambda value: {
            "supported": "observed_recovered",
            "supported_candidate": "observed_recovered",
            "supported_concordant": "observed_recovered",
            "no_prn_is_signal_detected": "observed_no_prn_signal",
            "tool_output_missing": "tool_output_missing",
            "unresolved": "unresolved",
        }.get(value, "unresolved_other")
    )
    merged["run_status"] = status_key.map(
        lambda value: {
            "supported": "completed",
            "supported_candidate": "completed",
            "supported_concordant": "completed",
            "no_prn_is_signal_detected": "completed",
            "tool_output_missing": "tool_output_missing",
            "unresolved": "unresolved",
        }.get(value, value or "missing")
    )
    for column in ["validation_method", "validator_version", "source_notes", "validation_notes"]:
        if column in merged.columns:
            merged[column] = merged[column].map(normalize_text)
        else:
            merged[column] = ""
    merged["count_in_total_denominator"] = 1
    merged["count_in_resolved_denominator"] = merged["status_bucket"].isin({"recovered", "true_nonrecovery"}).astype(int)
    merged["count_as_recovered"] = (merged["status_bucket"] == "recovered").astype(int)
    merged["count_as_true_nonrecovery"] = (merged["status_bucket"] == "true_nonrecovery").astype(int)
    merged["count_as_compatibility_excluded"] = (merged["status_bucket"] == "compatibility_excluded").astype(int)
    merged["source_file"] = str(source_path)
    merged["notes"] = merged.apply(
        lambda row: (
            f"parent_sample_id={normalize_text(row.get('parent_sample_id', ''))};"
            f"fraction={normalize_text(row.get('downsample_fraction_label', ''))};"
            f"replicate={normalize_text(row.get('downsample_replicate', ''))};"
            f"seed={normalize_text(row.get('downsample_seed', ''))};"
            f"read_validation_status={normalize_text(row.get('read_validation_status', ''))};"
            f"run_status={normalize_text(row.get('run_status', ''))}"
        ),
        axis=1,
    )
    return merged


def main() -> int:
    args = build_arg_parser().parse_args()
    if args.out is None:
        args.out = args.stress_root / "bp_prn_detectability_stress_results.tsv"

    if args.force and args.stress_root.exists():
        shutil.rmtree(args.stress_root)
    args.stress_root.mkdir(parents=True, exist_ok=True)

    subset_rows = load_tsv_rows(args.subset)
    exemplar_rows = choose_exemplar_rows(subset_rows)
    manifest_rows = build_replicate_rows(
        exemplar_rows,
        reads_root=args.reads_root,
        snippy_root=args.snippy_root,
        stress_root=args.stress_root,
        fractions=list(args.fractions),
        replicates=args.replicates,
        base_seed=args.seed,
    )

    manifest_fieldnames = []
    if manifest_rows:
        for key in manifest_rows[0].keys():
            if key not in manifest_fieldnames:
                manifest_fieldnames.append(key)
        for extra_key in [
            "parent_sample_id",
            "target_family_key",
            "target_family_label",
            "downsample_fraction",
            "downsample_fraction_label",
            "downsample_replicate",
            "downsample_seed",
            "downsample_n_read_pairs_total",
            "downsample_n_read_pairs_retained",
            "downsample_n_bam_records_retained",
            "downsample_reads_1_path",
            "downsample_reads_2_path",
            "downsample_bam_path",
        ]:
            if extra_key not in manifest_fieldnames:
                manifest_fieldnames.append(extra_key)
    subset_path = args.stress_root / "bp_prn_detectability_stress_subset.tsv"
    write_tsv(subset_path, manifest_fieldnames, manifest_rows)

    run_step4_validation(
        runner=args.step4_runner,
        parser=args.step4_parser,
        subset=subset_path,
        stress_root=args.stress_root,
        batch_label=args.batch_label,
        jobs=args.jobs,
        threads=args.threads,
        min_support=args.min_support,
        force=args.force,
        skip_tool_check=args.skip_tool_check,
        out_prefix=args.out,
    )

    validation_path = args.out.with_name(f"{args.out.stem}_read_validation.tsv")
    validation_rows = load_tsv_rows(validation_path)
    results = build_results_table(manifest_rows, validation_rows, source_path=args.out.with_name(f"{args.out.stem}_read_validation.tsv"))
    results.to_csv(args.out, sep="\t", index=False)

    print("[3/3] Wrote downsampling detectability results")
    print(f"  manifest rows: {len(manifest_rows)}")
    print(f"  validation rows: {len(validation_rows)}")
    print(f"  result rows: {len(results)}")
    print(f"  merged results: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
