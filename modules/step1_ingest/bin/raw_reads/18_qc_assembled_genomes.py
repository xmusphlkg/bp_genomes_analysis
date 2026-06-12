#!/usr/bin/env python3
"""Run QUAST/CheckM QC for assembled raw-read genomes and produce a filtered table.

Runtime environment:
  PROJECT_ENV_KEY: bio_tools
  PROJECT_ENV_NAME: pertussis-bio-tools
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
from pathlib import Path
from tempfile import NamedTemporaryFile

from raw_read_utils import project_module_data_root


OUTPUT_COLUMNS = [
    "run_accession",
    "sample_id_canonical",
    "biosample_accession",
    "analysis_cohort_id",
    "priority_tier",
    "priority_reason",
    "source_manifest",
    "run_source",
    "raw_read_link_status",
    "raw_read_run_count",
    "country",
    "year",
    "ena_library_layout",
    "ena_instrument_platform",
    "estimated_total_bytes",
    "assembly_server",
    "assembly_dir",
    "contigs_fasta",
    "latest_status",
    "latest_status_message",
    "latest_started_at",
    "latest_finished_at",
    "contig_count",
    "total_bases",
    "largest_contig",
    "contig_n50",
    "gc_percent",
    "n_bases",
    "quast_status",
    "quast_total_length",
    "quast_contig_count",
    "quast_largest_contig",
    "quast_gc_percent",
    "quast_n50",
    "quast_l50",
    "quast_n_per_100kbp",
    "checkm_status",
    "checkm_marker_lineage",
    "checkm_genome_count",
    "checkm_marker_count",
    "checkm_marker_set_count",
    "checkm_completeness",
    "checkm_contamination",
    "checkm_strain_heterogeneity",
    "qc_decision",
    "qc_reason",
]

QUAST_FIELD_MAP = {
    "Assembly": "assembly_name",
    "Total length": "quast_total_length",
    "# contigs": "quast_contig_count",
    "Largest contig": "quast_largest_contig",
    "GC (%)": "quast_gc_percent",
    "N50": "quast_n50",
    "L50": "quast_l50",
    "# N's per 100 kbp": "quast_n_per_100kbp",
}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[4]


STEP1_DATA_ROOT = project_module_data_root("step1_ingest")


def normalize_text(value: str | None) -> str:
    return (value or "").strip()


def parse_float(value: str | None) -> float | None:
    text = normalize_text(value)
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", newline="", encoding="utf-8", delete=False, dir=path.parent) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
        temp_path = Path(handle.name)
    temp_path.replace(path)


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"no rows found in {path}")
    return rows


def resolve_command(explicit: str | None, candidates: list[str]) -> str | None:
    if explicit:
        return explicit
    for candidate in candidates:
        found = shutil.which(candidate)
        if found:
            return found
    return None


def stage_inputs(rows: list[dict[str, str]], input_dir: Path) -> dict[str, str]:
    input_dir.mkdir(parents=True, exist_ok=True)
    labels: dict[str, str] = {}
    for row in rows:
        run_accession = row["run_accession"]
        source = Path(row["contigs_fasta"]).resolve()
        target = input_dir / f"{run_accession}.fa"
        if target.exists() or target.is_symlink():
            target.unlink()
        target.symlink_to(source)
        labels[f"{run_accession}.fa"] = run_accession
        labels[run_accession] = run_accession
    return labels


def run_command(command: list[str], *, cwd: Path | None = None) -> None:
    subprocess.run(command, cwd=str(cwd) if cwd else None, check=True)


def parse_quast_transposed(path: Path, labels: dict[str, str]) -> dict[str, dict[str, str]]:
    metrics: dict[str, dict[str, str]] = {}
    if not path.exists():
        return metrics
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            raw_name = normalize_text(row.get("Assembly", ""))
            if not raw_name:
                continue
            run_accession = labels.get(raw_name, Path(raw_name).stem)
            metrics.setdefault(run_accession, {})
            for source_name, target_name in QUAST_FIELD_MAP.items():
                if source_name == "Assembly":
                    continue
                metrics[run_accession][target_name] = normalize_text(row.get(source_name, ""))
    return metrics


def run_quast(
    rows: list[dict[str, str]],
    *,
    quast_cmd: str,
    workdir: Path,
    threads: int,
    batch_size: int,
) -> dict[str, dict[str, str]]:
    input_dir = workdir / "inputs"
    labels = stage_inputs(rows, input_dir)
    metrics: dict[str, dict[str, str]] = {}

    staged_paths = [input_dir / f"{row['run_accession']}.fa" for row in rows]
    for start in range(0, len(staged_paths), batch_size):
        chunk = staged_paths[start : start + batch_size]
        chunk_dir = workdir / "quast" / f"chunk_{start // batch_size:04d}"
        chunk_dir.mkdir(parents=True, exist_ok=True)
        command = [
            quast_cmd,
            "--threads",
            str(threads),
            "--output-dir",
            str(chunk_dir),
            "--min-contig",
            "500",
        ] + [str(path) for path in chunk]
        run_command(command)
        metrics.update(parse_quast_transposed(chunk_dir / "transposed_report.tsv", labels))
    return metrics


def parse_checkm_table(path: Path) -> dict[str, dict[str, str]]:
    metrics: dict[str, dict[str, str]] = {}
    if not path.exists():
        return metrics
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            run_accession = normalize_text(row.get("Bin Id", "")) or normalize_text(row.get("Bin", ""))
            if not run_accession:
                continue
            if run_accession.endswith(".fa"):
                run_accession = Path(run_accession).stem
            metrics[run_accession] = {
                "checkm_marker_lineage": normalize_text(row.get("Marker lineage", "")),
                "checkm_genome_count": normalize_text(row.get("# genomes", "")),
                "checkm_marker_count": normalize_text(row.get("# markers", "")),
                "checkm_marker_set_count": normalize_text(row.get("# marker sets", "")),
                "checkm_completeness": normalize_text(row.get("Completeness", "")),
                "checkm_contamination": normalize_text(row.get("Contamination", "")),
                "checkm_strain_heterogeneity": normalize_text(row.get("Strain heterogeneity", "")),
            }
    return metrics


def run_checkm(
    rows: list[dict[str, str]],
    *,
    checkm_cmd: str,
    workdir: Path,
    threads: int,
) -> dict[str, dict[str, str]]:
    input_dir = workdir / "inputs"
    stage_inputs(rows, input_dir)
    out_dir = workdir / "checkm"
    out_dir.mkdir(parents=True, exist_ok=True)
    lineage_ms = out_dir / "lineage.ms"
    qa_tsv = out_dir / "qa.tsv"

    run_command(
        [
            checkm_cmd,
            "lineage_wf",
            "-t",
            str(threads),
            "-x",
            "fa",
            str(input_dir),
            str(out_dir),
        ]
    )
    run_command(
        [
            checkm_cmd,
            "qa",
            "-t",
            str(threads),
            "-o",
            "2",
            "--tab_table",
            "-f",
            str(qa_tsv),
            str(lineage_ms),
            str(out_dir),
        ]
    )
    return parse_checkm_table(qa_tsv)


def build_output_rows(
    rows: list[dict[str, str]],
    *,
    quast_metrics: dict[str, dict[str, str]],
    checkm_metrics: dict[str, dict[str, str]],
    quast_ran: bool,
    checkm_ran: bool,
    min_completeness: float,
    max_contamination: float,
) -> list[dict[str, str]]:
    output_rows: list[dict[str, str]] = []
    for row in rows:
        run_accession = row["run_accession"]
        enriched = {column: normalize_text(row.get(column, "")) for column in OUTPUT_COLUMNS}
        quast = quast_metrics.get(run_accession, {})
        checkm = checkm_metrics.get(run_accession, {})

        enriched["quast_status"] = "completed" if quast_ran and quast else ("not_run" if not quast_ran else "missing")
        enriched["checkm_status"] = "completed" if checkm_ran and checkm else ("not_run" if not checkm_ran else "missing")

        for key, value in quast.items():
            enriched[key] = value
        for key, value in checkm.items():
            enriched[key] = value

        completeness = parse_float(enriched.get("checkm_completeness", ""))
        contamination = parse_float(enriched.get("checkm_contamination", ""))
        if completeness is None or contamination is None:
            enriched["qc_decision"] = "pending_checkm"
            enriched["qc_reason"] = "checkm_metrics_missing"
        elif completeness < min_completeness:
            enriched["qc_decision"] = "fail"
            enriched["qc_reason"] = f"completeness_lt_{min_completeness:g}"
        elif contamination > max_contamination:
            enriched["qc_decision"] = "fail"
            enriched["qc_reason"] = f"contamination_gt_{max_contamination:g}"
        else:
            enriched["qc_decision"] = "pass"
            enriched["qc_reason"] = "checkm_thresholds_passed"
        output_rows.append(enriched)
    output_rows.sort(key=lambda item: (item["qc_decision"] != "pass", item["assembly_server"], item["run_accession"]))
    return output_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run QUAST and CheckM QC for assembled raw-read genomes."
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_raw_read_assembly_manifest.tsv",
        help="Input manifest TSV produced by 17_collect_assembled_genomes.py.",
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "raw_read_qc_work",
        help="Work directory for staged FASTA, QUAST, and CheckM outputs.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_raw_read_assembly_qc.tsv",
        help="Output QC summary TSV.",
    )
    parser.add_argument(
        "--passed-output",
        type=Path,
        default=STEP1_DATA_ROOT / "outputs" / "bp_raw_read_assembly_qc_pass.tsv",
        help="Output TSV containing only QC-passed assemblies.",
    )
    parser.add_argument("--threads", type=int, default=12, help="Threads for QUAST/CheckM.")
    parser.add_argument("--quast-batch-size", type=int, default=400, help="Assemblies per QUAST batch.")
    parser.add_argument("--min-completeness", type=float, default=95.0, help="Minimum CheckM completeness.")
    parser.add_argument("--max-contamination", type=float, default=5.0, help="Maximum CheckM contamination.")
    parser.add_argument("--quast-cmd", default=None, help="Explicit QUAST executable path.")
    parser.add_argument("--checkm-cmd", default=None, help="Explicit CheckM executable path.")
    parser.add_argument("--skip-quast", action="store_true", help="Do not run QUAST even if available.")
    parser.add_argument("--skip-checkm", action="store_true", help="Do not run CheckM even if available.")
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    rows = load_manifest(args.manifest)
    args.workdir.mkdir(parents=True, exist_ok=True)

    quast_cmd = None if args.skip_quast else resolve_command(args.quast_cmd, ["quast.py", "quast"])
    checkm_cmd = None if args.skip_checkm else resolve_command(args.checkm_cmd, ["checkm"])

    quast_metrics: dict[str, dict[str, str]] = {}
    checkm_metrics: dict[str, dict[str, str]] = {}

    if quast_cmd:
        quast_metrics = run_quast(
            rows,
            quast_cmd=quast_cmd,
            workdir=args.workdir,
            threads=args.threads,
            batch_size=args.quast_batch_size,
        )
    if checkm_cmd:
        checkm_metrics = run_checkm(
            rows,
            checkm_cmd=checkm_cmd,
            workdir=args.workdir,
            threads=args.threads,
        )

    output_rows = build_output_rows(
        rows,
        quast_metrics=quast_metrics,
        checkm_metrics=checkm_metrics,
        quast_ran=bool(quast_cmd),
        checkm_ran=bool(checkm_cmd),
        min_completeness=args.min_completeness,
        max_contamination=args.max_contamination,
    )
    passed_rows = [row for row in output_rows if row["qc_decision"] == "pass"]

    write_tsv(args.output, output_rows, OUTPUT_COLUMNS)
    write_tsv(args.passed_output, passed_rows, OUTPUT_COLUMNS)

    print(f"Wrote QC summary: {args.output}")
    print(f"Wrote QC-passed assemblies: {args.passed_output}")
    print(f"Assemblies evaluated: {len(output_rows)}")
    print(f"QC passed: {len(passed_rows)}")
    if not quast_cmd:
        print("QUAST not run: command not found or --skip-quast specified")
    if not checkm_cmd:
        print("CheckM not run: command not found or --skip-checkm specified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
