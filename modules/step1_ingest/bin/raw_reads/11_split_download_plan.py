#!/usr/bin/env python3
"""Split a run-level download plan into deterministic shards for multiple servers."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from raw_read_utils import project_module_data_root


def repo_root() -> Path:
    return Path(__file__).resolve().parents[4]


STEP4_DATA_ROOT = project_module_data_root("step4_prn_validation")


def parse_servers(raw: str) -> list[str]:
    servers = [item.strip() for item in raw.split(",") if item.strip()]
    if not servers:
        raise ValueError("at least one server name must be provided")
    return servers


def load_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    return fieldnames, rows


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_int(value: str | None) -> int:
    text = (value or "").strip()
    if not text:
        return 0
    try:
        return int(text)
    except ValueError:
        return 0


def row_weight(row: dict[str, str]) -> int:
    estimated_total_bytes = parse_int(row.get("estimated_total_bytes", ""))
    if estimated_total_bytes > 0:
        return estimated_total_bytes
    fastq_bytes = sum(
        parse_int(part)
        for part in (row.get("ena_fastq_bytes", "") or "").split(";")
        if part.strip()
    )
    if fastq_bytes > 0:
        return fastq_bytes
    return 1


def write_run_list(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(f"{row['run_accession']}\n")


def write_server_launcher(path: Path, server_name: str, plan_tsv_path: Path, run_list_path: Path) -> None:
    content = "\n".join(
        [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            "",
            'ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"',
            '# shellcheck disable=SC1091',
            'source "${ROOT}/workflow/lib/runtime_envs.sh"',
            'project_env_load_config "$ROOT"',
            f'SERVER_NAME="{server_name}"',
            'PLAN_TSV="${PLAN_TSV:-' + str(plan_tsv_path) + '}"',
            'RUN_LIST="${RUN_LIST:-' + str(run_list_path) + '}"',
            'DATA_ROOT="${DATA_ROOT:-$(project_module_data_root step4_prn_validation)}"',
            'WORKDIR="${WORKDIR:-${DATA_ROOT}/work/${SERVER_NAME}}"',
            'OUTDIR="${OUTDIR:-${DATA_ROOT}/outputs/assemblies/${SERVER_NAME}}"',
            'THREADS="${THREADS:-12}"',
            'JOBS="${JOBS:-2}"',
            'MAX_RUNS="${MAX_RUNS:-}"',
            'CONDA_ENV="${CONDA_ENV:-}"',
            "",
            'CMD=("bash" "${ROOT}/modules/step1_ingest/bin/raw_reads/12_run_shard_to_assembly.sh"',
            '  --plan-tsv "${PLAN_TSV}"',
            '  --run-list "${RUN_LIST}"',
            '  --workdir "${WORKDIR}"',
            '  --outdir "${OUTDIR}"',
            '  --threads "${THREADS}"',
            '  --jobs "${JOBS}")',
            "",
            'if [[ -n "${MAX_RUNS}" ]]; then',
            '  CMD+=(--max-runs "${MAX_RUNS}")',
            "fi",
            "if [[ -n \"${CONDA_ENV:-}\" ]]; then",
            "  CMD+=(--conda-env \"${CONDA_ENV}\")",
            "fi",
            "",
            'echo "[Info] server=${SERVER_NAME} plan_tsv=${PLAN_TSV} run_list=${RUN_LIST}"',
            'echo "[Info] workdir=${WORKDIR} outdir=${OUTDIR} threads=${THREADS} jobs=${JOBS}"',
            '"${CMD[@]}"',
            "",
        ]
    )
    path.write_text(content, encoding="utf-8")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Split run-level raw-read plan across multiple servers."
    )
    parser.add_argument(
        "--plan",
        type=Path,
        default=STEP4_DATA_ROOT / "inputs" / "bp_raw_reads_download_plan.tsv",
        help="Run-level download plan TSV.",
    )
    parser.add_argument(
        "--servers",
        default="server1,server2,server3",
        help="Comma-separated server names.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=STEP4_DATA_ROOT / "inputs" / "shards",
        help="Output directory for plan shards and run lists.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    servers = parse_servers(args.servers)

    input_fields, rows = load_rows(args.plan)
    if not rows:
        raise ValueError(f"no plan rows found in {args.plan}")

    out_fields = list(input_fields)
    for extra in ["assigned_server", "assigned_rank"]:
        if extra not in out_fields:
            out_fields.append(extra)

    buckets: dict[str, list[dict[str, str]]] = {server: [] for server in servers}
    server_totals: dict[str, int] = {server: 0 for server in servers}

    ranked_rows = []
    for row in rows:
        run_accession = (row.get("run_accession") or "").strip()
        if not run_accession:
            continue
        ranked_rows.append((row_weight(row), row))

    ranked_rows.sort(
        key=lambda item: (
            parse_int(item[1].get("priority_tier", "")) or 999999,
            -item[0],
            item[1].get("plan_row_id", ""),
            item[1].get("run_accession", ""),
        )
    )

    for weight, row in ranked_rows:
        server_name = min(
            servers,
            key=lambda server: (
                server_totals[server],
                len(buckets[server]),
                server,
            ),
        )
        row_copy = dict(row)
        row_copy["assigned_server"] = server_name
        row_copy["assigned_rank"] = ""
        buckets[server_name].append(row_copy)
        server_totals[server_name] += weight

    args.outdir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, str]] = []
    for server_name in servers:
        server_rows = buckets[server_name]
        server_rows.sort(key=lambda row: row.get("plan_row_id", ""))
        for index, row in enumerate(server_rows, start=1):
            row["assigned_rank"] = str(index)

        shard_tsv = args.outdir / f"bp_raw_reads_download_plan.{server_name}.tsv"
        run_list = args.outdir / f"bp_raw_reads_runs.{server_name}.txt"
        launcher = args.outdir / f"run_{server_name}.sh"

        write_tsv(shard_tsv, out_fields, server_rows)
        write_run_list(run_list, server_rows)
        write_server_launcher(launcher, server_name, shard_tsv, run_list)

        summary_rows.append(
            {
                "server_name": server_name,
                "n_runs": str(len(server_rows)),
                "estimated_total_bytes": str(server_totals[server_name]),
                "plan_shard": str(shard_tsv),
                "run_list": str(run_list),
                "launcher": str(launcher),
            }
        )

    summary_path = args.outdir / "bp_raw_reads_shard_summary.tsv"
    write_tsv(
        summary_path,
        ["server_name", "n_runs", "estimated_total_bytes", "plan_shard", "run_list", "launcher"],
        summary_rows,
    )

    print(f"Wrote shard summary: {summary_path}")
    for row in summary_rows:
        print(f"{row['server_name']}: {row['n_runs']} runs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
