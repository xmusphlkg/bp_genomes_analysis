#!/usr/bin/env python3
"""Rebuild distributed raw-read shards from the unfinished subset of an existing shard set."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from datetime import datetime
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.lib.project_paths import project_module_data_root


DONE_STATUSES = {"assembled", "already_done", "skipped"}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def parse_servers(raw: str) -> list[str]:
    servers = [item.strip() for item in raw.split(",") if item.strip()]
    if not servers:
        raise ValueError("at least one server name must be provided")
    return servers


def parse_statuses(raw: str) -> set[str]:
    values = {item.strip() for item in raw.split(",") if item.strip()}
    if not values:
        raise ValueError("at least one done status must be provided")
    return values


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


def load_plan_rows(shard_dir: Path, servers: list[str]) -> tuple[list[str], list[dict[str, str]]]:
    fieldnames: list[str] = []
    rows: list[dict[str, str]] = []
    seen_runs: set[str] = set()
    for server in servers:
      plan_path = shard_dir / f"bp_raw_reads_download_plan.{server}.tsv"
      if not plan_path.exists():
          raise FileNotFoundError(f"missing plan shard: {plan_path}")
      with plan_path.open(newline="", encoding="utf-8") as handle:
          reader = csv.DictReader(handle, delimiter="\t")
          if not fieldnames:
              fieldnames = list(reader.fieldnames or [])
          for row in reader:
              run_accession = (row.get("run_accession") or "").strip()
              if not run_accession or run_accession in seen_runs:
                  continue
              row_copy = dict(row)
              row_copy["source_server"] = server
              rows.append(row_copy)
              seen_runs.add(run_accession)
    if not rows:
        raise ValueError(f"no plan rows found under {shard_dir}")
    return fieldnames, rows


def parse_timestamp(value: str | None) -> float:
    text = (value or "").strip()
    if not text:
        return float("-inf")
    try:
        return datetime.fromisoformat(text).timestamp()
    except ValueError:
        return float("-inf")


def load_status_records(assemblies_root: Path, servers: list[str]) -> dict[str, dict[str, object]]:
    latest: dict[str, dict[str, object]] = {}
    for server in servers:
        status_path = assemblies_root / server / "run_status.tsv"
        if not status_path.exists():
            continue
        with status_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for index, row in enumerate(reader):
                run_accession = (row.get("run_accession") or "").strip()
                if not run_accession:
                    continue
                status = (row.get("status") or "").strip()
                sort_key = (
                    parse_timestamp(row.get("finished_at")),
                    parse_timestamp(row.get("started_at")),
                    index,
                    server,
                )
                record = latest.setdefault(
                    run_accession,
                    {
                        "latest_status": "",
                        "latest_status_server": "",
                        "observed_statuses": set(),
                        "_sort_key": (float("-inf"), float("-inf"), -1, ""),
                    },
                )
                observed_statuses = record["observed_statuses"]
                assert isinstance(observed_statuses, set)
                observed_statuses.add(status)
                if sort_key >= record["_sort_key"]:
                    record["latest_status"] = status
                    record["latest_status_server"] = server
                    record["_sort_key"] = sort_key
    return latest


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


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
            'ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"',
            'source "${ROOT}/workflow/lib/runtime_envs.sh"',
            'project_env_load_config "${ROOT}"',
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
            'if [[ -n "${CONDA_ENV:-}" ]]; then',
            '  CMD+=(--conda-env "${CONDA_ENV}")',
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
        description="Rebalance unfinished distributed raw-read runs into a new shard directory."
    )
    parser.add_argument(
        "--input-shard-dir",
        type=Path,
        default=project_module_data_root("step4_prn_validation") / "inputs" / "external_gapfill_shards",
        help="Existing shard directory to rebalance from.",
    )
    parser.add_argument(
        "--output-shard-dir",
        type=Path,
        required=True,
        help="New shard directory for the rebalanced unfinished subset.",
    )
    parser.add_argument(
        "--assemblies-root",
        type=Path,
        default=project_module_data_root("step4_prn_validation") / "outputs" / "assemblies",
        help="Root directory containing per-server run_status.tsv files.",
    )
    parser.add_argument(
        "--servers",
        default="server1,server2,server3",
        help="Comma-separated server names for both input and output assignment.",
    )
    parser.add_argument(
        "--input-servers",
        default="",
        help="Comma-separated input server names. Defaults to --servers when omitted.",
    )
    parser.add_argument(
        "--output-servers",
        default="",
        help="Comma-separated output server names. Defaults to --servers when omitted.",
    )
    parser.add_argument(
        "--done-statuses",
        default="assembled,already_done,skipped",
        help="Comma-separated statuses to treat as finished.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    servers = parse_servers(args.servers)
    input_servers = parse_servers(args.input_servers) if args.input_servers else list(servers)
    output_servers = parse_servers(args.output_servers) if args.output_servers else list(servers)
    done_statuses = parse_statuses(args.done_statuses)

    input_fields, rows = load_plan_rows(args.input_shard_dir, input_servers)
    status_records = load_status_records(args.assemblies_root, input_servers)

    out_fields = list(input_fields)
    for extra in [
        "source_server",
        "latest_status",
        "latest_status_server",
        "observed_statuses",
        "assigned_server",
        "assigned_rank",
    ]:
        if extra not in out_fields:
            out_fields.append(extra)

    pending_rows: list[dict[str, str]] = []
    status_counter: dict[str, int] = defaultdict(int)
    for row in rows:
        run_accession = (row.get("run_accession") or "").strip()
        record = status_records.get(run_accession)
        latest_status = ""
        latest_status_server = ""
        observed_statuses: set[str] = set()
        if record:
            latest_status = str(record.get("latest_status", "") or "")
            latest_status_server = str(record.get("latest_status_server", "") or "")
            observed = record.get("observed_statuses", set())
            if isinstance(observed, set):
                observed_statuses = {str(item) for item in observed if str(item)}
        status_counter[latest_status or "missing"] += 1
        if observed_statuses & done_statuses:
            continue
        row_copy = dict(row)
        row_copy["latest_status"] = latest_status
        row_copy["latest_status_server"] = latest_status_server
        row_copy["observed_statuses"] = ",".join(sorted(observed_statuses))
        pending_rows.append(row_copy)

    if not pending_rows:
        raise SystemExit("No unfinished runs remain; nothing to rebalance.")

    buckets: dict[str, list[dict[str, str]]] = {server: [] for server in output_servers}
    server_totals: dict[str, int] = {server: 0 for server in output_servers}

    ranked_rows = [(row_weight(row), row) for row in pending_rows]
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
            output_servers,
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

    args.output_shard_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, str]] = []
    for server_name in output_servers:
        server_rows = buckets[server_name]
        server_rows.sort(key=lambda row: row.get("plan_row_id", ""))
        for index, row in enumerate(server_rows, start=1):
            row["assigned_rank"] = str(index)

        shard_tsv = args.output_shard_dir / f"bp_raw_reads_download_plan.{server_name}.tsv"
        run_list = args.output_shard_dir / f"bp_raw_reads_runs.{server_name}.txt"
        launcher = args.output_shard_dir / f"run_{server_name}.sh"

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

    summary_path = args.output_shard_dir / "bp_raw_reads_shard_summary.tsv"
    write_tsv(
        summary_path,
        ["server_name", "n_runs", "estimated_total_bytes", "plan_shard", "run_list", "launcher"],
        summary_rows,
    )

    print(f"Input shard dir: {args.input_shard_dir}")
    print(f"Assemblies root: {args.assemblies_root}")
    print(f"Input servers: {','.join(input_servers)}")
    print(f"Output servers: {','.join(output_servers)}")
    print(f"Done statuses: {','.join(sorted(done_statuses))}")
    print(f"Total original runs: {len(rows)}")
    print(f"Pending runs: {len(pending_rows)}")
    for status, count in sorted(status_counter.items()):
        print(f"status[{status}]={count}")
    print(f"Wrote shard summary: {summary_path}")
    for row in summary_rows:
        print(f"{row['server_name']}: {row['n_runs']} runs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
