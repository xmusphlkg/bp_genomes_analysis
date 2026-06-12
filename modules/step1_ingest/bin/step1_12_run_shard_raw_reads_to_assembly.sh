#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
exec "${ROOT}/modules/step1_ingest/bin/raw_reads/12_run_shard_to_assembly.sh" "$@"
