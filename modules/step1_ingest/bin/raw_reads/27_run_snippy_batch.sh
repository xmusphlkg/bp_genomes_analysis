#!/usr/bin/env bash
# T11 bootstrap: run Snippy in contig mode using the canonical M3 sample plan.
# Runtime environment:
#   PROJECT_ENV_KEY: bio_tools
#   PROJECT_ENV_NAME: pertussis-bio-tools

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/workflow/lib/runtime_envs.sh"
project_env_load_config "$ROOT"
project_env_prepend_path bio_tools
WORKFLOW_DATA_ROOT="$(project_workflow_root)"
PLAN="${WORKFLOW_DATA_ROOT}/snippy_ctg/snippy_ctg_plan.tsv"
REFERENCE="${ROOT}/pertussis_data/bp_genomes_qc/reference/tohama_i.fasta"
OUTDIR="${WORKFLOW_DATA_ROOT}/snippy_ctg"
BATCH_DIR="${OUTDIR}/batches"
LOGDIR="$(project_logs_root)/pipeline/snippy_ctg"
SUMMARY_TSV="${OUTDIR}/snippy_ctg_run_status.tsv"
SUMMARY_TXT="${OUTDIR}/snippy_ctg_run_summary.txt"
HISTORY_TSV="${OUTDIR}/snippy_ctg_run_history.tsv"
CPUS=0
JOBS=0
DRY_RUN=false
LIMIT=0
OFFSET=0
BATCH_LABEL=""
USE_DOCKER=false
DOCKER_IMAGE="${SNIPPY_DOCKER_IMAGE:-quay.io/biocontainers/snippy:4.6.0--hdfd78af_6}"
CPUS_SET=false
JOBS_SET=false
AUTO_PARALLEL_MODE="pending"
AUTO_CPUS_PER_SAMPLE="${SNIPPY_AUTO_CPUS_PER_SAMPLE:-2}"
AUTO_MEM_PER_SAMPLE_GB="${SNIPPY_AUTO_MEM_PER_SAMPLE_GB:-3}"
AUTO_RESERVE_CORES="${SNIPPY_AUTO_RESERVE_CORES:-4}"
AUTO_RESERVE_MEM_GB="${SNIPPY_AUTO_RESERVE_MEM_GB:-16}"
AUTO_MAX_JOBS="${SNIPPY_AUTO_MAX_JOBS:-0}"
AUTO_MAX_JOBS_DOCKER="${SNIPPY_AUTO_MAX_JOBS_DOCKER:-48}"
RESOLVE_PULL_DOCKER=true

usage() {
    cat <<'EOF'
Usage:
  bash modules/step1_ingest/bin/raw_reads/27_run_snippy_batch.sh \
        [--plan PATH] [--cpus N] [--jobs N] [--offset N] [--limit N] [--batch-label LABEL] [--dry-run]

Options:
  --plan    Path to snippy_ctg_plan.tsv. Default: workflow/snippy_ctg/snippy_ctg_plan.tsv
  --cpus    CPUs per Snippy sample job. Default: auto (2 unless overridden)
  --jobs    Parallel sample jobs. Default: auto (derived from host cores/RAM)
    --offset  Skip the first N included rows from the plan before selecting work.
  --limit   Restrict to the first N included rows from the plan.
    --batch-label  Optional stable label for batch audit files.
  --dry-run Report planned work without invoking Snippy.
EOF
}

detect_host_cores() {
    if command -v nproc >/dev/null 2>&1; then
        nproc
        return 0
    fi

    getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1
}

detect_host_mem_gb() {
    awk '/MemTotal/ {printf "%d\n", ($2 / 1024 / 1024)}' /proc/meminfo 2>/dev/null || echo 1
}

auto_tune_parallelism() {
    local selected_total="$1"
    local host_cores host_mem_gb usable_cores usable_mem_gb jobs_by_core jobs_by_mem jobs_cap

    if [[ "$CPUS_SET" == false ]]; then
        CPUS="$AUTO_CPUS_PER_SAMPLE"
    fi

    if [[ "$CPUS" -lt 1 ]]; then
        echo "ERROR: Computed --cpus must be >= 1" >&2
        exit 1
    fi

    if [[ "$JOBS_SET" == true ]]; then
        AUTO_PARALLEL_MODE="manual"
        return 0
    fi

    host_cores="$(detect_host_cores)"
    host_mem_gb="$(detect_host_mem_gb)"
    usable_cores=$((host_cores - AUTO_RESERVE_CORES))
    usable_mem_gb=$((host_mem_gb - AUTO_RESERVE_MEM_GB))

    if [[ "$usable_cores" -lt "$CPUS" ]]; then
        usable_cores="$CPUS"
    fi
    if [[ "$usable_mem_gb" -lt "$AUTO_MEM_PER_SAMPLE_GB" ]]; then
        usable_mem_gb="$AUTO_MEM_PER_SAMPLE_GB"
    fi

    jobs_by_core=$((usable_cores / CPUS))
    jobs_by_mem=$((usable_mem_gb / AUTO_MEM_PER_SAMPLE_GB))
    if [[ "$jobs_by_core" -lt 1 ]]; then
        jobs_by_core=1
    fi
    if [[ "$jobs_by_mem" -lt 1 ]]; then
        jobs_by_mem=1
    fi

    JOBS="$jobs_by_core"
    if [[ "$jobs_by_mem" -lt "$JOBS" ]]; then
        JOBS="$jobs_by_mem"
    fi

    if [[ "$AUTO_MAX_JOBS" -gt 0 && "$JOBS" -gt "$AUTO_MAX_JOBS" ]]; then
        JOBS="$AUTO_MAX_JOBS"
    fi

    if [[ "$USE_DOCKER" == true && "$AUTO_MAX_JOBS_DOCKER" -gt 0 && "$JOBS" -gt "$AUTO_MAX_JOBS_DOCKER" ]]; then
        JOBS="$AUTO_MAX_JOBS_DOCKER"
    fi

    if [[ "$selected_total" -gt 0 && "$JOBS" -gt "$selected_total" ]]; then
        JOBS="$selected_total"
    fi

    jobs_cap="none"
    if [[ "$AUTO_MAX_JOBS" -gt 0 ]]; then
        jobs_cap="$AUTO_MAX_JOBS"
    fi
    if [[ "$USE_DOCKER" == true && "$AUTO_MAX_JOBS_DOCKER" -gt 0 ]]; then
        jobs_cap="docker_cap=${AUTO_MAX_JOBS_DOCKER}${jobs_cap:+,global_cap=${jobs_cap}}"
    fi

    AUTO_PARALLEL_MODE="auto(host_cores=${host_cores},host_mem_gb=${host_mem_gb},reserve_cores=${AUTO_RESERVE_CORES},reserve_mem_gb=${AUTO_RESERVE_MEM_GB},cpus_per_sample=${CPUS},mem_per_sample_gb=${AUTO_MEM_PER_SAMPLE_GB},jobs_cap=${jobs_cap})"
}

resolve_snippy() {
    if command -v snippy >/dev/null 2>&1; then
        USE_DOCKER=false
        return 0
    fi

    if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
        if ! docker image inspect "$DOCKER_IMAGE" >/dev/null 2>&1; then
            if [[ "$RESOLVE_PULL_DOCKER" == true ]]; then
                echo "Pulling Docker image: $DOCKER_IMAGE"
                docker pull "$DOCKER_IMAGE" >/dev/null
            fi
        fi
        USE_DOCKER=true
        return 0
    fi

    echo "ERROR: Snippy is not available in PATH and Docker is unavailable." >&2
    echo "Configured bio_tools env: $(project_env_prefix bio_tools) ($(project_env_name bio_tools))" >&2
    echo "Run: bash workflow/bin/bootstrap_runtime_envs.sh --check" >&2
    echo "Or install Docker and allow pulling image: ${DOCKER_IMAGE}" >&2
    return 1
}

run_snippy() {
    if [[ "$USE_DOCKER" == true ]]; then
        local -a docker_args
        local pertussis_data_real=""
        docker_args=(
            --rm
            --user "$(id -u):$(id -g)"
            -e HOME=/tmp
            -v "${ROOT}:${ROOT}"
            -w "${ROOT}"
        )
        if [[ -L "${ROOT}/pertussis_data" ]]; then
            pertussis_data_real="$(readlink -f "${ROOT}/pertussis_data")"
            if [[ -n "$pertussis_data_real" && -e "$pertussis_data_real" && "$pertussis_data_real" != "${ROOT}/pertussis_data" ]]; then
                docker_args+=( -v "${pertussis_data_real}:${pertussis_data_real}" )
            fi
        fi
        docker run "${docker_args[@]}" \
            "$DOCKER_IMAGE" \
            snippy "$@"
    else
        snippy "$@"
    fi
}

run_one() {
    local acc="$1"
    local fasta="$2"
    local sample_outdir="${OUTDIR}/${acc}"
    local log_path="${LOGDIR}/${acc}.log"
    local status=""
    local run_rc=0

    if [[ ! -f "$fasta" ]]; then
        status="missing_fasta"
    elif [[ -f "${sample_outdir}/snps.aligned.fa" ]]; then
        status="skipped_existing"
    else
        rm -rf "$sample_outdir"
        mkdir -p "$sample_outdir"
        if run_snippy \
            --ctgs "$fasta" \
            --ref "$REFERENCE" \
            --outdir "$sample_outdir" \
            --cpus "$CPUS" \
            --force \
            --quiet >"$log_path" 2>&1; then
            run_rc=0
        else
            run_rc=$?
        fi

        if [[ -f "${sample_outdir}/snps.aligned.fa" ]]; then
            status="ok"
        else
            status="failed"
            if [[ "$run_rc" -ne 0 ]]; then
                printf '\n[wrapper] snippy exit code: %s\n' "$run_rc" >>"$log_path"
            fi
        fi
    fi

    {
        flock -x 200
        printf "%s\t%s\t%s\t%s\t%s\n" \
            "$acc" "$status" "$sample_outdir" "$log_path" "$(date -Iseconds)" >>"$SUMMARY_TSV"
        printf "%s\t%s\t%s\t%s\t%s\n" \
            "$acc" "$status" "$sample_outdir" "$log_path" "$(date -Iseconds)" >>"$BATCH_STATUS_TSV"
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$BATCH_LABEL" "$acc" "$status" "$sample_outdir" "$log_path" "$(date -Iseconds)" >>"$HISTORY_TSV"
    } 200>"$LOCKFILE"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --plan) PLAN="$2"; shift 2 ;;
        --cpus) CPUS="$2"; CPUS_SET=true; shift 2 ;;
        --jobs) JOBS="$2"; JOBS_SET=true; shift 2 ;;
        --offset) OFFSET="$2"; shift 2 ;;
        --limit) LIMIT="$2"; shift 2 ;;
        --batch-label) BATCH_LABEL="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        --help|-h) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ ! -f "$PLAN" ]]; then
    echo "ERROR: Snippy plan not found: $PLAN" >&2
    echo "Build it first with: python modules/step1_ingest/bin/raw_reads/28_build_snippy_ctg_plan.py" >&2
    exit 1
fi

if [[ ! -f "$REFERENCE" ]]; then
    echo "ERROR: Reference genome not found: $REFERENCE" >&2
    exit 1
fi

if [[ "$JOBS_SET" == true && "$JOBS" -lt 1 ]]; then
    echo "ERROR: --jobs must be >= 1" >&2
    exit 1
fi

if [[ "$CPUS_SET" == true && "$CPUS" -lt 1 ]]; then
    echo "ERROR: --cpus must be >= 1" >&2
    exit 1
fi

mkdir -p "$OUTDIR" "$LOGDIR" "$BATCH_DIR"

TASKS_FILE="$(mktemp)"
LOCKFILE="$(mktemp)"
cleanup() {
    rm -f "$TASKS_FILE" "$LOCKFILE"
}
trap cleanup EXIT

awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        idx[$i] = i
    }
    next
}
$idx["include_in_snippy_ctg"] == "True" {
    print $idx["assembly_accession"] "\t" $idx["fasta_path"]
}
' "$PLAN" >"$TASKS_FILE"

INCLUDED_TOTAL="$(wc -l < "$TASKS_FILE")"

if [[ "$OFFSET" -gt 0 ]]; then
    awk -v offset="$OFFSET" 'NR > offset' "$TASKS_FILE" >"${TASKS_FILE}.offset"
    mv "${TASKS_FILE}.offset" "$TASKS_FILE"
fi

if [[ "$LIMIT" -gt 0 ]]; then
    head -n "$LIMIT" "$TASKS_FILE" >"${TASKS_FILE}.limited"
    mv "${TASKS_FILE}.limited" "$TASKS_FILE"
fi

TOTAL="$(wc -l < "$TASKS_FILE")"
if [[ "$TOTAL" -eq 0 ]]; then
    echo "No rows selected from plan: $PLAN"
    exit 0
fi

if [[ -z "$BATCH_LABEL" ]]; then
    if [[ "$LIMIT" -gt 0 ]]; then
        BATCH_LABEL="offset${OFFSET}_limit${LIMIT}"
    else
        BATCH_LABEL="offset${OFFSET}_all"
    fi
fi

BATCH_STATUS_TSV="${BATCH_DIR}/${BATCH_LABEL}.run_status.tsv"
BATCH_SUMMARY_TXT="${BATCH_DIR}/${BATCH_LABEL}.run_summary.txt"
BATCH_SELECTION_TSV="${BATCH_DIR}/${BATCH_LABEL}.selection.tsv"

{
    echo -e "assembly_accession\tfasta_path"
    cat "$TASKS_FILE"
} >"$BATCH_SELECTION_TSV"

COMPLETED=0
PENDING=0
while IFS=$'\t' read -r acc _; do
    if [[ -f "${OUTDIR}/${acc}/snps.aligned.fa" ]]; then
        COMPLETED=$((COMPLETED + 1))
    else
        PENDING=$((PENDING + 1))
    fi
done <"$TASKS_FILE"

echo "=== Snippy-ctg Batch Runner ==="
echo "Plan: $PLAN"
echo "Batch label: $BATCH_LABEL"
echo "Included rows in full plan: $INCLUDED_TOTAL"
echo "Batch offset: $OFFSET"
echo "Reference: $REFERENCE"
echo "Selected rows this batch: $TOTAL"
echo "Already completed: $COMPLETED"
echo "Pending: $PENDING"
if [[ "$CPUS_SET" == true ]]; then
    echo "Requested CPUs per sample: $CPUS"
else
    echo "Requested CPUs per sample: auto"
fi
if [[ "$JOBS_SET" == true ]]; then
    echo "Requested parallel jobs: $JOBS"
else
    echo "Requested parallel jobs: auto"
fi
echo "Output dir: $OUTDIR"
echo "Batch selection: $BATCH_SELECTION_TSV"
echo "Docker image fallback: $DOCKER_IMAGE"
echo ""

if [[ "$DRY_RUN" == true ]]; then
    RESOLVE_PULL_DOCKER=false
fi
resolve_snippy
auto_tune_parallelism "$TOTAL"

echo "Execution backend: $( [[ "$USE_DOCKER" == true ]] && echo docker || echo path_or_configured_env )"
echo "Parallel policy: $AUTO_PARALLEL_MODE"
echo "Resolved CPUs per sample: $CPUS"
echo "Resolved parallel jobs: $JOBS"
echo ""

if [[ "$DRY_RUN" == true ]]; then
    echo "[DRY-RUN] No Snippy commands executed."
    exit 0
fi

printf "assembly_accession\tstatus\toutdir\tlog_path\ttimestamp\n" >"$SUMMARY_TSV"
printf "assembly_accession\tstatus\toutdir\tlog_path\ttimestamp\n" >"$BATCH_STATUS_TSV"
if [[ ! -f "$HISTORY_TSV" ]]; then
    printf "batch_label\tassembly_accession\tstatus\toutdir\tlog_path\ttimestamp\n" >"$HISTORY_TSV"
fi

export REFERENCE OUTDIR LOGDIR SUMMARY_TSV BATCH_STATUS_TSV HISTORY_TSV BATCH_LABEL LOCKFILE CPUS USE_DOCKER ROOT DOCKER_IMAGE
export -f run_snippy run_one

xargs -d '\n' -a "$TASKS_FILE" -P "$JOBS" -n 1 \
    bash -c 'IFS=$'"'"'\t'"'"' read -r acc fasta <<< "$1"; run_one "$acc" "$fasta"' _

OK=0
FAILED=0
SKIPPED=0
MISSING=0
while IFS=$'\t' read -r _ status _ _ _; do
    case "$status" in
        ok) OK=$((OK + 1)) ;;
        failed) FAILED=$((FAILED + 1)) ;;
        skipped_existing) SKIPPED=$((SKIPPED + 1)) ;;
        missing_fasta) MISSING=$((MISSING + 1)) ;;
    esac
done < <(tail -n +2 "$SUMMARY_TSV")

{
    echo "=== Snippy-ctg Run Summary ==="
    echo "Plan: $PLAN"
    echo "Batch label: $BATCH_LABEL"
    echo "Included rows in full plan: $INCLUDED_TOTAL"
    echo "Batch offset: $OFFSET"
    echo "Selected rows this batch: $TOTAL"
    echo "ok: $OK"
    echo "skipped_existing: $SKIPPED"
    echo "failed: $FAILED"
    echo "missing_fasta: $MISSING"
    echo "Parallel policy: $AUTO_PARALLEL_MODE"
    echo "Resolved CPUs per sample: $CPUS"
    echo "Resolved parallel jobs: $JOBS"
    echo "Status table: $SUMMARY_TSV"
    echo "Batch status table: $BATCH_STATUS_TSV"
    echo "Batch selection: $BATCH_SELECTION_TSV"
    echo "History table: $HISTORY_TSV"
} >"$SUMMARY_TXT"

cp "$SUMMARY_TXT" "$BATCH_SUMMARY_TXT"

cat "$SUMMARY_TXT"
