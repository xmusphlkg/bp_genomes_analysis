#!/usr/bin/env bash

if [[ -n "${STEP4_DISTRIBUTED_RAW_READS_LIB_LOADED:-}" ]]; then
  return 0
fi
STEP4_DISTRIBUTED_RAW_READS_LIB_LOADED=1

step4_init_distributed_env() {
  local env_file="$1"

  if [[ "$env_file" != */* && -f "./$env_file" ]]; then
    env_file="./$env_file"
  fi
  if [[ ! -f "$env_file" ]]; then
    echo "ERROR: env file not found: $env_file" >&2
    exit 1
  fi

  set -a
  source "$env_file"
  set +a

  if [[ -n "${SSH_KEY:-}" ]]; then
    SSH_KEY="${SSH_KEY/#\~/$HOME}"
  fi

  STEP4_REPO_DIR="${LOCAL_REPO_DIR:-$(pwd)}"
  STEP4_SSH_STRICT_HOST_KEY_CHECKING="${SSH_STRICT_HOST_KEY_CHECKING:-accept-new}"
  # shellcheck disable=SC1091
  source "${STEP4_REPO_DIR}/workflow/lib/runtime_envs.sh"
  project_env_load_config "${STEP4_REPO_DIR}"
  if [[ -n "${DATA_ROOT:-}" ]]; then
    STEP4_DATA_ROOT_DEFAULT="${DATA_ROOT}"
  else
    STEP4_DATA_ROOT_DEFAULT="$(project_module_data_root step4_prn_validation)"
  fi
}

step4_need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: required command not found: $1" >&2
    exit 1
  }
}

step4_server_value() {
  local idx="$1" key="$2"
  local primary="SERVER${idx}_${key}"
  local legacy="${key}${idx}"
  local value="${!primary:-}"
  if [[ -z "$value" ]]; then
    case "$key" in
      HOST|USER|PORT|PASS)
        value="${!legacy:-}"
        ;;
    esac
  fi
  printf '%s' "$value"
}

step4_server_has_remote() {
  local idx="$1"
  [[ -n "$(step4_server_value "$idx" HOST)" && -n "$(step4_server_value "$idx" USER)" ]]
}

step4_resolve_repo_for_server() {
  local idx="$1"
  local explicit user
  explicit="$(step4_server_value "$idx" REPO_DIR)"
  user="$(step4_server_value "$idx" USER)"
  if [[ -n "$explicit" ]]; then
    printf '%s' "$explicit"
  elif [[ -n "$user" ]]; then
    printf '/home/%s/pertussis/pertussis_gene' "$user"
  else
    printf '%s' "$STEP4_REPO_DIR"
  fi
}

step4_resolve_data_root_for_server() {
  local idx="$1"
  local explicit
  explicit="$(step4_server_value "$idx" DATA_ROOT)"
  if [[ -n "$explicit" ]]; then
    printf '%s' "$explicit"
  else
    printf '%s' "$STEP4_DATA_ROOT_DEFAULT"
  fi
}

step4_should_run_shard() {
  local shard="$1" filter="${2:-}"
  if [[ -z "$filter" ]]; then
    return 0
  fi
  [[ ",$filter," == *",$shard,"* ]]
}

step4_run_local() {
  local dry_run="${1:-0}"
  shift
  if [[ "$dry_run" -eq 1 ]]; then
    echo "[dry-run] $*"
  else
    eval "$@"
  fi
}

step4_run_remote() {
  local idx="$1" cmd="$2"
  local host user port pass
  host="$(step4_server_value "$idx" HOST)"
  user="$(step4_server_value "$idx" USER)"
  port="$(step4_server_value "$idx" PORT)"
  pass="$(step4_server_value "$idx" PASS)"
  [[ -n "$port" ]] || port="22"

  if [[ -n "${SSH_KEY:-}" ]]; then
    ssh -i "$SSH_KEY" -p "$port" -o BatchMode=yes -o ConnectTimeout=8 -o StrictHostKeyChecking="${STEP4_SSH_STRICT_HOST_KEY_CHECKING}" "$user@$host" "$cmd"
  elif [[ -n "$pass" ]]; then
    step4_need_cmd sshpass
    sshpass -p "$pass" ssh -p "$port" -o ConnectTimeout=8 -o StrictHostKeyChecking="${STEP4_SSH_STRICT_HOST_KEY_CHECKING}" "$user@$host" "$cmd"
  else
    ssh -p "$port" -o BatchMode=yes -o ConnectTimeout=8 -o StrictHostKeyChecking="${STEP4_SSH_STRICT_HOST_KEY_CHECKING}" "$user@$host" "$cmd"
  fi
}

step4_copy_to_remote() {
  local idx="$1" src="$2" dst="$3"
  local host user port pass
  host="$(step4_server_value "$idx" HOST)"
  user="$(step4_server_value "$idx" USER)"
  port="$(step4_server_value "$idx" PORT)"
  pass="$(step4_server_value "$idx" PASS)"
  [[ -n "$port" ]] || port="22"

  if [[ -n "${SSH_KEY:-}" ]]; then
    scp -i "$SSH_KEY" -P "$port" -o BatchMode=yes -o ConnectTimeout=8 -o StrictHostKeyChecking="${STEP4_SSH_STRICT_HOST_KEY_CHECKING}" "$src" "$user@$host:$dst"
  elif [[ -n "$pass" ]]; then
    step4_need_cmd sshpass
    sshpass -p "$pass" scp -P "$port" -o ConnectTimeout=8 -o StrictHostKeyChecking="${STEP4_SSH_STRICT_HOST_KEY_CHECKING}" "$src" "$user@$host:$dst"
  else
    scp -P "$port" -o BatchMode=yes -o ConnectTimeout=8 -o StrictHostKeyChecking="${STEP4_SSH_STRICT_HOST_KEY_CHECKING}" "$src" "$user@$host:$dst"
  fi
}

step4_build_ssh_cmd() {
  local idx="$1"
  local port
  port="$(step4_server_value "$idx" PORT)"
  [[ -n "$port" ]] || port="22"
  if [[ -n "${SSH_KEY:-}" ]]; then
    printf 'ssh -p %s -o BatchMode=yes -o ConnectTimeout=8 -o StrictHostKeyChecking=%s -i %q' \
      "$port" "$STEP4_SSH_STRICT_HOST_KEY_CHECKING" "$SSH_KEY"
  else
    printf 'ssh -p %s -o BatchMode=yes -o ConnectTimeout=8 -o StrictHostKeyChecking=%s' \
      "$port" "$STEP4_SSH_STRICT_HOST_KEY_CHECKING"
  fi
}

step4_path_run_list() {
  local repo="$1" shard_dir_rel="$2" shard="$3"
  printf '%s/%s/bp_raw_reads_runs.%s.txt' "$repo" "$shard_dir_rel" "$shard"
}

step4_path_plan_tsv() {
  local repo="$1" shard_dir_rel="$2" shard="$3"
  printf '%s/%s/bp_raw_reads_download_plan.%s.tsv' "$repo" "$shard_dir_rel" "$shard"
}

step4_path_workdir() {
  local idx="$1" shard="$2"
  local repo data_root
  repo="$(step4_resolve_repo_for_server "$idx")"
  data_root="$(step4_resolve_data_root_for_server "$idx")"
  if [[ -n "$data_root" ]]; then
    printf '%s/work/%s' "$data_root" "$shard"
  else
    printf '%s/modules/step4_prn_validation/work/%s' "$repo" "$shard"
  fi
}

step4_path_outdir() {
  local idx="$1" shard="$2"
  local repo data_root
  repo="$(step4_resolve_repo_for_server "$idx")"
  data_root="$(step4_resolve_data_root_for_server "$idx")"
  if [[ -n "$data_root" ]]; then
    printf '%s/outputs/assemblies/%s' "$data_root" "$shard"
  else
    printf '%s/modules/step4_prn_validation/outputs/assemblies/%s' "$repo" "$shard"
  fi
}

step4_path_launcher_log() {
  local idx="$1" shard="$2"
  printf '%s/launcher.log' "$(step4_path_workdir "$idx" "$shard")"
}

step4_path_launcher_pid() {
  local idx="$1" shard="$2"
  printf '%s/launcher.pid' "$(step4_path_workdir "$idx" "$shard")"
}

step4_path_launcher_command() {
  local idx="$1" shard="$2"
  printf '%s/launcher.command.sh' "$(step4_path_workdir "$idx" "$shard")"
}

step4_resolve_home_for_server() {
  local idx="$1"
  local explicit user
  explicit="$(step4_server_value "$idx" HOME_DIR)"
  user="$(step4_server_value "$idx" USER)"
  if [[ -n "$explicit" ]]; then
    printf '%s' "$explicit"
  elif [[ -n "$user" ]]; then
    printf '/home/%s' "$user"
  else
    printf '%s' "$HOME"
  fi
}

step4_path_tools_root() {
  local idx="$1"
  printf '%s/.local/share/bp_step4/aspera-cli' "$(step4_resolve_home_for_server "$idx")"
}

step4_path_aspera_bin_file() {
  local idx="$1" filename="$2"
  printf '%s/bin/%s' "$(step4_path_tools_root "$idx")" "$filename"
}

step4_path_aspera_etc_file() {
  local idx="$1" filename="$2"
  printf '%s/etc/%s' "$(step4_path_tools_root "$idx")" "$filename"
}

step4_path_aspera_bin() {
  local idx="$1"
  step4_path_aspera_bin_file "$idx" "ascp"
}

step4_path_aspera_key() {
  local idx="$1"
  step4_path_aspera_etc_file "$idx" "asperaweb_id_dsa.openssh"
}

step4_find_local_ascp() {
  local candidate
  if [[ -n "${ASPERA_ASCP:-}" && -x "${ASPERA_ASCP/#\~/$HOME}" ]]; then
    printf '%s' "${ASPERA_ASCP/#\~/$HOME}"
    return 0
  fi
  if command -v ascp >/dev/null 2>&1; then
    command -v ascp
    return 0
  fi
  for candidate in \
    "$HOME/miniconda3/bin/ascp" \
    "$HOME/miniforge3/bin/ascp" \
    "/opt/miniconda3/bin/ascp" \
    "/opt/miniconda/bin/ascp"; do
    if [[ -x "$candidate" ]]; then
      printf '%s' "$candidate"
      return 0
    fi
  done
  return 1
}

step4_find_local_aspera_key() {
  local ascp_bin="$1"
  local candidate key_root
  if [[ -n "${ASPERA_KEY:-}" && -f "${ASPERA_KEY/#\~/$HOME}" ]]; then
    printf '%s' "${ASPERA_KEY/#\~/$HOME}"
    return 0
  fi
  if [[ -n "$ascp_bin" ]]; then
    key_root="$(cd "$(dirname "$ascp_bin")/../etc" 2>/dev/null && pwd -P || true)"
    candidate="${key_root}/asperaweb_id_dsa.openssh"
    if [[ -f "$candidate" ]]; then
      printf '%s' "$candidate"
      return 0
    fi
  fi
  for candidate in \
    "$HOME/miniconda3/etc/asperaweb_id_dsa.openssh" \
    "$HOME/miniforge3/etc/asperaweb_id_dsa.openssh" \
    "/opt/miniconda3/etc/asperaweb_id_dsa.openssh" \
    "/opt/miniconda/etc/asperaweb_id_dsa.openssh"; do
    if [[ -f "$candidate" ]]; then
      printf '%s' "$candidate"
      return 0
    fi
  done
  return 1
}

step4_install_aspera_assets() {
  local idx="$1" source_ascp="$2" source_key="$3"
  local source_bin_dir source_etc_dir name
  source_bin_dir="$(cd "$(dirname "$source_ascp")" 2>/dev/null && pwd -P)"
  source_etc_dir="$(cd "$(dirname "$source_key")" 2>/dev/null && pwd -P)"
  mkdir -p "$(step4_path_tools_root "$idx")/bin" "$(step4_path_tools_root "$idx")/etc"

  for name in ascp ascp4 aspera .aspera_cli_conf; do
    if [[ -f "$source_bin_dir/$name" ]]; then
      install -m 755 "$source_bin_dir/$name" "$(step4_path_aspera_bin_file "$idx" "$name")"
    fi
  done

  for name in aspera-license aspera.conf aspera_tokenauth_id_rsa asperaweb_id_dsa.openssh asperaweb_id_dsa.putty; do
    if [[ -f "$source_etc_dir/$name" ]]; then
      install -m 644 "$source_etc_dir/$name" "$(step4_path_aspera_etc_file "$idx" "$name")"
    fi
  done
}
