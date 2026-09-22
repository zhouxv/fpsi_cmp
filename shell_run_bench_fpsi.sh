#!/usr/bin/env bash
set -euo pipefail

# Resolve paths relative to this script so it can run from any directory.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
FPSI_BIN="${SCRIPT_DIR}/build/fpsi"
NETWORK_SCRIPT="${SCRIPT_DIR}/shell_config_network.sh"

# Parse user options first, then fill in missing values from the preset.
# Explicit options override preset defaults regardless of argument order.
metrics=()
ns=()
dims=()
deltas=()
num_trials=""
interface="lo"
output_dir="${SCRIPT_DIR}"
dry_run=false
preset="quick"

print_help() {
  cat <<'EOF'
Usage:
  shell_run_bench_fpsi.sh [options]

The network is configured separately with shell_config_network.sh. This script
detects the current tcconfig settings and prints them before the benchmark.

Presets:
  full   Complete paper matrix, 81 cases and 3 trials per case
  quick  Representative Table 2 subset, 18 cases and 1 trial per case (default)

Full matrix:
  metric = 0 1 2          (0=Linf, 1=L1, 2=L2)
  nn     = 8 12 16        (set size N=2^nn)
  dim    = 2 6 10
  delta  = 10 60 250
  trials = 3

Quick matrix:
  metric = 0 1 2
  nn     = 12
  dim    = 2 6 10
  delta  = 10 250
  trials = 1

Options:
  --preset NAME                  full or quick (default: quick)
  --metric, -metric VALUES...   Metrics to run (0, 1, 2)
  --nn, -nn VALUES...           Base-2 set-size exponents
  --dim, -dim VALUES...         Dimensions
  --delta, -delta VALUES...     Distance thresholds
  --trials, -trait N            Trials averaged per case
  --interface IFACE             Interface inspected by tcconfig (default: lo)
  --output-dir DIR              Result directory (default: project root)
  --dry-run                     Print cases without executing them
  -h, --help, -help             Print this help message

Examples:
  ./shell_config_network.sh lan
  ./shell_run_bench_fpsi.sh --preset full

  ./shell_config_network.sh wan
  ./shell_run_bench_fpsi.sh --preset quick
EOF
}

# 1. Parse options. Experiment parameters accept lists; other options take one value.
while [[ $# -gt 0 ]]; do
  case "$1" in
    --metric|-metric|--nn|-nn|--dim|-dim|--delta|-delta)
      option="$1"
      shift
      values=()
      # Collect values until the next option, e.g. --dim 2 6 10 --delta 10.
      while [[ $# -gt 0 && "$1" != -* ]]; do
        values+=("$1")
        shift
      done
      if [[ ${#values[@]} -eq 0 ]]; then
        printf 'Error: %s requires at least one value.\n' "${option}" >&2
        exit 1
      fi
      case "${option}" in
        --metric|-metric) metrics=("${values[@]}") ;;
        --nn|-nn) ns=("${values[@]}") ;;
        --dim|-dim) dims=("${values[@]}") ;;
        --delta|-delta) deltas=("${values[@]}") ;;
      esac
      ;;
    --preset|--trials|-trait|--interface|--output-dir)
      if [[ $# -lt 2 || -z "$2" || "$2" == -* ]]; then
        printf 'Error: %s requires a value.\n' "$1" >&2
        exit 1
      fi
      case "$1" in
        --preset) preset="$2" ;;
        --trials|-trait) num_trials="$2" ;;
        --interface) interface="$2" ;;
        --output-dir) output_dir="$2" ;;
      esac
      shift 2
      ;;
    --dry-run)
      dry_run=true
      shift
      ;;
    -h|--help|-help)
      print_help
      exit 0
      ;;
    *)
      printf 'Error: unknown option: %s\n\n' "$1" >&2
      print_help >&2
      exit 1
      ;;
  esac
done

# 2. Select the preset and fill in parameters the user did not specify.
case "${preset}" in
  full)
    default_ns=(8 12 16)
    default_deltas=(10 60 250)
    default_trials=3
    ;;
  quick)
    default_ns=(12)
    default_deltas=(10 250)
    default_trials=1
    ;;
  *)
    printf 'Error: invalid preset %s (expected full or quick).\n' "${preset}" >&2
    exit 1
    ;;
esac

if [[ ${#metrics[@]} -eq 0 ]]; then metrics=(0 1 2); fi
if [[ ${#dims[@]} -eq 0 ]]; then dims=(2 6 10); fi
if [[ ${#ns[@]} -eq 0 ]]; then ns=("${default_ns[@]}"); fi
if [[ ${#deltas[@]} -eq 0 ]]; then deltas=("${default_deltas[@]}"); fi
if [[ -z "${num_trials}" ]]; then num_trials="${default_trials}"; fi

# Validate the experiment parameters before starting the benchmark.
for metric in "${metrics[@]}"; do
  if [[ ! "${metric}" =~ ^[012]$ ]]; then
    printf 'Error: invalid metric %s (expected 0, 1, or 2).\n' "${metric}" >&2
    exit 1
  fi
done
for nn in "${ns[@]}"; do
  if [[ ! "${nn}" =~ ^([1-9]|[1-5][0-9]|6[0-2])$ ]]; then
    printf 'Error: invalid nn %s (expected an integer from 1 to 62).\n' "${nn}" >&2
    exit 1
  fi
done
for value in "${dims[@]}" "${deltas[@]}" "${num_trials}"; do
  if [[ ! "${value}" =~ ^[1-9][0-9]*$ ]]; then
    printf 'Error: dimensions, deltas, and trials must be positive integers: %s\n' \
      "${value}" >&2
    exit 1
  fi
done

# 3. Read the current network settings for the startup message and CSV filename.
# Use shell_config_network.sh separately to change the network settings.
if [[ ! -x "${NETWORK_SCRIPT}" ]]; then
  printf 'Error: network helper not found or not executable: %s\n' \
    "${NETWORK_SCRIPT}" >&2
  exit 1
fi

network_metadata="$("${NETWORK_SCRIPT}" metadata --interface "${interface}")"
IFS=$'\t' read -r network_profile network_interface bandwidth rtt \
  <<<"${network_metadata}"

printf 'Detected network configuration:\n'
printf '  Profile   : %s\n' "${network_profile}"
printf '  Interface : %s\n' "${network_interface}"
printf '  Rate      : %s\n' "${bandwidth}"
printf '  RTT       : %s\n\n' "${rtt}"

case_count=$((${#metrics[@]} * ${#ns[@]} * ${#dims[@]} * ${#deltas[@]}))
printf 'Benchmark matrix: %d cases, %d trial(s) per case\n' \
  "${case_count}" "${num_trials}"
printf 'Preset: %s\n' "${preset}"

timestamp="$(date +%Y%m%d_%H%M%S)"
output_file="${output_dir}/fpsi_cmp_results_${network_profile}_${timestamp}.csv"

# A dry run only prints commands; it needs no binary and creates no output directory.
if [[ "${dry_run}" == false ]]; then
  if [[ ! -x "${FPSI_BIN}" ]]; then
    printf 'Error: benchmark binary not found or not executable: %s\n' \
      "${FPSI_BIN}" >&2
    printf 'Build it first with ./shell_build_cmd.sh\n' >&2
    exit 1
  fi
  mkdir -p "${output_dir}"
  printf '\n[Size] [Metric] [Dim] [Delta] [Offline_Com.(MB)] [Offline(s)] '
  printf '[Online_Com.(MB)] [Online(s)] [Total_Com.(MB)] [Total(s)]\n'
fi

# 4. Run each metric, set size, dimension, and threshold combination in order.
# fpsi averages trials and appends each result to one CSV. Stop on failure.
for metric in "${metrics[@]}"; do
  for nn in "${ns[@]}"; do
    for dim in "${dims[@]}"; do
      for delta in "${deltas[@]}"; do
        args=(-metric "${metric}" -nn "${nn}" -dim "${dim}"
              -delta "${delta}" -trait "${num_trials}" -out "${output_file}")
        if [[ "${dry_run}" == true ]]; then
          # Shell-escape arguments so the printed command can be copied directly.
          printf '%q ' "${FPSI_BIN}" "${args[@]}"
          printf '\n'
        else
          "${FPSI_BIN}" "${args[@]}"
        fi
      done
    done
  done
done

if [[ "${dry_run}" == true ]]; then
  printf 'Dry run complete; no benchmark was executed.\n'
else
  printf 'Results written to %s\n' "${output_file}"
fi
