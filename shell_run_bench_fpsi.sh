#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
FPSI_BIN="${SCRIPT_DIR}/build/fpsi"

metrics=(0 1 2)
ns=(8 12 16)
dims=(2 6 10 15)
deltas=(10 60 250)
num_trials=5

print_help() {
  cat <<EOF
Usage:
  ${0##*/} [-metric values...] [-nn values...] [-dim values...]
             [-delta values...] [-trait value]

Defaults: metric=(0 1 2), nn=(8 12 16), dim=(2 6 10 15),
          delta=(10 60 250), trait=5
Metric:   0=Linf, 1=L1, 2=L2

A timestamped CSV file is always created beside this script.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
  -metric) shift; metrics=(); while [[ $# -gt 0 && "$1" != -* ]]; do metrics+=("$1"); shift; done ;;
  -nn) shift; ns=(); while [[ $# -gt 0 && "$1" != -* ]]; do ns+=("$1"); shift; done ;;
  -dim) shift; dims=(); while [[ $# -gt 0 && "$1" != -* ]]; do dims+=("$1"); shift; done ;;
  -delta) shift; deltas=(); while [[ $# -gt 0 && "$1" != -* ]]; do deltas+=("$1"); shift; done ;;
  -trait) num_trials="$2"; shift 2 ;;
  -h | -help) print_help; exit 0 ;;
  *) shift ;;
  esac
done

if [[ ! -x "${FPSI_BIN}" ]]; then
  echo "Benchmark binary not found or not executable: ${FPSI_BIN}" >&2
  exit 1
fi

output_file="${SCRIPT_DIR}/fpsi_cmp_results_$(date +%Y%m%d_%H%M%S).csv"

printf "[Size] [Metric] [Dim] [Delta] [Online_Com.(MB)] [Online(s)] [Offline_Com.(MB)] [Offline(s)]\n"

run_case() {
  local metric="$1" nn="$2" dim="$3" delta="$4"
  local args=(-dim "${dim}" -delta "${delta}" -metric "${metric}"
              -nn "${nn}" -trait "${num_trials}" -out "${output_file}")
  "${FPSI_BIN}" "${args[@]}"
}

for metric in "${metrics[@]}"; do
  for nn in "${ns[@]}"; do
    for dim in "${dims[@]}"; do
      for delta in "${deltas[@]}"; do
        run_case "${metric}" "${nn}" "${dim}" "${delta}"
      done
      echo
    done
  done
done

printf "Results written to %s\n" "${output_file}"
