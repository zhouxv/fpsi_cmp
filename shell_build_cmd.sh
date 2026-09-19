#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"

if [[ "${1:-}" == "--clean" ]]; then
  rm -rf "${BUILD_DIR}"
elif [[ $# -gt 0 ]]; then
  printf 'Usage: %s [--clean]\n' "${0##*/}" >&2
  exit 1
fi

printf 'Building the project with CMake...\n'
cmake -S "${SCRIPT_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${BUILD_DIR}" --parallel "$(nproc)"
printf 'Build completed successfully: %s\n' "${BUILD_DIR}"
