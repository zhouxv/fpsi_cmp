#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
THIRDPARTY_DIR="${SCRIPT_DIR}/thirdparty"
INSTALL_DIR="${SCRIPT_DIR}/install"

VOLEPSI_REPOSITORY="https://github.com/ladnir/volepsi.git"
VOLEPSI_COMMIT="ed943f5f814591cdf864777c73b7bc9e7526c1a8"
TCCONFIG_REPOSITORY="https://github.com/thombashi/tcconfig.git"
TCCONFIG_COMMIT="ee55d275603f2cb8c10c3fbe84989ea71149f5de"
BOOST_URL="https://sourceforge.net/projects/boost/files/boost/1.86.0/boost_1_86_0.tar.bz2/download"
BOOST_SHA256="1bed88e40401b2cb7a1f76d4bab499e352fa4d0c5f31c0dbae64e24d34d7513b"

install_volepsi_requested=true
install_tcconfig_requested=true

case "${1:-}" in
  "") ;;
  --volepsi-only) install_tcconfig_requested=false ;;
  --tcconfig-only) install_volepsi_requested=false ;;
  -h|--help)
    cat <<'EOF'
Usage: shell_install_dependencies.sh [--volepsi-only|--tcconfig-only]

Without an option, installs both pinned dependencies under ./install.
EOF
    exit 0
    ;;
  *)
    printf 'Error: unknown option: %s\n' "$1" >&2
    exit 1
    ;;
esac

for command_name in git curl cmake python3 sha256sum; do
  if ! command -v "${command_name}" >/dev/null 2>&1; then
    printf 'Error: required command not found: %s\n' "${command_name}" >&2
    exit 1
  fi
done

if [[ "${install_tcconfig_requested}" == true ]] &&
   ! python3 -c 'import ensurepip' >/dev/null 2>&1; then
  cat >&2 <<'EOF'
Error: Python venv support is unavailable.
Install python3-venv (or the matching python3.X-venv package) and retry.
EOF
  exit 1
fi

install_volepsi() {
  local source_dir="${THIRDPARTY_DIR}/volepsi"
  local archive="${source_dir}/out/boost_1_86_0.tar.bz2"

  rm -rf "${source_dir}"
  printf 'Cloning volePSI at %s\n' "${VOLEPSI_COMMIT}"
  git clone "${VOLEPSI_REPOSITORY}" "${source_dir}"
  git -C "${source_dir}" checkout --detach "${VOLEPSI_COMMIT}"

  sed -i '/-DENABLE_SILENT_VOLE=ON/a\                       -DENABLE_FOLEAGE=ON' \
    "${source_dir}/thirdparty/getLibOTe.cmake"
  sed -i \
    's/set(libOTe_options silentot silent_vole circuits)/set(libOTe_options silentot silent_vole circuits foleage)/' \
    "${source_dir}/cmake/findDependancies.cmake"

  mkdir -p "${source_dir}/out"
  curl -fL --retry 3 -o "${archive}" "${BOOST_URL}"
  printf '%s  %s\n' "${BOOST_SHA256}" "${archive}" | sha256sum --check --status || {
    echo 'Error: Boost archive checksum verification failed.' >&2
    exit 1
  }

  (
    cd "${source_dir}"
    python3 build.py -DVOLE_PSI_ENABLE_BOOST=ON
    python3 build.py --install="${INSTALL_DIR}/volepsi"
  )
  cp "${source_dir}/out/build/linux/volePSI/config.h" \
    "${INSTALL_DIR}/volepsi/include/volePSI/config.h"
}

install_tcconfig() {
  local source_dir="${THIRDPARTY_DIR}/tcconfig"
  local venv_dir="${INSTALL_DIR}/tcconfig"

  rm -rf "${source_dir}" "${venv_dir}"
  printf 'Cloning tcconfig at %s\n' "${TCCONFIG_COMMIT}"
  git clone "${TCCONFIG_REPOSITORY}" "${source_dir}"
  git -C "${source_dir}" checkout --detach "${TCCONFIG_COMMIT}"

  python3 -m venv "${venv_dir}"
  "${venv_dir}/bin/python" -m pip install --no-cache-dir "${source_dir}"
}

mkdir -p "${THIRDPARTY_DIR}" "${INSTALL_DIR}"
if [[ "${install_volepsi_requested}" == true ]]; then
  install_volepsi
fi
if [[ "${install_tcconfig_requested}" == true ]]; then
  install_tcconfig
fi

printf '\nDependencies installed under %s\n' "${INSTALL_DIR}"
printf 'tcconfig commands are under %s/tcconfig/bin\n' "${INSTALL_DIR}"
