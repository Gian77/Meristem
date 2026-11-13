#!/usr/bin/env bash
# Install the BRAKER4 version that Meristem is tested with, unmodified, and the
# Snakemake environment that drives it. Run once by hand before the first run;
# Meristem itself never downloads BRAKER4.
#
# Usage: scripts/setup_braker4.sh [DEST] [SNAKEMAKE_ENV]
#   DEST           BRAKER4 folder (default: <Meristem>/vendor/BRAKER4-<short commit>)
#   SNAKEMAKE_ENV  Python virtual environment for Snakemake (default: <Meristem>/.snakemake-env)
# Environment:
#   PYTHON         Python 3.11 or newer used to build the environment (default: python3)
#
# Existing folders are reused, never overwritten or deleted.
set -euo pipefail

BRAKER4_COMMIT=6a5ee16591a0a15872aa838d9dee68e54f614ee0
BRAKER4_REPO=Gaius-Augustus/BRAKER4

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
dest=${1:-"${project_dir}/vendor/BRAKER4-${BRAKER4_COMMIT:0:7}"}
venv=${2:-"${project_dir}/.snakemake-env"}
python=${PYTHON:-python3}

# 1. BRAKER4 source at the pinned commit, as published (MIT licence)
if [[ -f "${dest}/MERISTEM_BRAKER4_SOURCE.txt" ]]; then
    echo "BRAKER4 already installed: ${dest}"
    cat "${dest}/MERISTEM_BRAKER4_SOURCE.txt"
elif [[ -e "${dest}" ]]; then
    echo "ERROR: ${dest} exists but was not installed by this script. Choose another DEST." >&2
    exit 1
else
    tmp=$(mktemp -d "${TMPDIR:-/tmp}/braker4-setup.XXXXXX")
    url="https://codeload.github.com/${BRAKER4_REPO}/tar.gz/${BRAKER4_COMMIT}"
    echo "Downloading BRAKER4 ${BRAKER4_COMMIT} from ${url}"
    curl -fsSL --retry 3 -o "${tmp}/braker4.tar.gz" "${url}"
    tar -xzf "${tmp}/braker4.tar.gz" -C "${tmp}"
    src="${tmp}/BRAKER4-${BRAKER4_COMMIT}"
    test -f "${src}/Snakefile" || { echo "ERROR: archive has no Snakefile" >&2; exit 1; }
    mkdir -p "$(dirname "${dest}")"
    mv "${src}" "${dest}"
    {
        printf 'repository\thttps://github.com/%s\n' "${BRAKER4_REPO}"
        printf 'commit\t%s\n' "${BRAKER4_COMMIT}"
        printf 'archive_sha256\t%s\n' "$(sha256sum "${tmp}/braker4.tar.gz" | cut -d' ' -f1)"
        printf 'installed\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'modified\tno\n'
    } > "${dest}/MERISTEM_BRAKER4_SOURCE.txt"
    echo "Installed BRAKER4 into ${dest}"
    echo "The downloaded archive is left in ${tmp}; delete it when convenient."
fi

# 2. Snakemake runtime, pinned in requirements-snakemake.txt
if [[ -x "${venv}/bin/snakemake" ]]; then
    echo "Snakemake environment already present: ${venv} ($("${venv}/bin/snakemake" --version))"
else
    "${python}" - <<'EOF'
import sys
if sys.version_info < (3, 11):
    sys.exit(f"ERROR: Python 3.11 or newer is needed for Snakemake 8, found {sys.version.split()[0]}. Set PYTHON.")
EOF
    "${python}" -m venv --copies "${venv}"
    "${venv}/bin/pip" install --no-cache-dir -r "${project_dir}/requirements-snakemake.txt"
    echo "Built Snakemake environment: ${venv} ($("${venv}/bin/snakemake" --version))"
fi

cat <<EOF

Done. Meristem uses these by default. If you chose other locations, pass:
  --braker4_source ${dest} --braker4_snakemake ${venv}/bin/snakemake
EOF
