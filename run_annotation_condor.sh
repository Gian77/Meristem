#!/usr/bin/env bash
# Submit-node launcher following the nf-broom execution convention.
# Keep Nextflow in the foreground so Condor supervises the head process.
set -euo pipefail

export HOME=/mnt/cephfs/linuxhome/benucci
CONDA_ENV="${HOME}/.conda/envs/nextflow"
export JAVA_HOME="${CONDA_ENV}"
export JAVA_CMD="${CONDA_ENV}/bin/java"

PROJECT_DIR=/mnt/cephfs/linuxhome/benucci/Meristem
export PATH="${PROJECT_DIR}/bin:${CONDA_ENV}/bin:${PATH}"
INPUT=${INPUT:-"${PROJECT_DIR}/assets/F10702.samplesheet.csv"}
PROTEINS=${PROTEINS:?Set PROTEINS to the protein-evidence FASTA}
BUSCO_DB=${BUSCO_DB:?Set BUSCO_DB to the local BUSCO database directory}
OUTDIR=${OUTDIR:-"${PROJECT_DIR}/results/F10702"}
WORKDIR=${WORKDIR:-"${PROJECT_DIR}/nf-work-F10702"}

EXTRA_ARGS=()
if [[ "${RUN_BRAKER3:-true}" == "false" ]]; then
    EXTRA_ARGS+=(--run_braker3 false)
fi
if [[ "${RUN_BRAKER4:-false}" == "true" ]]; then
    BRAKER4_SOURCE=${BRAKER4_SOURCE:-"${PROJECT_DIR}/vendor/BRAKER4"}
    EXTRA_ARGS+=(--run_braker4 true --braker4_source "${BRAKER4_SOURCE}")
fi
if [[ "${RUN_HELIXER:-false}" == "true" ]]; then
    EXTRA_ARGS+=(--run_helixer true)
fi
if [[ -n "${HELIXER_MODEL:-}" ]]; then
    EXTRA_ARGS+=(--helixer_model "${HELIXER_MODEL}")
fi

cd "${PROJECT_DIR}"
mkdir -p "${OUTDIR}" "${WORKDIR}" "${HOME}/.apptainer_cache"

exec "${CONDA_ENV}/bin/nextflow" run "${PROJECT_DIR}/main.nf" \
    -profile condor \
    --input "${INPUT}" \
    --proteins "${PROTEINS}" \
    --cluster_proteins "${CLUSTER_PROTEINS:-false}" \
    --cluster_min_seq_id "${CLUSTER_MIN_SEQ_ID:-0.90}" \
    --cluster_coverage "${CLUSTER_COVERAGE:-0.80}" \
    --busco_lineage "${BUSCO_LINEAGE:-poales_odb10}" \
    --busco_db "${BUSCO_DB}" \
    --outdir "${OUTDIR}" \
    -w "${WORKDIR}" \
    -resume \
    "${EXTRA_ARGS[@]}"
