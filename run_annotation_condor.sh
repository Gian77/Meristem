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

cd "${PROJECT_DIR}"
mkdir -p "${OUTDIR}" "${WORKDIR}" "${HOME}/.apptainer_cache"

exec "${CONDA_ENV}/bin/nextflow" run "${PROJECT_DIR}/main.nf" \
    -profile condor \
    --input "${INPUT}" \
    --proteins "${PROTEINS}" \
    --busco_lineage "${BUSCO_LINEAGE:-poales_odb10}" \
    --busco_db "${BUSCO_DB}" \
    --outdir "${OUTDIR}" \
    -w "${WORKDIR}" \
    -resume
