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
if [[ "${RUN_BRAKER4:-true}" == "true" ]]; then
    EXTRA_ARGS+=(--run_braker4 true)
    # Defaults come from nextflow.config (the pinned, unmodified BRAKER4 and poales_odb12).
    if [[ -n "${BRAKER4_SOURCE:-}" ]]; then
        EXTRA_ARGS+=(--braker4_source "${BRAKER4_SOURCE}")
    fi
    if [[ -n "${BRAKER4_BUSCO_LINEAGE:-}" ]]; then
        EXTRA_ARGS+=(--braker4_busco_lineage "${BRAKER4_BUSCO_LINEAGE}")
    fi
fi
if [[ "${RUN_HELIXER:-false}" == "true" ]]; then
    EXTRA_ARGS+=(--run_helixer true)
fi
if [[ -n "${HELIXER_MODEL:-}" ]]; then
    EXTRA_ARGS+=(--helixer_model "${HELIXER_MODEL}")
fi
# The gene origin check is report-only and reads these shared files without changing them.
if [[ "${RUN_GENE_ORIGIN:-true}" == "true" ]]; then
    SHARED_DB=/mnt/cephfs/linuxhome/compbio/nonbackup/databases
    EXTRA_ARGS+=(--gene_origin_db "${GENE_ORIGIN_DB:-${SHARED_DB}/funannotate_databases/funannotate_db/uniprot_sprot.fasta}")
    EXTRA_ARGS+=(--gene_origin_taxonomy "${GENE_ORIGIN_TAXONOMY:-${SHARED_DB}/PlusPFP_protozoa_fungi_plants/nodes.dmp}")
fi
# Resume this run's own session, never whichever session was launched last in
# this directory. Condor restarts after an eviction find it the same way. Set
# RESUME_SESSION to reuse another session, or NEW_RUN=true for a first run.
SESSION=${RESUME_SESSION:-$("${PROJECT_DIR}/scripts/resolve_resume_session.sh" \
    "${PROJECT_DIR}/.nextflow/history" "${OUTDIR}" "${WORKDIR}")}
RESUME_ARGS=()
if [[ -n "${SESSION}" ]]; then
    echo "Resuming Nextflow session ${SESSION}"
    RESUME_ARGS=(-resume "${SESSION}")
elif [[ "${NEW_RUN:-false}" != "true" ]]; then
    echo "No Nextflow session found for OUTDIR=${OUTDIR} and WORKDIR=${WORKDIR}." >&2
    echo "Set NEW_RUN=true to start a new run, or RESUME_SESSION to reuse a session." >&2
    exit 1
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
    "${RESUME_ARGS[@]}" \
    "${EXTRA_ARGS[@]}"
