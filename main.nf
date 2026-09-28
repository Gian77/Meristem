#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ANNOTATION } from './workflows/annotation'

workflow {
    def run_braker3 = params.run_braker3.toString().toBoolean()
    def run_braker4 = params.run_braker4.toString().toBoolean()
    def run_helixer = params.run_helixer.toString().toBoolean()
    def run_edta = params.run_edta.toString().toBoolean()
    def run_repeatmasker = params.run_repeatmasker.toString().toBoolean()
    def run_interproscan = params.interproscan.toString().toBoolean()
    def run_eggnog = params.eggnog.toString().toBoolean()
    if (!params.input) {
        error "Missing --input samplesheet"
    }
    if (!run_braker3 && !run_braker4 && !run_helixer) {
        error "Enable at least one annotation method: --run_braker3, --run_braker4, or --run_helixer"
    }
    if ((run_braker3 || run_braker4) && !params.proteins) {
        error "Missing --proteins protein FASTA"
    }
    if (run_braker4 && !params.braker4_source) {
        error "--run_braker4 requires --braker4_source pointing to a BRAKER4 checkout"
    }
    if (run_helixer && !params.helixer_model) {
        error "--run_helixer requires --helixer_model pointing to a local HDF5 model"
    }
    cluster_min_seq_id = params.cluster_min_seq_id.toString().toBigDecimal()
    cluster_coverage = params.cluster_coverage.toString().toBigDecimal()
    if (cluster_min_seq_id < 0 || cluster_min_seq_id > 1) {
        error "--cluster_min_seq_id must be between 0 and 1"
    }
    if (cluster_coverage <= 0 || cluster_coverage > 1) {
        error "--cluster_coverage must be greater than 0 and no more than 1"
    }
    if (run_interproscan && !params.interproscan_db) {
        error "--interproscan requires --interproscan_db"
    }
    if (run_eggnog && !params.eggnog_db) {
        error "--eggnog requires --eggnog_db"
    }
    if (!run_edta && !params.repeatmasker_library && !params.edta_db) {
        error "When --run_edta false, provide --repeatmasker_library or --edta_db"
    }

    samples = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !row.assembly) {
                error "Each samplesheet row requires sample and assembly columns"
            }
            tuple(row.sample.toString(), file(row.assembly.toString()))
        }

    proteins = params.proteins ? file(params.proteins, checkIfExists: true) : null
    braker4_source = params.braker4_source ? file(params.braker4_source, checkIfExists: true) : null
    helixer_model = run_helixer && params.helixer_model ? file(params.helixer_model, checkIfExists: true) : null

    ANNOTATION(samples, proteins, braker4_source, helixer_model)
}
