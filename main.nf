#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ANNOTATION } from './workflows/annotation'

workflow {
    if (!params.input) {
        error "Missing --input samplesheet"
    }
    if (!params.proteins) {
        error "Missing --proteins protein FASTA"
    }
    if (params.interproscan && !params.interproscan_db) {
        error "--interproscan requires --interproscan_db"
    }
    if (params.eggnog && !params.eggnog_db) {
        error "--eggnog requires --eggnog_db"
    }
    if (!params.run_edta && !params.repeatmasker_library && !params.edta_db) {
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

    ANNOTATION(
        samples,
        file(params.proteins, checkIfExists: true)
    )
}
