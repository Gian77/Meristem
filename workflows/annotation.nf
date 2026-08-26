nextflow.enable.dsl=2

include { VALIDATE_ASSEMBLY } from '../modules/local/validate_assembly'
include { EDTA_DISCOVERY } from '../modules/local/edta_discovery'
include { REPEATMASKER } from '../modules/local/repeatmasker'
include { BRAKER3 } from '../modules/local/braker3'
include { BUSCO } from '../modules/local/busco'
include { INTERPROSCAN } from '../modules/local/interproscan'
include { EGGNOG } from '../modules/local/eggnog'
include { FINAL_REPORT } from '../modules/local/final_report'

workflow ANNOTATION {
    take:
    samples
    proteins

    main:
    validated = VALIDATE_ASSEMBLY(samples)

    if (params.run_edta) {
        edta_out = EDTA_DISCOVERY(validated)
        repeats = edta_out.repeats
    } else {
        repeats = validated.map { sample, assembly, metadata -> tuple(sample, assembly, file(params.repeatmasker_library ?: params.edta_db), metadata) }
    }

    if (params.run_repeatmasker) {
        repeatmasker_out = REPEATMASKER(repeats)
        masked = repeatmasker_out.masked
    } else {
        masked = repeats.map { sample, assembly, repeat_library, metadata -> tuple(sample, assembly, repeat_library, metadata) }
    }
    braker_out = BRAKER3(masked, proteins)
    predictions = braker_out.predictions
    busco_out = BUSCO(predictions)
    qc = busco_out.qc

    if (params.interproscan) {
        INTERPROSCAN(predictions)
    }
    if (params.eggnog) {
        EGGNOG(predictions)
    }

    FINAL_REPORT(qc)
}
