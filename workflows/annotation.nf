nextflow.enable.dsl=2

include { VALIDATE_ASSEMBLY } from '../modules/local/validate_assembly'
include { EDTA_DISCOVERY } from '../modules/local/edta_discovery'
include { REPEATMASKER } from '../modules/local/repeatmasker'
include { BRAKER3 } from '../modules/local/braker3'
include { BRAKER4 } from '../modules/local/braker4'
include { HELIXER } from '../modules/local/helixer'
include { ANNOTATION_COMPARISON } from '../modules/local/annotation_comparison'
include { HELIXER_PROTEINS; ANNOTATION_BUSCO } from '../modules/local/annotation_busco'
include { CLUSTER_PROTEINS } from '../modules/local/cluster_proteins'
include { BUSCO } from '../modules/local/busco'
include { INTERPROSCAN } from '../modules/local/interproscan'
include { EGGNOG } from '../modules/local/eggnog'
include { FINAL_REPORT } from '../modules/local/final_report'
include { ANNOTATION_BUSCO_COMPARISON; ANNOTATION_FINAL_REPORT; TOOLS_REPORT } from '../modules/local/reports'

workflow ANNOTATION {
    take:
    samples
    proteins
    braker4_source
    helixer_model

    main:
    def run_braker3 = params.run_braker3.toString().toBoolean()
    def run_braker4 = params.run_braker4.toString().toBoolean()
    def run_helixer = params.run_helixer.toString().toBoolean()
    def run_edta = params.run_edta.toString().toBoolean()
    def run_repeatmasker = params.run_repeatmasker.toString().toBoolean()
    def run_interproscan = params.interproscan.toString().toBoolean()
    def run_eggnog = params.eggnog.toString().toBoolean()
    validated = VALIDATE_ASSEMBLY(samples)

    if (run_edta) {
        edta_out = EDTA_DISCOVERY(validated)
        repeats = edta_out.repeats
    } else {
        repeats = validated.map { sample, assembly, metadata -> tuple(sample, assembly, file(params.repeatmasker_library ?: params.edta_db), metadata) }
    }

    if (run_repeatmasker) {
        repeatmasker_out = REPEATMASKER(repeats)
        masked = repeatmasker_out.masked
    } else {
        masked = repeats.map { sample, assembly, repeat_library, metadata -> tuple(sample, assembly, repeat_library, metadata) }
    }
    if ((run_braker3 || run_braker4) && params.cluster_proteins) {
        clustered_proteins = CLUSTER_PROTEINS(proteins, params.cluster_min_seq_id, params.cluster_coverage)
        proteins_for_annotation = clustered_proteins.representatives
    } else {
        proteins_for_annotation = proteins
    }

    braker4_gff = null
    helixer_gff = null

    if (run_braker3) {
        braker_out = BRAKER3(masked, proteins_for_annotation)
        predictions = braker_out.predictions
        busco_out = BUSCO(predictions)
        FINAL_REPORT(busco_out.qc)
    }

    if (run_braker4) {
        braker4_out = BRAKER4(masked, proteins_for_annotation, braker4_source)
        braker4_gff = braker4_out.predictions.map { sample, gff, prediction_proteins, cds, metadata -> tuple(sample, gff) }
    }

    if (run_helixer) {
        helixer_out = HELIXER(masked, helixer_model)
        helixer_gff = helixer_out.predictions.map { sample, gff, metadata -> tuple(sample, gff) }
    }

    if (run_braker4 && run_helixer) {
        braker4_busco = braker4_out.predictions.map { sample, gff, prediction_proteins, cds, metadata ->
            tuple(sample, 'braker4', prediction_proteins, metadata)
        }
        helixer_with_genome = helixer_out.predictions.map { sample, gff, metadata -> tuple(sample, gff, metadata) }
            .join(masked.map { sample, genome, repeat_library, metadata -> tuple(sample, genome) })
            .map { sample, gff, metadata, genome -> tuple(sample, gff, genome, metadata) }
        helixer_proteins = HELIXER_PROTEINS(helixer_with_genome)
            .proteins
            .map { sample, protein_fasta, metadata -> tuple(sample, 'helixer', protein_fasta, metadata) }
        annotation_busco = ANNOTATION_BUSCO(braker4_busco.mix(helixer_proteins))

        braker4_summary = annotation_busco.qc
            .filter { sample, method, busco_dir, summary -> method == 'braker4' }
            .map { sample, method, busco_dir, summary -> tuple(sample, summary) }
        helixer_summary = annotation_busco.qc
            .filter { sample, method, busco_dir, summary -> method == 'helixer' }
            .map { sample, method, busco_dir, summary -> tuple(sample, summary) }
        busco_comparison = ANNOTATION_BUSCO_COMPARISON(braker4_summary.join(helixer_summary))
    }

    if (run_braker4 && run_helixer) {
        comparison = braker4_gff.join(helixer_gff)
        annotation_comparison = ANNOTATION_COMPARISON(comparison)

        report_inputs = braker4_summary
            .join(helixer_summary)
            .join(annotation_comparison.comparison)
            .join(busco_comparison.comparison)
            .join(braker4_out.predictions.map { sample, gff, prediction_proteins, cds, metadata -> tuple(sample, metadata) })
            .map { sample, braker4_summary_file, helixer_summary_file, annotation_comparison_file, busco_comparison_file, assembly_stats_file ->
                tuple(sample, braker4_summary_file, helixer_summary_file, annotation_comparison_file, busco_comparison_file, assembly_stats_file)
            }
        ANNOTATION_FINAL_REPORT(report_inputs)
    }

    TOOLS_REPORT()

    if (run_braker3 && run_interproscan) {
        INTERPROSCAN(predictions)
    }
    if (run_braker3 && run_eggnog) {
        EGGNOG(predictions)
    }
}
