nextflow.enable.dsl=2

include { VALIDATE_ASSEMBLY } from '../modules/local/validate_assembly'
include { EDTA_DISCOVERY } from '../modules/local/edta_discovery'
include { REPEATMASKER } from '../modules/local/repeatmasker'
include { BRAKER4 } from '../modules/local/braker4'
include { HELIXER } from '../modules/local/helixer'
include { ANNOTATION_COMPARISON } from '../modules/local/annotation_comparison'
include { HELIXER_PROTEINS; ANNOTATION_BUSCO } from '../modules/local/annotation_busco'
include { CLUSTER_PROTEINS } from '../modules/local/cluster_proteins'
include { INTERPROSCAN } from '../modules/local/interproscan'
include { EGGNOG } from '../modules/local/eggnog'
include { GENE_ORIGIN_SEARCH; GENE_ORIGIN_REPORT } from '../modules/local/gene_origin'
include { ANNOTATION_BUSCO_COMPARISON; ANNOTATION_FINAL_REPORT; ANNOTATION_METHOD_REPORT; TOOLS_REPORT } from '../modules/local/reports'

workflow ANNOTATION {
    take:
    samples
    proteins
    braker4_source
    braker4_snakemake
    helixer_model
    busco_db
    report_renderer
    gene_origin_db
    gene_origin_taxonomy
    gene_origin_script

    main:
    def run_braker4 = params.run_braker4.toString().toBoolean()
    def run_helixer = params.run_helixer.toString().toBoolean()
    def run_edta = params.run_edta.toString().toBoolean()
    def run_repeatmasker = params.run_repeatmasker.toString().toBoolean()
    def run_interproscan = params.interproscan.toString().toBoolean()
    def run_eggnog = params.eggnog.toString().toBoolean()
    tools_report = TOOLS_REPORT(report_renderer)
    tool_guides = tools_report.markdown.mix(tools_report.html, tools_report.docx, tools_report.pdf).collect()
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
    if (run_braker4 && params.cluster_proteins) {
        clustered_proteins = CLUSTER_PROTEINS(proteins, params.cluster_min_seq_id, params.cluster_coverage)
        proteins_for_annotation = clustered_proteins.representatives
    } else {
        proteins_for_annotation = proteins
    }

    braker4_gff = null
    helixer_gff = null
    method_reports = Channel.empty()

    if (run_braker4) {
        braker4_out = BRAKER4(masked, proteins_for_annotation, braker4_source, braker4_snakemake)
        braker4_gff = braker4_out.predictions.map { sample, gff, prediction_proteins, cds, metadata -> tuple(sample, gff) }
        method_reports = method_reports.mix(
            braker4_out.predictions.map { sample, gff, prediction_proteins, cds, metadata ->
                tuple(sample, 'braker4', gff, metadata)
            }
        )
    }

    if (run_helixer) {
        helixer_out = HELIXER(masked, helixer_model)
        helixer_gff = helixer_out.predictions.map { sample, gff, metadata -> tuple(sample, gff) }
        method_reports = method_reports.mix(
            helixer_out.predictions.map { sample, gff, metadata ->
                tuple(sample, 'helixer', gff, metadata)
            }
        )

        helixer_with_genome = helixer_out.predictions.map { sample, gff, metadata -> tuple(sample, gff, metadata) }
            .join(masked.map { sample, genome, repeat_library, metadata -> tuple(sample, genome) })
            .map { sample, gff, metadata, genome -> tuple(sample, gff, genome, metadata) }
        helixer_proteins = HELIXER_PROTEINS(helixer_with_genome)
            .proteins
            .map { sample, protein_fasta, coding_fasta, metadata -> tuple(sample, 'helixer', protein_fasta, metadata) }
    }

    ANNOTATION_METHOD_REPORT(method_reports, report_renderer)

    if (params.busco_db) {
        busco_inputs = Channel.empty()
        if (run_braker4) {
            busco_inputs = busco_inputs.mix(
                braker4_out.predictions.map { sample, gff, prediction_proteins, cds, metadata ->
                    tuple(sample, 'braker4', prediction_proteins, metadata)
                }
            )
        }
        if (run_helixer) {
            busco_inputs = busco_inputs.mix(helixer_proteins)
        }
        annotation_busco = ANNOTATION_BUSCO(busco_inputs, busco_db)
    }

    if (run_braker4 && run_helixer) {
        comparison = braker4_gff.join(helixer_gff)
        annotation_comparison = ANNOTATION_COMPARISON(comparison)
    }

    if (run_braker4 && run_helixer && params.busco_db) {
        braker4_summary = annotation_busco.qc
            .filter { sample, method, busco_dir, summary -> method == 'braker4' }
            .map { sample, method, busco_dir, summary -> tuple(sample, summary) }
        helixer_summary = annotation_busco.qc
            .filter { sample, method, busco_dir, summary -> method == 'helixer' }
            .map { sample, method, busco_dir, summary -> tuple(sample, summary) }
        busco_comparison = ANNOTATION_BUSCO_COMPARISON(braker4_summary.join(helixer_summary))
        report_inputs = braker4_summary
            .join(helixer_summary)
            .join(annotation_comparison.comparison)
            .join(busco_comparison.comparison)
            .join(braker4_out.predictions.map { sample, gff, prediction_proteins, cds, metadata -> tuple(sample, metadata) })
            .map { sample, braker4_summary_file, helixer_summary_file, annotation_comparison_file, busco_comparison_file, assembly_stats_file ->
                tuple(sample, braker4_summary_file, helixer_summary_file, annotation_comparison_file, busco_comparison_file, assembly_stats_file)
            }
        ANNOTATION_FINAL_REPORT(report_inputs, report_renderer, tool_guides)
    }

    // Report-only check of where predicted genes come from. Runs after gene prediction
    // because it needs the genes, and never changes the annotation or the assembly.
    if (params.gene_origin_db) {
        gene_origin_inputs = Channel.empty()
        if (run_braker4) {
            gene_origin_inputs = gene_origin_inputs.mix(
                braker4_out.predictions.map { sample, gff, prediction_proteins, cds, metadata ->
                    tuple(sample, 'braker4', gff, prediction_proteins)
                }
            )
        }
        if (run_helixer) {
            gene_origin_inputs = gene_origin_inputs.mix(
                helixer_gff
                    .join(helixer_proteins.map { sample, method, protein_fasta, metadata -> tuple(sample, protein_fasta) })
                    .map { sample, gff, protein_fasta -> tuple(sample, 'helixer', gff, protein_fasta) }
            )
        }
        gene_origin_id_maps = run_edta
            ? edta_out.id_map.map { id_map -> tuple(id_map.name - '.edta_id_map.tsv', id_map) }
            : samples.map { sample, assembly -> tuple(sample, []) }
        gene_origin_genomes = masked
            .map { sample, genome, repeat_library, metadata -> tuple(sample, genome) }
            .join(gene_origin_id_maps)

        gene_origin_search = GENE_ORIGIN_SEARCH(
            gene_origin_inputs.map { sample, method, gff, protein_fasta -> tuple(sample, method, protein_fasta) },
            gene_origin_db
        )
        GENE_ORIGIN_REPORT(
            gene_origin_inputs
                .map { sample, method, gff, protein_fasta -> tuple(sample, method, gff) }
                .join(gene_origin_search.hits, by: [0, 1])
                .combine(gene_origin_genomes, by: 0),
            gene_origin_taxonomy,
            gene_origin_script
        )
    }

    if (run_braker4 && run_interproscan) {
        INTERPROSCAN(braker4_out.predictions)
    }
    if (run_braker4 && run_eggnog) {
        EGGNOG(braker4_out.predictions)
    }
}
