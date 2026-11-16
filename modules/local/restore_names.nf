// Publish annotations and repeat files with the assembly's own sequence names.
// Every step before this uses EDTA's short ctg000000001 names; these processes
// translate them back with EDTA's id map, and the translated files are the ones
// published. Without EDTA (--run_edta false) no renaming happened, so the files
// are copied unchanged.

// Shell command that writes dst from src with original names, or copies it when there is no map.
def renameCommand(script, idMap, sample, kind, src, dst) {
    return idMap
        ? "python3 \"${script}\" rename --id-map \"${idMap}\" --sample \"${sample}\" --kind ${kind} \"${src}\" \"${dst}\""
        : "cp \"${src}\" \"${dst}\""
}

process RESTORE_PREDICTION_NAMES {
    tag "${sample}:${method}"
    label 'low'
    container params.multiqc_container
    // Files are named <sample>.<method>.*, so each goes to gene_prediction/<method>/
    publishDir "${params.outdir}/gene_prediction", mode: params.publish_dir_mode,
        saveAs: { fn -> (fn =~ /\.(braker4|helixer)\./) ? "${(fn =~ /\.(braker4|helixer)\./)[0][1]}/${fn}" : fn }

    input:
    tuple val(sample), val(method), path(gff, stageAs: 'input/*'), path(proteins, stageAs: 'input/*'), path(cds, stageAs: 'input/*'), path(id_map)
    path restore_script

    output:
    tuple val(sample), val(method), path("${sample}.${method}.gff3"), path("${sample}.${method}.proteins.faa"), path("${sample}.${method}.codingseq.fna"), emit: renamed
    tuple val(sample), val(method), path("${sample}.${method}.primary_proteins.faa"), emit: primary

    script:
    def out = "${sample}.${method}"
    def gff_cmd = renameCommand(restore_script, id_map, sample, 'gff', gff, "${out}.gff3")
    def protein_cmd = renameCommand(restore_script, id_map, sample, 'fasta', proteins, "${out}.proteins.faa")
    def cds_cmd = renameCommand(restore_script, id_map, sample, 'fasta', cds, "${out}.codingseq.fna")
    """
    ${gff_cmd}
    ${protein_cmd}
    ${cds_cmd}
    # One protein per gene (longest isoform), for functional annotation and for users
    python3 "$restore_script" primary --gff "${out}.gff3" --proteins "${out}.proteins.faa" "${out}.primary_proteins.faa"
    test -s "${out}.primary_proteins.faa"
    """

    stub:
    def out = "${sample}.${method}"
    """
    printf '##gff-version 3\\n' > "${out}.gff3"
    printf '>stub\\nM\\n' > "${out}.proteins.faa"
    printf '>stub\\nATG\\n' > "${out}.codingseq.fna"
    printf '>stub\\nM\\n' > "${out}.primary_proteins.faa"
    """
}

process RESTORE_REPEAT_NAMES {
    tag "$sample"
    label 'low'
    container params.multiqc_container
    publishDir "${params.outdir}/repeats", mode: params.publish_dir_mode,
        saveAs: { fn -> fn.endsWith('.EDTA.TEanno.gff3') ? "edta/${fn}" : "masked/${fn}" }

    input:
    tuple val(sample), path(te_gff, stageAs: 'input/*'), path(softmasked, stageAs: 'input/*'), path(rm_out, stageAs: 'input/*'), path(id_map)
    path restore_script

    output:
    path "${sample}.EDTA.TEanno.gff3", emit: te_gff
    path "${sample}.softmasked.fa", optional: true, emit: softmasked
    path "${sample}.repeatmasker.out", optional: true, emit: rm_out

    script:
    def te_cmd = renameCommand(restore_script, id_map, sample, 'gff', te_gff, "${sample}.EDTA.TEanno.gff3")
    def masked_cmd = softmasked ? renameCommand(restore_script, id_map, sample, 'fasta', softmasked, "${sample}.softmasked.fa") : ''
    def rm_cmd = rm_out ? renameCommand(restore_script, id_map, sample, 'rmout', rm_out, "${sample}.repeatmasker.out") : ''
    """
    ${te_cmd}
    ${masked_cmd}
    ${rm_cmd}
    """

    stub:
    """
    printf '##gff-version 3\\n' > "${sample}.EDTA.TEanno.gff3"
    """
}
