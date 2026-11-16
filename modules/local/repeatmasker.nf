process REPEATMASKER {
    tag "$sample"
    label 'repeat_masking'
    container params.repeatmasker_container
    // With EDTA, sequences carry short names here, so RESTORE_REPEAT_NAMES publishes the
    // masked genome and the RepeatMasker table with the assembly's names instead.
    publishDir "${params.outdir}/repeats/masked", mode: params.publish_dir_mode,
        saveAs: { fn -> params.run_edta.toString().toBoolean() && (fn.endsWith('.softmasked.fa') || fn.endsWith('.repeatmasker.out')) ? null : fn }

    input:
    tuple val(sample), path(assembly), path(repeat_library), path(metadata)

    output:
    tuple val(sample), path("${sample}.softmasked.fa"), path(repeat_library), path(metadata), emit: masked
    path "${sample}.repeatmasker.out", emit: report

    script:
    """
    test -s "$repeat_library"
    cp "$assembly" "${sample}.fa"
    RepeatMasker -pa $task.cpus -xsmall -gff -lib "$repeat_library" -dir . "${sample}.fa"
    masked="${sample}.fa.masked"
    if [[ ! -s "\$masked" ]]; then
        masked="${sample}.fa.softmasked"
    fi
    test -s "\$masked"
    cp "\$masked" "${sample}.softmasked.fa"
    test -s "${sample}.fa.out"
    cp "${sample}.fa.out" "${sample}.repeatmasker.out"
    """

    stub:
    """
    cp "$assembly" "${sample}.softmasked.fa"
    touch "${sample}.repeatmasker.out"
    """
}
