process EGGNOG {
    tag "$sample"
    label 'functional'
    container params.eggnog_container
    publishDir "${params.outdir}/functional/eggnog", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(gff), path(proteins), path(cds), path(metadata)

    output:
    tuple val(sample), path("${sample}.eggnog.emapper.annotations"), path(metadata), emit: annotations

    script:
    """
    emapper.py -i "$proteins" --itype proteins --output "${sample}.eggnog" --output_dir . \\
      --data_dir "${params.eggnog_db}" --cpu $task.cpus --override
    test -s "${sample}.eggnog.emapper.annotations"
    """

    stub:
    """
    touch "${sample}.eggnog.emapper.annotations"
    """
}
