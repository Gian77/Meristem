process INTERPROSCAN {
    tag "$sample"
    label 'functional'
    container params.interproscan_container
    publishDir "${params.outdir}/functional/interproscan", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(gff), path(proteins), path(cds), path(metadata)

    output:
    tuple val(sample), path("${sample}.interpro.tsv"), path(metadata), emit: annotations

    script:
    """
    interproscan.sh -i "$proteins" -f tsv,xml,gff3 -goterms -pa -cpu $task.cpus \\
      --output-file-base "${sample}.interpro" --datadir "${params.interproscan_db}"
    test -s "${sample}.interpro.tsv"
    """

    stub:
    """
    touch "${sample}.interpro.tsv"
    """
}
