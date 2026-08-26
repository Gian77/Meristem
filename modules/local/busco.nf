process BUSCO {
    tag "$sample"
    label 'busco'
    container params.busco_container
    publishDir "${params.outdir}/busco", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(gff), path(proteins), path(cds), path(metadata)

    output:
    tuple val(sample), path("busco_${sample}"), path(metadata), emit: qc

    script:
    def db_args = params.busco_db ? "--download_path '${params.busco_db}' --offline" : ""
    """
    busco -i "$proteins" -l "${params.busco_lineage}" -m proteins -o "busco_${sample}" --cpu $task.cpus $db_args
    test -s "busco_${sample}/short_summary.txt" || test -s "busco_${sample}"/short_summary*.txt
    """

    stub:
    """
    mkdir -p "busco_${sample}"
    printf 'C:0.0%%[S:0.0%%,D:0.0%%],F:0.0%%,M:100.0%%,n:1\\n' > "busco_${sample}/short_summary.txt"
    """
}
