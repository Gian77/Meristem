process HELIXER {
    tag "$sample"
    label 'helixer'
    container params.helixer_container
    containerOptions params.helixer_container_options
    publishDir "${params.outdir}/gene_prediction/helixer", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(genome), path(repeat_library), path(metadata)
    path model

    output:
    tuple val(sample), path("${sample}.helixer.gff3"), path(metadata), emit: predictions
    path "${sample}.helixer.log", emit: log

    script:
    def deterministic_arg = params.helixer_deterministic ? "--deterministic" : ""
    """
    test -s "$genome"
    Helixer.py \\
      --fasta-path "$genome" \\
      --lineage "${params.helixer_lineage}" \\
      --species "$sample" \\
      --gff-output-path "${sample}.helixer.gff3" \\
      --subsequence-length 64152 \\
      --model-filepath "$model" $deterministic_arg \\
      2>&1 | tee "${sample}.helixer.log"
    test -s "${sample}.helixer.gff3"
    """

    stub:
    """
    printf '##gff-version 3\\n' > "${sample}.helixer.gff3"
    printf 'Helixer stub\\n' > "${sample}.helixer.log"
    """
}
