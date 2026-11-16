process ANNOTATION_COMPARISON {
    tag "$sample"
    label 'low'
    publishDir "${params.outdir}/comparison", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(braker4_gff), path(helixer_gff)

    output:
    tuple val(sample), path("${sample}.annotation_comparison.tsv"), emit: comparison

    script:
    """
    count_features() {
      awk -v feature="\$1" '\$3 == feature { n++ } END { print n + 0 }' "\$2"
    }
    {
      printf 'sample\\tmetric\\tbraker4\\thelixer\\n'
      for feature in gene mRNA transcript CDS exon; do
        printf '%s\\t%s\\t%s\\t%s\\n' \\
          "$sample" "\$feature" "\$(count_features \"\$feature\" \"$braker4_gff\")" "\$(count_features \"\$feature\" \"$helixer_gff\")"
      done
    } > "${sample}.annotation_comparison.tsv"
    """

    stub:
    """
    printf 'sample\\tmetric\\tbraker4\\thelixer\\n%s\\tgene\\t0\\t0\\n' "$sample" > "${sample}.annotation_comparison.tsv"
    """
}
