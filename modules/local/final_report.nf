process FINAL_REPORT {
    tag "$sample"
    label 'low'
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(busco_dir), path(metadata)

    output:
    path "${sample}.annotation_summary.tsv", emit: summary

    script:
    """
    summary=\$(find "$busco_dir" -type f -name 'short_summary*.txt' | head -1)
    test -n "\$summary"
    {
      echo -e "sample\\tbusco_lineage\\tbusco_summary"
      printf '%s\\t%s\\t' "$sample" "${params.busco_lineage}"
      tr '\\n' ' ' < "\$summary"
      echo
    } > "${sample}.annotation_summary.tsv"
    cp "$metadata" "${sample}.assembly.tsv"
    """

    stub:
    """
    printf 'sample\\tbusco_lineage\\tbusco_summary\\n%s\\t%s\\tstub\\n' "$sample" "${params.busco_lineage}" > "${sample}.annotation_summary.tsv"
    """
}
