process GENE_ORIGIN_SEARCH {
    tag "${sample}:${method}"
    label 'protein_clustering'
    container params.mmseqs_container

    input:
    tuple val(sample), val(method), path(proteins)
    path protein_db

    output:
    tuple val(sample), val(method), path("${sample}.${method}.gene_origin_hits.tsv"), emit: hits

    script:
    """
    mmseqs easy-search "$proteins" "$protein_db" "${sample}.${method}.gene_origin_hits.tsv" mmseqs_tmp \\
      --threads ${task.cpus} -e 1e-5 --max-seqs 100 --max-accept 5 \\
      --format-output query,target,pident,qcov,tcov,evalue,bits,theader
    """

    stub:
    """
    touch "${sample}.${method}.gene_origin_hits.tsv"
    """
}

process GENE_ORIGIN_REPORT {
    tag "${sample}:${method}"
    label 'low'
    container params.multiqc_container
    publishDir "${params.outdir}/gene_origin", mode: params.publish_dir_mode

    input:
    tuple val(sample), val(method), path(gff), path(hits), path(genome), path(id_map)
    path taxonomy_nodes
    path gene_origin_script

    output:
    path "${sample}.${method}.gene_origin.tsv", emit: table
    path "${sample}.${method}.gene_origin.md", emit: summary

    script:
    def id_map_arg = id_map ? "--id-map \"${id_map}\"" : ''
    """
    python3 "$gene_origin_script" --sample "$sample" --method "$method" \\
      --gff "$gff" --genome "$genome" --hits "$hits" --nodes "$taxonomy_nodes" $id_map_arg
    """

    stub:
    """
    printf 'contig\\tcall\\n' > "${sample}.${method}.gene_origin.tsv"
    printf '# Gene origin check — %s — %s\\n\\nStub output.\\n' "$sample" "$method" > "${sample}.${method}.gene_origin.md"
    """
}
