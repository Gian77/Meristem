process CLUSTER_PROTEINS {
    tag "protein evidence"
    label 'protein_clustering'
    container params.mmseqs_container

    input:
    path proteins
    val min_seq_id
    val coverage

    output:
    path 'clustered_proteins.faa', emit: representatives

    script:
    """
    mmseqs easy-linclust "$proteins" clustered_proteins mmseqs_tmp \\
      --min-seq-id "$min_seq_id" \\
      -c "$coverage" \\
      --cov-mode 1
    test -s clustered_proteins_rep_seq.fasta
    cp clustered_proteins_rep_seq.fasta clustered_proteins.faa
    """

    stub:
    """
    cp "$proteins" clustered_proteins.faa
    """
}
