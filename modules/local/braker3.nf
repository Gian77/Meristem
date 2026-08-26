process BRAKER3 {
    tag "$sample"
    label 'braker'
    container params.braker_container
    publishDir "${params.outdir}/gene_prediction", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(genome), path(repeat_library), path(metadata)
    path proteins

    output:
    tuple val(sample), path("${sample}.braker.gff3"), path("${sample}.braker.proteins.faa"), path("${sample}.braker.codingseq.fna"), path(metadata), emit: predictions
    path "${sample}.braker.log", emit: log

    script:
    """
    mkdir -p braker_work augustus_config .home
    export HOME="\$PWD/.home"
    cp -r /opt/Augustus/config/. augustus_config/
    chmod -R u+rwX augustus_config
    # Condor may evict and requeue this long-running task.  BRAKER leaves its
    # species directory behind in that case, so allow the requeued attempt to
    # reuse it instead of failing before prediction starts.
    braker.pl --genome "$genome" --prot_seq "$proteins" --softmasking --gff3 --useexisting \\
      --species "meristem_${sample}" --threads $task.cpus --workingdir braker_work \\
      --AUGUSTUS_CONFIG_PATH "\$PWD/augustus_config" \\
      > "${sample}.braker.log" 2>&1
    gff=\$(find braker_work -maxdepth 2 -type f -name 'braker.gff3' | head -1)
    prot=\$(find braker_work -maxdepth 2 -type f \( -name 'augustus.hints.aa' -o -name 'braker.aa' \) | head -1)
    cds=\$(find braker_work -maxdepth 2 -type f \( -name 'augustus.hints.codingseq' -o -name 'braker.codingseq' \) | head -1)
    test -n "\$gff" && test -s "\$gff"
    test -n "\$prot" && test -s "\$prot"
    test -n "\$cds" && test -s "\$cds"
    cp "\$gff" "${sample}.braker.gff3"
    cp "\$prot" "${sample}.braker.proteins.faa"
    cp "\$cds" "${sample}.braker.codingseq.fna"
    """

    stub:
    """
    printf '##gff-version 3\\n' > "${sample}.braker.gff3"
    printf '>stub_protein\\nM\\n' > "${sample}.braker.proteins.faa"
    printf '>stub_cds\\nATG\\n' > "${sample}.braker.codingseq.fna"
    touch "${sample}.braker.log"
    """
}
