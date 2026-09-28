process BRAKER4 {
    tag "$sample"
    label 'braker4'
    container params.braker4_runner_container
    publishDir "${params.outdir}/gene_prediction/braker4", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(genome), path(repeat_library), path(metadata)
    path proteins
    path braker4_source

    output:
    tuple val(sample), path("${sample}.braker4.gff3"), path("${sample}.braker4.proteins.faa"), path("${sample}.braker4.codingseq.fna"), path(metadata), emit: predictions
    path "${sample}.braker4.log", emit: log

    script:
    """
    test -f "${braker4_source}/Snakefile"
    test -s "$genome"
    test -s "$proteins"
    mkdir -p augustus_config

    # BRAKER4 auto-detects EP mode from the protein_fasta column. The masked
    # assembly is supplied for both genome fields because repeat masking was
    # already performed by this parent workflow.
    printf '%s\\n' 'sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf' > samples.csv
    printf '%s,%s,%s,%s,,,,,,,,,%s,\\n' \\
      "$sample" "$genome" "$genome" "$proteins" "${params.busco_lineage}" >> samples.csv

    cat > config.ini <<'EOF'
[paths]
samples_file = samples.csv
augustus_config_path = augustus_config

[containers]
braker3_image = ${params.braker4_braker3_image}

[PARAMS]
min_contig = ${params.braker4_min_contig}
skip_optimize_augustus = 1
skip_busco = 1
masking_tool = repeatmasker
no_cleanup = 0
EOF

    test -x "${params.braker4_snakemake}"
      "${params.braker4_snakemake}" \\
        --snakefile "${braker4_source}/Snakefile" \\
        --directory . \\
        --cores "$task.cpus" \\
      --jobs 1 \\
      --use-singularity \\
      --singularity-prefix "${params.braker4_singularity_prefix}" \\
      --singularity-args "${params.braker4_singularity_args}" \\
      --latency-wait 120 \\
      --restart-times 3 \\
      --rerun-incomplete \\
      --keep-going 2>&1 | tee "${sample}.braker4.log"

    out="output/${sample}/results"
    test -s "\$out/braker.gff3.gz"
    test -s "\$out/braker.aa.gz"
    test -s "\$out/braker.codingseq.gz"
    gzip -cd "\$out/braker.gff3.gz" > "${sample}.braker4.gff3"
    gzip -cd "\$out/braker.aa.gz" > "${sample}.braker4.proteins.faa"
    gzip -cd "\$out/braker.codingseq.gz" > "${sample}.braker4.codingseq.fna"
    """

    stub:
    """
    printf '##gff-version 3\\n' > "${sample}.braker4.gff3"
    printf '>stub_protein\\nM\\n' > "${sample}.braker4.proteins.faa"
    printf '>stub_cds\\nATG\\n' > "${sample}.braker4.codingseq.fna"
    printf 'BRAKER4 stub\\n' > "${sample}.braker4.log"
    """
}
