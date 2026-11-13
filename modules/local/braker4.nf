process BRAKER4 {
    tag "$sample"
    label 'braker4'
    container params.braker4_runner_container
    publishDir "${params.outdir}/gene_prediction/braker4", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(genome), path(repeat_library), path(metadata)
    path proteins
    path braker4_source
    path braker4_snakemake

    output:
    tuple val(sample), path("${sample}.braker4.gff3"), path("${sample}.braker4.proteins.faa"), path("${sample}.braker4.codingseq.fna"), path(metadata), emit: predictions
    path "${sample}.braker4.log", emit: log
    path "${sample}.braker4.compleasm.txt", emit: compleasm

    script:
    // BRAKER4 reads its BUSCO data from <busco_db>/lineages/<lineage>. Point it at
    // the local copy so it never downloads at run time, and bind that folder into
    // BRAKER4's own containers.
    def busco_root = params.busco_db ? file(params.busco_db.toString()).toAbsolutePath().toString() : ''
    def lineage = params.braker4_busco_lineage
    """
    test -f "${braker4_source}/Snakefile"
    test -s "$genome"
    test -s "$proteins"
    test -d "${busco_root}/lineages/${lineage}" || {
      echo "BRAKER4 BUSCO lineage not found: ${busco_root}/lineages/${lineage}" >&2
      exit 1
    }
    mkdir -p augustus_config .home .cache tmp
    export HOME="\$PWD/.home"
    export XDG_CACHE_HOME="\$PWD/.cache"
    export TMPDIR="\$PWD/tmp"

    # BRAKER4 auto-detects EP mode from the protein_fasta column. The masked
    # assembly is supplied for both genome fields because repeat masking was
    # already performed by this parent workflow. BRAKER4 uses its BUSCO lineage
    # internally (compleasm hints and best_by_compleasm); it is odb12 because
    # BRAKER4 is built for odb12. Meristem's reported BUSCO uses --busco_lineage.
    printf '%s\\n' 'sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf' > samples.csv
    printf '%s,%s,%s,%s,,,,,,,,,%s,\\n' \\
      "$sample" "$genome" "$genome" "$proteins" "${lineage}" >> samples.csv

    cat > config.ini <<'EOF'
[paths]
samples_file = samples.csv
augustus_config_path = augustus_config
busco_download_path = ${busco_root}

[containers]
braker3_image = ${params.braker4_braker3_image}

[PARAMS]
min_contig = ${params.braker4_min_contig}
skip_optimize_augustus = 1
skip_busco = 1
masking_tool = repeatmasker
no_cleanup = 0
EOF

    test -x "${braker4_snakemake}"
    "./${braker4_snakemake}" \\
      --snakefile "${braker4_source}/Snakefile" \\
      --directory . \\
      --cores "$task.cpus" \\
      --jobs 1 \\
      --use-singularity \\
      --singularity-prefix "${params.braker4_singularity_prefix}" \\
      --singularity-args "${params.braker4_singularity_args} -B ${busco_root}" \\
      --latency-wait 120 \\
      --restart-times 3 \\
      --rerun-incomplete \\
      --keep-going 2>&1 | tee "${sample}.braker4.log"

    out="output/${sample}/results"
    test -s "\$out/braker.gff3.gz"
    test -s "\$out/braker.aa.gz"
    test -s "\$out/braker.codingseq.gz"

    # BRAKER4 lets compleasm fail without stopping (`|| true`) and then finishes
    # with fewer hints and no BUSCO rescue. Stop here instead, so a broken
    # compleasm never passes as a normal BRAKER4 result.
    comp_log="logs/${sample}/compleasm/compleasm.log"
    comp_summary="\$out/quality_control/compleasm_genome_out/summary.txt"
    bbc_log="logs/${sample}/best_by_compleasm/best_by_compleasm.log"
    n_hints=\$(sed -n 's/^\\[INFO\\] Generated \\([0-9]*\\) BUSCO-based hints.*/\\1/p' "\$comp_log" | tail -n 1)
    if [[ -z "\$n_hints" || "\$n_hints" -eq 0 ]]; then
      echo "BRAKER4 compleasm produced no BUSCO hints; see \$comp_log" >&2
      exit 1
    fi
    if ! grep -q '^S:' "\$comp_summary"; then
      echo "BRAKER4 compleasm wrote no BUSCO summary; see \$comp_summary" >&2
      exit 1
    fi
    if grep -q 'script may have failed' "\$bbc_log" 2>/dev/null; then
      echo "BRAKER4 best_by_compleasm failed; see \$bbc_log" >&2
      exit 1
    fi
    {
      echo "# BRAKER4 internal compleasm check (genome mode, lineage ${lineage})"
      echo "# BUSCO-based hints: \$n_hints"
      cat "\$comp_summary"
    } > "${sample}.braker4.compleasm.txt"

    # BRAKER numbers genes from g1 in every genome, so add the sample name to keep
    # IDs unique across genomes, consistently in the GFF3 and both FASTA files.
    gzip -cd "\$out/braker.gff3.gz" | awk -F '\\t' -v OFS='\\t' -v sample="$sample" '
      /^#/ || NF < 9 { print; next }
      {
        n = split(\$9, attrs, ";")
        for (i = 1; i <= n; i++) {
          eq = index(attrs[i], "=")
          key = substr(attrs[i], 1, eq - 1)
          if (key == "ID" || key == "Parent" || key == "gene_id" || key == "transcript_id") {
            m = split(substr(attrs[i], eq + 1), values, ",")
            value = sample "_" values[1]
            for (j = 2; j <= m; j++) value = value "," sample "_" values[j]
            attrs[i] = key "=" value
          }
        }
        \$9 = attrs[1]
        for (i = 2; i <= n; i++) \$9 = \$9 ";" attrs[i]
        print
      }
    ' > "${sample}.braker4.gff3"
    gzip -cd "\$out/braker.aa.gz" | sed "s/^>/>${sample}_/" > "${sample}.braker4.proteins.faa"
    gzip -cd "\$out/braker.codingseq.gz" | sed "s/^>/>${sample}_/" > "${sample}.braker4.codingseq.fna"
    """

    stub:
    """
    printf '##gff-version 3\\n' > "${sample}.braker4.gff3"
    printf '>stub_protein\\nM\\n' > "${sample}.braker4.proteins.faa"
    printf '>stub_cds\\nATG\\n' > "${sample}.braker4.codingseq.fna"
    printf 'BRAKER4 stub\\n' > "${sample}.braker4.log"
    printf 'S:0.00%%, 0\\n' > "${sample}.braker4.compleasm.txt"
    """
}
