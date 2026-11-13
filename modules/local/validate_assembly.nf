process VALIDATE_ASSEMBLY {
    tag "$sample"
    label 'low'
    publishDir "${params.outdir}/assembly_qc", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(assembly)

    output:
    tuple val(sample), path(assembly), path("${sample}.assembly.tsv"), emit: validated

    script:
    """
    test -s "$assembly"
    awk -v sample="$sample" -v minlen="${params.min_contig_size}" '
      /^>/ {
        if (seen && len >= minlen) { n++; total += len; if (len > max) max=len; a[n]=len }
        seen=1; len=0
      }
      !/^>/ { gsub(/[[:space:]]/, ""); len += length(\$0); gaps += gsub(/[Nn]/, "N") }
      END {
        if (len >= minlen) { n++; total += len; if (len > max) max=len; a[n]=len }
        if (n == 0) { print "No sequences passed minimum length" > "/dev/stderr"; exit 1 }
        asort(a)
        cum=0; n50=0
        for (i=n; i>=1; i--) { cum += a[i]; if (cum >= total/2) { n50=a[i]; break } }
        printf "sample\\t%s\\nsequences_ge_min\\t%d\\ntotal_length\\t%d\\nlargest\\t%d\\nN50\\t%d\\ngap_bases\\t%d\\n" , sample,n,total,max,n50,gaps
      }
    ' "$assembly" > "${sample}.assembly.tsv"
    """

    stub:
    """
    printf 'sample\\t%s\\nsequences_ge_min\\t1\\ntotal_length\\t1\\nlargest\\t1\\nN50\\t1\\ngap_bases\\t0\\n' "$sample" > "${sample}.assembly.tsv"
    """
}
