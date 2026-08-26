process EDTA_DISCOVERY {
    tag "$sample"
    label 'repeat_discovery'
    container params.edta_container
    publishDir "${params.outdir}/repeats/edta", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(assembly), path(metadata)

    output:
    tuple val(sample), path("${sample}.edta.fa"), path("${sample}.EDTA.TElib.fa"), path(metadata), emit: repeats
    path "${sample}.EDTA.anno.gff3", emit: edta_gff
    path "${sample}.edta_id_map.tsv", emit: id_map

    script:
    """
    awk -v map="${sample}.edta_id_map.tsv" '
      BEGIN { print "normalized_id\toriginal_id" > map }
      /^>/ {
        header=substr(\$0, 2)
        original=header
        sub(/[[:space:]].*\$/, "", original)
        n++
        normalized=sprintf("ctg%09d", n)
        print ">" normalized
        print normalized "\t" header >> map
        next
      }
      { print }
    ' "$assembly" > "${sample}.edta.fa"
    cp "${sample}.edta.fa" "${sample}.fa"
    EDTA.pl --genome "${sample}.fa" --species "${params.edta_species}" --step all --anno 1 --threads $task.cpus --overwrite 1
    lib=\$(find . -type f -name '*.mod.EDTA.TElib.fa' -o -type f -name '*.EDTA.TElib.fa' | head -1)
    test -n "\$lib" && test -s "\$lib"
    cp "\$lib" "${sample}.EDTA.TElib.fa"
    gff=\$(find . -type f -name '*.EDTA.anno.gff3' -o -type f -name '*.EDTA.anno.gff' | head -1 || true)
    gff_target="${sample}.EDTA.anno.gff3"
    if [[ -n "\$gff" ]]; then cp "\$gff" "\$gff_target"; else touch "\$gff_target"; fi
    """

    stub:
    """
    awk -v map="${sample}.edta_id_map.tsv" '
      BEGIN { print "normalized_id\toriginal_id" > map }
      /^>/ {
        header=substr(\$0, 2)
        original=header
        sub(/[[:space:]].*\$/, "", original)
        n++
        normalized=sprintf("ctg%09d", n)
        print ">" normalized
        print normalized "\t" header >> map
        next
      }
      { print }
    ' "$assembly" > "${sample}.edta.fa"
    cp "${sample}.edta.fa" "${sample}.EDTA.TElib.fa"
    touch "${sample}.EDTA.anno.gff3"
    """
}
