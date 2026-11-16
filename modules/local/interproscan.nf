// Protein domains and GO terms with InterProScan, for one protein per gene of each method.
// EBI's official image keeps InterProScan in /opt/interproscan and reads its data from
// /opt/interproscan/data, so the prepared data folder (--interproscan_db) is bound there.
// -dp disables the precalculated-match lookup, so nothing is sent to EBI during a run.
// No -pa (pathways): it copies every MetaCyc/Reactome pathway ever linked to a match's InterPro
// family, in any organism (551 entries on one Col-0 line), which made the output ~10 GB per
// proteome while saying little about the gene. No XML: it repeats the TSV and GFF3.
process INTERPROSCAN {
    tag "${sample}:${method}"
    label 'functional'
    container params.interproscan_container
    containerOptions params.interproscan_db ? "-B ${file(params.interproscan_db.toString()).toAbsolutePath()}:/opt/interproscan/data" : ''
    publishDir "${params.outdir}/functional/interproscan", mode: params.publish_dir_mode

    input:
    tuple val(sample), val(method), path(proteins)

    output:
    tuple val(sample), val(method), path("${sample}.${method}.interpro.tsv"), path("${sample}.${method}.interpro.go.tsv"), emit: annotations
    path "${sample}.${method}.interpro.gff3", emit: gff
    path "${sample}.${method}.interpro.log", emit: log

    script:
    def base = "${sample}.${method}.interpro"
    """
    test -s "$proteins"
    test -d /opt/interproscan/data/pfam || { echo "InterProScan data not found at /opt/interproscan/data" >&2; exit 1; }
    mkdir -p ipr_tmp
    /opt/interproscan/interproscan.sh -i "$proteins" -f TSV,GFF3 -goterms -dp \\
      -cpu ${task.cpus} -b "${base}" -T "\$PWD/ipr_tmp" 2>&1 | tee "${base}.log"
    test -f "${base}.tsv"
    test -s "${base}.gff3"

    # GO terms per protein, merged across member databases (column 14 of the TSV).
    awk -F '\\t' -v OFS='\\t' '
      NF >= 14 && \$14 != "" && \$14 != "-" {
        n = split(\$14, terms, "|")
        for (i = 1; i <= n; i++) { t = terms[i]; sub(/\\(.*\$/, "", t); if (t != "") go[\$1, t] = 1; ids[\$1] = 1 }
      }
      END {
        print "protein_id", "go_terms"
        for (key in go) { split(key, parts, SUBSEP); list[parts[1]] = list[parts[1]] (list[parts[1]] ? "," : "") parts[2] }
        for (id in list) print id, list[id]
      }
    ' "${base}.tsv" > "${base}.go.tsv"
    """

    stub:
    def base = "${sample}.${method}.interpro"
    """
    touch "${base}.tsv" "${base}.gff3" "${base}.log"
    printf 'protein_id\\tgo_terms\\n' > "${base}.go.tsv"
    """
}
