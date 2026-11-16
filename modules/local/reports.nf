process ANNOTATION_BUSCO_COMPARISON {
    tag "$sample"
    label 'low'
    publishDir "${params.outdir}/busco_comparison", mode: params.publish_dir_mode

    input:
    tuple val(sample),
          path(braker4_summary, stageAs: 'braker4_short_summary.txt'),
          path(helixer_summary, stageAs: 'helixer_short_summary.txt')

    output:
    tuple val(sample), path("${sample}.busco_annotation_comparison.tsv"), emit: comparison

    script:
    """
    metric() {
        awk -v key="\$1" '
            /C:[0-9.]+%/ {
                if (key == "complete") match(\$0, /C:([0-9.]+)%/, a)
                else if (key == "single") match(\$0, /S:([0-9.]+)%/, a)
                else if (key == "duplicated") match(\$0, /D:([0-9.]+)%/, a)
                else if (key == "fragmented") match(\$0, /F:([0-9.]+)%/, a)
                else if (key == "missing") match(\$0, /M:([0-9.]+)%/, a)
                else if (key == "total") match(\$0, /n[:=]([0-9]+)/, a)
                if (a[1] != "") { print a[1]; exit }
            }
        ' "\$2"
    }
    {
      printf 'sample\\tmethod\\tcomplete_pct\\tsingle_copy_pct\\tduplicated_pct\\tfragmented_pct\\tmissing_pct\\ttotal_buscos\\n'
      for method in braker4 helixer; do
        summary="\${method}_summary.txt"
        if [[ "\$method" == braker4 ]]; then cp "$braker4_summary" "\$summary"; else cp "$helixer_summary" "\$summary"; fi
        printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
          "$sample" "\$method" "\$(metric complete "\$summary")" "\$(metric single "\$summary")" \\
          "\$(metric duplicated "\$summary")" "\$(metric fragmented "\$summary")" \\
          "\$(metric missing "\$summary")" "\$(metric total "\$summary")"
      done
    } > "${sample}.busco_annotation_comparison.tsv"
    """

    stub:
    """
    printf 'sample\\tmethod\\tcomplete_pct\\tsingle_copy_pct\\tduplicated_pct\\tfragmented_pct\\tmissing_pct\\ttotal_buscos\\n%s\\tbraker4\\t0.0\\t0.0\\t0.0\\t0.0\\t100.0\\t1\\n%s\\thelixer\\t0.0\\t0.0\\t0.0\\t0.0\\t100.0\\t1\\n' \\
      "$sample" "$sample" > "${sample}.busco_annotation_comparison.tsv"
    """
}

process ANNOTATION_METHOD_REPORT {
    tag "${sample}:${method}"
    label 'low'
    publishDir { "${params.outdir}/reports/${sample}" }, mode: params.publish_dir_mode

    input:
    tuple val(sample),
          val(method),
          path(annotation_gff, stageAs: 'annotation.gff3'),
          path(assembly_stats, stageAs: 'assembly_stats.tsv')
    path report_renderer

    output:
    path "${sample}.${method}.annotation_summary.md", emit: summary
    path "${sample}.${method}.annotation_manifest.tsv", emit: manifest
    path "${sample}.${method}.analysis_report.html"

    script:
    """
    test -s annotation.gff3

    count_features() {
      awk -F '\\t' -v feature="\$1" '\$0 !~ /^#/ && \$3 == feature { n++ } END { print n + 0 }' annotation.gff3
    }

    case "${method}" in
      braker4)
        gff_output="${params.outdir}/gene_prediction/braker4/${sample}.braker4.gff3"
        protein_output="${params.outdir}/gene_prediction/braker4/${sample}.braker4.proteins.faa"
        coding_output="${params.outdir}/gene_prediction/braker4/${sample}.braker4.codingseq.fna"
        ;;
      helixer)
        gff_output="${params.outdir}/gene_prediction/helixer/${sample}.helixer.gff3"
        protein_output="${params.outdir}/gene_prediction/helixer/${sample}.helixer.proteins.faa"
        coding_output="${params.outdir}/gene_prediction/helixer/${sample}.helixer.codingseq.fna"
        ;;
      *)
        echo "Unsupported annotation method: ${method}" >&2
        exit 1
        ;;
    esac

    {
      printf 'sample\\tmethod\\tstatus\\tgff3\\tprotein_fasta\\tcoding_sequence\\tgene_count\\tmrna_count\\tcds_count\\texon_count\\n'
      printf '%s\\t%s\\tPASS\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "${sample}" "${method}" "\${gff_output}" "\${protein_output}" "\${coding_output}" \\
        "\$(count_features gene)" "\$(count_features mRNA)" "\$(count_features CDS)" "\$(count_features exon)"
    } > "${sample}.${method}.annotation_manifest.tsv"

    {
      echo "# Annotation result — ${sample} — ${method}"
      echo
      echo "## Run status"
      echo
      echo "**PASS** — the ${method} annotation GFF3 was produced and is non-empty."
      echo
      echo "## Annotation files"
      echo
      echo '| Product | Path |'
      echo '|---|---|'
      echo "| GFF3 annotation | \${gff_output} |"
      echo "| Protein FASTA | \${protein_output} |"
      echo "| Coding-sequence FASTA | \${coding_output} |"
      echo
      echo "## Feature counts"
      echo
      echo '| Feature | Count |'
      echo '|---|---:|'
      printf '| gene | %s |\\n' "\$(count_features gene)"
      printf '| mRNA | %s |\\n' "\$(count_features mRNA)"
      printf '| CDS | %s |\\n' "\$(count_features CDS)"
      printf '| exon | %s |\\n' "\$(count_features exon)"
      echo
      echo "## Assembly QC"
      echo
      echo '```text'
      cat assembly_stats.tsv
      echo '```'
      echo
      echo "The complete Nextflow run status is in \\`${params.outdir}/pipeline_info/RUN_STATUS.md\\`."
      echo "BUSCO and dual-method comparisons, when enabled, are in the sibling \\`busco_comparison/\\` and \\`comparison/\\` directories."
    } > "${sample}.${method}.annotation_summary.md"
    python3 "${report_renderer}" --sample "${sample}" --method "${method}" \\
      --analysis-markdown "${sample}.${method}.annotation_summary.md" \\
      --assembly assembly_stats.tsv --outdir .
    """

    stub:
    """
    printf 'sample\\tmethod\\tstatus\\tgff3\\tprotein_fasta\\tcoding_sequence\\tgene_count\\tmrna_count\\tcds_count\\texon_count\\n%s\\t%s\\tPASS\\t%s\\t%s\\t%s\\t0\\t0\\t0\\t0\\n' \\
      "${sample}" "${method}" "${params.outdir}/gene_prediction/${method}/${sample}.${method}.gff3" \\
      "${params.outdir}/gene_prediction/${method}/${sample}.${method}.proteins.faa" "not available in stub" \\
      > "${sample}.${method}.annotation_manifest.tsv"
    printf '# Annotation result — %s — %s\\n\\n**PASS** — stub output produced.\\n' \\
      "${sample}" "${method}" > "${sample}.${method}.annotation_summary.md"
    python3 "${report_renderer}" --sample "${sample}" --method "${method}" \\
      --analysis-markdown "${sample}.${method}.annotation_summary.md" \\
      --assembly assembly_stats.tsv --outdir .
    """
}

process ANNOTATION_FINAL_REPORT {
    tag "$sample"
    label 'low'
    container params.multiqc_container
    publishDir { "${params.outdir}/reports/${sample}" }, mode: params.publish_dir_mode

    input:
    tuple val(sample),
          path(braker4_summary, stageAs: 'busco_braker4_summary.txt'),
          path(helixer_summary, stageAs: 'busco_helixer_summary.txt'),
          path(annotation_comparison, stageAs: 'annotation_comparison.tsv'),
          path(busco_comparison, stageAs: 'busco_annotation_comparison.tsv'),
          path(assembly_stats, stageAs: 'assembly_stats.tsv')
    path report_renderer
    path tool_guides

    output:
    path "${sample}_annotation_summary.md", emit: summary
    path "${sample}_analysis_report.html"
    path "${sample}_analysis_report.md"
    path "${sample}_annotation_interpretation.md"
    path "${sample}_annotation_interpretation.html"
    path "${sample}_annotation_interpretation.pdf"
    path "${sample}_annotation_report.docx"
    path "${sample}_annotation_report.pdf"
    path 'multiqc_report.html'
    path 'multiqc_data'
    path "${sample}_annotation_reports.zip"

    script:
    """
    # These inputs are custom summary tables rather than native MultiQC logs.
    # Keep the HTML artifact even when no bundled MultiQC module recognizes
    # them; the custom report below still contains the pipeline summaries.
    if multiqc . --force --outdir multiqc_data_out; then
      :
    fi
    if [[ -s multiqc_data_out/multiqc_report.html ]]; then
      cp multiqc_data_out/multiqc_report.html multiqc_report.html
      if [[ -d multiqc_data_out/multiqc_data ]]; then
        mv multiqc_data_out/multiqc_data multiqc_data
      else
        mkdir -p multiqc_data
      fi
    else
      mkdir -p multiqc_data
      cat > multiqc_report.html <<'HTMLEOF'
<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>Meristem report</title></head>
<body><h1>Meristem annotation report</h1>
<p>MultiQC found no native module results for these custom summary inputs.</p>
<p>See the accompanying Markdown summary and TSV files for the pipeline results.</p>
</body></html>
HTMLEOF
    fi

    braker4_line="\$(grep -E 'C:[0-9.]+%' busco_braker4_summary.txt | head -1)"
    helixer_line="\$(grep -E 'C:[0-9.]+%' busco_helixer_summary.txt | head -1)"
    {
      echo "# Annotation summary — ${sample}"
      echo
      echo "_Generated by Meristem on \$(date -u '+%Y-%m-%d %H:%M UTC')_"
      echo
      echo "## Annotation completeness (BUSCO — ${params.busco_lineage})"
      echo
      echo '| Method | Complete | Single-copy | Duplicated | Fragmented | Missing | Total |'
      echo '|---|---:|---:|---:|---:|---:|---:|'
      awk -F '\\t' 'NR > 1 { printf "| %s | %s%% | %s%% | %s%% | %s%% | %s%% | %s |\\n", \$2, \$3, \$4, \$5, \$6, \$7, \$8 }' busco_annotation_comparison.tsv
      echo
      echo "## Annotation feature counts"
      echo
      echo '```text'
      cat annotation_comparison.tsv
      echo '```'
      echo
      echo "## Assembly input statistics"
      echo
      echo '```text'
      cat assembly_stats.tsv
      echo '```'
      echo
      echo "## Interpretation"
      echo
      echo "BRAKER4 BUSCO: \$braker4_line"
      echo
      echo "Helixer BUSCO: \$helixer_line"
      echo
      echo "See ${sample}_analysis_report.html for the analysis and pipeline_tools.html for the tools guide. The original BUSCO summaries retain the full marker-assessment details."
      echo
      echo "## Detailed interpretation and report formats"
      echo
      echo "See ${sample}_annotation_interpretation.md for step-by-step interpretation from assembly assessment through annotation quality. DOCX and PDF versions are provided as ${sample}_annotation_report.docx and ${sample}_annotation_report.pdf."
    } > "${sample}_annotation_summary.md"

    # The bundled renderer makes the HTML/data report useful even when the
    # stock MultiQC modules do not recognize these custom TSV inputs. It also
    # writes the detailed interpretation and portable DOCX/PDF versions.
    python3 render_annotation_reports.py \
      --sample "${sample}" \
      --lineage "${params.busco_lineage}" \
      --summary "${sample}_annotation_summary.md" \
      --busco busco_annotation_comparison.tsv \
      --comparison annotation_comparison.tsv \
      --assembly assembly_stats.tsv \
      --outdir . \
      --tools-link ../pipeline_tools/pipeline_tools.html \
      --run-edta "${params.run_edta}" \
      --run-repeatmasker "${params.run_repeatmasker}" \
      --interproscan "${params.interproscan}" \
      --eggnog "${params.eggnog}"

    mkdir -p package
    cp "${sample}_annotation_summary.md" package/
    cp "${sample}_analysis_report.html" "${sample}_analysis_report.md" package/
    cp pipeline_tools.* package/
    cp "${sample}_annotation_interpretation.md" package/
    cp "${sample}_annotation_report.docx" package/
    cp "${sample}_annotation_report.pdf" package/
    cp multiqc_report.html package/
    cp -a multiqc_data package/
    cp annotation_comparison.tsv busco_annotation_comparison.tsv package/
    python3 -m zipfile -c "${sample}_annotation_reports.zip" package
    """

    stub:
    """
    printf '# Annotation summary — %s\\n' "$sample" > "${sample}_annotation_summary.md"
    printf '<!doctype html><html><body><h1>Analysis report — %s</h1><p>Stub output.</p></body></html>\\n' "$sample" > "${sample}_analysis_report.html"
    printf '# Analysis report — %s\\n\\nStub output.\\n' "$sample" > "${sample}_analysis_report.md"
    printf '# Interpretation — %s\\n\\n**PASS** — stub report generated.\\n' "$sample" > "${sample}_annotation_interpretation.md"
    printf '<!doctype html><html><body><h1>Interpretation — %s</h1><p>Stub output.</p></body></html>\\n' "$sample" > "${sample}_annotation_interpretation.html"
    printf 'stub interpretation PDF report\\n' > "${sample}_annotation_interpretation.pdf"
    printf 'stub DOCX report\\n' > "${sample}_annotation_report.docx"
    printf 'stub PDF report\\n' > "${sample}_annotation_report.pdf"
    touch multiqc_report.html "${sample}_annotation_reports.zip"
    mkdir -p multiqc_data
    printf 'sample\\tsection\\tmetric\\tvalue\\n%s\\tstatus\\toverall\\tPASS\\n' "$sample" > multiqc_data/meristem_annotation_summary_mqc.tsv
    """
}

process TOOLS_REPORT {
    label 'low'
    container params.multiqc_container
    publishDir "${params.outdir}/reports/pipeline_tools", mode: params.publish_dir_mode

    input:
    path report_renderer

    output:
    path 'pipeline_tools.md', emit: markdown
    path 'pipeline_tools.html', emit: html
    path 'pipeline_tools.docx', emit: docx
    path 'pipeline_tools.pdf', emit: pdf

    script:
    // Only reporting choices and image references are serialized. Do not dump
    // arbitrary params/environment values into a distributable report.
    def reportConfig = groovy.json.JsonOutput.toJson([
        run_edta: params.run_edta,
        run_repeatmasker: params.run_repeatmasker,
        run_braker4: params.run_braker4,
        run_helixer: params.run_helixer,
        cluster_proteins: params.cluster_proteins,
        cluster_min_seq_id: params.cluster_min_seq_id,
        cluster_coverage: params.cluster_coverage,
        busco_db: params.busco_db ? true : false,
        busco_lineage: params.busco_lineage,
        interproscan: params.interproscan,
        eggnog: params.eggnog,
        edta_container: params.edta_container,
        repeatmasker_container: params.repeatmasker_container,
        mmseqs_container: params.mmseqs_container,
        braker4_runner_container: params.braker4_runner_container,
        helixer_container: params.helixer_container,
        busco_container: params.busco_container,
        multiqc_container: params.multiqc_container,
        interproscan_container: params.interproscan_container,
        eggnog_container: params.eggnog_container
    ])
    """
    cat > tools_config.json <<'CONFIGEOF'
${reportConfig}
CONFIGEOF
    python3 "${report_renderer}" --tools-only --tools-config tools_config.json --outdir .
    """
}
