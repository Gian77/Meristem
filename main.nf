#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ANNOTATION } from './workflows/annotation'

// Provenance helpers. They only read files: no git command is run.
def digestFile(path, algorithm) {
    def md = java.security.MessageDigest.getInstance(algorithm)
    path.withInputStream { stream ->
        stream.eachByte(1048576) { buffer, n -> md.update(buffer, 0, n) }
    }
    return md.digest().encodeHex().toString()
}

// Commit recorded in .git/HEAD, read as a file. It says nothing about edits
// made since that commit; PIPELINE_SOURCE.sha256 covers those.
def gitCommit(dir) {
    def head = new File(dir.toString(), '.git/HEAD')
    if (!head.exists()) {
        return 'not a git checkout'
    }
    def text = head.text.trim()
    if (!text.startsWith('ref: ')) {
        return text
    }
    def ref = text.substring(5)
    def branch = ref - 'refs/heads/'
    def refFile = new File(dir.toString(), ".git/${ref}")
    if (refFile.exists()) {
        return "${refFile.text.trim()} (branch ${branch})"
    }
    def packed = new File(dir.toString(), '.git/packed-refs')
    def line = packed.exists() ? packed.readLines().find { it.endsWith(" ${ref}") } : null
    return line ? "${line.split(' ')[0]} (branch ${branch})" : "unknown (branch ${branch})"
}

// sha256 of every file that defines the pipeline's behaviour, so a run made
// from edited, uncommitted code can be told apart from the committed version.
def pipelineSourceDigests(dir) {
    def root = new File(dir.toString())
    def files = ['main.nf', 'nextflow.config', 'nextflow_schema.json', 'requirements-snakemake.txt']
        .collect { new File(root, it) }
        .findAll { it.isFile() }
    ['conf', 'modules', 'workflows', 'scripts', 'bin'].each { sub ->
        def folder = new File(root, sub)
        if (folder.isDirectory()) {
            folder.eachFileRecurse { f ->
                if (f.isFile() && !f.path.contains('__pycache__') && !f.name.endsWith('.pyc')) {
                    files << f
                }
            }
        }
    }
    return files
        .collect { f -> [root.toPath().relativize(f.toPath()).toString(), digestFile(f.toPath(), 'SHA-256')] }
        .sort { it[0] }
}

workflow {
    def run_braker4 = params.run_braker4.toString().toBoolean()
    def run_helixer = params.run_helixer.toString().toBoolean()
    def run_edta = params.run_edta.toString().toBoolean()
    def run_repeatmasker = params.run_repeatmasker.toString().toBoolean()
    def run_interproscan = params.interproscan.toString().toBoolean()
    def run_eggnog = params.eggnog.toString().toBoolean()
    if (!params.input) {
        error "Missing --input samplesheet"
    }
    if (!run_braker4 && !run_helixer) {
        error "Enable at least one annotation method: --run_braker4 or --run_helixer"
    }
    if (run_braker4 && !params.proteins) {
        error "Missing --proteins protein FASTA"
    }
    if (run_braker4 && !params.braker4_source) {
        error "--run_braker4 requires --braker4_source pointing to a BRAKER4 checkout"
    }
    if (run_braker4 && !params.braker4_snakemake) {
        error "--run_braker4 requires --braker4_snakemake pointing to the Snakemake executable"
    }
    def braker4_installed_commit = null
    if (run_braker4) {
        def marker = file("${params.braker4_source}/MERISTEM_BRAKER4_SOURCE.txt")
        def commitLine = marker.exists() ? marker.readLines().find { it.startsWith('commit\t') } : null
        braker4_installed_commit = commitLine ? commitLine.split('\t')[1] : null
        def allow_unpinned = params.braker4_allow_unpinned.toString().toBoolean()
        if (braker4_installed_commit != params.braker4_commit && !allow_unpinned) {
            error "BRAKER4 at ${params.braker4_source} is not the tested version ${params.braker4_commit} " +
                "(found: ${braker4_installed_commit ?: 'no MERISTEM_BRAKER4_SOURCE.txt'}). " +
                "Install it with scripts/setup_braker4.sh, or pass --braker4_allow_unpinned true to use it anyway."
        }
        if (!workflow.stubRun) {
            if (!params.busco_db) {
                error "--run_braker4 requires --busco_db: BRAKER4 reads <busco_db>/lineages/${params.braker4_busco_lineage} and must not download it"
            }
            def braker4_lineage = file("${params.busco_db}/lineages/${params.braker4_busco_lineage}")
            if (!braker4_lineage.isDirectory()) {
                error "BRAKER4 BUSCO lineage not found: ${braker4_lineage}. Download ${params.braker4_busco_lineage} into ${params.busco_db}/lineages/"
            }
        }
    }
    if (run_helixer && !params.helixer_model) {
        error "--run_helixer requires --helixer_model pointing to a local HDF5 model"
    }
    cluster_min_seq_id = params.cluster_min_seq_id.toString().toBigDecimal()
    cluster_coverage = params.cluster_coverage.toString().toBigDecimal()
    if (cluster_min_seq_id < 0 || cluster_min_seq_id > 1) {
        error "--cluster_min_seq_id must be between 0 and 1"
    }
    if (cluster_coverage <= 0 || cluster_coverage > 1) {
        error "--cluster_coverage must be greater than 0 and no more than 1"
    }
    if (run_interproscan && !params.interproscan_db) {
        error "--interproscan requires --interproscan_db"
    }
    if (run_eggnog && !params.eggnog_db) {
        error "--eggnog requires --eggnog_db"
    }
    if (params.gene_origin_db && !params.gene_origin_taxonomy) {
        error "--gene_origin_db requires --gene_origin_taxonomy pointing to an NCBI nodes.dmp"
    }
    if (!run_edta && !params.repeatmasker_library && !params.edta_db) {
        error "When --run_edta false, provide --repeatmasker_library or --edta_db"
    }

    def seen_samples = [] as Set
    samples = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !row.assembly) {
                error "Each samplesheet row requires sample and assembly columns"
            }
            if (!seen_samples.add(row.sample.toString())) {
                error "Duplicate sample name in samplesheet: ${row.sample}"
            }
            tuple(row.sample.toString(), file(row.assembly.toString()))
        }

    proteins = params.proteins ? file(params.proteins, checkIfExists: true) : null
    braker4_source = run_braker4 && params.braker4_source ? file(params.braker4_source, checkIfExists: true) : null
    braker4_snakemake = run_braker4 && params.braker4_snakemake ? file(params.braker4_snakemake, checkIfExists: true) : null
    helixer_model = run_helixer && params.helixer_model ? file(params.helixer_model, checkIfExists: true) : null
    busco_db = params.busco_db ? file(params.busco_db, checkIfExists: true) : []
    report_renderer = file("${projectDir}/scripts/render_annotation_reports.py", checkIfExists: true)
    gene_origin_db = params.gene_origin_db ? file(params.gene_origin_db, checkIfExists: true) : []
    gene_origin_taxonomy = params.gene_origin_db ? file(params.gene_origin_taxonomy, checkIfExists: true) : []
    gene_origin_script = file("${projectDir}/scripts/summarize_gene_origin.py", checkIfExists: true)

    ANNOTATION(samples, proteins, braker4_source, braker4_snakemake, helixer_model, busco_db, report_renderer,
        gene_origin_db, gene_origin_taxonomy, gene_origin_script)

    // Record what this run used, before any task can change the inputs.
    def source_digests = pipelineSourceDigests(projectDir)
    def source_listing = source_digests.collect { "${it[1]}  ${it[0]}" }.join('\n') + '\n'
    def source_digest = java.security.MessageDigest.getInstance('SHA-256')
        .digest(source_listing.getBytes('UTF-8')).encodeHex().toString()
    def meristem_commit = gitCommit(projectDir)
    def input_checksums = file(params.input).splitCsv(header: true).collect { row ->
        def assembly = file(row.assembly.toString())
        [row.sample.toString(), assembly.toString(), assembly.exists() ? digestFile(assembly, 'MD5') : 'missing']
    }
    def helixer_model_md5 = helixer_model ? digestFile(helixer_model, 'MD5') : 'not used'
    def braker4_marker = run_braker4 ? file("${params.braker4_source}/MERISTEM_BRAKER4_SOURCE.txt") : null
    def braker4_record = !run_braker4 ? 'not used'
        : braker4_marker.exists() ? braker4_marker.text.trim().readLines().collect { "    ${it}" }.join('\n')
        : 'no MERISTEM_BRAKER4_SOURCE.txt (unpinned BRAKER4)'

    workflow.onComplete = {
        def status = workflow.success ? 'SUCCESS' : 'FAILED'
        def conclusion = workflow.success
            ? 'The workflow completed without failed tasks.'
            : "The workflow did not complete successfully. Review the Nextflow trace, `.nextflow.log`, and the failed task work directory before using any partial outputs.\n\nError: ${workflow.errorMessage ?: 'No error message was reported.'}"
        def pipelineInfo = new File(params.outdir.toString(), 'pipeline_info')
        pipelineInfo.mkdirs()
        new File(pipelineInfo, 'RUN_STATUS.md').text = """# Meristem pipeline run status

| Field | Value |
|---|---|
| Overall status | **${status}** |
| Run name | ${workflow.runName} |
| Started | ${workflow.start} |
| Completed | ${workflow.complete} |
| Duration | ${workflow.duration} |
| Output directory | `${params.outdir}` |
| Work directory | `${workflow.workDir}` |

## Where to look

- Annotation reports and per-method manifests: `${params.outdir}/reports/`
- Annotation files: `${params.outdir}/gene_prediction/`
- Assembly/repeat QC: `${params.outdir}/assembly_qc/` and `${params.outdir}/repeats/`
- Nextflow execution details: `${params.outdir}/pipeline_info/`

${conclusion}
"""
        def manifestRows = params
            .sort { left, right -> left.key.toString() <=> right.key.toString() }
            .collect { key, value -> "| `${key}` | `${value}` |" }
            .join('\n')
        new File(pipelineInfo, 'PIPELINE_SOURCE.sha256').text = source_listing
        new File(pipelineInfo, 'RUN_MANIFEST.md').text = """# Meristem run manifest

This file records the inputs and options used for this Nextflow invocation.

| Field | Value |
|---|---|
| Run name | `${workflow.runName}` |
| Nextflow version | `${nextflow.version}` |
| Command line | `${workflow.commandLine}` |
| Project directory | `${projectDir}` |
| Work directory | `${workflow.workDir}` |

## Provenance

| Item | Value |
|---|---|
| Meristem commit (from `.git/HEAD`) | `${meristem_commit}` |
| Meristem source digest | `${source_digest}` (sha256 of `PIPELINE_SOURCE.sha256`) |
| BRAKER4 commit | `${run_braker4 ? (braker4_installed_commit ?: 'unpinned') : 'not used'}` |
| Helixer model md5 | `${helixer_model_md5}` |

`PIPELINE_SOURCE.sha256` lists a checksum for each pipeline file. Run `sha256sum -c` on it from
a clean checkout of the commit above: any file that fails was edited after that commit.

BRAKER4 install record:

```
${braker4_record}
```

### Input assemblies

| Sample | Assembly | md5 |
|---|---|---|
${input_checksums.collect { "| ${it[0]} | `${it[1]}` | `${it[2]}` |" }.join('\n')}

## Parameters and required files

| Parameter | Value |
|---|---|
${manifestRows}

## Workflow graph and execution evidence

- `report.html`: Nextflow execution report.
- `timeline.html`: task timeline and resource usage.
- `trace.tsv`: machine-readable task trace.
- `dag.html` (when enabled by the Nextflow configuration): workflow DAG.
"""
        if (workflow.success) {
            def packageScript = new File(projectDir.toString(), 'scripts/package_deliverables.py')
            def packageProcess = new ProcessBuilder(
                'python3', packageScript.toString(), '--outdir', params.outdir.toString()
            ).redirectErrorStream(true).start()
            def packageOutput = packageProcess.inputStream.text
            def packageExit = packageProcess.waitFor()
            def packageStatus = packageExit == 0
                ? "\n\nDeliverables package created:\n${packageOutput.trim()}\n"
                : "\n\nDeliverables package creation failed (exit ${packageExit}):\n${packageOutput.trim()}\n"
            new File(pipelineInfo, 'RUN_STATUS.md').append(packageStatus)
        }
    }
}
