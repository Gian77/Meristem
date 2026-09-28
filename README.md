# Meristem plant genome annotation

This Nextflow DSL2 pipeline performs:

```text
assembly validation -> EDTA repeat discovery -> RepeatMasker soft masking
-> independent BRAKER3, BRAKER4/Snakemake, and/or Helixer prediction
-> comparison of alternative annotations
```

It is designed for the scaffolded plant assemblies produced by the Meristem project. It does not delete inputs or shared databases.

## Status

The workflow is implemented and ready for a small test run. Production use requires compatible containers and pre-installed databases. BRAKER3 and BRAKER4 protein-only modes do not require RNA-seq. BRAKER3 remains enabled by default; the alternative predictors are opt-in.

## Alternative annotation methods

BRAKER4 is run through the upstream Snakemake workflow. A checkout is staged
at `vendor/BRAKER4/` and the pinned Snakemake environment is staged at
`.snakemake-env/`; both are local runtime assets ignored by Git. The workflow
does not clone repositories or download databases at run time. Snakemake runs inside one Condor allocation with
`--restart-times 3`, `--rerun-incomplete`, and `--keep-going`, so an interrupted
run can be resumed by rerunning the same command. The nested BRAKER4 rules use
Apptainer for their tool containers and should use a pre-populated cache. The
checkout is currently at commit
`42be27bd08fd1a73bf769c21c2697ac66781d85f`, and the local environment contains
Snakemake 8.18.2:

```bash
nextflow run main.nf -profile condor \
  --input assets/F10702.samplesheet.csv \
  --proteins databases/plant_annotation/ensembl_plants_release63/processed/poales_reference_proteins.faa \
  --run_braker3 false \
  --run_braker4 true \
  --braker4_source vendor/BRAKER4 \
  --braker4_snakemake .snakemake-env/bin/snakemake \
  --braker4_singularity_prefix /path/to/shared/apptainer_cache \
  --busco_db /mnt/cephfs/linuxhome/compbio/nonbackup/busco \
  --outdir results/F10702-braker4
```

BRAKER4 receives the RepeatMasker soft-masked assembly and the configured
protein evidence, so it runs in protein-only EP mode. Its standardized outputs
are published under `gene_prediction/braker4/`.

Helixer is also opt-in and uses the staged local copy of
`gglyptodon/helixer-docker:helixer_v0.3.7_cuda_12.2.2-cudnn8_1` with the
upstream one-step interface and `land_plant` model by default:

```bash
nextflow run main.nf -profile condor \
  --input assets/F10702.samplesheet.csv \
  --run_braker3 false --run_helixer true \
  --helixer_container .apptainer_cache/helixer_v0.3.7_cuda_12.2.2-cudnn8_1.sif \
  --helixer_container_options '--nv' \
  --outdir results/F10702-helixer
```

The land-plant model is staged at
`helixer_models/land_plant/land_plant_v0.3_a_0080.h5`. Set `--helixer_model`
to another local HDF5 model when needed.
Use `--helixer_container_options '--nv'` only on a node with NVIDIA GPU
support; Helixer can run on CPU but will be substantially slower. When both
`--run_braker4 true` and `--run_helixer true` are selected, the workflow writes
feature-count comparisons to `comparison/<sample>.annotation_comparison.tsv`.

## Input samplesheet

```csv
sample,assembly
F10702,/home/glbrc.org/benucci/nf-broom/results_F10702_prev_2026-07-30/assembly/decontam/F10702/F10702_decontam.fasta
```

Protein evidence is supplied separately as one FASTA file. Protein IDs must be unique and sequence headers should be simple.

For a large, redundant protein set, clustering can optionally create a smaller
representative FASTA for BRAKER. The source database is never modified:

```bash
nextflow run main.nf \
  --input assets/F10702.samplesheet.csv \
  --proteins /path/to/plant_proteins.faa \
  --cluster_proteins true \
  --cluster_min_seq_id 0.90 \
  --cluster_coverage 0.80
```

The module uses MMseqs2 `easy-linclust`, 90% minimum sequence identity, 80%
coverage, and target-coverage mode by default. Lowering identity to 80% reduces
runtime more aggressively but can discard divergent homologs that are useful to
BRAKER. The standalone `scripts/cluster_braker_proteins.sh` helper remains
available for preparing evidence outside Nextflow.

The initial curated evidence set is available in the workspace at:

```text
databases/plant_annotation/ensembl_plants_release63/processed/poales_reference_proteins.faa
```

It contains 203,254 proteins from Sorghum bicolor, Zea mays, Setaria italica, and Oryza sativa. The source files and checksums are recorded in `databases/plant_annotation/ensembl_plants_release63/MANIFEST.tsv`.

A copy is also available at `/home/glbrc.org/benucci/plant_annotation/ensembl_plants_release63/processed/poales_reference_proteins.faa` for Condor runs.

## Example production command

```bash
nextflow run main.nf \
  -profile apptainer \
  --input assets/F10702.samplesheet.csv \
  --proteins /mnt/cephfs/linuxhome/compbio/nonbackup/proteins/poales_proteins.faa \
  --busco_lineage poales_odb10 \
  --busco_db /mnt/cephfs/linuxhome/compbio/nonbackup/busco \
  --outdir /path/to/meristem_annotation_results \
  -work-dir /path/to/meristem_annotation_work \
  -with-report /path/to/meristem_annotation_results/pipeline_info/report.html \
  -with-trace /path/to/meristem_annotation_results/pipeline_info/trace.tsv \
  -resume
```

The EDTA and RepeatMasker containers must be available to Apptainer. If a cluster cannot pull containers at runtime, download them into a cache first and override the container parameters in a site config.

## Functional annotation

Functional annotation is deliberately opt-in because InterProScan and eggNOG databases are large:

```bash
nextflow run main.nf -profile apptainer \
  --input samplesheet.csv \
  --proteins proteins.faa \
  --interproscan --interproscan_db /mnt/cephfs/linuxhome/compbio/nonbackup/interproscan \
  --eggnog --eggnog_db /mnt/cephfs/linuxhome/compbio/nonbackup/eggnog \
  --outdir results
```

## Repeat libraries

By default, EDTA discovers a library per genome and that library is passed to RepeatMasker. If EDTA has already been run, use `--run_edta false --repeatmasker_library /path/to/library.fa`.

The shared database directory `/mnt/cephfs/linuxhome/compbio/nonbackup/` is treated as read-only. Database preparation is intentionally outside this workflow.

## GLBRC HTCondor execution

This project follows the `nf-broom` execution pattern: launch the Nextflow head process on the submit node, use `-profile condor` for task submission, and run Apptainer containers on compute nodes. The bundled Nextflow environment should be selected explicitly:

Current Scarcity node availability is shown at [CRG Scarcity](https://crg.glbrc.org/#/scarcity). For large jobs, the Condor profile prefers `scarcity-12.glbrc.org`, followed by `scarcity-9.glbrc.org`; it can still fall back to other eligible nodes when those nodes are unavailable.

For HTCondor troubleshooting, see the [HTCondor troubleshooting guide](https://htcondor.readthedocs.io/en/main/users-manual/troubleshooting.html). For an idle job, start with `condor_q -better-analyze JOB_ID`, then inspect `condor_q -hold` and the job log if it is held or repeatedly returns to idle.

```bash
export HOME=/mnt/cephfs/linuxhome/benucci
export JAVA_HOME=/mnt/cephfs/linuxhome/benucci/.conda/envs/nextflow
export JAVA_CMD=/mnt/cephfs/linuxhome/benucci/.conda/envs/nextflow/bin/java
export PATH=/mnt/cephfs/linuxhome/benucci/.conda/envs/nextflow/bin:$PATH
```

For a supervised head process, set the required protein and BUSCO paths and submit:

```bash
mkdir -p condor_logs
PROTEINS=/mnt/cephfs/linuxhome/compbio/nonbackup/path/to/plant_proteins.faa \
BUSCO_DB=/mnt/cephfs/linuxhome/compbio/nonbackup/busco \
condor_submit run_annotation.submit
```

The launcher writes work files under `nf-work-F10702`, results under `results/F10702`, and Apptainer images under `~/.apptainer_cache`. It runs in the foreground under Condor, matching the working `nf-broom` launcher convention.

The Condor launcher prepends `bin/` to `PATH` so Nextflow worker submissions use the bundled `bin/condor_submit` wrapper. The wrapper routes submissions to the submit host schedd (`scarcity-ap-1.glbrc.org`); override this with `MERISTEM_CONDOR_SCHEDD` if the site scheduler changes.

To monitor a submitted head job every 10 minutes, pass the cluster ID printed by
`condor_submit` to the bundled monitor. It reports queue and history status and
stops when the job completes, is held, or is removed:

```bash
scripts/monitor_condor_job.sh 143159 600 | tee condor_logs/annotation-monitor.log
```

The second argument is optional and is the polling interval in seconds; `600`
is the default.

## Outputs

- `assembly_qc/`: assembly statistics
- `repeats/`: EDTA libraries, EDTA annotation, RepeatMasker reports, and soft-masked FASTA
- `gene_prediction/`: BRAKER3, BRAKER4, and/or Helixer annotations and logs
- `comparison/`: feature-count comparison between BRAKER4 and Helixer
- `busco/`: completeness reports
- `busco_comparison/`: side-by-side BRAKER4 and Helixer BUSCO summaries and TSV comparison
- `functional/`: optional InterProScan and eggNOG results
- `reports/`: per-sample Markdown summaries, MultiQC HTML/data, ZIP report bundles, and `pipeline_tools.md`
- `pipeline_info/`: Nextflow execution report, timeline, and trace

Nextflow reports should be written under `pipeline_info/` in the output directory.

## Safety

The workflow only writes to its Nextflow work directory and `--outdir`. It does not run cleanup commands, modify input assemblies, or write into the shared database directory. Configure filesystem ACLs and the host application's approval policy separately if stronger protection is required.
