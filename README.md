# Meristem plant genome annotation

This Nextflow DSL2 pipeline performs:

```text
assembly validation -> EDTA repeat discovery -> RepeatMasker soft masking
-> BRAKER3 protein-only prediction -> BUSCO -> optional InterProScan/eggNOG
```

It is designed for the scaffolded plant assemblies produced by the Meristem project. It does not delete inputs or shared databases.

## Status

The workflow is implemented and ready for a small test run. Production use requires compatible containers and pre-installed databases. BRAKER3 protein-only mode does not require RNA-seq.

## Input samplesheet

```csv
sample,assembly
F10702,/home/glbrc.org/benucci/nf-broom/results_F10702_prev_2026-07-30/assembly/decontam/F10702/F10702_decontam.fasta
```

Protein evidence is supplied separately as one FASTA file. Protein IDs must be unique and sequence headers should be simple.

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

## Outputs

- `assembly_qc/`: assembly statistics
- `repeats/`: EDTA libraries, EDTA annotation, RepeatMasker reports, and soft-masked FASTA
- `gene_prediction/`: BRAKER3 GFF3, proteins, CDS, and logs
- `busco/`: completeness reports
- `functional/`: optional InterProScan and eggNOG results
- `reports/`: per-sample summaries

Nextflow reports should be written under `pipeline_info/` in the output directory.

## Safety

The workflow only writes to its Nextflow work directory and `--outdir`. It does not run cleanup commands, modify input assemblies, or write into the shared database directory. Configure filesystem ACLs and the host application's approval policy separately if stronger protection is required.
