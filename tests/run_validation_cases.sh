#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
nextflow_bin=${NEXTFLOW_BIN:-nextflow}
scratch_dir=$(mktemp -d "${TMPDIR:-/tmp}/meristem-validation.XXXXXX")
trap 'rm -rf "$scratch_dir"' EXIT

run_expected_failure() {
    local name=$1
    shift
    if "$@" >"$scratch_dir/$name.log" 2>&1; then
        echo "FAIL: $name unexpectedly succeeded" >&2
        cat "$scratch_dir/$name.log" >&2
        exit 1
    fi
    echo "PASS: $name rejected invalid input"
}

cd "$project_dir"

"$nextflow_bin" run main.nf -profile test -stub-run \
    --input tests/data/stress.samplesheet.csv \
    --outdir "$scratch_dir/stress-results" \
    -w "$scratch_dir/stress-work" >/"$scratch_dir/stress.log" 2>&1
echo "PASS: five-sample stub stress run"

"$nextflow_bin" run main.nf -profile test -stub-run \
    --gene_origin_db tests/data/proteins.faa --gene_origin_taxonomy tests/data/gene_origin_nodes.dmp \
    --outdir "$scratch_dir/gene-origin-results" \
    -w "$scratch_dir/gene-origin-work" >"$scratch_dir/gene-origin.log" 2>&1
test -s "$scratch_dir/gene-origin-results/gene_origin/TEST.braker4.gene_origin.md"
echo "PASS: gene origin stub run"

run_expected_failure gene_origin_without_taxonomy \
    "$nextflow_bin" run main.nf -profile test -stub-run \
    --gene_origin_db tests/data/proteins.faa \
    --outdir "$scratch_dir/gene-origin-missing-results" -w "$scratch_dir/gene-origin-missing-work"

run_expected_failure missing_column \
    "$nextflow_bin" run main.nf -profile test -stub-run \
    --input tests/data/adversarial_missing_column.samplesheet.csv \
    --outdir "$scratch_dir/missing-column-results" -w "$scratch_dir/missing-column-work"

run_expected_failure missing_assembly \
    "$nextflow_bin" run main.nf -profile test -stub-run \
    --input tests/data/adversarial_missing_assembly.samplesheet.csv \
    --outdir "$scratch_dir/missing-assembly-results" -w "$scratch_dir/missing-assembly-work"

run_expected_failure duplicate_sample \
    "$nextflow_bin" run main.nf -profile test -stub-run \
    --input tests/data/adversarial_duplicate_sample.samplesheet.csv \
    --outdir "$scratch_dir/duplicate-sample-results" -w "$scratch_dir/duplicate-sample-work"

run_expected_failure no_methods \
    "$nextflow_bin" run main.nf -profile test -stub-run \
    --run_braker4 false --run_helixer false \
    --outdir "$scratch_dir/no-method-results" -w "$scratch_dir/no-method-work"

echo "Validation stress and adversarial cases passed."
