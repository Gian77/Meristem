#!/usr/bin/env bash
# Print the latest Nextflow session recorded for this output and work directory.
# Prints nothing when the run has no recorded session.
set -euo pipefail

history_file=$1
outdir=$2
workdir=$3

[[ -f "$history_file" ]] || exit 0

awk -F '\t' -v outdir="--outdir ${outdir} " -v workdir="-w ${workdir} " '
    index($7 " ", outdir) && index($7 " ", workdir) { session = $6 }
    END { if (session != "") print session }
' "$history_file"
