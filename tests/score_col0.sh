#!/usr/bin/env bash
# Score a Meristem Col-0 test run against the Araport11 reference annotation.
#
# Standalone check, not a pipeline step. It:
#   1. keeps the protein-coding genes of the Ensembl TAIR10 GFF3 on chromosomes 1-5,
#   2. lifts them from TAIR10.1 onto the Col-0 assembly with Liftoff, because the
#      assembly's coordinates differ from TAIR10,
#   3. renames the lifted sequence IDs to the pipeline's ctg IDs (EDTA id map),
#   4. compares coding structure (CDS only, UTRs ignored) of each method's
#      predictions with gffcompare, and
#   5. writes <results>/scoring/: per-method gffcompare output, a summary TSV and
#      a README.
#
# It reads the inputs without changing them and refuses to overwrite an existing
# scoring folder. The protein evidence included A. thaliana proteins, so
# BRAKER4's scores are optimistic.
#
# Usage: tests/score_col0.sh [RESULTS_DIR]   (THREADS=4 by default)
set -euo pipefail

PROJECT_DIR=/mnt/cephfs/linuxhome/benucci/Meristem
RESULTS=$(readlink -f "${1:-${PROJECT_DIR}/results/Col0-test-20261112}")
SAMPLE=${SAMPLE:-Col0}
THREADS=${THREADS:-4}
REF_FA=${REF_FA:-/mnt/cephfs/linuxhome/benucci/Canopy/benchmark_test/refs/arabidopsis/TAIR10.1_nuclear.fasta}
REF_GFF=${REF_GFF:-${PROJECT_DIR}/databases/test_col0/Arabidopsis_thaliana.TAIR10.63.gff3.gz}
TARGET_FA="${RESULTS}/assembly_qc/${SAMPLE}_final.fasta"
ID_MAP="${RESULTS}/repeats/edta/${SAMPLE}.edta_id_map.tsv"
OUT="${RESULTS}/scoring"
CACHE=${APPTAINER_IMAGE_CACHE:-/mnt/cephfs/linuxhome/benucci/.apptainer_cache}
LIFTOFF_URI=docker://quay.io/biocontainers/liftoff:1.6.3--pyhdfd78af_2
GFFCOMPARE_URI=docker://quay.io/biocontainers/gffcompare:0.12.10--h9948957_0

for f in "$REF_FA" "$REF_GFF" "$TARGET_FA" "$ID_MAP"; do
    [[ -s "$f" ]] || { echo "Missing input: $f" >&2; exit 1; }
done
declare -A PRED
for m in braker4 helixer; do
    g="${RESULTS}/gene_prediction/${m}/${SAMPLE}.${m}.gff3"
    [[ -s "$g" ]] && PRED[$m]="$g"
done
(( ${#PRED[@]} > 0 )) || { echo "No predictions found under ${RESULTS}/gene_prediction" >&2; exit 1; }
if [[ -e "$OUT" ]]; then
    echo "Scoring folder already exists: $OUT (move it away to rescore)" >&2
    exit 1
fi

image() {  # image <docker uri>: local .img path, pulled once into the shared cache
    local uri=$1 img
    img="${CACHE}/$(echo "${uri#docker://}" | tr '/:' '--').img"
    [[ -s "$img" ]] || apptainer pull "$img" "$uri" >&2
    echo "$img"
}
LIFTOFF_IMG=$(image "$LIFTOFF_URI")
GFFCOMPARE_IMG=$(image "$GFFCOMPARE_URI")

WORK="${OUT}/work"
mkdir -p "$WORK"
cd "$WORK"

# Copies, so the FASTA indexes Liftoff writes never land beside the inputs.
cp "$REF_FA" ref.fa
cp "$TARGET_FA" target.fa

# 1. Protein-coding reference genes on chromosomes 1-5, renamed to the TAIR10.1 IDs.
python3 -I - "$REF_GFF" > ref.protein_coding.gff3 <<'PY'
import gzip, sys
chroms = {"1": "NC_003070.9", "2": "NC_003071.7", "3": "NC_003074.8",
          "4": "NC_003075.7", "5": "NC_003076.8"}
def attrs(col):
    return dict(a.split("=", 1) for a in col.strip().split(";") if "=" in a)
with gzip.open(sys.argv[1], "rt") as fh:
    rows = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("#")]
rows = [r for r in rows if len(r) == 9 and r[0] in chroms]
genes = {attrs(r[8])["ID"] for r in rows
         if r[2] == "gene" and attrs(r[8]).get("biotype") == "protein_coding"}
tx = {attrs(r[8])["ID"] for r in rows
      if r[2] == "mRNA" and attrs(r[8]).get("Parent") in genes
      and attrs(r[8]).get("biotype") == "protein_coding"}
print("##gff-version 3")
for r in rows:
    a = attrs(r[8])
    keep = ((r[2] == "gene" and a.get("ID") in genes)
            or (r[2] == "mRNA" and a.get("ID") in tx)
            or (r[2] not in ("gene", "mRNA") and a.get("Parent") in tx))
    if keep:
        r[0] = chroms[r[0]]
        print("\t".join(r))
PY

# 2. Lift the reference genes onto the Col-0 assembly.
apptainer exec -B "$WORK" "$LIFTOFF_IMG" liftoff \
    -g ref.protein_coding.gff3 -o ref.lifted.gff3 -u ref.unmapped.txt \
    -p "$THREADS" -dir liftoff_intermediate target.fa ref.fa > liftoff.log 2>&1

# 3 and 4. Rename sequences to ctg IDs and keep coding structure only: one
# transcript per mRNA, with its CDS segments as exons.
cds_only() {  # cds_only <gff3> <rename: yes|no>
    python3 -I - "$1" "$ID_MAP" "$2" <<'PY'
import sys
gff, idmap, rename = sys.argv[1:4]
names = {}
if rename == "yes":
    with open(idmap) as fh:
        next(fh)
        for line in fh:
            new, old = line.rstrip("\n").split("\t")[:2]
            names[old] = new
print("##gff-version 3")
with open(gff) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        r = line.rstrip("\n").split("\t")
        if len(r) != 9 or r[2] not in ("mRNA", "CDS"):
            continue
        if rename == "yes":
            if r[0] not in names:
                sys.exit(f"sequence {r[0]} is not in the EDTA id map")
            r[0] = names[r[0]]
        if r[2] == "CDS":
            r[2] = "exon"
            r[8] = ";".join(a for a in r[8].split(";") if a.startswith("Parent="))
        print("\t".join(r))
PY
}
cds_only ref.lifted.gff3 yes > reference.cds.gff3
for m in "${!PRED[@]}"; do
    cds_only "${PRED[$m]}" no > "${m}.cds.gff3"
    apptainer exec -B "$OUT" "$GFFCOMPARE_IMG" gffcompare \
        -r reference.cds.gff3 -o "${OUT}/${SAMPLE}.${m}" "${m}.cds.gff3" > "gffcompare.${m}.log" 2>&1
done

# 5. Summary.
python3 -I - "$OUT" "$SAMPLE" "$WORK" "${!PRED[@]}" > "${OUT}/${SAMPLE}.scoring_summary.tsv" <<'PY'
import re, sys
out, sample, work, methods = sys.argv[1], sys.argv[2], sys.argv[3], sorted(sys.argv[4:])
levels = ["Base", "Exon", "Intron", "Intron chain", "Transcript", "Locus"]
print("sample\tmethod\tlevel\tsensitivity\tprecision")
for m in methods:
    text = open(f"{out}/{sample}.{m}").read()  # gffcompare's -o file holds the stats
    for lvl in levels:
        hit = re.search(rf"^\s*{lvl} level:\s+([\d.]+)\s+\|\s+([\d.]+)", text, re.M)
        if hit:
            print(f"{sample}\t{m}\t{lvl}\t{hit.group(1)}\t{hit.group(2)}")
    for key, label in [("Missed loci", "missed_ref_loci"), ("Novel loci", "novel_pred_loci")]:
        hit = re.search(rf"^\s*{key}:\s+(\d+)/(\d+)", text, re.M)
        if hit:
            print(f"{sample}\t{m}\t{label}\t{hit.group(1)}\t{hit.group(2)}")
PY

n_ref=$(awk -F'\t' '$3=="gene"' ref.protein_coding.gff3 | wc -l)
n_lift=$(awk -F'\t' '$3=="gene"' ref.lifted.gff3 | wc -l)
n_unmapped=$(grep -c . ref.unmapped.txt || true)
cat > "${OUT}/README.md" <<EOF
# ${SAMPLE} annotation scoring

Created by \`tests/score_col0.sh\` on $(date -u +%Y-%m-%dT%H:%MZ).

- Reference: protein-coding genes of \`$(basename "$REF_GFF")\` on chromosomes 1-5 (${n_ref} genes),
  lifted from TAIR10.1 onto \`$(basename "$TARGET_FA")\` with Liftoff 1.6.3: ${n_lift} lifted,
  ${n_unmapped} not lifted. Genes that were not lifted are absent from the reference here.
- Comparison: gffcompare 0.12.10 on coding structure only (CDS segments as exons). UTRs are
  ignored because BRAKER4 predicts almost none.
- \`${SAMPLE}.scoring_summary.tsv\`: sensitivity and precision per level. For
  \`missed_ref_loci\` and \`novel_pred_loci\`, the two columns are the count and the total.
- \`${SAMPLE}.<method>\`: full gffcompare report; \`.tracking\`, \`.loci\` and \`.annotated.gtf\`
  match each predicted transcript to the reference.
- \`work/\`: the reference before and after lifting, the CDS-only files and the logs.

Precision counts every predicted gene that overlaps no reference gene as wrong, including
genes on contigs where Liftoff placed nothing. The protein evidence included A. thaliana
proteins, so BRAKER4's scores are optimistic.
EOF
echo "Scoring written to ${OUT}"
