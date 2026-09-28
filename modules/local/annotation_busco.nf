process HELIXER_PROTEINS {
    tag "$sample"
    label 'busco'
    container params.busco_container

    input:
    tuple val(sample), path(gff), path(genome), path(metadata)

    output:
    tuple val(sample), path("${sample}.helixer.proteins.faa"), path(metadata), emit: proteins

    script:
    """
    python - "$gff" "$genome" "${sample}.helixer.proteins.faa" <<'PY'
import sys
from collections import defaultdict

gff_path, fasta_path, protein_path = sys.argv[1:]
complement = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
codon_table = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def fasta_records(path):
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        yield name, ''.join(chunks).upper()

genome = dict(fasta_records(fasta_path))
cds = defaultdict(list)
with open(gff_path) as handle:
    for line in handle:
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        if len(fields) != 9 or fields[2] != 'CDS':
            continue
        attrs = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
        parent = attrs.get('Parent')
        if not parent or fields[0] not in genome:
            continue
        cds[parent.split(',')[0]].append((fields[0], int(fields[3]), int(fields[4]), fields[6], int(fields[7]) if fields[7].isdigit() else 0))

with open(protein_path, 'w') as out:
    for transcript, parts in cds.items():
        strand = parts[0][3]
        parts.sort(key=lambda item: item[1], reverse=(strand == '-'))
        sequence = []
        for contig, start, end, strand, phase in parts:
            fragment = genome[contig][start - 1:end]
            if strand == '-':
                fragment = fragment.translate(complement)[::-1]
            sequence.append(fragment[phase:])
        coding = ''.join(sequence)
        protein = ''.join(codon_table.get(coding[i:i+3], 'X') for i in range(0, len(coding) - 2, 3)).rstrip('*')
        if protein:
            out.write('>' + transcript + '\\n')
            for i in range(0, len(protein), 80):
                out.write(protein[i:i+80] + '\\n')
PY
    test -s "${sample}.helixer.proteins.faa"
    """

    stub:
    """
    printf '>helixer_stub\\nM\\n' > "${sample}.helixer.proteins.faa"
    """
}

process ANNOTATION_BUSCO {
    tag "${sample}:${method}"
    label 'busco'
    container params.busco_container
    publishDir "${params.outdir}/busco_comparison", mode: params.publish_dir_mode

    input:
    tuple val(sample), val(method), path(proteins), path(metadata)

    output:
    tuple val(sample), val(method),
          path("busco_${method}_${sample}"),
          path("busco_${method}_${sample}/run_${params.busco_lineage}/short_summary.txt"), emit: qc

    script:
    def db_args = params.busco_db ? "--download_path '${params.busco_db}' --offline" : ""
    """
    busco -i "$proteins" -l "${params.busco_lineage}" -m proteins \\
      -o "busco_${method}_${sample}" --cpu $task.cpus $db_args
    test -n "\$(find "busco_${method}_${sample}" -type f -name 'short_summary*.txt' -size +0c -print -quit)"
    """

    stub:
    """
    mkdir -p "busco_${method}_${sample}/run_${params.busco_lineage}"
    printf 'C:0.0%%[S:0.0%%,D:0.0%%],F:0.0%%,M:100.0%%,n:1\\n' > \\
      "busco_${method}_${sample}/run_${params.busco_lineage}/short_summary.txt"
    """
}
