#!/usr/bin/env python3
"""Report the likely origin of predicted genes, per contig. Advisory only.

Each gene takes the kingdom of its best protein hit, then each contig is called from its
genes. This never edits an annotation or an assembly: flagged contigs are for a person to
review, and the before/after gene counts show what removing them would change.
"""
import argparse
from collections import Counter, defaultdict
from pathlib import Path
import re
import sys

KINGDOM_TAXIDS = {
    '33090': 'plant', '2': 'bacteria', '2157': 'archaea', '10239': 'virus',
    '4751': 'fungi', '33208': 'animal',
}
KINGDOMS = ['plant', 'bacteria', 'archaea', 'fungi', 'animal', 'other_eukaryote', 'virus',
            'unclassified']
STRONG_IDENTITY = 70.0
# Swiss-Prot protein existence level 5 ("uncertain") marks dubious ORFs. Some lie inside
# rRNA genes (ART2_YEAST in yeast 25S rDNA), so plant rDNA would otherwise look fungal.
UNCERTAIN = re.compile(r'\bPE=5\b')
# Transposable-element proteins. Swiss-Prot has few plant ones, so plant retrotransposon genes
# take their best hit from yeast, fly, worm or retrovirus Gag-Pol and look foreign. Their
# taxonomy says nothing about where the contig came from (gene origin control test, 2026-11-16).
TRANSPOSON = re.compile(r'transpos|retrovir|retroelement|gypsy|copia|Gag-Pol|Gag-Pro-Pol|Pro-Pol|'
                        r'Pol polyprotein|Gag polyprotein|(?<!Telomerase )reverse transcriptase|'
                        r'integrase|Ribonuclease HI\b|K02A2\.6', re.IGNORECASE)
# Canopy evidence labels (fungal references in its screen) that mark plant rDNA arrays:
# rRNA genes are conserved enough to align across kingdoms (Canopy release v6 notes).
RDNA_LABELS = {'Aniger', 'Rdelemar'}
# Calls for contigs whose genes mostly hit one non-plant group. Prokaryotes are checked first.
FOREIGN_CALLS = [('bacterial_like', ('bacteria', 'archaea')), ('fungal_like', ('fungi',)),
                 ('animal_like', ('animal',))]
FOREIGN = [call for call, _ in FOREIGN_CALLS]


def read_parents(path):
    parents = {}
    with open(path) as handle:
        for line in handle:
            fields = line.split('\t|\t')
            parents[fields[0]] = fields[1]
    return parents


def kingdom_of(taxid, parents, cache):
    if taxid in cache:
        return cache[taxid]
    node, result = taxid, 'unclassified'
    for _ in range(100):
        if node in KINGDOM_TAXIDS:
            result = KINGDOM_TAXIDS[node]
            break
        if node == '2759':
            result = 'other_eukaryote'
            break
        if node not in parents or parents[node] == node:
            break
        node = parents[node]
    cache[taxid] = result
    return result


def read_lengths(path):
    lengths = Counter()
    name = None
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def read_id_map(path):
    """Return normalized_id -> original sequence name (first word of the original header)."""
    if not path:
        return {}
    with open(path) as handle:
        next(handle)
        pairs = (line.rstrip('\n').split('\t')[:2] for line in handle if line.strip())
        return {short: original.split()[0] for short, original in pairs}


def read_canopy_labels(path):
    """Return (review, rdna) from the headers of Canopy's <sample>_full_flagged.fasta.

    Each header carries canopy_decision=...;rule=...;evidence_taxon=...
    review: name -> evidence taxon for sequences Canopy kept for review.
    rdna: name -> label for kept sequences whose evidence taxon is a fungal reference that
    Canopy's release notes identify as plant 45S rDNA (RDNA_LABELS).
    """
    if not path:
        return {}, {}
    review, rdna = {}, {}
    with open(path) as handle:
        for line in handle:
            if not line.startswith('>'):
                continue
            name, _, description = line[1:].rstrip('\n').partition(' ')
            decision = re.search(r'canopy_decision=([^;\s]+)', description)
            if not decision or decision.group(1) not in ('keep', 'review'):
                continue
            taxon = re.search(r'evidence_taxon=([^;\s]+)', description)
            taxon = taxon.group(1) if taxon else 'unknown'
            if decision.group(1) == 'review':
                review[name] = taxon
            if taxon in RDNA_LABELS:
                rdna[name] = taxon
    return review, rdna


def read_genes(path):
    """Return gene->contig, transcript->gene and gene->highest CDS count of any transcript."""
    gene_contig, transcript_gene = {}, {}
    cds_counts = Counter()
    with open(path) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            attrs = dict(x.split('=', 1) for x in fields[8].split(';') if '=' in x)
            if fields[2] == 'gene':
                gene_contig[attrs['ID']] = fields[0]
            elif fields[2] in ('mRNA', 'transcript'):
                transcript_gene[attrs['ID']] = attrs['Parent']
            elif fields[2] == 'CDS':
                for parent in attrs['Parent'].split(','):
                    cds_counts[parent] += 1
    gene_cds = Counter()
    for transcript, count in cds_counts.items():
        gene = transcript_gene.get(transcript)
        if gene is not None:
            gene_cds[gene] = max(gene_cds[gene], count)
    return gene_contig, transcript_gene, gene_cds


def read_best_hits(path, transcript_gene, parents):
    """Return (gene->(kingdom, identity), set_aside) from mmseqs rows.

    Each gene keeps its top-bitscore hit, ignoring hits to PE=5 (uncertain) entries and to
    transposable-element proteins. set_aside counts genes whose top hit overall was such an
    entry, as {'uncertain': n, 'transposon': n}.
    """
    best, top_bits, top_skip, cache = {}, {}, {}, {}
    rows = matched = 0
    with open(path) as handle:
        for line in handle:
            # query, target, pident, qcov, tcov, evalue, bits, theader
            fields = line.rstrip('\n').split('\t')
            rows += 1
            gene = transcript_gene.get(fields[0])
            if gene is None:
                continue
            matched += 1
            bits = float(fields[6])
            skip = ('uncertain' if UNCERTAIN.search(fields[7]) else
                    'transposon' if TRANSPOSON.search(fields[7]) else None)
            if bits > top_bits.get(gene, float('-inf')):
                top_bits[gene], top_skip[gene] = bits, skip
            if skip or (gene in best and bits <= best[gene][0]):
                continue
            taxid = re.search(r'\bOX=(\d+)', fields[7])
            kingdom = kingdom_of(taxid.group(1), parents, cache) if taxid else 'unclassified'
            best[gene] = (bits, kingdom, float(fields[2]))
    # A protein/GFF naming mismatch would otherwise report every gene as having no hit.
    if rows and not matched:
        sys.exit('No search hit matches a GFF3 transcript ID: proteins and GFF3 do not correspond')
    set_aside = Counter(reason for reason in top_skip.values() if reason)
    set_aside = {'uncertain': set_aside['uncertain'], 'transposon': set_aside['transposon']}
    return {gene: (kingdom, identity) for gene, (_, kingdom, identity) in best.items()}, set_aside


def call_contig(row, min_hits, min_fraction):
    if not row['genes']:
        return 'no_genes'
    if row['hits'] < min_hits:
        return 'too_few_hits'
    for call, groups in FOREIGN_CALLS:
        if sum(row[group] for group in groups) / row['hits'] >= min_fraction:
            return call
    if row['plant'] / row['hits'] >= min_fraction:
        return 'plant_like'
    return 'mixed'


def summarize(lengths, id_map, gene_contig, gene_cds, best, min_hits, min_fraction, review=None,
              rdna=None):
    review = review or {}
    rdna = rdna or {}
    rows = {contig: Counter() for contig in lengths}
    for gene, contig in gene_contig.items():
        if contig not in rows:
            sys.exit(f'GFF3 contig {contig} is not in the genome FASTA')
        row = rows[contig]
        row['genes'] += 1
        row['single_cds'] += gene_cds[gene] == 1
        if gene in best:
            kingdom, identity = best[gene]
            row['hits'] += 1
            row[kingdom] += 1
            row['strong_bacterial'] += kingdom in ('bacteria', 'archaea') and identity >= STRONG_IDENTITY
    table = []
    for contig, row in rows.items():
        original = id_map.get(contig, contig)
        table.append({'contig': contig, 'original_id': original,
                      'length': lengths[contig], **{k: row[k] for k in
                      ['genes', 'single_cds', 'hits', *KINGDOMS, 'strong_bacterial']},
                      'call': call_contig(row, min_hits, min_fraction),
                      'canopy_review': review.get(original, ''),
                      'canopy_rdna_label': rdna.get(original, '')})
    return table


def write_outputs(table, sample, method, outdir, review=None, set_aside=None):
    set_aside = set_aside or {'uncertain': 0, 'transposon': 0}
    columns = ['contig', 'original_id', 'length', 'genes', 'single_cds', 'hits', *KINGDOMS,
               'strong_bacterial', 'call', 'canopy_review', 'canopy_rdna_label']
    with open(Path(outdir) / f'{sample}.{method}.gene_origin.tsv', 'w') as out:
        out.write('\t'.join(columns) + '\n')
        for row in sorted(table, key=lambda r: -r['length']):
            out.write('\t'.join(str(row[c]) for c in columns) + '\n')

    groups = defaultdict(lambda: Counter())
    for row in table:
        group = groups[row['call']]
        group['contigs'] += 1
        for key in ['length', 'genes', 'single_cds', 'hits', *KINGDOMS]:
            group[key] += row[key]
    total = sum(row['genes'] for row in table)
    flagged = sorted((r for r in table if r['call'] in FOREIGN), key=lambda r: -r['length'])
    flagged_genes = sum(r['genes'] for r in flagged)
    strong = sum(r['strong_bacterial'] for r in table if r['call'] not in FOREIGN)
    foreign_label = ', '.join(f'`{call}`' for call in FOREIGN)

    lines = [
        f'# Gene origin check — {sample} — {method}', '',
        'Advisory only: nothing was removed from the annotation or the assembly. A best hit is '
        'evidence, not proof of origin, and "no hit" does not mean a gene is false.', '',
        '## Gene totals', '',
        '| | Genes |', '|---|---:|',
        f'| All predicted genes | {total} |',
        f'| On contigs called {foreign_label} | {flagged_genes} |',
        f'| Left if those contigs were removed | {total - flagged_genes} |', '',
        '## Contig calls', '',
        '| Call | Contigs | Mb | Genes | Single-CDS | Plant hit | Prokaryotic hit | Fungal hit | '
        'Animal hit | No hit |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for call in ['plant_like', 'mixed', 'too_few_hits', *FOREIGN, 'no_genes']:
        g = groups.get(call)
        if not g:
            continue
        genes = g['genes'] or 1
        lines.append(
            f"| {call} | {g['contigs']} | {g['length'] / 1e6:.2f} | {g['genes']} | "
            f"{100 * g['single_cds'] / genes:.1f}% | {100 * g['plant'] / genes:.1f}% | "
            f"{100 * (g['bacteria'] + g['archaea']) / genes:.1f}% | "
            f"{100 * g['fungi'] / genes:.1f}% | {100 * g['animal'] / genes:.1f}% | "
            f"{100 * (g['genes'] - g['hits']) / genes:.1f}% |")
    lines += ['', f'Genes outside {foreign_label} contigs with a prokaryotic best hit of at least '
              f'{STRONG_IDENTITY:.0f}% identity: {strong}. Conserved proteins can do this in clean '
              'plant sequence, so check whether they cluster before treating them as contamination.',
              '', 'Hits to Swiss-Prot entries of uncertain existence (PE=5) are ignored: some are '
              'dubious ORFs inside rRNA genes, such as ART2 in yeast 25S rDNA, and would make plant '
              f"rDNA look fungal. Genes whose top hit was such an entry: {set_aside['uncertain']}; they "
              'take their best remaining hit, or none.',
              '', 'Hits to transposable-element proteins (Gag-Pol, reverse transcriptase, integrase, '
              'transposase and similar) are also ignored. Swiss-Prot holds few plant ones, so genes '
              'predicted on plant retrotransposons take their best hit from yeast, insect, worm or '
              'retrovirus elements and would make clean plant contigs look fungal or animal. Genes '
              f"whose top hit was a transposable-element protein: {set_aside['transposon']}; they take "
              'their best remaining hit, or none. Many of them are likely transposon genes rather '
              'than host genes.',
              '', '## Contigs to review', '',
              'A contig is called foreign when at least half of its gene hits come from one non-plant '
              'group.', '']
    if flagged:
        lines += ['| Contig | Call | Length | Genes | Prokaryotic hit | Fungal hit | Animal hit | Plant hit | Note |',
                  '|---|---|---:|---:|---:|---:|---:|---:|---|']
        lines += [f"| {r['original_id']} | {r['call']} | {r['length']} | {r['genes']} | "
                  f"{r['bacteria'] + r['archaea']} | {r['fungi']} | {r['animal']} | {r['plant']} | "
                  f"{'likely rDNA (Canopy: ' + r['canopy_rdna_label'] + ')' if r['canopy_rdna_label'] else ''} |"
                  for r in flagged]
        rdna_flagged = [r for r in flagged if r['canopy_rdna_label']]
        if rdna_flagged:
            labels = ' or '.join(f'`{label}`' for label in sorted({r['canopy_rdna_label'] for r in rdna_flagged}))
            n = len(rdna_flagged)
            lines += ['', f"**{n} of the {len(flagged)} flagged contig{'s' if len(flagged) > 1 else ''} "
                      f"{'are' if n > 1 else 'is'} likely rDNA, not contamination** "
                      f"({sum(r['genes'] for r in rdna_flagged)} genes). Canopy labelled "
                      f"{'them' if n > 1 else 'it'} {labels}, the name of the fungal reference "
                      'it partly aligns to. Canopy\'s release notes show that this label usually marks plant '
                      '45S rDNA arrays, because rRNA genes are conserved across kingdoms. Gene predictions '
                      'in rDNA are not real protein-coding genes, and they hit small yeast ORFs that lie in '
                      'or act on yeast rDNA (such as ART2 and RRT15), so the contig looks fungal. These '
                      'contigs are reported, not removed: screen them with your own criteria.']
    else:
        lines.append(f'No contig was called {foreign_label}.')
    if review is not None:
        lines += ['', '## Contigs Canopy kept for review', '',
                  'Canopy keeps a contig in the final genome when its contamination evidence is '
                  'unresolved. Genes on these contigs deserve a look even when the call above is plant.', '']
        reviewed = sorted((r for r in table if r['canopy_review']), key=lambda r: -r['length'])
        missing = len(review) - len(reviewed)
        if reviewed:
            lines += ['| Contig | Canopy evidence | Call | Length | Genes | Hits | Plant | Prokaryotic | Fungal | Animal |',
                      '|---|---|---|---:|---:|---:|---:|---:|---:|---:|']
            lines += [f"| {r['original_id']} | {r['canopy_review']} | {r['call']} | {r['length']} | "
                      f"{r['genes']} | {r['hits']} | {r['plant']} | {r['bacteria'] + r['archaea']} | "
                      f"{r['fungi']} | {r['animal']} |" for r in reviewed]
        else:
            lines.append('Canopy kept no contig of this genome for review.')
        if missing:
            lines += ['', f'{missing} review contig(s) named by Canopy are not in this genome; check that '
                      'the flagged FASTA belongs to the same release as the assembly.']
    Path(outdir, f'{sample}.{method}.gene_origin.md').write_text('\n'.join(lines) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--method', required=True)
    parser.add_argument('--gff', required=True)
    parser.add_argument('--genome', required=True)
    parser.add_argument('--hits', required=True)
    parser.add_argument('--nodes', required=True, help='NCBI taxonomy nodes.dmp')
    parser.add_argument('--id-map', default=None, help='EDTA normalized_id/original_id table')
    parser.add_argument('--review-fasta', default=None,
                        help="Canopy <sample>_full_flagged.fasta; review contigs are listed and "
                             "contigs labelled Aniger/Rdelemar are marked as likely rDNA")
    parser.add_argument('--min-hits', type=int, default=3)
    parser.add_argument('--min-fraction', type=float, default=0.5)
    parser.add_argument('--outdir', default='.')
    args = parser.parse_args()

    gene_contig, transcript_gene, gene_cds = read_genes(args.gff)
    if not gene_contig or not transcript_gene:
        sys.exit(f'No gene and transcript features found in {args.gff}')
    best, set_aside = read_best_hits(args.hits, transcript_gene, read_parents(args.nodes))
    review, rdna = read_canopy_labels(args.review_fasta) if args.review_fasta else (None, {})
    table = summarize(read_lengths(args.genome), read_id_map(args.id_map), gene_contig, gene_cds,
                      best, args.min_hits, args.min_fraction, review, rdna)
    write_outputs(table, args.sample, args.method, args.outdir, review, set_aside)


if __name__ == '__main__':
    main()
