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
    if not path:
        return {}
    with open(path) as handle:
        next(handle)
        return dict(line.rstrip('\n').split('\t')[:2] for line in handle if line.strip())


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
    """Return gene->(kingdom, identity) from mmseqs rows, keeping the top bitscore per gene."""
    best, cache = {}, {}
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
            if gene in best and bits <= best[gene][0]:
                continue
            taxid = re.search(r'\bOX=(\d+)', fields[7])
            kingdom = kingdom_of(taxid.group(1), parents, cache) if taxid else 'unclassified'
            best[gene] = (bits, kingdom, float(fields[2]))
    # A protein/GFF naming mismatch would otherwise report every gene as having no hit.
    if rows and not matched:
        sys.exit('No search hit matches a GFF3 transcript ID: proteins and GFF3 do not correspond')
    return {gene: (kingdom, identity) for gene, (_, kingdom, identity) in best.items()}


def call_contig(row, min_hits, min_fraction):
    if not row['genes']:
        return 'no_genes'
    if row['hits'] < min_hits:
        return 'too_few_hits'
    if (row['bacteria'] + row['archaea']) / row['hits'] >= min_fraction:
        return 'bacterial_like'
    if row['plant'] / row['hits'] >= min_fraction:
        return 'plant_like'
    return 'mixed'


def summarize(lengths, id_map, gene_contig, gene_cds, best, min_hits, min_fraction):
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
        table.append({'contig': contig, 'original_id': id_map.get(contig, contig),
                      'length': lengths[contig], **{k: row[k] for k in
                      ['genes', 'single_cds', 'hits', *KINGDOMS, 'strong_bacterial']},
                      'call': call_contig(row, min_hits, min_fraction)})
    return table


def write_outputs(table, sample, method, outdir):
    columns = ['contig', 'original_id', 'length', 'genes', 'single_cds', 'hits', *KINGDOMS,
               'strong_bacterial', 'call']
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
    flagged = sorted((r for r in table if r['call'] == 'bacterial_like'), key=lambda r: -r['length'])
    flagged_genes = sum(r['genes'] for r in flagged)
    strong = sum(r['strong_bacterial'] for r in table if r['call'] != 'bacterial_like')

    lines = [
        f'# Gene origin check — {sample} — {method}', '',
        'Advisory only: nothing was removed from the annotation or the assembly. A best hit is '
        'evidence, not proof of origin, and "no hit" does not mean a gene is false.', '',
        '## Gene totals', '',
        '| | Genes |', '|---|---:|',
        f'| All predicted genes | {total} |',
        f'| On contigs called `bacterial_like` | {flagged_genes} |',
        f'| Left if those contigs were removed | {total - flagged_genes} |', '',
        '## Contig calls', '',
        '| Call | Contigs | Mb | Genes | Single-CDS | Plant hit | Bacterial hit | No hit |',
        '|---|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for call in ['plant_like', 'mixed', 'too_few_hits', 'bacterial_like', 'no_genes']:
        g = groups.get(call)
        if not g:
            continue
        genes = g['genes'] or 1
        lines.append(
            f"| {call} | {g['contigs']} | {g['length'] / 1e6:.2f} | {g['genes']} | "
            f"{100 * g['single_cds'] / genes:.1f}% | {100 * g['plant'] / genes:.1f}% | "
            f"{100 * (g['bacteria'] + g['archaea']) / genes:.1f}% | "
            f"{100 * (g['genes'] - g['hits']) / genes:.1f}% |")
    lines += ['', f'Genes outside `bacterial_like` contigs with a prokaryotic best hit of at least '
              f'{STRONG_IDENTITY:.0f}% identity: {strong}. Conserved proteins can do this in clean '
              'plant sequence, so check whether they cluster before treating them as contamination.',
              '', '## Contigs to review', '']
    if flagged:
        lines += ['| Contig | Length | Genes | Prokaryotic hit | Plant hit |', '|---|---:|---:|---:|---:|']
        lines += [f"| {r['original_id']} | {r['length']} | {r['genes']} | "
                  f"{r['bacteria'] + r['archaea']} | {r['plant']} |" for r in flagged]
    else:
        lines.append('No contig was called `bacterial_like`.')
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
    parser.add_argument('--min-hits', type=int, default=3)
    parser.add_argument('--min-fraction', type=float, default=0.5)
    parser.add_argument('--outdir', default='.')
    args = parser.parse_args()

    gene_contig, transcript_gene, gene_cds = read_genes(args.gff)
    if not gene_contig or not transcript_gene:
        sys.exit(f'No gene and transcript features found in {args.gff}')
    best = read_best_hits(args.hits, transcript_gene, read_parents(args.nodes))
    table = summarize(read_lengths(args.genome), read_id_map(args.id_map), gene_contig, gene_cds,
                      best, args.min_hits, args.min_fraction)
    write_outputs(table, args.sample, args.method, args.outdir)


if __name__ == '__main__':
    main()
