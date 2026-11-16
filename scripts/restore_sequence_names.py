#!/usr/bin/env python3
"""Put the assembly's own sequence names back into Meristem's published files.

EDTA needs short sequence names, so Meristem renames every sequence to ctg000000001,
ctg000000002, ... before repeat discovery, and every later step uses those names. This script
translates them back with EDTA's id map (normalized_id <TAB> original header), for delivery:

  rename   GFF3, FASTA (headers only) or RepeatMasker .out: every ctg<9 digits> token becomes
           the original sequence name. Gene IDs that contain the short name change with it, and
           a doubled sample prefix (F25101_F25101_Chr01_000001, from Helixer IDs) is collapsed.
  primary  One protein per gene, the longest isoform, using the GFF3 mRNA -> gene links.
           Trailing stop codons (*) are removed and internal ones become X, as InterProScan
           requires.

The script refuses to write a file in which a short name has no entry in the map, so a partial
translation can never be published.
"""
import argparse
import re
import sys

SHORT_NAME = re.compile(r'ctg\d{9}')


def read_id_map(path):
    """Return normalized_id -> original sequence name (first word of the original header)."""
    names = {}
    with open(path) as handle:
        header = next(handle, '')
        if not header.startswith('normalized_id'):
            sys.exit(f'{path}: expected a normalized_id/original_id header line')
        for line in handle:
            if not line.strip():
                continue
            short, original = line.rstrip('\n').split('\t')[:2]
            names[short] = original.split()[0]
    if not names:
        sys.exit(f'{path}: the id map is empty')
    return names


def make_translator(names, sample):
    doubled = f'{sample}_{sample}_' if sample else None

    def translate(text, where):
        def swap(match):
            short = match.group(0)
            if short not in names:
                sys.exit(f'{where}: {short} is not in the id map')
            return names[short]
        text = SHORT_NAME.sub(swap, text)
        if doubled:
            text = text.replace(doubled, f'{sample}_')
        return text
    return translate


def rename(path_in, path_out, kind, names, sample):
    translate = make_translator(names, sample)
    changed = 0
    with open(path_in) as src, open(path_out, 'w') as out:
        for number, line in enumerate(src, 1):
            if kind == 'fasta' and not line.startswith('>'):
                out.write(line)
                continue
            new = translate(line, f'{path_in}:{number}')
            changed += new != line
            out.write(new)
    print(f'{path_out}: {changed} lines renamed')


def read_fasta(path):
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        yield name, ''.join(chunks)


def transcript_genes(gff):
    links = {}
    with open(gff) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] not in ('mRNA', 'transcript'):
                continue
            attrs = dict(x.split('=', 1) for x in fields[8].split(';') if '=' in x)
            if 'ID' in attrs and 'Parent' in attrs:
                links[attrs['ID']] = attrs['Parent'].split(',')[0]
    return links


def primary(gff, proteins, path_out):
    links = transcript_genes(gff)
    if not links:
        sys.exit(f'{gff}: no mRNA or transcript features with ID and Parent')
    best, order, unlinked = {}, [], 0
    for name, seq in read_fasta(proteins):
        gene = links.get(name)
        if gene is None:
            unlinked += 1
            continue
        if gene not in best:
            order.append(gene)
        if gene not in best or len(seq) > len(best[gene][1]):
            best[gene] = (name, seq)
    if not best:
        sys.exit(f'{proteins}: no protein ID matches a GFF3 transcript ID')
    if unlinked:
        sys.exit(f'{proteins}: {unlinked} protein(s) have no GFF3 transcript; files do not correspond')
    with open(path_out, 'w') as out:
        for gene in order:
            name, seq = best[gene]
            seq = seq.rstrip('*').replace('*', 'X')
            out.write(f'>{name} gene={gene}\n')
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + '\n')
    print(f'{path_out}: {len(order)} genes, one protein each')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command')
    r = sub.add_parser('rename')
    r.add_argument('--id-map', required=True)
    r.add_argument('--sample', default='')
    r.add_argument('--kind', choices=['gff', 'fasta', 'rmout'], required=True)
    r.add_argument('input')
    r.add_argument('output')
    p = sub.add_parser('primary')
    p.add_argument('--gff', required=True)
    p.add_argument('--proteins', required=True)
    p.add_argument('output')
    args = parser.parse_args()
    if args.command is None:
        parser.error('choose a command: rename or primary')
    if args.command == 'rename':
        rename(args.input, args.output, args.kind, read_id_map(args.id_map), args.sample)
    else:
        primary(args.gff, args.proteins, args.output)


if __name__ == '__main__':
    main()
