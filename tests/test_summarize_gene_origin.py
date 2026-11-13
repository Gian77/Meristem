"""Gene origin checks: contig calls, gene totals, and refusal of mismatched inputs."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


SPEC = importlib.util.spec_from_file_location(
    'gene_origin', Path(__file__).resolve().parents[1] / 'scripts/summarize_gene_origin.py')
origin = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(origin)

NODES = ''.join(f'{child}\t|\t{parent}\t|\trank\t|\n' for child, parent in [
    ('1', '1'), ('2', '1'), ('2759', '1'), ('33090', '2759'), ('3702', '33090'),
    ('562', '2'), ('9606', '33208'), ('33208', '2759'), ('5052', '2759')])


def gene(contig, name, start, cds=1):
    rows = [f'{contig}\tt\tgene\t{start}\t{start + 9}\t.\t+\t.\tID={name}',
            f'{contig}\tt\tmRNA\t{start}\t{start + 9}\t.\t+\t.\tID={name}.t1;Parent={name}']
    rows += [f'{contig}\tt\tCDS\t{start}\t{start + 9}\t.\t+\t0\tID=c{name}{i};Parent={name}.t1' for i in range(cds)]
    return '\n'.join(rows) + '\n'


def hit(name, taxid, bits=100, identity=60):
    return f'{name}.t1\tsp\t{identity}\t1\t1\t1e-20\t{bits}\tProtein OS=Thing OX={taxid} PE=1\n'


class GeneOriginTests(unittest.TestCase):
    def run_check(self, gff, hits, id_map=None):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / 'g.gff3').write_text(gff)
            (tmp / 'h.tsv').write_text(hits)
            (tmp / 'n.dmp').write_text(NODES)
            (tmp / 'g.fa').write_text('>bac\nACGTACGT\n>chr\nACGTACGTACGT\n>empty\nAC\n')
            argv = ['x', '--sample', 'S', '--method', 'm', '--gff', str(tmp / 'g.gff3'),
                    '--genome', str(tmp / 'g.fa'), '--hits', str(tmp / 'h.tsv'),
                    '--nodes', str(tmp / 'n.dmp'), '--outdir', str(tmp)]
            if id_map:
                (tmp / 'm.tsv').write_text(id_map)
                argv += ['--id-map', str(tmp / 'm.tsv')]
            with mock.patch.object(sys, 'argv', argv):
                origin.main()
            lines = (tmp / 'S.m.gene_origin.tsv').read_text().splitlines()
            header = lines[0].split('\t')
            table = {r.split('\t')[0]: dict(zip(header, r.split('\t'))) for r in lines[1:]}
            return table, (tmp / 'S.m.gene_origin.md').read_text()

    def example(self):
        gff = ''.join(gene('bac', f'b{i}', 1 + 20 * i) for i in range(4))
        gff += ''.join(gene('chr', f'p{i}', 1 + 20 * i, cds=3) for i in range(3))
        hits = hit('b0', 562) + hit('b1', 562) + hit('b2', 3702)
        hits += hit('p0', 3702) + hit('p1', 3702) + hit('p2', 562, identity=80)
        return gff, hits

    def test_contigs_called_from_their_genes(self):
        table, _ = self.run_check(*self.example())
        self.assertEqual(table['bac']['call'], 'bacterial_like')
        self.assertEqual(table['chr']['call'], 'plant_like')
        self.assertEqual(table['empty']['call'], 'no_genes')
        self.assertEqual(table['bac']['single_cds'], '4')
        self.assertEqual(table['chr']['single_cds'], '0')
        self.assertEqual(table['chr']['strong_bacterial'], '1')

    def test_gene_totals_show_what_removal_would_change(self):
        _, summary = self.run_check(*self.example())
        self.assertIn('| All predicted genes | 7 |', summary)
        self.assertIn('| On contigs called `bacterial_like` | 4 |', summary)
        self.assertIn('| Left if those contigs were removed | 3 |', summary)
        self.assertIn('nothing was removed', summary)

    def test_few_hits_are_not_called(self):
        table, _ = self.run_check(gene('bac', 'b0', 1) + gene('chr', 'p0', 1), hit('b0', 562))
        self.assertEqual(table['bac']['call'], 'too_few_hits')
        self.assertEqual(table['chr']['call'], 'too_few_hits')

    def test_top_bitscore_wins_across_rows(self):
        gff = ''.join(gene('chr', f'p{i}', 1 + 20 * i) for i in range(3))
        hits = hit('p0', 562, bits=50) + hit('p0', 3702, bits=90) + hit('p1', 3702) + hit('p2', 3702)
        table, _ = self.run_check(gff, hits)
        self.assertEqual(table['chr']['bacteria'], '0')
        self.assertEqual(table['chr']['plant'], '3')

    def test_unknown_taxon_is_unclassified_not_guessed(self):
        table, _ = self.run_check(gene('chr', 'p0', 1), hit('p0', 999999))
        self.assertEqual(table['chr']['unclassified'], '1')

    def test_original_contig_names_reach_the_report(self):
        gff, hits = self.example()
        table, summary = self.run_check(gff, hits, 'normalized_id\toriginal_id\nbac\tcontig_77\n')
        self.assertEqual(table['bac']['original_id'], 'contig_77')
        self.assertIn('| contig_77 |', summary)

    def test_mismatched_protein_ids_fail_instead_of_reporting_no_hits(self):
        with self.assertRaises(SystemExit):
            self.run_check(gene('chr', 'p0', 1), hit('other', 3702))

    def test_gff_contig_missing_from_genome_fails(self):
        with self.assertRaises(SystemExit):
            self.run_check(gene('absent', 'p0', 1), '')


if __name__ == '__main__':
    unittest.main()
