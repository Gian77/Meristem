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
    ('562', '2'), ('9606', '33208'), ('33208', '2759'), ('5052', '2759'),
    ('4751', '2759'), ('5061', '4751')])


def gene(contig, name, start, cds=1):
    rows = [f'{contig}\tt\tgene\t{start}\t{start + 9}\t.\t+\t.\tID={name}',
            f'{contig}\tt\tmRNA\t{start}\t{start + 9}\t.\t+\t.\tID={name}.t1;Parent={name}']
    rows += [f'{contig}\tt\tCDS\t{start}\t{start + 9}\t.\t+\t0\tID=c{name}{i};Parent={name}.t1' for i in range(cds)]
    return '\n'.join(rows) + '\n'


def hit(name, taxid, bits=100, identity=60, pe=1, desc='Protein'):
    return f'{name}.t1\tsp\t{identity}\t1\t1\t1e-20\t{bits}\t{desc} OS=Thing OX={taxid} PE={pe} SV=1\n'


class GeneOriginTests(unittest.TestCase):
    def run_check(self, gff, hits, id_map=None, review=None):
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
            if review is not None:
                (tmp / 'flagged.fa').write_text(review)
                argv += ['--review-fasta', str(tmp / 'flagged.fa')]
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
        self.assertIn('| On contigs called `bacterial_like`, `fungal_like`, `animal_like` | 4 |', summary)
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

    def test_fungal_and_animal_contigs_are_called(self):
        gff = ''.join(gene('bac', f'f{i}', 1 + 20 * i) for i in range(3))
        gff += ''.join(gene('chr', f'a{i}', 1 + 20 * i) for i in range(3))
        hits = ''.join(hit(f'f{i}', 5061) for i in range(3)) + ''.join(hit(f'a{i}', 9606) for i in range(3))
        table, summary = self.run_check(gff, hits)
        self.assertEqual(table['bac']['call'], 'fungal_like')
        self.assertEqual(table['chr']['call'], 'animal_like')
        self.assertIn('| On contigs called `bacterial_like`, `fungal_like`, `animal_like` | 6 |', summary)
        self.assertIn('| bac | fungal_like |', summary)

    def test_uncertain_entries_are_ignored(self):
        # Plant rDNA genes whose top hit is a dubious yeast ORF (PE=5, like ART2_YEAST)
        # must not make the contig fungal. A gene with a weaker certain hit keeps that hit.
        gff = ''.join(gene('bac', f'r{i}', 1 + 20 * i) for i in range(4))
        hits = ''.join(hit(f'r{i}', 5061, bits=200, pe=5) for i in range(4))
        hits += hit('r0', 3702, bits=50)
        table, summary = self.run_check(gff, hits)
        self.assertEqual(table['bac']['fungi'], '0')
        self.assertEqual(table['bac']['plant'], '1')
        self.assertEqual(table['bac']['hits'], '1')
        self.assertEqual(table['bac']['call'], 'too_few_hits')
        self.assertIn('Genes whose top hit was such an entry: 4;', summary)

    def test_transposon_hits_are_ignored(self):
        # Plant retrotransposon genes whose top hits are yeast or fly Gag-Pol proteins must not
        # make the contig look foreign. Telomerase and RNase HII are host genes and still count.
        gff = ''.join(gene('bac', f't{i}', 1 + 20 * i) for i in range(5))
        hits = hit('t0', 5061, bits=200, desc='Transposon Ty3-I Gag-Pol polyprotein')
        hits += hit('t1', 33208, bits=200, desc='Retrovirus-related Pol polyprotein from transposon opus')
        hits += hit('t2', 33208, bits=200, desc='Copia protein')
        hits += hit('t0', 3702, bits=50)
        hits += hit('t3', 562, desc='Ribonuclease HII') + hit('t4', 5061, desc='Telomerase reverse transcriptase')
        table, summary = self.run_check(gff, hits)
        self.assertEqual(table['bac']['fungi'], '1')
        self.assertEqual(table['bac']['animal'], '0')
        self.assertEqual(table['bac']['bacteria'], '1')
        self.assertEqual(table['bac']['plant'], '1')
        self.assertEqual(table['bac']['call'], 'mixed')
        self.assertIn('Genes whose top hit was a transposable-element protein: 3;', summary)

    def test_canopy_rdna_label_is_marked_not_removed(self):
        gff = ''.join(gene('bac', f'f{i}', 1 + 20 * i) for i in range(3))
        hits = ''.join(hit(f'f{i}', 5061) for i in range(3))
        flagged = ('>bac canopy_decision=keep;rule=none;evidence_taxon=Aniger\nACGT\n'
                   '>chr canopy_decision=keep;rule=none;evidence_taxon=Streptophyta\nACGT\n'
                   '>gone canopy_decision=remove;rule=R1;evidence_taxon=Rdelemar\nACGT\n')
        table, summary = self.run_check(gff, hits, review=flagged)
        self.assertEqual(table['bac']['call'], 'fungal_like')
        self.assertEqual(table['bac']['canopy_rdna_label'], 'Aniger')
        self.assertEqual(table['chr']['canopy_rdna_label'], '')
        self.assertIn('| likely rDNA (Canopy: Aniger) |', summary)
        self.assertIn('**1 of the 1 flagged contig is likely rDNA, not contamination** (3 genes)', summary)
        self.assertIn('Canopy labelled it `Aniger`,', summary)
        self.assertIn('| On contigs called `bacterial_like`, `fungal_like`, `animal_like` | 3 |', summary)

    def test_canopy_review_contigs_are_listed(self):
        gff, hits = self.example()
        review = ('>contig_77 canopy_decision=review;rule=none;evidence_taxon=Aspergillus_niger\nACGT\n'
                  '>chr canopy_decision=keep;rule=none\nACGT\n'
                  '>gone canopy_decision=review;rule=none;evidence_taxon=no-hit\nACGT\n')
        table, summary = self.run_check(gff, hits, 'normalized_id\toriginal_id\nbac\tcontig_77 desc\n', review)
        self.assertEqual(table['bac']['canopy_review'], 'Aspergillus_niger')
        self.assertEqual(table['chr']['canopy_review'], '')
        self.assertIn('## Contigs Canopy kept for review', summary)
        self.assertIn('| contig_77 | Aspergillus_niger | bacterial_like |', summary)
        self.assertIn('1 review contig(s) named by Canopy are not in this genome', summary)

    def test_review_section_absent_without_flagged_fasta(self):
        _, summary = self.run_check(*self.example())
        self.assertNotIn('Canopy kept', summary)

    def test_mismatched_protein_ids_fail_instead_of_reporting_no_hits(self):
        with self.assertRaises(SystemExit):
            self.run_check(gene('chr', 'p0', 1), hit('other', 3702))

    def test_gff_contig_missing_from_genome_fails(self):
        with self.assertRaises(SystemExit):
            self.run_check(gene('absent', 'p0', 1), '')


if __name__ == '__main__':
    unittest.main()
