"""Restoring assembly sequence names, and choosing one protein per gene."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


SPEC = importlib.util.spec_from_file_location(
    'restore', Path(__file__).resolve().parents[1] / 'scripts/restore_sequence_names.py')
restore = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(restore)

ID_MAP = 'normalized_id\toriginal_id\nctg000000001\tF1_Chr01 some description\nctg000000002\tF1_ptg000123l\n'


class RestoreNameTests(unittest.TestCase):
    def run_cli(self, files, argv):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            for name, text in files.items():
                (tmp / name).write_text(text)
            args = ['x'] + [str(tmp / a) if a in files or a.startswith('out') else a for a in argv]
            with mock.patch.object(sys, 'argv', args):
                restore.main()
            return {p.name: p.read_text() for p in tmp.iterdir() if p.name.startswith('out')}

    def rename(self, kind, text, sample='F1'):
        out = self.run_cli({'map.tsv': ID_MAP, 'in.txt': text},
                           ['rename', '--id-map', 'map.tsv', '--sample', sample, '--kind', kind,
                            'in.txt', 'out.txt'])
        return out['out.txt']

    def test_gff_seqids_and_helixer_ids_use_original_names(self):
        gff = ('##sequence-region ctg000000001 1 100\n'
               'ctg000000001\tHelixer\tgene\t1\t90\t.\t+\t.\tID=F1_ctg000000001_000001\n')
        out = self.rename('gff', gff)
        self.assertIn('##sequence-region F1_Chr01 1 100', out)
        # the doubled sample prefix from Helixer IDs is collapsed
        self.assertIn('F1_Chr01\tHelixer\tgene\t1\t90\t.\t+\t.\tID=F1_Chr01_000001', out)
        self.assertNotIn('ctg0', out)

    def test_braker_ids_without_short_names_are_unchanged(self):
        out = self.rename('gff', 'ctg000000002\tAUGUSTUS\tmRNA\t1\t9\t.\t+\t.\tID=F1_g1.t1;Parent=F1_g1\n')
        self.assertEqual(out, 'F1_ptg000123l\tAUGUSTUS\tmRNA\t1\t9\t.\t+\t.\tID=F1_g1.t1;Parent=F1_g1\n')

    def test_fasta_changes_headers_only(self):
        out = self.rename('fasta', '>F1_ctg000000001_000001.1\nMctgKL\n')
        self.assertEqual(out, '>F1_Chr01_000001.1\nMctgKL\n')

    def test_unknown_short_name_fails_instead_of_half_renaming(self):
        with self.assertRaises(SystemExit):
            self.rename('gff', 'ctg000000009\tx\tgene\t1\t9\t.\t+\t.\tID=g\n')

    def test_primary_keeps_longest_isoform_and_cleans_stops(self):
        gff = ('c\tx\tmRNA\t1\t9\t.\t+\t.\tID=g1.t1;Parent=g1\n'
               'c\tx\tmRNA\t1\t9\t.\t+\t.\tID=g1.t2;Parent=g1\n'
               'c\tx\tmRNA\t1\t9\t.\t+\t.\tID=g2.t1;Parent=g2\n')
        proteins = '>g1.t1\nMKV*\n>g1.t2\nMKVLL*\n>g2.t1\nMA*Q*\n'
        out = self.run_cli({'a.gff3': gff, 'p.faa': proteins},
                           ['primary', '--gff', 'a.gff3', '--proteins', 'p.faa', 'out.faa'])['out.faa']
        self.assertEqual(out, '>g1.t2 gene=g1\nMKVLL\n>g2.t1 gene=g2\nMAXQ\n')

    def test_primary_refuses_proteins_missing_from_gff(self):
        with self.assertRaises(SystemExit):
            self.run_cli({'a.gff3': 'c\tx\tmRNA\t1\t9\t.\t+\t.\tID=g1.t1;Parent=g1\n',
                          'p.faa': '>g1.t1\nM\n>other\nM\n'},
                         ['primary', '--gff', 'a.gff3', '--proteins', 'p.faa', 'out.faa'])


if __name__ == '__main__':
    unittest.main()
