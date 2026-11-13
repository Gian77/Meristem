"""Reporting checks: interpretation, provenance, and safe offline rendering."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


SPEC = importlib.util.spec_from_file_location(
    'report_renderer', Path(__file__).resolve().parents[1] / 'scripts/render_annotation_reports.py')
report = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(report)


class ReportTests(unittest.TestCase):
    def interpretation(self, busco):
        return report.interpretation('TEST', 'poales_odb10', {}, busco,
                                     {'gene': {'braker4': '10', 'helixer': '12'}},
                                     False, False, False, False)

    def test_interpretation_recognizes_either_method_advantage(self):
        busco = {'braker4': {'complete_pct': '91.2', 'missing_pct': '6.7'},
                 'helixer': {'complete_pct': '98.5', 'missing_pct': '0.9'}}
        text = self.interpretation(busco)
        self.assertIn('Helixer has the stronger BUSCO result', text)
        busco['braker4'], busco['helixer'] = busco['helixer'], busco['braker4']
        self.assertIn('BRAKER4 has the stronger BUSCO result', self.interpretation(busco))

    def test_incomplete_measurements_do_not_become_zero_or_a_ranking(self):
        text = self.interpretation({'braker4': {'complete_pct': '', 'missing_pct': 'nan'}})
        self.assertIn('values are incomplete', text)
        self.assertNotIn('0.0% complete', text)
        page = report.render_analysis_page('TEST', 'lineage', text, {},
                                          {'braker4': {'complete_pct': 'nan'}}, 'pipeline_tools.html')
        self.assertIn('Not measured', page)

    def test_markdown_escapes_html_and_rejects_unsafe_links(self):
        body, toc = report.markdown_sections(
            '# Same\n\n## Same\n\n<script>alert(1)</script>\n\n'
            '[bad](javascript:alert)\n\n| Name | Value |\n|---|---|\n| <img> | **ok** |')
        self.assertNotIn('<script>', body)
        self.assertNotIn('href="javascript:', body)
        self.assertIn('&lt;img&gt;', body)
        self.assertEqual(len({key for key, _ in toc}), 2)
        self.assertIn('<table>', body)

    def test_tool_selection_is_distinct_from_execution(self):
        config = {'run_helixer': True, 'run_braker4': True, 'interproscan': False}
        self.assertIn('not verified', report.tool_status('run_helixer', 'HELIXER', config, set()))
        self.assertEqual(report.tool_status('run_helixer', 'HELIXER', config, {'HELIXER'}),
                         'Observed completed/cached')
        self.assertEqual(report.tool_status('interproscan', 'INTERPROSCAN', config, set()), 'Not selected')
        self.assertEqual(report.tool_status('run_edta', 'EDTA_DISCOVERY', {}, set()), 'Selection not recorded')

    def test_tools_guide_uses_trace_and_does_not_dump_unrelated_config(self):
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            trace = out / 'trace.tsv'
            trace.write_text('name\tstatus\nANNOTATION:HELIXER (TEST)\tCOMPLETED\n'
                             'ANNOTATION:EDTA_DISCOVERY (TEST)\tFAILED\n')
            report.render_tools_report(out, {'run_helixer': True, 'run_edta': True,
                                            'access_token': 'UNRELATED_SECRET_VALUE'}, trace,
                                       'TEST_analysis_report.html')
            page = (out / 'pipeline_tools.html').read_text()
            text = (out / 'pipeline_tools.md').read_text()
            self.assertNotIn('UNRELATED_SECRET_VALUE', page + text)
            self.assertIn('Observed completed/cached', page)
            self.assertIn('Selected; completion not verified here', page)
            self.assertIn('href="TEST_analysis_report.html"', page)
            self.assertIn('Find a tool', page)
            self.assertNotIn('<script src=', page)
            for suffix in ('md', 'html', 'pdf', 'docx'):
                self.assertGreater((out / f'pipeline_tools.{suffix}').stat().st_size, 0)

    def test_report_leaves_multiqc_output_alone_and_names_no_genome(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            (out / 'multiqc_report.html').write_text('written by MultiQC')
            inputs = {
                'summary': '# Summary\n',
                'busco': 'sample\tmethod\tcomplete_pct\tsingle_copy_pct\tduplicated_pct\tfragmented_pct\tmissing_pct\ttotal_buscos\n'
                         'TEST\tbraker4\t90\t88\t2\t5\t5\t10\n',
                'comparison': 'feature\tbraker4\thelixer\ngene\t10\t12\n',
                'assembly': 'sample\tTEST\ntotal_length\t100\n',
            }
            argv = ['x', '--sample', 'TEST', '--outdir', str(out)]
            for name, text in inputs.items():
                (out / f'{name}.tsv').write_text(text)
                argv += [f'--{name}', str(out / f'{name}.tsv')]
            with mock.patch.object(sys, 'argv', argv):
                report.main()
            self.assertEqual((out / 'multiqc_report.html').read_text(), 'written by MultiQC')
            self.assertTrue((out / 'multiqc_data/meristem_annotation_summary.html').exists())
            self.assertTrue((out / 'TEST_analysis_report.html').exists())
            text = (out / 'TEST_annotation_interpretation.md').read_text()
            self.assertNotIn('release candidate', text)
            self.assertNotIn('F10702', text)


if __name__ == '__main__':
    unittest.main()
