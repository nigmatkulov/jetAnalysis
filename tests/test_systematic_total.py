"""Contract checks for the legacy CSV adapter; synthetic values are test-only."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from hist_analysis.python.systematic_total import (
    combine_systematics, archive_systematic_inputs, write_total_systematic_csv,
)


class TotalSystematicsTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.metadata = dict(observable='full', frame='CM', direction='combined',
            selection='jetId', eta_cut=1.9, units='fraction', convention='symmetric_relative',
            target_sample='pPb', binning='normal', production_sample='data', extraction='fit')
        self.config = {}
        for name, value in [('a', .03), ('b', .04)]:
            filename = f'{name}_full_relative_jetId_etaCM_19_ptave_50_60.csv'
            np.savetxt(self.directory / filename, [[-1, 1, 1, value], [1, 1, 1, value]], delimiter=',')
            self.config[name] = dict(directory=self.directory, pattern=filename, metadata=dict(self.metadata))

    def combine(self):
        return combine_systematics(self.config, [(50, 60)], binning='normal')

    def test_quadrature_and_hashes(self):
        tables, provenance = self.combine()
        np.testing.assert_allclose(tables[(50, 60)].total, .05)
        self.assertEqual(len(provenance), 2)
        self.assertTrue(provenance.sha256.str.len().eq(64).all())

    def test_archive_and_total_csv(self):
        tables, provenance = self.combine()
        output = self.directory / 'output'
        archived = archive_systematic_inputs(provenance, output)
        for record in archived.to_dict('records'):
            self.assertEqual(Path(record['csv']).read_bytes(),
                             (output / record['copied_csv']).read_bytes())
        path = write_total_systematic_csv(tables[(50, 60)], output / 'total.csv')
        np.testing.assert_allclose(np.loadtxt(path, delimiter=','),
                                   [[-1, 1, 1, .05], [1, 1, 1, .05]])
        self.assertNotIn('copied_csv', provenance)

    def test_archive_rejects_changed_input(self):
        _, provenance = self.combine()
        Path(provenance.iloc[0]['csv']).write_text('changed after combination')
        with self.assertRaisesRegex(ValueError, 'Input changed'):
            archive_systematic_inputs(provenance, self.directory / 'output')

    def test_incompatible_units_direction_binning(self):
        for key, value in [('units', 'percent'), ('direction', 'pgoing'), ('binning', 'test')]:
            with self.subTest(key=key):
                self.config['b']['metadata'][key] = value
                with self.assertRaises(ValueError):
                    self.combine()
                self.config['b']['metadata'][key] = self.metadata[key]

    def test_duplicate_file(self):
        self.config['b'] = self.config['a']
        with self.assertRaisesRegex(ValueError, 'counted more than once'):
            self.combine()

    def test_missing_and_ambiguous(self):
        original = self.config['a']['pattern']
        self.config['a']['pattern'] = 'absent.csv'
        with self.assertRaisesRegex(FileNotFoundError, 'a, ptAve'):
            self.combine()
        self.config['a']['pattern'] = '*.csv'
        with self.assertRaisesRegex(ValueError, 'Ambiguous'):
            self.combine()
        self.config['a']['pattern'] = original

    def test_invalid_values_and_x_bins(self):
        path = self.directory / self.config['b']['pattern']
        for values in (
            [[-1, 1, 1, np.nan], [1, 1, 1, .04]],
            [[-1, 1, 1, -.04], [1, 1, 1, .04]],
            [[-1, 1, 1, .04], [.9, 1, 1, .04]],
            [[-1, 1, 1, .04], [1.1, 1, 1.1, .04]],
        ):
            np.savetxt(path, values, delimiter=',')
            with self.assertRaises(ValueError):
                self.combine()

    def test_test_mode_canonical_interval(self):
        for config in self.config.values():
            old = self.directory / config['pattern']
            config['pattern'] = config['pattern'].replace('50_60', '50_100')
            old.rename(self.directory / config['pattern'])
            config['metadata']['binning'] = 'test'
        tables, _ = combine_systematics(self.config, [(50, 100)], binning='test')
        np.testing.assert_allclose(tables[(50, 100)].total, .05)


if __name__ == '__main__':
    unittest.main()
