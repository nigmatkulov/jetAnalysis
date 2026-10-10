"""Regression checks for reusable ROOT inputs and consumer compatibility."""
import os
import json
import re
import tempfile
import unittest
from array import array
from pathlib import Path
from unittest.mock import patch

from hist_analysis.python.notebook_setup import load_root
ROOT = load_root(batch=True)
from hist_analysis.python.histogram_production import produce_source, projection_key, resolve_catalog_file, load_eta_projection, source_catalog, ETA_KEY_PATTERN
from hist_analysis.python.histogram_io import resolve_unfolding_mc_file, resolve_unfolding_truth_file
from hist_analysis.python.data_distributions import project_data_histogram, forward_backward_from_full


class HistogramProductionTest(unittest.TestCase):
    def test_systematic_projection_metadata_and_sparse_copy(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / 'source.root'
            output = directory / 'produced_histograms' / 'output.root'
            file = ROOT.TFile.Open(str(source), 'RECREATE')
            key = 'hRecoDijetPtEtaCMJerDefExtra_5'
            histogram = ROOT.TH2D(key, ';pT;eta', 100, 0, 1000, 72, -3.6, 3.6)
            histogram.Sumw2()
            histogram.Fill(65, .15, 2.5)
            histogram.Write()
            # The response name also matches the eta pattern: it must be copied,
            # never treated as a TH2 projection source.
            sparse = ROOT.THnSparseD('hGenDijetPtEtaCMVsRecoPtEtaCM_5', '', 4,
                                    array('i', [2]*4), array('d', [0]*4), array('d', [1]*4))
            sparse.Sumw2()
            sparse.Fill(array('d', [.5]*4), 2.5)
            sparse.Write()
            ROOT.TNamed('jerSmearingMode', '0').Write()
            file.Close()
            result = produce_source(dict(source=source, output=output,
                category='official_mc', dataset='fixture', direction='pp', trigger=None))
            self.assertEqual(result['copied'], 3)
            file = ROOT.TFile.Open(str(output))
            try:
                projected = file.Get(projection_key(key, 'standard', (60, 70)))
                self.assertEqual(projected.Integral(), 2.5)
                self.assertEqual(projected.GetBinError(projected.FindBin(.15)), 2.5)
                self.assertEqual(file.Get('jerSmearingMode').GetTitle(), '0')
                metadata = json.loads(file.Get('histogramProductionMetadata').GetTitle())
                self.assertEqual(metadata['schema_version'], 3)
                self.assertEqual(metadata['eta_cut_index_to_value']['4'], 1.9)
                self.assertEqual(file.Get(sparse.GetName()).GetBinContent(0, array('i', [0]*4)), 2.5)
            finally:
                file.Close()
            custom = project_data_histogram(output, key, 'eta', (60, 80), name='custom')
            self.assertEqual(custom.Integral(), 2.5)
            # A legacy produced file without stored TH1 projections remains usable.
            file = ROOT.TFile.Open(str(output), 'UPDATE')
            file.Delete(projection_key(key, 'standard', (60, 70)) + ';*')
            file.Close()
            legacy = project_data_histogram(output, key, 'eta', (60, 70), name='legacy')
            self.assertEqual(legacy.Integral(), 2.5)

    def test_forward_backward_uses_raw_weighted_independent_errors(self):
        full = ROOT.TH1D('weighted_fb_fixture', '', 4, -2, 2)
        full.SetDirectory(0)
        full.Sumw2()
        full.Fill(.5, 2.5)
        full.Fill(-.5, 4.)
        full.Fill(1.5, 1.)  # Second backward bin is deliberately empty.
        ratio = forward_backward_from_full(full, name='weighted_fb_ratio')
        self.assertEqual(ratio.GetBinContent(1), 2.5 / 4.)
        self.assertAlmostEqual(ratio.GetBinError(1), (2. * (2.5 / 4.)**2)**.5)
        self.assertEqual(ratio.GetBinContent(2), 0.)
        self.assertEqual(ratio.GetBinError(2), 0.)
        self.assertEqual(full.Integral(), 7.5)

    def test_physical_cut_keys_match_upstream_cut_order(self):
        from hist_analysis.config.histograms import DIJET_ETA_CUTS
        cpp = (Path(__file__).resolve().parents[1] / 'processForestSimple.C').read_text()
        values = re.search(r'etaCuts\{([^}]+)\}', cpp).group(1)
        self.assertEqual(tuple(float(value.strip()) for value in values.split(',')), DIJET_ETA_CUTS)
        for index, label in enumerate(('1p5', '1p6', '1p7', '1p8', '1p9', '2p3', '2p4', '3p0')):
            self.assertEqual(projection_key(f'hGenDijetPtEtaCM_{index}', 'standard', (60, 70)),
                             f'hGenDijetEtaCM_etaCut{label}_standard_ptAve60to70')
        self.assertIsNone(ETA_KEY_PATTERN.match('hRecoDijetPtEtaCMRun_285480'))
        with self.assertRaisesRegex(ValueError, 'Unknown eta-cut index'):
            projection_key('hGenDijetPtEtaCM_8', 'standard', (60, 70))

    def test_legacy_index_projection_is_read_by_both_consumers(self):
        from hist_analysis.python.dijet_closures import _project_eta
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'produced_histograms' / 'legacy.root'
            path.parent.mkdir()
            file = ROOT.TFile.Open(str(path), 'RECREATE')
            key = 'hRecoDijetEtaCM_etaCut5_standard_ptAve60to70'
            histogram = ROOT.TH1D(key, '', 72, -3.6, 3.6)
            histogram.Sumw2()
            histogram.Fill(.15, 2.5)
            histogram.Write()
            file.Close()
            # No copied TH2: success proves that the old persisted key is read.
            actual = project_data_histogram(path, 'hRecoDijetPtEtaCM_5', 'eta', (60, 70), name='data_legacy')
            closure = _project_eta(path, 'hRecoDijetPtEtaCM_5', (60, 70), name='closure_legacy', rebin_eta=1)
            self.assertEqual(actual.Integral(), 2.5)
            self.assertEqual(closure.Integral(), 2.5)

    def test_canonical_and_legacy_bundle_filenames(self):
        from hist_analysis.python.histogram_io import find_produced_histogram_file
        with tempfile.TemporaryDirectory() as directory:
            canonical = Path(directory) / 'official_mc_pythia_combined.root'
            legacy = Path(directory) / 'official_mc__pythia__combined.root'
            legacy.touch()
            with patch('hist_analysis.python.histogram_io.produced_histogram_path', return_value=canonical):
                self.assertEqual(find_produced_histogram_file('official_mc', 'pythia', 'combined'), legacy)
                canonical.touch()
                self.assertEqual(find_produced_histogram_file('official_mc', 'pythia', 'combined'), canonical)
        catalog = source_catalog()
        self.assertIn('official_mc_pythia_cp5_pp', catalog)
        self.assertIn('standalone_cp5_pp', catalog)

    def test_projection_interval_boundaries_and_no_extra_eta_cut(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / 'input.root'
            output = directory / 'produced_histograms' / 'output.root'
            file = ROOT.TFile.Open(str(source), 'RECREATE')
            h = ROOT.TH2D('hRefDijetPtEtaCM_4', '', 100, 0, 1000, 72, -3.6, 3.6)
            h.Sumw2()
            for pt, weight in [(59.9, 1.), (60., 2.), (69.9, 3.), (70., 4.)]:
                h.Fill(pt, 2.15, weight)  # Ref eta may differ from selected reco eta.
            h.Write()
            run = h.Clone('hRecoDijetPtEtaCMRun_285480')
            run.Write()
            file.Close()
            result = produce_source(dict(source=source, output=output,
                category='official_mc', dataset='fixture', direction='pp', trigger=None))
            self.assertEqual(result['projected'], 22)  # Only the eta-cut-index TH2.
            actual = project_data_histogram(output, 'hRefDijetPtEtaCM_4', 'eta', (60, 70), name='boundary')
            self.assertEqual(actual.Integral(), 5.)  # [60,70), without eta truncation.
            self.assertAlmostEqual(actual.GetBinError(actual.FindBin(2.15)), 13.**.5)
            file = ROOT.TFile.Open(str(output))
            try:
                self.assertTrue(file.Get('hRecoDijetPtEtaCMRun_285480'))
                keys = [k.GetName() for k in file.GetListOfKeys()]
                self.assertFalse(any('etaCut285480' in k for k in keys))
            finally:
                file.Close()
            with self.assertRaisesRegex(ValueError, 'overwrite its source'):
                produce_source(dict(source=source, output=source))

    def test_old_source_cut_order_is_rejected_and_metadata_protects_consumers(self):
        from hist_analysis.python.dijet_closures import _project_eta
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / 'old_data.root'
            output = directory / 'produced_histograms' / 'old_data.root'
            file = ROOT.TFile.Open(str(source), 'RECREATE')
            old = (1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.4, 3.0)
            for index, cut in enumerate(old):
                h = ROOT.TH2D(f'hRecoDijetPtEtaCM_{index}',
                             f'Reco (|#eta| < {cut})', 100, 0, 1000, 72, -3.6, 3.6)
                h.Fill(65, .15, index + 1.)
                h.Write()
            file.Close()
            entry = dict(source=source, output=output, category='data',
                         dataset='fixture', direction='combined', trigger='Jet80')
            with self.assertRaisesRegex(ValueError, 'differs from the current'):
                produce_source(entry)
            self.assertFalse(output.exists())
            # An archived bundle can still select physical cuts by its metadata.
            output.parent.mkdir(parents=True, exist_ok=True)
            original = ROOT.TFile.Open(str(source))
            file = ROOT.TFile.Open(str(output), 'RECREATE')
            for index, cut in enumerate(old):
                h = original.Get(f'hRecoDijetPtEtaCM_{index}').Clone()
                file.cd()
                h.Write()
                physical = projection_key(h.GetName(), 'standard', (60, 70), eta_cut_values=old)
                projected = h.ProjectionY(physical, 7, 7, 'e')
                projected.Write()
            ROOT.TNamed('histogramProductionMetadata', json.dumps({
                'schema_version': 2, 'eta_cut_index_to_value': dict(enumerate(old))})).Write()
            file.Close()
            original.Close()
            # Current index 6 means 2.4; archived source metadata determines its index.
            actual = project_data_histogram(output, 'hRecoDijetPtEtaCM_6', 'eta', (60, 70), name='old_2p4')
            self.assertEqual(actual.Integral(), 7.)
            closure = _project_eta(output, 'hRecoDijetPtEtaCM_7', (60, 70), name='old_3p0', rebin_eta=1)
            self.assertEqual(closure.Integral(), 8.)
            custom = project_data_histogram(output, 'hRecoDijetPtEtaCM_6', 'eta', (60, 80), name='custom_old')
            self.assertEqual(custom.Integral(), 7.)
            with self.assertRaisesRegex(ValueError, 'no unique source cut 2.3'):
                project_data_histogram(output, 'hRecoDijetPtEtaCM_5', 'eta', (60, 70), name='absent_2p3')

    def test_unfolding_overrides_take_precedence(self):
        with patch.object(Path, 'is_file', return_value=True):
            with patch.dict(os.environ, {'CP5_UNFOLDING_FILE': '/tmp/custom-cp5.root'}):
                self.assertEqual(resolve_unfolding_mc_file(Path('/source'), 'pythia_cp5', 'pp'),
                                 Path('/tmp/custom-cp5.root'))
            with patch.dict(os.environ, {'CP5_GEN_FILE': '/tmp/custom-truth.root'}):
                self.assertEqual(resolve_unfolding_truth_file(Path('/source'), 'pythia_cp5_standalone', 'pp'),
                                 Path('/tmp/custom-truth.root'))

    def test_catalog_keeps_source_fallback(self):
        entry = {'output': Path('/produced.root'), 'source': Path('/original.root'),
                 'category': 'fixture', 'dataset': 'fixture', 'direction': 'pp', 'trigger': None}
        with patch('hist_analysis.python.histogram_production.source_catalog', return_value={'fixture': entry}):
            with patch.object(Path, 'is_file', return_value=False):
                self.assertEqual(resolve_catalog_file('fixture'), entry['source'])
            with patch.object(Path, 'is_file', return_value=True):
                self.assertEqual(resolve_catalog_file('fixture'), entry['output'])

    def test_non_nominal_stem_does_not_use_nominal_production(self):
        with patch.object(Path, 'is_file', return_value=True):
            self.assertEqual(resolve_unfolding_mc_file(Path('/source'), 'embedding', 'combined', 'trkMax'),
                             Path('/source/embedding/embedding_trkMax.root'))


if __name__ == '__main__':
    unittest.main()
