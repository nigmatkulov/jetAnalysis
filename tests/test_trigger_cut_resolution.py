"""The unindexed nominal histogram is a 1.9 selection, never a wider cut."""
import unittest
from unittest.mock import patch
import json
import tempfile
from pathlib import Path
from hist_analysis.python.notebook_setup import load_root
ROOT = load_root()
from hist_analysis.python.data_trigger_stitching import dijet_histogram_key


class TriggerCutResolutionTest(unittest.TestCase):
    def test_nominal_fallback_only_for_physical_1p9_cut(self):
        for frame, base in [('cm', 'hRecoDijetPtEtaCM'), ('lab', 'hRecoDijetPtEtaLab')]:
            def load(filename, key):
                if key != base: raise KeyError(key)
                return object()
            with patch('hist_analysis.python.data_trigger_stitching.resolve_eta_source_key', side_effect=lambda f, k: k), \
                 patch('hist_analysis.python.data_trigger_stitching.load_histogram', side_effect=load):
                self.assertEqual(dijet_histogram_key('fixture.root', frame, 4), base)
                for index in [0, 5, 6, 7]:
                    with self.assertRaises(KeyError):
                        dijet_histogram_key('fixture.root', frame, index)

    def test_existing_indexed_inputs_are_unchanged(self):
        with patch('hist_analysis.python.data_trigger_stitching.resolve_eta_source_key', side_effect=lambda f, k: k), \
             patch('hist_analysis.python.data_trigger_stitching.load_histogram', return_value=object()):
            for index in range(8):
                self.assertEqual(dijet_histogram_key('fixture.root', 'cm', index),
                                 f'hRecoDijetPtEtaCM_{index}')

    def test_archived_metadata_resolves_physical_cut(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'archive.root'
            file = ROOT.TFile(str(path), 'RECREATE')
            histogram = ROOT.TH2D('hRecoDijetPtEtaCM_6', 'Reco (|#eta| < 2.4)', 2, 0, 100, 2, -3, 3)
            histogram.Write()
            ROOT.TNamed('histogramProductionMetadata', json.dumps({
                'eta_cut_index_to_value': {'5': 2.3, '6': 2.4}})).Write()
            file.Close()
            self.assertEqual(dijet_histogram_key(path, 'cm', 5), 'hRecoDijetPtEtaCM_6')
            with self.assertRaisesRegex(ValueError, 'no unique source cut 2.5'):
                dijet_histogram_key(path, 'cm', 6)
