"""Exact stored EPPS21 inputs remain unchanged after the notebook split."""
import unittest
from pathlib import Path
import numpy as np
from hist_analysis.python.notebook_setup import load_root
ROOT = load_root()
from hist_analysis.python.npdf_comparison import prediction_registry, load_prediction


class Epps21SelectionTests(unittest.TestCase):
    def test_stored_central_values_and_asymmetric_totals(self):
        definition = prediction_registry(Path.cwd())['epps21_ct18']
        if not definition['path'].is_file():
            self.skipTest('EPPS21 ROOT input is unavailable')
        source = ROOT.TFile.Open(str(definition['path']))
        try:
            for interval in definition['intervals']:
                directory = f'pt{interval[0]:03d}_{interval[1]:03d}_eta1p9'
                for observable, (central_key, band_key) in definition['keys'].items():
                    central, band = load_prediction(definition, interval, observable, 1.9)
                    original = source.Get(f'{directory}/{central_key}')
                    total = source.Get(f'{directory}/{band_key}')
                    np.testing.assert_array_equal(
                        [(central.GetBinContent(i),central.GetBinError(i)) for i in range(central.GetNbinsX()+2)],
                        [(original.GetBinContent(i),original.GetBinError(i)) for i in range(original.GetNbinsX()+2)])
                    np.testing.assert_array_equal(
                        [(band.GetPointX(i),band.GetPointY(i),band.GetErrorXlow(i),band.GetErrorXhigh(i),band.GetErrorYlow(i),band.GetErrorYhigh(i)) for i in range(band.GetN())],
                        [(total.GetPointX(i),total.GetPointY(i),total.GetErrorXlow(i),total.GetErrorXhigh(i),total.GetErrorYlow(i),total.GetErrorYhigh(i)) for i in range(total.GetN())])
        finally:
            source.Close()
