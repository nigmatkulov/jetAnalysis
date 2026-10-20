"""Checks for JER construction and its existing downstream estimators."""
import math
import unittest
from hist_analysis.python.notebook_setup import load_root
ROOT = load_root(batch=True)
from hist_analysis.python.jer_systematic_methods import construct_jer_variations
from hist_analysis.python.histogram_ops import ratio_to_nominal
from hist_analysis.python.systematic_fits import calculate_bin_by_bin_systematic


class JERMethodsTest(unittest.TestCase):
    def histogram(self, name, values):
        h = ROOT.TH1D(name, '', len(values), 0, len(values))
        h.SetDirectory(0)
        for i, value in enumerate(values, 1):
            h.SetBinContent(i, value)
            h.SetBinError(i, 0.02)
        return h

    def test_reflection_and_raw_estimators(self):
        up = self.histogram('up_method_test', [0.95, 0.98, 1, 1.02, 1.05, 1-1e-12])
        down = self.histogram('down_method_test', [1.01, 0.96, 1, 0.97, 1.03, 1+1e-12])
        raw = {'Up': up, 'Down': down}
        self.assertIs(construct_jer_variations(raw, method='direct_ratio'), raw)
        folded = construct_jer_variations(raw, method='folded_deviation')
        for label, h in raw.items():
            for i in range(1, h.GetNbinsX()+1):
                self.assertAlmostEqual(folded[label].GetBinContent(i), 1+abs(1-h.GetBinContent(i)))
                self.assertEqual(folded[label].GetBinError(i), h.GetBinError(i))
                self.assertGreaterEqual(folded[label].GetBinContent(i), 1)
        for combination in ('average', 'maximum'):
            a = calculate_bin_by_bin_systematic(up, down, name='raw_test', combination=combination)
            b = calculate_bin_by_bin_systematic(folded['Up'], folded['Down'], name='fold_test', combination=combination)
            for i in range(1, up.GetNbinsX()+1):
                self.assertAlmostEqual(a.GetBinContent(i), b.GetBinContent(i), places=15)

    def test_root_division_errors_and_empty_bins(self):
        numerator = self.histogram('num_method_test', [0.95, 1, 0, 1.05])
        denominator = self.histogram('den_method_test', [1, 1, 0, 1])
        for option in ('', 'B'):
            ratio = ratio_to_nominal(numerator, denominator, option=option)
            folded = construct_jer_variations({'ratio': ratio}, method='folded_deviation')['ratio']
            for i in range(1, 5):
                self.assertEqual(folded.GetBinError(i), ratio.GetBinError(i))
            self.assertEqual(folded.GetBinContent(3), 0)
            self.assertEqual(folded.GetBinError(3), 0)
            self.assertEqual(folded.GetBinContent(2), 1)

    def test_invalid_inputs(self):
        h = self.histogram('invalid_method_test', [1])
        with self.assertRaises(ValueError):
            construct_jer_variations({'Up': h}, method='unknown')
        for value in (math.nan, math.inf):
            h.SetBinContent(1, value)
            with self.assertRaises(ValueError):
                construct_jer_variations({'Up': h}, method='folded_deviation')
        h.SetBinContent(1, 1)
        h.SetBinErrorOption(ROOT.TH1.kPoisson)
        with self.assertRaises(ValueError):
            construct_jer_variations({'Up': h}, method='folded_deviation')


if __name__ == '__main__':
    unittest.main()
