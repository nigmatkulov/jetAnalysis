"""Eta-window normalization and complete-bin protection checks."""
import unittest
import math
from hist_analysis.python.notebook_setup import load_root
ROOT = load_root()
from hist_analysis.python.jet_selection import aligned_eta_bins, normalize_eta_curves, ratio_display_range, selection_ratio_to_gen


class EtaWindowTests(unittest.TestCase):
    def setUp(self):
        self.g = ROOT.TH1D('window_gen', '', 160, -4, 4)
        self.g.SetDirectory(0)
        self.g.Sumw2()
        self.s = self.g.Clone('window_selected')
        self.s.SetDirectory(0)
        for i in range(1, 161):
            self.g.SetBinContent(i, 10)
            self.g.SetBinError(i, 2)
            self.s.SetBinContent(i, 5 if abs(self.s.GetBinCenter(i)) < .4 else 8)
        self.g.SetBinContent(0, 1000)

    def test_shared_window_and_unchanged_inputs(self):
        curves = normalize_eta_curves({'Gen': self.g, 'RefSel': self.s}, 'integral', (-.4, .4))
        first, last = aligned_eta_bins(self.g, (-.4, .4))
        self.assertAlmostEqual(curves['Gen'].Integral(first, last), 1)
        self.assertAlmostEqual(curves['RefSel'].Integral(first, last), 1)
        self.assertAlmostEqual(curves['RefSel'].GetBinContent(100)/curves['Gen'].GetBinContent(100), 1.6)
        self.assertEqual(self.g.GetBinContent(100), 10)
        self.assertEqual(self.g.GetBinError(100), 2)
        density = normalize_eta_curves({'Gen': self.g}, 'bin_width', (-.4, .4))['Gen']
        self.assertAlmostEqual(density.Integral(first, last, 'width'), 1)

    def test_boundaries(self):
        self.g.Rebin(8)
        with self.assertRaises(ValueError):
            aligned_eta_bins(self.g, (-1, 1))
        bounds = ratio_display_range(self.g, (-3, 3), (-.335, 1.265))
        self.assertAlmostEqual(bounds[0], 0.)
        self.assertAlmostEqual(bounds[1], 1.2)
        with self.assertRaises(ValueError):
            ratio_display_range(self.g, (-3, 3), (-.8, .8), (-1.2, 1.2))
        with self.assertRaises(ValueError):
            normalize_eta_curves({'Gen': self.g}, 'integral', (.2, .1))


class BinomialRatioTests(unittest.TestCase):
    def test_analytic_weighted_subset_and_equal_yields(self):
        gen = ROOT.TH1D('weighted_gen', '', 2, 0, 2)
        selected = ROOT.TH1D('weighted_selected', '', 2, 0, 2)
        for h in (gen, selected):
            h.SetDirectory(0)
            h.Sumw2()
        # Bin 1: successes weights 2,3; failures weight 4. Bin 2: all pass.
        for weight in (2, 3, 4): gen.Fill(.5, weight)
        for weight in (2, 3): selected.Fill(.5, weight)
        for weight in (1, 2):
            gen.Fill(1.5, weight)
            selected.Fill(1.5, weight)
        efficiency = 5 / 9
        variance = ((1-efficiency)**2 * 13 + efficiency**2 * 16) / 81
        for mode in ('none', 'integral', 'bin_width'):
            ratio = selection_ratio_to_gen(selected, gen, mode, (0, 2), option='B')
            scale = 1 if mode == 'none' else 12 / 8
            self.assertAlmostEqual(ratio.GetBinContent(1), efficiency * scale)
            # Explicit ROOT weighted-B formula on the inputs AFTER normalization.
            numerator_scale = 1 if mode == 'none' else 1 / 8
            gen_scale = 1 if mode == 'none' else 1 / 12
            for bin_index, raw_num, raw_gen, num_var, gen_var in [(1, 5, 9, 13, 29), (2, 3, 3, 5, 5)]:
                n, d = raw_num*numerator_scale, raw_gen*gen_scale
                p = n / d
                expected_variance = (0 if n == d else
                    abs(((1-2*p)*num_var*numerator_scale**2 + p*p*gen_var*gen_scale**2) / d**2))
                self.assertAlmostEqual(ratio.GetBinError(bin_index), math.sqrt(expected_variance))
            self.assertAlmostEqual(ratio.GetBinContent(2), scale)
            if mode == 'none':
                self.assertEqual(ratio.GetBinError(2), 0)
            else:
                self.assertGreater(ratio.GetBinError(2), 0)
        self.assertEqual(selected.GetBinContent(1), 5)
        self.assertEqual(gen.GetBinContent(1), 9)

    def test_independent_errors_unchanged(self):
        gen = ROOT.TH1D('independent_gen', '', 2, 0, 2)
        selected = ROOT.TH1D('independent_selected', '', 2, 0, 2)
        for h in (gen, selected):
            h.SetDirectory(0)
            h.Sumw2()
        gen.Fill(.5, 2)
        gen.Fill(1.5, 4)
        selected.Fill(.5, 1)
        selected.Fill(1.5, 2)
        ratio = selection_ratio_to_gen(selected, gen, 'integral', (0, 2), option='')
        self.assertAlmostEqual(ratio.GetBinContent(1), 1)
        self.assertAlmostEqual(ratio.GetBinError(1), math.sqrt(2))


if __name__ == '__main__':
    unittest.main()
