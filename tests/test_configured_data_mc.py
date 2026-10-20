"""Configured ratios must retain the existing independent-error convention."""
import math
import unittest
from hist_analysis.python.notebook_setup import load_root
ROOT = load_root()
from hist_analysis.python.data_mc_plotting import prepare_configured_ratios


class ConfiguredRatioTests(unittest.TestCase):
    def test_required_ratios_and_errors(self):
        histograms = {}
        for key, value in [('data', 6.), ('official_mc_gen', 3.), ('official_mc_reco', 4.), ('pythia_cp5', 2.)]:
            h = ROOT.TH1D('configured_' + key, '', 2, -1, 1)
            h.SetDirectory(0)
            h.SetBinContent(1, value)
            h.SetBinError(1, value / 10)
            histograms[key] = h
        pairs = [('data', 'official_mc_gen'), ('official_mc_reco', 'official_mc_gen'), ('data', 'pythia_cp5')]
        styles = {key: (ROOT.kRed, 20) for key in histograms}
        results = prepare_configured_ratios(histograms, pairs, styles=styles, name='required')
        for pair, expected in zip(pairs, (2., 4/3, 3.)):
            self.assertAlmostEqual(results[pair].GetBinContent(1), expected)
            self.assertAlmostEqual(results[pair].GetBinError(1), expected * math.sqrt(.02))
            self.assertEqual(results[pair]._valid_ratio_bins, (1,))
        with self.assertRaisesRegex(ValueError, 'EPPS21'):
            prepare_configured_ratios(histograms, [('data', 'epps21_ct18')], styles=styles, name='bad')
        with self.assertWarnsRegex(UserWarning, 'unavailable'):
            self.assertEqual(prepare_configured_ratios(histograms, [('data', 'missing')], styles=styles, name='missing'), {})
        bad = ROOT.TH1D('configured_bad', '', 3, -1, 1)
        bad.SetDirectory(0)
        with self.assertRaisesRegex(ValueError, 'Incompatible binning'):
            prepare_configured_ratios({**histograms, 'bad': bad}, [('data', 'bad')], styles=styles, name='bad')

    def test_all_three_y_axis_categories_accept_exact_limits_and_auto(self):
        from hist_analysis.python.data_mc_plotting import draw_comparison_figure
        import numpy as np
        for category, limits in [('full', (0., .15)), ('ratio', (.7, 1.3)), ('fb', (.6, 1.4))]:
            h = ROOT.TH1D('range_' + category, '', 2, -1, 1)
            h.SetDirectory(0)
            for i in (1, 2):
                h.SetBinContent(i, .1 if category == 'full' else 1.)
                h.SetBinError(i, .01)
            total = dict(edges=np.array([[-1., 0.], [0., 1.]]),
                         fraction=np.array([.02, .03]), csv='synthetic',
                         eta_cut=1., observable='full')
            before = [(h.GetBinContent(i), h.GetBinError(i)) for i in (1, 2)]
            for requested in (limits, None):
                with self.subTest(category=category, requested=requested):
                    canvas = draw_comparison_figure({'Data': h},
                        name=f'limits_{category}_{requested}', x_title='eta',
                        y_title=category, x_range=(-1., 1.), styles={'Data': (ROOT.kBlack, 21)},
                        annotations=(), totals=(('Data', total),),
                        y_range=requested, reference=category != 'full')
                    actual = (canvas.GetUymin(), canvas.GetUymax())
                    self.assertTrue(all(math.isfinite(v) for v in actual))
                    self.assertLess(*actual)
                    if requested is not None:
                        self.assertAlmostEqual(actual[0], requested[0])
                        self.assertAlmostEqual(actual[1], requested[1])
                    else:
                        self.assertLessEqual(actual[0], h.GetBinContent(1)-.01)
                        self.assertGreater(actual[1], h.GetBinContent(1)+.01)
                    self.assertEqual(before, [(h.GetBinContent(i), h.GetBinError(i)) for i in (1, 2)])
                    canvas.Close()

    def test_missing_ratio_bins_keep_finite_axes_and_are_omitted(self):
        from hist_analysis.python.data_mc_plotting import data_over_mc, draw_comparison_figure
        for populated in (True, False):
            numerator = ROOT.TH1D(f'missing_n_{populated}', '', 2, -1, 1)
            denominator = ROOT.TH1D(f'missing_d_{populated}', '', 2, -1, 1)
            numerator.SetDirectory(0)
            denominator.SetDirectory(0)
            if populated:
                numerator.SetBinContent(1, 6.)
                numerator.SetBinError(1, .6)
                denominator.SetBinContent(1, 3.)
                denominator.SetBinError(1, .3)
            ratio = data_over_mc(numerator, denominator,
                name=f'missing_ratio_{populated}', style=(ROOT.kRed, 20))
            for requested in ((.5, 2.5), None):
                with self.subTest(populated=populated, limits=requested):
                    canvas = draw_comparison_figure({'ratio': ratio},
                        name=f'missing_canvas_{populated}_{requested}',
                        x_title='eta', y_title='ratio', x_range=(-1., 1.),
                        styles={'ratio': (ROOT.kRed, 20)}, annotations=(),
                        y_range=requested, reference=True)
                    low, high = canvas.GetUymin(), canvas.GetUymax()
                    self.assertTrue(math.isfinite(low) and math.isfinite(high))
                    self.assertLess(low, high)
                    if requested:
                        self.assertAlmostEqual(low, requested[0])
                        self.assertAlmostEqual(high, requested[1])
                    displayed = next(h for h in canvas._comparison_figure_objects
                                     if h and h.InheritsFrom('TH1')
                                     and h.GetName().endswith(ratio.GetName()))
                    self.assertTrue(math.isnan(displayed.GetBinContent(2)))
                    self.assertEqual(ratio.GetBinContent(2), 0.)
                    canvas.Close()
