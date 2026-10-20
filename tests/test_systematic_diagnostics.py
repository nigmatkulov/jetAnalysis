"""ROOT-error semantics and isolation of optional systematic comparisons."""
import math
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from hist_analysis.python.notebook_setup import load_root
ROOT = load_root(batch=True)
from hist_analysis.python.histogram_ops import ratio_to_nominal
from hist_analysis.python.systematic_fits import smooth_systematic_running_max
from hist_analysis.python.systematic_diagnostics import (
    estimate_ratio_error_variant, make_systematic_diagnostics,
)


def histogram(name, contents, errors):
    h = ROOT.TH1D(name, '', len(contents), -2, 2)
    h.SetDirectory(0)
    for i, (value, error) in enumerate(zip(contents, errors), 1):
        h.SetBinContent(i, value)
        h.SetBinError(i, error)
    return h


def snapshot(h):
    return tuple((h.GetBinContent(i), h.GetBinError(i)) for i in range(h.GetNbinsX()+2))


class SystematicDiagnosticsTest(unittest.TestCase):
    def test_installed_root_weighted_division_formula(self):
        # Includes a weighted ratio above one, equal contents with unequal
        # variances, and a zero denominator. Option B is not an efficiency here.
        nums = (80, 120, 100, 0, 25)
        dens = (100, 100, 100, 0, 100)
        ne = (12, 14, 13, 0, 8)
        de = (15, 15, 15, 0, 15)
        n = histogram('diag_formula_n', nums, ne)
        d = histogram('diag_formula_d', dens, de)
        standard = ratio_to_nominal(n, d, name='diag_formula_standard', option='')
        binomial = ratio_to_nominal(n, d, name='diag_formula_binomial', option='B')
        for i, (a, b, ea, eb) in enumerate(zip(nums, dens, ne, de), 1):
            self.assertEqual(standard.GetBinContent(i), binomial.GetBinContent(i))
            if b == 0:
                self.assertEqual(standard.GetBinError(i), 0)
                self.assertEqual(binomial.GetBinError(i), 0)
                continue
            r = a/b
            std_var = (ea*ea + r*r*eb*eb)/(b*b)
            b_var = abs(((1-2*r)*ea*ea + r*r*eb*eb)/(b*b)) if a != b else 0
            self.assertAlmostEqual(standard.GetBinError(i)**2, std_var, places=14)
            self.assertAlmostEqual(binomial.GetBinError(i)**2, b_var, places=14)

    def test_variants_preserve_inputs_and_estimator_configuration(self):
        d = histogram('diag_pipeline_d', [1000]*8, [math.sqrt(1000)]*8)
        up_values = [940, 980, 1010, 1000, 1040, 970, 1070, 1030]
        down_values = [1030, 1010, 980, 1000, 960, 1050, 970, 1020]
        up = histogram('diag_pipeline_up', up_values, [math.sqrt(v) for v in up_values])
        down = histogram('diag_pipeline_down', down_values, [math.sqrt(v) for v in down_values])
        before = [snapshot(h) for h in (up, down, d)]
        fit_difference_seen = False
        for method in ('direct_ratio', 'folded_deviation'):
            for extraction in ('fit', 'bin_by_bin'):
                for weights in ('', 'W'):
                    for two_sided in (False, True):
                        for combination in (('average', 'maximum') if two_sided else (None,)):
                            numerators = {'Up': up, 'Down': down} if two_sided else {'Reco / Ref': up}
                            for smoothing in (False, True):
                                variants = [estimate_ratio_error_variant(numerators, d,
                                    ratio_option=option, method=method, extraction=extraction,
                                    formula='pol2', fit_options='RQS0'+weights,
                                    initial_values={k: (1., 0., 0.) for k in numerators},
                                    evaluation_range=(-1.9, 1.9), apply_smoothing=smoothing,
                                    smoothing_origin=-0.465+0.00001,
                                    name_prefix=f'diag_unit_{method}_{extraction}_{weights}_{two_sided}_{combination}_{smoothing}_{tag}',
                                    combination=combination)
                                    for tag, option in (('std', ''), ('b', 'B'))]
                                for key in numerators:
                                    for i in range(10):
                                        self.assertEqual(variants[0]['raw_ratios'][key].GetBinContent(i), variants[1]['raw_ratios'][key].GetBinContent(i))
                                        self.assertEqual(variants[0]['ratios'][key].GetBinContent(i), variants[1]['ratios'][key].GetBinContent(i))
                                delta = max(abs(a[0]-b[0]) for a,b in zip(snapshot(variants[0]['systematic']),snapshot(variants[1]['systematic'])))
                                if extraction=='bin_by_bin' or weights=='W':
                                    self.assertLess(delta, 1e-9)
                                else:
                                    fit_difference_seen |= delta > 1e-6
                                for v in variants:
                                    if extraction=='bin_by_bin':self.assertEqual(v['functions'], {})
                                    for summary in v['fit_summaries'].values():
                                        self.assertEqual(summary['formula'], 'pol2')
                                        self.assertEqual(summary['range'], (-1.9, 1.9))
                                        self.assertEqual(summary['status'], 0)
        self.assertTrue(fit_difference_seen)
        self.assertEqual(before, [snapshot(h) for h in (up, down, d)])
        self.assertTrue(all(h.GetListOfFunctions().GetSize()==0 for h in (up, down, d)))

    def test_smoothing_reuses_nominal_and_scales_only_display_clones(self):
        raw = histogram('diag_reuse_raw', [.01, .03, .02, .01, .02, .01, .04, .03], [0]*8)
        smoothed = smooth_systematic_running_max(raw, name='diag_reuse_smooth', evaluation_range=(-1.9,1.9), smoothing_origin=0)
        before = (snapshot(raw), snapshot(smoothed))
        renderer = Mock(return_value=object())
        kwargs = dict(unsmoothed=raw, nominal=smoothed, numerators=None, denominator=None,
            plot_smoothed_vs_unsmoothed=True, plot_binomial_vs_standard_errors=False,
            method='folded_deviation', extraction='fit', formula='pol2', fit_options='RQS0',
            initial_values=None, evaluation_range=(-1.9,1.9), apply_smoothing=True,
            smoothing_origin=0, renderer=renderer, nominal_csv=Path('/tmp/nominal.csv'),
            source='test', observable='CM full', x_title='eta', x_range=(-2,2),
            y_range=(0,5), annotations=('sample',), style=None)
        with patch('hist_analysis.python.systematic_diagnostics.estimate_ratio_error_variant', side_effect=AssertionError('must not refit')):
            result = make_systematic_diagnostics(**kwargs)['smoothed_vs_unsmoothed']
        self.assertIs(result['curves']['Unsmoothed'], raw)
        self.assertIs(result['curves']['Smoothed'], smoothed)
        for label,h in result['curves'].items():
            percent = result['percent_curves'][label]
            self.assertIsNot(percent,h)
            for i in range(1,9):self.assertAlmostEqual(percent.GetBinContent(i),100*h.GetBinContent(i))
        self.assertTrue(renderer.call_args.kwargs['mplhep'])
        self.assertEqual(result['output'].name,'nominal_diagnostic_smoothed_vs_unsmoothed.pdf')
        self.assertEqual(before, (snapshot(raw),snapshot(smoothed)))
        kwargs.update(apply_smoothing=False, nominal=raw)
        result = make_systematic_diagnostics(**kwargs)['smoothed_vs_unsmoothed']
        self.assertEqual(snapshot(result['curves']['Smoothed']),snapshot(smoothed))
        self.assertIsNot(result['curves']['Smoothed'], raw)
        kwargs.update(plot_smoothed_vs_unsmoothed=False)
        renderer.reset_mock()
        self.assertEqual(make_systematic_diagnostics(**kwargs), {})
        renderer.assert_not_called()


if __name__ == '__main__':
    unittest.main()
