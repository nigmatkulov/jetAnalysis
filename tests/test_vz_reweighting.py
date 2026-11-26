"""Analytic checks of the v_z reweighting helpers (no input files needed)."""
import math
import unittest

from hist_analysis.python import vz_reweighting as vz

EDGES = [-15.0 + 0.5 * i for i in range(61)]
CENTRES = [0.5 * (a + b) for a, b in zip(EDGES, EDGES[1:])]


def gaussian_counts(mean, sigma, total):
    values = [math.exp(-0.5 * ((x - mean) / sigma) ** 2) for x in CENTRES]
    scale = total / sum(values)
    return [scale * v for v in values]


class NormalizationAndRatioTests(unittest.TestCase):
    def test_unit_area_excludes_flow_and_scales_errors(self):
        h = vz.histogram_from_values('t_norm', EDGES, gaussian_counts(0.0, 5.0, 1e6))
        h.SetBinContent(0, 1000.0)
        h.SetBinContent(61, 2000.0)
        norm, total = vz.normalize_unit_area(h, 't_norm_n', norm_range=(-15.0, 15.0))
        self.assertAlmostEqual(total, 1e6, delta=1e-6)
        self.assertAlmostEqual(norm.Integral(1, 60), 1.0, places=12)
        for i in (1, 30, 60):
            self.assertAlmostEqual(norm.GetBinError(i), h.GetBinError(i) / total, places=15)
        self.assertEqual(vz.bin_range(h, (-15.0, 15.0)), (1, 60))

    def test_ratio_values_errors_and_invalid_bins(self):
        a = vz.histogram_from_values('t_a', EDGES, [4.0] * 60)
        b = vz.histogram_from_values('t_b', EDGES, [16.0] * 60)
        b.SetBinContent(7, 0.0)
        b.SetBinError(7, 0.0)
        r, invalid = vz.ratio_histogram(a, b, 't_r')
        self.assertEqual(invalid, [7])
        self.assertEqual(r.GetBinContent(7), 0.0)
        self.assertAlmostEqual(r.GetBinContent(1), 0.25, places=15)
        expected = 0.25 * math.hypot(2.0 / 4.0, 4.0 / 16.0)
        self.assertAlmostEqual(r.GetBinError(1), expected, places=15)

    def test_density_is_fraction_over_width(self):
        h = vz.histogram_from_values('t_d', EDGES, [1.0] * 60)
        norm, _ = vz.normalize_unit_area(h, 't_d_n', norm_range=(-15.0, 15.0))
        density = vz.density_from_fraction(norm, 't_d_density')
        self.assertAlmostEqual(density.GetBinContent(5), (1 / 60) / 0.5, places=14)
        self.assertAlmostEqual(density.Integral(1, 60, 'width'), 1.0, places=12)


class FitAndClosureTests(unittest.TestCase):
    def setUp(self):
        # Data: Gaussian; MC = data * g with a known positive polynomial g.
        self.g = (0.8, 0.01, 0.004, 0.0, 1e-5)
        data_counts = gaussian_counts(-1.0, 4.7, 1e9)
        mc_counts = [d * sum(c * x ** k for k, c in enumerate(self.g)) for d, x in zip(data_counts, CENTRES)]
        self.data = vz.histogram_from_values('t_data', EDGES, data_counts)
        self.mc = vz.histogram_from_values('t_mc', EDGES, mc_counts)
        self.data_n, _ = vz.normalize_unit_area(self.data, 't_data_n', norm_range=(-15, 15))
        self.mc_n, _ = vz.normalize_unit_area(self.mc, 't_mc_n', norm_range=(-15, 15))
        self.ratio, _ = vz.ratio_histogram(self.mc_n, self.data_n, 't_ratio')

    def test_fit_recovers_shape_up_to_normalization(self):
        fit = vz.fit_polynomial(self.ratio, name='t_fit', order=4, fit_range=(-15.1, 15.1),
                                start=vz.PRODUCTION_POL8_PARAMETERS)
        self.assertEqual(fit.status, 0)
        # The ratio of unit-normalized shapes equals g / <g>_data: constant factor only.
        factor = fit.parameters[0] / self.g[0]
        for k, c in enumerate(self.g):
            self.assertAlmostEqual(fit.parameters[k], factor * c, delta=1e-9 + 1e-6 * abs(factor * c))
        self.assertLess(fit.chi2, 1e-6)
        self.assertLess(fit.max_linear_relative_difference, 1e-6)

    def test_weight_is_reciprocal_and_closure_is_exact(self):
        fit = vz.fit_polynomial(self.ratio, name='t_fit2', order=4, fit_range=(-15.1, 15.1))
        weight = vz.weight_function(fit.function, 't_weight')
        for x in (-14.9, -3.0, 0.0, 7.5, 15.0):
            self.assertAlmostEqual(weight.Eval(x) * fit.function.Eval(x), 1.0, places=12)
        rw = vz.reweight_histogram(self.mc_n, fit.function, 't_rw', mode='centre')
        rw, _ = vz.normalize_unit_area(rw, 't_rw_n', norm_range=(-15, 15))
        for i in range(1, 61):
            self.assertAlmostEqual(rw.GetBinContent(i) / self.data_n.GetBinContent(i), 1.0, places=6)
        before = vz.agreement(self.mc_n, self.data_n, bins=(1, 60))
        after = vz.agreement(rw, self.data_n, bins=(1, 60))
        self.assertLess(after['chi2'], before['chi2'])
        self.assertEqual(after['ndf'], 59)

    def test_bin_average_of_reciprocal(self):
        f = vz.ROOT.TF1('t_lin', 'pol1', -15, 15)
        f.SetParameters(2.0, 0.1)
        expected = math.log((2.0 + 0.1 * 1.0) / 2.0) / 0.1  # integral of 1/(2+0.1x) over [0,1]
        self.assertAlmostEqual(vz.bin_averaged_weight(f, 0.0, 1.0), expected, places=10)

    def test_stability_detects_non_positive_function(self):
        f = vz.ROOT.TF1('t_neg', 'pol2', -15, 15)
        f.SetParameters(-0.1, 0.0, 0.01)
        checks = vz.stability_checks(f)
        self.assertFalse(checks['all_positive'])
        self.assertEqual(checks['extremum_count'], 1)


class ExportAndMixtureTests(unittest.TestCase):
    def test_cpp_snippet_round_trip_and_convention(self):
        parameters = list(vz.PRODUCTION_POL8_PARAMETERS)
        snippet = vz.cpp_snippet(parameters, order=8, fit_range=(-15.1, 15.1))
        self.assertEqual(vz.parse_cpp_parameters(snippet), parameters)
        self.assertIn('1. / fVzWeight.Eval(vz)', snippet)
        self.assertIn('"pol8", -15.1, 15.1', snippet)

    def test_mixture_requires_unit_fractions(self):
        a = vz.histogram_from_values('t_m1', EDGES, [1.0] * 60)
        b = vz.histogram_from_values('t_m2', EDGES, [3.0] * 60)
        mixed = vz.mixture({'a': a, 'b': b}, {'a': 0.25, 'b': 0.75}, 't_mix')
        self.assertAlmostEqual(mixed.GetBinContent(3), 2.5, places=14)
        with self.assertRaises(ValueError):
            vz.mixture({'a': a, 'b': b}, {'a': 0.5, 'b': 0.6}, 't_mix_bad')

    def test_f_test_monotonic(self):
        self.assertLess(vz.f_test(200.0, 55, 60.0, 54)['p_value'], 1e-6)
        self.assertGreater(vz.f_test(60.5, 55, 60.0, 54)['p_value'], 0.3)


if __name__ == '__main__':
    unittest.main()
