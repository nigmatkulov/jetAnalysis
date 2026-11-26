"""Analytic weighted examples for the diagnostic unfolding chain."""

import unittest

import numpy as np

from hist_analysis.python.notebook_setup import load_root, load_roounfold

ROOT = load_root(batch=True)

from hist_analysis.python.unfolding import FlattenedBinning
from hist_analysis.python.unfolding_validation import (
    Spectrum, _bayes_derivative, forward_backward, histogram, make_training,
    project, ratio, run_unfolding, values, variances,
)


class MatchedValidationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        load_roounfold(ROOT, project_root=".")

    def setUp(self):
        self.template = ROOT.TH1D("eta", "", 2, -1, 1)
        self.template.SetDirectory(0)
        self.flat = ROOT.TH1D("flat", "", 2, .5, 2.5)
        self.flat.SetDirectory(0)
        self.layout = FlattenedBinning(((50., 60.),), 2)
        self.matrix = ROOT.TH2D("toy_matrix", "", 2, .5, 2.5, 2, .5, 2.5)
        self.matrix.SetDirectory(0)
        self.matrix.Sumw2()
        # Rows = Reco, columns = Gen; deliberately asymmetric.
        for i, j, weight, count in ((1, 1, 2., 4), (2, 1, 1., 2),
                                     (1, 2, 1., 2), (2, 2, 3., 2)):
            for _ in range(count):
                self.matrix.Fill(i, j, weight)

    def train(self):
        reco = histogram(self.flat, "total_r", [12., 16.], [22., 36.])
        gen = histogram(self.flat, "total_g", [20., 10.], [38., 22.])
        return make_training(self.matrix, reco, gen, self.layout, [self.template], {})

    def test_weighted_subset_and_bayesian_closure(self):
        t = self.train()
        np.testing.assert_array_equal(values(t.matched_gen), [10, 8])
        np.testing.assert_array_equal(values(t.matched_reco), [10, 8])
        np.testing.assert_allclose(variances(t.matched_gen), [18, 20])
        np.testing.assert_allclose(variances(t.matched_reco), [18, 20])
        np.testing.assert_allclose(values(t.purity), [5/6, .5])
        np.testing.assert_allclose(values(t.inverse_efficiency), [2., 1.25])
        expected_p2 = ((1-.5)**2 * 20 + .5**2 * 16)/16**2
        expected_e1 = ((1-.5)**2 * 18 + .5**2 * 20)/20**2
        self.assertAlmostEqual(variances(t.purity)[1], expected_p2)
        self.assertAlmostEqual(variances(t.inverse_efficiency)[0], expected_e1/.5**4)
        nominal = run_unfolding(t, t.reco, iterations=3, scale=1.)
        large = run_unfolding(t, t.reco, iterations=3, scale=1.e12)
        np.testing.assert_allclose(values(nominal.unfolded.histogram), [20, 10])
        np.testing.assert_allclose(values(nominal.refolded.histogram), [12, 16])
        np.testing.assert_allclose(large.unfolded.covariance, nominal.unfolded.covariance)
        np.testing.assert_allclose(values(nominal.bundle.response.Hmeasured()), [10, 8])

    def test_first_iteration_and_correlated_refold(self):
        t = self.train()
        reco = histogram(self.flat, "pseudo", [18., 8.], [16., 9.])
        vr = np.array([[16., 2.], [2., 9.]])
        result = run_unfolding(t, reco, iterations=1, measured_covariance=vr)
        # Bayes for one iteration with its training prior: B[j,i]=M[i,j]/row_sum[i].
        posterior = np.array([[.8, .25], [.2, .75]])
        purity = np.diag([5/6, .5])
        efficiency_inverse = np.diag([2., 1.25])
        jacobian = efficiency_inverse @ posterior @ purity
        expected = jacobian @ np.array([18., 8.])
        np.testing.assert_allclose(values(result.unfolded.histogram), expected)
        np.testing.assert_allclose(result.unfolded.covariance, jacobian @ vr @ jacobian.T)
        a = np.array([[.8, .25], [.2, .75]])
        fold = np.diag([6/5, 2.]) @ a @ posterior @ purity
        np.testing.assert_allclose(result.refolded_cross, fold @ vr)
        rr = ratio(result.refolded, result.reco, name="rr", cross=result.refolded_cross)
        y = np.array([18., 8.])
        derivative = np.diag(1/y) @ fold - np.diag((fold @ y)/y**2)
        np.testing.assert_allclose(rr.covariance, derivative @ vr @ derivative.T)

    def test_iterative_derivative_finite_difference(self):
        a = np.array([[.8, .25], [.2, .75]])
        measured = np.array([15., 4.])
        prior = np.array([10., 8.])
        _, jac = _bayes_derivative(a, measured, prior, 4)
        for i in range(2):
            change = np.eye(2)[i] * 1.e-4
            up, _ = _bayes_derivative(a, measured+change, prior, 4)
            down, _ = _bayes_derivative(a, measured-change, prior, 4)
            np.testing.assert_allclose((up-down)/2.e-4, jac[:, i], rtol=1.e-8)

    def test_projection_and_fb_covariance(self):
        h = histogram(self.template, "eta_spectrum", [5., 10.], [4., 9.])
        v = np.array([[4., 2.], [2., 9.]])
        s = Spectrum(h, v, np.ones(2, dtype=bool))
        total = ROOT.TH1D("sum_eta", "", 1, -1, 1)
        total.SetDirectory(0)
        p = project(s, [[1., 1.]], total, name="sum_eta_projected")
        np.testing.assert_array_equal(values(p.histogram), [15.])
        np.testing.assert_array_equal(p.covariance, [[17.]])
        fb = forward_backward(s, name="fb")
        self.assertAlmostEqual(values(fb.histogram)[0], 2.)
        self.assertAlmostEqual(variances(fb.histogram)[0], 17./25.)

    def test_unsupported_and_flow(self):
        self.matrix.Reset()
        self.matrix.Fill(1, 1, 1.)
        t = self.train()
        self.assertTrue(np.isnan(values(t.inverse_efficiency)[1]))
        result = run_unfolding(t, t.reco, iterations=1, scale=1.)
        self.assertEqual(result.unfolded.supported.tolist(), [True, False])
        self.assertEqual(result.refolded.supported.tolist(), [True, False])
        self.matrix.Fill(0, 1, 1.)
        with self.assertRaisesRegex(ValueError, "flow"):
            self.train()


if __name__ == "__main__":
    unittest.main()
