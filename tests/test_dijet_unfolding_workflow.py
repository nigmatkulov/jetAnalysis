"""Independent examples for forward folding and the conditional toy ensemble."""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import ROOT

from hist_analysis.python.dijet_unfolding_workflow import (
    bayes_toy_trajectory, fold_expected, poisson_toys, risk_metrics,
)
from hist_analysis.python.histogram_io import load_object
from hist_analysis.python.notebook_setup import load_roounfold
from hist_analysis.python.unfolding import FlattenedBinning, write_unfolding_output
from hist_analysis.python.unfolding_validation import (
    histogram, make_training, matrix_values, run_unfolding, validation_views, values, variances,
)


class DijetStudyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ROOT.gROOT.SetBatch(True)
        ROOT.TH1.AddDirectory(False)
        load_roounfold(ROOT, project_root=".")

    def training(self, migration, reco, gen):
        matrix = ROOT.TH2D("study_migration", "", 4, .5, 4.5, 4, .5, 4.5)
        matrix.SetDirectory(0)
        matrix.Sumw2()
        for i, row in enumerate(migration, 1):
            for j, count in enumerate(row, 1):
                matrix.SetBinContent(i, j, count)
                matrix.SetBinError(i, j, np.sqrt(count))
        flat = ROOT.TH1D("study_flat", "", 4, .5, 4.5)
        eta = ROOT.TH1D("study_eta", "", 2, -1., 1.)
        return make_training(matrix, histogram(flat, "reco", reco, reco),
                             histogram(flat, "gen", gen, gen),
                             FlattenedBinning(((50., 60.), (60., 70.)), 2), [eta, eta], {})

    def test_forward_fold_across_both_axes_and_root_round_trip(self):
        t = self.training([[4, 0, 0, 2], [0, 4, 0, 0], [2, 0, 6, 0], [0, 0, 0, 6]],
                          [12, 4, 8, 9], [8, 4, 12, 10])
        # Gen bin 4 migrates in both pT and eta to Reco bin 1. Explicit
        # probabilities include losses. Reco rows 1 and 4 add fakes.
        np.testing.assert_allclose(fold_expected(t, np.array([16, 8, 60, 30])),
                                   [(4/8*16 + 2/10*30)*2, 8., 2/8*16+6/12*60, 27.])
        np.testing.assert_allclose(fold_expected(t, values(t.gen)), values(t.reco))
        np.testing.assert_array_equal(values(t.matched_reco), [6, 4, 8, 6])
        np.testing.assert_array_equal(values(t.matched_gen), [6, 4, 6, 8])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "response.root"
            write_unfolding_output(path, histograms=(t.matrix, t.reco), metadata={})
            restored = load_object(str(path), t.matrix.GetName())
            for errors in (False, True):
                np.testing.assert_array_equal(matrix_values(restored, errors=errors),
                                              matrix_values(t.matrix, errors=errors))

    def test_identity_with_fixed_fake_efficiency_and_missing_support(self):
        t = self.training(np.diag([4, 4, 4, 4]), [8, 8, 8, 8], [16, 16, 16, 16])
        truth = np.array([1., 2., 3., 4.])
        expected = fold_expected(t, truth)
        for step in bayes_toy_trajectory(t, expected[None, :], 4):
            np.testing.assert_allclose(step[0], truth)
        result = run_unfolding(t, histogram(t.reco, "toy", expected, expected), iterations=4)
        np.testing.assert_allclose(values(result.unfolded.histogram), truth)
        np.testing.assert_allclose(values(result.refolded.histogram), expected)
        self.assertGreaterEqual(np.linalg.eigvalsh(result.unfolded.covariance).min(), 0.)
        np.testing.assert_allclose(variances(result.unfolded.histogram), 4*expected)

        empty = self.training(np.diag([4, 4, 4, 0]), [4, 4, 4, 0], [4, 4, 4, 4])
        with self.assertRaisesRegex(ValueError, "unmeasured Gen"):
            fold_expected(empty, truth)
        np.testing.assert_array_equal(fold_expected(empty, np.array([1., 0., 0., 0.])),
                                      [1., 0., 0., 0.])

    def test_scaled_poisson_moments_and_reproducibility(self):
        expected = np.array([4., 8., 0.])
        weight = np.array([.5, 2., .25])
        first = poisson_toys(expected, weight, 100000, np.random.default_rng(737))
        second = poisson_toys(expected, weight, 100000, np.random.default_rng(737))
        np.testing.assert_array_equal(first, second)
        np.testing.assert_allclose(first.mean(axis=0), expected, rtol=.012)
        np.testing.assert_allclose(first.var(axis=0), expected*weight, rtol=.025)
        np.testing.assert_array_equal(first[:, 2], 0)
        self.assertTrue((first >= 0).all())

    def test_shared_event_cross_covariance_cancels_identity_ratio_error(self):
        counts = [4., 9., 16., 25.]
        t = self.training(np.diag(counts), counts, counts)
        result = run_unfolding(t, t.reco, iterations=3)
        cross = np.diag(counts)
        flat, _, checks = validation_views(t, result, t.gen, reco_truth_cross=cross,
                                           unfolded_truth_cross=cross)
        for curve in flat["to_gen"].values():
            np.testing.assert_allclose(values(curve.histogram), 1.)
            np.testing.assert_allclose(curve.covariance, 0., atol=1e-15)
        for i in range(2):
            np.testing.assert_allclose(checks[f"fb_closure_pt{i}"]["residual"], 0., atol=1e-15)

    def test_bias_variance_mse(self):
        mean, bias, variance, mse = risk_metrics(np.array([[1., 6.], [5., 6.]]),
                                                 np.array([2., 8.]))
        for actual, expected in zip((mean, bias, variance, mse), ([3, 6], [1, -2], [4, 0], [5, 4])):
            np.testing.assert_array_equal(actual, expected)


if __name__ == "__main__":
    unittest.main()
