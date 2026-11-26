"""Numerical tests for the eta-dependent same-event inference."""

from pathlib import Path
import tempfile
import unittest

import numpy as np

from hist_analysis.python.dijet_eta_unfolding import Response, iterate, same_event_fiducial
from hist_analysis.python.dijet_eta_unfolding_study import bottom_line


class EtaDijetInferenceTests(unittest.TestCase):
    def setUp(self):
        migration = np.array([[60., 10., 4.], [20., 70., 6.]])
        truth = np.array([100., 100., 25.])
        reco = migration.sum(axis=1) + np.array([2., 1.])
        fake = reco - migration.sum(axis=1)
        self.response = Response(migration, migration, truth, truth, reco, reco, fake, fake,
                                 ("feed_in", "fake"), 2)

    def test_extended_denominator_and_equivalent_fixed_models(self):
        b = self.response
        bridge = same_event_fiducial(b)
        np.testing.assert_allclose(b.efficiency, [.8, .8, .4])
        np.testing.assert_allclose(b.background(), [6., 7.])
        for model in (b, bridge):
            np.testing.assert_allclose(model.fold(model.truth[:2], background=model.background()),
                                       model.measured)
            np.testing.assert_allclose(
                iterate(model.kernel[:, :2], model.measured, model.background(), [100., 100.], 4),
                [100., 100.], rtol=1e-12,
            )

    def test_known_fixed_kernel_iterates(self):
        model = self.response
        k = model.kernel[:, :2]
        estimate = iterate(k, [70., 90.], [0., 0.], [50., 150.], 1)
        self.assertAlmostEqual(estimate[0], 70.*(30./45.)/.8 + 90.*(10./115.)/.8)
        self.assertAlmostEqual(estimate.sum(), 200.)

    def test_data_and_rate_jacobians_against_finite_differences(self):
        model = self.response
        k, y, background = model.kernel[:, :2], np.array([45., 110.]), model.background()
        start = [80., 120.]
        templates = model.components()
        _, data_jac, rate_jac = iterate(k, y, background, start, 4, derivative=True,
                                       background_templates=templates)
        step = 1e-4
        for i in range(len(y)):
            delta = np.eye(len(y))[i]*step
            expected = (iterate(k, y + delta, background, start, 4)
                        - iterate(k, y - delta, background, start, 4))/(2*step)
            np.testing.assert_allclose(data_jac[:, i], expected, rtol=1e-8, atol=1e-9)
        for i, template in enumerate(templates.T):
            delta = template*step
            expected = (iterate(k, y, background + delta, start, 4)
                        - iterate(k, y, background - delta, start, 4))/(2*step)
            np.testing.assert_allclose(rate_jac[:, i], expected, rtol=1e-8, atol=1e-9)

    def test_sparse_and_background_dominated_rows(self):
        result = iterate([[.5, 0.], [0., 0.]], [1., 0.], [10., 0.], [10., 0.], 1)
        np.testing.assert_allclose(result, [2/3, 0.])
        with self.assertRaisesRegex(ValueError, "outside"):
            iterate([[.5, 0.], [0., 0.]], [1., 1.], [10., 0.], [10., 0.], 1)

    def test_bottom_line_separate_calibration_and_test_toys(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder)
            result = bottom_line(self.response, [[100., 100.], [120., 80.]], [1., 1.],
                                 output, iterations=2, count=20, seed=17, observed=[80., 80.])
            with np.load(output / "bottom_line.npz") as archive:
                self.assertEqual(archive["reco"].shape, (2, 121))
            self.assertEqual(result["unfolded_covariance_ranks"], [2, 2])
            self.assertGreater(result["statistics"]["reco"][0]["observed_pvalue"], 0)
