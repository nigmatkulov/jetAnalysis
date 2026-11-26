"""Analytic cross-pT migration and correction test using RooUnfold."""
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

from hist_analysis.python.notebook_setup import load_root, load_roounfold

ROOT = load_root()

from hist_analysis.python.unfolding import (
    FlattenedBinning, apply_efficiency_correction, build_roounfold_response,
    calculate_response_diagnostics, prepare_factorized_corrections, unfold_bayes,
)
from hist_analysis.python.unfolding_studies import run_correction_study


class CorrectionStudyTest(unittest.TestCase):
    def test_first_iteration_and_covariance_with_cross_pt_migrations(self):
        load_roounfold(ROOT, project_root=Path.cwd())
        ROOT.TH1.AddDirectory(False)
        matrix = ROOT.TH2D("migration", "", 4, 0.5, 4.5, 4, 0.5, 4.5)
        truth, reco, miss, fake = [
            ROOT.TH1D(name, "", 4, 0.5, 4.5) for name in ("truth", "reco", "miss", "fake")
        ]
        for hist in (matrix, truth, reco, miss, fake):
            hist.Sumw2()

        # Outer bins are in different pT and eta bins. Each matched cell and
        # explicit loss has one weighted entry with variance equal to w^2.
        counts = np.array([[60., 0., 0., 10.], [0., 30., 0., 0.],
                           [0., 0., 30., 0.], [20., 0., 0., 40.]])
        for i, j in zip(*np.nonzero(counts)):
            w = counts[i, j]
            matrix.Fill(i + 1, j + 1, w)
            truth.Fill(j + 1, w)
            reco.Fill(i + 1, w)
        for i, m, f in ((1, 20., 30.), (4, 50., 40.)):
            truth.Fill(i, m)
            miss.Fill(i, m)
            reco.Fill(i, f)
            fake.Fill(i, f)
        diagnostics = calculate_response_diagnostics(
            matrix, truth, reco, explicit_miss=miss, explicit_fake=fake,
        )
        self.assertAlmostEqual(diagnostics.effective_miss.GetBinError(4), 50.)
        self.assertAlmostEqual(diagnostics.effective_fake.GetBinError(1), 30.)
        factors = prepare_factorized_corrections(
            reco, reco, diagnostics.matched_measured, truth, diagnostics.matched_truth,
        )
        response = build_roounfold_response(
            ROOT, diagnostics.matched_truth, diagnostics.matched_measured, matrix,
            require_fakes=False, matched_only=True, scale=1e3,
        )
        matched = unfold_bayes(ROOT, response, factors.measured_signal, iterations=1,
                               handle_fakes=False)
        full, covariance = apply_efficiency_correction(
            matched.histogram, matched.covariance, factors.efficiency,
            name="full",
        )
        template = ROOT.TH1D("eta", "", 2, -1, 1)
        with TemporaryDirectory() as directory, patch(
            "hist_analysis.python.unfolding_studies.save_comparison"
        ) as plots:
            chains = ("purity_only", "efficiency_only", "uncorrected", "fully_corrected")
            results = run_correction_study(
                ROOT, response, diagnostics, factors, matched, full, covariance,
                measured=reco, training_reco=reco, training_truth=truth,
                reference_truth=truth, layout=FlattenedBinning(((50, 60), (60, 70)), 2),
                eta_templates=(template, template), iterations=1, chains=chains,
                output_dir=directory, output_tag="toy", annotations=(),
                cms_header_left="CMS Simulation", cms_header_right="pp",
                explicit_miss=miss, explicit_fake=fake,
            )
            self.assertEqual(plots.call_count, 34)
            output = ROOT.TFile.Open(str(Path(directory) / "toy_fully_corrected.root"))
            self.assertTrue(output.Get("hUnfoldedCovariance"))
            output.Close()
            with ROOT.TFile.Open(str(Path(directory) / "toy_forward_backward.root")) as output:
                for chain, (hist, cov) in results.items():
                    for pt_index in range(2):
                        back, fwd = 2 * pt_index, 2 * pt_index + 1
                        f = hist.GetBinContent(fwd + 1)
                        b = hist.GetBinContent(back + 1)
                        key = f"{chain}_fb_pt{pt_index}"
                        ratio = output.Get(key)
                        self.assertAlmostEqual(output.Get(f"{key}_forward").GetBinContent(1), f)
                        self.assertAlmostEqual(output.Get(f"{key}_backward").GetBinContent(1), b)
                        self.assertEqual(ratio.GetXaxis().GetBinLowEdge(1), 0.)
                        self.assertEqual(ratio.GetXaxis().GetBinUpEdge(1), 1.)
                        self.assertAlmostEqual(ratio.GetBinContent(1), f / b)
                        self.assertAlmostEqual(
                            ratio.GetBinError(1)**2,
                            cov[fwd, fwd] / b**2 + f**2 * cov[back, back] / b**4,
                        )

        # With the matched training prior the first-iteration Jacobian is the
        # posterior from Bayes' rule, counts[i,j] / row_sum[i].
        posterior = (counts / counts.sum(axis=1, keepdims=True)).T
        measured = np.array([100., 30., 30., 100.])
        measured_variance = np.array([4600., 900., 900., 3600.])
        purity = np.array([.7, 1., 1., .6])
        efficiency = np.array([.8, 1., 1., .5])
        for chain, pre, post in zip(chains, (purity, 1, 1, purity),
                                    (1, 1 / efficiency, 1, 1 / efficiency)):
            hist, cov = results[chain]
            jacobian = posterior * pre * np.broadcast_to(post, (4,))[:, None]
            np.testing.assert_allclose(
                [hist.GetBinContent(i) for i in range(1, 5)], jacobian @ measured,
            )
            np.testing.assert_allclose(cov, (jacobian * measured_variance) @ jacobian.T)
