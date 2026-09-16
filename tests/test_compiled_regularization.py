import unittest
from pathlib import Path

from hist_analysis.python.notebook_setup import load_root, load_roounfold

ROOT = load_root(batch=True)
PROJECT_ROOT = Path(__file__).resolve().parents[1]
load_roounfold(ROOT, project_root=PROJECT_ROOT)
import RooUnfold

from hist_analysis.python.compiled_regularization import (
    resolve_compiled_scan_library,
    scan_compiled_bayes_moments,
)
from hist_analysis.python.unfolding import build_roounfold_response, unfold_bayes


class CompiledRegularizationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            resolve_compiled_scan_library(PROJECT_ROOT)
        except FileNotFoundError as error:
            raise unittest.SkipTest(str(error))
        ROOT.TH1.AddDirectory(False)

    def test_one_pass_moments_match_roounfold_at_every_iteration(self):
        truth = ROOT.TH1D("compiledTruth", "", 2, 0.5, 2.5)
        measured = ROOT.TH1D("compiledMeasured", "", 2, 0.5, 2.5)
        matrix = ROOT.TH2D("compiledMatrix", "", 2, 0.5, 2.5, 2, 0.5, 2.5)
        for measured_bin, truth_bin, value in (
            (1, 1, 120.0), (2, 1, 30.0),
            (1, 2, 10.0), (2, 2, 40.0),
        ):
            matrix.SetBinContent(measured_bin, truth_bin, value)
        truth_projection = matrix.ProjectionY("compiledTruthProjection")
        measured_projection = matrix.ProjectionX("compiledMeasuredProjection")
        for bin_index in (1, 2):
            truth.SetBinContent(bin_index, truth_projection.GetBinContent(bin_index))
            measured.SetBinContent(
                bin_index, measured_projection.GetBinContent(bin_index)
            )
        efficiency = truth.Clone("compiledEfficiency")
        efficiency.SetBinContent(1, 0.8)
        efficiency.SetBinContent(2, 0.6)
        toys = []
        for toy_index, values in enumerate(((100.0, 75.0), (140.0, 45.0),
                                             (90.0, 95.0))):
            toy = measured.Clone(f"compiledToy{toy_index}")
            for bin_index, value in enumerate(values, 1):
                toy.SetBinContent(bin_index, value)
            toys.append(toy)

        bundle = build_roounfold_response(
            RooUnfold, truth, measured, matrix, scale=1.0, require_fakes=False,
        )
        compiled = scan_compiled_bayes_moments(
            matrix, truth, efficiency, toys, max_iterations=5,
            project_root=PROJECT_ROOT, worker_count=2,
        )
        single_worker = scan_compiled_bayes_moments(
            matrix, truth, efficiency, toys, max_iterations=5,
            project_root=PROJECT_ROOT, worker_count=1,
        )
        for iteration in range(1, 6):
            expected_sums = [0.0, 0.0]
            expected_squares = [0.0, 0.0]
            for toy_index, toy in enumerate(toys):
                unfolded = unfold_bayes(
                    RooUnfold, bundle, toy, iterations=iteration,
                    handle_fakes=False, calculate_covariance=False,
                    name=f"compiledReference_iter{iteration}_toy{toy_index}",
                ).histogram
                for truth_bin in range(2):
                    value = (unfolded.GetBinContent(truth_bin + 1)
                             / efficiency.GetBinContent(truth_bin + 1))
                    expected_sums[truth_bin] += value
                    expected_squares[truth_bin] += value * value
            for truth_bin in range(2):
                self.assertAlmostEqual(
                    compiled.sums[iteration - 1][truth_bin],
                    expected_sums[truth_bin], places=11,
                )
                self.assertAlmostEqual(
                    compiled.sums_squared[iteration - 1][truth_bin],
                    expected_squares[truth_bin], places=9,
                )
                self.assertAlmostEqual(
                    compiled.sums[iteration - 1][truth_bin],
                    single_worker.sums[iteration - 1][truth_bin], places=12,
                )


if __name__ == "__main__":
    unittest.main()
