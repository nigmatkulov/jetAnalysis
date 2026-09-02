import tempfile
import unittest
from array import array
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from hist_analysis.python.notebook_setup import load_root

ROOT = load_root(batch=True)

from hist_analysis.python.unfolding import (
    apply_efficiency_correction,
    FlattenedBinning, RegularizationDistributionWriter,
    RegularizationScanResult, ResponseDiagnostics,
    RooUnfoldResponseBundle,
    as_pt_intervals, calculate_response_diagnostics,
    calculate_closure_metrics,
    calculate_toy_metrics, copy_histogram_errors, histogram_fingerprint,
    load_regularization_scan_cache, load_toy_cache, make_gaussian_toys,
    regularization_distribution_file_matches,
    match_pt_block_yields, match_truth_to_reco_pt_yields,
    scale_pt_blocks,
    extract_eta_block, flatten_pt_eta_projections, project_eta_by_pt,
    flatten_sparse_response, scan_bayes_regularization,
    unflatten_to_eta_projections,
    prepare_factorized_corrections, restrict_histogram_range, validate_response_accounting,
    write_regularization_scan_cache, write_toy_cache, write_unfolding_output,
)


class UnfoldingHelpersTest(unittest.TestCase):
    def setUp(self):
        ROOT.TH1.AddDirectory(False)

    def test_interval_forms_and_global_mapping(self):
        self.assertEqual(as_pt_intervals((0, 40, 80)), ((0.0, 40.0), (40.0, 80.0)))
        self.assertEqual(as_pt_intervals(((50, 60), (80, 100))), ((50.0, 60.0), (80.0, 100.0)))
        layout = FlattenedBinning(((0.0, 40.0), (40.0, 80.0)), 3)
        self.assertEqual(layout.global_bin(1, 2), 5)
        self.assertEqual(layout.indices(5), (1, 2))

    def test_project_flatten_and_unflatten_round_trip(self):
        source = ROOT.TH2D("source", "", 4, 0, 80, 3, -1.5, 1.5)
        source.Sumw2()
        for x_bin in range(1, 5):
            for y_bin in range(1, 4):
                source.SetBinContent(x_bin, y_bin, 10 * x_bin + y_bin)
                source.SetBinError(x_bin, y_bin, 0.1 * (10 * x_bin + y_bin))
        projections = project_eta_by_pt(source, (0, 40, 80), name_prefix="eta")
        flattened, layout = flatten_pt_eta_projections(
            projections, name="flat", pt_bins=(0, 40, 80),
        )
        restored = unflatten_to_eta_projections(
            flattened, projections, layout, name_prefix="restored",
        )
        for original, result in zip(projections, restored):
            for eta_bin in range(1, 4):
                self.assertEqual(result.GetBinContent(eta_bin), original.GetBinContent(eta_bin))
                self.assertEqual(result.GetBinError(eta_bin), original.GetBinError(eta_bin))
        block = extract_eta_block(flattened, projections[1], layout, 1, name="block")
        self.assertEqual(block.GetBinContent(2), projections[1].GetBinContent(2))

    def test_eta_restriction_is_shared_by_projection_and_response(self):
        source = ROOT.TH2D("restricted_source", "", 2, 0, 80, 6, -3, 3)
        source.Sumw2()
        for x_bin in range(1, 3):
            for y_bin in range(1, 7):
                source.SetBinContent(x_bin, y_bin, 10 * x_bin + y_bin)
        projections = project_eta_by_pt(
            source, (0, 40, 80), name_prefix="restricted_eta", eta_range=(-2, 2),
        )
        self.assertEqual(projections[0].GetNbinsX(), 4)
        self.assertEqual(projections[0].GetXaxis().GetXmin(), -2.0)
        self.assertEqual(projections[0].GetXaxis().GetXmax(), 2.0)

        bins = array("i", [2, 6, 2, 6])
        low = array("d", [0.0, -3.0, 0.0, -3.0])
        high = array("d", [80.0, 3.0, 80.0, 3.0])
        sparse = ROOT.THnSparseD("restricted_sparse", "", 4, bins, low, high)
        sparse.Fill(array("d", [20.0, -1.5, 60.0, 1.5]), 3.0)
        response, layout = flatten_sparse_response(
            sparse, (0, 40, 80), eta_range=(-2, 2),
        )
        self.assertEqual(layout.n_eta_bins, 4)
        self.assertEqual(response.GetNbinsX(), 8)
        self.assertEqual(response.GetBinContent(layout.global_bin(1, 4),
                                                layout.global_bin(0, 1)), 3.0)

    def test_eta_cut_keeps_rebinned_boundary_bin(self):
        # Width 0.2 gives a final accepted bin [1.8, 2.0] centered at 1.9.
        source = ROOT.TH1D("rebinned_eta", "", 22, -2.2, 2.2)
        restricted = restrict_histogram_range(
            source, (-1.9, 1.9), name="restricted_rebinned_eta",
        )
        self.assertEqual(restricted.GetNbinsX(), 20)
        self.assertAlmostEqual(restricted.GetXaxis().GetXmin(), -2.0)
        self.assertAlmostEqual(restricted.GetXaxis().GetXmax(), 2.0)
        self.assertAlmostEqual(restricted.GetXaxis().GetBinCenter(1), -1.9)
        self.assertAlmostEqual(restricted.GetXaxis().GetBinCenter(20), 1.9)

    def test_response_diagnostics(self):
        truth = ROOT.TH1D("truth", "", 2, 0.5, 2.5)
        measured = ROOT.TH1D("measured", "", 2, 0.5, 2.5)
        response = ROOT.TH2D("matrix", "", 2, 0.5, 2.5, 2, 0.5, 2.5)
        miss = ROOT.TH1D("miss", "", 2, 0.5, 2.5)
        fake = ROOT.TH1D("fake", "", 2, 0.5, 2.5)
        response.SetBinContent(1, 1, 8.0)
        response.SetBinContent(2, 2, 15.0)
        truth.SetBinContent(1, 10.0)
        truth.SetBinContent(2, 20.0)
        measured.SetBinContent(1, 9.0)
        measured.SetBinContent(2, 18.0)
        miss.SetBinContent(1, 1.0)
        miss.SetBinContent(2, 4.0)
        fake.SetBinContent(1, 0.5)
        fake.SetBinContent(2, 2.0)
        diagnostics = calculate_response_diagnostics(
            response, truth, measured, explicit_miss=miss, explicit_fake=fake,
        )
        self.assertEqual(diagnostics.effective_miss.GetBinContent(1), 2.0)
        self.assertEqual(diagnostics.effective_fake.GetBinContent(2), 3.0)
        self.assertEqual(diagnostics.boundary_miss.GetBinContent(1), 1.0)
        self.assertEqual(diagnostics.boundary_fake.GetBinContent(2), 1.0)
        accounting = validate_response_accounting(
            response, truth, measured, diagnostics,
        )
        self.assertEqual(accounting.max_truth_residual, 0.0)
        self.assertEqual(accounting.max_measured_residual, 0.0)
        self.assertAlmostEqual(accounting.global_efficiency, 23.0 / 30.0)
        self.assertAlmostEqual(accounting.global_fake_fraction, 4.0 / 27.0)

    def test_factorized_fake_and_efficiency_corrections(self):
        measured = ROOT.TH1D("factor_measured", "", 2, 0.5, 2.5)
        training_measured = ROOT.TH1D("factor_training_measured", "", 2, 0.5, 2.5)
        matched_measured = ROOT.TH1D("factor_matched_measured", "", 2, 0.5, 2.5)
        training_truth = ROOT.TH1D("factor_training_truth", "", 2, 0.5, 2.5)
        matched_truth = ROOT.TH1D("factor_matched_truth", "", 2, 0.5, 2.5)
        for histogram, values in (
            (measured, (120.0, 90.0)),
            (training_measured, (100.0, 100.0)),
            (matched_measured, (80.0, 90.0)),
            (training_truth, (100.0, 80.0)),
            (matched_truth, (50.0, 60.0)),
        ):
            for index, value in enumerate(values, 1):
                histogram.SetBinContent(index, value)
        measured.SetBinError(1, 12.0)
        corrections = prepare_factorized_corrections(
            measured, training_measured, matched_measured,
            training_truth, matched_truth,
        )
        self.assertAlmostEqual(corrections.purity.GetBinContent(1), 0.8)
        self.assertAlmostEqual(corrections.measured_signal.GetBinContent(1), 96.0)
        self.assertAlmostEqual(corrections.measured_signal.GetBinError(1), 9.6)
        self.assertAlmostEqual(corrections.efficiency.GetBinContent(2), 0.75)

        unfolded_matched = ROOT.TH1D("factor_unfolded_matched", "", 2, 0.5, 2.5)
        unfolded_matched.SetBinContent(1, 50.0)
        unfolded_matched.SetBinError(1, 5.0)
        unfolded_matched.SetBinContent(2, 60.0)
        covariance = ROOT.TMatrixD(2, 2)
        covariance[0][0] = 25.0
        covariance[0][1] = covariance[1][0] = 3.0
        covariance[1][1] = 16.0
        corrected, corrected_covariance = apply_efficiency_correction(
            unfolded_matched, covariance, corrections.efficiency,
        )
        self.assertAlmostEqual(corrected.GetBinContent(1), 100.0)
        self.assertAlmostEqual(corrected.GetBinError(1), 10.0)
        self.assertAlmostEqual(corrected_covariance[0][0], 100.0)
        self.assertAlmostEqual(corrected_covariance[0][1], 8.0)

    def test_closure_metrics_ignore_tiny_bins_only_in_relative_mean(self):
        reference = ROOT.TH1D("metric_reference", "", 3, 0.5, 3.5)
        estimate = ROOT.TH1D("metric_estimate", "", 3, 0.5, 3.5)
        for index, (expected, observed) in enumerate(
            ((100.0, 90.0), (50.0, 55.0), (0.01, 1.0)), 1
        ):
            reference.SetBinContent(index, expected)
            estimate.SetBinContent(index, observed)
        metrics = calculate_closure_metrics(
            estimate, reference, populated_fraction=1.0e-3,
        )
        self.assertAlmostEqual(metrics.integral_ratio, 146.0 / 150.01)
        self.assertAlmostEqual(metrics.relative_l1, 15.99 / 150.01)
        self.assertAlmostEqual(metrics.mean_absolute_relative, 0.1)
        self.assertEqual(metrics.compared_bins, 2)

    def test_sparse_response_keeps_off_diagonal_pt_migration(self):
        bins = array("i", [2, 2, 2, 2])
        low = array("d", [0.0, -1.0, 0.0, -1.0])
        high = array("d", [80.0, 1.0, 80.0, 1.0])
        sparse = ROOT.THnSparseD("sparse", "", 4, bins, low, high)
        # Gen in the first pT block and reco in the second pT block.
        sparse.Fill(array("d", [20.0, -0.5, 60.0, 0.5]), 3.0)
        response, layout = flatten_sparse_response(sparse, (0, 40, 80))
        global_gen = layout.global_bin(0, 1)
        global_reco = layout.global_bin(1, 2)
        self.assertEqual(response.GetBinContent(global_reco, global_gen), 3.0)
        self.assertEqual(sparse.GetAxis(0).GetFirst(), 1)
        self.assertEqual(sparse.GetAxis(0).GetLast(), 2)
        self.assertEqual(sparse.GetAxis(2).GetFirst(), 1)
        self.assertEqual(sparse.GetAxis(2).GetLast(), 2)

    def test_write_output_metadata_as_json(self):
        histogram = ROOT.TH1D("saved", "", 1, 0, 1)
        histogram.SetBinContent(1, 2.0)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.root"
            write_unfolding_output(path, histograms=[histogram], metadata={"iterations": 4})
            root_file = ROOT.TFile.Open(str(path))
            try:
                self.assertEqual(root_file.Get("saved").GetBinContent(1), 2.0)
                self.assertEqual(root_file.Get("unfoldingConfiguration").GetTitle(), '{"iterations":4}')
            finally:
                root_file.Close()

    def test_match_pt_yields_and_preserve_data_errors(self):
        layout = FlattenedBinning(((0.0, 10.0), (10.0, 20.0)), 2)
        source = ROOT.TH1D("yield_source", "", 4, 0.5, 4.5)
        target = ROOT.TH1D("yield_target", "", 4, 0.5, 4.5)
        for index, (source_value, target_value, error) in enumerate(
            ((1.0, 4.0, 0.4), (1.0, 6.0, 0.6),
             (2.0, 3.0, 0.3), (2.0, 5.0, 0.5)), 1
        ):
            source.SetBinContent(index, source_value)
            target.SetBinContent(index, target_value)
            target.SetBinError(index, error)
        matched, factors = match_pt_block_yields(source, target, layout)
        self.assertEqual(factors, (5.0, 2.0))
        self.assertEqual(matched.Integral(1, 2), target.Integral(1, 2))
        self.assertEqual(matched.Integral(3, 4), target.Integral(3, 4))
        self.assertEqual(matched.GetBinError(4), 0.5)
        scaled = scale_pt_blocks(source, layout, factors, name="scaled_truth")
        self.assertEqual(scaled.GetBinContent(1), 5.0)
        self.assertEqual(scaled.GetBinContent(3), 4.0)

    def test_truth_yield_matching_folds_the_scaled_truth(self):
        layout = FlattenedBinning(((0.0, 10.0), (10.0, 20.0)), 1)
        truth = ROOT.TH1D("iterative_truth", "", 2, 0.5, 2.5)
        target = ROOT.TH1D("iterative_target", "", 2, 0.5, 2.5)
        truth.SetBinContent(1, 10.0)
        truth.SetBinContent(2, 10.0)
        target.SetBinContent(1, 20.0)
        target.SetBinContent(2, 5.0)

        class IdentityResponse:
            @staticmethod
            def ApplyToTruth(histogram):
                return histogram.Clone("identity_fold")

        zero = truth.Clone("iterative_zero")
        zero.Reset("ICES")
        diagnostics = ResponseDiagnostics(
            truth.Clone("iterative_matched_truth"),
            truth.Clone("iterative_matched_reco"), zero, zero.Clone(), None, None,
        )
        bundle = RooUnfoldResponseBundle(
            IdentityResponse(), truth.Clone(), truth.Clone(), None, 1.0,
        )
        matched = match_truth_to_reco_pt_yields(
            bundle, truth, target, diagnostics, layout,
            tolerance=1.0e-10, relaxation=1.0,
        )
        self.assertAlmostEqual(matched.truth.GetBinContent(1), 20.0)
        self.assertAlmostEqual(matched.truth.GetBinContent(2), 5.0)
        self.assertAlmostEqual(matched.folded.GetBinContent(1), 20.0)
        self.assertAlmostEqual(matched.folded.GetBinContent(2), 5.0)

        target.SetBinError(1, 3.0)
        expectation = copy_histogram_errors(
            matched.folded, target, name="iterative_expectation",
        )
        self.assertAlmostEqual(expectation.GetBinContent(1), 20.0)
        self.assertAlmostEqual(expectation.GetBinError(1), 3.0)

    def test_gaussian_toys_are_seeded_and_keep_errors(self):
        expectation = ROOT.TH1D("toy_expectation", "", 2, 0.5, 2.5)
        expectation.SetBinContent(1, 100.0)
        expectation.SetBinError(1, 2.0)
        expectation.SetBinContent(2, 50.0)
        expectation.SetBinError(2, 3.0)
        left = make_gaussian_toys(expectation, n_toys=3, random_seed=17)
        right = make_gaussian_toys(expectation, n_toys=3, random_seed=17)
        for left_toy, right_toy in zip(left, right):
            self.assertEqual(left_toy.GetBinContent(1), right_toy.GetBinContent(1))
            self.assertEqual(left_toy.GetBinError(1), 2.0)
            self.assertEqual(left_toy.GetBinError(2), 3.0)

    def test_seeded_intermediate_caches_validate_metadata(self):
        expectation = ROOT.TH1D("cache_expectation", "", 2, 0.5, 2.5)
        expectation.SetBinContent(1, 10.0)
        expectation.SetBinContent(2, 20.0)
        toys = make_gaussian_toys(expectation, n_toys=2, random_seed=7)
        metadata = {"seed": 7, "expectation": histogram_fingerprint(expectation)}
        with tempfile.TemporaryDirectory() as directory:
            toy_path = Path(directory) / "toys_seed_7.root"
            write_toy_cache(toy_path, toys, metadata)
            loaded_toys = load_toy_cache(
                toy_path, n_toys=2, metadata=metadata,
            )
            self.assertIsNotNone(loaded_toys)
            self.assertEqual(loaded_toys[1].GetBinContent(2), toys[1].GetBinContent(2))
            self.assertIsNone(load_toy_cache(
                toy_path, n_toys=2, metadata={**metadata, "seed": 8},
            ))

            metric = calculate_toy_metrics(toys, expectation, iterations=1)
            absolute = ROOT.TH2D(
                "hRegularizationAbsoluteSummary", "", 1, 0.5, 1.5, 3, 0.5, 3.5,
            )
            relative = ROOT.TH2D(
                "hRegularizationRelativeSummary", "", 1, 0.5, 1.5, 3, 0.5, 3.5,
            )
            for index, key in enumerate(("bias_squared", "variance", "mse"), 1):
                absolute.SetBinContent(1, index, metric.absolute[key])
                relative.SetBinContent(1, index, metric.relative[key])
            scan = RegularizationScanResult((metric,), absolute, relative, 1)
            scan_path = Path(directory) / "scan_seed_7.root"
            write_regularization_scan_cache(scan_path, scan, metadata)
            loaded_scan = load_regularization_scan_cache(
                scan_path, max_iterations=1, metadata=metadata,
            )
            self.assertIsNotNone(loaded_scan)
            self.assertEqual(loaded_scan.selected_iterations, 1)
            self.assertAlmostEqual(
                loaded_scan.metrics[0].relative["mse"], metric.relative["mse"],
            )

            distributions_path = Path(directory) / "distributions_seed_7.root"
            writer = RegularizationDistributionWriter(
                distributions_path,
                histogram_groups={"flattened_inputs": [expectation], "toys": toys},
                metadata=metadata,
            )
            writer.write_unfolded(1, 0, toys[0])
            writer.write_unfolded(1, 1, toys[1])
            writer.finalize(scan)
            self.assertTrue(regularization_distribution_file_matches(
                distributions_path, metadata=metadata,
                max_iterations=1, n_toys=2,
            ))
            self.assertFalse(regularization_distribution_file_matches(
                distributions_path, metadata={**metadata, "seed": 8},
                max_iterations=1, n_toys=2,
            ))
            root_file = ROOT.TFile.Open(str(distributions_path))
            try:
                self.assertTrue(root_file.Get("flattened_inputs/cache_expectation"))
                self.assertTrue(root_file.Get(
                    "iteration_001/unfolded_toys/hUnfoldedToy_000001"
                ))
                self.assertTrue(root_file.Get(
                    "iteration_001/metrics/hRegularizationMSE_iter1"
                ))
            finally:
                root_file.Close()

    def test_negative_gaussian_draw_is_reported(self):
        expectation = ROOT.TH1D("negative_expectation", "", 1, 0.5, 1.5)
        expectation.SetBinContent(1, 0.0)
        expectation.SetBinError(1, 1.0)
        with self.assertRaisesRegex(ValueError, "Negative Gaussian toy content"):
            make_gaussian_toys(expectation, n_toys=2, random_seed=5)

    def test_toy_metrics_use_per_bin_bias_and_unbiased_variance(self):
        truth = ROOT.TH1D("metric_truth", "", 2, 0.5, 2.5)
        truth.SetBinContent(1, 10.0)
        truth.SetBinContent(2, 10.0)
        toy_a = truth.Clone("metric_toy_a")
        toy_b = truth.Clone("metric_toy_b")
        # Mean biases are +1 and -1: they must not cancel before squaring.
        toy_a.SetBinContent(1, 10.0)
        toy_b.SetBinContent(1, 12.0)
        toy_a.SetBinContent(2, 8.0)
        toy_b.SetBinContent(2, 10.0)
        result = calculate_toy_metrics([toy_a, toy_b], truth, iterations=1)
        self.assertAlmostEqual(result.absolute["bias_squared"], 1.0)
        self.assertAlmostEqual(result.absolute["variance"], 2.0)
        self.assertAlmostEqual(result.absolute["mse"], 3.0)
        self.assertAlmostEqual(result.relative["mse"], 0.03)

    def test_absolute_metric_preserves_larger_absolute_bin_contribution(self):
        truth = ROOT.TH1D("weighted_metric_truth", "", 2, 0.5, 2.5)
        truth.SetBinContent(1, 1000.0)
        truth.SetBinContent(2, 100.0)
        toy_a = truth.Clone("weighted_metric_toy_a")
        toy_b = truth.Clone("weighted_metric_toy_b")
        # Both bins are biased upward by 10%, despite their unequal contents.
        for toy in (toy_a, toy_b):
            toy.SetBinContent(1, 1100.0)
            toy.SetBinContent(2, 110.0)
        result = calculate_toy_metrics([toy_a, toy_b], truth, iterations=1)
        self.assertAlmostEqual(result.absolute["bias_squared"], 5050.0)
        self.assertAlmostEqual(result.relative["bias_squared"], 0.01)

    def test_regularization_selection_mode_is_validated(self):
        with self.assertRaisesRegex(ValueError, "selection_mode"):
            scan_bayes_regularization(
                None, None, (), None, max_iterations=1,
                selection_mode="unsupported",
            )

    def test_both_selection_modes_record_their_own_mse_minima(self):
        truth = ROOT.TH1D("selection_truth", "", 2, 0.5, 2.5)
        truth.SetBinContent(1, 1000.0)
        truth.SetBinContent(2, 10.0)
        toys = [truth.Clone("selection_toy_0"), truth.Clone("selection_toy_1")]

        def fake_unfold(*args, **kwargs):
            iteration = int(kwargs["name"].rsplit("iter", 1)[1])
            histogram = truth.Clone(f"selection_unfolded_{iteration}")
            # Iteration 1 misses only the small bin; iteration 2 misses only the
            # populated bin. Absolute and fractional MSE therefore prefer 1 and 2.
            histogram.SetBinContent(1, 1000.0 if iteration == 1 else 1100.0)
            histogram.SetBinContent(2, 20.0 if iteration == 1 else 10.0)
            return SimpleNamespace(histogram=histogram)

        with patch("hist_analysis.python.unfolding.unfold_bayes", fake_unfold):
            scan = scan_bayes_regularization(
                None, None, toys, truth, max_iterations=2,
                selection_mode="both",
            )
        self.assertEqual(scan.absolute_selected_iterations, 1)
        self.assertEqual(scan.relative_selected_iterations, 2)
        self.assertEqual(scan.selected_iterations, 1)


if __name__ == "__main__":
    unittest.main()
