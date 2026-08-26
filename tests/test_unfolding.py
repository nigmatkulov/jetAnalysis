import tempfile
import unittest
from array import array
from pathlib import Path

from hist_analysis.python.notebook_setup import load_root

ROOT = load_root(batch=True)

from hist_analysis.python.unfolding import (
    FlattenedBinning, as_pt_intervals, calculate_response_diagnostics,
    extract_eta_block, flatten_pt_eta_projections, project_eta_by_pt,
    flatten_sparse_response, unflatten_to_eta_projections,
    write_unfolding_output,
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


if __name__ == "__main__":
    unittest.main()
