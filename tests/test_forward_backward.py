"""Physics checks for folding full eta yields, weighted errors and legacy files."""
from array import array
from pathlib import Path
import tempfile
import unittest

from hist_analysis.python.notebook_setup import load_root
ROOT = load_root()
from hist_analysis.python.forward_backward import fold_eta_histogram, forward_backward_ratio
from hist_analysis.python.dijet_closures import DijetClosureCurve, build_dijet_gen_comparisons
from hist_analysis.python.histogram_production import produce_source


class ForwardBackwardTest(unittest.TestCase):
    def test_weighted_signed_yields_match_direct_filling_away_from_edges(self):
        full = ROOT.TH1D('weighted_full', '', 72, -3.6, 3.6)
        expected = [ROOT.TH1D('expected_'+s, '', 36, 0., 3.6) for s in ('f', 'b')]
        for h in [full, *expected]:
            h.SetDirectory(0)
            h.Sumw2()
        for eta, weight in [(0.15, 2.), (0.15, -0.5), (-0.15, 3.),
                            (-0.15, 4.), (1.25, 2.5), (-1.25, -1.)]:
            full.Fill(eta, weight)
            expected[eta < 0].Fill(abs(eta), weight)
        folded = fold_eta_histogram(full, name='weighted_fold')
        for actual, reference in zip(folded, expected):
            for i in range(38):
                self.assertAlmostEqual(actual.GetBinContent(i), reference.GetBinContent(i))
                self.assertAlmostEqual(actual.GetBinError(i), reference.GetBinError(i))
        ratio = forward_backward_ratio(*folded, name='weighted_ratio')
        i = ratio.FindBin(.15)
        self.assertAlmostEqual(ratio.GetBinContent(i), 1.5/7.)
        self.assertAlmostEqual(ratio.GetBinError(i)**2, 4.25/49. + 1.5**2*25./7.**4)
        self.assertEqual(ratio.GetBinContent(ratio.FindBin(1.05)), 0.)
        self.assertEqual(ratio.GetBinError(ratio.FindBin(1.05)), 0.)

    def test_exact_edges_follow_full_bins_instead_of_legacy_abs_eta(self):
        full = ROOT.TH1D('edge_full', '', 72, -3.6, 3.6)
        old = ROOT.TH1D('edge_old_b', '', 36, 0., 3.6)
        full.Sumw2(); old.Sumw2()
        full.Fill(-1.5, 2.5); old.Fill(1.5, 2.5)
        _, backward = fold_eta_histogram(full, name='edge_fold')
        source = full.FindBin(-1.5)
        mapped = 37-source
        self.assertEqual(backward.GetBinContent(mapped), 2.5)
        self.assertEqual(backward.GetBinError(mapped), 2.5)
        self.assertNotEqual(mapped, old.FindBin(1.5))
        self.assertEqual(backward.Integral(), old.Integral())

    def test_2d_flow_bins_and_cpp_agree(self):
        full = ROOT.TH2D('two_full', ';pT;eta', 2, array('d', [0., 50., 100.]),
                         4, array('d', [-2., -1., 0., 1., 2.]))
        full.Sumw2()
        for pt, eta, weight in [(-1., -.25, 2.), (150., .25, 3.),
                                (25., -3., 4.), (25., 3., 5.)]:
            full.Fill(pt, eta, weight)
        forward, backward = fold_eta_histogram(full, name='flow_fold')
        self.assertEqual(backward.GetBinContent(0, 1), 2.)
        self.assertEqual(forward.GetBinContent(3, 1), 3.)
        self.assertEqual(backward.GetBinContent(1, 3), 4.)
        self.assertEqual(forward.GetBinContent(1, 3), 5.)
        header = Path(__file__).resolve().parents[1]/'hist_analysis/cpp/eta_folding.h'
        ROOT.gInterpreter.Declare(f'#include "{header}"')
        ROOT.gInterpreter.Declare('TH1* testFoldCpp(const TH1& h, bool b) { return eta_folding::fold(h,b,"cpp_fold").release(); }')
        for side, expected in [(False, forward), (True, backward)]:
            actual = ROOT.testFoldCpp(full, side)
            ROOT.SetOwnership(actual, True)
            for i in range(expected.GetNcells()):
                self.assertEqual(actual.GetBinContent(i), expected.GetBinContent(i))
                self.assertEqual(actual.GetBinError(i), expected.GetBinError(i))

    def test_rejects_straddling_or_asymmetric_bins(self):
        for bins in [[-2., -.5, .5, 2.], [-2., -1., 0., .5, 2.]]:
            full = ROOT.TH1D('invalid_full', '', len(bins)-1, array('d', bins))
            with self.assertRaises(ValueError): fold_eta_histogram(full, name='invalid')

    def test_closure_uses_full_only_and_producer_discards_legacy_sides(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/'source.root'
            output = Path(directory)/'produced_histograms/bundle.root'
            f = ROOT.TFile.Open(str(source), 'RECREATE')
            full = ROOT.TH2D('hGenDijetPtEtaCM_4', ';pT;eta', 100, 0., 1000., 72, -3.6, 3.6)
            full.Sumw2(); full.Fill(65., .15, 2.); full.Fill(65., -.15, 4.)
            full.Write()
            for side in ['Forward', 'Backward']:
                legacy = ROOT.TH2D('hGenDijetPtEta'+side+'_4', '', 100, 0., 1000., 36, 0., 3.6)
                legacy.Fill(65., .15, 999.)
                legacy.Write()
            f.Close()
            curve = DijetClosureCurve('Gen', 'hGenDijetPtEtaCM_{eta_cut_index}')
            args = dict(eta_cut_index=4, ptave_range=(60, 70), rebin_eta=2, normalization='integral')
            before = build_dijet_gen_comparisons(source, (curve,), **args)
            ratio = before[1]['Gen']
            self.assertEqual(ratio.GetBinContent(ratio.FindBin(.15)), .5)
            produce_source(dict(source=source, output=output, category='official_mc', dataset='pythia', direction='combined', trigger=None))
            check = ROOT.TFile.Open(str(output))
            self.assertFalse(any('Forward' in k.GetName() or 'Backward' in k.GetName() for k in check.GetListOfKeys()))
            check.Close()
            after = build_dijet_gen_comparisons(output, (curve,), **args)
            for a, b in zip(before[:2], after[:2]):
                for i in range(a['Gen'].GetNcells()):
                    self.assertAlmostEqual(a['Gen'].GetBinContent(i), b['Gen'].GetBinContent(i))
                    self.assertAlmostEqual(a['Gen'].GetBinError(i), b['Gen'].GetBinError(i))
            with self.assertRaises(ValueError):
                build_dijet_gen_comparisons(source, (curve,), **args, ratio_option='B')


if __name__ == '__main__': unittest.main()
