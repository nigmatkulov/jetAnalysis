"""Stored-envelope validation must reject silent changes to theory point sets."""
import contextlib
import io
import unittest
from hist_analysis.python.notebook_setup import load_root
from hist_analysis.python.theory_plotting import validate_total_band

ROOT = load_root()


class TotalBandValidation(unittest.TestCase):
    def setUp(self):
        self.h = ROOT.TH1D('central_validation', '', 3, 0., 3.)
        self.h.SetDirectory(0)
        for i, y in enumerate((1., float('nan'), 2.), 1):
            self.h.SetBinContent(i, y)
        self.g = ROOT.TGraphAsymmErrors(2)
        for i, (x, y) in enumerate(((.5, 1.), (2.5, 2.))):
            self.g.SetPoint(i, x, y)
            self.g.SetPointError(i, .5, .5, .1, .2)

    def check_band(self):
        with contextlib.redirect_stdout(io.StringIO()):
            return validate_total_band(self.h, self.g, context='test')

    def test_gap_and_asymmetric_errors_are_preserved(self):
        self.assertEqual(self.check_band(), [(.5, 1., .1, .2), (2.5, 2., .1, .2)])

    def test_missing_point_is_rejected(self):
        self.g.RemovePoint(1)
        with self.assertRaisesRegex(ValueError, 'valid-point mismatch'):
            self.check_band()

    def test_central_value_mismatch_is_rejected(self):
        self.g.SetPoint(1, 2.5, 2.01)
        with self.assertRaisesRegex(ValueError, 'alignment mismatch'):
            self.check_band()

    def test_bin_edge_mismatch_is_rejected(self):
        self.g.SetPointError(1, .4, .5, .1, .2)
        with self.assertRaisesRegex(ValueError, 'alignment mismatch'):
            self.check_band()

    def test_negative_error_is_rejected(self):
        self.g.SetPointError(1, .5, .5, -.1, .2)
        with self.assertRaisesRegex(ValueError, 'invalid point'):
            self.check_band()


if __name__ == '__main__':
    unittest.main()
