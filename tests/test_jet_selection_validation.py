"""Analytical expectations independent of ROOT's division implementation."""
import math
import unittest
from array import array
from hist_analysis.python.notebook_setup import load_root
ROOT = load_root()
from hist_analysis.python.jet_selection import (
    project_selection_eta, independent_shape_ratio, ratio_variance,
)


class SelectionValidation(unittest.TestCase):
    def histogram(self, name, values, variances):
        h = ROOT.TH1D(name, '', len(values), 0, len(values))
        h.SetDirectory(0)
        h.Sumw2()
        for i, (v, var) in enumerate(zip(values, variances), 1):
            h.SetBinContent(i, v)
            h.SetBinError(i, math.sqrt(var))
        return h

    def test_raw_correlated_and_independent(self):
        # Weighted subset: successes 2,3; failure 4. Cov(N,D)=13.
        expected = ((1-5/9)**2*13 + (5/9)**2*16)/81
        self.assertAlmostEqual(ratio_variance(5, 9, 13, 29, 13), expected)
        self.assertAlmostEqual(ratio_variance(5, 9, 13, 29), 13/81+25*29/9**4)
        self.assertEqual(ratio_variance(9, 9, 29, 29, 29), 0)
        self.assertAlmostEqual(ratio_variance(4, 10, 4, 10, 4), .4*.6/10)
        with self.assertRaises(ValueError): ratio_variance(0, 0, 0, 0)
        with self.assertRaises(ValueError): ratio_variance(1, 1, 1, 1, 2)

    def test_random_normalization_analytic_derivatives(self):
        n = self.histogram('shape_n', [2, 3], [4, 9])
        d = self.histogram('shape_d', [4, 6], [16, 36])
        # r=(n1/(n1+n2))/(d1/(d1+d2)); direct partial derivatives.
        # dr/d(n1,n2,d1,d2)=(.3,-.2,-.15,.1)
        expected = .3**2*4 + .2**2*9 + .15**2*16 + .1**2*36
        for mode in ('integral', 'bin_width'):
            r = independent_shape_ratio(n, d, mode, (0, 2))
            self.assertAlmostEqual(r.GetBinContent(1), 1)
            self.assertAlmostEqual(r.GetBinError(1)**2, expected)
        one = independent_shape_ratio(n, d, 'integral', (0, 1))
        self.assertAlmostEqual(one.GetBinError(1), 0)
        # A bin outside the normalization window includes sum uncertainty.
        self.assertAlmostEqual(one.GetBinError(2)**2, 4.)

    def test_empty_bins(self):
        n = self.histogram('empty_n', [0, 2], [0, 2])
        d = self.histogram('empty_d', [0, 4], [0, 4])
        r = independent_shape_ratio(n, d, 'integral', (0, 2))
        self.assertEqual(r.GetBinContent(1), 0)
        self.assertEqual(r.GetBinError(1), 0)

    def test_projection_boundaries_flows_sumw2_and_transpose(self):
        for transpose in (False, True):
            args = (2, -1., 1., 3, 0., 300.) if transpose else (3, 0., 300., 2, -1., 1.)
            h = ROOT.TH2D(f'projection_{transpose}', '', *args)
            h.SetDirectory(0)
            h.Sumw2()
            for pt, eta, weight in [(0,-.5,2),(100,-.5,3),(199,.5,4),
                                     (200,.5,5),(-1,-.5,7),(301,-.5,8),
                                     (150,-2,9),(150,2,10)]:
                h.Fill(eta,pt,weight) if transpose else h.Fill(pt,eta,weight)
            h.GetXaxis().SetRange(1,1)
            r = project_selection_eta(h, (100,200), name=f'projected_{transpose}')
            self.assertEqual([r.GetBinContent(i) for i in range(4)], [9,3,4,10])
            self.assertEqual([round(r.GetBinError(i)**2) for i in range(4)], [81,9,16,100])
            self.assertEqual(h.GetXaxis().GetLast(), 1)
            rebinned = project_selection_eta(h,(100,200),rebin=2)
            self.assertEqual(rebinned.GetBinContent(1),7)
            self.assertAlmostEqual(rebinned.GetBinError(1)**2,25)
            for bounds in ((101,200),(200,100),(400,500),(100,100)):
                with self.assertRaises(ValueError): project_selection_eta(h,bounds)
            with self.assertRaises(ValueError): project_selection_eta(h,(100,200),rebin=3)

    def test_unweighted_variable_eta_width(self):
        h=ROOT.TH2D('variable_eta','',2,array('d',[0,100,200]),
                    2,array('d',[-2,0,1]))
        h.SetDirectory(0)
        for pt,eta in [(100,-1),(100,-1),(100,.5),(200,.5)]: h.Fill(pt,eta)
        r=project_selection_eta(h,(100,200))
        self.assertEqual(r.GetBinContent(1),2)
        self.assertAlmostEqual(r.GetBinError(1)**2,2)
        self.assertEqual(r.GetBinContent(2),1)


if __name__ == '__main__': unittest.main()
