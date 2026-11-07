"""Whole-bin alignment, weighted F/B, and protected theory selection."""
import math
import unittest
from unittest.mock import patch
from hist_analysis.python.notebook_setup import load_root
ROOT = load_root()
from hist_analysis.python.pdf_cp5_comparison import prepare_pdf_cp5, load_pdf_central


class PdfCp5Tests(unittest.TestCase):
    def fixtures(self):
        pdf = ROOT.TH1D('pdf_fixture', '', 4, -2, 2)
        pdf.SetDirectory(0); pdf.Sumw2()
        for i, (value, error) in enumerate(((0, 0), (10, 3), (15, 4), (2, 1)), 1):
            pdf.SetBinContent(i, value); pdf.SetBinError(i, error)
        cp5 = ROOT.TH2D('cp5_fixture', '', 4, 40, 80, 8, -2, 2)
        cp5.SetDirectory(0); cp5.Sumw2()
        cp5.GetXaxis().SetTitle('p_{T}^{ave}'); cp5.GetYaxis().SetTitle('#eta_{CM}')
        cp5.Fill(55, -.75, 3); cp5.Fill(55, -.25, 4)
        cp5.Fill(55, .25, 6); cp5.Fill(55, .75, 8)
        cp5.Fill(65, .25, 100)  # excluded upper pTave interval
        return pdf, cp5

    def test_raw_weighted_alignment_and_zero_denominator(self):
        pdf, cp5 = self.fixtures()
        with patch('hist_analysis.python.pdf_cp5_comparison.load_pdf_central', return_value=pdf):
            result = prepare_pdf_cp5('CT18', (50, 60), cp5)
        h = result['raw']['PYTHIA CP5']
        self.assertEqual(h.Integral(), 21)
        self.assertEqual(h.GetBinContent(2), 7)
        self.assertEqual(h.GetBinError(2), 5)
        ratio = result['fb']['PYTHIA CP5']
        self.assertEqual(ratio.GetBinContent(1), 2)
        self.assertAlmostEqual(ratio.GetBinError(1), math.sqrt(200)/7)
        self.assertEqual(ratio.GetBinContent(2), 0)
        self.assertEqual(ratio.GetBinError(2), 0)
        self.assertEqual(ratio._valid_ratio_bins, (1,))
        self.assertEqual(pdf.Integral(), 27)  # input not mutated
        self.assertEqual(result['shapes']['CT18'].Integral(), 1)
        full_ratio = result['comparison_ratios']['shapes']
        cp5_shape, pdf_shape = result['shapes']['PYTHIA CP5'], result['shapes']['CT18']
        self.assertAlmostEqual(full_ratio.GetBinContent(3), (14/21)/(15/27))
        self.assertAlmostEqual(full_ratio.GetBinError(3), math.hypot(
            cp5_shape.GetBinError(3)/pdf_shape.GetBinContent(3),
            cp5_shape.GetBinContent(3)*pdf_shape.GetBinError(3)/pdf_shape.GetBinContent(3)**2))
        self.assertNotIn(1, full_ratio._valid_ratio_bins)
        fb_ratio = result['comparison_ratios']['fb']
        self.assertAlmostEqual(fb_ratio.GetBinContent(1), 2/1.5)
        theory_fb = result['fb']['CT18']
        self.assertAlmostEqual(fb_ratio.GetBinError(1), math.hypot(
            ratio.GetBinError(1)/1.5, 2*theory_fb.GetBinError(1)/1.5**2))
        self.assertEqual(fb_ratio._valid_ratio_bins, (1,))

    def test_undefined_cp5_fb_is_masked_even_with_nonzero_pdf(self):
        pdf, cp5 = self.fixtures()
        pdf.SetBinContent(1, 2)  # PDF F/B outer bin is now defined
        with patch('hist_analysis.python.pdf_cp5_comparison.load_pdf_central', return_value=pdf):
            result = prepare_pdf_cp5('CT18', (50, 60), cp5)
        ratio = result['comparison_ratios']['fb']
        self.assertEqual(ratio.GetBinContent(2), 0)  # ROOT numerical convention
        self.assertNotIn(2, ratio._valid_ratio_bins)  # undefined CP5 input

    def test_density_normalization_and_exact_ptave(self):
        pdf, cp5 = self.fixtures()
        with patch('hist_analysis.python.pdf_cp5_comparison.load_pdf_central', return_value=pdf):
            result = prepare_pdf_cp5('CT18', (50, 60), cp5, normalization='bin_width')
            self.assertAlmostEqual(result['shapes']['PYTHIA CP5'].Integral('width'), 1)
            with self.assertRaisesRegex(ValueError, 'exactly cover'):
                prepare_pdf_cp5('CT18', (51, 60), cp5)
            with self.assertRaisesRegex(ValueError, 'normalization'):
                prepare_pdf_cp5('CT18', (50, 60), cp5, normalization='none')

    def test_unavailable_theory_selections_fail(self):
        for family, interval, cut in [('EPPS21', (50, 60), 1.9),
                                      ('CT18', (90, 100), 1.9),
                                      ('nCTEQ15HQ', (50, 60), 2.4)]:
            with self.assertRaises(ValueError):
                load_pdf_central(family, interval, eta_cut=cut)

    def test_selected_key_is_propagated_and_recorded(self):
        pdf, cp5 = self.fixtures()
        for key in ('pdf_central', 'pdf_central_mirrored'):
            with patch('hist_analysis.python.pdf_cp5_comparison.load_pdf_central', return_value=pdf) as load:
                result = prepare_pdf_cp5('CT18', (50, 60), cp5, pdf_key=key)
            self.assertEqual(load.call_args.kwargs['pdf_key'], key)
            self.assertEqual(result['validation']['pdf_key'], key)
        with patch('hist_analysis.python.pdf_cp5_comparison.load_pdf_central', return_value=pdf) as load:
            prepare_pdf_cp5('CT18', (50, 60), cp5)
        self.assertEqual(load.call_args.kwargs['pdf_key'], 'pdf_central_mirrored')
        with self.assertRaisesRegex(ValueError, 'central key'):
            load_pdf_central('CT18', (50, 60), pdf_key='npdf_central')
