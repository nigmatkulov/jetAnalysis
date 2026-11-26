"""Selection text survives ROOT-to-MPLHEP export without changing data."""
import unittest
from hist_analysis.python.notebook_setup import load_root
ROOT = load_root()
from hist_analysis.python.root_style import draw_text_block, draw_cms_header
from hist_analysis.python.mplhep_plotting import _pad_ptave_annotations


class AnnotationTests(unittest.TestCase):
    def test_all_ndc_selections_survive_and_headers_are_not_duplicated(self):
        canvas = ROOT.TCanvas('selection_annotation_test', '', 600, 600)
        canvas.SetLeftMargin(.2)
        canvas.SetRightMargin(.1)
        canvas.SetBottomMargin(.15)
        canvas.SetTopMargin(.1)
        labels = draw_text_block(canvas, ('Reco jet selection: trkMax',
                                          '-1 < #eta_{lab}^{ref} < 1',
                                          'Embedding MC'), x=.235, y=.85)
        headers = draw_cms_header(canvas, 'CMS Simulation', 'pPb (8.16 TeV)')
        annotations = _pad_ptave_annotations(canvas)
        self.assertEqual(len(annotations), 3)
        self.assertTrue(any('trkMax' in line for _, _, line in annotations))
        self.assertTrue(any('Embedding MC' in line for _, _, line in annotations))
        self.assertAlmostEqual(annotations[0][0], .05, places=6)
        self.assertAlmostEqual(annotations[0][1], .7/.75, places=6)
        canvas.Close()


    def test_comparison_legend_avoids_high_eta_points(self):
        from hist_analysis.python.data_mc_plotting import draw_comparison_figure
        h = ROOT.TH1D('legend_geometry_test', '', 10, 0, 2)
        h.SetDirectory(0)
        for i in range(1, 11):
            h.SetBinContent(i, 1. if i < 8 else 1.18)
            h.SetBinError(i, .02)
        before = [(h.GetBinContent(i), h.GetBinError(i)) for i in range(12)]
        canvas = draw_comparison_figure(
            {'MC': h}, name='legend_geometry_canvas', x_title='eta', y_title='F/B',
            x_range=(0, 2), y_range=(.75, 1.25), styles={'MC': (ROOT.kRed, 24)},
            annotations=('Reco input jet selection: jetId',), output=None)
        self.assertEqual(canvas._comparison_legend_location, 'lower right')
        self.assertEqual(before, [(h.GetBinContent(i), h.GetBinError(i)) for i in range(12)])
        canvas.Close()
        band = ROOT.TGraphAsymmErrors(1)
        band.SetPoint(0, 1.8, .80)
        band.SetPointError(0, .1, .1, .03, .03)
        canvas2 = draw_comparison_figure(
            {'MC': h}, name='legend_band_geometry_canvas', x_title='eta', y_title='F/B',
            x_range=(0, 2), y_range=(.75, 1.25), styles={'MC': (ROOT.kRed, 24)},
            annotations=('Reco input jet selection: jetId',), bands=((band, 'Theory'),), output=None)
        self.assertEqual(canvas2._comparison_legend_location, 'lower left')
        self.assertEqual(before, [(h.GetBinContent(i), h.GetBinError(i)) for i in range(12)])
        canvas2.Close()


if __name__ == '__main__':
    unittest.main()
