"""Empty 2D cells reveal the background without changing stored yields."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from hist_analysis.python.notebook_setup import load_root
ROOT = load_root()
from hist_analysis.python import mplhep_plotting as plotting
from hist_analysis.python.histogram_io import save_histograms


class Empty2DBinTests(unittest.TestCase):
    def fixture(self):
        h = ROOT.TH2D('empty_2d_fixture', '', 3, 0, 3, 2, 0, 2)
        h.SetDirectory(0)
        h.Sumw2()
        h.SetBinContent(1, 1, .0001)  # Nonzero, below the default color minimum.
        h.SetBinError(1, 1, .0002)
        h.SetBinContent(2, 2, 2.)
        h.SetBinError(2, 2, .5)
        h.SetBinContent(0, 1, 7.)  # Flow bins are storage, not displayed cells.
        h.SetBinError(0, 1, 3.)
        return h

    def snapshot(self, h):
        return [(h.GetBinContent(i), h.GetBinError(i)) for i in range(h.GetNcells())]

    def test_root_roundtrip_preserves_zero_bins_errors_and_flows(self):
        h = self.fixture()
        expected = self.snapshot(h)
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'fixture.root'
            save_histograms(str(path), [h])
            source = ROOT.TFile.Open(str(path))
            try:
                self.assertEqual(self.snapshot(source.Get(h.GetName())), expected)
            finally:
                source.Close()

    def test_root_negative_range_empty_cell_reveals_background(self):
        from PIL import Image
        h = self.fixture()
        h.SetBinContent(3, 2, -1.)
        canvas = ROOT.TCanvas('empty_native_negative', '', 500, 500)
        canvas.SetFrameFillColor(ROOT.kMagenta)
        h.Draw('COLZ1')
        canvas.Update()
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'native.png'
            canvas.SaveAs(str(path))
            with Image.open(path) as image:
                pixel = image.convert('RGB').getpixel(
                    (canvas.XtoPixel(1.5), canvas.YtoPixel(.5)))
                self.assertEqual(pixel, (255, 0, 255))
        canvas.Close()

    def test_jes_jer_export_has_transparent_zero_faces(self):
        from matplotlib.figure import Figure
        from hist_analysis.python.jes_jer_plotting import save_jes_jer_canvas
        h = self.fixture()
        expected = self.snapshot(h)
        values, _, _ = plotting._th2_arrays(h)
        empty = (values.T == 0).ravel()
        canvas = ROOT.TCanvas('empty_jes_jer_canvas', '', 500, 500)
        h.Draw('COLZ1')
        canvas.Update()
        original_save = Figure.savefig
        def checked_save(fig, *args, **kwargs):
            fig.axes[0].set_facecolor('#ff00ff')
            fig.canvas.draw()
            colors = fig.axes[0].collections[0].get_facecolors()
            np.testing.assert_array_equal(colors[empty, 3], 0.)
            np.testing.assert_array_equal(colors[~empty, 3], 1.)
            return original_save(fig, *args, **kwargs)
        with tempfile.TemporaryDirectory() as d:
            output = Path(d) / 'jes_jer.pdf'
            with patch.object(Figure, 'savefig', checked_save):
                save_jes_jer_canvas(canvas, output=output, save_png=True)
            self.assertTrue(output.is_file())
            self.assertTrue(output.with_suffix('.png').is_file())
        canvas.Close()
        self.assertEqual(self.snapshot(h), expected)

    def test_both_mplhep_paths_have_transparent_zero_faces(self):
        h = self.fixture()
        expected = self.snapshot(h)
        values, _, _ = plotting._th2_arrays(h)
        empty = (values.T == 0).ravel()
        original_save = plotting._save
        with tempfile.TemporaryDirectory() as d:
            for renderer in ('standalone', 'canvas'):
                def inspect_and_save(fig, output, **kwargs):
                    axis = fig.axes[0]
                    axis.set_facecolor('#ff00ff')
                    fig.canvas.draw()
                    rgba = axis.collections[0].get_facecolors()
                    np.testing.assert_array_equal(rgba[empty, 3], 0.)
                    np.testing.assert_array_equal(rgba[~empty, 3], 1.)
                    if renderer == 'standalone':
                        self.assertEqual(axis.collections[0].norm.vmin, .001)
                    original_save(fig, output, **kwargs)
                output = Path(d) / (renderer + '.pdf')
                with patch.object(plotting, '_save', inspect_and_save):
                    if renderer == 'standalone':
                        plotting.save_2d(h, output=output, x_label='x', y_label='y', save_png=True)
                    else:
                        c = ROOT.TCanvas('empty_2d_canvas', '', 500, 500)
                        h.Draw('COLZ1')
                        c.Update()
                        plotting.save_canvas(c, output=output, save_png=True)
                        c.Close()
                self.assertTrue(output.is_file())
                self.assertTrue(output.with_suffix('.png').is_file())
                self.assertEqual(self.snapshot(h), expected)


if __name__ == '__main__':
    unittest.main()
