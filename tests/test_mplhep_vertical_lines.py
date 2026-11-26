"""ROOT vertical-line and legend preservation in mplhep mirrors."""
import ast
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from hist_analysis.python.notebook_setup import load_root
ROOT = load_root()
from hist_analysis.python import mplhep_plotting as plotting
from hist_analysis.python.plotting import draw_overlay


class VerticalLineExportTests(unittest.TestCase):
    def test_full_export_preserves_boundary_lines_and_labels(self):
        histogram = ROOT.TH1D('line_export_input', '', 16, -4, 4)
        histogram.SetDirectory(0)
        for i in range(1, 17):
            histogram.SetBinContent(i, 1)
            histogram.SetBinError(i, .1)
        canvas = draw_overlay({'Gen': histogram}, title='', x_title='eta', y_title='yield',
                              x_range=(-3, 3), y_range=(.1, 2))
        legend = next(o for o in canvas.GetListOfPrimitives() if o.InheritsFrom('TLegend'))
        lines = []
        canvas.cd()
        for eta, style in ((2.4, 1), (1.9, 2)):
            pair = []
            for sign in (-1, 1):
                line = ROOT.TLine(sign*eta, .1, sign*eta, 2)
                line.SetLineColor(ROOT.kBlack)
                line.SetLineStyle(style)
                line.SetLineWidth(2)
                line.Draw('SAME')
                lines.append(line)
                pair.append(line)
            legend.AddEntry(pair[0], f'|#eta_{{CM}}| = {eta}', 'l')
        canvas.Update()
        old_contents = [histogram.GetBinContent(i) for i in range(histogram.GetNcells())]
        seen = []
        def inspect(fig, output, **kwargs):
            axis = fig.axes[0]
            boundaries = [line for line in axis.lines if line.get_linestyle() in ('-', '--')
                          and len(line.get_xdata()) == 2
                          and line.get_xdata()[0] == line.get_xdata()[1]]
            self.assertEqual(len(boundaries), 4)
            self.assertEqual(sorted(float(line.get_xdata()[0]) for line in boundaries), [-2.4,-1.9,1.9,2.4])
            for line in boundaries:
                self.assertEqual(line.get_linestyle(), '-' if abs(line.get_xdata()[0]) == 2.4 else '--')
                self.assertEqual(line.get_linewidth(), 2)
                self.assertEqual(list(line.get_ydata()), [0, 1])
                self.assertEqual(line.get_color(), '0.6')
                self.assertEqual(line.get_transform(), axis.get_xaxis_transform())
            labels = [text.get_text() for text in axis.get_legend().get_texts()]
            self.assertEqual(labels, ['Gen', plotting._root_label('|#eta_{CM}| = 2.4'), plotting._root_label('|#eta_{CM}| = 1.9')])
            self.assertEqual(labels.count(plotting._root_label('|#eta_{CM}| = 2.4')), 1)
            self.assertEqual(labels.count(plotting._root_label('|#eta_{CM}| = 1.9')), 1)
            self.assertEqual(tuple(axis.get_xlim()), (-3, 3))
            self.assertEqual(tuple(axis.get_ylim()), (.1, 2))
            seen.append(True)
        with tempfile.TemporaryDirectory() as directory, patch('matplotlib.figure.Figure.savefig', inspect):
            plotting.save_canvas(canvas, output=Path(directory)/'mirror.pdf')
        self.assertEqual(seen, [True])
        self.assertEqual(old_contents, [histogram.GetBinContent(i) for i in range(histogram.GetNcells())])

    def test_notebook_boundaries_export_for_every_selection_and_plot_type(self):
        notebook_path = Path(__file__).resolve().parents[1] / 'hist_analysis/notebooks/02_mc_jet_selection.ipynb'
        notebook = json.loads(notebook_path.read_text())
        source = next(''.join(cell['source']) for cell in notebook['cells']
                      if cell.get('id') == 'draw-projections')
        function = next(node for node in ast.parse(source).body
                        if isinstance(node, ast.FunctionDef) and node.name == 'add_cm_eta_boundaries')
        exported = []
        def inspect(fig, output, **kwargs):
            axis = fig.axes[0]
            lines = [line for line in axis.lines if line.get_color() == '0.6']
            self.assertEqual(sorted(float(line.get_xdata()[0]) for line in lines), [-2.4, -1.9, 1.9, 2.4])
            self.assertTrue(all(list(line.get_ydata()) == [0, 1] for line in lines))
            labels = [text.get_text() for text in axis.get_legend().get_texts()]
            self.assertIn(plotting._root_label('|#eta_{CM}| = 2.4'), labels)
            self.assertIn(plotting._root_label('|#eta_{CM}| = 1.9'), labels)
            exported.append(Path(output).stem)
        namespace = dict(ROOT=ROOT, SAVE_PNG=False,
                         save_canvas=lambda canvas, output, **kwargs:
                         plotting.save_canvas(canvas, output=output))
        exec(compile(ast.Module(body=[function], type_ignores=[]), '<notebook-boundaries>', 'exec'), namespace)
        with tempfile.TemporaryDirectory() as directory, patch('matplotlib.figure.Figure.savefig', inspect):
            for selection in ('noSel', 'trkMax', 'jetId'):
                for kind in ('distribution', 'ratio'):
                    with self.subTest(selection=selection, kind=kind):
                        histogram = ROOT.TH1D(f'boundary_{selection}_{kind}', '', 16, -4, 4)
                        histogram.SetDirectory(0)
                        for i in range(1, 17):
                            histogram.SetBinContent(i, 1)
                            histogram.SetBinError(i, .01)
                        canvas = draw_overlay({'Gen': histogram}, title='', x_title='eta', y_title='yield',
                                              x_range=(-3.6, 3.6), y_range=(.1, 2))
                        namespace['add_cm_eta_boundaries'](canvas, Path(directory) / f'{selection}_{kind}.pdf')
        self.assertEqual(len(exported), 6)

    def test_lines_span_frame_after_log_scale_and_range_changes(self):
        plt, _ = plotting._mplhep()
        canvas = ROOT.TCanvas('line_full_height_log', '', 400, 400)
        line = ROOT.TLine(1.9, .4, 1.9, .6)
        line.SetLineColor(ROOT.kRed)
        line.Draw()
        fig, axis = plt.subplots()
        try:
            plotting._draw_vertical_lines(axis, canvas)
            axis.set_yscale('log')
            axis.set_ylim(.01, 100)
            fig.canvas.draw()
            artist = axis.lines[0]
            pixels = artist.get_transform().transform(artist.get_xydata())
            self.assertAlmostEqual(pixels[0][1], axis.bbox.ymin)
            self.assertAlmostEqual(pixels[1][1], axis.bbox.ymax)
            self.assertEqual(artist.get_color(), '0.6')
            self.assertEqual(line.GetLineColor(), ROOT.kRed)
            self.assertEqual((line.GetY1(), line.GetY2()), (.4, .6))
        finally:
            plt.close(fig)

    def test_horizontal_and_ndc_lines_are_not_vertical_boundaries(self):
        plt, _ = plotting._mplhep()
        canvas = ROOT.TCanvas('line_ndc_filter', '', 400, 400)
        horizontal = ROOT.TLine(0, 1, 2, 1)
        horizontal.Draw()
        ndc = ROOT.TLine(.5, .1, .5, .9)
        ndc.SetNDC(True)
        ndc.Draw()
        fig, axis = plt.subplots()
        try:
            plotting._draw_vertical_lines(axis, canvas)
            self.assertEqual(len(axis.lines), 0)
        finally:
            plt.close(fig)


if __name__ == '__main__':
    unittest.main()
