"""Presentation scaling must leave ROOT analysis and systematic exports intact."""
import json
from pathlib import Path
import tempfile
import unittest

from hist_analysis.python.notebook_setup import load_root
ROOT = load_root(batch=True)
from hist_analysis.python.systematic_plotting import (
    statistical_plot_copies, validate_stat_uncertainty_scale,
    draw_systematic_overlay, draw_systematic_variation_overlay, build_systematic_band,
)
from hist_analysis.python.systematic_fits import (
    fit_histogram_variations, calculate_one_sided_systematic,
    smooth_systematic_running_max, write_systematic_csv,
)


def bins(h):
    return tuple((h.GetBinContent(i), h.GetBinError(i)) for i in range(h.GetNcells()))


class StatisticalDisplayTest(unittest.TestCase):
    def test_invalid_configuration(self):
        for value in (0, -1, float('nan'), float('inf'), True, '0.5', None, 1j):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'finite positive'):
                validate_stat_uncertainty_scale(value)

    def test_fits_systematics_and_csv_invariant(self):
        original = ROOT.TH1D('stat_scaling_original', '', 12, -2, 2)
        original.SetDirectory(0)
        for i in range(original.GetNcells()):
            x = original.GetXaxis().GetBinCenter(i)
            original.SetBinContent(i, 1.02 + 0.01*x*x)
            original.SetBinError(i, 0.01 + 0.001*i)
        original_bins = bins(original)
        snapshots = []
        with tempfile.TemporaryDirectory() as directory:
            for factor in (1.0, 0.5, 2.0):
                fits, summaries = fit_histogram_variations(
                    {'variation': original}, formula='pol2', fit_range=(-2, 2),
                    name_prefix=f'stat_scale_fit_{factor}')
                systematic = calculate_one_sided_systematic(
                    original, name=f'stat_scale_syst_{factor}',
                    variation_function=fits['variation'], evaluation_range=(-2, 2))
                smoothed = smooth_systematic_running_max(
                    systematic, name=f'stat_scale_smooth_{factor}',
                    evaluation_range=(-2, 2), smoothing_origin=0)
                before = (bins(systematic), bins(smoothed))
                band = build_systematic_band(smoothed, evaluation_range=(-2, 2), name='stat_band')
                # Exercise both ROOT and MPLHEP with the real final renderer.
                for repeat in range(2):
                    display = statistical_plot_copies({'variation': original}, factor)['variation']
                    for (value, error), (old_value, old_error) in zip(bins(display), original_bins):
                        self.assertEqual(value, old_value)
                        self.assertEqual(error, factor*old_error)
                    canvas = draw_systematic_variation_overlay(
                        {'variation': original}, method='direct_ratio',
                        renderer=draw_systematic_overlay, stat_uncertainty_scale=factor,
                        title='', x_title='eta', y_title='Density',
                        canvas_name=f'stat_scale_{factor}_{repeat}',
                        output=Path(directory)/f'plot_{factor}_{repeat}.pdf')
                    # Match the notebook band stage: redraw only retained display points.
                    band.Draw('E2 SAME')
                    for point_histogram in canvas._systematic_plot_histograms.values():
                        point_histogram.Draw('E1 SAME')
                    plotted = [o for o in canvas.GetListOfPrimitives() if o.InheritsFrom('TH1')]
                    drawn = next(iter(canvas._systematic_plot_histograms.values()))
                    self.assertTrue(all(ROOT.addressof(o) == ROOT.addressof(drawn) for o in plotted))
                    for (_, error), (_, old_error) in zip(bins(drawn), original_bins):
                        self.assertEqual(error, factor*old_error)
                    self.assertEqual(bins(original), original_bins)
                    self.assertEqual((bins(systematic), bins(smoothed)), before)
                csv = write_systematic_csv(smoothed, Path(directory)/f'{factor}.csv', evaluation_range=(-2, 2))
                snapshots.append((summaries['variation']['parameters'], summaries['variation']['parameter_errors'],
                                  before, csv.read_bytes(), tuple(band.GetErrorYhigh(i) for i in range(band.GetN()))))
            self.assertEqual(snapshots[0], snapshots[1])
            self.assertEqual(snapshots[0], snapshots[2])

    def test_folded_display_only(self):
        from hist_analysis.python.systematic_methods import construct_systematic_variations
        from hist_analysis.python.systematic_plotting import draw_systematic_variation_overlay
        h = ROOT.TH1D('folded_scale_test', '', 4, 0, 4)
        h.SetDirectory(0)
        for i, value in enumerate((0.9, 1.0, 1.1, 1.2), 1):
            h.SetBinContent(i, value)
            h.SetBinError(i, 0.02)
        folded = construct_systematic_variations({'Up / Default': h}, method='folded_deviation')
        original = bins(folded['Up / Default'])
        def renderer(histograms, **kwargs):
            return next(iter(histograms.values()))
        for scale in (1.0, 0.5, 2.0):
            for repeat in range(2):
                display = draw_systematic_variation_overlay(
                    folded, method='folded_deviation', renderer=renderer,
                    stat_uncertainty_scale=scale)
                for (value, error), (old_value, old_error) in zip(bins(display), original):
                    self.assertEqual(value, old_value)
                    self.assertEqual(error, scale * old_error)
                self.assertEqual(bins(folded['Up / Default']), original)
            direct = draw_systematic_variation_overlay(
                {'Up / Default': h}, method='direct_ratio', renderer=renderer,
                stat_uncertainty_scale=scale)
            for (value, error), (old_value, old_error) in zip(bins(direct), bins(h)):
                self.assertEqual(value, old_value)
                self.assertEqual(error, scale * old_error)
            self.assertEqual(h.GetBinError(1), 0.02)

    def test_default_preserves_existing_render_path(self):
        from unittest.mock import patch
        histograms = {'original': object()}
        with patch('hist_analysis.python.plotting.draw_overlay') as renderer, \
             patch('hist_analysis.python.systematic_plotting.save_systematic_canvas'):
            canvas = draw_systematic_overlay(histograms, mplhep=False)
        self.assertIs(renderer.call_args.args[0], histograms)
        self.assertIs(canvas, renderer.return_value)
        self.assertNotIn('stat_uncertainty_scale', renderer.call_args.kwargs)

    def test_notebook_opt_in_scope_and_syntax(self):
        root = Path(__file__).resolve().parents[1]
        for name, count in [('JER', 3), ('JEU', 1), ('pileup', 1), ('pointing_resolution', 1)]:
            nb = json.loads((root/f'hist_analysis/notebooks/04_systematics_{name}.ipynb').read_text())
            source = '\n'.join(''.join(c['source']) for c in nb['cells'] if c['cell_type'] == 'code')
            self.assertIn('STAT_UNCERTAINTY_SCALE = ', source)
            self.assertNotIn('eta_shapes, stat_uncertainty_scale=', source)
            self.assertNotIn('fb_ratios, stat_uncertainty_scale=', source)
            self.assertIn('canvas._systematic_plot_histograms.values()', source)
            self.assertNotIn('in ratios.values():', source)
            self.assertNotIn("if SYSTEMATIC_METHOD == 'folded_deviation' else 1.0", source)
            self.assertEqual(source.count('stat_uncertainty_scale=STAT_UNCERTAINTY_SCALE'), count)
            for i, cell in enumerate(nb['cells']):
                if cell['cell_type'] == 'code':
                    code = '\n'.join(line for line in ''.join(cell['source']).splitlines() if not line.startswith('%'))
                    compile(code, f'{name}:cell{i}', 'exec')


if __name__ == '__main__':
    unittest.main()
