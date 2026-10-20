"""One-sided estimators and shared mathematical plot labels."""
import unittest
from hist_analysis.python.notebook_setup import load_root
ROOT = load_root(batch=True)
from hist_analysis.python.systematic_methods import construct_systematic_variations
from hist_analysis.python.jer_systematic_methods import construct_jer_variations
from hist_analysis.python.systematic_fits import calculate_one_sided_systematic
from hist_analysis.python.systematic_plotting import draw_systematic_variation_overlay


class SharedSystematicMethodsTest(unittest.TestCase):
    def test_one_sided_raw_invariance_and_compatibility(self):
        self.assertIs(construct_systematic_variations, construct_jer_variations)
        h = ROOT.TH1D('one_sided_methods_test', '', 6, 0, 6)
        h.SetDirectory(0)
        for i, value in enumerate((0.95, 0.98, 1, 1.02, 1.05), 1):
            h.SetBinContent(i, value)
            h.SetBinError(i, 0.01)
        folded = construct_systematic_variations({'Reco / Ref': h}, method='folded_deviation')['Reco / Ref']
        direct_s = calculate_one_sided_systematic(h, name='one_sided_direct_test')
        folded_s = calculate_one_sided_systematic(folded, name='one_sided_folded_test')
        for i in range(1, 7):
            self.assertAlmostEqual(direct_s.GetBinContent(i), folded_s.GetBinContent(i), places=15)
        self.assertEqual(folded.GetBinContent(6), 0)

    def test_renderer_uses_matching_labels_colors_and_fits(self):
        def renderer(histograms, **kwargs):
            return histograms, kwargs
        values = {'Reco / Ref': object()}
        options = dict(y_title='Reco / Ref', overlay_functions={'Reco / Ref': object()},
                       style_indices={'Reco / Ref': 0},
                       overlay_text={'Reco / Ref': ('Reco / Ref: p0=1', 'chi2/ndf=0/1')})
        direct, direct_options = draw_systematic_variation_overlay(values, method='direct_ratio', renderer=renderer, **options)
        self.assertIs(direct, values)
        self.assertEqual(direct_options, options)
        folded, folded_options = draw_systematic_variation_overlay(values, method='folded_deviation', renderer=renderer, **options)
        label = '1 + |1 - Reco / Ref|'
        self.assertEqual(set(folded), {label})
        self.assertEqual(set(folded_options['overlay_functions']), {label})
        self.assertEqual(folded_options['style_indices'], {label: 0})
        self.assertEqual(folded_options['y_title'], label)
        self.assertTrue(folded_options['overlay_text'][label][0].startswith('Folded Reco / Ref:'))
        self.assertEqual(options['y_title'], 'Reco / Ref')
        with self.assertRaises(ValueError):
            draw_systematic_variation_overlay(values, method='unsupported', renderer=renderer)


if __name__ == '__main__':
    unittest.main()
