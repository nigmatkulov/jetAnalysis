"""Validate the consumer of stored totals, including metadata and ratio scaling."""
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from hist_analysis.python.notebook_setup import load_root
ROOT = load_root()
from hist_analysis.python.data_mc_plotting import load_data_total, data_total_band, selection_lines, decorate_pad, set_acceptance_frame, accepted_eta_range, check_distinct_generator_inputs, check_distinct_generator_predictions, load_comparison_source, data_over_mc, draw_root_systematic_brackets


class StoredDataTotalTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        folder = self.base / 'full'; folder.mkdir()
        self.stem = folder / 'example_ptave_50_60'
        self.record = dict(source='JEU', observable='full', frame='CM', direction='combined',
                           selection='jetId', binning='normal', eta_cut=1.9, units='fraction',
                           convention='symmetric_relative', target_sample='pPb 8.16 TeV dijets',
                           pt_low=50, pt_high=60, production_sample='data', csv='minimumbias_jeu.csv')
        self.manifest = Path(str(self.stem)+'_provenance.json')
        self.manifest.write_text(json.dumps([self.record]))
        np.savetxt(str(self.stem)+'_total.csv', [[-.5,1,.5,.02],[.5,1,.5,.03]], delimiter=',')
        Path(str(self.stem)+'.csv').write_text('x_low,x_high,total\n-1,0,.02\n0,1,.03\n')

    def load(self, **kwargs):
        options = dict(observable='full', direction='combined', selection='jetId', binning='normal',
                       eta_cut=1.9, sample='mb')
        options.update(kwargs)
        return load_data_total(self.base, (50,60), **options)

    def test_fraction_conversion_and_statistical_errors_unchanged(self):
        h=ROOT.TH1D('consumer_test','',2,-1,1);h.SetDirectory(0)
        for i,y in enumerate((2.,4.),1):h.SetBinContent(i,y);h.SetBinError(i,.4)
        before=[(h.GetBinContent(i),h.GetBinError(i)) for i in (1,2)]
        g=data_total_band(h,self.load(),name='test_band')
        self.assertAlmostEqual(g.GetErrorYlow(0),.04)
        self.assertAlmostEqual(g.GetErrorYhigh(1),.12)
        self.assertEqual(before,[(h.GetBinContent(i),h.GetBinError(i)) for i in (1,2)])
        # Ratios use their own central y, never absolute-distribution errors.
        h.SetBinContent(1,1.2)
        g=data_total_band(h,self.load(),name='test_ratio_band')
        self.assertAlmostEqual(g.GetErrorYlow(0),.024)

    def test_wrong_trigger_or_configuration_fails(self):
        for kwargs in [dict(sample='jet60'),dict(binning='test'),dict(direction='pgoing'),dict(observable='fb')]:
            with self.assertRaisesRegex(ValueError,'expected one total'):self.load(**kwargs)

    def test_ambiguous_variant_fails(self):
        (self.base/'full'/'other_ptave_50_60_provenance.json').write_text(self.manifest.read_text())
        with self.assertRaisesRegex(ValueError,'found 2'):self.load()

    def test_percent_metadata_is_not_treated_as_fraction(self):
        self.record['units']='percent';self.manifest.write_text(json.dumps([self.record]))
        with self.assertRaises(ValueError):self.load()

    def test_audit_disagreement_fails(self):
        Path(str(self.stem)+'.csv').write_text('x_low,x_high,total\n-1,0,.2\n0,1,.03\n')
        with self.assertRaisesRegex(ValueError,'alignment mismatch'):self.load()

    def test_missing_systematic_bin_fails_without_interpolation(self):
        h=ROOT.TH1D('consumer_mismatch','',2,-.9,1.1);h.SetDirectory(0);h.SetBinContent(1,1)
        with self.assertRaisesRegex(ValueError,'no unique systematic bin'):
            data_total_band(h,self.load(),name='bad_band')

    def test_missing_data_point_is_omitted(self):
        h=ROOT.TH1D('consumer_gap','',2,-1,1);h.SetDirectory(0)
        h.SetBinContent(1,float('nan'));h.SetBinContent(2,1)
        self.assertEqual(data_total_band(h,self.load(),name='gap_band').GetN(),1)

    def test_exact_acceptance_frame_preserves_original_bin_edges(self):
        h = ROOT.TH1D('acceptance_test', '', 20, -2, 2)
        h.SetDirectory(0)
        h.SetBinContent(1, 3.)
        h.SetBinError(1, .2)
        before = [(h.GetBinContent(i), h.GetBinError(i)) for i in range(22)]
        c = ROOT.TCanvas('acceptance_canvas', '', 600, 600)
        h.Draw('E1'); c.Update()
        frame = set_acceptance_frame(c, h, (-1.9, 1.9))
        self.assertAlmostEqual(c.GetUxmin(), -1.9)
        self.assertAlmostEqual(c.GetUxmax(), 1.9)
        self.assertEqual(h.GetXaxis().GetXmin(), -2.)
        self.assertEqual(before, [(h.GetBinContent(i), h.GetBinError(i)) for i in range(22)])
        c.Close()

    def test_systematic_caps_follow_sample_color_and_toggle_only_representation(self):
        h = ROOT.TH1D('caps_test', '', 2, -1, 1)
        h.SetDirectory(0)
        h.SetBinContent(1, 2.); h.SetBinError(1, .4)
        h.SetMarkerColor(ROOT.kRed+1)
        before = (h.GetBinContent(1), h.GetBinError(1))
        for enabled in (True, False):
            c = ROOT.TCanvas('caps_canvas_' + str(enabled), '', 600, 600)
            h.Draw('E1')
            objects = decorate_pad(c, {'ratio': h},
                totals=(('ratio', self.load()),) if enabled else ())
            graphs = [o for o in objects if o.InheritsFrom('TGraphAsymmErrors')]
            self.assertEqual(len(graphs), int(enabled))
            if graphs:
                self.assertEqual(graphs[0].GetLineColor(), ROOT.kRed+1)
                self.assertEqual(graphs[0].GetFillStyle(), 0)
                self.assertAlmostEqual(graphs[0].GetErrorYhigh(0), .04)
            self.assertEqual(before, (h.GetBinContent(1), h.GetBinError(1)))
            c.Close()

    def test_jet_acceptance_retains_straddling_bins_and_exact_edges(self):
        h = ROOT.TH1D('edge_acceptance', '', 36, -3.6, 3.6)
        h.SetDirectory(0)
        for cut, expected in ((1.9, 2.0), (1.8, 1.8), (1.5, 1.6)):
            low, high = accepted_eta_range(h, cut)
            self.assertAlmostEqual(low, -expected)
            self.assertAlmostEqual(high, expected)
        fb = ROOT.TH1D('edge_acceptance_fb', '', 18, 0, 3.6)
        fb.SetDirectory(0)
        low, high = accepted_eta_range(fb, 1.9, observable='fb')
        self.assertEqual(low, 0.)
        self.assertAlmostEqual(high, 2.)

    def test_duplicate_generator_sources_and_predictions_require_decision(self):
        with self.assertRaisesRegex(ValueError, 'aliases'):
            check_distinct_generator_inputs((('official', self.base/'same.root', 'gen'),
                                            ('optional', self.base/'same.root', 'gen')), 'official')
        h = ROOT.TH1D('duplicate_gen', '', 2, -1, 1)
        h.SetDirectory(0); h.SetBinContent(1, 1.)
        copy = h.Clone('duplicate_gen_copy'); copy.SetDirectory(0)
        with self.assertRaisesRegex(ValueError, 'numerically identical'):
            check_distinct_generator_predictions({'official': h, 'optional': copy},
                                                  'official', ('optional',))
        copy.SetBinContent(1, 2.)
        check_distinct_generator_predictions({'official': h, 'optional': copy},
                                             'official', ('optional',))

    def test_cached_sources_survive_file_close_and_projection(self):
        path = self.base / 'source.root'
        f = ROOT.TFile(str(path), 'RECREATE')
        h = ROOT.TH2D('source', '', 2, 0, 2, 2, -1, 1)
        h.SetBinContent(1, 1, 3.); h.SetBinError(1, 1, .4)
        h.Write(); f.Close()
        cache = {}
        source = load_comparison_source(path, 'source', cache)
        projection = source.ProjectionY('cached_projection', 1, 1, 'e')
        projection.SetDirectory(0); projection.Scale(2.)
        self.assertIs(source, load_comparison_source(path, 'source', cache))
        self.assertFalse(source.GetDirectory())
        self.assertEqual(source.GetBinContent(1, 1), 3.)
        self.assertAlmostEqual(source.GetBinError(1, 1), .4)
        self.assertEqual(len(cache), 1)

    def test_bracket_paths_keep_uncertainty_endpoints_and_color(self):
        h=ROOT.TH1D('bracket_frame_test','',2,-1,1); h.SetDirectory(0)
        h.SetMaximum(4.); h.SetMinimum(0.)
        c=ROOT.TCanvas('bracket_geometry_test','',800,800); h.Draw(); c.Update()
        g=ROOT.TGraphAsymmErrors(1)
        g.SetPoint(0,0.,2.); g.SetPointError(0,.1,.1,.3,.4); g.SetLineColor(ROOT.kBlue+1)
        before=(g.GetPointY(0),g.GetErrorYlow(0),g.GetErrorYhigh(0))
        paths=draw_root_systematic_brackets(c,g)
        self.assertEqual(len(paths),2)
        for path,endpoint in zip(paths,(1.7,2.4)):
            self.assertEqual(path.GetN(),4)
            self.assertAlmostEqual(path.GetY()[1],endpoint)
            self.assertAlmostEqual(path.GetY()[2],endpoint)
            self.assertNotEqual(path.GetY()[0],endpoint)
            self.assertEqual(path.GetLineColor(),ROOT.kBlue+1)
        self.assertEqual(before,(g.GetPointY(0),g.GetErrorYlow(0),g.GetErrorYhigh(0)))
        c.Close()

    def test_data_over_mc_keeps_independent_statistical_errors(self):
        d=ROOT.TH1D('d_ratio_test','',2,-1,1); d.SetDirectory(0)
        m=ROOT.TH1D('m_ratio_test','',2,-1,1); m.SetDirectory(0)
        d.SetBinContent(1,6.); d.SetBinError(1,.6)
        m.SetBinContent(1,3.); m.SetBinError(1,.3)
        d.SetBinContent(2,2.)  # zero MC denominator stays numerical zero, omitted in display
        before=(d.GetBinContent(1),d.GetBinError(1),m.GetBinContent(1),m.GetBinError(1))
        r=data_over_mc(d,m,name='data_over_mc_test',style=(ROOT.kRed+1,20))
        self.assertEqual(r.GetBinContent(1),2.)
        self.assertAlmostEqual(r.GetBinError(1),np.sqrt(.6**2/3**2+6**2*.3**2/3**4))
        self.assertEqual(r.GetMarkerColor(),ROOT.kRed+1)
        self.assertEqual(r._valid_ratio_bins,(1,))
        self.assertEqual(before,(d.GetBinContent(1),d.GetBinError(1),m.GetBinContent(1),m.GetBinError(1)))

    def test_selection_uses_current_contract(self):
        root=Path(__file__).resolve().parents[1]
        lines=selection_lines(root,(50,60),eta_cut=1.9,direction='combined',selection='jetId')
        self.assertIn('p_{T,1} #geq 50, p_{T,2} #geq 40 GeV',lines)
        self.assertIn('2#pi/3',lines[2])
        self.assertEqual(lines[0], '#it{anti}-k_{T} R = 0.4')
        self.assertFalse(any(token in ' '.join(lines) for token in ('combined','jetId','Jet60','Minimum bias')))


if __name__=='__main__':unittest.main()
