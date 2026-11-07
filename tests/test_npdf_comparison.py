"""nPDF configuration, statistical-band identity, and systematic ownership."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from hist_analysis.python.notebook_setup import load_root
ROOT = load_root()
from hist_analysis.config.histograms import (
    NPDF_TRIGGER_BINS, SYSTEMATIC_PTAVE_BINS, with_npdf_trigger_bins,
)
from hist_analysis.python.npdf_comparison import (
    statistical_band, prediction_registry, select_predictions, draw_npdf_comparisons,
    load_prediction,
)


class NpdfComparisonTests(unittest.TestCase):
    def test_exact_trigger_extension(self):
        original = {'mb': [(50, 60)], 'jet60': [(90, 100)],
                    'jet80': [(100, 110)], 'jet100': [(160, 180)]}
        extended = with_npdf_trigger_bins(original)
        self.assertEqual(original['jet60'], [(90, 100)])
        self.assertEqual(extended['jet80'], original['jet80'])
        for sample, intervals in NPDF_TRIGGER_BINS.items():
            self.assertTrue(set(intervals) <= set(extended[sample]))
            self.assertTrue(set(intervals) <= set(SYSTEMATIC_PTAVE_BINS))
        self.assertEqual(NPDF_TRIGGER_BINS['jet60'], ((80,90),(90,120)))
        self.assertEqual(NPDF_TRIGGER_BINS['jet100'], ((120,160),(200,300)))

    def test_unavailable_prediction(self):
        with tempfile.TemporaryDirectory() as directory:
            registry = prediction_registry(directory)
            with self.assertWarnsRegex(UserWarning, 'ncteq15hq'):
                self.assertEqual(select_predictions(registry, ('ncteq15hq',)), {})
            with self.assertRaises(ValueError):
                select_predictions(registry, ('unknown',))

    def test_statistical_band_and_data_only_systematics(self):
        h = ROOT.TH1D('npdf_stat_fixture', '', 2, -1, 1); h.SetDirectory(0)
        f = ROOT.TH1D('npdf_fb_fixture', '', 1, 0, 1); f.SetDirectory(0)
        for i in (1,2):h.SetBinContent(i, 1); h.SetBinError(i, .04)
        f.SetBinContent(1, 1.2); f.SetBinError(1, .15)
        graph = statistical_band(f, name='fixture', color=2)
        self.assertEqual(graph.GetPointY(0), 1.2)
        self.assertEqual(graph.GetErrorYlow(0), .15)
        self.assertEqual(graph.GetErrorYhigh(0), .15)
        self.assertEqual(graph.GetMarkerSize(), 0)
        entry = dict(sample='jet60', ptave_range=(90,120), result=dict(
            histograms={'Data':h}, ratios={'CP5':h, 'CUETP8M1':h},
            forward_backward={'Data':f, 'CP5':f, 'CUETP8M1':f}))
        for models in (('cp5',), ('cue',), ('cp5','cue')):
            for enabled in (False,True):
                total = {'fraction':np.array([.05])} if enabled else None
                calls=[]
                config=dict(DATA_LABEL='Data', ETA_CUT_VALUE=1., COMPARISON_MODELS=models,
                    MODEL_REGISTRY={'cp5':{'label':'CP5'}, 'cue':{'label':'CUETP8M1'}},
                    CURVE_STYLES={'Data':(1,21), 'CP5':(2,20), 'CUETP8M1':(4,24)},
                    data_totals=lambda *a:total, display_label=lambda x:x,
                    selection_lines=lambda *a,**k:(), PROJECT_ROOT=Path.cwd(),
                    DATA_DIRECTION='combined', DATA_SELECTION='jetId',
                    NPDF_Y_RANGES={'full':(.7,1.3), 'fb':(.7,1.3)},
                    DRAW_GRID=True, SAVE_PDF=False, SAVE_PNG=False)
                def draw(plotted,**kwargs):calls.append((plotted,kwargs));return kwargs
                with patch('hist_analysis.python.npdf_comparison.draw_comparison_figure',draw):
                    draw_npdf_comparisons(entry,{},config)
                full,fb=calls
                self.assertEqual(len(full[0]),len(models))
                self.assertEqual(tuple(fb[0]),('Data',))
                self.assertEqual(len(fb[1]['bands']),len(models))
                self.assertEqual(tuple(k for k,t in fb[1]['totals']), ('Data',) if enabled else ())
                self.assertEqual(len(full[1]['totals']),len(models) if enabled else 0)
                for band,label in fb[1]['bands']:
                    self.assertEqual(band.GetPointY(0),1.2)
                    self.assertEqual(band.GetErrorYlow(0),.15)
                    self.assertEqual(band.GetMarkerSize(),0)

    def test_supplied_predictions_match_final_root_totals(self):
        import math
        registry = prediction_registry(Path(__file__).resolve().parents[1])
        for selected in (('epps21_ct18',), ('ncteq15hq',), tuple(registry)):
            self.assertEqual(tuple(select_predictions(registry, selected)), selected)
        edges = {}
        for model, definition in registry.items():
            source = ROOT.TFile.Open(str(definition['path']), 'READ')
            try:
                self.assertEqual(source.Get('pdf_confidence_level').GetVal(), 68)
                for interval in definition['intervals']:
                    directory = f'pt{interval[0]:03d}_{interval[1]:03d}_eta1p9'
                    for observable, prefix in (('full', ''), ('fb', 'fb_')):
                        central, band = load_prediction(definition, interval, observable, 1.9)
                        original = source.Get(directory+'/'+definition['keys'][observable][0])
                        down = source.Get(directory+'/'+prefix+'total_down')
                        up = source.Get(directory+'/'+prefix+'total_up')
                        finite = [i for i in range(1, central.GetNbinsX()+1)
                                  if math.isfinite(central.GetBinContent(i))]
                        for point, i in enumerate(finite):
                            self.assertEqual(central.GetBinContent(i), original.GetBinContent(i))
                            self.assertEqual(band.GetPointY(point), original.GetBinContent(i))
                            self.assertEqual(band.GetErrorYlow(point), down.GetBinContent(i))
                            self.assertEqual(band.GetErrorYhigh(point), up.GetBinContent(i))
                        axis = central.GetXaxis()
                        key = (interval, observable)
                        current = tuple(axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins()+2))
                        if key in edges:
                            self.assertEqual(current, edges[key])
                        edges[key] = current
            finally:
                source.Close()

    def test_standalone_projection_cache_is_additive(self):
        from hist_analysis.python.histogram_production import (
            produce_source, standalone_projection_key, load_cached_eta_projection,
            ensure_npdf_projections,
        )
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/'source.root'; output = Path(directory)/'produced_histograms/cache.root'
            file = ROOT.TFile.Open(str(source),'RECREATE')
            raw = ROOT.TH2D('hGenDijetPtEtaCM','',100,0,1000,4,-2,2); raw.Sumw2()
            raw.Fill(90.,.5,3.); raw.Fill(119.,.5,4.); raw.Fill(120.,.5,12.)
            raw.Write();file.Close()
            entry=dict(source=source,output=output,category='standalone',dataset='cp5',direction='pp',trigger=None,selection='jetId')
            produce_source(entry)
            key=standalone_projection_key('hGenDijetPtEtaCM','standard',(90,120))
            self.assertNotIn('etaCut',key)
            cached=load_cached_eta_projection(output,'hGenDijetPtEtaCM','standard',(90,120))
            self.assertEqual(cached.Integral(),7.)
            self.assertEqual(cached.GetBinError(cached.FindBin(.5)),5.)
            self.assertEqual(ensure_npdf_projections(entry)['added'],[])
            file=ROOT.TFile.Open(str(output),'UPDATE');file.Delete(key+';*');file.Close()
            extension=ensure_npdf_projections(entry)
            self.assertEqual(extension['added'],[key])
            self.assertNotEqual(extension['output'],output)
            restored=load_cached_eta_projection(extension['output'],'hGenDijetPtEtaCM','standard',(90,120))
            self.assertEqual(restored.Integral(),7.)
            self.assertEqual(restored.GetBinError(restored.FindBin(.5)),5.)
