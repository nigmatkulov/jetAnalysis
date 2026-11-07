"""Exact selection names, read compatibility, and preservation of old artifacts."""
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from hist_analysis.python.output_naming import (
    OutputDirectory, OutputSelection, output_filename, number_token,
    output_record, find_outputs, sibling_output, related_output,
)


class OutputNamingTests(unittest.TestCase):
    def test_dated_renderer_layout_and_single_underscores(self):
        from hist_analysis.python.output_naming import output_date, mplhep_output_path
        with tempfile.TemporaryDirectory() as d:
            for method in ('noSel', 'trkMax', 'jetId'):
                directory = OutputDirectory(d, OutputSelection(method, 1.9, [(80,90)]))
                native = directory / 'plot___diagnostic.pdf'
                rendered = directory / 'mplhep' / 'plot___diagnostic.pdf'
                self.assertEqual(native.parent, Path(d)/output_date()/method/'root')
                self.assertEqual(rendered.parent, Path(d)/output_date()/method/'mplhep')
                self.assertEqual(mplhep_output_path(native), rendered)
                self.assertNotIn('__', native.name)
                self.assertNotIn('__', rendered.name)

    def test_custom_workflow_layout_and_legacy_reader_path(self):
        with tempfile.TemporaryDirectory() as d:
            base = Path(d)/'custom'
            original = base/'jetId'/'normal'
            context = OutputSelection('jetId',1.9,[(80,90)])
            output = OutputDirectory(original, context, workflow_root=base)
            path = output/'plot.pdf'
            path.write_bytes(b'fixture')
            self.assertEqual(path.parent.parent.parent.parent, base)
            self.assertEqual(find_outputs(original,'plot.pdf',selection='jetId',eta_cut=1.9,interval=(80,90)),[path])

    def test_all_selection_dimensions_are_distinct(self):
        names=set()
        for selection in ('noSel','trkMax','jetId'):
            for cut in (1.5,1.9,2.0,1.9000000000000001):
                for interval in ((80,90),(90,120),(120,160),(200,300),(80.00000000000001,90)):
                    name=output_filename('dijet_eta','data_vs_mc',selection=selection,
                                         eta_cut=cut,pt_interval=interval)
                    self.assertIn(selection,name)
                    self.assertIn('etaCM'+number_token(cut),name)
                    self.assertIn('ptAve'+number_token(interval[0])+'to'+number_token(interval[1]),name)
                    self.assertNotIn(name,names);names.add(name)
        self.assertEqual(len(names),60)
        self.assertNotEqual(number_token(1.9000000000000001),number_token(1.9))
        self.assertEqual(number_token(Decimal('1.9000000000000000000001')),'1p9000000000000000000001')

    def test_complete_binset_identity_and_manifest(self):
        a=OutputSelection('jetId',1.9,[(80,90),(90,120)])
        b=OutputSelection('jetId',1.9,[(80,100),(100,120)])
        self.assertNotEqual(a.token(),b.token())
        with tempfile.TemporaryDirectory() as d:
            output=OutputDirectory(d,a)
            path=output/'multiple_intervals.root';path.write_bytes(b'ROOT fixture')
            self.assertIn('ptAveBins',path.name)
            self.assertEqual(output_record(path)['identity']['pt_bins'],[['80','90'],['90','120']])

    def test_existing_outputs_and_runs_are_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            old=Path(d)/'plot.pdf';old.write_bytes(b'validated')
            context=OutputSelection('trkMax',2.0,[(90,120)])
            a=OutputDirectory(d,context)/'plot.pdf';a.write_bytes(b'new run')
            b=OutputDirectory(d,context)/'plot.pdf';b.write_bytes(b'another run')
            self.assertNotEqual(a,b)
            self.assertEqual(old.read_bytes(),b'validated')
            self.assertEqual(a.read_bytes(),b'new run')
            found=find_outputs(d,'plot.pdf',selection='trkMax',eta_cut=2.,interval=(90,120))
            self.assertEqual(found,[b])
            self.assertEqual(find_outputs(d,'plot.pdf',selection='noSel'),[old])

    def test_current_interval_and_sibling_lookup(self):
        with tempfile.TemporaryDirectory() as d:
            output=OutputDirectory(d,OutputSelection('noSel',1.9,[(80,90),(90,120)]))
            stem='total_ptave_90_120'
            paths={suffix:output/(stem+suffix) for suffix in ('_provenance.json','_total.csv','.csv')}
            for p in paths.values():p.write_text('fixture')
            for path in paths.values():
                self.assertIn('ptAve90to120',path.name)
                self.assertNotIn('ptAveBins',path.name)
            self.assertEqual(sibling_output(paths['_provenance.json'],'_provenance.json','_total.csv'),paths['_total.csv'])
            diagnostic=related_output(paths['_total.csv'],'_diagnostic')
            self.assertIn('noSel_etaCM1p9_ptAve90to120',diagnostic.name)

    def test_rounded_legacy_stem_is_never_ambiguous(self):
        context=OutputSelection('jetId',1.9,[(80.00000000000001,90),(80.00000000000003,90)])
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError,'Ambiguous rounded interval'):
                OutputDirectory(d,context)/'full_ptave_80_90.pdf'

    def test_binset_digest_collision_is_detected(self):
        with tempfile.TemporaryDirectory() as d, patch('hist_analysis.python.output_naming._digest',return_value='a'*64):
            a=OutputDirectory(d,OutputSelection('jetId',1.9,[(80,90),(90,120)]),new_run=False,run_id='test')
            b=OutputDirectory(d,OutputSelection('jetId',1.9,[(80,100),(100,120)]),new_run=False,run_id='test')
            p=a/'summary.root';p.write_bytes(b'keep')
            with self.assertRaisesRegex(ValueError,'identity collision'):b/'summary.root'
            self.assertEqual(p.read_bytes(),b'keep')

    def test_long_systematic_names_and_frame_acceptance(self):
        with tempfile.TemporaryDirectory() as d:
            directory=OutputDirectory(d,OutputSelection('jetId',{'CM':1.9,'Lab':2.4},[(200,300)]))
            p=directory/('pointing_'+'long_method_'*20+'etaLab_24_ptave_200_300.csv')
            self.assertLessEqual(len(p.name.encode()),220)
            self.assertIn('etaLab2p4',p.name)
            self.assertNotIn('etaCM1p9',p.name)
            self.assertIn('ptAve200to300',p.name)

    def test_invalid_boundaries_fail(self):
        for value in (float('nan'),float('inf'),True):
            with self.assertRaises(ValueError):number_token(value)
        with self.assertRaises(ValueError):OutputSelection('guess',1.9,[(80,90)]).identity()
        with self.assertRaises(ValueError):OutputSelection('jetId',1.9,[(90,80)]).identity()

    def test_root_and_plot_writers_for_each_selection(self):
        from hist_analysis.python.notebook_setup import load_root
        ROOT=load_root()
        from hist_analysis.python.histogram_io import save_histograms
        from hist_analysis.python.root_style import save_canvas, configure_mplhep_export
        configure_mplhep_export(False)
        with tempfile.TemporaryDirectory() as d:
            for method in ('noSel','trkMax','jetId'):
                for cut in (1.5,1.9):
                    for interval in ((80,90),(120,160)):
                        directory=OutputDirectory(d,OutputSelection(method,cut,[interval]))
                        root=directory/'fixture.root'
                        h=ROOT.TH1D('hNamingFixture','',2,-1,1);h.SetDirectory(0);h.SetBinContent(1,3.)
                        save_histograms(str(root),[h])
                        file=ROOT.TFile.Open(str(root));self.assertEqual(file.Get('hNamingFixture').GetBinContent(1),3.);file.Close()
                        self.assertEqual(output_record(root)['identity']['selection'],[method])
                        canvas=ROOT.TCanvas('naming_canvas','',250,250);h.Draw('E1')
                        pdf=directory/'fixture.pdf';save_canvas(canvas,pdf,save_png=True)
                        self.assertTrue(pdf.is_file());self.assertTrue(pdf.with_suffix('.png').is_file())
                        canvas.Close()

    def test_cache_without_interval_provenance_is_not_reused(self):
        from hist_analysis.python.notebook_setup import load_root
        ROOT=load_root()
        from hist_analysis.python.histogram_production import load_cached_eta_projection
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'unverified.root'
            file=ROOT.TFile.Open(str(path),'RECREATE')
            ROOT.TH1D('hGenDijetEtaCM_standard_ptAve80to90','',2,-1,1).Write();file.Close()
            with self.assertRaisesRegex(KeyError,'no declared .*provenance'):
                load_cached_eta_projection(path,'hGenDijetPtEtaCM','standard',(80,90))

    def test_selection_roles_and_multi_eta_outputs(self):
        data=OutputSelection({'data':'jetId','response':'trkMax'},1.9,[(80,90)])
        swapped=OutputSelection({'data':'trkMax','response':'jetId'},1.9,[(80,90)])
        self.assertNotEqual(data.token(),swapped.token())
        with tempfile.TemporaryDirectory() as d:
            directory=OutputDirectory(d,OutputSelection('jetId',(1.5,1.9),[(80,90)]))
            first=directory/'full_etaCM_15_ptave_80_90.pdf'
            second=directory/'full_etaCM_19_ptave_80_90.pdf'
            self.assertIn('etaCM1p5_',first.name);self.assertNotIn('etaCM1p9',first.name)
            self.assertIn('etaCM1p9_',second.name);self.assertNotIn('etaCM1p5',second.name)

    def test_reference_acceptance_is_verified_instead_of_guessed(self):
        from hist_analysis.python.notebook_setup import load_root
        ROOT=load_root()
        from hist_analysis.python.output_naming import verified_reference_eta_cuts
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'misleading_etaCM2p4.root'
            file=ROOT.TFile.Open(str(path),'RECREATE')
            h=ROOT.TH2D('hGenDijetPtEtaCM','',2,0,200,2,-2,2);h.Sumw2();h.Fill(80,.5,3.);h.Write()
            a=h.Clone('hGenDijetPtEtaCM_0');a.Write()
            b=h.Clone('hGenDijetPtEtaCM_1');b.Fill(80,.5,2.);b.Write()
            ROOT.TNamed('histogramProductionMetadata',json.dumps({'eta_cut_index_to_value':{'0':1.9,'1':2.4}})).Write()
            file.Close()
            self.assertEqual(verified_reference_eta_cuts([('CP5GenCM',path,'hGenDijetPtEtaCM')]),{'CP5GenCM':1.9})

    def test_reference_without_acceptance_evidence_is_rejected(self):
        from hist_analysis.python.notebook_setup import load_root
        ROOT=load_root()
        from hist_analysis.python.output_naming import verified_reference_eta_cuts
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'looks_like_etaCM1p9.root'
            file=ROOT.TFile.Open(str(path),'RECREATE')
            h=ROOT.TH2D('hGenDijetPtEtaCM','',2,0,200,2,-2,2);h.Write();file.Close()
            with self.assertRaisesRegex(ValueError,'acceptance is undeclared'):
                verified_reference_eta_cuts([('CP5GenCM',path,'hGenDijetPtEtaCM')])

    def test_output_paths_cannot_escape_the_new_run(self):
        with tempfile.TemporaryDirectory() as d:
            directory=OutputDirectory(d,OutputSelection('jetId',1.9,[(80,90)]))
            with self.assertRaisesRegex(ValueError,'beneath'):directory/'../old.pdf'
            with self.assertRaisesRegex(ValueError,'beneath'):directory.file('../old.pdf')
