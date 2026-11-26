"""Native ROOT round-trip and notebook wiring checks (no production inputs)."""
from array import array
from dataclasses import dataclass
import json
from pathlib import Path
import tempfile
import unittest

from hist_analysis.python.notebook_setup import load_root
ROOT = load_root()
from hist_analysis.python.notebook_output import (
    NotebookRootOutput, capture_canvas, systematic_table_histograms,
)


def keys(directory, prefix=''):
    result = {}
    for key in directory.GetListOfKeys():
        name = prefix + key.GetName()
        result[name] = key.GetCycle()
        if key.GetClassName() == 'TDirectoryFile':
            result.update(keys(directory.GetDirectory(key.GetName()), name + '/'))
    return result


class NotebookOutputTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.archive = NotebookRootOutput(self.directory, 'analysis')

    def histogram(self):
        hist = ROOT.TH1D('original', 'weighted;x;yield', 3, array('d', [-2, 0, .5, 4]))
        hist.SetDirectory(0)
        hist.Sumw2()
        for x, weight in [(-3, 2), (-1, 3), (.2, -1), (.2, 4), (6, 5)]:
            hist.Fill(x, weight)
        return hist

    def assert_histogram(self, expected, actual):
        self.assertTrue(actual)
        self.assertEqual(expected.ClassName(), actual.ClassName())
        self.assertEqual(expected.GetDimension(), actual.GetDimension())
        self.assertEqual(expected.GetNcells(), actual.GetNcells())
        self.assertEqual(expected.GetSumw2N(), actual.GetSumw2N())
        for i in range(expected.GetNcells()):
            self.assertEqual(expected.GetBinContent(i), actual.GetBinContent(i))
            self.assertEqual(expected.GetBinError(i), actual.GetBinError(i))
        for getter in ['GetXaxis', 'GetYaxis', 'GetZaxis'][:expected.GetDimension()]:
            left, right = getattr(expected, getter)(), getattr(actual, getter)()
            self.assertEqual(left.GetNbins(), right.GetNbins())
            for i in range(1, left.GetNbins()+2):
                self.assertEqual(left.GetBinLowEdge(i), right.GetBinLowEdge(i))
        self.assertEqual(expected.Integral(), actual.Integral())
        self.assertEqual(expected.GetEntries(), actual.GetEntries())

    def test_native_roundtrip_and_lifetime(self):
        hist = self.histogram()
        empty = ROOT.TH2D('empty', '', 2, 0, 1, 3, -1, 2)
        empty.SetDirectory(0)
        h3 = ROOT.TH3D('three', '', 2, 0, 1, 3, 0, 3, 2, 0, 2)
        h3.SetDirectory(0); h3.Sumw2(); h3.Fill(.2, 2, .5, 3)
        @dataclass
        class Result:
            histogram: object
        directory_before = ROOT.gDirectory.GetPath()
        self.archive.collect({'weighted': hist, 'empty': empty, 'result': Result(h3)})
        self.archive.finish({'SELECTION': 'jetId'})
        self.assertEqual(directory_before, ROOT.gDirectory.GetPath())
        f = ROOT.TFile.Open(str(self.archive.path))
        try:
            for expected, name in [(hist,'weighted'),(empty,'empty'),(h3,'result/histogram')]:
                self.assert_histogram(expected, f.Get('results/'+name))
            self.assertEqual(f.Get('results/weighted').GetName(), 'results__weighted')
            self.assertEqual(hist.GetName(), 'original')
            self.assertFalse(hist.GetDirectory())
            metadata = json.loads(f.Get('notebookConfiguration').GetTitle())
            self.assertEqual(metadata['configuration']['SELECTION'], 'jetId')
        finally:
            f.Close()
        hist.Fill(1, 2)  # original still alive after output closure

    def test_plot_snapshot_and_key_overwrite(self):
        hist = self.histogram()
        canvas = ROOT.TCanvas('random_name', '', 400, 400)
        canvas.Divide(1, 2); canvas.cd(1); hist.Draw('hist')
        capture_canvas(canvas, self.directory/'figure.pdf')
        hist.SetBinContent(1, 999)
        f = ROOT.TFile.Open(str(self.archive.path))
        try:
            histogram_keys = [p for p in keys(f) if f.Get(p).InheritsFrom('TH1')]
            self.assertEqual(len(histogram_keys), 1)
            self.assertEqual(f.Get(histogram_keys[0]).GetBinContent(1), 3)
        finally: f.Close()
        capture_canvas(canvas, self.directory/'figure.pdf')
        self.archive.finish({'histogram':hist})
        self.archive.finish({'histogram':hist})
        f = ROOT.TFile.Open(str(self.archive.path))
        try:
            self.assertTrue(all(cycle == 1 for cycle in keys(f).values()))
            self.assertEqual(f.Get(histogram_keys[0]).GetBinContent(1), 999)
        finally: f.Close()
        # A new execution removes stale objects, not just their newest cycles.
        self.archive = NotebookRootOutput(self.directory, 'analysis')
        self.archive.collect({'new':hist})
        f = ROOT.TFile.Open(str(self.archive.path))
        try:
            self.assertFalse(f.Get('results/histogram'))
            self.assertTrue(f.Get('results/new'))
        finally: f.Close()
        canvas.Close()

    def test_explicit_exports_and_independent_ratio(self):
        forward, backward = self.histogram(), self.histogram()
        backward.Scale(2)
        ratio = forward.Clone('fb'); ratio.SetDirectory(0)
        ratio.Divide(forward, backward, 1., 1., '')
        self.archive.write_unfolding_output(self.directory/'legacy.root', histograms=[ratio],
                                            metadata={'error': 'independent'})
        self.assertFalse((self.directory/'legacy.root').exists())
        f = ROOT.TFile.Open(str(self.archive.path))
        try:
            saved = f.Get('exports/legacy/fb')
            self.assert_histogram(ratio, saved)
            expected_error = ((forward.GetBinError(1)/backward.GetBinContent(1))**2 +
                              (forward.GetBinContent(1)*backward.GetBinError(1)/
                               backward.GetBinContent(1)**2)**2)**.5
            self.assertAlmostEqual(saved.GetBinError(1),expected_error)
        finally: f.Close()

    def test_systematic_table(self):
        import pandas as pd
        table = pd.DataFrame(dict(x_low=[0,.5],x_high=[.5,2],source=[.02,.03],total=[.04,.05]))
        histograms = systematic_table_histograms(table)
        self.assertEqual(histograms['total_percent'].GetBinContent(2),5)
        self.assertEqual(histograms['total_fraction'].GetBinError(1),0)
        self.archive.collect(histograms)
        f = ROOT.TFile.Open(str(self.archive.path))
        try:self.assert_histogram(histograms['source_percent'], f.Get('results/source_percent'))
        finally:f.Close()

    def test_sparse_graph_profile_and_producer_import(self):
        sparse = ROOT.THnSparseD('sparse', '', 2, array('i',[2,3]),
                                  array('d',[0,0]), array('d',[2,3]))
        sparse.Sumw2(); sparse.Fill(array('d',[.2,1.2]),2.5)
        graph = ROOT.TGraphErrors(2)
        graph.SetPoint(0,.25,3); graph.SetPointError(0,.25,.7)
        graph.SetPoint(1,1,5); graph.SetPointError(1,.5,.8)
        profile = ROOT.TProfile('profile','',2,0,2)
        profile.SetDirectory(0); profile.Fill(.2,3); profile.Fill(.2,5)
        self.archive.collect({'sparse':sparse,'graph':graph,'profile':profile,
                              'imported_class':ROOT.TH1D})
        source = self.directory/'producer.root'
        f = ROOT.TFile.Open(str(source),'RECREATE')
        try:
            f.mkdir('nested').cd(); profile.Write()
        finally:f.Close()
        self.archive.import_file(source,prefix='sample')
        f = ROOT.TFile.Open(str(self.archive.path))
        try:
            saved = f.Get('results/sparse')
            coordinates = array('i',[1,2])
            self.assertEqual(saved.GetBinContent(coordinates),2.5)
            self.assertEqual(saved.GetBinError(coordinates),2.5)
            self.assertEqual(f.Get('results/graph').GetErrorY(0),.7)
            self.assert_histogram(profile,f.Get('results/profile'))
            self.assert_histogram(profile,f.Get('production/sample/nested/profile'))
            self.assertFalse(f.Get('results/imported_class'))
        finally:f.Close()

    def test_all_notebooks_wired(self):
        import nbformat
        from IPython.core.interactiveshell import InteractiveShell
        transform = InteractiveShell.instance().input_transformer_manager.transform_cell
        paths = sorted(Path('hist_analysis/notebooks').glob('*.ipynb'))
        standalone = {
            '01_inventory',
            '04_systematics_old_vs_new',
            '05_unfold2D_validation',
            '06_data_jet_corrections_comparison',
        }
        self.assertTrue(paths)
        for path in paths:
            with self.subTest(notebook=path.name):
                notebook = nbformat.read(path, as_version=4)
                nbformat.validate(notebook)
                code = '\n'.join(c.source for c in notebook.cells if c.cell_type == 'code')
                for index, cell in enumerate(notebook.cells):
                    if cell.cell_type == 'code':
                        compile(transform(cell.source), f'{path}:{index}', 'exec')
                # The local study notebooks write their own per-run artifacts.
                if path.stem in standalone or path.stem.startswith(('05_unfold2D_eta_', '05_unfold2D_study_')):
                    self.assertNotIn('NotebookRootOutput', code)
                else:
                    self.assertEqual(code.count('ROOT_OUTPUT = NotebookRootOutput('), 1)
                    self.assertIn(repr(path.stem), code)
                    self.assertEqual(code.count('ROOT_OUTPUT.finish('), 1)


if __name__ == '__main__': unittest.main()
