"""Generic sparse projections preserve input state and existing numerics."""
import unittest
from array import array
from unittest.mock import patch

from hist_analysis.python.notebook_setup import load_root
ROOT = load_root(batch=True)
from hist_analysis.python.projections import project_histogram, _range_to_bins


class SparseProjectionTest(unittest.TestCase):
    def setUp(self):
        self.add_directory = ROOT.TH1.AddDirectoryStatus()
        ROOT.TH1.AddDirectory(False)

    def tearDown(self):
        ROOT.TH1.AddDirectory(self.add_directory)

    def make_sparse(self):
        obj = ROOT.THnSparseD('range_fixture', ';x;y;z', 3,
                             array('i', [4]*3), array('d', [0]*3), array('d', [4]*3))
        obj.Sumw2()
        for coords, weight in [([.5, .5, .5], 2.5), ([2.5, 2.5, 2.5], -1.25),
                               ([-.5, .5, .5], 3.), ([4.5, .5, .5], 4.)]:
            obj.Fill(array('d', coords), weight)
        return obj

    def state(self, obj):
        return [(obj.GetAxis(i).GetFirst(), obj.GetAxis(i).GetLast(),
                 bool(obj.GetAxis(i).TestBit(ROOT.TAxis.kAxisRange)))
                for i in range(obj.GetNdimensions())]

    def assert_hist_same(self, a, b):
        self.assertEqual(a.GetNcells(), b.GetNcells())
        self.assertEqual(a.GetEntries(), b.GetEntries())
        for i in range(a.GetNcells()):
            self.assertEqual(a.GetBinContent(i), b.GetBinContent(i))
            self.assertEqual(a.GetBinError(i), b.GetBinError(i))
        for i in range(1, a.GetNbinsX()+2):
            self.assertEqual(a.GetXaxis().GetBinLowEdge(i), b.GetXaxis().GetBinLowEdge(i))
        self.assertFalse(b.GetDirectory())

    def test_existing_projection_results_and_all_range_states(self):
        for bounds in [None, (1, 4), (2, 3), (-1, 0), (0, 5), (5, 5)]:
            for spec in [{'axis_index': 1},
                         {'axis_index': 1, 'axis_ranges': {0: (0., 1.)}},
                         {'axis_index': 0, 'axis_ranges': {0: (1., 3.), 2: None}}]:
                with self.subTest(bounds=bounds, spec=spec):
                    source = self.make_sparse()
                    if bounds: source.GetAxis(0).SetRange(*bounds)
                    source.GetAxis(2).SetRange(1, 3)
                    before = self.state(source)
                    legacy = source.Clone('legacy')
                    for idx, interval in spec.get('axis_ranges', {}).items():
                        legacy.GetAxis(idx).SetRange(*_range_to_bins(legacy.GetAxis(idx), interval))
                    expected = legacy.Projection(spec['axis_index'])
                    expected.SetDirectory(0)
                    actual = project_histogram(source, spec)
                    self.assert_hist_same(expected, actual)
                    self.assertEqual(before, self.state(source))

    def test_repeated_calls_do_not_inherit_previous_selection(self):
        source = self.make_sparse()
        full = project_histogram(source, {'axis_index': 1})
        restricted = project_histogram(source, {'axis_index': 1, 'axis_ranges': {0: (0., 1.)}})
        self.assertEqual(restricted.Integral(), 2.5)
        self.assert_hist_same(full, project_histogram(source, {'axis_index': 1}))

    def test_ranges_restored_after_selection_and_projection_failures(self):
        for clone_failure in [False, True]:
            source = self.make_sparse()
            source.GetAxis(0).SetRange(-1, 0)
            before = self.state(source)
            spec = {'axis_index': 1, 'axis_ranges': {0: (1., 2.)}}
            if clone_failure:
                with patch('hist_analysis.python.projections._clone_detached', side_effect=RuntimeError('test')):
                    with self.assertRaises(RuntimeError): project_histogram(source, spec)
            else:
                spec['axis_ranges'][2] = (1.,)  # fail after the first axis was changed
                with self.assertRaises(ValueError): project_histogram(source, spec)
            self.assertEqual(before, self.state(source))


if __name__ == '__main__':
    unittest.main()
