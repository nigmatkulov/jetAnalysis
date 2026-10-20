"""Stable output identity, archive isolation, and total-plot presentation."""
import tempfile
import unittest
from pathlib import Path
import hashlib

import pandas as pd
from hist_analysis.python.systematic_total import (
    total_systematic_configuration_tag, total_systematic_output_stem,
    archive_systematic_inputs, plot_systematics,
)


class TotalOutputTest(unittest.TestCase):
    def setUp(self):
        self.manifest = pd.DataFrame([dict(source=s, csv=f'/inputs/{s}/source.csv',
            sha256='old', pt_low=50, pt_high=60, observable='full', frame='CM',
            direction='combined', selection='jetId', eta_cut=1.9,
            binning='normal', production_sample=sample, extraction='fit',
            full_fit='pol2', fb_fit='pol1', method='folded_deviation',
            combination='maximum', smoothing='smoothed', weights='Std',
            reference='refWithRecoPtAve' if s=='pointing_resolution' else None)
            for s,sample in [('JER','embedding'),('JEU','data'),
                             ('pileup','data'),('pointing_resolution','embedding')]])

    def test_identity_ignores_rerun_bytes_order_and_bin_subset(self):
        original = total_systematic_output_stem(self.manifest,(50,60))
        changed = self.manifest.iloc[::-1].copy()
        changed['sha256']='refreshed bytes'
        changed['pt_low']=80
        changed['pt_high']=90
        changed['copied_csv']='new archive path'
        self.assertEqual(original,total_systematic_output_stem(changed,(50,60)))
        self.assertNotEqual(original,total_systematic_output_stem(changed,(80,90)))
        self.assertNotIn('__',original)
        self.assertLess(len(original+'_provenance.json'),256)
        self.assertIn('_normal_full_',original)
        self.assertIn('_fullFit_pol2_fbFit_pol1_',original)
        self.assertIn('_method_folded_deviation_',original)

    def test_meaningful_variants_do_not_overwrite(self):
        baseline=total_systematic_output_stem(self.manifest,(50,60))
        for field,value in [('observable','fb'),('direction','Pbgoing'),
             ('selection','trkMax'),('binning','test'),('eta_cut',1.5),
             ('method','direct_ratio'),('extraction','bin_by_bin'),
             ('full_fit','pol4'),('fb_fit','pol2'),('combination','average'),
             ('smoothing','nonsmoothed'),('weights','1'),
             ('reference','refWithRefPtAve'),('csv','/another-production/source.csv')]:
            with self.subTest(field=field):
                m=self.manifest.copy();m.loc[0,field]=value
                self.assertNotEqual(baseline,total_systematic_output_stem(m,(50,60)))
        self.assertNotEqual(baseline,total_systematic_output_stem(self.manifest.iloc[:2],(50,60)))

    def test_archives_with_same_input_filename_remain_separate(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);a=root/'a'/'input.csv';b=root/'b'/'input.csv'
            a.parent.mkdir();b.parent.mkdir();a.write_text('first production');b.write_text('second production')
            outputs=root/'outputs'
            archived=[]
            for path in (a,b):
                m=self.manifest.iloc[:1].copy()
                m['csv']=str(path);m['sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
                prefix=total_systematic_configuration_tag(m)+'_'
                first=archive_systematic_inputs(m,outputs,filename_prefix=prefix)
                second=archive_systematic_inputs(m,outputs,filename_prefix=prefix)
                self.assertEqual(first.copied_csv.tolist(),second.copied_csv.tolist())
                archived.append(first.iloc[0].copied_csv)
            self.assertNotEqual(*archived)
            self.assertEqual((outputs/archived[0]).read_bytes(),a.read_bytes())
            self.assertEqual((outputs/archived[1]).read_bytes(),b.read_bytes())

    def test_plot_has_shared_headers_complete_legend_and_unchanged_arrays(self):
        table=pd.DataFrame({'x_low':[-2.,0.], 'x_high':[0.,2.],
            'JER':[.01,.02], 'JEU':[.02,.03], 'pileup':[.005,.01],
            'pointing_resolution':[.015,.025], 'total':[.03,.045]})
        before=table.copy(deep=True)
        fig=plot_systematics(table,self.manifest,(50,60),y_range=(0,5))
        try:
            fig.canvas.draw()  # Includes mathtext rendering.
            pd.testing.assert_frame_equal(table,before)
            ax=fig.axes[0]
            labels=set(t.get_text() for t in ax.get_legend().get_texts())
            self.assertEqual(labels,{'JER','JES','Pileup','Pointing resolution','Total'})
            texts=' '.join(t.get_text() for t in ax.texts)
            self.assertNotIn('Preliminary',texts)
            self.assertNotIn('Simulation',texts)
            self.assertIn('Data',texts)
            self.assertIn('Embedding MC',texts)
            self.assertIn('Combined',texts)
            self.assertEqual(ax.get_ylim(),(0,5))
            self.assertEqual(tuple(fig.get_size_inches()),(8,8))
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)
        m=self.manifest[self.manifest.production_sample=='embedding']
        fig=plot_systematics(table,m,(50,60))
        try:
            self.assertIn('Simulation',' '.join(t.get_text() for t in fig.axes[0].texts))
        finally:
            plt.close(fig)


if __name__=='__main__':
    unittest.main()
