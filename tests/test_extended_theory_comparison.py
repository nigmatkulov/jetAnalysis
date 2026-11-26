"""Direct source-bin checks, independent of the shared projection/alignment."""
import math
import unittest
import ROOT
from hist_analysis.config.files import BASE_DIR
from hist_analysis.config.histograms import NPDF_TRIGGER_BINS
from hist_analysis.python.histogram_io import resolve_data_file, load_histogram
from hist_analysis.python.histogram_production import resolve_eta_source_key
from hist_analysis.python.extended_theory_comparison import prepare_comparison, MODELS

THEORY = 'hist_analysis/root_files/extended-comparison.root'


class ExtendedComparisonTest(unittest.TestCase):
    def test_direct_source_sums_and_ratios(self):
        for trigger, intervals in NPDF_TRIGGER_BINS.items():
            path = resolve_data_file(BASE_DIR/'exp',trigger,'combined','jetId')
            source = load_histogram(str(path),resolve_eta_source_key(path,'hRecoDijetPtEtaCM_4'))
            for interval in intervals:
                result = prepare_comparison(path, THEORY, interval)
                data = result['raw']['Data']
                for i in range(1,data.GetNbinsX()+1):
                    low, high = data.GetXaxis().GetBinLowEdge(i),data.GetXaxis().GetBinUpEdge(i)
                    values=[]; variances=[]
                    for x in range(1,source.GetNbinsX()+1):
                        if not interval[0] <= source.GetXaxis().GetBinCenter(x) < interval[1]:
                            continue
                        for y in range(1,source.GetNbinsY()+1):
                            if low <= source.GetYaxis().GetBinCenter(y) < high:
                                values.append(source.GetBinContent(x,y));variances.append(source.GetBinError(x,y)**2)
                    self.assertAlmostEqual(data.GetBinContent(i),math.fsum(values),delta=max(1e-8,abs(math.fsum(values))*1e-12))
                    self.assertAlmostEqual(data.GetBinError(i)**2,math.fsum(variances),delta=max(1e-8,math.fsum(variances)*1e-12))
                for label,prefix in MODELS.items():
                    theory=load_histogram(THEORY,f'pt{interval[0]:03d}_{interval[1]:03d}_eta1p9/{prefix}_combined_propagated')
                    raw=result['raw'][label]
                    for i in range(1,raw.GetNbinsX()+1):
                        j=theory.FindBin(raw.GetBinCenter(i))
                        self.assertAlmostEqual(raw.GetBinContent(i),theory.GetBinContent(j))
                        self.assertAlmostEqual(raw.GetBinError(i),theory.GetBinError(j))
                    Dtot=data.Integral();Ttot=raw.Integral()
                    for i in result['ratios'][label]._valid_ratio_bins:
                        d=data.GetBinContent(i)/Dtot;t=raw.GetBinContent(i)/Ttot
                        err=math.hypot(data.GetBinError(i)/Dtot/t,d*raw.GetBinError(i)/Ttot/t**2)
                        self.assertAlmostEqual(result['ratios'][label].GetBinContent(i),d/t)
                        self.assertAlmostEqual(result['ratios'][label].GetBinError(i),err)

    def test_rejects_invalid_configuration(self):
        with self.assertRaises(ValueError):
            prepare_comparison('',THEORY,(50,60),campaign='unknown')

if __name__=='__main__':
    unittest.main()
