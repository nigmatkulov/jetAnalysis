"""Compare old/new producer outputs; report accepted full-only folding differences.

Run from the repository root with py-env/bin/python. Inputs are read-only.
"""
import argparse
import json
from pathlib import Path

from hist_analysis.python.notebook_setup import load_root
ROOT = load_root()
from hist_analysis.config.histograms import DIJET_PTAVE_BINS
from hist_analysis.python.forward_backward import fold_eta_histogram, forward_backward_ratio

ROOT.gInterpreter.Declare(r'''
bool sameFullEtaObject(TObject* a, TObject* b) {
    if (std::string(a->ClassName()) != b->ClassName()) return false;
    if (auto *ha = dynamic_cast<TH1*>(a)) {
        auto *hb = dynamic_cast<TH1*>(b);
        if (ha->GetNcells() != hb->GetNcells() || ha->GetEntries() != hb->GetEntries()) return false;
        for (int axis = 0; axis < ha->GetDimension(); ++axis) {
            auto *aa = axis == 0 ? ha->GetXaxis() : axis == 1 ? ha->GetYaxis() : ha->GetZaxis();
            auto *ab = axis == 0 ? hb->GetXaxis() : axis == 1 ? hb->GetYaxis() : hb->GetZaxis();
            if (aa->GetNbins() != ab->GetNbins()) return false;
            for (int i = 1; i <= aa->GetNbins()+1; ++i)
                if (aa->GetBinLowEdge(i) != ab->GetBinLowEdge(i)) return false;
        }
        for (int i = 0; i < ha->GetNcells(); ++i)
            if (ha->GetBinContent(i) != hb->GetBinContent(i) || ha->GetBinError(i) != hb->GetBinError(i)) return false;
    } else if (auto *sa = dynamic_cast<THnSparse*>(a)) {
        auto *sb = dynamic_cast<THnSparse*>(b);
        if (sa->GetNbins() != sb->GetNbins() || sa->GetEntries() != sb->GetEntries() || sa->GetNdimensions() != sb->GetNdimensions()) return false;
        std::vector<int> ca(sa->GetNdimensions()), cb(sa->GetNdimensions());
        for (int axis = 0; axis < sa->GetNdimensions(); ++axis) {
            auto *aa = sa->GetAxis(axis), *ab = sb->GetAxis(axis);
            if (aa->GetNbins() != ab->GetNbins()) return false;
            for (int i = 1; i <= aa->GetNbins()+1; ++i)
                if (aa->GetBinLowEdge(i) != ab->GetBinLowEdge(i)) return false;
        }
        for (Long64_t i = 0; i < sa->GetNbins(); ++i)
            if (sa->GetBinContent(i,ca.data()) != sb->GetBinContent(i,cb.data()) || ca != cb || sa->GetBinError2(i) != sb->GetBinError2(i)) return false;
    } else if (auto *na = dynamic_cast<TNamed*>(a)) {
        if (std::string(na->GetTitle()) != dynamic_cast<TNamed*>(b)->GetTitle()) return false;
    }
    return true;
}
''')


def compare(before, after, *, audit_existing=False):
    old = ROOT.TFile.Open(str(before))
    new = ROOT.TFile.Open(str(after))
    if not old or old.IsZombie() or not new or new.IsZombie():
        raise OSError('Cannot open comparison inputs')
    try:
        keys = {k.GetName() for k in old.GetListOfKeys()}
        revised = {k.GetName() for k in new.GetListOfKeys()}
        hemisphere = {k for k in keys if 'Dijet' in k and ('Forward' in k or 'Backward' in k)}
        if audit_existing:
            revised = keys-hemisphere
        if not audit_existing and revised != keys-hemisphere:
            raise AssertionError(f'Unexpected key changes: missing={keys-hemisphere-revised}, extra={revised-(keys-hemisphere)}')
        changed = [] if audit_existing else [k for k in sorted(revised) if not ROOT.sameFullEtaObject(old.Get(k), new.Get(k))]
        if changed: raise AssertionError(f'Retained objects changed: {changed}')
        report = dict(before=str(before), after=str(after), unchanged_objects=len(revised), removed_hemisphere_objects=len(hemisphere), sides={})
        if audit_existing:
            report.pop('unchanged_objects')
            report.pop('removed_hemisphere_objects')
            report['audit_only'] = True
        for label in ['Forward', 'Backward']:
            cells = errors = 0
            largest = 0.
            for key in sorted(k for k in hemisphere if label in k):
                full_key = key.replace('PtEta'+label, 'PtEtaCM').replace('EtaCM'+label+'Shifted', 'EtaCMShifted')
                source = new.Get(full_key)
                if not source: raise AssertionError(f'No full partner for {key}')
                derived = fold_eta_histogram(source, name='compare_'+key)[label == 'Backward']
                stored = old.Get(key)
                for i in range(stored.GetNcells()):
                    for attr, is_error in [('GetBinContent',False), ('GetBinError',True)]:
                        a,b = getattr(derived,attr)(i), getattr(stored,attr)(i)
                        if abs(a-b) > 1e-10*max(abs(a),abs(b),1e-20):
                            if is_error: errors += 1
                            else:
                                cells += 1
                                largest = max(largest, abs(a-b)/max(abs(a),abs(b),1e-20))
            report['sides'][label] = dict(content_difference_bins=cells,error_difference_bins=errors,max_relative_content_difference=largest)
        ratio_differences = ratio_error_differences = ratios_checked = 0
        largest_ratio_difference = 0.
        for key in sorted(k for k in hemisphere if 'Forward' in k):
            backward_key = key.replace('Forward', 'Backward')
            full_key = key.replace('PtEtaForward', 'PtEtaCM').replace('EtaCMForwardShifted', 'EtaCMShifted')
            stored_f, stored_b = old.Get(key), old.Get(backward_key)
            folded_f, folded_b = fold_eta_histogram(new.Get(full_key), name='ratio_'+key)
            intervals = DIJET_PTAVE_BINS if stored_f.GetDimension() == 2 else [None]
            for interval in intervals:
                inputs = []
                for index, h in enumerate([stored_f, stored_b, folded_f, folded_b]):
                    if interval is None:
                        projection = h.Clone(f'audit_{index}')
                    else:
                        x = h.GetXaxis()
                        low, high = x.FindBin(interval[0]+1e-6), x.FindBin(interval[1]-1e-6)
                        projection = h.ProjectionY(f'audit_{index}', low, high, 'e')
                    projection.SetDirectory(0)
                    projection.Rebin(2)
                    inputs.append(projection)
                expected = forward_backward_ratio(*inputs[:2], name='expected_ratio')
                actual = forward_backward_ratio(*inputs[2:], name='actual_ratio')
                ratios_checked += 1
                for i in range(actual.GetNcells()):
                    for attr, is_error in [('GetBinContent', False), ('GetBinError', True)]:
                        a,b = getattr(actual, attr)(i), getattr(expected, attr)(i)
                        if abs(a-b) > 1e-10*max(abs(a),abs(b),1e-20):
                            if is_error: ratio_error_differences += 1
                            else:
                                ratio_differences += 1
                                largest_ratio_difference = max(largest_ratio_difference, abs(a-b)/max(abs(a),abs(b),1e-20))
        report['fb_ratios_rebin2'] = dict(projections_checked=ratios_checked,
            content_difference_bins=ratio_differences, error_difference_bins=ratio_error_differences,
            max_relative_content_difference=largest_ratio_difference)
        return report
    finally:
        old.Close(); new.Close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('before', type=Path)
    parser.add_argument('after', type=Path, nargs='?')
    parser.add_argument('--audit-existing', action='store_true', help='Compare folded and stored hemispheres in one archived file')
    args = parser.parse_args()
    if not args.audit_existing and args.after is None: parser.error('after is required unless --audit-existing is used')
    print(json.dumps(compare(args.before, args.before if args.audit_existing else args.after, audit_existing=args.audit_existing), indent=2))
