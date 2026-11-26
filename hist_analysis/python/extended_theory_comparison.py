"""Normalized etaCM data shapes versus supplied NLOJet central spectra."""
from array import array
import math

from hist_analysis.python.histogram_io import load_histogram
from hist_analysis.python.histogram_ops import align_common_binning, normalize_histogram, ratio_to_nominal
from hist_analysis.python.data_distributions import project_data_histogram

MODELS = {'EPPS21 x CT18': 'npdf_central', 'CT18 (pp reference)': 'pdf_central'}


def prepare_comparison(data_file, theory_file, interval, *, eta_cut=1.9,
                       campaign='combined', normalization='integral'):
    """Keep whole accepted edge bins; exclude flows and invalid exterior bins.

    Uniform theory eta widths cancel in shape normalization. Stored integration
    errors and data Sumw2 are scaled with fixed normalization constants; no
    normalization covariance, PDF, scale or experimental systematics are added.
    """
    import ROOT
    from hist_analysis.config.histograms import eta_cut_index
    if campaign not in {'combined', 'nov09', 'nov18'}:
        raise ValueError('Unknown campaign')
    if eta_cut not in (1.5, 1.9, 2.4):
        raise ValueError('Acceptance not supplied')
    if normalization not in {'integral', 'bin_width'}:
        raise ValueError('Shape normalization must be integral or bin_width')
    key = f'hRecoDijetPtEtaCM_{eta_cut_index(eta_cut)}'
    from hist_analysis.python.histogram_production import resolve_eta_source_key
    source = load_histogram(str(data_file), resolve_eta_source_key(data_file, key))
    axis = source.GetXaxis()
    pt_edges = [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins()+2)]
    if not interval[0] < interval[1] or not all(
            any(math.isclose(edge, bound, abs_tol=1e-8, rel_tol=0) for edge in pt_edges)
            for bound in interval):
        raise ValueError(f'Data bins do not exactly cover {interval}')
    data = project_data_histogram(data_file, key, 'eta', interval, name='data_eta')
    directory = f'pt{interval[0]:03d}_{interval[1]:03d}_eta{eta_cut:.1f}'.replace('.', 'p')
    originals = {'Data': data}
    for label, prefix in MODELS.items():
        h = load_histogram(str(theory_file), f'{directory}/{prefix}_{campaign}_propagated')
        if h.GetDimension() != 1:
            raise ValueError('Expected a one-dimensional central spectrum')
        originals[label] = h
    # Accepted theory bins may straddle the jet cut; preserve their full edges.
    theory = originals[next(iter(MODELS))]
    selected = [i for i in range(1, theory.GetNbinsX()+1)
                if abs(theory.GetBinCenter(i)) <= eta_cut + 1e-8]
    first, last = min(selected), max(selected)
    edges = [theory.GetXaxis().GetBinLowEdge(i) for i in range(first, last+2)]
    if not all(math.isclose(edges[i+1]-edges[i], edges[1]-edges[0], abs_tol=1e-9)
               for i in range(len(edges)-1)):
        raise ValueError('Nonuniform theory widths require an explicit units convention')
    cleaned = {'Data': data}
    for label in MODELS:
        h = originals[label]
        if any(not math.isclose(h.GetXaxis().GetBinLowEdge(i), theory.GetXaxis().GetBinLowEdge(i), abs_tol=1e-9)
               for i in range(1, theory.GetNbinsX()+2)):
            raise ValueError('Theory axes differ')
        c = ROOT.TH1D(f'{h.GetName()}_accepted', '', len(edges)-1, array('d', edges))
        c.SetDirectory(0); c.Sumw2()
        for j, i in enumerate(selected, 1):
            y, e = h.GetBinContent(i), h.GetBinError(i)
            if not math.isfinite(y) or not math.isfinite(e) or e < 0:
                raise ValueError(f'Invalid accepted theory bin: {label}/{interval}/{i}')
            c.SetBinContent(j, y); c.SetBinError(j, e)
        cleaned[label] = c
    raw = align_common_binning(cleaned)
    for label, h in raw.items():
        source = cleaned[label]
        # No finite nonzero yield or uncertainty may be discarded, including flows.
        outside = [i for i in range(source.GetNbinsX()+2)
                   if i in (0, source.GetNbinsX()+1) or
                   source.GetBinCenter(i) < edges[0] or source.GetBinCenter(i) > edges[-1]]
        if any(source.GetBinContent(i) != 0 or source.GetBinError(i) != 0 for i in outside):
            raise ValueError(f'Alignment discards yield/errors: {label}')
        if not math.isclose(h.Integral(), source.Integral(), rel_tol=1e-10, abs_tol=1e-8):
            raise ValueError(f'Integral changed: {label}')
        if any(not math.isfinite(h.GetBinContent(i)) or not math.isfinite(h.GetBinError(i))
               for i in range(1, h.GetNbinsX()+1)):
            raise ValueError(f'Invalid aligned data: {label}')
    shapes = {label: normalize_histogram(h, normalization) for label, h in raw.items()}
    ratios = {}
    for label in MODELS:
        r = ratio_to_nominal(shapes['Data'], shapes[label], name='Data_over_'+label)
        r._valid_ratio_bins = tuple(i for i in range(1,r.GetNbinsX()+1) if shapes[label].GetBinContent(i) != 0)
        ratios[label] = r
    return dict(raw=raw, shapes=shapes, ratios=ratios, x_range=(edges[0],edges[-1]))
