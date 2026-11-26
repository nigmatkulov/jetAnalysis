"""Validated eta windows for inclusive-jet selection comparisons."""
from __future__ import annotations

import math

from hist_analysis.python.histogram_ops import normalize_histogram


def aligned_eta_bins(histogram, window):
    """Return in-range bins for an exact half-open, edge-aligned eta window."""
    low, high = window
    if not (math.isfinite(low) and math.isfinite(high) and low < high):
        raise ValueError('Eta window must have finite, ordered boundaries')
    axis = histogram.GetXaxis()
    edges = [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 2)]
    indices = []
    for value in (low, high):
        matches = [i for i, edge in enumerate(edges) if math.isclose(value, edge, rel_tol=0, abs_tol=1e-9)]
        if not matches:
            raise ValueError(f'Eta boundary {value:g} is not an in-range rebinned edge of {histogram.GetName()}')
        indices.append(matches[0])
    return indices[0] + 1, indices[1]


def normalize_eta_curves(curves, mode, window):
    """Normalize each curve in the same exact window; never include flows.

    Bin-width mode uses this same window and additionally divides by bin width.
    Stored-yield mode leaves contents and errors unchanged in detached clones.
    """
    if mode not in {'none', 'integral', 'bin_width'}:
        raise ValueError(f'Unsupported normalization mode: {mode!r}')
    normalized = {}
    for label, histogram in curves.items():
        if mode == 'none' or window is None:
            normalized[label] = normalize_histogram(histogram, mode)
            continue
        first, last = aligned_eta_bins(histogram, window)
        integral = histogram.Integral(first, last)
        if not math.isfinite(integral) or integral <= 0:
            raise ValueError(f'Nonpositive or nonfinite normalization for {label}: {integral}')
        clone = normalize_histogram(histogram, 'none')
        clone.Scale(1 / integral, 'width' if mode == 'bin_width' else '')
        normalized[label] = clone
    return normalized


def ratio_display_range(histogram, display_range, validity_range, normalization_range=None):
    """Intersect display/validity ranges, retaining only complete valid bins."""
    if validity_range is None:
        return display_range
    low, high = validity_range
    if not (math.isfinite(low) and math.isfinite(high) and low < high):
        raise ValueError('Validity range must have finite, ordered boundaries')
    axis = histogram.GetXaxis()
    valid = [i for i in range(1, axis.GetNbins() + 1)
             if axis.GetBinLowEdge(i) >= low - 1e-9
             and axis.GetBinUpEdge(i) <= high + 1e-9]
    if not valid:
        raise ValueError('Validity range contains no complete rebinned eta bins')
    first, last = valid[0], valid[-1]
    if normalization_range is not None:
        aligned_eta_bins(histogram, normalization_range)
        if not validity_range[0] <= normalization_range[0] < normalization_range[1] <= validity_range[1]:
            raise ValueError('Normalization window must lie inside the jetId validity range')
    axis = histogram.GetXaxis()
    low, high = axis.GetBinLowEdge(first), axis.GetBinUpEdge(last)
    if display_range is not None:
        low, high = max(low, display_range[0]), min(high, display_range[1])
    if low >= high:
        raise ValueError('Display range does not overlap the jetId validity range')
    return low, high


def selection_ratio_to_gen(numerator, gen, mode, window, *, option='', name=None):
    """Normalize rebinned inputs first, then divide with the selected ROOT option.

    Mode 'none' skips normalization. B is applied directly to the normalized
    histograms when normalization is enabled; no ratio-level scaling follows.
    Normalization-integral uncertainty/covariance is not included.
    """
    from hist_analysis.python.histogram_ops import ratio_to_nominal

    normalized = normalize_eta_curves({'numerator': numerator, 'Gen': gen}, mode, window)
    return ratio_to_nominal(normalized['numerator'], normalized['Gen'], name=name, option=option)


def project_selection_eta(histogram, pt_range, *, rebin=1, name=None):
    """Project inclusive-jet eta over exact, half-open pT bins, excluding pT flows.

    Eta flow contents/errors are retained, but never used by normalization.
    Reject partial or out-of-axis intervals instead of silently selecting bins.
    """
    from hist_analysis.python.projections import _semantic_axis_index

    if isinstance(rebin, bool) or not isinstance(rebin, int) or rebin < 1:
        raise ValueError('rebin must be a positive integer')
    target = _semantic_axis_index(histogram, 'eta')
    pt_axis = histogram.GetXaxis() if target == 1 else histogram.GetYaxis()
    # Reuse exact-edge validation on the integrated coordinate.
    class AxisView:
        def GetXaxis(self): return pt_axis
        def GetName(self): return histogram.GetName()
    first, last = aligned_eta_bins(AxisView(), pt_range)
    eta_axis = histogram.GetYaxis() if target == 1 else histogram.GetXaxis()
    if eta_axis.GetNbins() % rebin:
        raise ValueError('rebin must divide the eta-bin count')
    # Ignore display ranges without mutating the source histogram.
    source = histogram.Clone(f'{name or histogram.GetName()}_projection_source')
    source.SetDirectory(0)
    source.GetXaxis().SetRange(0, 0)
    source.GetYaxis().SetRange(0, 0)
    project = source.ProjectionY if target == 1 else source.ProjectionX
    result = project(name or f'{histogram.GetName()}_eta', first, last, 'e')
    result.SetDirectory(0)
    if rebin > 1:
        result.Rebin(rebin)
    return result


def independent_shape_ratio(numerator, gen, mode, window, *, name=None):
    """Delta-method shape errors under independent stored-bin/sample assumptions.

    Includes each random normalization sum and its covariance with its own bin.
    Omits shared-event inter-bin and numerator/Gen covariance, unavailable in TH1.
    This is an explicit alternative approximation, not a measured MC covariance.
    Zero-denominator bins follow ROOT's zero content/error convention.
    """
    result = selection_ratio_to_gen(numerator, gen, mode, window, name=name)
    if mode == 'none':
        return result
    def relative_variance(hist, i):
        first, last = ((1, hist.GetNbinsX()) if window is None
                       else aligned_eta_bins(hist, window))
        total = hist.Integral(first, last)
        variance = sum(hist.GetBinError(j)**2 for j in range(first, last + 1))
        value, error2 = hist.GetBinContent(i), hist.GetBinError(i)**2
        return error2/value**2 + variance/total**2 - (
            2*error2/(value*total) if first <= i <= last else 0)
    for i in range(result.GetNbinsX()+2):
        n, d = numerator.GetBinContent(i), gen.GetBinContent(i)
        if d == 0:
            continue
        if n == 0:
            # Derivative wrt numerator remains finite, even for cancellation.
            continue
        variance = relative_variance(numerator, i) + relative_variance(gen, i)
        result.SetBinError(i, abs(result.GetBinContent(i))*math.sqrt(max(0., variance)))
    return result


def ratio_variance(n, d, vn, vd, covariance=0.):
    """First-order raw ratio variance with an explicitly supplied covariance.

    Undefined denominator raises; callers must explicitly mask that bin.
    Covariance must obey the Cauchy bound. No covariance is inferred from yields.
    """
    if not all(math.isfinite(x) for x in (n, d, vn, vd, covariance)):
        raise ValueError('Ratio inputs must be finite')
    if d == 0 or vn < 0 or vd < 0:
        raise ValueError('Require nonzero denominator and nonnegative variances')
    if abs(covariance) > math.sqrt(vn*vd) + 1e-12:
        raise ValueError('Covariance exceeds the Cauchy bound')
    return max(0., (vn + (n/d)**2*vd - 2*(n/d)*covariance)/d**2)


def selection_diagnostics(raw_curves, mode, window):
    """Binwise legacy-B/independent/shape comparisons and weighted statistics."""
    gen = raw_curves['Gen']
    rows = []
    for label, numerator in raw_curves.items():
        if label == 'Gen':
            continue
        b = selection_ratio_to_gen(numerator, gen, mode, window, option='B')
        independent = selection_ratio_to_gen(numerator, gen, mode, window)
        shape = independent_shape_ratio(numerator, gen, mode, window)
        for i in range(1, gen.GetNbinsX()+1):
            n, d = numerator.GetBinContent(i), gen.GetBinContent(i)
            vn, vd = numerator.GetBinError(i)**2, gen.GetBinError(i)**2
            rows.append(dict(level=label, bin=i, eta=gen.GetBinCenter(i),
                numerator=n, denominator=d, numerator_sumw2=vn, gen_sumw2=vd,
                numerator_neff=n*n/vn if vn else 0., gen_neff=d*d/vd if vd else 0.,
                valid_denominator=d != 0, ratio=b.GetBinContent(i),
                error_B=b.GetBinError(i), error_independent=independent.GetBinError(i),
                error_shape_independent=shape.GetBinError(i)))
    return rows
