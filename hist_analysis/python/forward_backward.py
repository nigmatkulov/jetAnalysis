"""Derive unnormalized eta hemispheres from symmetric full ROOT histograms."""

from array import array
import math

import ROOT


def fold_eta_histogram(full, *, name: str):
    """Return detached (Forward, Backward) TH1/TH2 with copied Sumw2 errors.

    TH2 inputs have pT on x and eta on y. Fold before eta rebinning. Zero must
    be an edge, with mirrored bins on both sides; no interpolation is allowed.
    Eta flow bins become the corresponding hemisphere overflow. Interior
    exact-edge entries follow the full histogram's ROOT bin assignment.
    """
    # Folding is a permutation/copy of whole bins, so stored weighted errors
    # survive even for signed weights. SetBinContent does not reconstruct the
    # original event count or exact within-bin moments; use integrals for yields.
    dimension = full.GetDimension()
    if dimension not in (1, 2):
        raise TypeError("Eta folding requires a TH1 or TH2")
    axis = full.GetYaxis() if dimension == 2 else full.GetXaxis()
    n = axis.GetNbins()
    half = n // 2
    edges = [axis.GetBinLowEdge(i) for i in range(1, n + 2)]
    if (n % 2 or not math.isclose(edges[half], 0., abs_tol=1e-9)
            or any(not math.isclose(edges[i], -edges[n-i], rel_tol=1e-10,
                                    abs_tol=1e-9) for i in range(n+1))):
        raise ValueError(f"{full.GetName()}: eta bins must be symmetric with zero as an edge")
    positive_edges = array('d', [0., *edges[half+1:]])
    sides = []
    for label in ('forward', 'backward'):
        if dimension == 1:
            result = ROOT.TH1D(f'{name}_{label}', full.GetTitle(), half, positive_edges)
        else:
            x = full.GetXaxis()
            x_edges = array('d', [x.GetBinLowEdge(i) for i in range(1, x.GetNbins()+2)])
            result = ROOT.TH2D(f'{name}_{label}', full.GetTitle(), x.GetNbins(), x_edges,
                               half, positive_edges)
        result.SetDirectory(0)
        result.Sumw2()
        result.GetXaxis().SetTitle(full.GetXaxis().GetTitle() if dimension == 2 else '|#eta_{CM}|')
        result.GetYaxis().SetTitle('|#eta_{CM}|' if dimension == 2 else full.GetYaxis().GetTitle())
        sides.append(result)
    forward, backward = sides
    for output_bin in range(1, half+2):
        fbin = half + output_bin if output_bin <= half else n+1
        bbin = half - output_bin + 1 if output_bin <= half else 0
        for result, source_bin in ((forward, fbin), (backward, bbin)):
            if dimension == 1:
                result.SetBinContent(output_bin, full.GetBinContent(source_bin))
                result.SetBinError(output_bin, full.GetBinError(source_bin))
            else:
                for xbin in range(full.GetNbinsX()+2):
                    result.SetBinContent(xbin, output_bin, full.GetBinContent(xbin, source_bin))
                    result.SetBinError(xbin, output_bin, full.GetBinError(xbin, source_bin))
    return forward, backward


def forward_backward_ratio(forward, backward, *, name: str):
    """Divide raw hemispheres using ROOT independent-error propagation."""
    ratio = forward.Clone(name)
    ratio.SetDirectory(0)
    # F and B are disjoint hemispheres, not pass/total. ROOT's default assumes
    # independent bins; a zero B produces content=error=0 (undefined physics
    # ratio, not evidence of a measured zero). No cross-bin covariance is stored.
    ratio.Divide(forward, backward, 1., 1., '')
    ratio._forward_backward_inputs = (forward, backward)
    return ratio
