"""Compare inclusive data jets using each JEC stage's own pT coordinate."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import ROOT

from hist_analysis.python.histogram_io import resolve_data_file
from hist_analysis.python.mplhep_plotting import (
    DEFAULT_MPLHEP_STYLE, _draw_cms_header, _mplhep,
)
from hist_analysis.python.plot_labels import CMS_PRELIMINARY_HEADER, beam_orientation_label

RAW = "hRecoInclusiveJetRawPtEtaLabUnflipped"
CORRECTED = "hRecoInclusiveJetPtEtaLabUnflipped"
REFERENCE = "l2_relative"
VARIANTS = {
    REFERENCE: "_noResiduals",
    "l2l3_residual_old": "_oldResiduals",
    "l2l3_residual_new": "",
}
STYLES = {
    "raw": ("Raw", "0.5", "o"),
    REFERENCE: ("L2-Rel.", "black", "s"),
    "l2l3_residual_old": ("L2L3-Res.Old", "tab:blue", "^"),
    "l2l3_residual_new": ("L2L3-Res.New", "tab:red", "D"),
}


@dataclass
class Histogram:
    edges: tuple
    values: np.ndarray
    variance: np.ndarray
    entries: float

    def matches(self, other):
        return (
            self.entries == other.entries
            and all(np.array_equal(a, b) for a, b in zip(self.edges, other.edges))
            and np.array_equal(self.values, other.values)
            and np.array_equal(self.variance, other.variance)
        )


def _read(file, key, dimension):
    h = file.Get(key)
    context = f"{file.GetName()}:{key}"
    if not h or not h.InheritsFrom("TH1") or h.GetDimension() != dimension:
        raise ValueError(f"{context}: expected a {dimension}D histogram")
    if h.GetSumw2N() != h.GetNcells():
        raise ValueError(f"{context}: expected stored Sumw2")
    axes = (h.GetXaxis(), h.GetYaxis())[:dimension]
    edges = tuple(np.array([axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 2)])
                  for axis in axes)
    # ROOT global bins are x-fastest, including underflow and overflow.
    shape = tuple(len(edge) + 1 for edge in edges)
    values = np.array([h.GetBinContent(i) for i in range(h.GetNcells())]).reshape(shape, order="F")
    sumw2 = h.GetSumw2()
    variance = np.array([sumw2.At(i) for i in range(h.GetNcells())]).reshape(shape, order="F")
    if (not np.all(np.isfinite(values)) or not np.all(np.isfinite(variance))
            or np.any(values < 0) or np.any(variance < 0)):
        raise ValueError(f"{context}: invalid data counts or variances")
    return Histogram(edges, values, variance, h.GetEntries())


def load_corrections(data_dir, trigger, direction):
    """Read original campaign files and exclude an old sample if its controls differ.

    Equal aggregate controls are necessary checks, not event-by-event identity.
    The variant-to-payload mapping comes from the producer and campaign naming.
    """
    base = resolve_data_file(data_dir, trigger, direction, "jetId", prefer_produced=False)
    inputs = {}
    for level, suffix in VARIANTS.items():
        path = base.with_stem(base.stem + suffix)
        with ROOT.TFile.Open(str(path), "READ") as file:
            inputs[level] = [_read(file, key, dim) for key, dim in (
                (RAW, 2), (CORRECTED, 2), ("hVz", 1),
            )]
        raw, corrected, _ = inputs[level]
        if raw.entries != corrected.entries:
            raise ValueError(f"{path}: raw and corrected inclusive entries differ")

    raw, ref, vz = inputs[REFERENCE]
    result = {"raw": raw, REFERENCE: ref}
    excluded = {}
    for level in VARIANTS:
        if level == REFERENCE:
            continue
        candidate_raw, corrected, candidate_vz = inputs[level]
        mismatches = [key for key, a, b in (
            (RAW, raw, candidate_raw), ("hVz", vz, candidate_vz),
        ) if not a.matches(b)]
        if mismatches:
            reason = f"{level}: differs from L2-Rel. in {', '.join(mismatches)}"
            if level != "l2l3_residual_old":
                raise ValueError(reason)
            excluded[level] = reason
        else:
            result[level] = corrected

    if any(not np.array_equal(h.edges[1], ref.edges[1]) for h in result.values()):
        raise ValueError(f"{trigger}/{direction}: incompatible eta axes")
    return result, excluded


def _bin_slice(edges, interval):
    indices = []
    for boundary in interval:
        matches = np.flatnonzero(np.isclose(edges, boundary, rtol=0, atol=1e-10))
        if len(matches) != 1:
            raise ValueError(f"Boundary {boundary} is not a histogram edge")
        indices.append(int(matches[0]))
    if indices[0] >= indices[1]:
        raise ValueError(f"Invalid interval {interval}")
    return slice(*indices)


def project_eta(histogram, pt_range, eta_range):
    """Project regular bins only; reject cuts that would split an input bin."""
    pt_edges, eta_edges = histogram.edges
    pt = _bin_slice(pt_edges, pt_range)
    eta = _bin_slice(eta_edges, eta_range)
    values = histogram.values[1:-1, 1:-1][pt, eta].sum(axis=0)
    variance = histogram.variance[1:-1, 1:-1][pt, eta].sum(axis=0)
    return eta_edges[eta.start:eta.stop + 1], values, variance


def ratio_to_l2(numerator, reference):
    """Return the ratio and delta-method variance with covariance set to zero."""
    n, vn = numerator
    d, vd = reference
    ratio, variance = np.full_like(d, np.nan), np.full_like(d, np.nan)
    valid = d > 0
    ratio[valid] = n[valid] / d[valid]
    variance[valid] = (vn[valid] + ratio[valid]**2 * vd[valid]) / d[valid]**2
    return ratio, variance


def diagnostics(edges, counts, ratios):
    """Return yields, relative changes, and eta dependence without significance claims."""
    eta = (edges[:-1] + edges[1:]) / 2
    ref = counts[REFERENCE][0]
    rows = []
    for level, (n, _) in counts.items():
        row = {"level": level, "yield": float(n.sum())}
        for label, mask in (
            ("change", np.ones_like(eta, dtype=bool)),
            ("central_change", abs(eta) < 1.9),
            ("outer_change", abs(eta) >= 1.9),
        ):
            denom = ref[mask].sum()
            row[label] = float(n[mask].sum() / denom - 1) if denom else None
        if level in ratios:
            ratio = ratios[level][0]
            if np.isfinite(ratio).any():
                peak = np.nanargmax(abs(ratio - 1))
                row.update(peak_eta_low=float(edges[peak]), peak_eta_high=float(edges[peak + 1]),
                           peak_change=float(ratio[peak] - 1), peak_l2_yield=float(ref[peak]))
        rows.append(row)
    return rows


def plot_comparison(edges, counts, ratios, *, trigger, direction, pt_range, excluded,
                    output_dir, eta_ylim, ratio_ylim, formats):
    """Save standalone count and ratio figures with consistent CMS styling."""
    plt, hep = _mplhep()
    eta = (edges[1:] + edges[:-1]) / 2
    half_width = np.diff(edges) / 2
    low, high = pt_range
    orientation = beam_orientation_label(direction) or "Combined beams"
    stem = (f"{trigger}_{direction}_jetId_etaLabUnflipped_{edges[0]:g}_{edges[-1]:g}"
            f"_pt_{low:g}_{high:g}").replace(".", "p")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    with plt.style.context(hep.style.CMS):
        for observable, series, ylabel, ylim in (
            ("eta_overlay", counts, "Jets / bin", eta_ylim),
            ("ratio_to_l2", ratios, "Jets / L2-Rel.", ratio_ylim),
        ):
            is_ratio = observable == "ratio_to_l2"
            fig, ax = plt.subplots(figsize=(8, 7), layout="constrained")
            try:
                _draw_cms_header(hep, ax, left=CMS_PRELIMINARY_HEADER,
                                 right="pPb (8.16 TeV)", style=DEFAULT_MPLHEP_STYLE)
                if is_ratio:
                    ax.axhline(1, color="0.45", ls="--", lw=1)
                for level, (values, variance) in series.items():
                    label, color, marker = STYLES[level]
                    valid = np.isfinite(values) & (values != 0)
                    ax.errorbar(eta[valid], values[valid], xerr=half_width[valid],
                                yerr=np.sqrt(variance[valid]), label=label, color=color,
                                fmt=marker, markersize=4, lw=1)
                detail = (f"{trigger} • {orientation}\n"
                          f"AK4 PF, TightLepVeto jet ID, "
                          rf"${edges[0]:g}\leq\eta_{{\mathrm{{lab}}}}<{edges[-1]:g}$"
                          f"\n${low:g}\\leq p_{{\\mathrm{{T}}}}<{high:g}$ GeV")
                if excluded:
                    detail += "\nOld residual omitted: sample mismatch"
                if is_ratio:
                    detail += "\nStat. bars: covariance set to zero"
                # Put selection text above the panel to avoid covering data.
                ax.set_title(detail, loc="left", fontsize=12, pad=40)
                ax.set(xlim=(edges[0], edges[-1]), ylim=ylim,
                       xlabel=r"$\eta_{\mathrm{lab}}^{\mathrm{jet}}$ (unflipped)", ylabel=ylabel)
                ax.legend(loc="best", fontsize=13)
                if not is_ratio and ylim is None:
                    ax.set_ylim(bottom=0)
                for fmt in formats:
                    path = output_dir / f"{stem}_{observable}.{fmt}"
                    fig.savefig(path)
                    paths.append(path)
            finally:
                plt.close(fig)
    return paths
