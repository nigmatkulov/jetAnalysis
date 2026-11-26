"""CMS diagnostics for the joint pTave–eta validation, with no yield rescaling."""

from dataclasses import replace
from pathlib import Path

import numpy as np

from hist_analysis.python.dijet_data_plots import _interval_label, _limits
from hist_analysis.python.mplhep_plotting import (
    DEFAULT_MPLHEP_STYLE, _curve_style, _draw_cms_header, _mplhep,
)
from hist_analysis.python.unfolding_validation import edges, matrix_values, values

ETA = r"\eta_{\mathrm{CM}}^{\mathrm{dijet}}"
KINDS = {
    "eta": ("Weighted yield / bin", False),
    "to_gen": ("Ratio to Gen", True),
    "refolded": ("Weighted yield / bin", False),
    "refolded_ratio": ("Refolded / Reco", True),
    "matched_refolded": ("Weighted yield / bin", False),
    "matched_refolded_ratio": ("Matched refolded / corrected Reco", True),
    "fb": (rf"Forward (${ETA}>0$) / Backward (${ETA}<0$)", True),
    "fb_double": ("Ratio of Forward / Backward", True),
    "method_comparison": ("Method B / Method A", True),
    "purity": ("Reco matched fraction", False),
    "inverse_efficiency": ("Inverse reconstruction efficiency", False),
    "efficiency": ("Reconstruction efficiency", False),
    "bin_migration": ("Diagonal fraction of matched yield", False),
    "unfolded_over_reco": ("Unfolded / Reco", True),
    "original": ("Weighted yield / bin", False),
}


class ValidationPlots:
    def __init__(self, directory, formats, *, left, right, annotations, limits=None,
                 style=DEFAULT_MPLHEP_STYLE, figsize=None):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.formats = tuple(dict.fromkeys(formats))
        if not self.formats or set(self.formats) - {"pdf", "png"}:
            raise ValueError("Choose PDF and/or PNG")
        self.left, self.right, self.annotations = left, right, tuple(annotations)
        self.limits, self.style = limits or {}, style
        self.figsize = figsize
        self.plt, self.hep = _mplhep()
        self.paths = []

    def _export(self, fig, stem):
        try:
            for ext in self.formats:
                path = self.directory / f"{stem}.{ext}"
                fig.savefig(path, dpi=155)
                self.paths.append(path)
        finally:
            self.plt.close(fig)

    def _header(self, fig, *, left=None, top=.90, right=.96, small=False):
        left = .16 if left is None else left
        axis = fig.add_axes((left, top, right-left, .008), frameon=False)
        axis.set_axis_off()
        _draw_cms_header(self.hep, axis, left=self.left, right=self.right,
                         style=replace(self.style, cms_label_size=15 if small else 20))

    @staticmethod
    def _arrays(s):
        h = s.histogram
        y = values(h)
        error = np.sqrt(np.maximum(s.covariance.diagonal(), 0))
        return np.where(s.supported & (y != 0), y, np.nan), edges(h.GetXaxis()), error

    def _draw(self, axis, curves, *, unity, grid):
        if unity:
            axis.axhline(1., color="0.45", lw=1, ls="--")
        for index, (label, spectrum) in enumerate(curves.items()):
            y, bins, error = self._arrays(spectrum)
            source = label.split(" / ")[0]
            if source in {"Reco", "Corrected Reco"}:
                color, marker, face = "tab:red", "o", "tab:red"
            elif source == "Gen":
                color, marker, face = "tab:blue", "o", "none"
            elif source == "Unfolded":
                color, marker, face = "black", "s", "black"
            else:
                _, color, marker, face = _curve_style(label, index + 1)
            widths = np.diff(bins)/2
            axis.errorbar(bins[:-1]+widths, y, xerr=widths, yerr=error, fmt=marker,
                          color=color, markerfacecolor=face, label=label, lw=.9,
                          markersize=self.style.marker_size*(.46 if grid else .7),
                          capsize=self.style.capsize*(.55 if grid else 1))
        axis.tick_params(labelsize=11 if grid else 13, top=False, right=False)

    def comparisons(self, kind, flat, intervals, *, prefix, grids=True):
        ylabel, unity = KINDS[kind]
        plt = self.plt
        xlabel = rf"$|{ETA}|$" if kind in {"fb", "fb_double"} else rf"${ETA}$"
        with plt.style.context(self.hep.style.CMS):
            cases = [(f"{prefix}_{kind}_pt{lo:g}_{hi:g}", curves, _interval_label((lo, hi)),
                      xlabel) for (lo, hi), curves in intervals.items()]
            if flat is not None:
                cases.insert(0, (f"{prefix}_{kind}_flat", flat, None, "Flattened ROOT bin"))
            for stem, curves, interval, x_title in cases:
                fig, axis = plt.subplots(figsize=self.figsize or (7, 7))
                fig.subplots_adjust(left=.17, right=.96, bottom=.13, top=.74)
                self._header(fig, left=.17, top=.915, small=True)
                self._draw(axis, curves, unity=unity, grid=False)
                axis.set_xlabel(x_title, fontsize=15)
                axis.set_ylabel(ylabel, fontsize=13)
                if kind in {"eta", "refolded", "matched_refolded", "original"}:
                    axis.yaxis.get_offset_text().set_size(11)
                axis.set_ylim(*(self.limits.get(kind) or
                                _limits([self._arrays(s) for s in curves.values()],
                                        fb=unity, headroom=.12)))
                if interval is None and intervals:
                    width = next(iter(next(iter(intervals.values())).values())).histogram.GetNbinsX()
                    for i in range(1, len(intervals)):
                        axis.axvline(i*width+.5, color=".75", lw=.6, zorder=0)
                notes = (*self.annotations, interval) if interval else self.annotations
                fig.text(.17, .90, "\n".join(notes), va="top", fontsize=11)
                handles, labels = axis.get_legend_handles_labels()
                fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(.97,.90),
                           frameon=False, fontsize=11)
                self._export(fig, stem)
            if not intervals or not grids:
                return
            limit = self.limits.get(kind) or _limits(
                [self._arrays(s) for curves in intervals.values() for s in curves.values()],
                fb=unity, headroom=.25)
            for start in range(0, len(intervals), 16):
                page = list(intervals.items())[start:start+16]
                fig, axes = plt.subplots(4, 4, figsize=(13, 13), sharex=True, sharey=True)
                fig.subplots_adjust(left=.10, bottom=.08, top=.90, right=.92, wspace=0, hspace=0)
                self._header(fig, left=.10, top=.96, right=.92)
                fig.text(.10, .953, "   ".join(self.annotations), va="top", fontsize=13)
                for i, axis in enumerate(axes.flat):
                    if i >= len(page):
                        axis.set_visible(False)
                        continue
                    interval, curves = page[i]
                    self._draw(axis, curves, unity=unity, grid=True)
                    axis.set_ylim(*limit)
                    axis.xaxis.set_major_locator(plt.MaxNLocator(4, prune="both"))
                    axis.yaxis.set_major_locator(plt.MaxNLocator(4, prune="both"))
                    axis.tick_params(labelleft=i % 4 == 0, labelbottom=i+4 >= len(page))
                    axis.yaxis.get_offset_text().set_visible(False)
                    axis.text(.5, .98, _interval_label(interval), transform=axis.transAxes,
                              va="top", ha="center", fontsize=12)
                handles, labels = axes.flat[0].get_legend_handles_labels()
                fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(.51, .912),
                           ncol=len(labels), frameon=False, fontsize=13)
                fig.supxlabel(xlabel, y=.03, fontsize=18)
                if kind in {"eta", "refolded", "matched_refolded", "original"}:
                    axes.flat[0].yaxis.get_offset_text().set_visible(True)
                    axes.flat[0].yaxis.get_offset_text().set_size(11)
                fig.supylabel(ylabel, x=.02, fontsize=16)
                self._export(fig, f"{prefix}_{kind}_grid{start//16+1}")

    def response(self, training):
        from matplotlib.colors import LogNorm

        layout = training.layout
        flat = np.arange(layout.n_global_bins+1)+.5
        eta = edges(training.templates[0].GetXaxis())
        pt = np.arange(layout.n_pt_bins+1)+.5
        m = matrix_values(training.matrix)
        a = training.probabilities
        block = m.reshape(layout.n_pt_bins, layout.n_eta_bins,
                          layout.n_pt_bins, layout.n_eta_bins).sum(axis=(1, 3))

        def draw(z, bins, xlabel, ylabel, colorlabel, stem, note="", boundaries=(), ticks=()):
            with self.plt.style.context(self.hep.style.CMS):
                fig, ax = self.plt.subplots(figsize=self.figsize or (8, 7))
                fig.subplots_adjust(left=.15, bottom=.13, top=.78, right=.86)
                self._header(fig, left=.15, top=.935, right=.96, small=True)
                cmap = self.plt.get_cmap("viridis").with_extremes(bad="none")
                nonzero = z[z > 0]
                norm = LogNorm(vmin=nonzero.min(), vmax=max(nonzero.max(), 1.01*nonzero.min())) if nonzero.size else None
                mesh = ax.pcolormesh(bins, bins, np.ma.masked_less_equal(z.T, 0),
                                     cmap=cmap, norm=norm, rasterized=True)
                if nonzero.size:
                    fig.colorbar(mesh, ax=ax, pad=.025).set_label(colorlabel, fontsize=13)
                ax.set(xlabel=xlabel, ylabel=ylabel)
                ax.xaxis.label.set_size(15)
                ax.yaxis.label.set_size(15)
                for b in boundaries:
                    ax.axvline(b, color="white", lw=.45, alpha=.65)
                    ax.axhline(b, color="white", lw=.45, alpha=.65)
                if len(ticks):
                    ax.set_xticks(ticks)
                    ax.set_yticks(ticks)
                ax.tick_params(labelsize=10)
                fig.text(.15, .905, "\n".join((*self.annotations, note)), va="top", fontsize=9)
                self._export(fig, stem)

        boundaries = flat[layout.n_eta_bins:-1:layout.n_eta_bins]
        draw(m, flat, "Reco flattened ROOT bin", "Gen flattened ROOT bin", "Matched weighted yield",
             "response_weighted", "Lines mark pTave block boundaries", boundaries)
        draw(a, flat, "Reco flattened ROOT bin", "Gen flattened ROOT bin",
             "Probability given matched Gen bin", "response_probability",
             "Sum over Reco = 1 on supported Gen bins", boundaries)
        draw(block, pt, "Reco pTave interval index", "Gen pTave interval index",
             "Matched weighted yield", "response_pt_migration",
             "Includes off-diagonal pTave migrations", ticks=np.arange(1, len(pt)))
        for i, (lo, hi) in enumerate(layout.pt_intervals):
            s = slice(i*layout.n_eta_bins, (i+1)*layout.n_eta_bins)
            draw(m[s, s], eta, rf"Reco ${ETA}$", rf"Gen ${ETA}$", "Matched weighted yield",
                 f"response_pt{lo:g}_{hi:g}", f"Reco and Gen pTave: [{lo:g}, {hi:g}) GeV; diagonal block only")
