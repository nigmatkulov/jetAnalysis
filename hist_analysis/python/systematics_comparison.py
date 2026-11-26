"""Read and compare the stored relative dijet systematic uncertainties."""

import json
import re

import numpy as np
import pandas as pd

from hist_analysis.python.mplhep_plotting import (
    DEFAULT_MPLHEP_STYLE, _draw_cms_header, _graph_arrays, _mplhep,
)

OBSERVABLES = {"full": "full_eta", "fb": "forward_backward"}
SOURCES = ("JER", "JEU", "pileup", "pointing_resolution")
EDGE_ATOL = 1e-9  # Float roundoff only; no tolerance for physically different bins.


def _json(obj):
    return json.loads(str(obj.GetString() if obj.InheritsFrom("TObjString") else obj.GetTitle()))


def _bins(x, low, high, fraction):
    table = pd.DataFrame(dict(eta=x, low=low, high=high, percent=100 * fraction))
    if (table.empty or not np.isfinite(table).all().all()
            or (table.percent < 0).any()
            or not ((table.low < table.eta) & (table.eta < table.high)).all()
            or not np.allclose(table.high.to_numpy()[:-1], table.low.to_numpy()[1:],
                               atol=EDGE_ATOL, rtol=0)):
        raise ValueError("Invalid systematic uncertainty bins")
    return table


def read_published(archive, lab_to_cm_shift):
    """Use the measurement-band fractions; translate only full-eta coordinates."""
    provenance = _json(archive.Get("metadata/provenance"))
    result = {}
    for observable, directory in OBSERVABLES.items():
        for key in archive.GetDirectory(directory).GetListOfKeys():
            panel = archive.GetDirectory(f"{directory}/{key.GetName()}")
            pt = tuple(_json(panel.Get("provenance"))["pt_gev"])
            graph = panel.Get("data_syst_relative")
            if not graph or not graph.InheritsFrom("TGraphAsymmErrors"):
                raise ValueError(f"Missing systematic graph: {panel.GetPath()}")
            arrays = _graph_arrays(graph)
            central = _graph_arrays(panel.Get("data"))["y"]
            absolute = _graph_arrays(panel.Get("data_syst"))
            np.testing.assert_allclose(absolute["ylow"], absolute["yhigh"], rtol=1e-12)
            fraction = np.asarray(arrays["y"])
            np.testing.assert_allclose(np.abs(central) * fraction, absolute["ylow"], rtol=1e-12)
            eta = np.asarray(arrays["x"]) + (lab_to_cm_shift if observable == "full" else 0)
            result[observable, pt] = _bins(eta, eta - arrays["xlow"], eta + arrays["xhigh"], fraction)
    return result, provenance


def read_current(archive, expected):
    """Read notebook snapshots and verify their configuration, axes and totals."""
    metadata = _json(archive.Get("notebookConfiguration"))
    cfg = metadata["configuration"]
    for field, value in expected.items():
        if cfg[field] != value:
            raise ValueError(f"Current archive {field}={cfg[field]!r}, expected {value!r}")
    result = {}
    for key in archive.GetListOfKeys():
        name = key.GetName()
        match = re.fullmatch(r"systematics_(full|fb)_ptave__(\d+)_(\d+)_[0-9a-f]+", name)
        if not match:
            continue
        observable, low, high = match.groups()
        pt = int(low), int(high)
        if (observable, pt) in result or list(pt) not in cfg["PTAVE_BINS"]:
            raise ValueError(f"Unexpected or duplicate archive interval: {name}")
        directory = archive.GetDirectory(name)
        total = directory.Get("total_fraction")
        if not total or not total.InheritsFrom("TH1") or total.GetDimension() != 1:
            raise ValueError(f"Missing 1D total in {name}")
        bins = range(1, total.GetNbinsX() + 1)
        edges = np.array([total.GetBinLowEdge(i) for i in (*bins, total.GetNbinsX() + 1)])
        def contents(hist):
            if not hist or hist.GetNbinsX() != len(bins):
                raise ValueError(f"Missing or incompatible component in {name}")
            np.testing.assert_allclose([hist.GetBinLowEdge(i) for i in range(1, len(edges)+1)],
                                       edges, rtol=0, atol=EDGE_ATOL)
            return np.array([hist.GetBinContent(i) for i in bins])

        fraction = contents(total)
        components = [contents(directory.Get(f"{s}_fraction")) for s in SOURCES]
        np.testing.assert_allclose(fraction, np.linalg.norm(components, axis=0), rtol=1e-12, atol=1e-15)
        np.testing.assert_allclose(100 * fraction, contents(directory.Get("total_percent")),
                                   rtol=1e-12, atol=1e-13)
        eta_cut = cfg["ETA_CUT"]
        # The fiducial cut can fall inside the last histogram bin, just as in
        # systematic_total.load_source. Preserve the nominal histogram edges.
        left = 0 if observable == "fb" else -eta_cut
        if (edges[0] > left + EDGE_ATOL or edges[-1] < eta_cut - EDGE_ATOL
                or (edges[1:] <= left).any() or (edges[:-1] >= eta_cut).any()
                or (observable == "fb" and abs(edges[0]) > EDGE_ATOL)):
            raise ValueError(f"Histogram does not cover its declared acceptance: {name}")
        result[observable, pt] = _bins((edges[:-1] + edges[1:]) / 2, edges[:-1], edges[1:], fraction)
    return result, cfg


def compare_bins(old, new):
    """Return differences in percentage points only where both physical edges agree."""
    i, j = np.where(np.isclose(old.low.to_numpy()[:, None], new.low.to_numpy(), rtol=0, atol=EDGE_ATOL)
                    & np.isclose(old.high.to_numpy()[:, None], new.high.to_numpy(), rtol=0, atol=EDGE_ATOL))
    matched = old.iloc[i].reset_index(drop=True).rename(columns={"percent": "old_percent"})
    matched["new_percent"] = new.percent.to_numpy()[j]
    matched["delta_pp"] = matched.new_percent - matched.old_percent
    matched["abs_delta_pp"] = matched.delta_pp.abs()
    return matched


def plot_comparisons(observable, panels, *, grid, xlim, ylim, styles, marker_offset,
                     cuts, header):
    """Draw one overlay or a row-major page of up to 16 shared-axis overlays."""
    from dataclasses import replace
    from matplotlib.ticker import MaxNLocator
    from hist_analysis.python.systematic_total import (
        _place_total_annotations, _place_total_legend,
    )

    plt, hep = _mplhep()
    style = DEFAULT_MPLHEP_STYLE
    xlabel = (r'$|\eta_{\mathrm{CM}}^{\mathrm{dijet}}|$' if observable == "fb"
              else r'$\eta_{\mathrm{CM}}^{\mathrm{dijet}}$')
    ylabel = 'Relative total systematic uncertainty (%)'
    with plt.style.context(hep.style.CMS):
        if grid:
            fig, axes = plt.subplots(4, 4, figsize=(16, 16), sharex=True, sharey=True)
            left, right, bottom, top = .10, .99, .065, .955
            fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top,
                                wspace=0, hspace=0)
            notes = {}
            for index, group in zip((4, 5), cuts):
                index = min(index if index < len(panels) else index - 4, len(panels) - 1)
                notes[index] = (*notes.get(index, ()), *group)
        else:
            fig, axes = plt.subplots(figsize=(style.figure_size, style.figure_size),
                                     squeeze=False, layout='constrained')
        first = axes.flat[0]
        for index, (ax, (pt, old, new)) in enumerate(zip(axes.flat, panels)):
            # Match centers only for visual displacement; physical tables stay untouched.
            coincident = np.isclose(old.eta.to_numpy()[:, None], new.eta.to_numpy(),
                                    rtol=0, atol=EDGE_ATOL)
            for side, (label, table, match) in enumerate(
                    zip(("Old", "New"), (old, new),
                        (coincident.any(axis=1), coincident.any(axis=0)))):
                offset = marker_offset * (2 * side - 1) if observable == "full" else 0
                ax.plot(table.eta + offset * match, table.percent, label=label,
                        linestyle="none", markersize=style.marker_size if grid else 6,
                        **styles[side])
            ax.set(xlim=xlim, ylim=ylim)
            ax.grid(alpha=.18)
            ax.tick_params(labelsize=16)
            label = rf'${pt[0]:g} \leq p_{{\mathrm{{T}}}}^{{\mathrm{{ave}}}} < {pt[1]:g}$ GeV'
            if grid:
                ax.tick_params(axis='x', which='both', top=False,
                               labelbottom=index + 4 >= len(panels))
                ax.text(.5, .975, label, transform=ax.transAxes,
                        ha='center', va='top', fontsize=16)
                for line, cut in enumerate(notes.get(index, ())):
                    ax.text(.045, .87 - line * .105, cut, transform=ax.transAxes,
                            ha='left', va='top', fontsize=15)
            else:
                ax.set_xlabel(xlabel, fontsize=20)
                ax.set_ylabel(ylabel, fontsize=18)
                _draw_cms_header(hep, ax, left=header[0], right=header[1], style=style)
                annotations = [ax.text(.04, style.annotation_y - i * .045, text,
                                       transform=ax.transAxes, ha='left', va='top', fontsize=12)
                               for i, text in enumerate((label, *cuts[0], *cuts[1]))]
                _place_total_annotations(ax, annotations)
                _place_total_legend(ax)
        for ax in list(axes.flat)[len(panels):]:
            ax.set_visible(False)
        if grid:
            first.xaxis.set_major_locator(MaxNLocator(nbins=5, prune='both'))
            first.yaxis.set_major_locator(MaxNLocator(nbins=5, prune='both'))
            visible_bottom = top - ((len(panels) + 3) // 4) * (top - bottom) / 4
            legend_ax = axes.flat[min(1, len(panels) - 1)]
            legend_ax.legend(loc='upper right', bbox_to_anchor=(.97, .89),
                             ncols=2, fontsize=16, frameon=False)
            fig.supxlabel(xlabel, y=visible_bottom - .056, fontsize=27)
            fig.supylabel(ylabel, x=.04, y=(top + visible_bottom) / 2, fontsize=26)
            header_axis = fig.add_axes((left, top + .004, right - left, .001), frameon=False)
            header_axis.set_axis_off()
            _draw_cms_header(hep, header_axis, left=header[0], right=header[1],
                             style=replace(style, cms_label_size=27))
        return fig


def save_mplhep_comparison(fig, path, *, grid, save_png, formats, partial_grid=False):
    """Export through the project MPLHEP PDF/PNG convention."""
    from hist_analysis.python.mplhep_plotting import _save

    path.parent.mkdir(parents=True, exist_ok=True)
    if not grid and "pdf" in formats:
        _save(fig, path, save_png=save_png and "png" in formats, style=DEFAULT_MPLHEP_STYLE)
        return
    options = dict(dpi=160 if grid else 200, bbox_inches="tight" if partial_grid else None)
    if "pdf" in formats:
        fig.savefig(path, **options)
    if save_png and "png" in formats:
        fig.savefig(path.with_suffix(".png"), **options)
