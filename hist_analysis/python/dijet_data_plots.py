"""Standalone and paginated CMS figures for existing reconstructed dijet TH1s."""

from dataclasses import replace
from pathlib import Path

import numpy as np

from hist_analysis.python.mplhep_plotting import (
    DEFAULT_MPLHEP_STYLE, _curve_style, _draw_cms_header, _graph_arrays, _mplhep, _th1_arrays,
)
from hist_analysis.python.output_naming import output_filename


_ETA = r"\eta_{\mathrm{CM}}^{\mathrm{dijet}}"
_Y_LABELS = {
    "none": "Dijet yield per bin",
    "integral": "Fraction of dijets per bin",
    "bin_width": rf"$(1/N)\,dN/d{_ETA}$",
    "fb": rf"Forward (${_ETA}>0$) / Backward (${_ETA}<0$)",
}


def _interval_label(interval):
    low, high = interval
    return rf"${low:g}\leq p_{{\mathrm{{T}}}}^{{\mathrm{{ave}}}}<{high:g}\ \mathrm{{GeV}}$"


def _limits(arrays, *, fb, headroom):
    # Include both ends of every displayed statistical error and the F/B line.
    lower, upper = (1.0, 1.0) if fb else (0.0, 0.0)
    for values, _, errors in arrays:
        valid = np.isfinite(values) & np.isfinite(errors)
        if valid.any():
            values, errors = np.asarray(values)[valid], np.asarray(errors)[valid]
            lower = min(lower, np.min(values - errors))
            upper = max(upper, np.max(values + errors))
    span = upper - lower or 1.0
    return lower - (0.05 if fb else 0.025) * span, upper + headroom * span


def save_combined_dijet_figures(
    histograms, *, observable, normalization, output_dir, selection, jet_eta_cut,
    annotations, cms_header_left, cms_header_right, individual_ylim=None, grid_notes=(),
    grid_ylim=None, totals=None, save=True, formats=("pdf", "png"), style=DEFAULT_MPLHEP_STYLE,
):
    """Render interval-to-TH1 mappings without changing the source histograms."""
    if observable not in {"eta", "fb"}:
        raise ValueError(f"Unknown dijet observable {observable!r}")
    if not histograms:
        raise ValueError("At least one dijet interval is required")
    if not save:
        return []
    if not formats or set(formats) - {"pdf", "png"}:
        raise ValueError("Choose PDF and/or PNG for the combined dijet figures")

    from hist_analysis.python.data_mc_plotting import accepted_eta_range

    for name, limits in (("individual_ylim", individual_ylim), ("grid_ylim", grid_ylim)):
        if limits is not None and (len(limits) != 2 or not np.isfinite(limits).all()
                                   or limits[0] >= limits[1]):
            raise ValueError(f"{name} must be a finite (lower, upper) pair")

    plt, hep = _mplhep()
    fb = observable == "fb"
    xlabel = rf"$|{_ETA}|$" if fb else rf"${_ETA}$"
    ylabel = _Y_LABELS["fb" if fb else normalization]
    intervals = sorted(histograms)
    arrays = {key: _th1_arrays(histograms[key]) for key in intervals}
    x_ranges = [accepted_eta_range(histograms[key], jet_eta_cut,
                                   observable="fb" if fb else "full") for key in intervals]
    xlim = (min(low for low, _ in x_ranges), max(high for _, high in x_ranges))
    brackets = {}
    if totals is not None:
        from hist_analysis.python.data_mc_plotting import data_total_band, draw_mpl_systematic_brackets

        brackets = {
            key: _graph_arrays(data_total_band(histograms[key], totals[key], name=f"data_syst_{observable}_{index}"))
            for index, key in enumerate(intervals)
        }

    def automatic_ylim(keys):
        inputs = [arrays[key] for key in keys]
        inputs.extend((brackets[key]["y"], None, brackets[key]["ylow"])
                      for key in keys if key in brackets)
        return _limits(inputs, fb=fb, headroom=0.10)

    common_grid_ylim = grid_ylim if grid_ylim is not None else (
        automatic_ylim(intervals)
    )
    label, color, marker, fill = _curve_style("Data", 0)
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    written = []

    def draw(axis, interval, *, grid):
        values, edges, errors = arrays[interval]
        widths = np.diff(edges)
        axis.set_xlim(*xlim)
        if fb:
            axis.axhline(1.0, ls="--", lw=1, color="0.45", zorder=1)
        handle = axis.errorbar(
            edges[:-1] + widths / 2, values, yerr=errors, xerr=widths / 2,
            fmt=marker, color=color, markerfacecolor=fill, label=f"{label} (stat.)",
            markersize=style.marker_size * (0.52 if grid else 0.73),
            capsize=style.capsize * (0.65 if grid else 1.0), elinewidth=1.1, zorder=2,
        )
        handles = [handle]
        if interval in brackets:
            draw_mpl_systematic_brackets(axis, brackets[interval], color)
            systematic = axis.lines[-1]
            systematic.set_label("Syst. uncertainty")
            handles.append(systematic)
        axis.tick_params(labelsize=style.tick_label_size * (0.80 if grid else 1))
        axis.xaxis.set_major_locator(plt.MaxNLocator(nbins=4 if grid else 7))
        axis.yaxis.set_major_locator(plt.MaxNLocator(nbins=4 if grid else 6))
        axis.minorticks_off()
        return handles

    def export(fig, suffix, **identity):
        try:
            for extension in dict.fromkeys(formats):
                name = output_filename(
                    f"combined_cm_dijet_{observable}", suffix,
                    selection=selection, eta_cut=jet_eta_cut,
                    extension=extension, **identity,
                )
                path = directory / name
                fig.savefig(path, dpi=180)
                written.append(path)
        finally:
            plt.close(fig)

    with plt.style.context(hep.style.CMS):
        for interval in intervals:
            fig, axis = plt.subplots(figsize=(style.figure_size,) * 2)
            fig.subplots_adjust(left=0.15, right=0.96, bottom=0.13, top=0.77)
            handles = draw(axis, interval, grid=False)
            axis.set_ylim(*(individual_ylim if individual_ylim is not None else
                            automatic_ylim([interval])))
            axis.set_xlabel(xlabel, fontsize=style.axis_title_size)
            axis.set_ylabel(ylabel, fontsize=style.axis_title_size * (0.8 if fb else 1))
            header = fig.add_axes((0.15, 0.918, 0.81, 0.025), frameon=False)
            header.set_axis_off()
            _draw_cms_header(hep, header, left=cms_header_left, right=cms_header_right,
                             style=replace(style, cms_label_size=style.cms_label_size * 0.85))
            fig.text(0.15, 0.915, "\n".join((_interval_label(interval), *annotations)),
                     va="top", fontsize=style.annotation_text_size * 0.86)
            fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.96, 0.915),
                       fontsize=style.legend_text_size * 0.83, frameon=False)
            export(fig, "individual", pt_interval=interval)

        for start in range(0, len(intervals), 16):
            page = intervals[start:start + 16]
            fig, axes = plt.subplots(4, 4, figsize=(style.figure_size * 1.85,) * 2,
                                     sharex=True, sharey=True)
            left, bottom, span = 0.105, 0.08, 0.84
            top = bottom + span
            fig.subplots_adjust(left=left, right=left + span, bottom=bottom,
                                top=top, wspace=0, hspace=0)
            # Keep the header just above the occupied top row.
            header_width = 0.42 if len(page) == 1 else span * min(len(page), 4) / 4
            header = fig.add_axes((left, top + 0.004, header_width, 0.006), frameon=False)
            header.set_axis_off()
            _draw_cms_header(hep, header, left=cms_header_left, right=cms_header_right,
                             style=style)
            for index, axis in enumerate(axes.flat):
                if index >= len(page):
                    axis.set_visible(False)
                    continue
                interval = page[index]
                handles = draw(axis, interval, grid=True)
                axis.set_ylim(*common_grid_ylim)
                axis.text(0.5, 0.99 if len(page) == 1 else 0.975,
                          _interval_label(interval), transform=axis.transAxes,
                          ha="center", va="top", fontsize=style.annotation_text_size)
                axis.tick_params(top=False, labelbottom=index + 4 >= len(page),
                                 labelleft=index % 4 == 0)
                axis.set_box_aspect(1)
                if index == min(1, len(page) - 1):
                    single = len(page) == 1
                    axis.legend(
                        handles=handles, loc="upper right" if single else "upper left",
                        bbox_to_anchor=(0.97, 0.72) if single else (0.03, 0.85),
                        fontsize=style.legend_text_size * 0.80, frameon=False,
                        borderaxespad=0,
                    )
            # Give a lone continuation panel room for both cut blocks and the legend.
            if len(page) == 1:
                axes.flat[0].set_position((left, top - 0.42, 0.42, 0.42))
            # Keep endpoint tick labels from touching across panel borders.
            axes[0, 0].xaxis.set_major_locator(plt.MaxNLocator(nbins=4, prune="both"))
            axes[0, 0].yaxis.set_major_locator(plt.MaxNLocator(nbins=4, prune="both"))
            note_panels = (4, 5) if len(page) >= 6 else (0, min(2, len(page) - 1))
            for group, (panel, notes) in enumerate(zip(note_panels, grid_notes)):
                right = (len(page) == 1 and group == 0) or (len(page) == 2 and group == 1)
                axis = axes.flat[panel]
                axis.text(
                    0.95 if right else 0.05, 0.90 if len(page) == 1 else 0.83,
                    "\n".join(notes), transform=axis.transAxes,
                    ha="right" if right else "left",
                    va="top",
                    fontsize=style.annotation_text_size * 0.90,
                )
            # Align the titles with the occupied frames on a partial last page.
            last = axes.flat[len(page) - 1].get_position()
            x_label_y = last.y0 - 0.055
            x_center = (axes[0, 0].get_position().x0 +
                        axes[0, min(len(page), 4) - 1].get_position().x1) / 2
            y_center = (axes[0, 0].get_position().y1 + last.y0) / 2
            fig.supxlabel(xlabel, x=x_center, y=x_label_y, fontsize=style.axis_title_size)
            short_page = len(page) <= 4
            grid_ylabel = ylabel.replace(" / ", " /\n") if fb and short_page else ylabel
            fig.supylabel(grid_ylabel, x=0.010 if fb and short_page else 0.026, y=y_center,
                          fontsize=style.axis_title_size * (0.78 if short_page else 1))
            suffix = "4x4" if len(intervals) <= 16 else f"4x4_page{start // 16 + 1:02}"
            export(fig, suffix, pt_bins=page)

    return written
