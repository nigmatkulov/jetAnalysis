"""MPLHEP grids for the existing percent systematic-uncertainty displays."""
from dataclasses import replace
from pathlib import Path

import numpy as np

from hist_analysis.config.histograms import DIJET_PTAVE_BINS
from hist_analysis.python.mplhep_plotting import (
    DEFAULT_MPLHEP_STYLE, _curve_style, _draw_cms_header, _mplhep,
    _root_label, _root_pad_range, _th1_arrays,
)
from hist_analysis.python.plot_labels import CMS_SIMULATION_HEADER, PPB_MPLHEP_HEADER_RIGHT


def _page_y_range(curves, x_range, y_range, *, annotation_lines):
    extrema = [0.0]
    for values, edges, errors in curves:
        edges = np.asarray(edges)
        centers = (edges[1:] + edges[:-1]) / 2
        values = np.asarray(values)
        visible = (centers >= x_range[0]) & (centers <= x_range[1])
        lower, upper = np.broadcast_to(0 if errors is None else errors, (2, len(values)))
        for bound in (values - lower, values + upper):
            extrema.extend(bound[visible & np.isfinite(bound)])
    low, high = min(extrema), max(extrema)
    padding = (high - low or 1.0) * 0.08
    # Reserve space above the error bars for one or two lines of panel text.
    headroom = (high - low or 1.0) * (0.15 + 0.09 * annotation_lines)
    natural_range = (low - padding if low < 0 else 0, high + headroom)
    if y_range is None:
        return natural_range
    return min(y_range[0], natural_range[0]), max(y_range[1], natural_range[1])


def save_systematic_grid(source, panels, dijet_ptave_bins, *, output, x_title,
                         observable='dijet shape', x_range=None, y_range=None,
                         sample_label=None, eta_phi_cuts=(), secondary_cuts=(), panel_labels=None,
                         cms_header_left=CMS_SIMULATION_HEADER, cms_header_right=PPB_MPLHEP_HEADER_RIGHT,
                         save_png=False, stepped=False, canvas_view=None, fits=None, raw_ratios=None,
                         ratio_title=None, stat_scale=1.0, contributions=None):
    """Save row-major 4×4 pages from TH1s or (values, edges, errors) tuples.

    Inputs are already in display units. Errors may be None or symmetric or
    asymmetric arrays accepted by Matplotlib. The configured y range is
    expanded per page when necessary to show all visible points and their errors.
    """
    bins = sorted(map(tuple, dijet_ptave_bins))
    if not bins or len(bins) != len(set(bins)) or set(bins) != set(panels):
        raise ValueError('Grid panels must match unique configured ptAve intervals')

    # Match the primary grid in the dijet closure notebook. Keep every other
    # configured interval on supplemental pages, including overlapping nPDF bins.
    standard = DIJET_PTAVE_BINS[:16]
    primary = standard if set(standard).issubset(panels) else []
    remaining = [interval for interval in bins if interval not in primary]
    pages = [primary] if primary else []
    pages.extend(remaining[start:start + 16] for start in range(0, len(remaining), 16))

    if canvas_view is None:
        curves = {
            interval: _th1_arrays(panel) if hasattr(panel, 'GetNbinsX') else panel
            for interval, panel in panels.items()
        }
        if x_range is None:
            x_range = (min(edges[0] for _, edges, _ in curves.values()),
                       max(edges[-1] for _, edges, _ in curves.values()))
    else:
        from hist_analysis.python.systematic_plotting import _draw_pad

        if x_range is None:
            limits = [_root_pad_range(canvas, 'x') for canvas in panels.values()]
            x_range = min(low for low, _ in limits), max(high for _, high in limits)
    if contributions:
        from hist_analysis.python.systematic_total import draw_systematic_contributions
    plt, hep = _mplhep()
    style = replace(DEFAULT_MPLHEP_STYLE, cms_label_size=27)
    _, color, marker, face = _curve_style(source, 0)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    paths = []
    with plt.style.context(hep.style.CMS):
        for page_number, page in enumerate(pages, 1):
            left, right, bottom, top = .095, .985, .070, .960
            if canvas_view is not None or contributions:
                left, top = .145, .910  # Space for the shared legend above square panels.
            fig, axes = plt.subplots(4, 4, figsize=(16, 16), sharex=True, sharey=True)
            try:
                fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top,
                                    wspace=0, hspace=0)
                rows = (len(page) + 3) // 4
                visible_bottom = top - rows * (top - bottom) / 4
                notes = {0: (sample_label,) if sample_label else ()}
                for index, cuts in ((4, eta_phi_cuts), (5, secondary_cuts)):
                    if cuts:
                        if index >= len(page):
                            index -= 3
                        index = min(index, len(page) - 1)
                        notes[index] = (*notes.get(index, ()), *cuts)
                legend = {}
                for index, ax in enumerate(axes.flat):
                    if index >= len(page):
                        ax.set_visible(False)
                        continue
                    interval = page[index]
                    if canvas_view is not None:
                        x_title, canvas_y_title = _draw_pad(
                            ax, panels[interval], hep, cms_header_left=cms_header_left,
                            grid_view=canvas_view, fits=(fits or {}).get(interval, {}),
                            raw_ratios=None if raw_ratios is None else raw_ratios[interval],
                            stat_scale=stat_scale,
                        )
                    else:
                        values, edges, errors = curves[interval]
                        if stepped:
                            if contributions:
                                draw_systematic_contributions(ax, hep, contributions[interval],
                                                              edges, style=style)
                            ax.stairs(values, edges, color='black', linewidth=2, zorder=4)
                            centers = (np.asarray(edges[1:]) + edges[:-1]) / 2
                            ax.errorbar(centers, values, yerr=errors, fmt='ks', markersize=5,
                                        capsize=style.capsize, zorder=5,
                                        label='Total' if contributions else '_nolegend_')
                        else:
                            hep.histplot(values, bins=edges, yerr=False if errors is None else errors,
                                         xerr=True, histtype='errorbar', ax=ax, color=color,
                                         marker=marker, markerfacecolor=face, markeredgecolor=color,
                                         markersize=style.marker_size, capsize=style.capsize)
                    if canvas_view is not None or contributions:
                        handles, labels = ax.get_legend_handles_labels()
                        for label, handle in zip(labels, handles):
                            legend.setdefault(label, handle)
                    low_pt, high_pt = interval
                    label = rf'${low_pt:g} \leq p_{{\mathrm{{T}}}}^{{\mathrm{{ave}}}} < {high_pt:g}$ GeV'
                    if panel_labels and interval in panel_labels:
                        label += '\n' + panel_labels[interval]
                    ax.text(.5, .975, label, transform=ax.transAxes, ha='center', va='top',
                            fontsize=16, bbox=dict(facecolor='white', edgecolor='none', pad=2))
                    note_y = .78 if panel_labels and interval in panel_labels else .87
                    note_spacing = .105 if canvas_view is not None or contributions else .09
                    for line, note in enumerate(notes.get(index, ())):
                        ax.text(.045, note_y - note_spacing * line, _root_label(note),
                                transform=ax.transAxes, ha='left', va='top', fontsize=16)
                    ax.tick_params(axis='x', which='both', top=False,
                                   labelsize=16, labelbottom=index + 4 >= len(page))
                    ax.tick_params(axis='y', labelsize=16)
                    ax.grid(alpha=.18)
                if canvas_view is None:
                    page_limits = _page_y_range((curves[interval] for interval in page),
                                                x_range, y_range,
                                                annotation_lines=max(map(len, notes.values())) + bool(panel_labels))
                else:
                    limits = [_root_pad_range(panels[interval], 'y') for interval in page]
                    page_limits = min(low for low, _ in limits), max(high for _, high in limits)
                if legend:
                    fig.legend(legend.values(), legend.keys(), loc='upper center',
                               bbox_to_anchor=((left + right) / 2, .956),
                               ncols=len(legend) if contributions else min(4, len(legend)),
                               fontsize=16, frameon=False)
                axes[0, 0].set(xlim=x_range, ylim=page_limits)
                from matplotlib.ticker import MaxNLocator
                axes[0, 0].xaxis.set_major_locator(MaxNLocator(nbins=5, prune='both'))
                axes[0, 0].yaxis.set_major_locator(MaxNLocator(nbins=5, prune='both'))
                fig.supxlabel(x_title, fontsize=29, y=visible_bottom - .055)
                y_title = f'{source} relative systematic uncertainty of {observable} (%)'
                if canvas_view is not None:
                    prefix = f'{source} fits' if canvas_view == 'fits' else source
                    separator = '\n' if rows == 1 else ' — '
                    y_title = prefix + separator + (ratio_title if raw_ratios is not None else canvas_y_title)
                elif rows == 1:
                    y_title = f'{source} relative\nsystematic uncertainty\nof {observable} (%)'
                fig.supylabel(y_title, fontsize=19 if rows == 1 else 28,
                              x=.016 if rows == 1 else .022, y=(top + visible_bottom) / 2)
                header = fig.add_axes((left, .960, right - left, .001), frameon=False)
                header.set_axis_off()
                _draw_cms_header(hep, header, left=cms_header_left, right=cms_header_right,
                                 style=style)
                path = output if page_number == 1 else output.with_stem(
                    f'{output.stem}_page{page_number:02d}')
                save_options = dict(dpi=160, bbox_inches='tight' if rows < 4 else None)
                fig.savefig(path, **save_options)
                paths.append(path)
                if save_png:
                    fig.savefig(path.with_suffix('.png'), **save_options)
                    paths.append(path.with_suffix('.png'))
            finally:
                plt.close(fig)
    return paths


def save_systematic_canvas_grid(source, canvases, dijet_ptave_bins, **options):
    """Reuse the percent TH1 and eta range of each individual source plot."""
    panels = {
        interval: next(item for item in canvas.GetListOfPrimitives()
                       if item.InheritsFrom('TH1') and not item.GetName().startswith('hframe'))
        for interval, canvas in canvases.items()
    }
    first = next(iter(panels.values()))
    return save_systematic_grid(
        source, panels, dijet_ptave_bins, x_title=_root_label(first.GetXaxis().GetTitle()),
        x_range=_root_pad_range(next(iter(canvases.values())), 'x'),
        **options,
    )


def save_systematic_diagnostic_grids(source, results, intervals, *, overlay_key, band_key,
                                     fit_key, raw_ratio_key, ratio_title, stat_scale, output, **options):
    """Reuse the variation canvases, their uncertainty rectangles, and the fitted TF1s."""
    from hist_analysis.python.systematic_methods import systematic_variation_labels

    overlays = {interval: result['canvases'][overlay_key] for interval, result in results.items()}
    bands = {interval: result['canvases'][band_key] if band_key in result['canvases']
             else result[band_key] for interval, result in results.items()}
    fits = {
        interval: systematic_variation_labels(result[fit_key], method=result['systematic_method'])
        for interval, result in results.items()
    }
    raw_ratios = {interval: result[raw_ratio_key] for interval, result in results.items()}
    output = Path(output)
    for view, tag, canvases in (
        ('overlay', 'overlay', overlays),
        ('ratios', 'ratio_with_band', bands),
        ('fits', 'fits_with_band', bands),
    ):
        if view == 'fits' and not any(fits.values()):
            continue
        stem = output.stem.replace('_all_ptave_4x4', f'_{tag}_all_ptave_4x4')
        save_systematic_grid(source, canvases, intervals, x_title='', canvas_view=view, fits=fits,
                             raw_ratios=raw_ratios if view == 'ratios' else None,
                             ratio_title=ratio_title, stat_scale=stat_scale,
                             output=output.with_stem(stem), **options)
