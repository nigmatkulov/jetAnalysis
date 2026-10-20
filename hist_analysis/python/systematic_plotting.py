"""Shared CMS presentation for systematic notebooks; ROOT owns all numerics."""
from pathlib import Path
import math
import ROOT

from hist_analysis.python.mplhep_plotting import (
    _mplhep, _th1_arrays, _root_label, _root_color, _curve_style,
    _draw_cms_header, DEFAULT_MPLHEP_STYLE, _leaf_pads, _root_pad_range,
)


def validate_stat_uncertainty_scale(value):
    """Require a finite positive real factor; booleans and strings are errors."""
    from numbers import Real
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value) or value <= 0:
        raise ValueError('STAT_UNCERTAINTY_SCALE must be a finite positive number')
    return float(value)


def statistical_plot_copies(histograms, scale):
    """Detached TH1 display copies; call only after extraction and smoothing.

    These workflows store symmetric statistical errors in TH1 Sumw2. Analysis
    histograms (including their errors, entries and statistics) remain untouched.
    Always pass original analysis histograms, never a previous display copy.
    """
    scale = validate_stat_uncertainty_scale(scale)
    copies = {}
    from uuid import uuid4
    for label, histogram in histograms.items():
        if not histogram.InheritsFrom('TH1') or histogram.GetDimension() != 1:
            raise TypeError('Statistical display scaling requires one-dimensional TH1 inputs')
        copy = histogram.Clone(f'{histogram.GetName()}_stat_display_{uuid4().hex}')
        copy.SetDirectory(0)
        for i in range(histogram.GetNcells()):
            copy.SetBinError(i, scale * histogram.GetBinError(i))
        copies[label] = copy
    return copies


def _label(text):
    return _root_label(text.replace('#chi', r'\chi').replace('#pm', r'\pm'))


def _draw_pad(ax, pad, hep, *, cms_header_left, show_header=True,
              x_range=None, y_range=None):
    import numpy as np
    # ROOT may list a histogram more than once after E1 SAME redraws.
    objects = list({int(ROOT.addressof(p)): p
                    for p in pad.GetListOfPrimitives()}.values())
    histograms = [p for p in objects if p.InheritsFrom('TH1')
                  and not p.GetName().startswith('hframe')]
    frame = next((p for p in objects if p.InheritsFrom('TH1')
                  and p.GetName().startswith('hframe')), None)
    reference = frame if frame is not None else histograms[0]
    labels = {}
    for p in objects:
        if p.InheritsFrom('TLegend'):
            for entry in p.GetListOfPrimitives():
                if entry.GetObject():
                    labels[entry.GetObject().GetName()] = entry.GetLabel()
    for graph in objects:
        if not graph.InheritsFrom('TGraphAsymmErrors'):
            continue
        # E2 represents independent rectangular bins, not interpolation between centers.
        for i in range(graph.GetN()):
            x, y = graph.GetPointX(i), graph.GetPointY(i)
            ax.fill_between([x-graph.GetErrorXlow(i), x+graph.GetErrorXhigh(i)],
                            y-graph.GetErrorYlow(i), y+graph.GetErrorYhigh(i),
                            color='lightpink', alpha=0.45, linewidth=0, zorder=1,
                            label='Systematic uncertainty' if i == 0 else '_nolegend_')
    colors = {}
    for index, h in enumerate(histograms):
        label = labels.get(h.GetName(), '_nolegend_')
        _, color, marker, marker_face = _curve_style(label, index)
        colors[h.GetLineColor()] = color
        values, edges, errors = _th1_arrays(h)
        hep.histplot(values, bins=edges, yerr=errors, xerr=True, histtype='errorbar',
                     ax=ax, color=color, marker=marker, markerfacecolor=marker_face,
                     markeredgecolor=color, markeredgewidth=1.5,
                     markersize=DEFAULT_MPLHEP_STYLE.marker_size,
                     capsize=DEFAULT_MPLHEP_STYLE.capsize, zorder=3,
                     label=labels.get(h.GetName(), '_nolegend_'))
    for p in objects:
        if p.InheritsFrom('TF1'):
            x = np.linspace(p.GetXmin(), p.GetXmax(), 500)
            y = [p.Eval(float(v)) for v in x]
            if not all(math.isfinite(v) for v in y):
                raise ValueError(f'Nonfinite fit curve: {p.GetName()}')
            ax.plot(x, y, color=colors.get(p.GetLineColor(), _root_color(p.GetLineColor())),
                    linewidth=1.8, zorder=4, label='_nolegend_')
        elif p.InheritsFrom('TLine'):
            ax.plot([p.GetX1(), p.GetX2()], [p.GetY1(), p.GetY2()],
                    color='gray', linestyle='--', linewidth=1, zorder=2)
    annotations = [p for p in objects if p.InheritsFrom('TLatex')
                   and p.GetTitle() and 'CMS' not in p.GetTitle()
                   and '8.16' not in p.GetTitle()]
    context = [p for p in annotations if '#chi' not in p.GetTitle()]
    summaries = [p for p in annotations if '#chi' in p.GetTitle()]
    for i, p in enumerate(context):
        ax.text(0.04, 0.955 - i * 0.045, _label(p.GetTitle()), transform=ax.transAxes,
                ha='left', va='top', fontsize=12)
    # Beam diagnostic summaries contain several parameters: split into short lines.
    for i, p in enumerate(summaries):
        text = '\n'.join(_label(part) for part in p.GetTitle().split(', '))
        ax.text(0.46, 0.32 - i * 0.15, text, transform=ax.transAxes,
                va='top', fontsize=8,
                color=colors.get(p.GetTextColor(), _root_color(p.GetTextColor())))
    # Keep the optional numerical fit summary, with the same curve colors.
    for p in objects:
        if p.InheritsFrom('TPaveText'):
            for i, line in enumerate(p.GetListOfLines()):
                ax.text(0.46, 0.32-i*0.025, _label(line.GetTitle()),
                        color=colors.get(line.GetTextColor(), _root_color(line.GetTextColor())), fontsize=8,
                        transform=ax.transAxes)
    handles, names = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, names, loc='upper right', fontsize=12, frameon=False)
    ax.set_xlabel(_label(reference.GetXaxis().GetTitle()), fontsize=20)
    ax.set_ylabel(_label(reference.GetYaxis().GetTitle()), fontsize=18)
    ax.set_xlim(*(x_range or _root_pad_range(pad, 'x')))
    ax.set_ylim(*(y_range or _root_pad_range(pad, 'y')))
    ax.set_xscale('log' if pad.GetLogx() else 'linear')
    ax.set_yscale('log' if pad.GetLogy() else 'linear')
    ax.tick_params(labelsize=16)
    ax.grid(bool(pad.GetGridx() or pad.GetGridy()), alpha=0.2)
    if show_header:
        _draw_cms_header(hep, ax, left=cms_header_left,
                         right='pPb (8.16 TeV)', style=DEFAULT_MPLHEP_STYLE)


def save_systematic_mplhep(canvas, output, *, save_png=False,
                           cms_header_left='CMS Simulation', x_range=None, y_range=None):
    """Export one panel, a closure with ratio, or a side-by-side diagnostic."""
    plt, hep = _mplhep()
    hep.style.use('CMS')
    canvas.Update()
    pads = [p for p in _leaf_pads(canvas)
            if any(o.InheritsFrom('TH1') for o in p.GetListOfPrimitives())]
    if not pads:
        raise ValueError(f'No histogram panels in {canvas.GetName()}')
    side_by_side = len(pads) > 1 and len({round(p.GetXlowNDC(), 4) for p in pads}) > 1
    ratio_layout = len(pads) == 2 and not side_by_side
    if ratio_layout:
        pads.sort(key=lambda p: p.GetYlowNDC(), reverse=True)
    rows, columns = (1, len(pads)) if side_by_side else (len(pads), 1)
    fig, axes = plt.subplots(rows, columns,
        figsize=(8 * columns, 8 if len(pads) <= 2 else 6 * rows),
        layout='constrained', squeeze=False,
        gridspec_kw={'height_ratios': (3, 1)} if ratio_layout else None)
    for i, (ax, pad) in enumerate(zip(axes.flat, pads)):
        _draw_pad(ax, pad, hep, cms_header_left=cms_header_left,
                  show_header=not ratio_layout or i == 0,
                  x_range=x_range, y_range=y_range)
        if ratio_layout and i == 1:
            # Reuse top-panel semantic colors for named direction ratios.
            ax.set_ylabel(ax.get_ylabel(), fontsize=15)
            ax.tick_params(labelsize=14)
        if ratio_layout and i == 0:
            ax.set_xlabel('')
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    if save_png:
        fig.savefig(output.with_suffix('.png'), dpi=160)
    from IPython import get_ipython
    if get_ipython() is not None:
        from io import BytesIO
        from IPython.display import Image, display
        preview = BytesIO()
        fig.savefig(preview, format='png', dpi=100)
        display(Image(data=preview.getvalue()))
    plt.close(fig)


def save_systematic_canvas(canvas, output, *, save_png=False, mplhep=True,
                           mplhep_png=False, cms_header_left='CMS Simulation',
                           mplhep_x_range=None, mplhep_y_range=None):
    from hist_analysis.python.root_style import save_canvas
    save_canvas(canvas, output, save_png=save_png)
    if mplhep and output is not None:
        path = Path(output)
        save_systematic_mplhep(canvas, path.parent / 'mplhep' / path.name,
            save_png=mplhep_png, cms_header_left=cms_header_left,
            x_range=mplhep_x_range, y_range=mplhep_y_range)


def draw_systematic_overlay(histograms, *args, stat_uncertainty_scale=1.0,
                            mplhep=True, mplhep_png=False,
                            cms_header_left='CMS Simulation', **kwargs):
    from hist_analysis.python.plotting import draw_overlay
    output = kwargs.pop('output', None)
    scale = validate_stat_uncertainty_scale(stat_uncertainty_scale)
    # Default follows the existing renderer exactly; non-default affects only
    # explicitly opted-in final distributions, never variation/default inputs.
    if scale != 1.0:
        histograms = statistical_plot_copies(histograms, scale)
    canvas = draw_overlay(histograms, *args, output=None, cms_header_left=cms_header_left,
                          cms_header_right='pPb (8.16 TeV)', **kwargs)
    # Retain exactly the plotted points for redraws above systematic bands.
    canvas._systematic_plot_histograms = histograms
    save_systematic_canvas(canvas, output, save_png=kwargs.get('save_png', False),
        mplhep=mplhep, mplhep_png=mplhep_png, cms_header_left=cms_header_left)
    return canvas


def draw_systematic_closure(*args, mplhep=True, mplhep_png=False,
                            cms_header_left='CMS Simulation', **kwargs):
    from hist_analysis.python.plotting import draw_closure
    output = kwargs.pop('output', None)
    canvas, ratios = draw_closure(*args, output=None, cms_header_left=cms_header_left,
                                  cms_header_right='pPb (8.16 TeV)', **kwargs)
    save_systematic_canvas(canvas, output, save_png=kwargs.get('save_png', False),
        mplhep=mplhep, mplhep_png=mplhep_png, cms_header_left=cms_header_left)
    return canvas, ratios


def build_systematic_band(systematic, *, evaluation_range, name):
    """Build pink 1±delta rectangles from fractional bins within acceptance."""
    low, high = evaluation_range
    graph = ROOT.TGraphAsymmErrors()
    graph.SetName(name)
    axis = systematic.GetXaxis()
    for i in range(1, systematic.GetNbinsX() + 1):
        left, right = max(low, axis.GetBinLowEdge(i)), min(high, axis.GetBinUpEdge(i))
        if left >= right:
            continue
        index = graph.GetN()
        graph.SetPoint(index, (left + right) / 2, 1.0)
        half_width, delta = (right - left) / 2, systematic.GetBinContent(i)
        graph.SetPointError(index, half_width, half_width, delta, delta)
    graph.SetFillColorAlpha(ROOT.kPink - 9, 0.45)
    graph.SetLineColor(ROOT.kPink - 9)
    return graph


def draw_systematic_variation_overlay(histograms, *, method, renderer,
                                      stat_uncertainty_scale=1.0, **kwargs):
    """Render finalized direct or folded ratios with scaled display errors."""
    from hist_analysis.python.systematic_methods import (
        validate_systematic_method, systematic_variation_labels,
    )
    validate_systematic_method(method)
    scale = validate_stat_uncertainty_scale(stat_uncertainty_scale)
    if scale != 1.0:
        histograms = statistical_plot_copies(histograms, scale)
    if method == 'direct_ratio':
        return renderer(histograms, **kwargs)
    for key in ('overlay_functions', 'style_indices'):
        if kwargs.get(key) is not None:
            kwargs[key] = systematic_variation_labels(kwargs[key], method=method)
    if kwargs.get('overlay_text') is not None:
        kwargs['overlay_text'] = {
            f'1 + |1 - {label}|': tuple(
                line.replace(f'{label}:', f'Folded {label}:', 1) for line in lines
            )
            for label, lines in kwargs['overlay_text'].items()
        }
    kwargs['y_title'] = f"1 + |1 - {kwargs.get('y_title', 'Variation / default')}|"
    return renderer(systematic_variation_labels(histograms, method=method), **kwargs)
