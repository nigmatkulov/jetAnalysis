"""MPLHEP presentation of the JES/JER notebook's ROOT canvases.

Read drawable objects only: projection, fit, and histogram values stay in ROOT.
"""
from pathlib import Path
import re

from hist_analysis.python.mplhep_plotting import (
    _mplhep, _th1_arrays, _th2_arrays, _root_color, _curve_style, _draw_cms_header,
    DEFAULT_MPLHEP_STYLE,
)
from hist_analysis.python.plot_labels import CMS_SIMULATION_HEADER


def physics_label(text):
    """Translate the notebook's small TLatex vocabulary into standard mathtext."""
    tokens = {'eta': r'\eta', 'sigma': r'\sigma', 'mu': r'\mu'}
    text = re.sub(r'#(eta|sigma|mu)(?![A-Za-z])', lambda m: tokens[m[1]], text)
    text = re.sub(r'\{(lab|ref|truth|reco|smeared|gen|T)\}', r'{\\mathrm{\1}}', text)
    if '\\' in text or '_{' in text or '^{' in text:
        # Preserve prose in mixed labels, including the eta-dependent variation.
        if text.startswith('gen (nominal, '):
            return r'gen (nominal, $\eta$-dep)'
        return '$' + text.replace(' ', r'\ ') + '$'
    return text


def save_jes_jer_canvas(canvas, *, output, save_png=False):
    """Export all selections, legends, and TF1 curves with one common layout."""
    import numpy as np

    plt, hep = _mplhep()
    hep.style.use('CMS')
    primitives = list(canvas.GetListOfPrimitives())
    histograms = [p for p in primitives if p.InheritsFrom('TH1')]
    reference = histograms[0]
    labels = {}
    for primitive in primitives:
        if primitive.InheritsFrom('TLegend'):
            for entry in primitive.GetListOfPrimitives():
                if entry.GetObject():
                    labels[entry.GetObject().GetName()] = physics_label(entry.GetLabel())
    annotations = [physics_label(p.GetTitle()) for p in primitives
                   if p.InheritsFrom('TLatex') and 'CMS' not in p.GetTitle()
                   and '8.16' not in p.GetTitle()]
    fig, axis = plt.subplots(figsize=(8, 8), layout='constrained')
    # Keep common context in the upper left and curve selections in the upper right.
    # Split the mixed reco/gen response definition to keep the left block compact.
    annotation_text = '\n'.join(annotations).replace(r';\ ', '$\n$').replace('; ', '\n')
    axis.text(0.04, 0.955, annotation_text, transform=axis.transAxes,
              ha='left', va='top', fontsize=16, linespacing=1.15,
              bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=3)
              if reference.InheritsFrom('TH2') else None)
    if reference.InheritsFrom('TH2'):
        # The axes-divider colorbar needs explicit space with this square layout.
        fig.set_layout_engine(None)
        fig.subplots_adjust(left=0.14, right=0.84, bottom=0.12, top=0.93)
        values, x_edges, y_edges = _th2_arrays(reference)
        hep.hist2dplot(values, xbins=x_edges, ybins=y_edges, ax=axis,
                       cbar=True, cbarsize='4%', cbarpad=0.04)
    else:
        fit_colors = {}
        for index, histogram in enumerate(histograms):
            values, edges, errors = _th1_arrays(histogram)
            label = labels.get(histogram.GetName(), histogram.GetName())
            _, color, marker, marker_face = _curve_style(label, index)
            fit_colors[histogram.GetMarkerColor()] = color
            hep.histplot(values, bins=edges, yerr=errors, xerr=True,
                         histtype='errorbar', ax=axis,
                         color=color, marker=marker, markerfacecolor=marker_face,
                         markeredgecolor=color, markeredgewidth=1.5,
                         markersize=DEFAULT_MPLHEP_STYLE.marker_size,
                         capsize=DEFAULT_MPLHEP_STYLE.capsize,
                         label=labels.get(histogram.GetName()))
        for fit in primitives:
            if fit.InheritsFrom('TF1'):
                x = np.linspace(fit.GetXmin(), fit.GetXmax(), 400)
                axis.plot(x, [fit.Eval(float(v)) for v in x],
                          color=fit_colors.get(fit.GetLineColor(), _root_color(fit.GetLineColor())),
                          linewidth=1.5)
        if labels:
            axis.legend(loc='upper right',
                        fontsize=11 if len(labels) > 7 else 17, frameon=False,
                        labelspacing=0.15 if len(labels) > 7 else 0.5)
    axis.set_xlabel(physics_label(reference.GetXaxis().GetTitle()), fontsize=22)
    axis.set_ylabel(physics_label(reference.GetYaxis().GetTitle()), fontsize=22)
    axis.set_xlim(canvas.GetUxmin(), canvas.GetUxmax())
    axis.set_ylim(canvas.GetUymin(), canvas.GetUymax())
    axis.tick_params(labelsize=17)
    axis.grid(False)
    fig.canvas.draw()
    _draw_cms_header(hep, axis, left=CMS_SIMULATION_HEADER,
                     right='pPb (8.16 TeV)', style=DEFAULT_MPLHEP_STYLE)
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    if save_png:
        fig.savefig(path.with_suffix('.png'), dpi=150)
    # Display the same figure that was saved, including in a batch ROOT notebook.
    from IPython import get_ipython
    if get_ipython() is not None:
        from io import BytesIO
        from IPython.display import Image, display
        preview = BytesIO()
        fig.savefig(preview, format='png', dpi=100)
        display(Image(data=preview.getvalue()))
    plt.close(fig)
