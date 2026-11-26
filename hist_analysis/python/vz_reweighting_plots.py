"""ROOT and native mplhep figures for the v_z reweighting notebook.

A figure is described once as panels of drawable items, then rendered by
both backends.  Native mplhep rendering is used because the generic
ROOT-canvas mirror in ``mplhep_plotting.save_canvas`` does not draw TF1
curves, which are the central objects here.  Displayed contents are the
histograms passed in; nothing is renormalized or rebinned while plotting.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence
from uuid import uuid4

import ROOT

from hist_analysis.python.mplhep_plotting import (
    DEFAULT_MPLHEP_STYLE, _draw_cms_header, _mplhep, _root_label,
)
from hist_analysis.python.root_style import (
    DEFAULT_PLOT_STYLE, draw_cms_header, draw_text_block, save_canvas,
    set_1d_style, set_legend_style,
)


OPEN_MARKERS = {24, 25, 26, 27, 28, 30, 32}
# Filled <-> open counterpart, used for reweighted MC.
MARKER_COMPLEMENT = {20: 24, 21: 25, 22: 26, 23: 32, 33: 27, 34: 28, 29: 30}
MARKER_COMPLEMENT.update({v: k for k, v in list(MARKER_COMPLEMENT.items())})
MPL_MARKERS = {20: "o", 24: "o", 21: "s", 25: "s", 22: "^", 26: "^", 23: "v", 32: "v",
               33: "D", 27: "D", 34: "P", 28: "P", 29: "*", 30: "*"}


@dataclass(frozen=True)
class CurveStyle:
    """ROOT colour and marker; open/filled follows the ROOT marker code."""
    root_color: int
    root_marker: int

    @property
    def filled(self) -> bool:
        return self.root_marker not in OPEN_MARKERS

    @property
    def hex(self) -> str:
        return ROOT.gROOT.GetColor(self.root_color).AsHexString()


# Defaults; notebooks override them through ``styles_from_config``.
DEFAULT_STYLE_CONFIG = {
    "data": (ROOT.kBlack, 20),
    "embedding_Pbgoing": (ROOT.kRed, 21),
    "embedding_pgoing": (ROOT.kBlue, 22),
    "pythia_Pbgoing": (ROOT.kOrange + 7, 25),
    "pythia_pgoing": (ROOT.kMagenta, 26),
    "cp5_Pbgoing": (ROOT.kPink + 1, 33),
}


def styles_from_config(config) -> dict:
    """Map {key: (ROOT colour, ROOT marker)} to CurveStyle; validate markers."""
    styles = {}
    for key, (color, marker) in config.items():
        if marker not in MPL_MARKERS:
            raise ValueError(f"{key}: marker {marker} has no mplhep equivalent; use {sorted(MPL_MARKERS)}")
        if not ROOT.gROOT.GetColor(int(color)):
            raise ValueError(f"{key}: unknown ROOT colour {color}")
        styles[key] = CurveStyle(int(color), int(marker))
    return styles


_DEFAULT_STYLES = styles_from_config(DEFAULT_STYLE_CONFIG)
DATA_STYLE = _DEFAULT_STYLES["data"]
SAMPLE_STYLES = {k: v for k, v in _DEFAULT_STYLES.items() if k != "data"}


@dataclass
class Item:
    kind: str                     # 'hist', 'func', 'hline'
    obj: object = None
    label: str | None = None
    style: CurveStyle | None = None
    line_style: int = 1
    value: float | None = None


@dataclass
class Panel:
    items: list
    y_title: str
    y_range: tuple | None = None
    log_y: bool = False
    text: Sequence[str] = ()
    legend: bool = True
    legend_bounds: tuple = (0.55, 0.62, 0.93, 0.88)
    legend_columns: int = 1
    legend_header: str | None = None
    text_scale: float = 1.0


@dataclass
class Figure:
    panels: list                  # one or two (main + ratio-like lower panel)
    x_title: str
    x_range: tuple
    name: str
    extra: dict = field(default_factory=dict)


def _mpl_label(text: str | None) -> str:
    """Translate only the TLatex sub/superscript tokens to inline mathtext.

    ``v_{z}`` -> ``$v_{z}$``, ``#chi^{2}`` -> ``$\\chi^{2}$``; ordinary words stay
    upright text.
    """
    if not text:
        return ""
    import re
    token = re.compile(r"#?[A-Za-z]+(?:[_^]\{[^}]*\})+")
    text = text.replace("#pm", "\u00b1")
    return token.sub(lambda m: "$" + m.group(0).replace("#", "\\") + "$", text)


def open_style(style: CurveStyle) -> CurveStyle:
    """Complementary (open <-> filled) marker of the same colour."""
    return CurveStyle(style.root_color, MARKER_COMPLEMENT.get(style.root_marker, style.root_marker))


# ------------------------------------------------------------------
# ROOT backend
# ------------------------------------------------------------------
def _root_pad_draw(pad, panel: Panel, *, x_title: str | None, x_range, lower: bool,
                   keep: list):
    pad.cd()
    pad.SetLogy(panel.log_y)
    frame_drawn = False
    n_entries = sum(1 for item in panel.items if item.label and item.kind in {"hist", "func"})
    x1, y1, x2, y2 = panel.legend_bounds
    # Height follows the number of entries; the top edge stays fixed.
    rows = -(-max(n_entries, 1) // panel.legend_columns) + (1 if panel.legend_header else 0)
    y1 = max(y1, y2 - (0.13 if lower else 0.055) * rows)
    legend = ROOT.TLegend(x1, y1, x2, y2)
    set_legend_style(legend)
    legend.SetFillStyle(0)
    legend.SetNColumns(panel.legend_columns)
    if panel.legend_columns > 1:
        legend.SetTextSize(0.024)
    if panel.legend_header:
        legend.SetHeader(panel.legend_header)
    if lower:
        legend.SetTextSize(0.075)
    style = DEFAULT_PLOT_STYLE
    for item in panel.items:
        if item.kind != "hist":
            continue
        histogram = item.obj.Clone(f"{item.obj.GetName()}_draw_{uuid4().hex[:8]}")
        histogram.SetDirectory(0)
        set_1d_style(histogram)
        histogram.SetLineColor(item.style.root_color)
        histogram.SetMarkerColor(item.style.root_color)
        histogram.SetMarkerStyle(item.style.root_marker)
        histogram.SetMarkerSize(1.0)
        histogram.SetTitle("")
        histogram.GetXaxis().SetRangeUser(*x_range)
        histogram.GetYaxis().SetTitle(panel.y_title)
        histogram.GetXaxis().SetTitle(x_title or "")
        if panel.y_range:
            histogram.SetMinimum(panel.y_range[0])
            histogram.SetMaximum(panel.y_range[1])
        if lower:
            for axis in (histogram.GetXaxis(), histogram.GetYaxis()):
                axis.SetTitleSize(style.ratio_axis_title_size)
                axis.SetLabelSize(style.ratio_axis_label_size)
            histogram.GetYaxis().SetTitleOffset(style.ratio_y_title_offset)
            histogram.GetXaxis().SetTitleOffset(style.ratio_x_title_offset)
            histogram.GetYaxis().SetNdivisions(style.ratio_axis_divisions)
        else:
            if x_title is None:  # upper pad of a two-pad figure: no x labels
                histogram.GetXaxis().SetLabelSize(0.0)
                histogram.GetXaxis().SetTitleSize(0.0)
            histogram.GetYaxis().SetTitleSize(style.ratio_top_axis_title_size)
            histogram.GetYaxis().SetLabelSize(style.ratio_top_axis_label_size)
            histogram.GetYaxis().SetTitleOffset(style.ratio_top_y_title_offset)
        histogram.Draw("E1 X0 SAME" if frame_drawn else "E1 X0")
        frame_drawn = True
        keep.append(histogram)
        if item.label:
            legend.AddEntry(histogram, item.label, "p")
    if not frame_drawn:
        raise ValueError("panel needs at least one histogram to define its frame")
    for item in panel.items:
        if item.kind == "func":
            function = item.obj
            function.SetLineColor(item.style.root_color)
            function.SetLineStyle(item.line_style)
            function.SetLineWidth(2)
            function.Draw("L SAME")
            keep.append(function)
            if item.label:
                legend.AddEntry(function, item.label, "l")
        elif item.kind == "hline":
            line = ROOT.TLine(x_range[0], item.value, x_range[1], item.value)
            line.SetLineStyle(2)
            line.SetLineColor(ROOT.kGray + 2)
            line.Draw()
            keep.append(line)
    if panel.legend and legend.GetNRows():
        legend.Draw()
        keep.append(legend)
    if panel.text:
        size = (0.032 if not lower else 0.07) * panel.text_scale
        keep.extend(draw_text_block(pad, panel.text, x=pad.GetLeftMargin() + 0.04,
                                    y=0.86 if not lower else 0.90, text_size=size,
                                    line_spacing=1.45 * size))


def draw_root(figure: Figure, output: Path | None, *, cms_left: str, cms_right: str,
              save_png: bool = False):
    style = DEFAULT_PLOT_STYLE
    suffix = uuid4().hex
    canvas = ROOT.TCanvas(f"c_{figure.name}_{suffix}", figure.name,
                          style.canvas_width, style.canvas_height)
    keep = []
    if len(figure.panels) == 1:
        canvas.SetLeftMargin(style.single_panel_left_margin)
        canvas.SetBottomMargin(style.single_panel_bottom_margin)
        canvas.SetRightMargin(style.right_margin)
        canvas.SetTopMargin(style.top_margin)
        _root_pad_draw(canvas, figure.panels[0], x_title=figure.x_title,
                       x_range=figure.x_range, lower=False, keep=keep)
        header_pad = canvas
    else:
        top = ROOT.TPad(f"top_{suffix}", "", 0.0, style.ratio_fraction, 1.0, 1.0)
        bottom = ROOT.TPad(f"bottom_{suffix}", "", 0.0, 0.0, 1.0, style.ratio_fraction)
        for pad in (top, bottom):
            pad.SetLeftMargin(style.ratio_left_margin)
            pad.SetRightMargin(style.right_margin)
        top.SetTopMargin(style.top_margin)
        top.SetBottomMargin(style.ratio_top_bottom_margin)
        bottom.SetTopMargin(style.ratio_bottom_top_margin)
        bottom.SetBottomMargin(style.ratio_bottom_margin)
        canvas.cd(); top.Draw(); bottom.Draw()
        _root_pad_draw(top, figure.panels[0], x_title=None, x_range=figure.x_range,
                       lower=False, keep=keep)
        _root_pad_draw(bottom, figure.panels[1], x_title=figure.x_title,
                       x_range=figure.x_range, lower=True, keep=keep)
        keep.extend((top, bottom))
        header_pad = top
    keep.extend(draw_cms_header(header_pad, cms_left, cms_right))
    canvas._cms_header_objects = keep[-2:]
    canvas._vz_objects = keep
    canvas.Modified()
    canvas.Update()
    if output is not None:
        save_canvas(canvas, output, save_png=save_png)
    return canvas


# ------------------------------------------------------------------
# mplhep backend
# ------------------------------------------------------------------
def _hist_arrays(histogram, x_range):
    centres, values, errors, widths = [], [], [], []
    axis = histogram.GetXaxis()
    for i in range(1, histogram.GetNbinsX() + 1):
        centre = axis.GetBinCenter(i)
        if not x_range[0] <= centre <= x_range[1]:
            continue
        value, error = histogram.GetBinContent(i), histogram.GetBinError(i)
        if value == 0.0 and error == 0.0:
            continue  # empty/invalid bin: suppress point and error bars
        centres.append(centre); values.append(value); errors.append(error)
        widths.append(0.5 * axis.GetBinWidth(i))
    return centres, values, errors, widths


def _mpl_panel(axis, panel: Panel, x_range, *, lower: bool, style):
    for item in panel.items:
        if item.kind == "hist":
            x, y, ey, ex = _hist_arrays(item.obj, x_range)
            axis.errorbar(x, y, yerr=ey, xerr=ex, linestyle="none",
                          marker=MPL_MARKERS.get(item.style.root_marker, "o"),
                          color=item.style.hex,
                          markerfacecolor=item.style.hex if item.style.filled else "none",
                          markeredgecolor=item.style.hex, capsize=style.capsize,
                          markersize=(style.ratio_marker_size if lower else style.marker_size) * 0.7,
                          label=_mpl_label(item.label) if item.label else "_nolegend_")
        elif item.kind == "func":
            function = item.obj
            low = max(x_range[0], function.GetXmin())
            high = min(x_range[1], function.GetXmax())
            xs = [low + (high - low) * k / 600 for k in range(601)]
            axis.plot(xs, [function.Eval(v) for v in xs], color=item.style.hex,
                      linestyle="-" if item.line_style == 1 else "--", linewidth=2,
                      label=_mpl_label(item.label) if item.label else "_nolegend_")
        elif item.kind == "hline":
            axis.axhline(item.value, color="0.4", linestyle="--", linewidth=1)
    axis.set_xlim(*x_range)
    if panel.y_range:
        axis.set_ylim(*panel.y_range)
    if panel.log_y:
        axis.set_yscale("log")
    axis.set_ylabel(_mpl_label(panel.y_title),
                    fontsize=style.ratio_axis_title_size if lower else style.axis_title_size)
    axis.tick_params(labelsize=style.ratio_tick_label_size if lower else style.tick_label_size)
    if panel.legend and any(item.label for item in panel.items):
        columns = panel.legend_columns
        axis.legend(fontsize=style.legend_text_size * (0.75 if lower else (0.85 if columns == 1 else 0.6)),
                    loc="upper right", frameon=False, ncol=columns,
                    title=_mpl_label(panel.legend_header) if panel.legend_header else None,
                    title_fontsize=style.legend_text_size * 0.7)
    if panel.text:
        axis.text(0.04, 0.95, "\n".join(_mpl_label(t) for t in panel.text),
                  transform=axis.transAxes, va="top", ha="left",
                  fontsize=style.annotation_text_size * (0.75 if lower else 0.85) * panel.text_scale)


def draw_mplhep(figure: Figure, output: Path | None, *, cms_left: str, cms_right: str,
                save_png: bool = False):
    plt, hep = _mplhep()
    hep.style.use("CMS")
    style = DEFAULT_MPLHEP_STYLE
    if len(figure.panels) == 1:
        fig, axes = plt.subplots(1, 1, figsize=(style.figure_size, style.figure_size))
        axes = [axes]
    else:
        fig, axes = plt.subplots(2, 1, figsize=(style.figure_size, style.figure_size),
                                 sharex=True, gridspec_kw=dict(height_ratios=(3, 1), hspace=0.05))
    for index, (axis, panel) in enumerate(zip(axes, figure.panels)):
        _mpl_panel(axis, panel, figure.x_range, lower=index > 0, style=style)
    axes[-1].set_xlabel(_mpl_label(figure.x_title), fontsize=style.axis_title_size)
    _draw_cms_header(hep, axes[0], left=cms_left, right=cms_right, style=style)
    if output is not None:
        path = Path(output)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=200, bbox_inches="tight")
        if save_png:
            fig.savefig(path.with_suffix(".png"), dpi=200, bbox_inches="tight")
    return fig
