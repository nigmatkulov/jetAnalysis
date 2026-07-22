"""Two-pad closure overlays and ratios to an explicit nominal curve."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping

import ROOT

from hist_analysis.python.histogram_ops import ratio_to_nominal, style_histograms
from hist_analysis.python.root_style import (
    DEFAULT_PLOT_STYLE, PlotStyle, draw_text_block, save_canvas, set_legend_style,
)


def figure_title(generator: str, direction: str, observable: str) -> str:
    return f"{generator} {direction} {observable}"


def draw_closure(histograms: Mapping[str, object], nominal: str, *, title: str,
                 x_title: str, y_title: str = "Arbitrary units",
                 ratio_range: tuple[float, float] = (0.5, 1.5),
                 log_y: bool = False, grid: bool = True,
                 output: str | Path | None = None, save_png: bool = False,
                 draw_nominal_ratio: bool = True,
                 reference_line_x_range: tuple[float, float] | None = None,
                 headroom: float = 1.3, canvas_name: str | None = None,
                 annotations: Iterable[str] = (),
                 x_range: tuple[float, float] | None = None,
                 style: PlotStyle = DEFAULT_PLOT_STYLE):
    """Draw an overlay and ratios to an explicitly named nominal histogram.

    ``draw_nominal_ratio=False`` suppresses the trivial nominal/self ratio.
    ``x_range`` is applied to both panels, while ``reference_line_x_range`` can
    independently constrain the horizontal ratio reference line. ``annotations``
    are drawn as a shared-style text block in the upper pad. ``headroom`` reserves
    vertical space for legends and annotations, and ``style`` controls the common
    canvas, margin, font, axis, and text geometry.
    """

    if nominal not in histograms:
        raise KeyError(f"Nominal {nominal!r} is not among {list(histograms)}")
    if not histograms:
        raise ValueError("No histograms supplied")

    style_histograms(histograms)
    if headroom <= 1.0:
        raise ValueError("headroom must be greater than 1")
    if x_range is not None and x_range[0] >= x_range[1]:
        raise ValueError("x_range must satisfy low < high")
    annotation_lines = tuple(annotations)
    suffix = str(abs(hash(canvas_name or title)))
    canvas = ROOT.TCanvas(
        f"c_{suffix}", title, style.canvas_width, style.canvas_height,
    )
    top = ROOT.TPad(
        f"top_{suffix}", "top", 0.0, style.ratio_fraction, 1.0, 1.0,
    )
    bottom = ROOT.TPad(
        f"bottom_{suffix}", "bottom", 0.0, 0.0, 1.0, style.ratio_fraction,
    )
    top.SetBottomMargin(style.ratio_top_bottom_margin)
    top.SetLeftMargin(style.ratio_left_margin)
    bottom.SetTopMargin(style.ratio_bottom_top_margin)
    bottom.SetBottomMargin(style.ratio_bottom_margin)
    bottom.SetLeftMargin(style.ratio_left_margin)
    top.SetGridx(grid)
    top.SetGridy(grid)
    bottom.SetGridx(grid)
    bottom.SetGridy(grid)
    top.Draw()
    bottom.Draw()

    top.cd()
    top.SetLogy(log_y)
    maximum = max(
        hist.GetBinContent(bin_index) + hist.GetBinError(bin_index)
        for hist in histograms.values()
        for bin_index in range(1, hist.GetNbinsX() + 1)
    )
    legend = ROOT.TLegend(0.68, 0.72, 0.88, 0.88)
    set_legend_style(legend, style=style)
    for index, (label, histogram) in enumerate(histograms.items()):
        histogram.SetTitle(title)
        histogram.GetYaxis().SetTitle(y_title)
        histogram.GetXaxis().SetTitleFont(style.font)
        histogram.GetXaxis().SetLabelFont(style.font)
        histogram.GetYaxis().SetTitleFont(style.font)
        histogram.GetYaxis().SetLabelFont(style.font)
        histogram.GetYaxis().SetTitleSize(style.ratio_top_axis_title_size)
        histogram.GetYaxis().SetLabelSize(style.ratio_top_axis_label_size)
        histogram.GetYaxis().SetTitleOffset(style.ratio_top_y_title_offset)
        histogram.GetXaxis().SetNdivisions(style.ratio_axis_divisions)
        histogram.GetYaxis().SetNdivisions(style.ratio_axis_divisions)
        histogram.GetXaxis().CenterTitle(True)
        histogram.GetYaxis().CenterTitle(True)
        histogram.GetXaxis().SetLabelSize(0)
        if x_range is not None:
            histogram.GetXaxis().SetRangeUser(*x_range)
        histogram.SetMaximum(maximum * (20.0 if log_y else headroom))
        histogram.Draw("E1" if index == 0 else "E1 SAME")
        legend.AddEntry(histogram, label, "p")
    legend.Draw()
    annotation_objects = draw_text_block(top, annotation_lines, style=style)

    bottom.cd()
    ratios = {}
    for index, (label, histogram) in enumerate(histograms.items()):
        if label == nominal and not draw_nominal_ratio:
            continue
        ratio = ratio_to_nominal(histogram, histograms[nominal], name=f"ratio_{index}")
        ratios[label] = ratio
        ratio.SetTitle("")
        ratio.GetYaxis().SetTitle(f"Ratio to {nominal}")
        ratio.GetXaxis().SetTitle(x_title)
        ratio.GetYaxis().SetRangeUser(*ratio_range)
        ratio.GetXaxis().SetTitleFont(style.font)
        ratio.GetXaxis().SetLabelFont(style.font)
        ratio.GetYaxis().SetTitleFont(style.font)
        ratio.GetYaxis().SetLabelFont(style.font)
        ratio.GetXaxis().SetNdivisions(style.ratio_axis_divisions)
        ratio.GetYaxis().SetNdivisions(style.ratio_axis_divisions)
        ratio.GetXaxis().CenterTitle(True)
        ratio.GetYaxis().CenterTitle(True)
        ratio.GetXaxis().SetTitleSize(style.ratio_axis_title_size)
        ratio.GetYaxis().SetTitleSize(style.ratio_axis_title_size)
        ratio.GetXaxis().SetLabelSize(style.ratio_axis_label_size)
        ratio.GetYaxis().SetLabelSize(style.ratio_axis_label_size)
        ratio.GetXaxis().SetTitleOffset(style.ratio_x_title_offset)
        ratio.GetYaxis().SetTitleOffset(style.ratio_y_title_offset)
        if x_range is not None:
            ratio.GetXaxis().SetRangeUser(*x_range)
        ratio.Draw("E1" if len(ratios) == 1 else "E1 SAME")
    if not ratios:
        raise ValueError("No ratios remain after applying the ratio-selection options")
    ratio_axis = next(iter(ratios.values())).GetXaxis()
    line_x_low, line_x_high = reference_line_x_range or x_range or (
        ratio_axis.GetXmin(), ratio_axis.GetXmax()
    )
    line = ROOT.TLine(line_x_low, 1.0, line_x_high, 1.0)
    line.SetLineStyle(2)
    line.Draw()

    canvas.cd()
    canvas.Update()
    save_canvas(canvas, output, save_png=save_png)
    # Keep drawn ROOT objects alive for notebook display.
    canvas._closure_objects = [
        top, bottom, legend, line, *annotation_objects,
        *histograms.values(), *ratios.values(),
    ]
    return canvas, ratios
