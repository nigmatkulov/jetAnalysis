"""Two-pad closure overlays and ratios to an explicit nominal curve."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import ROOT

from hist_analysis.python.histogram_ops import ratio_to_nominal, style_histograms
from hist_analysis.python.root_style import save_canvas, set_legend_style


def figure_title(generator: str, direction: str, observable: str) -> str:
    return f"{generator} {direction} {observable}"


def draw_closure(histograms: Mapping[str, object], nominal: str, *, title: str,
                 x_title: str, y_title: str = "Arbitrary units",
                 ratio_range: tuple[float, float] = (0.5, 1.5),
                 log_y: bool = False, grid: bool = True,
                 output: str | Path | None = None, save_png: bool = False):
    """Draw an overlay and ratios to an explicitly named nominal histogram."""

    if nominal not in histograms:
        raise KeyError(f"Nominal {nominal!r} is not among {list(histograms)}")
    if not histograms:
        raise ValueError("No histograms supplied")

    style_histograms(histograms)
    suffix = str(abs(hash(title)))
    canvas = ROOT.TCanvas(f"c_{suffix}", title, 800, 800)
    top = ROOT.TPad(f"top_{suffix}", "top", 0.0, 0.30, 1.0, 1.0)
    bottom = ROOT.TPad(f"bottom_{suffix}", "bottom", 0.0, 0.0, 1.0, 0.30)
    top.SetBottomMargin(0.02)
    bottom.SetTopMargin(0.03)
    bottom.SetBottomMargin(0.32)
    top.SetGridx(grid)
    top.SetGridy(grid)
    bottom.SetGridx(grid)
    bottom.SetGridy(grid)
    top.Draw()
    bottom.Draw()

    top.cd()
    top.SetLogy(log_y)
    maximum = max(hist.GetMaximum() for hist in histograms.values())
    legend = ROOT.TLegend(0.7, 0.7, 0.88, 0.88)
    set_legend_style(legend)
    for index, (label, histogram) in enumerate(histograms.items()):
        histogram.SetTitle(title)
        histogram.GetYaxis().SetTitle(y_title)
        histogram.GetXaxis().SetLabelSize(0)
        histogram.SetMaximum(maximum * (20.0 if log_y else 1.3))
        histogram.Draw("E1" if index == 0 else "E1 SAME")
        legend.AddEntry(histogram, label, "p")
    legend.Draw()

    bottom.cd()
    ratios = {}
    for index, (label, histogram) in enumerate(histograms.items()):
        ratio = ratio_to_nominal(histogram, histograms[nominal], name=f"ratio_{index}")
        ratios[label] = ratio
        ratio.SetTitle("")
        ratio.GetYaxis().SetTitle(f"Ratio to {nominal}")
        ratio.GetXaxis().SetTitle(x_title)
        ratio.GetYaxis().SetRangeUser(*ratio_range)
        ratio.GetYaxis().SetNdivisions(505)
        ratio.GetYaxis().SetTitleSize(0.10)
        ratio.GetYaxis().SetLabelSize(0.09)
        ratio.GetYaxis().SetTitleOffset(0.45)
        ratio.GetXaxis().SetTitleSize(0.12)
        ratio.GetXaxis().SetLabelSize(0.10)
        ratio.Draw("E1" if index == 0 else "E1 SAME")
    line = ROOT.TLine(next(iter(ratios.values())).GetXaxis().GetXmin(), 1.0,
                      next(iter(ratios.values())).GetXaxis().GetXmax(), 1.0)
    line.SetLineStyle(2)
    line.Draw()

    canvas.cd()
    canvas.Update()
    save_canvas(canvas, output, save_png=save_png)
    # Keep drawn ROOT objects alive for notebook display.
    canvas._closure_objects = [top, bottom, legend, line, *histograms.values(), *ratios.values()]
    return canvas, ratios
