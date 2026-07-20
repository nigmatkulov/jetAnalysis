"""Event-level overlays and dijet overweight-protection diagnostics."""

from __future__ import annotations

from array import array
from pathlib import Path
from typing import Mapping

import ROOT

from hist_analysis.python.histogram_ops import normalize_histogram, style_histograms
from hist_analysis.python.root_style import (
    save_canvas, set_1d_style, set_2d_style, set_legend_style, set_pad_style,
    set_palette_style,
)


def _format_percentage(fraction: float) -> str:
    """Format a fractional probability as a percentage without fixed rounding."""

    return f"{100.0 * fraction:.6g}"


def _draw_fit_annotation(fits: Mapping[str, object]):
    """Draw the fit expression and all fitted parameters in NDC coordinates."""

    text_objects = []
    equation = ROOT.TLatex()
    equation.SetNDC(True)
    equation.SetTextFont(42)
    equation.SetTextSize(0.026)
    text_objects.append(equation.DrawLatex(
        0.43, 0.86,
        "f(x) = p_{0} + p_{1}e^{-p_{2}x} + p_{3}e^{-p_{4}x}",
    ))

    y_position = 0.815
    for label, fit in fits.items():
        values = [fit.GetParameter(index) for index in range(fit.GetNpar())]
        split = (len(values) + 1) // 2
        for start, stop in ((0, split), (split, len(values))):
            if start == stop:
                continue
            parameters = ", ".join(
                f"p_{{{index}}}={values[index]:.4g}" for index in range(start, stop)
            )
            prefix = f"{label}: " if start == 0 else ""
            text = ROOT.TLatex()
            text.SetNDC(True)
            text.SetTextFont(42)
            text.SetTextSize(0.023)
            text.SetTextColor(fit.GetLineColor())
            text_objects.append(text.DrawLatex(0.43, y_position, prefix + parameters))
            y_position -= 0.035
        y_position -= 0.012
    return text_objects


def upper_tail_threshold(histogram, cut_fraction: float = 0.0001, *, name: str | None = None):
    """Return the per-X-bin upper-tail threshold of a TH2 Y distribution.

    This mirrors ``plotOverweightProtection`` in ``macro/plotMcClosures.C``:
    underflow and overflow are excluded, and the threshold is the center of the
    first Y bin for which the integral accumulated from high Y reaches
    ``cut_fraction`` of the in-range integral.
    """

    if not histogram.InheritsFrom("TH2"):
        raise TypeError(f"Expected TH2, got {histogram.ClassName()}")
    if not 0.0 < cut_fraction <= 1.0:
        raise ValueError("cut_fraction must be in (0, 1]")

    x_axis = histogram.GetXaxis()
    edges = [x_axis.GetBinLowEdge(index) for index in range(1, x_axis.GetNbins() + 2)]

    threshold = ROOT.TH1D(
        name or f"{histogram.GetName()}_upperTailThreshold",
        histogram.GetTitle(), x_axis.GetNbins(), array("d", edges),
    )
    threshold.SetDirectory(0)
    threshold.SetStats(False)

    y_axis = histogram.GetYaxis()
    for x_bin in range(1, x_axis.GetNbins() + 1):
        total = sum(
            histogram.GetBinContent(x_bin, y_bin)
            for y_bin in range(1, y_axis.GetNbins() + 1)
        )
        if total <= 0.0:
            continue

        target = cut_fraction * total
        accumulated = 0.0
        for y_bin in range(y_axis.GetNbins(), 0, -1):
            accumulated += histogram.GetBinContent(x_bin, y_bin)
            if accumulated >= target:
                threshold.SetBinContent(x_bin, y_axis.GetBinCenter(y_bin))
                threshold.SetBinError(x_bin, 0.05)
                break

    return threshold


def draw_event_histograms(histograms: Mapping[str, object], *, title: str,
                          x_title: str, normalization: str = "integral",
                          log_y: bool = False, grid: bool = True,
                          output: str | Path | None = None, save_png: bool = False):
    """Draw event-level TH1 histograms with the requested normalization."""

    if not histograms:
        raise ValueError("No histograms supplied")
    plotted = {
        label: normalize_histogram(histogram, normalization)
        for label, histogram in histograms.items()
    }
    style_histograms(plotted)

    suffix = str(abs(hash((title, tuple(plotted)))))
    canvas = ROOT.TCanvas(f"c_event_{suffix}", title, 800, 800)
    set_pad_style(canvas, grid_x=grid, grid_y=grid)
    canvas.SetLogy(log_y)
    legend = ROOT.TLegend(0.7, 0.7, 0.88, 0.88)
    set_legend_style(legend)
    maximum = max(histogram.GetMaximum() for histogram in plotted.values())
    for index, (label, histogram) in enumerate(plotted.items()):
        histogram.SetTitle(title)
        histogram.GetXaxis().SetTitle(x_title)
        histogram.GetYaxis().SetTitle(
            "Unit-normalized events" if normalization == "integral" else "Events"
        )
        histogram.SetMaximum(maximum * (20.0 if log_y else 1.3))
        histogram.Draw("E1" if index == 0 else "E1 SAME")
        legend.AddEntry(histogram, label, "p")
    legend.Draw()
    canvas.Update()

    save_canvas(canvas, output, save_png=save_png)
    canvas._event_objects = [legend, *plotted.values()]
    return canvas, plotted


def draw_overweight_protection(gen_histogram, reco_histogram, *, cut_fraction: float = 0.0001,
                               grid: bool = True, output: str | Path | None = None,
                               save_png: bool = False):
    """Plot gen/reco upper-tail thresholds and two-exponential empirical fits.

    Both fits use ``p0 + p1*exp(-p2*x) + p3*exp(-p4*x)`` over 15–950 GeV.
    The threshold-bin definition follows ``plotOverweightProtection`` in
    ``macro/plotMcClosures.C``.
    """

    thresholds = {
        "Gen": upper_tail_threshold(gen_histogram, cut_fraction, name="hGenUpperTailThreshold"),
        "Reco": upper_tail_threshold(reco_histogram, cut_fraction, name="hRecoUpperTailThreshold"),
    }
    set_1d_style(thresholds["Gen"], 1)
    set_1d_style(thresholds["Reco"], 0)

    fits = {
        "Gen": ROOT.TF1("fGenUpperTailThreshold", "[0]+[1]*exp(-[2]*x)+[3]*exp(-[4]*x)", 15.0, 950.0),
        "Reco": ROOT.TF1("fRecoUpperTailThreshold", "[0]+[1]*exp(-[2]*x)+[3]*exp(-[4]*x)", 15.0, 950.0),
        # "Gen": ROOT.TF1("fGenUpperTailThreshold", "[0]+[1]*x+[2]*exp(-[3]*x)", 15.0, 950.0),
        # "Reco": ROOT.TF1("fRecoUpperTailThreshold", "[0]+[1]*x+[2]*exp(-[3]*x)", 15.0, 950.0),
    }
    for label, fit in fits.items():
        fit.SetParameters(1.0, 1.0, 0.01, 1.0, 0.01)
        # fit.SetParameters(1.0, 1.0, 1.0, 0.1)
        fit.SetLineColor(thresholds[label].GetLineColor())
        fit.SetLineWidth(3)
        fit.SetNpx(1000)
        thresholds[label].Fit(fit, "MRE0")

    canvas = ROOT.TCanvas("c_overweight_threshold", "Overweight protection", 800, 800)
    set_pad_style(canvas, grid_x=grid, grid_y=grid)
    maximum = max(histogram.GetMaximum() for histogram in thresholds.values()) or 1.0
    first = thresholds["Gen"]
    first.SetTitle("Overweight protection upper-tail threshold")
    first.GetXaxis().SetTitle("#hat{p}_{T} (GeV)")
    first.GetYaxis().SetTitle(
        f"Tail threshold ({_format_percentage(cut_fraction)}%)"
    )
    first.GetYaxis().SetTitleSize(0.045)
    first.SetMinimum(0.75)
    first.SetMaximum(1.2 * maximum)
    first.Draw("P")
    thresholds["Reco"].Draw("P SAME")
    fits["Gen"].Draw("SAME")
    fits["Reco"].Draw("SAME")
    legend = ROOT.TLegend(0.18, 0.72, 0.42, 0.84)
    set_legend_style(legend)
    for label in ("Gen", "Reco"):
        legend.AddEntry(thresholds[label], label, "p")
    legend.Draw()
    fit_text = _draw_fit_annotation(fits)
    canvas.Update()

    save_canvas(canvas, output, save_png=save_png)
    canvas._overweight_objects = [legend, *fit_text, *thresholds.values(), *fits.values()]
    return canvas, thresholds, fits


def draw_overweight_map(histogram, *, title: str, fit=None, grid: bool = True,
                        log_z: bool = True, output: str | Path | None = None,
                        save_png: bool = False, palette_width: float = 0.02):
    """Draw an overweight TH2, optionally overlaying its fitted threshold."""

    suffix = str(abs(hash((histogram.GetName(), title, bool(fit)))))
    canvas = ROOT.TCanvas(f"c_overweight_map_{suffix}", title, 800, 800)
    set_pad_style(canvas, grid_x=grid, grid_y=grid)
    canvas.SetRightMargin(0.12)
    canvas.SetLogz(log_z)
    histogram.SetTitle(title)
    histogram.GetXaxis().SetTitle("#hat{p}_{T} (GeV)")
    histogram.GetYaxis().SetTitle("p_{T}^{ave}/#hat{p}_{T}")
    set_2d_style(histogram)
    histogram.Draw("COLZ")
    canvas.Update()
    palette = set_palette_style(histogram, width=palette_width)

    objects = [histogram, palette]
    if fit is not None:
        fit.Draw("SAME")
        objects.append(fit)
    canvas.Update()
    save_canvas(canvas, output, save_png=save_png)
    canvas._overweight_map_objects = objects
    return canvas
