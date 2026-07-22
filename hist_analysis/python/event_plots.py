"""Event-level overlays and jet overweight-protection diagnostics."""

from __future__ import annotations

from array import array
import math
from pathlib import Path
from typing import Mapping

import ROOT

from hist_analysis.python.histogram_ops import normalize_histogram, style_histograms
from hist_analysis.python.root_style import (
    DEFAULT_PLOT_STYLE, draw_text_block, save_canvas, set_1d_style, set_2d_style,
    set_legend_style, set_pad_style, set_palette_style,
)


def _format_percentage(fraction: float) -> str:
    """Format a fractional probability as a percentage without fixed rounding."""

    return f"{100.0 * fraction:.6g}"


def _draw_plot_context(generator_label: str | None, orientation_label: str | None):
    """Draw generator and beam-orientation labels in the upper-left corner."""

    return draw_text_block(
        ROOT.gPad, (generator_label, orientation_label), x=0.18, y=0.86,
    )


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


def print_fit_summary(fits: Mapping[str, object]) -> None:
    """Print every parameter with its error, followed by chi2/NDF."""

    for level, fit in fits.items():
        print(f"{level} fit:")
        for index in range(fit.GetNpar()):
            print(
                f"  p{index} = {fit.GetParameter(index):.8g} "
                f"+/- {fit.GetParError(index):.3g}"
            )
        chi2 = fit.GetChisquare()
        ndf = fit.GetNDF()
        chi2_ndf = chi2 / ndf if ndf > 0 else float("nan")
        print(f"  chi2 = {chi2:.6g}, NDF = {ndf}, chi2/NDF = {chi2_ndf:.6g}")


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


def rejected_fraction_by_x(total_histogram, passed_histogram, *, name: str | None = None):
    """Return ``(all - pass) / all`` for every TH2 X bin.

    Y underflow and overflow are excluded. Errors account for pass being a
    subset of all: the failed-event sum of weights squared is obtained from
    ``sumw2(all) - sumw2(pass)`` and propagated with its covariance with all.
    """

    for histogram in (total_histogram, passed_histogram):
        if not histogram.InheritsFrom("TH2"):
            raise TypeError(f"Expected TH2, got {histogram.ClassName()}")
    if (total_histogram.GetNbinsX() != passed_histogram.GetNbinsX()
            or total_histogram.GetNbinsY() != passed_histogram.GetNbinsY()):
        raise ValueError("Total and pass histograms have incompatible dimensions")

    total_x = total_histogram.GetXaxis()
    pass_x = passed_histogram.GetXaxis()
    for x_bin in range(1, total_x.GetNbins() + 2):
        if not math.isclose(
            total_x.GetBinLowEdge(x_bin), pass_x.GetBinLowEdge(x_bin),
            rel_tol=0.0, abs_tol=1e-9,
        ):
            raise ValueError("Total and pass histograms have incompatible X binning")
    total_y = total_histogram.GetYaxis()
    pass_y = passed_histogram.GetYaxis()
    for y_bin in range(1, total_y.GetNbins() + 2):
        if not math.isclose(
            total_y.GetBinLowEdge(y_bin), pass_y.GetBinLowEdge(y_bin),
            rel_tol=0.0, abs_tol=1e-9,
        ):
            raise ValueError("Total and pass histograms have incompatible Y binning")

    edges = [total_x.GetBinLowEdge(index) for index in range(1, total_x.GetNbins() + 2)]
    result = ROOT.TH1D(
        name or f"{total_histogram.GetName()}_rejectedFraction",
        "", total_x.GetNbins(), array("d", edges),
    )
    result.SetDirectory(0)
    result.Sumw2()

    n_y_bins = total_histogram.GetNbinsY()
    for x_bin in range(1, total_x.GetNbins() + 1):
        total = sum(
            total_histogram.GetBinContent(x_bin, y_bin)
            for y_bin in range(1, n_y_bins + 1)
        )
        passed = sum(
            passed_histogram.GetBinContent(x_bin, y_bin)
            for y_bin in range(1, n_y_bins + 1)
        )
        if total <= 0.0:
            continue
        tolerance = 1e-10 * max(abs(total), abs(passed), 1e-300)
        if passed > total + tolerance:
            raise ValueError(
                f"Pass integral exceeds total in X bin {x_bin}: {passed} > {total}"
            )

        failed = max(0.0, total - passed)
        total_variance = sum(
            total_histogram.GetBinError(x_bin, y_bin) ** 2
            for y_bin in range(1, n_y_bins + 1)
        )
        pass_variance = sum(
            passed_histogram.GetBinError(x_bin, y_bin) ** 2
            for y_bin in range(1, n_y_bins + 1)
        )
        failed_variance = max(0.0, total_variance - pass_variance)
        fraction = failed / total
        variance = (
            failed_variance / total**2
            + failed**2 * total_variance / total**4
            - 2.0 * failed * failed_variance / total**3
        )
        result.SetBinContent(x_bin, fraction)
        result.SetBinError(x_bin, math.sqrt(max(0.0, variance)))

    return result


def draw_event_histograms(histograms: Mapping[str, object], *, title: str,
                          x_title: str, normalization: str = "integral",
                          log_y: bool = False, grid: bool = True,
                          output: str | Path | None = None, save_png: bool = False,
                          y_title: str | None = None,
                          y_range: tuple[float, float] | None = None,
                          generator_label: str | None = None,
                          orientation_label: str | None = None):
    """Draw event-level TH1 histograms with the requested normalization."""

    if not histograms:
        raise ValueError("No histograms supplied")
    plotted = {
        label: normalize_histogram(histogram, normalization)
        for label, histogram in histograms.items()
    }
    style_histograms(plotted)

    suffix = str(abs(hash((title, tuple(plotted)))))
    canvas = ROOT.TCanvas(
        f"c_event_{suffix}", title,
        DEFAULT_PLOT_STYLE.canvas_width, DEFAULT_PLOT_STYLE.canvas_height,
    )
    set_pad_style(canvas, grid_x=grid, grid_y=grid)
    canvas.SetLogy(log_y)
    legend = ROOT.TLegend(0.7, 0.7, 0.88, 0.88)
    set_legend_style(legend)
    maximum = max(histogram.GetMaximum() for histogram in plotted.values())
    for index, (label, histogram) in enumerate(plotted.items()):
        histogram.SetTitle(title)
        histogram.GetXaxis().SetTitle(x_title)
        histogram.GetYaxis().SetTitle(
            y_title or ("Unit-normalized events" if normalization == "integral" else "Events")
        )
        histogram.SetMaximum(maximum * (20.0 if log_y else 1.3))
        if y_range is not None:
            histogram.GetYaxis().SetRangeUser(*y_range)
        histogram.Draw("E1" if index == 0 else "E1 SAME")
        legend.AddEntry(histogram, label, "p")
    legend.Draw()
    context_text = _draw_plot_context(generator_label, orientation_label)
    canvas.Update()

    save_canvas(canvas, output, save_png=save_png)
    canvas._event_objects = [legend, *context_text, *plotted.values()]
    return canvas, plotted


def draw_overweight_protection(gen_histogram, reco_histogram, *, cut_fraction: float = 0.0001,
                               grid: bool = True, output: str | Path | None = None,
                               save_png: bool = False, observable_label: str = "Dijet",
                               name_prefix: str = "Dijet",
                               generator_label: str | None = None,
                               orientation_label: str | None = None):
    """Plot gen/reco upper-tail thresholds and two-exponential empirical fits.

    Both fits use ``p0 + p1*exp(-p2*x) + p3*exp(-p4*x)`` over 15–950 GeV.
    The threshold-bin definition follows ``plotOverweightProtection`` in
    ``macro/plotMcClosures.C``.
    """

    thresholds = {
        "Gen": upper_tail_threshold(
            gen_histogram, cut_fraction, name=f"hGen{name_prefix}UpperTailThreshold",
        ),
        "Reco": upper_tail_threshold(
            reco_histogram, cut_fraction, name=f"hReco{name_prefix}UpperTailThreshold",
        ),
    }
    set_1d_style(thresholds["Gen"], 1)
    set_1d_style(thresholds["Reco"], 0)

    fits = {
        "Gen": ROOT.TF1(f"fGen{name_prefix}UpperTailThreshold", "[0]+[1]*exp(-[2]*x)+[3]*exp(-[4]*x)", 15.0, 950.0),
        "Reco": ROOT.TF1(f"fReco{name_prefix}UpperTailThreshold", "[0]+[1]*exp(-[2]*x)+[3]*exp(-[4]*x)", 15.0, 950.0),
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

    canvas = ROOT.TCanvas(
        f"c_{name_prefix}_overweight_threshold",
        f"{observable_label} overweight protection",
        DEFAULT_PLOT_STYLE.canvas_width, DEFAULT_PLOT_STYLE.canvas_height,
    )
    set_pad_style(canvas, grid_x=grid, grid_y=grid)
    maximum = max(histogram.GetMaximum() for histogram in thresholds.values()) or 1.0
    first = thresholds["Gen"]
    first.SetTitle(f"{observable_label} overweight protection upper-tail threshold")
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
    legend = ROOT.TLegend(0.18, 0.64, 0.42, 0.74)
    set_legend_style(legend)
    for label in ("Gen", "Reco"):
        legend.AddEntry(thresholds[label], label, "p")
    legend.Draw()
    context_text = _draw_plot_context(generator_label, orientation_label)
    fit_text = _draw_fit_annotation(fits)
    canvas.Update()

    save_canvas(canvas, output, save_png=save_png)
    canvas._overweight_objects = [
        legend, *context_text, *fit_text, *thresholds.values(), *fits.values(),
    ]
    return canvas, thresholds, fits


def draw_overweight_map(histogram, *, title: str, fit=None, grid: bool = True,
                        log_z: bool = True, output: str | Path | None = None,
                        save_png: bool = False, palette_width: float | None = None,
                        y_title: str = "p_{T}^{ave}/#hat{p}_{T}",
                        generator_label: str | None = None,
                        orientation_label: str | None = None):
    """Draw an overweight TH2, optionally overlaying its fitted threshold."""

    suffix = str(abs(hash((histogram.GetName(), title, bool(fit)))))
    canvas = ROOT.TCanvas(
        f"c_overweight_map_{suffix}", title,
        DEFAULT_PLOT_STYLE.canvas_width, DEFAULT_PLOT_STYLE.canvas_height,
    )
    set_pad_style(canvas, grid_x=grid, grid_y=grid)
    canvas.SetRightMargin(DEFAULT_PLOT_STYLE.palette_right_margin)
    canvas.SetLogz(log_z)
    histogram.SetTitle(title)
    histogram.GetXaxis().SetTitle("#hat{p}_{T} (GeV)")
    histogram.GetYaxis().SetTitle(y_title)
    set_2d_style(histogram)
    histogram.Draw("COLZ")
    canvas.Update()
    palette = set_palette_style(histogram, width=palette_width)

    context_text = _draw_plot_context(generator_label, orientation_label)
    objects = [histogram, palette, *context_text]
    if fit is not None:
        fit.Draw("SAME")
        objects.append(fit)
    canvas.Update()
    save_canvas(canvas, output, save_png=save_png)
    canvas._overweight_map_objects = objects
    return canvas
