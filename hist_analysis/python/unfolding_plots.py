"""Plots shared by flattened dijet unfolding notebooks."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping, Sequence

import ROOT

from hist_analysis.python.root_style import (
    DEFAULT_PLOT_STYLE, PlotStyle, draw_text_block, save_canvas, set_2d_style,
    set_legend_style, set_pad_style, set_unfolding_1d_style,
)
from hist_analysis.python.unfolding import (
    FlattenedBinning, extract_eta_block,
)


def draw_unfolding_classification(
    histogram,
    *,
    annotations: Iterable[str] = (),
    output: str | Path | None = None,
    save_png: bool = False,
    grid: bool = False,
    canvas_name: str = "canvas_unfolding_classification",
    style: PlotStyle = DEFAULT_PLOT_STYLE,
):
    result = histogram.Clone(f"{histogram.GetName()}Display")
    result.SetDirectory(0)
    integral = result.Integral()
    if integral <= 0.0:
        raise ValueError(f"Cannot normalize classification with integral {integral}")
    result.Scale(1.0 / integral)
    positive = [
        result.GetBinContent(index) for index in range(1, result.GetNbinsX() + 1)
        if result.GetBinContent(index) > 0.0
    ]
    if not positive:
        raise ValueError("Classification has no positive bins")
    set_unfolding_1d_style(result, "reco", style=style)
    result.SetTitle(";Pair category;Normalized weighted events")
    result.SetMinimum(0.5 * min(positive))
    result.SetMaximum(10.0 * max(positive))
    result.GetXaxis().LabelsOption("v")
    result.GetXaxis().SetLabelSize(0.035)
    canvas = ROOT.TCanvas(canvas_name, "Unfolding pair classification", 1000, 750)
    set_pad_style(canvas, grid_x=grid, grid_y=grid, style=style)
    canvas.SetBottomMargin(0.30)
    canvas.SetLogy(True)
    result.Draw("E1")
    labels = draw_text_block(canvas, annotations, style=style)
    canvas.Modified()
    canvas.Update()
    save_canvas(canvas, output, save_png=save_png)
    canvas._unfolding_objects = [result, *labels]
    return canvas, result


def draw_projection_response(
    spectra: Mapping[str, tuple[str, object]],
    response,
    *,
    eta_range: tuple[float, float],
    response_titles: tuple[str, str],
    annotations: Iterable[str] = (),
    plot_miss_and_fakes: bool = True,
    output: str | Path | None = None,
    save_png: bool = False,
    grid: bool = True,
    canvas_name: str = "canvas_projection_response",
    style: PlotStyle = DEFAULT_PLOT_STYLE,
):
    displayed = {
        role: value for role, value in spectra.items()
        if plot_miss_and_fakes or role not in {"miss", "fake"}
    }
    if not displayed:
        raise ValueError("No spectra supplied")
    canvas = ROOT.TCanvas(
        canvas_name, "Projection and response",
        style.side_by_side_canvas_width, style.side_by_side_canvas_height,
    )
    canvas.Divide(2, 1)
    canvas.cd(1)
    set_pad_style(ROOT.gPad, grid_x=grid, grid_y=grid, style=style)
    maximum = max(hist.GetMaximum() for _, hist in displayed.values())
    legend = ROOT.TLegend(0.62, 0.67, 0.88, 0.88)
    set_legend_style(legend, style=style)
    for index, (role, (label, histogram)) in enumerate(displayed.items()):
        set_unfolding_1d_style(histogram, role, style=style)
        histogram.SetTitle(";#eta_{CM};Entries")
        histogram.GetXaxis().SetRangeUser(*eta_range)
        if index == 0:
            histogram.SetMaximum(1.25 * maximum if maximum > 0.0 else 1.0)
        histogram.Draw("E1" if index == 0 else "E1 SAME")
        legend.AddEntry(histogram, label, "p")
    legend.Draw()
    left_labels = draw_text_block(ROOT.gPad, annotations, style=style)
    canvas.cd(2)
    set_pad_style(ROOT.gPad, grid_x=grid, grid_y=grid, style=style)
    ROOT.gPad.SetRightMargin(style.palette_right_margin)
    set_2d_style(response, style=style)
    response.SetTitle(f";{response_titles[0]};{response_titles[1]}")
    response.GetXaxis().SetRangeUser(*eta_range)
    response.GetYaxis().SetRangeUser(*eta_range)
    response.Draw("COLZ")
    right_labels = draw_text_block(ROOT.gPad, annotations, style=style)
    canvas.Modified()
    canvas.Update()
    save_canvas(canvas, output, save_png=save_png)
    canvas._unfolding_objects = [
        legend, *left_labels, *right_labels, response,
        *(hist for _, hist in displayed.values()),
    ]
    return canvas


def draw_flattened_response(
    spectra: Mapping[str, tuple[str, object]],
    response,
    *,
    annotations: Iterable[str] = (),
    plot_miss_and_fakes: bool = True,
    output: str | Path | None = None,
    save_png: bool = False,
    grid: bool = True,
    canvas_name: str = "canvas_flattened_response",
    style: PlotStyle = DEFAULT_PLOT_STYLE,
):
    displayed = {
        role: value for role, value in spectra.items()
        if plot_miss_and_fakes or role not in {"miss", "fake"}
    }
    if not displayed:
        raise ValueError("No flattened spectra supplied")
    canvas = ROOT.TCanvas(
        canvas_name, "Flattened histograms and response",
        style.side_by_side_canvas_width, style.side_by_side_canvas_height,
    )
    canvas.Divide(2, 1)
    canvas.cd(1)
    set_pad_style(ROOT.gPad, grid_x=grid, grid_y=grid, style=style)
    maximum = max(hist.GetMaximum() for _, hist in displayed.values())
    legend = ROOT.TLegend(0.55, 0.67, 0.88, 0.88)
    set_legend_style(legend, style=style)
    for index, (role, (label, histogram)) in enumerate(displayed.items()):
        set_unfolding_1d_style(histogram, role, style=style)
        histogram.SetTitle(";global #eta_{CM} bin;Entries")
        if index == 0:
            histogram.SetMaximum(1.25 * maximum if maximum > 0.0 else 1.0)
        histogram.Draw("E1" if index == 0 else "E1 SAME")
        legend.AddEntry(histogram, label, "p")
    legend.Draw()
    left_labels = draw_text_block(ROOT.gPad, annotations, style=style)
    canvas.cd(2)
    set_pad_style(ROOT.gPad, grid_x=grid, grid_y=grid, style=style)
    ROOT.gPad.SetRightMargin(style.palette_right_margin)
    set_2d_style(response, style=style)
    response.Draw("COLZ")
    right_labels = draw_text_block(ROOT.gPad, annotations, style=style)
    canvas.Modified()
    canvas.Update()
    save_canvas(canvas, output, save_png=save_png)
    canvas._unfolding_objects = [
        legend, *left_labels, *right_labels, response,
        *(hist for _, hist in displayed.values()),
    ]
    return canvas


def draw_unfolding_closure(
    target,
    measured,
    unfolded,
    *,
    target_label: str = "Gen",
    target_role: str = "gen",
    measured_label: str = "Reco",
    x_title: str = "#eta_{CM}",
    ratio_range: tuple[float, float] = (0.5, 1.5),
    x_range: tuple[float, float] | None = None,
    annotations: Iterable[str] = (),
    output: str | Path | None = None,
    save_png: bool = False,
    grid: bool = False,
    canvas_name: str = "canvas_unfolding_closure",
    ratio_name_prefix: str = "hClosure",
    style: PlotStyle = DEFAULT_PLOT_STYLE,
):
    if ratio_range[0] >= ratio_range[1]:
        raise ValueError("ratio_range must be increasing")
    measured_ratio = measured.Clone(f"{ratio_name_prefix}MeasuredToTarget")
    unfolded_ratio = unfolded.Clone(f"{ratio_name_prefix}UnfoldedToTarget")
    measured_ratio.SetDirectory(0)
    unfolded_ratio.SetDirectory(0)
    measured_ratio.Divide(target)
    unfolded_ratio.Divide(target)
    for histogram, role in (
        (target, target_role), (measured, "reco"), (unfolded, "unfolded"),
        (measured_ratio, "reco"), (unfolded_ratio, "unfolded"),
    ):
        set_unfolding_1d_style(histogram, role, style=style)

    canvas = ROOT.TCanvas(canvas_name, "Unfolding closure", 800, 800)
    top = ROOT.TPad(f"{canvas_name}_top", "", 0.0, style.ratio_fraction, 1.0, 1.0)
    bottom = ROOT.TPad(f"{canvas_name}_bottom", "", 0.0, 0.0, 1.0, style.ratio_fraction)
    top.Draw()
    bottom.Draw()
    top.cd()
    set_pad_style(top, grid_x=grid, grid_y=grid, style=style)
    top.SetBottomMargin(style.ratio_top_bottom_margin)
    maximum = max(target.GetMaximum(), measured.GetMaximum(), unfolded.GetMaximum())
    target.SetMaximum(1.30 * maximum if maximum > 0.0 else 1.0)
    target.SetTitle(f";{x_title};Entries")
    target.GetXaxis().SetLabelSize(0.0)
    if x_range is not None:
        target.GetXaxis().SetRangeUser(*x_range)
    target.Draw("E1")
    measured.Draw("E1 SAME")
    unfolded.Draw("E1 SAME")
    legend = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    set_legend_style(legend, style=style)
    legend.AddEntry(target, target_label, "p")
    legend.AddEntry(measured, measured_label, "p")
    legend.AddEntry(unfolded, "Unfolded", "p")
    legend.Draw()
    labels = draw_text_block(top, annotations, style=style)
    bottom.cd()
    set_pad_style(bottom, grid_x=grid, grid_y=grid, style=style)
    bottom.SetTopMargin(style.ratio_bottom_top_margin)
    bottom.SetBottomMargin(style.ratio_bottom_margin)
    for ratio in (measured_ratio, unfolded_ratio):
        ratio.GetXaxis().SetTitleSize(style.ratio_axis_title_size)
        ratio.GetXaxis().SetLabelSize(style.ratio_axis_label_size)
        ratio.GetYaxis().SetTitleSize(style.ratio_axis_title_size)
        ratio.GetYaxis().SetLabelSize(style.ratio_axis_label_size)
        ratio.GetYaxis().SetTitleOffset(style.ratio_y_title_offset)
        ratio.GetYaxis().SetNdivisions(style.ratio_axis_divisions)
        if x_range is not None:
            ratio.GetXaxis().SetRangeUser(*x_range)
    measured_ratio.SetTitle(f";{x_title};Ratio to {target_label}")
    measured_ratio.GetYaxis().SetRangeUser(*ratio_range)
    measured_ratio.Draw("E1")
    unfolded_ratio.Draw("E1 SAME")
    line_low, line_high = x_range or (
        measured_ratio.GetXaxis().GetXmin(), measured_ratio.GetXaxis().GetXmax(),
    )
    line = ROOT.TLine(line_low, 1.0, line_high, 1.0)
    line.SetLineStyle(2)
    line.Draw()
    canvas.Modified()
    canvas.Update()
    save_canvas(canvas, output, save_png=save_png)
    canvas._unfolding_objects = [
        top, bottom, legend, *labels, line, target, measured, unfolded,
        measured_ratio, unfolded_ratio,
    ]
    return canvas, measured_ratio, unfolded_ratio


def draw_unfolding_closure_by_pt(
    flattened_unfolded,
    targets: Sequence,
    measured: Sequence,
    layout: FlattenedBinning,
    *,
    output_dir: str | Path,
    output_tag: str,
    target_label: str = "Gen",
    target_role: str = "gen",
    measured_label: str = "Reco",
    eta_range: tuple[float, float] | None = None,
    ratio_range: tuple[float, float] = (0.5, 1.5),
    annotation_prefix: Iterable[str] = (),
    save_png: bool = False,
    grid: bool = False,
):
    if len(targets) != layout.n_pt_bins or len(measured) != layout.n_pt_bins:
        raise ValueError("Per-pT inputs do not match flattening layout")
    unfolded_by_pt = []
    measured_ratios = []
    unfolded_ratios = []
    canvases = []
    for pt_index, ((low, high), target, reco) in enumerate(
        zip(layout.pt_intervals, targets, measured)
    ):
        unfolded = extract_eta_block(
            flattened_unfolded, target, layout, pt_index,
            name=f"hUnfoldedEtaCM_pt{pt_index}",
        )
        annotations = (
            *tuple(annotation_prefix),
            f"{low:g} < p_{{T}}^{{ave}} < {high:g} GeV",
        )
        canvas, reco_ratio, unfolded_ratio = draw_unfolding_closure(
            target, reco, unfolded,
            target_label=target_label, target_role=target_role,
            measured_label=measured_label, ratio_range=ratio_range,
            x_range=eta_range, annotations=annotations,
            output=Path(output_dir) / f"{output_tag}_eta_ratio_pt_{low:g}_{high:g}.pdf",
            save_png=save_png, grid=grid,
            canvas_name=f"canvas_unfolding_closure_pt{pt_index}",
            ratio_name_prefix=f"hClosurePt{pt_index}",
        )
        unfolded_by_pt.append(unfolded)
        measured_ratios.append(reco_ratio)
        unfolded_ratios.append(unfolded_ratio)
        canvases.append(canvas)
    return unfolded_by_pt, measured_ratios, unfolded_ratios, canvases
