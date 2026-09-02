"""Plots shared by flattened dijet unfolding notebooks.

Plotting functions clone ratio/display objects when necessary and retain them on
their canvas to satisfy PyROOT lifetime rules.  They never renormalize the
histograms used as response or unfolding inputs.  Unless ``normalize=True`` is
explicitly requested, closure ratios compare stored weighted yields bin by bin.
"""

from __future__ import annotations

import math
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


def range_covering_histogram_bins(histogram, requested_range):
    """Expand limits inside bins to bin edges while preserving exact edges."""

    if requested_range is None:
        return None
    low, high = requested_range
    axis = histogram.GetXaxis()
    edges = [axis.GetBinLowEdge(index) for index in range(1, axis.GetNbins() + 1)]
    edges.append(axis.GetBinUpEdge(axis.GetNbins()))
    tolerance = 1.0e-10 * max(1.0, abs(low), abs(high), abs(edges[-1] - edges[0]))

    def snapped_edge(value):
        return next(
            (edge for edge in edges if math.isclose(value, edge, abs_tol=tolerance)),
            None,
        )

    snapped_low = snapped_edge(low)
    snapped_high = snapped_edge(high)
    if snapped_low is None:
        low_bin = min(max(axis.FindFixBin(low), 1), axis.GetNbins())
        snapped_low = axis.GetBinLowEdge(low_bin)
    if snapped_high is None:
        high_bin = min(max(axis.FindFixBin(high), 1), axis.GetNbins())
        snapped_high = axis.GetBinUpEdge(high_bin)
    return snapped_low, snapped_high


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
    """Draw normalized producer-level match/miss/fake category weights."""
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
    """Show one pTave projection beside its diagonal eta-response block.

    This response panel is diagnostic only.  The actual flattened unfolding
    matrix also contains off-diagonal pTave migrations.
    """
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
    """Show flattened spectra and the complete two-dimensional response."""
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


def draw_response_matrix(
    response,
    *,
    annotations: Iterable[str] = (),
    output: str | Path | None = None,
    save_png: bool = False,
    grid: bool = False,
    log_z: bool = True,
    canvas_name: str = "canvas_response_matrix",
    title: str = ";global reco #eta_{CM} bin;global gen #eta_{CM} bin",
    style: PlotStyle = DEFAULT_PLOT_STYLE,
):
    """Draw the complete measured-bin versus truth-bin migration matrix."""
    """Draw and save a flattened 2D response matrix on its own canvas."""

    if not response.InheritsFrom("TH2"):
        raise TypeError(f"{response.GetName()} is not a TH2 response matrix")
    canvas = ROOT.TCanvas(canvas_name, "Flattened response matrix", 900, 800)
    set_pad_style(canvas, grid_x=grid, grid_y=grid, style=style)
    canvas.SetRightMargin(style.palette_right_margin)
    canvas.SetLogz(log_z)
    set_2d_style(response, style=style)
    response.SetTitle(title)
    if log_z:
        positive = [
            response.GetBinContent(x_bin, y_bin)
            for x_bin in range(1, response.GetNbinsX() + 1)
            for y_bin in range(1, response.GetNbinsY() + 1)
            if response.GetBinContent(x_bin, y_bin) > 0.0
        ]
        if not positive:
            raise ValueError("Cannot draw an empty response matrix with log_z=True")
        response.SetMinimum(min(positive))
    response.Draw("COLZ")
    labels = draw_text_block(canvas, annotations, style=style)
    canvas.Modified()
    canvas.Update()
    save_canvas(canvas, output, save_png=save_png)
    canvas._unfolding_objects = [response, *labels]
    return canvas


def draw_unfolding_closure(
    target,
    measured,
    unfolded,
    *,
    target_label: str = "Gen",
    target_role: str = "gen",
    measured_label: str = "Reco",
    unfolded_label: str = "Unfolded",
    ratio_target_label: str | None = None,
    x_title: str = "#eta_{CM}",
    y_title: str = "Entries",
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
    """Compare measured and unfolded yields with a truth/reference histogram.

    Ratios use ROOT's standard independent-error propagation for display.  The
    unfolded bins are correlated, so quantitative inference should use the
    saved unfolding covariance rather than treating these ratio points as
    statistically independent.
    """
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
    target.SetTitle(f";{x_title};{y_title}")
    target.GetXaxis().SetLabelSize(0.0)
    displayed_x_range = range_covering_histogram_bins(target, x_range)
    if displayed_x_range is not None:
        target.GetXaxis().SetRangeUser(*displayed_x_range)
    target.Draw("E1")
    measured.Draw("E1 SAME")
    unfolded.Draw("E1 SAME")
    legend = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    set_legend_style(legend, style=style)
    legend.AddEntry(target, target_label, "p")
    legend.AddEntry(measured, measured_label, "p")
    legend.AddEntry(unfolded, unfolded_label, "p")
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
        if displayed_x_range is not None:
            ratio.GetXaxis().SetRangeUser(*displayed_x_range)
    ratio_axis_label = ratio_target_label or target_label
    measured_ratio.SetTitle(f";{x_title};Ratio to {ratio_axis_label}")
    measured_ratio.GetYaxis().SetRangeUser(*ratio_range)
    measured_ratio.Draw("E1")
    unfolded_ratio.Draw("E1 SAME")
    line_low, line_high = displayed_x_range or (
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
    ratio_target_label: str | None = None,
    unfolded_label: str = "Unfolded",
    normalize: bool = False,
    output_suffix: str = "eta_ratio",
    object_name_prefix: str = "ClosurePt",
    eta_range: tuple[float, float] | None = None,
    ratio_range: tuple[float, float] = (0.5, 1.5),
    annotation_prefix: Iterable[str] = (),
    save_png: bool = False,
    grid: bool = False,
):
    """Extract each pTave block and draw its eta closure.

    ``normalize=False`` tests yield closure.  ``normalize=True`` clones and
    independently unit-normalizes the three displayed shapes, intentionally
    removing normalization closure from that diagnostic.
    """
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
            name=f"h{object_name_prefix}UnfoldedEtaCM_pt{pt_index}",
        )
        displayed_target = target
        displayed_reco = reco
        displayed_unfolded = unfolded
        if normalize:
            normalized = []
            for role, histogram in (
                ("Target", target), ("Reco", reco), ("Unfolded", unfolded),
            ):
                clone = histogram.Clone(
                    f"h{object_name_prefix}{role}NormalizedPt{pt_index}"
                )
                clone.SetDirectory(0)
                integral = clone.Integral("width")
                if integral <= 0.0:
                    raise ValueError(
                        f"Cannot normalize {histogram.GetName()} with integral {integral}"
                    )
                clone.Scale(1.0 / integral)
                normalized.append(clone)
            displayed_target, displayed_reco, displayed_unfolded = normalized
        annotations = (
            *tuple(annotation_prefix),
            f"{low:g} < p_{{T}}^{{ave}} < {high:g} GeV",
        )
        canvas, reco_ratio, unfolded_ratio = draw_unfolding_closure(
            displayed_target, displayed_reco, displayed_unfolded,
            target_label=target_label, target_role=target_role,
            measured_label=measured_label, unfolded_label=unfolded_label,
            y_title="1/N dN/d#eta_{CM}" if normalize else "Entries",
            ratio_range=ratio_range,
            ratio_target_label=ratio_target_label,
            x_range=eta_range, annotations=annotations,
            output=Path(output_dir) / f"{output_tag}_{output_suffix}_pt_{low:g}_{high:g}.pdf",
            save_png=save_png, grid=grid,
            canvas_name=f"canvas_{object_name_prefix}_pt{pt_index}",
            ratio_name_prefix=f"h{object_name_prefix}{pt_index}",
        )
        unfolded_by_pt.append(displayed_unfolded)
        measured_ratios.append(reco_ratio)
        unfolded_ratios.append(unfolded_ratio)
        canvases.append(canvas)
    return unfolded_by_pt, measured_ratios, unfolded_ratios, canvases


def draw_response_components(
    *,
    explicit_miss,
    boundary_miss,
    effective_miss,
    explicit_fake,
    boundary_fake,
    effective_fake,
    annotations: Iterable[str] = (),
    output: str | Path | None = None,
    save_png: bool = False,
    grid: bool = True,
    canvas_name: str = "canvas_response_components",
    style: PlotStyle = DEFAULT_PLOT_STYLE,
):
    """Draw explicit, boundary, and effective miss/fake spectra separately."""

    canvas = ROOT.TCanvas(
        canvas_name, "Response miss and fake components",
        style.side_by_side_canvas_width, style.side_by_side_canvas_height,
    )
    canvas.Divide(2, 1)
    retained = []
    for pad_index, (title, role, components) in enumerate((
        ("Miss", "miss", (
            ("Explicit selection/matching", explicit_miss, 2),
            ("p_{T}^{ave} boundary migration", boundary_miss, 3),
            ("Effective (used by RooUnfold)", effective_miss, 1),
        )),
        ("Fake", "fake", (
            ("Explicit selection/matching", explicit_fake, 2),
            ("p_{T}^{ave} boundary migration", boundary_fake, 3),
            ("Effective (used by RooUnfold)", effective_fake, 1),
        )),
    ), 1):
        canvas.cd(pad_index)
        set_pad_style(ROOT.gPad, grid_x=grid, grid_y=grid, style=style)
        maximum = max(histogram.GetMaximum() for _, histogram, _ in components)
        legend = ROOT.TLegend(0.48, 0.68, 0.88, 0.88)
        set_legend_style(legend, style=style)
        for index, (label, histogram, line_style) in enumerate(components):
            set_unfolding_1d_style(histogram, role, style=style)
            histogram.SetLineStyle(line_style)
            histogram.SetLineWidth(2)
            histogram.SetTitle(f";global #eta_{{CM}} bin;{title} weighted events")
            if index == 0:
                histogram.SetMaximum(1.25 * maximum if maximum > 0.0 else 1.0)
            histogram.Draw("HIST" if index == 0 else "HIST SAME")
            legend.AddEntry(histogram, label, "l")
        legend.Draw()
        labels = draw_text_block(ROOT.gPad, annotations, style=style)
        retained.extend((legend, *labels, *(hist for _, hist, _ in components)))
    canvas.Modified()
    canvas.Update()
    save_canvas(canvas, output, save_png=save_png)
    canvas._unfolding_objects = retained
    return canvas


def draw_controlled_studied_shape_comparisons(
    controlled_truth: Sequence,
    studied_truth: Sequence,
    controlled_reco: Sequence,
    studied_reco: Sequence,
    pt_intervals: Sequence[tuple[float, float]],
    *,
    output_dir: str | Path,
    output_tag: str,
    eta_range: tuple[float, float] | None = None,
    ratio_range: tuple[float, float] = (0.5, 1.5),
    annotations: Iterable[str] = (),
    save_png: bool = False,
    grid: bool = True,
    style: PlotStyle = DEFAULT_PLOT_STYLE,
):
    """Compare normalized shapes and studied/controlled ratios by pT interval."""

    if ratio_range[0] >= ratio_range[1]:
        raise ValueError("ratio_range must be increasing")

    collections = (
        tuple(controlled_truth), tuple(studied_truth),
        tuple(controlled_reco), tuple(studied_reco),
    )
    n_intervals = len(tuple(pt_intervals))
    if any(len(collection) != n_intervals for collection in collections):
        raise ValueError("Shape-comparison inputs do not match pT intervals")
    output_dir = Path(output_dir)
    canvases = []
    normalized_histograms = []
    ratio_histograms = []
    for pt_index, (low, high) in enumerate(pt_intervals):
        canvas = ROOT.TCanvas(
            f"canvas_controlled_studied_shapes_pt{pt_index}", "",
            800, 800,
        )
        canvas.Divide(2, 2, 0.001, 0.001)
        retained = []
        for column, (role, legend_role, controlled, studied) in enumerate((
            ("gen", "Gen", controlled_truth[pt_index], studied_truth[pt_index]),
            ("reco", "Reco", controlled_reco[pt_index], studied_reco[pt_index]),
        )):
            top_pad_index = column + 1
            ratio_pad_index = column + 3
            canvas.cd(top_pad_index)
            set_pad_style(ROOT.gPad, grid_x=grid, grid_y=grid, style=style)
            displayed = []
            for sample, source, line_style, marker_style in (
                ("Controlled", controlled, 1, 20),
                ("Studied", studied, 1, 24),
            ):
                histogram = source.Clone(
                    f"{source.GetName()}Normalized{sample}Pt{pt_index}"
                )
                histogram.SetDirectory(0)
                integral = histogram.Integral("width")
                if integral <= 0.0:
                    raise ValueError(
                        f"Cannot normalize {source.GetName()} with integral {integral}"
                    )
                histogram.Scale(1.0 / integral)
                set_unfolding_1d_style(histogram, role, style=style)
                histogram.SetLineStyle(line_style)
                histogram.SetMarkerStyle(marker_style)
                histogram.SetTitle(";#eta_{CM};1/N dN/d#eta_{CM}")
                displayed.append((sample, histogram))
            displayed_eta_range = range_covering_histogram_bins(
                displayed[0][1], eta_range,
            )
            if displayed_eta_range is not None:
                for _, histogram in displayed:
                    histogram.GetXaxis().SetRangeUser(*displayed_eta_range)
            maximum = max(histogram.GetMaximum() for _, histogram in displayed)
            displayed[0][1].SetMaximum(1.25 * maximum if maximum > 0.0 else 1.0)
            legend = ROOT.TLegend(0.62, 0.75, 0.88, 0.88)
            set_legend_style(legend, style=style)
            for index, (sample, histogram) in enumerate(displayed):
                histogram.Draw("E1" if index == 0 else "E1 SAME")
                legend.AddEntry(histogram, f"{sample} {legend_role}", "pl")
            legend.Draw()
            labels = draw_text_block(
                ROOT.gPad,
                (*tuple(annotations), f"{low:g} < p_{{T}}^{{ave}} < {high:g} GeV"),
                style=style,
            )
            retained.extend((legend, *labels, *(hist for _, hist in displayed)))
            normalized_histograms.extend(hist for _, hist in displayed)

            controlled_normalized = displayed[0][1]
            studied_normalized = displayed[1][1]
            ratio = studied_normalized.Clone(
                f"hStudiedToControlled{legend_role}Pt{pt_index}"
            )
            ratio.SetDirectory(0)
            ratio.Divide(controlled_normalized)
            set_unfolding_1d_style(ratio, "unfolded", style=style)
            ratio.SetLineStyle(1)
            ratio.SetMarkerStyle(20)
            ratio.SetTitle(";#eta_{CM};Studied / Controlled")
            ratio.GetYaxis().SetRangeUser(*ratio_range)
            if displayed_eta_range is not None:
                ratio.GetXaxis().SetRangeUser(*displayed_eta_range)
            canvas.cd(ratio_pad_index)
            set_pad_style(ROOT.gPad, grid_x=grid, grid_y=grid, style=style)
            ratio.Draw("E1")
            line_low, line_high = displayed_eta_range or (
                ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax(),
            )
            unity = ROOT.TLine(line_low, 1.0, line_high, 1.0)
            unity.SetLineStyle(2)
            unity.Draw()
            retained.extend((ratio, unity))
            ratio_histograms.append(ratio)
        canvas.Modified()
        canvas.Update()
        save_canvas(
            canvas,
            output_dir / f"{output_tag}_controlled_studied_shapes_pt_{low:g}_{high:g}.pdf",
            save_png=save_png,
        )
        canvas._unfolding_objects = retained
        canvases.append(canvas)
    return canvases, normalized_histograms, ratio_histograms


def draw_regularization_metrics(
    scan_result,
    *,
    mode: str = "relative",
    output: str | Path | None = None,
    save_png: bool = False,
    grid: bool = True,
    log_y: bool = False,
    canvas_name: str | None = None,
    style: PlotStyle = DEFAULT_PLOT_STYLE,
):
    """Plot bin-averaged bias squared, variance, and MSE versus iterations.

    Absolute values are the arithmetic mean of the per-bin metric and retain
    yield-squared units.  Relative values first divide each valid bin's metric
    by that bin's truth squared and then average, so they are dimensionless.
    """

    if mode not in {"absolute", "relative"}:
        raise ValueError("mode must be 'absolute' or 'relative'")
    names = (("bias_squared", "Bias^{2}"), ("variance", "Variance"), ("mse", "MSE"))
    roles = ("reco", "gen", "unfolded")
    histograms = []
    for index, ((key, label), role) in enumerate(zip(names, roles)):
        histogram = ROOT.TH1D(
            f"h{mode.capitalize()}{key.title().replace('_', '')}Plot", "",
            len(scan_result.metrics), 0.5, len(scan_result.metrics) + 0.5,
        )
        histogram.SetDirectory(0)
        for result in scan_result.metrics:
            histogram.SetBinContent(result.iterations, getattr(result, mode)[key])
        y_title = ("Mean fractional per-bin metric"
                   if mode == "relative" else "Mean per-bin metric [yield^{2}]")
        histogram.SetTitle(f";Bayesian iterations;{y_title}")
        set_unfolding_1d_style(histogram, role, style=style)
        histogram.SetMarkerStyle(20)
        histograms.append((label, histogram))
    canvas_name = canvas_name or f"canvas_regularization_{mode}"
    canvas = ROOT.TCanvas(canvas_name, f"{mode} regularization metrics", 800, 700)
    set_pad_style(canvas, grid_x=grid, grid_y=grid, style=style)
    canvas.SetLogy(log_y)
    maximum = max(histogram.GetMaximum() for _, histogram in histograms)
    positive_values = [
        histogram.GetBinContent(bin_index)
        for _, histogram in histograms
        for bin_index in range(1, histogram.GetNbinsX() + 1)
        if histogram.GetBinContent(bin_index) > 0.0
    ]
    if log_y and not positive_values:
        raise ValueError("Logarithmic metric plot has no positive values")
    legend = ROOT.TLegend(0.62, 0.70, 0.88, 0.88)
    set_legend_style(legend, style=style)
    for index, (label, histogram) in enumerate(histograms):
        if index == 0:
            histogram.SetMaximum(1.25 * maximum if maximum > 0.0 else 1.0)
            if log_y:
                histogram.SetMinimum(0.5 * min(positive_values))
        histogram.Draw("LP" if index == 0 else "LP SAME")
        legend.AddEntry(histogram, label, "lp")
    legend.Draw()
    # Each canvas marks the minimum belonging to the metric it displays. This
    # remains meaningful when the scan reports both absolute and relative optima.
    selected_iterations = (
        scan_result.absolute_selected_iterations if mode == "absolute"
        else scan_result.relative_selected_iterations
    ) or scan_result.selected_iterations
    selected = ROOT.TLine(selected_iterations, 0.0,
                          selected_iterations,
                          1.25 * maximum if maximum > 0.0 else 1.0)
    selected.SetLineStyle(2)
    selected.Draw()
    labels = draw_text_block(
        canvas, (f"{mode.capitalize()} MSE minimum: {selected_iterations}",),
        style=style,
    )
    canvas.Modified()
    canvas.Update()
    save_canvas(canvas, output, save_png=save_png)
    canvas._unfolding_objects = [legend, selected, *labels, *(h for _, h in histograms)]
    return canvas, [histogram for _, histogram in histograms]


def draw_toy_input_diagnostic(
    prior,
    folded,
    data_errors,
    toy,
    *,
    output: str | Path | None = None,
    save_png: bool = False,
    grid: bool = True,
    canvas_name: str = "canvas_toy_input_diagnostic",
    style: PlotStyle = DEFAULT_PLOT_STYLE,
):
    """Show truth and measured toy inputs in separate coordinate panels.

    Truth and reco use the same flattened bin numbering but are different
    observables.  Keeping them on separate pads prevents a visual overlay from
    being mistaken for a closure comparison.
    """

    canvas = ROOT.TCanvas(canvas_name, "Regularization toy inputs", 1400, 650)
    canvas.Divide(2, 1)
    canvas.cd(1)
    set_pad_style(ROOT.gPad, grid_x=grid, grid_y=grid, style=style)
    set_unfolding_1d_style(prior, "gen", style=style)
    prior.SetTitle(";global truth bin;Entries")
    prior.Draw("E1")
    truth_legend = ROOT.TLegend(0.58, 0.78, 0.88, 0.88)
    set_legend_style(truth_legend, style=style)
    truth_legend.AddEntry(prior, "Target distorted truth", "p")
    truth_legend.Draw()

    canvas.cd(2)
    set_pad_style(ROOT.gPad, grid_x=grid, grid_y=grid, style=style)
    reco_spectra = (
        (folded, "reco", "Inclusive-reco toy mean"),
        (data_errors, "fake", "Trigger data (error source)"),
        (toy, "unfolded", "Representative inclusive-reco toy"),
    )
    maximum = max(histogram.GetMaximum() for histogram, _, _ in reco_spectra)
    reco_legend = ROOT.TLegend(0.52, 0.70, 0.88, 0.88)
    set_legend_style(reco_legend, style=style)
    for index, (histogram, role, label) in enumerate(reco_spectra):
        set_unfolding_1d_style(histogram, role, style=style)
        histogram.SetTitle(";global measured bin;Entries")
        if index == 0:
            histogram.SetMaximum(1.25 * maximum if maximum > 0.0 else 1.0)
        histogram.Draw("E1" if index == 0 else "E1 SAME")
        reco_legend.AddEntry(histogram, label, "p")
    reco_legend.Draw()
    canvas.Modified()
    canvas.Update()
    save_canvas(canvas, output, save_png=save_png)
    canvas._unfolding_objects = [
        truth_legend, reco_legend, prior,
        *(histogram for histogram, _, _ in reco_spectra),
    ]
    return canvas


def draw_regularization_metrics_by_pt(
    scan_result,
    layout: FlattenedBinning,
    *,
    output_dir: str | Path,
    output_tag: str,
    mode: str = "relative",
    save_png: bool = False,
    grid: bool = True,
    log_y: bool = False,
    style: PlotStyle = DEFAULT_PLOT_STYLE,
):
    """Plot per-pT averages using the same absolute/relative definitions."""

    if mode not in {"absolute", "relative"}:
        raise ValueError("mode must be 'absolute' or 'relative'")
    output_dir = Path(output_dir)
    canvases = []
    all_histograms = []
    for pt_index, (low, high) in enumerate(layout.pt_intervals):
        metric_histograms = []
        first = layout.global_bin(pt_index, 1)
        last = layout.global_bin(pt_index, layout.n_eta_bins)
        for key, source_name, role, label in (
            ("bias_squared", "bias_squared", "reco", "Bias^{2}"),
            ("variance", "variance", "gen", "Variance"),
            ("mse", "mse", "unfolded", "MSE"),
        ):
            y_title = ("Mean fractional per-bin metric"
                       if mode == "relative" else "Mean per-bin metric [yield^{2}]")
            histogram = ROOT.TH1D(
                f"h{mode.capitalize()}{key.title().replace('_', '')}Pt{pt_index}",
                f";Bayesian iterations;{y_title}",
                len(scan_result.metrics), 0.5, len(scan_result.metrics) + 0.5,
            )
            histogram.SetDirectory(0)
            for result in scan_result.metrics:
                source = getattr(result, source_name)
                values = []
                for bin_index in range(first, last + 1):
                    value = source.GetBinContent(bin_index)
                    if mode == "relative":
                        truth_value = result.mean.GetBinContent(bin_index) - result.bias.GetBinContent(bin_index)
                        if bin_index not in result.valid_relative_bins:
                            continue
                        value /= truth_value * truth_value
                    values.append(value)
                histogram.SetBinContent(result.iterations, sum(values) / len(values))
            set_unfolding_1d_style(histogram, role, style=style)
            histogram.SetMarkerStyle(20)
            metric_histograms.append((label, histogram))
        canvas = ROOT.TCanvas(
            f"canvas_regularization_{mode}_pt{pt_index}", "", 800, 700,
        )
        set_pad_style(canvas, grid_x=grid, grid_y=grid, style=style)
        canvas.SetLogy(log_y)
        maximum = max(histogram.GetMaximum() for _, histogram in metric_histograms)
        positive_values = [
            histogram.GetBinContent(bin_index)
            for _, histogram in metric_histograms
            for bin_index in range(1, histogram.GetNbinsX() + 1)
            if histogram.GetBinContent(bin_index) > 0.0
        ]
        if log_y and not positive_values:
            raise ValueError(
                f"Logarithmic pT-block metric plot {pt_index} has no positive values"
            )
        legend = ROOT.TLegend(0.62, 0.70, 0.88, 0.88)
        set_legend_style(legend, style=style)
        for index, (label, histogram) in enumerate(metric_histograms):
            if index == 0:
                histogram.SetMaximum(1.25 * maximum if maximum > 0.0 else 1.0)
                if log_y:
                    histogram.SetMinimum(0.5 * min(positive_values))
            histogram.Draw("LP" if index == 0 else "LP SAME")
            legend.AddEntry(histogram, label, "lp")
        legend.Draw()
        labels = draw_text_block(
            canvas, (f"{low:g} < p_{{T}}^{{ave}} < {high:g} GeV",), style=style,
        )
        canvas.Modified()
        canvas.Update()
        save_canvas(
            canvas,
            output_dir / f"{output_tag}_{mode}_metrics_pt_{low:g}_{high:g}.pdf",
            save_png=save_png,
        )
        canvas._unfolding_objects = [
            legend, *labels, *(histogram for _, histogram in metric_histograms),
        ]
        canvases.append(canvas)
        all_histograms.extend(histogram for _, histogram in metric_histograms)
    return canvases, all_histograms
