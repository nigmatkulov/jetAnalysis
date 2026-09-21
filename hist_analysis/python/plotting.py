"""Two-pad closure overlays and ratios to an explicit nominal curve."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping

import ROOT

from hist_analysis.python.histogram_ops import ratio_to_nominal, style_histograms
from hist_analysis.python.root_style import (
    DEFAULT_PLOT_STYLE, PlotStyle, apply_curve_style, draw_cms_header,
    draw_text_block, save_canvas, set_1d_style, set_legend_style, set_pad_style,
    style_single_panel_axes,
)


def figure_title(generator: str, direction: str, observable: str) -> str:
    return f"{generator} {direction} {observable}"


def _remember_display_ranges(pad, histograms: Iterable[object]) -> None:
    """Attach the final ROOT pad ranges for later histogram-aware exports."""

    x_range = (float(pad.GetUxmin()), float(pad.GetUxmax()))
    y_range = (float(pad.GetUymin()), float(pad.GetUymax()))
    if pad.GetLogx():
        x_range = (10.0 ** x_range[0], 10.0 ** x_range[1])
    if pad.GetLogy():
        y_range = (10.0 ** y_range[0], 10.0 ** y_range[1])
    for histogram in histograms:
        histogram._root_display_x_range = x_range
        histogram._root_display_y_range = y_range


def draw_overlay(histograms: Mapping[str, object], *, title: str,
                 x_title: str, y_title: str,
                 grid: bool = True, log_y: bool = False,
                 output: str | Path | None = None, save_png: bool = False,
                 headroom: float = 1.3, canvas_name: str | None = None,
                 annotations: Iterable[str] = (),
                 x_range: tuple[float, float] | None = None,
                 y_range: tuple[float, float] | None = None,
                 legend_bounds: tuple[float, float, float, float] | None = None,
                 reference_y: float | None = None,
                 overlay_functions: Mapping[str, object] | None = None,
                 overlay_text: Mapping[str, Iterable[str]] | None = None,
                 overlay_text_bounds: tuple[float, float, float, float] = (
                     0.43, 0.18, 0.88, 0.40,
                 ),
                 overlay_text_size: float = 0.018,
                 style_indices: Mapping[str, int] | None = None,
                 show_legend: bool = True,
                 legend_columns: int = 1,
                 legend_text_size: float | None = None,
                 right_margin: float | None = None,
                 curve_styles: Mapping[str, tuple[int, int]] | None = None,
                 cms_header_left: str | None = None,
                 cms_header_right: str | None = None,
                 style: PlotStyle = DEFAULT_PLOT_STYLE):
    """Draw histograms with optional same-label functions and text box."""

    if not histograms:
        raise ValueError("No histograms supplied")
    if headroom <= 1.0:
        raise ValueError("headroom must be greater than 1")
    if x_range is not None and x_range[0] >= x_range[1]:
        raise ValueError("x_range must satisfy low < high")
    if y_range is not None and y_range[0] >= y_range[1]:
        raise ValueError("y_range must satisfy low < high")
    if isinstance(legend_columns, bool) or not isinstance(legend_columns, int) or legend_columns < 1:
        raise ValueError("legend_columns must be a positive integer")
    if legend_text_size is not None and legend_text_size <= 0.0:
        raise ValueError("legend_text_size must be positive")
    if right_margin is not None and not 0.0 < right_margin < 0.8:
        raise ValueError("right_margin must lie between 0 and 0.8")
    if overlay_functions is not None:
        unknown_labels = set(overlay_functions) - set(histograms)
        if unknown_labels:
            raise KeyError(
                f"Functions have no matching histograms: {sorted(unknown_labels)}"
            )
    if overlay_text is not None:
        unknown_labels = set(overlay_text) - set(histograms)
        if unknown_labels:
            raise KeyError(
                f"Text has no matching histograms: {sorted(unknown_labels)}"
            )
        x1, y1, x2, y2 = overlay_text_bounds
        if not (0.0 <= x1 < x2 <= 1.0 and 0.0 <= y1 < y2 <= 1.0):
            raise ValueError("overlay_text_bounds must be ordered NDC coordinates")
        if overlay_text_size <= 0.0:
            raise ValueError("overlay_text_size must be positive")

    suffix = str(abs(hash(canvas_name or title)))
    canvas = ROOT.TCanvas(
        f"c_overlay_{suffix}", title,
        style.canvas_width, style.canvas_height,
    )
    set_pad_style(canvas, grid_x=grid, grid_y=grid, style=style)
    canvas.SetLeftMargin(style.single_panel_left_margin)
    canvas.SetBottomMargin(style.single_panel_bottom_margin)
    if right_margin is not None:
        canvas.SetRightMargin(right_margin)
    canvas.SetLogy(log_y)

    maximum = max(
        histogram.GetBinContent(bin_index) + histogram.GetBinError(bin_index)
        for histogram in histograms.values()
        for bin_index in range(1, histogram.GetNbinsX() + 1)
    )
    legend = ROOT.TLegend(*(legend_bounds or (0.65, 0.75, 0.88, 0.85)))
    set_legend_style(legend, style=style)
    legend.SetNColumns(legend_columns)
    if legend_text_size is not None:
        legend.SetTextSize(legend_text_size)
    for index, (label, histogram) in enumerate(histograms.items()):
        style_index = index if style_indices is None else style_indices.get(label, index)
        if curve_styles is not None and label in curve_styles:
            apply_curve_style(histogram, *curve_styles[label], style=style)
        else:
            set_1d_style(histogram, style_index, style=style)
        style_single_panel_axes(histogram, style=style)
        histogram.SetTitle(title)
        histogram.GetXaxis().SetTitle(x_title)
        histogram.GetYaxis().SetTitle(y_title)
        if x_range is not None:
            histogram.GetXaxis().SetRangeUser(*x_range)
        if y_range is not None:
            histogram.SetMinimum(y_range[0])
            histogram.SetMaximum(y_range[1])
        else:
            histogram.SetMaximum(maximum * (20.0 if log_y else headroom))
        histogram.Draw("E1" if index == 0 else "E1 SAME")
        legend.AddEntry(histogram, label, "p")
    drawn_functions = []
    for label, function in (overlay_functions or {}).items():
        function.SetLineColor(histograms[label].GetLineColor())
        function.Draw("SAME")
        drawn_functions.append(function)
    if show_legend:
        legend.Draw()
    annotation_objects = draw_text_block(canvas, annotations, style=style)
    header_objects = draw_cms_header(
        canvas, cms_header_left, cms_header_right, style=style,
    )

    line = None
    if reference_y is not None:
        axis = next(iter(histograms.values())).GetXaxis()
        line_low, line_high = x_range or (axis.GetXmin(), axis.GetXmax())
        line = ROOT.TLine(line_low, reference_y, line_high, reference_y)
        line.SetLineStyle(2)
        line.Draw()

    text_box = None
    text_lines = []
    if overlay_text:
        text_box = ROOT.TPaveText(*overlay_text_bounds, "NDC NB")
        text_box.SetFillColor(ROOT.kWhite)
        text_box.SetFillStyle(1001)
        text_box.SetBorderSize(0)
        text_box.SetTextAlign(12)
        text_box.SetTextFont(style.font)
        text_box.SetTextSize(overlay_text_size)
        for label, lines in overlay_text.items():
            for text in lines:
                rendered = text_box.AddText(text)
                rendered.SetTextColor(histograms[label].GetLineColor())
                text_lines.append(rendered)
        text_box.Draw()

    canvas.Modified()
    canvas.Update()
    _remember_display_ranges(canvas, histograms.values())
    save_canvas(canvas, output, save_png=save_png)
    canvas._overlay_objects = [
        *([legend] if show_legend else []),
        *annotation_objects, *header_objects, *histograms.values(), *drawn_functions,
        *([] if line is None else [line]),
        *([] if text_box is None else [text_box, *text_lines]),
    ]
    return canvas


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
                 y_range: tuple[float, float] | None = None,
                 ratio_option: str = "",
                 ratio_pairs: Mapping[str, tuple[str, str]] | None = None,
                 ratio_y_title: str | None = None,
                 legend_bounds: tuple[float, float, float, float] | None = None,
                 annotation_x: float | None = None,
                 annotation_y: float | None = None,
                 right_margin: float | None = None,
                 style_indices: Mapping[str, int] | None = None,
                 curve_styles: Mapping[str, tuple[int, int]] | None = None,
                 cms_header_left: str | None = None,
                 cms_header_right: str | None = None,
                 style: PlotStyle = DEFAULT_PLOT_STYLE):
    """Draw an overlay and ratios to an explicitly named nominal histogram.

    ``draw_nominal_ratio=False`` suppresses the trivial nominal/self ratio.
    ``ratio_pairs`` optionally maps displayed ratio labels to explicit
    ``(numerator, denominator)`` histogram labels instead of dividing every
    curve by ``nominal``.
    ``ratio_y_title`` overrides the automatically generated lower-panel title.
    ``legend_bounds`` and the annotation coordinates allow crowded overlays to
    reserve separate non-overlapping regions for labels.
    ``right_margin`` can reserve a sidebar for those labels on both pads.
    ``x_range`` is applied to both panels, while ``reference_line_x_range`` can
    independently constrain the horizontal ratio reference line. ``y_range``
    constrains the upper-panel distributions. ``ratio_option`` is passed to
    ROOT histogram division: ``''`` selects regular independent errors and
    ``'B'`` selects binomial errors. ``style_indices`` can assign
    explicit ``set_1d_style`` indices by curve label. ``annotations`` are drawn
    as a shared-style text block in the upper pad. ``headroom`` reserves vertical
    space for legends and annotations, and ``style`` controls the common canvas,
    margin, font, axis, and text geometry.
    """

    if nominal not in histograms:
        raise KeyError(f"Nominal {nominal!r} is not among {list(histograms)}")
    if not histograms:
        raise ValueError("No histograms supplied")
    if ratio_pairs is not None:
        missing_ratio_labels = {
            label
            for pair in ratio_pairs.values()
            for label in pair
            if label not in histograms
        }
        if missing_ratio_labels:
            raise KeyError(
                "Ratio pairs reference unknown histograms: "
                f"{sorted(missing_ratio_labels)}"
            )

    if curve_styles is not None:
        for index, (label, histogram) in enumerate(histograms.items()):
            if label in curve_styles:
                apply_curve_style(histogram, *curve_styles[label], style=style)
            else:
                set_1d_style(histogram, index, style=style)
    elif style_indices is None:
        style_histograms(histograms)
    else:
        for index, (label, histogram) in enumerate(histograms.items()):
            set_1d_style(
                histogram, style_indices.get(label, index), style=style,
            )
    if headroom <= 1.0:
        raise ValueError("headroom must be greater than 1")
    if x_range is not None and x_range[0] >= x_range[1]:
        raise ValueError("x_range must satisfy low < high")
    if y_range is not None and y_range[0] >= y_range[1]:
        raise ValueError("y_range must satisfy low < high")
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
    top.SetTopMargin(style.top_margin)
    top.SetLeftMargin(style.ratio_left_margin)
    top.SetRightMargin(style.right_margin if right_margin is None else right_margin)
    bottom.SetTopMargin(style.ratio_bottom_top_margin)
    bottom.SetBottomMargin(style.ratio_bottom_margin)
    bottom.SetLeftMargin(style.ratio_left_margin)
    bottom.SetRightMargin(style.right_margin if right_margin is None else right_margin)
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
    legend = ROOT.TLegend(*(legend_bounds or (0.65, 0.75, 0.88, 0.85)))
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
        if y_range is None:
            histogram.SetMaximum(maximum * (20.0 if log_y else headroom))
        else:
            histogram.SetMinimum(y_range[0])
            histogram.SetMaximum(y_range[1])
        histogram.Draw("E1" if index == 0 else "E1 SAME")
        legend.AddEntry(histogram, label, "p")
    legend.Draw()
    annotation_objects = draw_text_block(
        top, annotation_lines, x=annotation_x, y=annotation_y, style=style,
    )
    canvas.cd()
    canvas.SetLeftMargin(top.GetLeftMargin())
    canvas.SetRightMargin(top.GetRightMargin())
    header_objects = draw_cms_header(
        canvas, cms_header_left, cms_header_right, style=style,
    )

    bottom.cd()
    ratios = {}
    configured_ratios = tuple(
        ratio_pairs.items() if ratio_pairs is not None else (
            (label, (label, nominal))
            for label in histograms
            if label != nominal or draw_nominal_ratio
        )
    )
    for index, (label, (numerator_label, denominator_label)) in enumerate(
        configured_ratios
    ):
        ratio = ratio_to_nominal(
            histograms[numerator_label], histograms[denominator_label],
            name=f"ratio_{index}",
            option=ratio_option,
        )
        if curve_styles is not None and numerator_label in curve_styles:
            apply_curve_style(
                ratio, *curve_styles[numerator_label], style=style,
            )
        if ratio_pairs is not None:
            repeated_numerator = sum(
                pair[0] == numerator_label
                for _, pair in configured_ratios
            ) > 1
            style_label = denominator_label if repeated_numerator else numerator_label
            ratio_style_index = (
                style_indices.get(style_label, index)
                if style_indices is not None else index
            )
            set_1d_style(ratio, ratio_style_index, style=style)
        ratios[label] = ratio
        ratio.SetTitle("")
        ratio.GetYaxis().SetTitle(
            ratio_y_title or (
                "Ratio" if ratio_pairs is not None else f"Ratio to {nominal}"
            )
        )
        ratio.GetXaxis().SetTitle(x_title)
        ratio.SetMinimum(ratio_range[0])
        ratio.SetMaximum(ratio_range[1])
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
    _remember_display_ranges(top, histograms.values())
    _remember_display_ranges(bottom, ratios.values())
    save_canvas(canvas, output, save_png=save_png)
    # Keep drawn ROOT objects alive for notebook display.
    canvas._closure_objects = [
        top, bottom, legend, line, *annotation_objects, *header_objects,
        *histograms.values(), *ratios.values(),
    ]
    return canvas, ratios
