"""Standalone ROOT and mplhep figures for dijet closure comparisons."""

from dataclasses import replace
from uuid import uuid4

import ROOT

from hist_analysis.python.histogram_ops import ratio_to_nominal
from hist_analysis.python.plotting import draw_overlay
from hist_analysis.python.root_style import (
    COLORS, DEFAULT_PLOT_STYLE, MARKERS, draw_cms_header, save_canvas,
)


# Shared palette index and ROOT marker code for each physics curve.
CURVE_STYLES = {
    "Reco": (0, 20),  # red filled circle
    "Ref": (1, 24),   # blue open circle
    "Gen": (5, 22),   # green filled upward triangle
    "RefSel": (3, 23),
}
FIGURE_STYLE = replace(
    DEFAULT_PLOT_STYLE,
    single_panel_left_margin=0.15,
    annotation_x=0.18,
    annotation_y=0.875,
    annotation_text_size=0.026,
    annotation_line_spacing=0.038,
    legend_text_size=0.026,
)


def closure_styles(labels):
    """Assign styles once in the configured order, including extra curves."""
    used = {color for color, _ in CURVE_STYLES.values()}
    available = iter(index for index in range(len(COLORS)) if index not in used)
    styles = {}
    for label in labels:
        curve_style = CURVE_STYLES.get(label)
        if curve_style is None:
            index = next(available, 0)
            curve_style = index, MARKERS[index]
        color, marker = curve_style
        styles[label] = COLORS[color], marker
    return styles


def signed_eta_plot_range(histogram, eta_cut):
    """Use whole bins that overlap the selected signed eta acceptance."""
    axis = histogram.GetXaxis()
    # Probe just inside the cuts so roundoff at an aligned bin edge cannot
    # select the adjacent bin. Keep the whole bin when a cut falls inside it.
    edge_tolerance = 1e-9
    first = max(1, axis.FindFixBin(-eta_cut + edge_tolerance))
    last = min(axis.GetNbins(), axis.FindFixBin(eta_cut - edge_tolerance))
    return axis.GetBinLowEdge(first), axis.GetBinUpEdge(last)


def _display_copies(histograms):
    copies = {}
    for label, histogram in histograms.items():
        copy = histogram.Clone()
        copy.SetDirectory(0)
        copies[label] = copy
    return copies


def draw_dijet_figure(histograms, *, output_dir, output_tag, x_title, y_title,
                      annotations, styles, grid, save_png, x_range=None,
                      y_range=None, reference_y=None, line_styles=None):
    """Draw selection text at upper left and the legend at upper right."""
    legend_top = FIGURE_STYLE.annotation_y
    legend_bottom = legend_top - 0.039 * len(histograms)
    return draw_overlay(
        histograms, title="", x_title=x_title, y_title=y_title,
        annotations=annotations, curve_styles=styles, curve_line_styles=line_styles,
        legend_bounds=(0.66, legend_bottom, 0.93, legend_top),
        x_range=x_range, y_range=y_range, reference_y=reference_y,
        grid=grid, headroom=1.6, style=FIGURE_STYLE, mirror_root_style=True,
        output=output_dir / f"{output_tag}.pdf", save_png=save_png,
        canvas_name=output_tag,
    )


def draw_dijet_comparison(histograms, nominal, *, output_tag, y_title, ratio_y_title,
                          ratio_range, ratio_option, styles, fb=False, y_range=None,
                          **figure_options):
    """Draw a distribution and its nominal comparisons as two figures."""
    numerators = [label for label in histograms if label != nominal]
    if not numerators:
        raise ValueError("At least one non-nominal curve is required")
    ratios = {
        label: ratio_to_nominal(histograms[label], histograms[nominal],
                                name=f"ratio_{index}", option=ratio_option)
        for index, label in enumerate(numerators)
    }

    x_title = "|#eta_{CM}^{dijet}|" if fb else "#eta_{CM}^{dijet}"
    # ROOT's axis selection also affects Integral() and histogram statistics.
    display_histograms = histograms if fb else _display_copies(histograms)
    display_ratios = ratios if fb else _display_copies(ratios)
    distribution = draw_dijet_figure(
        display_histograms, output_tag=f"{output_tag}_distributions",
        x_title=x_title, y_title=y_title, styles=styles, y_range=y_range, **figure_options,
    )
    ratio_labels = {
        label: f"({label} F/B) / ({nominal} F/B)" if fb else f"{label} / {nominal}"
        for label in ratios
    }
    ratio_canvas = draw_dijet_figure(
        {ratio_labels[label]: ratio for label, ratio in display_ratios.items()},
        output_tag=f"{output_tag}_ratios_to_gen",
        x_title=x_title, y_title=ratio_y_title, reference_y=1.0,
        styles={display: styles[label] for label, display in ratio_labels.items()},
        y_range=ratio_range, **figure_options,
    )
    return distribution, ratio_canvas, ratios


def save_dijet_grid(panels, *, output, cms_header_left, cms_header_right, save_png=False):
    """Save a 4×4 grid of existing closure canvases in row-major order."""
    if len(panels) != 16:
        raise ValueError("A dijet grid requires exactly 16 panels")

    canvas = ROOT.TCanvas(f"c_dijet_grid_{uuid4().hex}", "", 2400, 2400)
    canvas.SetCanvasSize(2400, 2400)
    canvas._mirror_root_style = True
    canvas._shared_axes_grid = True
    left, bottom, span = 0.085, 0.075, 0.895
    step = span / 4
    canvas.SetLeftMargin(left)
    canvas.SetRightMargin(1 - left - span)
    canvas._grid_pads = []
    x_range = min(p.GetUxmin() for p in panels), max(p.GetUxmax() for p in panels)
    y_range = min(p.GetUymin() for p in panels), max(p.GetUymax() for p in panels)
    legend = next(item for item in panels[1].GetListOfPrimitives()
                  if item.InheritsFrom("TLegend"))
    legend_rows = legend.GetListOfPrimitives().GetSize()
    text_top, line_step = 0.885, 0.10
    for index, source in enumerate(panels):
        row, column = divmod(index, 4)
        frame_left = left + column * step
        frame_bottom = bottom + (3 - row) * step
        pad_left = 0 if column == 0 else frame_left
        pad_bottom = 0 if row == 3 else frame_bottom
        canvas.cd()
        pad = ROOT.TPad(f"{canvas.GetName()}_{index}", "", pad_left, pad_bottom,
                        frame_left + step, frame_bottom + step)
        pad.Draw()
        pad.cd()
        source.DrawClonePad()
        pad.SetLeftMargin((frame_left - pad_left) / (frame_left + step - pad_left))
        pad.SetBottomMargin((frame_bottom - pad_bottom) / (frame_bottom + step - pad_bottom))
        pad.SetRightMargin(0)
        pad.SetTopMargin(0)
        left_margin, bottom_margin = pad.GetLeftMargin(), pad.GetBottomMargin()

        def map_x(x):
            return left_margin + x * (1 - left_margin)

        def map_y(y):
            return bottom_margin + y * (1 - bottom_margin)

        header_titles = {item.GetTitle() for item in source._cms_header_objects}
        first_annotation = True
        for item in list(pad.GetListOfPrimitives()):
            if item.InheritsFrom("TH1"):
                x_axis, y_axis = item.GetXaxis(), item.GetYaxis()
                if index == 0:
                    canvas._grid_axis_titles = x_axis.GetTitle(), y_axis.GetTitle()
                x_axis.SetRangeUser(*x_range)
                item.SetMinimum(y_range[0])
                item.SetMaximum(y_range[1])
                for axis, show_labels in ((x_axis, row == 3), (y_axis, column == 0)):
                    axis.SetTitle("")
                    axis.SetLabelSize(0.062 if show_labels else 0)
                    # Omit the endpoint labels where adjacent panels meet.
                    axis.ChangeLabel(1, -1, 0)
                    axis.ChangeLabel(-1, -1, 0)
            elif item.InheritsFrom("TLatex"):
                if item.GetTitle() in header_titles:
                    pad.GetListOfPrimitives().Remove(item)
                    continue
                title = item.GetTitle()
                item.SetTextFont(43)
                item.SetTextSize(36)
                item.SetX(map_x(0.045))
                item.SetY(map_y(text_top))
                if "p_{T}^{ave}" in title:
                    item.SetX(map_x(0.5))
                    item.SetY(map_y(0.985))
                    item.SetTextAlign(23)
                elif first_annotation and index == 0:
                    pass
                elif index == 4 and title.startswith(("|#eta_", "|#Delta#phi|")):
                    if title.startswith("|#Delta#phi|"):
                        item.SetY(map_y(text_top - line_step))
                elif index == 5 and "p_{T}^{Lead}" in title:
                    leading, subleading = title.split(", ", 1)
                    item.SetTitle(leading)
                    item.DrawLatex(map_x(0.045), map_y(text_top - line_step), subleading)
                else:
                    pad.GetListOfPrimitives().Remove(item)
                first_annotation = False
            elif item.InheritsFrom("TLegend"):
                if index == 1:
                    item.SetX1NDC(map_x(0.035))
                    item.SetX2NDC(map_x(0.97))
                    item.SetY1NDC(map_y(text_top - line_step * legend_rows))
                    item.SetY2NDC(map_y(text_top))
                    item.SetTextSize(0.067)
                else:
                    pad.GetListOfPrimitives().Remove(item)
        pad.Modified()
        canvas._grid_pads.append(pad)

    canvas.cd()
    header_style = replace(FIGURE_STYLE, cms_header_y=0.997, cms_header_text_size=0.023)
    draw_cms_header(canvas, cms_header_left, cms_header_right, style=header_style)
    x_title, y_title = canvas._grid_axis_titles
    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextAlign(22)
    label.SetTextFont(42)
    label.SetTextSize(0.026)
    label.DrawLatex(left + span / 2, 0.025, x_title)
    label.SetTextAngle(90)
    label.DrawLatex(0.025, bottom + span / 2, y_title)
    canvas.Modified()
    canvas.Update()
    save_canvas(canvas, output, save_png=save_png)
    return canvas
