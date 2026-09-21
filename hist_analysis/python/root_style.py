"""ROOT styling and PDF/PNG canvas-output helpers shared by notebooks."""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
from typing import Iterable

import ROOT

from hist_analysis.python.plot_labels import (
    CMS_HEADER_CHOICES, CMS_SIMULATION_HEADER, PPB_HEADER_RIGHT,
)


# ROOT equivalents of the ordering used by plotMcClosures.C::set1DStyle.
COLORS = (
    ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kMagenta,
    ROOT.kOrange + 7, ROOT.kGreen + 2, ROOT.kAzure + 1,
    ROOT.kPink + 1, ROOT.kCyan + 1, ROOT.kTeal + 1,
    ROOT.kGray + 1, ROOT.kSpring + 5, ROOT.kViolet + 1,
    ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kOrange + 7,
)
MARKERS = (20, 21, 20, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 43, 44)

# Semantic colors for gen/reco unfolding diagnostics. Keep these mappings in
# one place so a physics role has the same appearance in every notebook panel.
UNFOLDING_COLORS = {
    "gen": ROOT.kBlue,
    "ref": ROOT.kMagenta + 1,
    "reco": ROOT.kRed,
    "miss": ROOT.kOrange + 7,
    "fake": ROOT.kGreen + 2,
    "unfolded": ROOT.kBlack,
}

UNFOLDING_MARKERS = {
    "gen": 21,
    "ref": 23,
    "reco": 20,
    "miss": 24,
    "fake": 25,
    "unfolded": 22,
}


@dataclass(frozen=True)
class PlotStyle:
    """Shared geometry and typography for ROOT figures."""

    canvas_width: int = 800
    canvas_height: int = 800
    side_by_side_canvas_width: int = 1400
    side_by_side_canvas_height: int = 700
    font: int = 42
    annotation_x: float = 0.22
    annotation_y: float = 0.87
    annotation_text_size: float = 0.035
    annotation_line_spacing: float = 0.052
    legend_text_size: float = 0.035
    top_margin: float = 0.10
    bottom_margin: float = 0.15
    right_margin: float = 0.05
    left_margin: float = 0.15
    axis_1d_title_size: float = 0.05
    axis_1d_label_size: float = 0.05
    axis_1d_divisions: int = 208
    axis_1d_x_title_offset: float = 1.20
    axis_1d_y_title_offset: float = 1.50
    axis_2d_title_size: float = 0.06
    axis_2d_label_size: float = 0.06
    axis_2d_divisions: int = 205
    axis_2d_y_title_offset: float = 1.00
    single_axis_title_size: float = 0.045
    single_axis_label_size: float = 0.040
    single_axis_divisions: int = 208
    single_x_title_offset: float = 1.20
    single_y_title_offset: float = 1.65
    single_panel_left_margin: float = 0.20
    single_panel_bottom_margin: float = 0.17
    ratio_fraction: float = 0.30
    ratio_left_margin: float = 0.20
    ratio_top_bottom_margin: float = 0.02
    ratio_bottom_top_margin: float = 0.03
    ratio_bottom_margin: float = 0.35
    ratio_top_axis_title_size: float = 0.050
    ratio_top_axis_label_size: float = 0.045
    ratio_top_y_title_offset: float = 1.55
    ratio_axis_title_size: float = 0.115
    ratio_axis_label_size: float = 0.105
    ratio_x_title_offset: float = 1.05
    ratio_y_title_offset: float = 0.55
    ratio_axis_divisions: int = 505
    palette_x1: float = 0.89
    palette_width: float = 0.02
    palette_y1: float = 0.15
    palette_y2: float = 0.90
    palette_label_size: float = 0.035
    palette_right_margin: float = 0.12
    cms_header_text_size: float = 0.040
    cms_header_y: float = 0.975


DEFAULT_PLOT_STYLE = PlotStyle()

_DEFAULT_CMS_HEADER_LEFT: str | None = CMS_SIMULATION_HEADER
_DEFAULT_CMS_HEADER_RIGHT: str | None = PPB_HEADER_RIGHT
_SAVE_MPLHEP = os.environ.get("HIST_ANALYSIS_SAVE_MPLHEP", "0") == "1"
_SAVE_MPLHEP_PNG = os.environ.get("HIST_ANALYSIS_SAVE_MPLHEP_PNG", "0") == "1"


def configure_cms_header(left: str | None, right: str | None) -> None:
    """Configure the header automatically added when canvases are saved."""

    if left is not None and left not in CMS_HEADER_CHOICES:
        raise ValueError(
            f"left CMS header must be one of {CMS_HEADER_CHOICES}, got {left!r}"
        )
    if right is not None and (not isinstance(right, str) or not right.strip()):
        raise ValueError("right CMS header must be None or a non-empty string")
    global _DEFAULT_CMS_HEADER_LEFT, _DEFAULT_CMS_HEADER_RIGHT
    _DEFAULT_CMS_HEADER_LEFT = left
    _DEFAULT_CMS_HEADER_RIGHT = right


def configure_mplhep_export(enabled: bool, *, save_png: bool = False) -> None:
    """Enable optional Python re-rendering for every canvas saved thereafter."""

    if not isinstance(enabled, bool) or not isinstance(save_png, bool):
        raise TypeError("enabled and save_png must be bool")
    global _SAVE_MPLHEP, _SAVE_MPLHEP_PNG
    _SAVE_MPLHEP = enabled
    _SAVE_MPLHEP_PNG = save_png


def draw_cms_header(pad, left: str | None, right: str | None, *,
                    style: PlotStyle = DEFAULT_PLOT_STYLE):
    """Draw an optional CMS header above the ROOT plotting frame."""

    if not left and not right:
        return []
    left_labels = {
        "CMS Simulation": "#bf{CMS} #it{Simulation}",
        "CMS Preliminary": "#bf{CMS} #it{Preliminary}",
        "CMS": "#bf{CMS}",
    }
    if left not in left_labels:
        raise ValueError(
            "left CMS header must be CMS Simulation, CMS Preliminary, or CMS"
        )
    pad.cd()
    objects = []
    left_text = ROOT.TLatex(
        pad.GetLeftMargin(), style.cms_header_y, left_labels[left],
    )
    left_text.SetNDC(True)
    left_text.SetTextFont(style.font)
    left_text.SetTextSize(style.cms_header_text_size)
    left_text.SetTextAlign(13)
    left_text.Draw()
    objects.append(left_text)
    if right:
        right_text = ROOT.TLatex(
            1.0 - pad.GetRightMargin(), style.cms_header_y,
            right.replace("$", ""),
        )
        right_text.SetNDC(True)
        right_text.SetTextFont(style.font)
        right_text.SetTextSize(style.cms_header_text_size)
        right_text.SetTextAlign(33)
        right_text.Draw()
        objects.append(right_text)
    pad._cms_header_objects = objects
    return objects


def apply_curve_style(histogram, color: int, marker: int, *,
                      style: PlotStyle = DEFAULT_PLOT_STYLE) -> None:
    """Apply an explicit ROOT color and marker while retaining shared axes."""

    set_1d_style(histogram, style=style)
    histogram.SetLineColor(color)
    histogram.SetMarkerColor(color)
    histogram.SetMarkerStyle(marker)


def set_pad_style(pad, *, grid_x: bool = True, grid_y: bool = True,
                  style: PlotStyle = DEFAULT_PLOT_STYLE) -> None:
    """Apply ``plotMcClosures.C::setPadStyle`` with configurable grids."""

    pad.SetTopMargin(style.top_margin)
    pad.SetBottomMargin(style.bottom_margin)
    pad.SetRightMargin(style.right_margin)
    pad.SetLeftMargin(style.left_margin)
    pad.SetGridx(grid_x)
    pad.SetGridy(grid_y)


def set_1d_style(histogram, style_index: int = 0, *,
                 style: PlotStyle = DEFAULT_PLOT_STYLE) -> None:
    """Apply the 1D marker, line, and axis style from ``set1DStyle``."""

    if style_index < 0:
        raise ValueError("style_index must be non-negative")
    if style_index < len(COLORS):
        color = COLORS[style_index]
        marker = MARKERS[style_index]
    else:
        color = ROOT.kMagenta
        marker = 45

    histogram.SetLineColor(color)
    histogram.SetMarkerColor(color)
    histogram.SetMarkerStyle(marker)
    histogram.SetLineWidth(2)
    histogram.SetMarkerSize(1.3)
    for axis in (histogram.GetXaxis(), histogram.GetYaxis()):
        axis.SetTitleFont(style.font)
        axis.SetLabelFont(style.font)
        axis.SetTitleSize(style.axis_1d_title_size)
        axis.SetLabelSize(style.axis_1d_label_size)
        axis.SetNdivisions(style.axis_1d_divisions)
    histogram.GetXaxis().SetTitleOffset(style.axis_1d_x_title_offset)
    histogram.GetYaxis().SetTitleOffset(style.axis_1d_y_title_offset)


def set_unfolding_1d_style(histogram, role: str, *,
                           style: PlotStyle = DEFAULT_PLOT_STYLE) -> None:
    """Style a 1D unfolding histogram according to its physics role."""

    if role not in UNFOLDING_COLORS:
        raise ValueError(
            f"Unknown unfolding role {role!r}; expected one of "
            f"{tuple(UNFOLDING_COLORS)}"
        )
    set_1d_style(histogram, style=style)
    color = UNFOLDING_COLORS[role]
    histogram.SetLineColor(color)
    histogram.SetMarkerColor(color)
    histogram.SetMarkerStyle(UNFOLDING_MARKERS[role])


def set_2d_style(histogram, *, style: PlotStyle = DEFAULT_PLOT_STYLE) -> None:
    """Apply ``plotMcClosures.C::set2DStyle``."""

    for axis in (histogram.GetXaxis(), histogram.GetYaxis()):
        axis.SetTitleFont(style.font)
        axis.SetLabelFont(style.font)
        axis.SetTitleSize(style.axis_2d_title_size)
        axis.SetLabelSize(style.axis_2d_label_size)
        axis.SetNdivisions(style.axis_2d_divisions)
    histogram.GetYaxis().SetTitleOffset(style.axis_2d_y_title_offset)


def style_single_panel_axes(histogram, *, style: PlotStyle = DEFAULT_PLOT_STYLE) -> None:
    """Apply common font, sizing, divisions, and offsets to TH1 axes."""

    for axis in (histogram.GetXaxis(), histogram.GetYaxis()):
        axis.SetTitleFont(style.font)
        axis.SetLabelFont(style.font)
        axis.SetTitleSize(style.single_axis_title_size)
        axis.SetLabelSize(style.single_axis_label_size)
        axis.SetNdivisions(style.single_axis_divisions)
        axis.CenterTitle(True)
    histogram.GetXaxis().SetTitleOffset(style.single_x_title_offset)
    histogram.GetYaxis().SetTitleOffset(style.single_y_title_offset)


def draw_text_block(pad, lines: Iterable[str], *, x: float | None = None,
                    y: float | None = None, text_size: float | None = None,
                    line_spacing: float | None = None,
                    style: PlotStyle = DEFAULT_PLOT_STYLE):
    """Draw non-empty TLatex lines with spacing safe for ROOT scripts."""

    rendered_lines = [line for line in lines if line]
    if not rendered_lines:
        return []
    size = style.annotation_text_size if text_size is None else text_size
    if size <= 0.0:
        raise ValueError("text_size must be positive")
    spacing = max(
        style.annotation_line_spacing if line_spacing is None else line_spacing,
        1.45 * size,
    )
    x_position = style.annotation_x if x is None else x
    y_position = style.annotation_y if y is None else y
    if y_position - spacing * (len(rendered_lines) - 1) <= pad.GetBottomMargin():
        raise ValueError("Text block does not fit above the pad bottom margin")

    pad.cd()
    labels = []
    for line_index, line in enumerate(rendered_lines):
        label = ROOT.TLatex(x_position, y_position - spacing * line_index, line)
        label.SetNDC(True)
        label.SetTextFont(style.font)
        label.SetTextSize(size)
        label.SetTextAlign(13)
        label.Draw()
        labels.append(label)
    return labels


def set_palette_style(histogram, *, x1: float | None = None,
                      width: float | None = None,
                      style: PlotStyle = DEFAULT_PLOT_STYLE):
    """Narrow and position the COLZ palette inside the canvas right margin."""

    palette_x1 = style.palette_x1 if x1 is None else x1
    palette_width = style.palette_width if width is None else width
    if palette_width <= 0.0:
        raise ValueError("Palette width must be positive")
    palette = histogram.GetListOfFunctions().FindObject("palette")
    if palette is None:
        raise RuntimeError("ROOT palette is unavailable; update the canvas after drawing COLZ")
    palette.SetX1NDC(palette_x1)
    palette.SetX2NDC(palette_x1 + palette_width)
    palette.SetY1NDC(style.palette_y1)
    palette.SetY2NDC(style.palette_y2)
    palette.SetLabelFont(style.font)
    palette.SetLabelSize(style.palette_label_size)
    return palette


def set_legend_style(legend, *, style: PlotStyle = DEFAULT_PLOT_STYLE) -> None:
    legend.SetFillColor(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(style.font)
    legend.SetTextSize(style.legend_text_size)


def _clear_canvas_titles(pad) -> None:
    """Remove ROOT title boxes that would overlap a configured CMS header."""

    primitives = pad.GetListOfPrimitives()
    if primitives is None:
        return
    for primitive in tuple(primitives):
        if primitive.InheritsFrom("TH1"):
            primitive.SetTitle("")
        elif primitive.InheritsFrom("TPad"):
            _clear_canvas_titles(primitive)
    title_box = primitives.FindObject("title")
    if title_box is not None:
        primitives.Remove(title_box)


def save_canvas(canvas, pdf_path: str | Path | None, *, save_png: bool = False,
                mplhep_x_range: tuple[float, float] | None = None,
                mplhep_y_range: tuple[float, float] | None = None) -> None:
    """Save PDF and, when requested, a same-named PNG."""

    if pdf_path is None:
        return
    if _DEFAULT_CMS_HEADER_LEFT or _DEFAULT_CMS_HEADER_RIGHT:
        _clear_canvas_titles(canvas)
    if (
        (_DEFAULT_CMS_HEADER_LEFT or _DEFAULT_CMS_HEADER_RIGHT)
        and not getattr(canvas, "_cms_header_objects", None)
    ):
        draw_cms_header(
            canvas, _DEFAULT_CMS_HEADER_LEFT, _DEFAULT_CMS_HEADER_RIGHT,
        )
        canvas.Modified()
        canvas.Update()
    path = Path(pdf_path)
    if path.suffix.lower() != ".pdf":
        raise ValueError(f"Canvas output must have a .pdf suffix: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(path))
    if save_png:
        canvas.SaveAs(str(path.with_suffix(".png")))
    if _SAVE_MPLHEP:
        from hist_analysis.python.mplhep_plotting import save_canvas as save_mplhep_canvas
        mplhep_path = path.parent / "mplhep" / path.name
        save_mplhep_canvas(
            canvas, output=mplhep_path,
            cms_header_left=_DEFAULT_CMS_HEADER_LEFT or CMS_SIMULATION_HEADER,
            cms_header_right=_DEFAULT_CMS_HEADER_RIGHT or PPB_HEADER_RIGHT,
            save_png=_SAVE_MPLHEP_PNG,
            x_range=mplhep_x_range, y_range=mplhep_y_range,
        )
