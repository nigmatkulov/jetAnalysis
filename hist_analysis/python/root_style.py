"""ROOT styling and PDF/PNG canvas-output helpers shared by notebooks."""

from __future__ import annotations

from pathlib import Path

import ROOT


# ROOT equivalents of the ordering used by plotMcClosures.C::set1DStyle.
COLORS = (
    ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kMagenta,
    ROOT.kOrange + 7, ROOT.kGreen + 2, ROOT.kAzure + 1,
    ROOT.kPink + 1, ROOT.kCyan + 1, ROOT.kTeal + 1,
    ROOT.kGray + 1, ROOT.kSpring + 5, ROOT.kViolet + 1,
)
MARKERS = (20, 21, 20, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32)


def set_pad_style(pad, *, grid_x: bool = True, grid_y: bool = True) -> None:
    """Apply ``plotMcClosures.C::setPadStyle`` with configurable grids."""

    pad.SetTopMargin(0.10)
    pad.SetBottomMargin(0.15)
    pad.SetRightMargin(0.05)
    pad.SetLeftMargin(0.15)
    pad.SetGridx(grid_x)
    pad.SetGridy(grid_y)


def set_1d_style(histogram, style_index: int = 0) -> None:
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
    histogram.GetYaxis().SetTitleSize(0.05)
    histogram.GetYaxis().SetLabelSize(0.05)
    histogram.GetXaxis().SetTitleSize(0.05)
    histogram.GetXaxis().SetLabelSize(0.05)
    histogram.GetXaxis().SetTitleOffset(1.2)
    histogram.GetXaxis().SetNdivisions(208)
    histogram.GetYaxis().SetNdivisions(208)
    histogram.GetYaxis().SetTitleOffset(1.5)


def set_2d_style(histogram) -> None:
    """Apply ``plotMcClosures.C::set2DStyle``."""

    histogram.GetYaxis().SetTitleSize(0.06)
    histogram.GetYaxis().SetLabelSize(0.06)
    histogram.GetXaxis().SetTitleSize(0.06)
    histogram.GetXaxis().SetLabelSize(0.06)
    histogram.GetXaxis().SetNdivisions(205)
    histogram.GetYaxis().SetNdivisions(205)
    histogram.GetYaxis().SetTitleOffset(1.0)


def set_palette_style(histogram, *, x1: float = 0.89, width: float = 0.02):
    """Narrow and position the COLZ palette inside the canvas right margin."""

    if width <= 0.0:
        raise ValueError("Palette width must be positive")
    palette = histogram.GetListOfFunctions().FindObject("palette")
    if palette is None:
        raise RuntimeError("ROOT palette is unavailable; update the canvas after drawing COLZ")
    palette.SetX1NDC(x1)
    palette.SetX2NDC(x1 + width)
    palette.SetY1NDC(0.15)
    palette.SetY2NDC(0.90)
    palette.SetLabelSize(0.035)
    return palette


def set_legend_style(legend) -> None:
    legend.SetFillColor(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)


def save_canvas(canvas, pdf_path: str | Path | None, *, save_png: bool = False) -> None:
    """Save PDF and, when requested, a same-named PNG."""

    if pdf_path is None:
        return
    path = Path(pdf_path)
    if path.suffix.lower() != ".pdf":
        raise ValueError(f"Canvas output must have a .pdf suffix: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(path))
    if save_png:
        canvas.SaveAs(str(path.with_suffix(".png")))
