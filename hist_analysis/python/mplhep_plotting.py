"""Optional mplhep rendering of detached ROOT histograms."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import re
from typing import Iterable, Mapping

from hist_analysis.python.plot_labels import (
    CMS_HEADER_CHOICES, CMS_SIMULATION_HEADER, PPB_MPLHEP_HEADER_RIGHT,
)


@dataclass(frozen=True)
class MplhepStyle:
    """User-facing mplhep figure geometry and typography."""

    figure_size: float = 8.0
    axis_title_size: float = 22.0
    ratio_axis_title_size: float = 17.0
    tick_label_size: float = 17.0
    ratio_tick_label_size: float = 15.0
    annotation_text_size: float = 16.0
    annotation_y: float = 0.955
    legend_text_size: float = 17.0
    cms_label_size: float = 20.0
    marker_size: float = 10.0
    ratio_marker_size: float = 9.0
    capsize: float = 2.0
    colorbar_tick_size: float = 14.0

    def __post_init__(self) -> None:
        for field_name, value in vars(self).items():
            if field_name == "annotation_y" and not 0.0 < value < 1.0:
                raise ValueError(f"annotation_y must lie between 0 and 1, got {value}")
            if field_name == "capsize" and value < 0:
                raise ValueError(f"capsize must be non-negative, got {value}")
            if field_name not in {"annotation_y", "capsize"} and value <= 0:
                raise ValueError(f"{field_name} must be positive, got {value}")


DEFAULT_MPLHEP_STYLE = MplhepStyle()
CMS_HEADER_LEFT_CHOICES = CMS_HEADER_CHOICES
MPLHEP_CURVE_STYLES = (
    ("tab:red", "o", True),
    ("tab:blue", "o", False),
    ("tab:green", "^", True),
    ("tab:orange", "^", False),
    ("tab:purple", "s", True),
    ("tab:brown", "s", False),
    ("tab:pink", "D", True),
    ("tab:gray", "D", False),
)


def _curve_style(label: str, index: int):
    """Return semantic or deterministic color, marker, and fill settings."""

    normalized = label.upper().replace(" ", "")
    if normalized == "DATA":
        return label, "black", "s", "black"
    if "CUETP8M1" in normalized or "CUEP8M1" in normalized:
        display_label = (
            label
            if "/" in label or any(token in normalized for token in ("OFFICIAL", "STANDALONE"))
            else "CUEP8M1"
        )
        marker_face = "tab:blue" if "OFFICIAL" in normalized else "none"
        return display_label, "tab:blue", "o", marker_face
    if "CP5" in normalized:
        marker_face = "none" if "STANDALONE" in normalized else "tab:red"
        return label, "tab:red", "o", marker_face
    if "EPPS21" in normalized:
        return label, "tab:green", "D", "tab:green"
    color, marker, filled = MPLHEP_CURVE_STYLES[index % len(MPLHEP_CURVE_STYLES)]
    return label, color, marker, color if filled else "none"


def _mplhep():
    """Import optional plotting dependencies only when export is requested."""

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mplhep as hep

    return plt, hep


def _draw_cms_header(hep, axis, *, left: str, right: str,
                     style: MplhepStyle) -> None:
    """Draw one of the supported CMS left labels with a custom right label."""

    if left not in CMS_HEADER_LEFT_CHOICES:
        raise ValueError(
            f"cms_header_left must be one of {CMS_HEADER_LEFT_CHOICES}, got {left!r}"
        )
    kwargs = {
        "rlabel": right,
        "loc": 0,
        "fontsize": style.cms_label_size,
        "ax": axis,
    }
    if left == "CMS Simulation":
        hep.cms.label(data=False, **kwargs)
    elif left == "CMS Preliminary":
        hep.cms.label("Preliminary", data=True, **kwargs)
    else:
        hep.cms.label(data=True, **kwargs)


def _th1_arrays(histogram, *, mask_zero: bool = True):
    axis = histogram.GetXaxis()
    edges = [axis.GetBinLowEdge(index) for index in range(1, histogram.GetNbinsX() + 2)]
    values = [histogram.GetBinContent(index) for index in range(1, histogram.GetNbinsX() + 1)]
    errors = [histogram.GetBinError(index) for index in range(1, histogram.GetNbinsX() + 1)]
    if mask_zero:
        for index, value in enumerate(values):
            if value == 0.0:
                values[index] = float("nan")
                errors[index] = float("nan")
    return values, edges, errors


def _th2_arrays(histogram):
    x_axis = histogram.GetXaxis()
    y_axis = histogram.GetYaxis()
    x_edges = [x_axis.GetBinLowEdge(index) for index in range(1, histogram.GetNbinsX() + 2)]
    y_edges = [y_axis.GetBinLowEdge(index) for index in range(1, histogram.GetNbinsY() + 2)]
    values = [
        [histogram.GetBinContent(x_index, y_index) for y_index in range(1, histogram.GetNbinsY() + 1)]
        for x_index in range(1, histogram.GetNbinsX() + 1)
    ]
    return values, x_edges, y_edges


def _save(fig, output: str | Path, *, save_png: bool,
          style: MplhepStyle) -> None:
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.set_size_inches(style.figure_size, style.figure_size, forward=True)
    fig.savefig(path)
    if save_png:
        fig.savefig(path.with_suffix(".png"), dpi=200)


def _root_label(text: str) -> str:
    """Translate the ROOT TLatex subset used by this project to math text."""

    if not text:
        return ""
    translated = text.replace("#geq", r"\geq").replace("#leq", r"\leq")
    translated = translated.replace("#Delta", r"\Delta").replace("#phi", r"\phi")
    translated = translated.replace("#eta", r"\eta").replace("#hat", r"\hat")
    translated = translated.replace("#sqrt", r"\sqrt").replace("#pi", r"\pi")
    translated = translated.replace("#bf", r"\mathbf").replace("#it", r"\mathit")
    if "\\" in translated or "_{" in translated or "^{" in translated:
        translated = re.sub(r"\s+", r"\\ ", translated)
        return f"${translated}$"
    return translated


def _leaf_pads(container):
    pads = [item for item in container.GetListOfPrimitives() if item.InheritsFrom("TPad")]
    if not pads:
        return [container]
    leaves = []
    for pad in pads:
        leaves.extend(_leaf_pads(pad))
    return leaves


def _pad_histograms(pad):
    return [
        item for item in pad.GetListOfPrimitives()
        if (
            item.InheritsFrom("TH1")
            and not item.InheritsFrom("THStack")
            and not item.GetName().startswith("hframe")
        )
    ]


def _pad_graphs(pad):
    return [
        item for item in pad.GetListOfPrimitives()
        if item.InheritsFrom("TGraphAsymmErrors")
    ]


def _pad_frame(pad):
    for item in pad.GetListOfPrimitives():
        if item.InheritsFrom("TH1") and item.GetName().startswith("hframe"):
            return item
    return None


def _graph_arrays(graph):
    count = graph.GetN()
    return {
        "x": [float(graph.GetPointX(index)) for index in range(count)],
        "y": [float(graph.GetPointY(index)) for index in range(count)],
        "xlow": [float(graph.GetErrorXlow(index)) for index in range(count)],
        "xhigh": [float(graph.GetErrorXhigh(index)) for index in range(count)],
        "ylow": [float(graph.GetErrorYlow(index)) for index in range(count)],
        "yhigh": [float(graph.GetErrorYhigh(index)) for index in range(count)],
    }


def _root_color(color_index: int) -> str:
    import ROOT

    color = ROOT.gROOT.GetColor(color_index)
    return color.AsHexString() if color else "black"


def _graph_marker(marker_style: int) -> tuple[str, bool]:
    markers = {20: "o", 21: "s", 22: "^", 23: "v", 24: "o", 25: "s", 26: "^", 33: "D"}
    return markers.get(marker_style, "D"), marker_style not in {24, 25, 26, 32}


def _root_pad_range(pad, axis: str) -> tuple[float, float]:
    """Return the final displayed ROOT pad range in linear coordinates."""

    if axis == "x":
        low, high, logarithmic = pad.GetUxmin(), pad.GetUxmax(), pad.GetLogx()
    elif axis == "y":
        low, high, logarithmic = pad.GetUymin(), pad.GetUymax(), pad.GetLogy()
    else:
        raise ValueError("axis must be 'x' or 'y'")
    if logarithmic:
        low, high = 10.0 ** low, 10.0 ** high
    return float(low), float(high)


def _hist_axis_range(histogram, axis: str) -> tuple[float, float]:
    root_axis = histogram.GetXaxis() if axis == "x" else histogram.GetYaxis()
    first, last = root_axis.GetFirst(), root_axis.GetLast()
    return root_axis.GetBinLowEdge(first), root_axis.GetBinUpEdge(last)


def _legend_labels(pad):
    labels = {}
    for primitive in pad.GetListOfPrimitives():
        if not primitive.InheritsFrom("TLegend"):
            continue
        for entry in primitive.GetListOfPrimitives():
            target = entry.GetObject()
            if target:
                labels[target.GetName()] = entry.GetLabel()
    return labels


def _pad_ptave_annotations(pad):
    """Return NDC pTave-selection labels drawn on a ROOT pad."""

    annotations = []
    for primitive in pad.GetListOfPrimitives():
        if not primitive.InheritsFrom("TLatex"):
            continue
        label = primitive.GetTitle()
        if not any(token in label for token in (
                "p_{T}^{ave}", "|#eta_{CM}^{jet}|", "anti-k_{T}", "#it{anti}", "p_{T,1}", "#Delta#phi",
                "Minimum bias", "Jet60", "Jet80", "Jet100")):
            continue
        annotations.append((
            float(primitive.GetX()), float(primitive.GetY()), _root_label(label),
        ))
    return annotations


def save_canvas(
    canvas,
    *,
    output: str | Path,
    cms_header_left: str = CMS_SIMULATION_HEADER,
    cms_header_right: str = PPB_MPLHEP_HEADER_RIGHT,
    save_png: bool = False,
    style: MplhepStyle = DEFAULT_MPLHEP_STYLE,
    x_range: tuple[float, float] | None = None,
    y_range: tuple[float, float] | None = None,
) -> None:
    """Re-render a ROOT canvas, including nested pads, with mplhep."""

    plt, hep = _mplhep()
    hep.style.use("CMS")
    pads = [
        pad for pad in _leaf_pads(canvas)
        if _pad_histograms(pad) or _pad_graphs(pad)
    ]
    if not pads:
        raise ValueError(f"Canvas {canvas.GetName()} contains no drawable histograms")
    count = len(pads)
    columns = 1 if count <= 2 else math.ceil(math.sqrt(count))
    rows = math.ceil(count / columns)
    panel_size = 4.0
    figure_width = style.figure_size if count <= 2 else max(style.figure_size, panel_size * columns)
    figure_height = style.figure_size if count <= 2 else max(style.figure_size, panel_size * rows)
    font_scale = 1.0 if count <= 2 else max(0.45, 1.0 / math.sqrt(count / 2.0))
    comparison_canvas = getattr(canvas, '_data_mc_plot', False) or any(p.InheritsFrom('TLatex') and p.GetTitle().startswith('anti-k_{T}')
                            for pad in pads for p in pad.GetListOfPrimitives())
    ratio_layout = comparison_canvas and count == 2 and len({round(p.GetXlowNDC(), 4) for p in pads}) == 1
    fig, axes = plt.subplots(
        rows, columns, figsize=(figure_width, figure_height),
        gridspec_kw={'height_ratios': (3, 1)} if ratio_layout else None,
        squeeze=False, layout="constrained" if count <= 2 else None,
    )
    flat_axes = list(axes.flat)
    for axis, pad in zip(flat_axes, pads):
        histograms = _pad_histograms(pad)
        graphs = _pad_graphs(pad)
        legend_labels = _legend_labels(pad)
        epps21_panel = any(g.GetName().startswith("g_epps21_") for g in graphs)
        comparison_panel = comparison_canvas or epps21_panel
        if comparison_panel:
            histograms = list({h.GetName(): h for h in histograms}.values())
        two_dimensional = [hist for hist in histograms if hist.InheritsFrom("TH2")]
        if two_dimensional:
            histogram = two_dimensional[-1]
            values, x_edges, y_edges = _th2_arrays(histogram)
            artists = hep.hist2dplot(
                values, xbins=x_edges, ybins=y_edges, ax=axis, cbar=True,
                cbarsize="4%", cbarpad=0.04,
            )
            artists.cbar.ax.tick_params(labelsize=style.colorbar_tick_size)
        else:
            for graph in graphs:
                if graph.GetFillStyle() == 0 and graph.GetMarkerSize() > 0.0:
                    continue
                arrays = _graph_arrays(graph)
                normalized_name = graph.GetName().lower()
                if "pdf" in normalized_name:
                    color, hatch, alpha = "tab:green", None, 0.20
                elif "scale" in normalized_name:
                    color, hatch, alpha = "tab:orange", "////", 0.14
                else:
                    color = _root_color(graph.GetFillColor())
                    hatch = "////" if graph.GetFillStyle() >= 3000 else None
                    alpha = 0.28 if hatch is None else 0.12
                if graph.GetName().startswith('g_data_total_'):
                    from hist_analysis.python.data_mc_plotting import draw_mpl_systematic_brackets
                    draw_mpl_systematic_brackets(axis, arrays, _root_color(graph.GetLineColor()))
                elif graph.GetName().startswith("g_epps21_") and "_total_" in graph.GetName():
                    import ROOT
                    color = _root_color(graph.GetLineColor())
                    alpha = ROOT.gROOT.GetColor(graph.GetFillColor()).GetAlpha()
                    # Independent stored bin rectangles preserve gaps and outer half bins.
                    for i, (x, y) in enumerate(zip(arrays["x"], arrays["y"])):
                        bounds = [x-arrays["xlow"][i], x+arrays["xhigh"][i]]
                        axis.fill_between(bounds, y-arrays["ylow"][i], y+arrays["yhigh"][i],
                                          color=color, alpha=alpha, linewidth=0, zorder=1,
                                          label=legend_labels.get(graph.GetName()) if i == 0 else "_nolegend_")
                else:
                    axis.fill_between(
                        arrays["x"],
                        [value - error for value, error in zip(arrays["y"], arrays["ylow"])],
                        [value + error for value, error in zip(arrays["y"], arrays["yhigh"])],
                        color=color, alpha=alpha,
                        hatch=hatch, edgecolor=color, linewidth=0.8,
                        label=legend_labels.get(graph.GetName()),
                    )
            for index, histogram in enumerate(histograms):
                values, edges, errors = _th1_arrays(histogram)
                label = legend_labels.get(histogram.GetName(), histogram.GetName())
                display_label, color, marker, marker_face = _curve_style(label, index)
                if comparison_panel:
                    color = _root_color(histogram.GetMarkerColor())
                    marker, filled = _graph_marker(histogram.GetMarkerStyle())
                    marker_face = color if filled else "none"
                    display_label = legend_labels.get(histogram.GetName(), "_nolegend_")
                hep.histplot(
                    values, bins=edges, yerr=errors, xerr=True,
                    histtype="errorbar", ax=axis, label=display_label, color=color,
                    marker=marker, markerfacecolor=marker_face,
                    markeredgecolor=color, markeredgewidth=1.5,
                    markersize=style.marker_size * font_scale,
                    capsize=style.capsize,
                    zorder=(6 if comparison_panel and histogram.GetMarkerStyle() == 21
                            else 4 if comparison_panel else 2),
                )
            for graph in graphs:
                if graph.GetFillStyle() != 0 or graph.GetMarkerSize() <= 0.0:
                    continue
                arrays = _graph_arrays(graph)
                color = _root_color(graph.GetMarkerColor())
                marker, filled = _graph_marker(graph.GetMarkerStyle())
                axis.errorbar(
                    arrays["x"], arrays["y"],
                    xerr=(arrays["xlow"], arrays["xhigh"]),
                    yerr=(arrays["ylow"], arrays["yhigh"]),
                    linestyle="none", marker=marker, color=color,
                    markerfacecolor=color if filled else "none",
                    markeredgecolor=color,
                    markersize=style.marker_size * font_scale,
                    capsize=style.capsize,
                    label=("_nolegend_" if graph.GetName().startswith("g_epps21_central_")
                           else legend_labels.get(graph.GetName(), graph.GetName())),
                )
            if any(label for label in legend_labels.values()):
                handles, names = axis.get_legend_handles_labels()
                if getattr(canvas, '_marker_only_legend', comparison_panel):
                    from matplotlib.container import ErrorbarContainer
                    from matplotlib.lines import Line2D
                    # Apply to all notebook plot families, independently of theory content.
                    for i, handle in enumerate(handles):
                        if isinstance(handle, ErrorbarContainer):
                            point = handle.lines[0]
                            handles[i] = Line2D([], [], linestyle='none',
                                marker=point.get_marker(), markersize=point.get_markersize(),
                                markerfacecolor=point.get_markerfacecolor(),
                                markeredgecolor=point.get_markeredgecolor(),
                                markeredgewidth=point.get_markeredgewidth())
                if comparison_panel:
                    axis.legend(handles, names, fontsize=12 * font_scale, loc="upper right",
                                bbox_to_anchor=(0.98, 0.955), borderaxespad=0,
                                borderpad=0, frameon=False)
                else:
                    axis.legend(handles, names, fontsize=style.legend_text_size * font_scale)
        frame = _pad_frame(pad)
        reference = frame if frame is not None else histograms[-1]
        axis.set_xlabel(
            _root_label(reference.GetXaxis().GetTitle()),
            fontsize=(20 if comparison_panel else style.axis_title_size) * font_scale,
        )
        axis.set_ylabel(
            _root_label(reference.GetYaxis().GetTitle()),
            fontsize=(18 if comparison_panel else style.axis_title_size) * font_scale,
        )
        axis.tick_params(axis="both", labelsize=(16 if comparison_panel else style.tick_label_size) * font_scale)
        axis.set_xscale("log" if pad.GetLogx() else "linear")
        axis.set_yscale("log" if pad.GetLogy() else "linear")
        axis.set_xlim(*(x_range or _root_pad_range(pad, "x")))
        axis.set_ylim(*(y_range or _root_pad_range(pad, "y")))
        if getattr(pad, '_comparison_reference', False):
            axis.axhline(1., color='gray', linestyle='--', linewidth=1, zorder=2)
        if pad.GetGridx() or pad.GetGridy():
            axis.grid(True, alpha=0.25)
        else:
            axis.grid(False)
        for annotation_index, (x_position, y_position, label) in enumerate(_pad_ptave_annotations(pad)):
            axis.text(
                0.04 if comparison_panel else x_position,
                (0.955 - annotation_index * 0.045) if comparison_panel else y_position, label,
                transform=axis.transAxes, ha="left", va="top",
                fontsize=(12 if comparison_panel else style.annotation_text_size) * font_scale,
            )
    for axis in flat_axes[count:]:
        axis.remove()
    rendered_right = cms_header_right.replace("nb^{-1}", r"nb$^{-1}$")
    if count <= 2:
        _draw_cms_header(
            hep, flat_axes[0], left=cms_header_left,
            right=rendered_right, style=style,
        )
    else:
        left_label = {
            "CMS Simulation": r"$\bf{CMS}$ $\it{Simulation}$",
            "CMS Preliminary": r"$\bf{CMS}$ $\it{Preliminary}$",
            "CMS": r"$\bf{CMS}$",
        }[cms_header_left]
        fig.text(0.01, 0.995, left_label, ha="left", va="top",
                 fontsize=style.cms_label_size)
        fig.text(0.99, 0.995, rendered_right, ha="right", va="top",
                 fontsize=style.cms_label_size)
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.965))
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    if save_png:
        fig.savefig(path.with_suffix(".png"), dpi=200)
    plt.close(fig)


def save_comparison(
    histograms: Mapping[str, object],
    *,
    output: str | Path,
    x_label: str,
    y_label: str,
    ratios: Mapping[str, object] | None = None,
    ratio_label: str = "CUEP8M1 / CP5",
    x_range: tuple[float, float] | None = None,
    y_range: tuple[float, float] | None = None,
    ratio_range: tuple[float, float] | None = None,
    log_y: bool = False,
    annotations: Iterable[str] = (),
    save_png: bool = False,
    style: MplhepStyle = DEFAULT_MPLHEP_STYLE,
    cms_header_left: str = CMS_SIMULATION_HEADER,
    cms_header_right: str = PPB_MPLHEP_HEADER_RIGHT,
) -> None:
    """Save an 8x8-inch CMS-style overlay, optionally with a ratio panel."""

    if not histograms:
        raise ValueError("No histograms supplied")
    plt, hep = _mplhep()
    hep.style.use("CMS")
    if ratios:
        fig, (axis, ratio_axis) = plt.subplots(
            2, 1, figsize=(style.figure_size, style.figure_size),
            sharex=True, layout="constrained",
            gridspec_kw={"height_ratios": (3, 1)},
        )
    else:
        fig, axis = plt.subplots(
            figsize=(style.figure_size, style.figure_size), layout="constrained",
        )
        ratio_axis = None

    for index, (label, histogram) in enumerate(histograms.items()):
        values, edges, errors = _th1_arrays(histogram)
        display_label, color, marker, marker_face = _curve_style(label, index)
        hep.histplot(
            values, bins=edges, yerr=errors, xerr=True, histtype="errorbar",
            label=display_label, color=color, marker=marker,
            markerfacecolor=marker_face, markeredgecolor=color,
            markeredgewidth=1.5, markersize=style.marker_size,
            elinewidth=1.5, capsize=style.capsize,
            ax=axis,
        )
    axis.set_ylabel(y_label, fontsize=style.axis_title_size, labelpad=8)
    axis.tick_params(axis="both", labelsize=style.tick_label_size)
    axis.set_yscale("log" if log_y else "linear")
    if x_range is not None:
        axis.set_xlim(*x_range)
    else:
        reference = next(iter(histograms.values()))
        axis.set_xlim(*getattr(
            reference, "_root_display_x_range", _hist_axis_range(reference, "x"),
        ))
    if y_range is not None:
        axis.set_ylim(*y_range)
    else:
        reference = next(iter(histograms.values()))
        remembered = getattr(reference, "_root_display_y_range", None)
        if remembered is not None:
            axis.set_ylim(*remembered)
        else:
            stored_low = reference.GetMinimumStored()
            stored_high = reference.GetMaximumStored()
            automatic_low, automatic_high = axis.get_ylim()
            axis.set_ylim(
                automatic_low if stored_low == -1111.0 else stored_low,
                automatic_high if stored_high == -1111.0 else stored_high,
            )
    axis.legend(
        loc="upper right", fontsize=style.legend_text_size, markerscale=1.2,
    )
    axis.grid(True, alpha=0.25)
    _draw_cms_header(
        hep, axis, left=cms_header_left, right=cms_header_right, style=style,
    )
    annotation_lines = tuple(line for line in annotations if line)
    if annotation_lines:
        axis.text(
            0.04, style.annotation_y, "\n".join(annotation_lines),
            transform=axis.transAxes, va="top",
            fontsize=style.annotation_text_size,
        )

    if ratio_axis is None:
        axis.set_xlabel(x_label, fontsize=style.axis_title_size, labelpad=8)
    else:
        for index, (label, histogram) in enumerate(ratios.items()):
            values, edges, errors = _th1_arrays(histogram)
            display_label, color, marker, marker_face = _curve_style(label, index + 1)
            hep.histplot(
                values, bins=edges, yerr=errors, xerr=True, histtype="errorbar",
                label=display_label, color=color, marker=marker,
                markerfacecolor=marker_face, markeredgecolor=color,
                markeredgewidth=1.5, markersize=style.ratio_marker_size,
                elinewidth=1.5,
                capsize=style.capsize, ax=ratio_axis,
            )
        ratio_axis.axhline(1.0, color="black", linestyle="--", linewidth=1)
        ratio_axis.set_xlabel(
            x_label, fontsize=style.axis_title_size, labelpad=8,
        )
        ratio_axis.set_ylabel(
            ratio_label, fontsize=style.ratio_axis_title_size, labelpad=6,
        )
        ratio_axis.tick_params(
            axis="both", labelsize=style.ratio_tick_label_size,
        )
        if x_range is not None:
            ratio_axis.set_xlim(*x_range)
        if ratio_range is not None:
            ratio_axis.set_ylim(*ratio_range)
        else:
            ratio_reference = next(iter(ratios.values()))
            remembered_ratio = getattr(
                ratio_reference, "_root_display_y_range", None,
            )
            if remembered_ratio is not None:
                ratio_axis.set_ylim(*remembered_ratio)
        ratio_axis.grid(True, alpha=0.25)

    _save(fig, output, save_png=save_png, style=style)
    plt.close(fig)


def save_2d(
    histogram,
    *,
    output: str | Path,
    x_label: str,
    y_label: str,
    x_range: tuple[float, float] | None = None,
    y_range: tuple[float, float] | None = None,
    annotation: str = "",
    save_png: bool = False,
    style: MplhepStyle = DEFAULT_MPLHEP_STYLE,
    cms_header_left: str = CMS_SIMULATION_HEADER,
    cms_header_right: str = PPB_MPLHEP_HEADER_RIGHT,
) -> None:
    """Save one ROOT TH2 as an 8x8-inch CMS-style mplhep figure."""

    if not histogram.InheritsFrom("TH2"):
        raise TypeError(f"Expected TH2, got {histogram.ClassName()}")
    plt, hep = _mplhep()
    hep.style.use("CMS")
    fig, axis = plt.subplots(
        figsize=(style.figure_size, style.figure_size), layout="constrained",
    )
    values, x_edges, y_edges = _th2_arrays(histogram)
    artists = hep.hist2dplot(
        values, xbins=x_edges, ybins=y_edges, ax=axis, cbar=True,
        cbarsize="4%", cbarpad=0.08,
    )
    artists.cbar.ax.tick_params(labelsize=style.colorbar_tick_size)
    axis.set_xlabel(x_label, fontsize=style.axis_title_size, labelpad=8)
    axis.set_ylabel(y_label, fontsize=style.axis_title_size, labelpad=8)
    axis.tick_params(axis="both", labelsize=style.tick_label_size)
    if x_range is not None:
        axis.set_xlim(*x_range)
    else:
        axis.set_xlim(*_hist_axis_range(histogram, "x"))
    if y_range is not None:
        axis.set_ylim(*y_range)
    else:
        axis.set_ylim(*_hist_axis_range(histogram, "y"))
    _draw_cms_header(
        hep, axis, left=cms_header_left, right=cms_header_right, style=style,
    )
    if annotation:
        axis.text(
            0.04, style.annotation_y, annotation,
            transform=axis.transAxes, va="top",
            fontsize=style.annotation_text_size,
        )
    _save(fig, output, save_png=save_png, style=style)
    plt.close(fig)
