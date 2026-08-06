"""ROOT projections and forward/backward folding for data checks."""

from __future__ import annotations

from array import array
from pathlib import Path

import ROOT
from IPython.display import display

from hist_analysis.python.histogram_io import load_histogram
from hist_analysis.python.histogram_ops import normalize_histogram, ratio_to_nominal
from hist_analysis.python.plotting import draw_closure, draw_overlay
from hist_analysis.python.projections import project_semantic_th2
from hist_analysis.python.root_style import (
    draw_text_block, save_canvas, set_2d_style, set_pad_style,
    set_palette_style,
)


def project_data_histogram(
    filename: str | Path,
    key: str,
    observable: str,
    selection_range: tuple[float, float] | None,
    *,
    name: str,
    rebin: int = 1,
    normalization: str = "none",
):
    """Project a stored (pT, eta) TH2 using the shared half-open bin convention."""

    if observable not in {"pt", "ptave", "eta"}:
        raise ValueError(f"Unsupported observable: {observable!r}")
    if isinstance(rebin, bool) or not isinstance(rebin, int) or rebin < 1:
        raise ValueError(f"rebin must be a positive integer, got {rebin!r}")
    source = load_histogram(str(filename), key)
    if not source.InheritsFrom("TH2"):
        raise TypeError(f"{key!r} in {filename} is not a TH2")
    source.Rebin2D(
        rebin if observable in {"pt", "ptave"} else 1,
        rebin if observable == "eta" else 1,
    )
    projection = project_semantic_th2(
        source, observable, selection_range, name=name,
    )
    return normalize_histogram(projection, normalization)


def forward_backward_from_full(full, *, name: str):
    """Fold a full CM eta histogram and form F/B with independent ROOT errors."""

    axis = full.GetXaxis()
    middle_bin = axis.FindBin(0.0 + 0.001)
    positive_bins = full.GetNbinsX() - middle_bin + 1
    negative_bins = middle_bin - 1
    output_bins = min(positive_bins, negative_bins)
    if output_bins < 1:
        raise ValueError(f"{full.GetName()} has no bins on both sides of eta=0")

    edges = [axis.GetBinLowEdge(middle_bin)]
    edges.extend(
        axis.GetBinUpEdge(middle_bin + offset)
        for offset in range(output_bins)
    )
    forward = ROOT.TH1D(f"{name}_forward", "", output_bins, array("d", edges))
    backward = ROOT.TH1D(f"{name}_backward", "", output_bins, array("d", edges))
    forward.Sumw2()
    backward.Sumw2()
    forward.SetDirectory(0)
    backward.SetDirectory(0)

    for output_bin in range(1, output_bins + 1):
        forward_bin = middle_bin + output_bin - 1
        backward_bin = middle_bin - output_bin
        forward.SetBinContent(output_bin, full.GetBinContent(forward_bin))
        forward.SetBinError(output_bin, full.GetBinError(forward_bin))
        backward.SetBinContent(output_bin, full.GetBinContent(backward_bin))
        backward.SetBinError(output_bin, full.GetBinError(backward_bin))

    ratio = forward.Clone(name)
    ratio.SetDirectory(0)
    ratio.Divide(forward, backward, 1.0, 1.0, "")
    ratio._forward_backward_inputs = (forward, backward)
    return ratio


def _tag(text: str) -> str:
    return "".join(
        character.lower() if character.isalnum() else "_" for character in text
    ).strip("_")


def draw_data_distributions(
    trigger: str,
    filename: str | Path,
    frames,
    pt_bins,
    *,
    jet_kind: str,
    output_dir: str | Path,
    rebin_pt: int,
    rebin_eta: int,
    pt_display_range,
    pt_eta_range,
    fb_y_range,
    save_png: bool,
    grid: bool,
    style,
):
    """Draw independent pT/eta frame checks and CM F/B for one trigger."""

    if jet_kind not in {"jet", "dijet"}:
        raise ValueError("jet_kind must be 'jet' or 'dijet'")
    path = Path(filename)
    output = Path(output_dir)
    trigger_tag = _tag(trigger)
    momentum = "p_{T}^{jet}" if jet_kind == "jet" else "p_{T}^{ave}"
    object_label = "Inclusive jets" if jet_kind == "jet" else "Dijets"
    pt_observable = "pt" if jet_kind == "jet" else "ptave"
    results = {}
    for frame, (frame_label, key, eta_title) in frames.items():
        frame_tag = _tag(frame)
        pt_selection = pt_eta_range if jet_kind == "jet" else None
        pt = project_data_histogram(
            path, key, pt_observable, pt_selection,
            name=f"h_{trigger_tag}_{frame_tag}_{jet_kind}_pt",
            rebin=rebin_pt,
        )
        pt_canvas = draw_overlay(
            {trigger: pt}, title="", x_title=f"{momentum} (GeV)",
            y_title=(
                "dN/dp_{T}^{jet}" if jet_kind == "jet"
                else "dN/dp_{T}^{ave}"
            ), x_range=pt_display_range,
            log_y=True, annotations=(
                "pPb 8.16 TeV data", trigger, frame_label, object_label,
            ), style_indices={trigger: 2}, grid=grid, style=style,
            output=output / f"{trigger_tag}_{frame_tag}_{jet_kind}_pt.pdf",
            save_png=save_png, headroom=2.2,
            canvas_name=f"{trigger_tag}_{frame_tag}_{jet_kind}_pt",
        )
        display(pt_canvas)
        eta_results = {}
        for pt_range in pt_bins:
            low, high = pt_range
            pt_tag = f"{low:g}_{high:g}".replace(".", "p")
            raw_eta = project_data_histogram(
                path, key, "eta", pt_range,
                name=f"h_{trigger_tag}_{frame_tag}_{jet_kind}_eta_{pt_tag}",
                rebin=rebin_eta,
            )
            eta = normalize_histogram(raw_eta, "integral")
            annotations = (
                "pPb 8.16 TeV data", trigger, frame_label, object_label,
                f"{low:g} #leq {momentum} < {high:g} GeV",
            )
            eta_canvas = draw_overlay(
                {trigger: eta}, title="", x_title=eta_title,
                y_title=f"1/N dN/d{eta_title}", annotations=annotations,
                style_indices={trigger: 2}, grid=grid, style=style,
                headroom=2.2,
                output=output / (
                    f"{trigger_tag}_{frame_tag}_{jet_kind}_eta_pt_{pt_tag}.pdf"
                ), save_png=save_png,
                canvas_name=(
                    f"{trigger_tag}_{frame_tag}_{jet_kind}_eta_{pt_tag}"
                ),
            )
            display(eta_canvas)
            entry = {"raw_eta": raw_eta, "eta": eta, "eta_canvas": eta_canvas}
            if frame == "cm":
                fb = forward_backward_from_full(
                    raw_eta,
                    name=f"h_{trigger_tag}_{jet_kind}_fb_{pt_tag}",
                )
                fb_canvas = draw_overlay(
                    {trigger: fb}, title="",
                    x_title=f"|{eta_title}|", y_title="Forward / Backward",
                    y_range=fb_y_range, annotations=annotations,
                    style_indices={trigger: 2}, grid=grid, style=style,
                    reference_y=1.0,
                    legend_bounds=(0.22, 0.20, 0.48, 0.29),
                    output=output / (
                        f"{trigger_tag}_cm_{jet_kind}_fb_pt_{pt_tag}.pdf"
                    ), save_png=save_png,
                    canvas_name=f"{trigger_tag}_cm_{jet_kind}_fb_{pt_tag}",
                )
                display(fb_canvas)
                entry.update({"forward_backward": fb, "fb_canvas": fb_canvas})
            eta_results[pt_range] = entry
        results[frame] = {
            "pt": pt, "pt_canvas": pt_canvas, "eta": eta_results,
        }
    return results


def draw_orientation_comparisons(
    trigger: str,
    direction_files,
    frames,
    pt_bins,
    *,
    jet_kind: str,
    output_dir: str | Path,
    rebin_eta: int,
    rebin_pt: int,
    pt_display_range,
    pt_eta_range,
    pt_normalization_range,
    ratio_range,
    selection: str,
    save_png: bool,
    grid: bool,
    style,
    eta_display_range: tuple[float, float] | None = None,
    eta_y_range: tuple[float, float] | None = None,
    fb_y_range: tuple[float, float] | None = None,
    fb_ratio_range: tuple[float, float] | None = None,
):
    """Overlay orientation shapes and CM F/B ratios relative to combined."""

    if jet_kind not in {"jet", "dijet"}:
        raise ValueError("jet_kind must be 'jet' or 'dijet'")
    output = Path(output_dir)
    trigger_tag = _tag(trigger)
    default_momentum = "p_{T}^{jet}" if jet_kind == "jet" else "p_{T}^{ave}"
    object_label = "Inclusive jets" if jet_kind == "jet" else "Dijets"
    styles = {"Pb-going": 1, "p-going": 0, "combined": 2}
    results = {}
    for frame, frame_config in frames.items():
        if len(frame_config) == 3:
            frame_label, key, eta_title = frame_config
            momentum = default_momentum
        elif len(frame_config) == 4:
            frame_label, key, eta_title, momentum = frame_config
        else:
            raise ValueError(
                f"Frame {frame!r} must contain 3 or 4 configuration values"
            )
        maps = {
            direction: load_histogram(str(filename), key)
            for direction, filename in direction_files.items()
        }
        for histogram in maps.values():
            if not histogram.InheritsFrom("TH2"):
                raise TypeError(f"{histogram.GetName()!r} is not a TH2")
        positive_values = [
            histogram.GetBinContent(x_bin, y_bin)
            for histogram in maps.values()
            for x_bin in range(1, histogram.GetNbinsX() + 1)
            for y_bin in range(1, histogram.GetNbinsY() + 1)
            if histogram.GetBinContent(x_bin, y_bin) > 0.0
        ]
        if not positive_values:
            raise ValueError(f"All {frame_label} orientation maps are empty")
        common_z_range = (min(positive_values), max(positive_values))
        map_canvases = {}
        for direction, histogram in maps.items():
            direction_tag = _tag(direction)
            canvas = ROOT.TCanvas(
                f"c_{trigger_tag}_{frame}_{jet_kind}_{direction_tag}_2d", "",
                style.canvas_width, style.canvas_height,
            )
            set_pad_style(canvas, grid_x=False, grid_y=False, style=style)
            canvas.SetRightMargin(style.palette_right_margin)
            canvas.SetLogz(True)
            histogram.SetTitle("")
            histogram.GetXaxis().SetTitle(f"{momentum} (GeV)")
            histogram.GetYaxis().SetTitle(eta_title)
            histogram.GetZaxis().SetTitle(f"{object_label} / bin")
            histogram.SetMinimum(common_z_range[0])
            histogram.SetMaximum(common_z_range[1])
            set_2d_style(histogram, style=style)
            histogram.Draw("COLZ")
            canvas.Update()
            palette = set_palette_style(histogram, style=style)
            annotations_2d = draw_text_block(
                canvas, (
                    "pPb 8.16 TeV data", trigger, direction, frame_label,
                    object_label, selection,
                ), style=style,
            )
            canvas.Modified()
            canvas.Update()
            map_output = output / (
                f"{trigger_tag}_{frame}_{jet_kind}_{direction_tag}_2d.pdf"
            )
            save_canvas(canvas, map_output, save_png=save_png)
            canvas._data_map_objects = [
                histogram, palette, *annotations_2d,
            ]
            display(canvas)
            map_canvases[direction] = canvas
        results[(frame, "2d")] = {
            "histograms": maps, "canvases": map_canvases,
            "common_z_range": common_z_range,
        }
        pt_selection = pt_eta_range if jet_kind == "jet" else None
        pt_histograms = {}
        for direction, filename in direction_files.items():
            raw_pt = project_data_histogram(
                filename, key, "pt" if jet_kind == "jet" else "ptave",
                pt_selection,
                name=(
                    f"h_{trigger_tag}_{frame}_{_tag(direction)}_"
                    f"{jet_kind}_pt"
                ), rebin=rebin_pt,
            )
            pt_histograms[direction] = normalize_histogram(
                raw_pt, "integral", integral_range=pt_normalization_range,
            )
        pt_canvas, pt_ratios = draw_closure(
            pt_histograms, "combined", title="",
            x_title=f"{momentum} (GeV)",
            y_title=(
                "1/N dN/dp_{T}^{jet}" if jet_kind == "jet"
                else "1/N dN/dp_{T}^{ave}"
            ), ratio_range=ratio_range, ratio_option="",
            draw_nominal_ratio=False, style_indices=styles,
            annotations=(
                "pPb 8.16 TeV data", trigger, frame_label, object_label,
                selection,
                (
                    f"Normalized in {pt_normalization_range[0]:g} #leq "
                    f"{momentum} < {pt_normalization_range[1]:g} GeV"
                ),
            ), x_range=pt_display_range, log_y=True, headroom=2.2,
            grid=grid, style=style,
            output=output / (
                f"{trigger_tag}_{frame}_{jet_kind}_orientation_pt.pdf"
            ), save_png=save_png,
            canvas_name=f"{trigger_tag}_{frame}_{jet_kind}_orientation_pt",
        )
        display(pt_canvas)
        results[(frame, "pt")] = {
            "histograms": pt_histograms,
            "ratios_to_combined": pt_ratios,
            "canvas": pt_canvas,
        }
        for pt_range in pt_bins:
            low, high = pt_range
            pt_tag = f"{low:g}_{high:g}".replace(".", "p")
            raw_histograms = {
                direction: project_data_histogram(
                    filename, key, "eta", pt_range,
                    name=(
                        f"h_{trigger_tag}_{frame}_{_tag(direction)}_"
                        f"{jet_kind}_{pt_tag}"
                    ), rebin=rebin_eta,
                )
                for direction, filename in direction_files.items()
            }
            histograms = {
                direction: normalize_histogram(histogram, "integral")
                for direction, histogram in raw_histograms.items()
            }
            canvas, ratios = draw_closure(
                histograms, "combined", title="", x_title=eta_title,
                y_title=f"1/N dN/d{eta_title}", ratio_range=ratio_range,
                ratio_option="", draw_nominal_ratio=False,
                style_indices=styles, annotations=(
                    "pPb 8.16 TeV data", trigger, frame_label, object_label,
                    selection,
                    f"{low:g} #leq {momentum} < {high:g} GeV",
                ), x_range=eta_display_range, y_range=eta_y_range,
                grid=grid, style=style, headroom=2.2,
                output=output / (
                    f"{trigger_tag}_{frame}_{jet_kind}_orientation_pt_{pt_tag}.pdf"
                ), save_png=save_png,
                canvas_name=(
                    f"{trigger_tag}_{frame}_{jet_kind}_orientation_{pt_tag}"
                ),
            )
            display(canvas)
            results[(frame, pt_range)] = {
                "raw_histograms": raw_histograms,
                "histograms": histograms,
                "ratios_to_combined": ratios,
                "canvas": canvas,
            }
            if frame == "cm" and jet_kind == "dijet":
                forward_backward = {
                    direction: forward_backward_from_full(
                        histogram,
                        name=(
                            f"h_{trigger_tag}_{_tag(direction)}_"
                            f"{jet_kind}_fb_{pt_tag}"
                        ),
                    )
                    for direction, histogram in raw_histograms.items()
                }
                fb_canvas, fb_ratios = draw_closure(
                    forward_backward, "combined", title="",
                    x_title=f"|{eta_title}|",
                    y_title="Forward / Backward",
                    y_range=fb_y_range,
                    ratio_range=fb_ratio_range or ratio_range,
                    ratio_option="", draw_nominal_ratio=False,
                    style_indices=styles, annotations=(
                        "pPb 8.16 TeV data", trigger, frame_label,
                        object_label, selection,
                        f"{low:g} #leq {momentum} < {high:g} GeV",
                    ), grid=grid, style=style, headroom=2.2,
                    output=output / (
                        f"{trigger_tag}_{frame}_{jet_kind}_"
                        f"orientation_fb_pt_{pt_tag}.pdf"
                    ), save_png=save_png,
                    canvas_name=(
                        f"{trigger_tag}_{frame}_{jet_kind}_"
                        f"orientation_fb_{pt_tag}"
                    ),
                )
                display(fb_canvas)
                results[(frame, pt_range)].update({
                    "forward_backward": forward_backward,
                    "fb_ratios_to_combined": fb_ratios,
                    "fb_canvas": fb_canvas,
                })
    return results


def draw_corrected_raw_eta_comparisons(
    trigger: str,
    direction_files,
    pt_bins,
    *,
    corrected_key: str,
    raw_key: str,
    output_dir: str | Path,
    rebin_eta: int,
    eta_display_range: tuple[float, float] | None,
    ratio_range: tuple[float, float] | None,
    direction_ratio_range: tuple[float, float],
    selection: str,
    save_png: bool,
    grid: bool,
    style,
):
    """Compare unflipped-Lab corrected/raw eta ratios by beam direction."""

    directions = ("Pb-going", "p-going")
    missing = [
        direction for direction in directions
        if direction not in direction_files
    ]
    if missing:
        raise KeyError(f"Missing direction files for: {', '.join(missing)}")

    output = Path(output_dir)
    trigger_tag = _tag(trigger)
    styles = {"Pb-going": 1, "p-going": 0}
    results = {}
    for low, high in pt_bins:
        pt_range = (low, high)
        pt_tag = f"{low:g}_{high:g}".replace(".", "p")
        ratios = {}
        inputs = {}
        for direction in directions:
            filename = direction_files[direction]
            corrected = project_data_histogram(
                filename, corrected_key, "eta", pt_range,
                name=(
                    f"h_{trigger_tag}_{_tag(direction)}_"
                    f"corrected_eta_{pt_tag}"
                ),
                rebin=rebin_eta,
            )
            raw = project_data_histogram(
                filename, raw_key, "eta", pt_range,
                name=f"h_{trigger_tag}_{_tag(direction)}_raw_eta_{pt_tag}",
                rebin=rebin_eta,
            )
            ratios[direction] = ratio_to_nominal(
                corrected, raw,
                name=(
                    f"h_{trigger_tag}_{_tag(direction)}_"
                    f"corrected_over_raw_eta_{pt_tag}"
                ),
                option="",
            )
            inputs[direction] = {"corrected": corrected, "raw": raw}

        canvas, ratios_to_pgoing = draw_closure(
            ratios, "p-going", title="",
            x_title="#eta_{Lab,unflipped}^{jet}",
            y_title="Corrected / raw", x_range=eta_display_range,
            y_range=ratio_range, ratio_range=direction_ratio_range,
            ratio_option="", draw_nominal_ratio=False,
            annotations=(
                "pPb 8.16 TeV data", trigger, "Lab unflipped",
                "Inclusive jets", selection,
                f"{low:g} #leq p_{{T}}^{{jet}} < {high:g} GeV",
            ),
            style_indices=styles, grid=grid, style=style,
            output=output / (
                f"{trigger_tag}_lab_unflipped_jet_"
                f"corrected_over_raw_eta_pt_{pt_tag}.pdf"
            ),
            save_png=save_png,
            canvas_name=(
                f"{trigger_tag}_lab_unflipped_jet_"
                f"corrected_over_raw_eta_{pt_tag}"
            ),
        )
        display(canvas)
        results[pt_range] = {
            "inputs": inputs,
            "ratios": ratios,
            "corrected_over_raw": ratios,
            "ratios_to_pgoing": ratios_to_pgoing,
            "canvas": canvas,
        }
    return results
