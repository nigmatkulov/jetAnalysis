"""Data-trigger pT projections, scale extraction, and stitching plots."""

from __future__ import annotations

import ctypes
from dataclasses import dataclass
import json
from math import hypot
from pathlib import Path

import ROOT

from .histogram_io import load_histogram
from .projections import project_semantic_th2
from .root_style import (
    DEFAULT_PLOT_STYLE,
    draw_text_block,
    save_canvas,
    set_1d_style,
    set_legend_style,
)


@dataclass(frozen=True)
class StitchingScale:
    """Multiplicative candidate-to-reference trigger normalization."""

    reference: str
    candidate: str
    pt_range: tuple[float, float]
    local_scale: float
    local_uncertainty: float
    cumulative_scale: float
    cumulative_uncertainty: float
    reference_integral: float
    candidate_integral: float


@dataclass(frozen=True)
class TriggerStitchingConfig:
    """Validated cumulative trigger weights and their provenance."""

    source: Path
    tag: str
    weights: dict[str, float]
    uncertainties: dict[str, float]
    histogram_keys: dict[str, str]
    metadata: dict


def load_trigger_stitching_config(
        filename, *, direction: str, selection: str,
        eta_range: tuple[float, float], eta_cut_index: int,
        observable: str = "dijet", frame: str = "cm",
        tolerance: float = 1.0e-9) -> TriggerStitchingConfig:
    """Load and validate a trigger-stitching JSON artifact.

    The returned weights multiply raw trigger histograms onto the MinimumBias
    scale. A provenance mismatch fails explicitly rather than applying weights
    extracted with a different observable or acceptance.
    """

    path = Path(filename).expanduser().resolve()
    if path.suffix.lower() != ".json":
        raise ValueError(f"Stitching weights must be a JSON file: {path}")
    if not path.is_file():
        raise FileNotFoundError(path)
    metadata = json.loads(path.read_text())
    expected_scalars = {
        "direction": direction,
        "selection": selection,
        "observable": observable,
        "frame": frame,
        "dijet_eta_cut_index": eta_cut_index,
    }
    mismatches = [
        f"{key}={metadata.get(key)!r}, expected {expected!r}"
        for key, expected in expected_scalars.items()
        if metadata.get(key) != expected
    ]
    stored_eta_range = metadata.get("projection_eta_range")
    if (not isinstance(stored_eta_range, list)
            or len(stored_eta_range) != 2
            or any(abs(float(actual) - expected) > tolerance
                   for actual, expected in zip(stored_eta_range, eta_range))):
        mismatches.append(
            f"projection_eta_range={stored_eta_range!r}, expected {eta_range!r}"
        )
    if mismatches:
        raise ValueError(
            f"Incompatible trigger-stitching artifact {path}: "
            + "; ".join(mismatches)
        )

    weights = {"MinimumBias": 1.0}
    uncertainties = {"MinimumBias": 0.0}
    for record in metadata.get("scales", ()):
        trigger = record.get("candidate")
        if trigger in weights:
            raise ValueError(f"Duplicate stitching weight for {trigger!r} in {path}")
        scale = float(record["cumulative_scale"])
        uncertainty = float(record["cumulative_uncertainty"])
        if not (scale > 0.0 and uncertainty >= 0.0):
            raise ValueError(
                f"Invalid cumulative weight for {trigger!r}: "
                f"{scale} +/- {uncertainty}"
            )
        weights[trigger] = scale
        uncertainties[trigger] = uncertainty
    required = {"MinimumBias", "Jet60", "Jet80", "Jet100"}
    missing = sorted(required - set(weights))
    if missing:
        raise ValueError(f"Missing stitching weights in {path}: {missing}")

    histogram_keys = metadata.get("histogram_keys", {})
    missing_keys = sorted(required - set(histogram_keys))
    if missing_keys:
        raise ValueError(f"Missing histogram keys in {path}: {missing_keys}")
    return TriggerStitchingConfig(
        source=path,
        tag=str(metadata.get("stitching_tag", path.stem)),
        weights=weights,
        uncertainties=uncertainties,
        histogram_keys={key: str(value) for key, value in histogram_keys.items()},
        metadata=metadata,
    )


def dijet_histogram_key(filename, frame: str, eta_cut_index: int) -> str:
    """Resolve an indexed dijet key, with the merged 1.9 fallback."""

    if frame not in {"lab", "cm"}:
        raise ValueError("frame must be 'lab' or 'cm'")
    base = f"hRecoDijetPtEta{'Lab' if frame == 'lab' else 'CM'}"
    indexed = f"{base}_{eta_cut_index}"
    try:
        load_histogram(str(filename), indexed)
    except KeyError:
        pass
    else:
        return indexed
    if eta_cut_index == 5:
        try:
            load_histogram(str(filename), base)
        except KeyError:
            pass
        else:
            return base
    raise KeyError(
        f"Neither {indexed!r} nor an allowed standard-cut fallback exists in "
        f"{filename}"
    )


def histogram_key(filename, observable: str, frame: str,
                  eta_cut_index: int = 5) -> str:
    """Return the configured inclusive-jet or dijet TH2 key."""

    if observable == "inclusive":
        if frame not in {"lab", "cm"}:
            raise ValueError("frame must be 'lab' or 'cm'")
        return f"hRecoInclusiveJetPtEta{'Lab' if frame == 'lab' else 'CM'}"
    if observable == "dijet":
        return dijet_histogram_key(filename, frame, eta_cut_index)
    raise ValueError("observable must be 'inclusive' or 'dijet'")


def load_pt_projection(filename, *, observable: str, frame: str,
                       eta_range: tuple[float, float], eta_cut_index: int = 5,
                       name: str):
    """Load a TH2 and project its pT axis in a half-open eta interval."""

    key = histogram_key(filename, observable, frame, eta_cut_index)
    source = load_histogram(str(filename), key)
    projection = project_semantic_th2(
        source, "pt" if observable == "inclusive" else "ptave",
        eta_range, name=name,
    )
    if projection.GetSumw2N() == 0:
        projection.Sumw2()
    return projection, key


def _integral_and_error(histogram, pt_range: tuple[float, float]):
    axis = histogram.GetXaxis()
    low, high = pt_range
    first = max(1, axis.FindBin(low + 0.001))
    last = min(axis.GetNbins(), axis.FindBin(high - 0.001))
    error = ctypes.c_double(0.0)
    value = histogram.IntegralAndError(first, last, error)
    return float(value), error.value


def derive_scale(reference, candidate, *, reference_label: str,
                 candidate_label: str, pt_range: tuple[float, float],
                 reference_cumulative_scale: float = 1.0,
                 reference_cumulative_uncertainty: float = 0.0) -> StitchingScale:
    """Derive an integral-ratio scale with independent statistical errors.

    The local factor is integral(reference)/integral(candidate). The cumulative
    candidate-to-MB factor includes the already-applied reference scale.
    """

    reference_integral, reference_error = _integral_and_error(reference, pt_range)
    candidate_integral, candidate_error = _integral_and_error(candidate, pt_range)
    if reference_integral <= 0.0 or candidate_integral <= 0.0:
        raise ValueError(
            f"Non-positive matching integral in {pt_range}: "
            f"{reference_integral=}, {candidate_integral=}"
        )
    local = reference_integral / candidate_integral
    local_error = local * hypot(
        reference_error / reference_integral,
        candidate_error / candidate_integral,
    )
    cumulative = reference_cumulative_scale * local
    cumulative_error = cumulative * hypot(
        reference_cumulative_uncertainty / reference_cumulative_scale
        if reference_cumulative_scale else 0.0,
        local_error / local,
    )
    return StitchingScale(
        reference_label, candidate_label, pt_range, local, local_error,
        cumulative, cumulative_error, reference_integral, candidate_integral,
    )


def scaled_clone(histogram, factor: float, name: str):
    clone = histogram.Clone(name)
    clone.SetDirectory(0)
    clone.Scale(factor)
    return clone


def make_ratio(numerator, denominator, name: str):
    """Form a trigger ratio using independent-error propagation."""

    ratio = numerator.Clone(name)
    ratio.SetDirectory(0)
    ratio.Divide(numerator, denominator, 1.0, 1.0, "")
    return ratio


def _configure_axes(histogram, *, x_title: str, y_title: str,
                    display_range: tuple[float, float], y_range=None):
    histogram.SetTitle("")
    histogram.GetXaxis().SetTitle(x_title)
    histogram.GetYaxis().SetTitle(y_title)
    histogram.GetXaxis().SetRangeUser(*display_range)
    if y_range is not None:
        histogram.SetMinimum(y_range[0])
        histogram.SetMaximum(y_range[1])


def _short_trigger_label(label: str) -> str:
    """Return a compact trigger name suitable for a ratio-axis title."""

    trigger = label.split()[0]
    if trigger in {"MinimumBias", "MB"}:
        return "MB"
    if trigger.startswith("Jet"):
        return trigger.replace("Jet", "J", 1)
    return trigger


def draw_comparison(reference, candidate, *, reference_label: str,
                    candidate_label: str, annotation_lines, output_path,
                    display_range=(40.0, 500.0), ratio_range=(0.0, 2.0),
                    log_y=True, save_png=False, grid=True,
                    style=DEFAULT_PLOT_STYLE):
    """Draw an unnormalized overlay and candidate/reference ratio."""

    canvas = ROOT.TCanvas(
        f"c_{candidate.GetName()}_over_{reference.GetName()}", "",
        style.canvas_width, style.canvas_height,
    )
    top = ROOT.TPad("top", "", 0.0, style.ratio_fraction, 1.0, 1.0)
    bottom = ROOT.TPad("bottom", "", 0.0, 0.0, 1.0, style.ratio_fraction)
    top.SetLeftMargin(style.ratio_left_margin)
    top.SetBottomMargin(style.ratio_top_bottom_margin)
    top.SetTopMargin(style.top_margin)
    top.SetRightMargin(style.right_margin)
    bottom.SetLeftMargin(style.ratio_left_margin)
    bottom.SetBottomMargin(style.ratio_bottom_margin)
    bottom.SetTopMargin(style.ratio_bottom_top_margin)
    bottom.SetRightMargin(style.right_margin)
    top.SetGridx(grid); top.SetGridy(grid)
    bottom.SetGridx(grid); bottom.SetGridy(grid)
    top.SetLogy(log_y)
    top.Draw(); bottom.Draw()

    set_1d_style(reference, 2, style=style)
    set_1d_style(candidate, 0, style=style)
    _configure_axes(reference, x_title="", y_title="Stored yield",
                    display_range=display_range)
    reference.GetXaxis().SetLabelSize(0.0)
    reference.GetXaxis().SetTickLength(0.0)
    maximum = max(reference.GetMaximum(), candidate.GetMaximum())
    positive = [h.GetBinContent(i) for h in (reference, candidate)
                for i in range(1, h.GetNbinsX() + 1) if h.GetBinContent(i) > 0]
    reference.SetMaximum((5.0 if log_y else 1.25) * maximum)
    reference.SetMinimum(0.5 * min(positive) if log_y and positive else 0.0)
    top.cd(); reference.Draw("E1"); candidate.Draw("E1 SAME")
    legend = ROOT.TLegend(0.61, 0.71, 0.92, 0.88)
    set_legend_style(legend, style=style)
    legend.AddEntry(reference, reference_label, "p")
    legend.AddEntry(candidate, candidate_label, "p")
    legend.Draw()
    labels = draw_text_block(top, annotation_lines, x=0.23, y=0.38,
                             text_size=0.029, style=style)

    ratio = make_ratio(candidate, reference,
                       f"hRatio_{candidate.GetName()}_over_{reference.GetName()}")
    set_1d_style(ratio, 0, style=style)
    _configure_axes(ratio, x_title="p_{T} (GeV)",
                    y_title=(f"{_short_trigger_label(candidate_label)} / "
                             f"{_short_trigger_label(reference_label)}"),
                    display_range=display_range, y_range=ratio_range)
    ratio.GetXaxis().SetTitleSize(0.10)
    ratio.GetXaxis().SetLabelSize(0.085)
    ratio.GetYaxis().SetTitleSize(0.095)
    ratio.GetYaxis().SetLabelSize(0.080)
    ratio.GetYaxis().SetTitleOffset(0.70)
    ratio.GetYaxis().SetNdivisions(style.ratio_axis_divisions)
    bottom.cd(); ratio.Draw("E1")
    unity = ROOT.TLine(display_range[0], 1.0, display_range[1], 1.0)
    unity.SetLineStyle(2); unity.Draw("SAME")
    canvas.Update()
    save_canvas(canvas, output_path, save_png=save_png)
    canvas._stitching_objects = (top, bottom, legend, *labels, ratio, unity,
                                 reference, candidate)
    return {"canvas": canvas, "ratio": ratio}


def draw_final_overlay(histograms: dict[str, object], *, annotation_lines,
                       output_path, display_range=(40.0, 500.0),
                       ratio_range=(0.0, 2.0), log_y=True, save_png=False,
                       grid=True, style=DEFAULT_PLOT_STYLE):
    """Draw all MB-normalized trigger spectra and ratios to MB."""

    reference = histograms["MinimumBias"]
    canvas = ROOT.TCanvas("c_data_trigger_stitching_final", "", 800, 900)
    top = ROOT.TPad("final_top", "", 0.0, 0.36, 1.0, 1.0)
    bottom = ROOT.TPad("final_bottom", "", 0.0, 0.0, 1.0, 0.36)
    for pad in (top, bottom):
        pad.SetLeftMargin(0.18); pad.SetRightMargin(0.05)
        pad.SetGridx(grid); pad.SetGridy(grid)
    top.SetBottomMargin(0.02); top.SetTopMargin(0.08); top.SetLogy(log_y)
    bottom.SetTopMargin(0.03); bottom.SetBottomMargin(0.30)
    top.Draw(); bottom.Draw()
    retained = [top, bottom]
    legend = ROOT.TLegend(0.62, 0.66, 0.93, 0.90)
    set_legend_style(legend, style=style)
    maximum_label, maximum_histogram = max(
        histograms.items(), key=lambda item: item[1].GetMaximum(),
    )
    maximum_bin = maximum_histogram.GetMaximumBin()
    maximum = maximum_histogram.GetBinContent(maximum_bin)
    maximum_axis = maximum_histogram.GetXaxis()
    maximum_bin_center = maximum_axis.GetBinCenter(maximum_bin)
    maximum_bin_low = maximum_axis.GetBinLowEdge(maximum_bin)
    maximum_bin_high = maximum_axis.GetBinUpEdge(maximum_bin)
    print(
        "Final-overlay global maximum: "
        f"trigger={maximum_label}, value={maximum:.8g}, "
        f"x_bin={maximum_bin}, x_center={maximum_bin_center:.8g} GeV, "
        f"x_edges=[{maximum_bin_low:.8g}, {maximum_bin_high:.8g}) GeV"
    )
    positive = [hist.GetBinContent(i) for hist in histograms.values()
                for i in range(1, hist.GetNbinsX() + 1)
                if hist.GetBinContent(i) > 0]
    top.cd()
    for index, (label, histogram) in enumerate(histograms.items()):
        set_1d_style(histogram, (2, 0, 1, 3)[index], style=style)
        _configure_axes(histogram, x_title="", y_title="Stored yield (MB scale)",
                        display_range=display_range)
        histogram.GetXaxis().SetLabelSize(0.0)
        histogram.GetXaxis().SetTickLength(0.0)
        histogram.SetMaximum(5.0 * maximum if log_y else 1.25 * maximum)
        histogram.SetMinimum(0.5 * min(positive) if log_y and positive else 0.0)
        histogram.Draw("E1" if index == 0 else "E1 SAME")
        legend.AddEntry(histogram, label, "p")
        retained.append(histogram)
    legend.Draw()
    retained.extend((legend, *draw_text_block(
        top, annotation_lines, x=0.22, y=0.42, text_size=0.028, style=style,
    )))

    ratios = {}
    bottom.cd()
    for index, (label, histogram) in enumerate(tuple(histograms.items())[1:], 1):
        ratio = make_ratio(histogram, reference,
                           f"hRatio_{label}_over_MinimumBias_final")
        set_1d_style(ratio, (0, 1, 3)[index - 1], style=style)
        _configure_axes(ratio, x_title="p_{T} (GeV)",
                        y_title="Trigger / MB", display_range=display_range,
                        y_range=ratio_range)
        ratio.GetXaxis().SetTitleSize(0.09); ratio.GetXaxis().SetLabelSize(0.08)
        ratio.GetYaxis().SetTitleSize(0.085); ratio.GetYaxis().SetLabelSize(0.075)
        ratio.GetYaxis().SetTitleOffset(0.80); ratio.GetYaxis().SetNdivisions(505)
        ratio.Draw("E1" if index == 1 else "E1 SAME")
        ratios[label] = ratio; retained.append(ratio)
    unity = ROOT.TLine(display_range[0], 1.0, display_range[1], 1.0)
    unity.SetLineStyle(2); unity.Draw("SAME"); retained.append(unity)
    canvas.Update(); save_canvas(canvas, output_path, save_png=save_png)
    canvas._stitching_objects = tuple(retained)
    return {"canvas": canvas, "ratios": ratios}
