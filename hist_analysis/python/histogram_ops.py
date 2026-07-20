"""Normalization, common-bin alignment, ratios, and stable ROOT styles."""

from __future__ import annotations

from dataclasses import dataclass
from array import array
import math
from typing import Iterable, Mapping

import ROOT

from hist_analysis.python.root_style import set_1d_style


@dataclass(frozen=True)
class HistogramSpec:
    name: str
    title: str = ""


def normalize_labels(labels: Iterable[str]) -> list[str]:
    return [label.strip() for label in labels if label.strip()]


def normalize_histogram(histogram, mode: str = "none"):
    """Return a detached clone normalized according to ``mode``.

    Supported modes are ``none``, ``integral``, and ``bin_width``.  The latter
    divides by the integral and bin width, producing a unit-area density.
    """

    result = histogram.Clone(f"{histogram.GetName()}_{mode}")
    result.SetDirectory(0)
    if mode == "none":
        return result
    if mode not in {"integral", "bin_width"}:
        raise ValueError(f"Unsupported normalization mode: {mode!r}")

    integral = result.Integral(1, result.GetNbinsX())
    if integral <= 0:
        raise ValueError(f"Cannot normalize {histogram.GetName()!r}: integral is {integral}")
    result.Scale(1.0 / integral, "width" if mode == "bin_width" else "")
    return result


def _same_binning(left, right) -> bool:
    if left.GetNbinsX() != right.GetNbinsX():
        return False
    for bin_index in range(1, left.GetNbinsX() + 2):
        if abs(left.GetXaxis().GetBinLowEdge(bin_index) - right.GetXaxis().GetBinLowEdge(bin_index)) > 1e-9:
            return False
    return True


def align_common_binning(histograms: Mapping[str, object]) -> dict[str, object]:
    """Restrict 1D histograms to bin edges shared by every input.

    Input bins between consecutive shared edges are summed. This handles files
    with the same physical binning but different outer ranges, and fails rather
    than interpolating when no meaningful common bin intervals exist.
    """

    if not histograms:
        raise ValueError("No histograms supplied")

    edge_sets = []
    for histogram in histograms.values():
        axis = histogram.GetXaxis()
        edge_sets.append({round(axis.GetBinLowEdge(i), 9)
                          for i in range(1, histogram.GetNbinsX() + 2)})
    common_edges = sorted(set.intersection(*edge_sets))
    if len(common_edges) < 2:
        raise ValueError("Histograms have no common bin intervals")

    aligned = {}
    for index, (label, histogram) in enumerate(histograms.items()):
        result = ROOT.TH1D(
            f"{histogram.GetName()}_aligned_{index}", histogram.GetTitle(),
            len(common_edges) - 1, array("d", common_edges),
        )
        result.SetDirectory(0)
        result.Sumw2()
        axis = histogram.GetXaxis()
        for output_bin, (low, high) in enumerate(zip(common_edges, common_edges[1:]), 1):
            first = axis.FindFixBin(low)
            last = axis.FindFixBin(math.nextafter(high, -math.inf))
            content = sum(histogram.GetBinContent(i) for i in range(first, last + 1))
            error2 = sum(histogram.GetBinError(i) ** 2 for i in range(first, last + 1))
            result.SetBinContent(output_bin, content)
            result.SetBinError(output_bin, math.sqrt(error2))
        aligned[label] = result
    return aligned


def ratio_to_nominal(histogram, nominal, *, name: str | None = None):
    """Return ``histogram / nominal`` with standard independent errors."""

    if not _same_binning(histogram, nominal):
        raise ValueError(
            f"Incompatible binning: {histogram.GetName()!r} and {nominal.GetName()!r}"
        )
    ratio = histogram.Clone(name or f"{histogram.GetName()}_over_{nominal.GetName()}")
    ratio.SetDirectory(0)
    if histogram is nominal:
        for bin_index in range(0, ratio.GetNbinsX() + 2):
            ratio.SetBinContent(bin_index, 1.0 if nominal.GetBinContent(bin_index) != 0 else 0.0)
            ratio.SetBinError(bin_index, 0.0)
        return ratio
    ratio.Divide(nominal)
    return ratio


def style_histograms(histograms: Mapping[str, object]) -> None:
    """Apply the macro's stable 1D styles in mapping order."""

    for index, histogram in enumerate(histograms.values()):
        set_1d_style(histogram, index)
