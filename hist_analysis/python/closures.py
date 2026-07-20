"""Curve definitions and histogram preparation for closure comparisons."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from hist_analysis.python.histogram_io import load_histogram
from hist_analysis.python.histogram_ops import align_common_binning, normalize_histogram
from hist_analysis.python.projections import project_semantic_th2


@dataclass(frozen=True)
class ClosureCurve:
    """One curve in a closure plot."""

    label: str
    filename: Path
    histogram_names: tuple[str, ...]


def first_available_histogram(curve: ClosureCurve):
    """Load the first available key and return ``(histogram, key)``."""

    errors = []
    for name in curve.histogram_names:
        try:
            return load_histogram(str(curve.filename), name), name
        except KeyError as error:
            errors.append(str(error))
    raise KeyError(
        f"None of {curve.histogram_names} found in {curve.filename}: {'; '.join(errors)}"
    )


def build_closure_histograms(curves: Iterable[ClosureCurve], *, observable: str,
                             selection_range=None, normalization: str = "integral"):
    """Load, semantically project, align, and optionally normalize closure curves.

    The requested selection is interpreted as ``[low, high)`` by the projection
    helper. Returns the histogram mapping and a mapping recording the selected
    ROOT key for each curve.
    """

    projected_histograms = {}
    selected_keys = {}
    for index, curve in enumerate(curves):
        if curve.label in projected_histograms:
            raise ValueError(f"Duplicate closure label: {curve.label!r}")
        source, key = first_available_histogram(curve)
        projected = project_semantic_th2(
            source, observable, selection_range,
            name=f"closure_{index}_{observable}",
        )
        projected_histograms[curve.label] = projected
        selected_keys[curve.label] = key
    aligned = align_common_binning(projected_histograms)
    histograms = {
        label: normalize_histogram(histogram, normalization)
        for label, histogram in aligned.items()
    }
    return histograms, selected_keys
