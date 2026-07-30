"""Range-aware ROOT projections, including semantic jet-axis detection."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import ROOT


def _semantic_axis_index(obj, semantic: str) -> int:
    """Locate pT or eta using ranges, then validate against axis titles.

    Some existing reference-jet files have swapped/misleading axis titles, so
    title-only identification is unsafe.  Jet pT axes extend well above 20 GeV,
    whereas eta axes in these outputs are bounded within |eta| < 10.
    """

    if not obj.InheritsFrom("TH2"):
        raise TypeError("Semantic projection currently requires a TH2 histogram")
    axes = (obj.GetXaxis(), obj.GetYaxis())
    if semantic in {"pt", "ptave"}:
        candidates = [i for i, axis in enumerate(axes) if axis.GetXmax() > 20]
    elif semantic == "eta":
        candidates = [i for i, axis in enumerate(axes)
                      if axis.GetXmin() >= -10 and axis.GetXmax() <= 10]
    else:
        raise ValueError(f"Unsupported semantic axis: {semantic!r}")
    if len(candidates) != 1:
        ranges = [(axis.GetXmin(), axis.GetXmax(), axis.GetTitle()) for axis in axes]
        raise ValueError(f"Cannot identify {semantic!r} axis from {ranges}")
    return candidates[0]


def project_semantic_th2(obj, observable: str, selection_range=None, *, name: str | None = None):
    """Project a jet TH2 onto eta, pT, or pTave with a half-open other-axis selection."""

    target = _semantic_axis_index(obj, observable)
    selection = 1 - target
    selection_axis = obj.GetXaxis() if selection == 0 else obj.GetYaxis()
    first, last = _range_to_bins(selection_axis, selection_range)
    projection_name = name or f"{obj.GetName()}_{observable}_projection"
    if target == 0:
        result = obj.ProjectionX(projection_name, first, last, "e")
    else:
        result = obj.ProjectionY(projection_name, first, last, "e")
    result.SetDirectory(0)
    return result


def _clone_detached(obj):
    clone = obj.Clone()
    if hasattr(clone, "SetDirectory"):
        clone.SetDirectory(0)
    return clone


def _range_to_bins(axis, value_range):
    if value_range is None:
        return 1, axis.GetNbins()

    low, high = value_range
    # Interpret configured intervals as [low, high), with a small inward
    # offset matching the analysis notebooks and ROOT macros.
    low_bin = max(1, axis.FindBin(low + 0.001))
    high_bin = min(axis.GetNbins(), axis.FindBin(high - 0.001))
    if low_bin > high_bin:
        low_bin, high_bin = high_bin, low_bin
    return low_bin, high_bin


def project_histogram(obj, spec: Mapping[str, Any] | None = None, *, name: str | None = None):
    """
    Convert a ROOT object to a 1D histogram when a projection is requested.

    Supported inputs:
    - TH1-like objects are cloned and returned directly.
    - TH2 and TH3 objects can be projected on x, y, or z with optional axis ranges.
    - THnSparse objects can be projected on a named axis index after applying axis ranges.

    The caller provides `spec` only for objects that need a projection.
    """

    if obj.InheritsFrom("THnSparse"):
        if spec is None:
            return _clone_detached(obj)

        axis_index = spec.get("axis_index")
        if axis_index is None:
            raise ValueError("THnSparse projection requires 'axis_index' in spec")

        axis_ranges = spec.get("axis_ranges", {})
        for idx, value_range in axis_ranges.items():
            axis = obj.GetAxis(idx)
            low_bin, high_bin = _range_to_bins(axis, value_range)
            axis.SetRange(low_bin, high_bin)

        return _clone_detached(obj.Projection(axis_index))

    if obj.InheritsFrom("TH3"):
        if spec is None:
            return _clone_detached(obj)

        axis = spec.get("axis", "x").lower()
        x_range = spec.get("x_range")
        y_range = spec.get("y_range")
        z_range = spec.get("z_range")

        xbin1, xbin2 = _range_to_bins(obj.GetXaxis(), x_range)
        ybin1, ybin2 = _range_to_bins(obj.GetYaxis(), y_range)
        zbin1, zbin2 = _range_to_bins(obj.GetZaxis(), z_range)

        proj_name = name or f"{obj.GetName()}_proj"
        if axis == "x":
            return _clone_detached(obj.ProjectionX(proj_name, ybin1, ybin2, zbin1, zbin2))
        if axis == "y":
            return _clone_detached(obj.ProjectionY(proj_name, xbin1, xbin2, zbin1, zbin2))
        if axis == "z":
            return _clone_detached(obj.ProjectionZ(proj_name, xbin1, xbin2, ybin1, ybin2))
        raise ValueError(f"Unsupported TH3 projection axis: {axis!r}")

    if obj.InheritsFrom("TH2"):
        if spec is None:
            return _clone_detached(obj)

        axis = spec.get("axis", "x").lower()
        x_range = spec.get("x_range")
        y_range = spec.get("y_range")

        xbin1, xbin2 = _range_to_bins(obj.GetXaxis(), x_range)
        ybin1, ybin2 = _range_to_bins(obj.GetYaxis(), y_range)

        proj_name = name or f"{obj.GetName()}_proj"
        if axis == "x":
            return _clone_detached(obj.ProjectionX(proj_name, ybin1, ybin2))
        if axis == "y":
            return _clone_detached(obj.ProjectionY(proj_name, xbin1, xbin2))
        raise ValueError(f"Unsupported TH2 projection axis: {axis!r}")

    return _clone_detached(obj)
