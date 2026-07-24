"""Dijet eta-shape and forward/backward closure preparation with ROOT."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from hist_analysis.python.histogram_io import load_histogram
from hist_analysis.python.histogram_ops import normalize_histogram
from hist_analysis.python.projections import project_semantic_th2


@dataclass(frozen=True)
class DijetClosureCurve:
    """ROOT keys and display label for one distribution compared with gen."""

    label: str
    cm_key: str
    forward_key: str
    backward_key: str


def _key(template: str, eta_cut_index: int) -> str:
    try:
        return template.format(eta_cut_index=eta_cut_index)
    except (IndexError, KeyError) as error:
        raise ValueError(
            "Histogram-key templates must use the {eta_cut_index} field"
        ) from error


def _project_eta(filename: Path, key: str, ptave_range: tuple[float, float],
                 *, name: str, rebin_eta: int):
    source = load_histogram(str(filename), key)
    if not source.InheritsFrom("TH2"):
        raise TypeError(f"{key!r} in {filename} is not a TH2")
    source.Rebin2D(1, rebin_eta)
    return project_semantic_th2(
        source, "eta", ptave_range, name=name,
    )


def _root_divide(numerator, denominator, *, name: str, option: str):
    """Clone and divide through ROOT.TH1::Divide."""

    result = numerator.Clone(name)
    result.SetDirectory(0)
    result.Divide(numerator, denominator, 1.0, 1.0, option)
    return result


def _validate_common_binning(histograms, *, observable: str) -> None:
    """Require identical ROOT bin edges without rebuilding the histograms."""

    reference_label, reference = next(iter(histograms.items()))
    reference_axis = reference.GetXaxis()
    for label, histogram in histograms.items():
        axis = histogram.GetXaxis()
        same_edges = histogram.GetNbinsX() == reference.GetNbinsX() and all(
            abs(
                axis.GetBinLowEdge(bin_index)
                - reference_axis.GetBinLowEdge(bin_index)
            ) < 1e-9
            for bin_index in range(1, histogram.GetNbinsX() + 2)
        )
        if not same_edges:
            raise ValueError(
                f"Incompatible {observable} binning for {label!r} and "
                f"{reference_label!r}; ROOT projections were not altered"
            )


def build_dijet_gen_comparisons(
    filename: str | Path,
    curves: Iterable[DijetClosureCurve],
    *,
    eta_cut_index: int,
    ptave_range: tuple[float, float],
    nominal: str = "Gen",
    rebin_eta: int = 2,
    normalization: str = "integral",
    ratio_option: str = "B",
):
    """Build normalized eta shapes and unnormalized forward/backward ratios.

    The pTave interval is interpreted as half-open ``[low, high)``. Full eta
    projections are normalized with ``TH1::Scale`` according to
    ``normalization`` (``none``, ``integral``, or ``bin_width``).
    Forward/backward ratios are never normalized and use ``TH1::Divide``
    directly; ``ratio_option='B'`` reproduces the macro's default
    binomial-error setting, while an empty string requests independent errors.

    Returns ``(eta_shapes, fb_ratios, selected_keys)``. Both mappings retain the
    supplied curve order and are ready for ``plotting.draw_closure``, which uses
    ROOT division again for each curve's ratio to the nominal gen distribution.
    """

    if isinstance(rebin_eta, bool) or not isinstance(rebin_eta, int) or rebin_eta < 1:
        raise ValueError(f"rebin_eta must be a positive integer, got {rebin_eta!r}")
    if normalization not in {"none", "integral", "bin_width"}:
        raise ValueError(
            "normalization must be 'none', 'integral', or 'bin_width'"
        )
    if ratio_option not in {"", "B"}:
        raise ValueError("ratio_option must be '' or 'B'")
    low, high = ptave_range
    if low >= high:
        raise ValueError("ptave_range must satisfy low < high")

    path = Path(filename)
    curve_list = tuple(curves)
    labels = [curve.label for curve in curve_list]
    if not curve_list:
        raise ValueError("At least one closure curve is required")
    if len(set(labels)) != len(labels):
        raise ValueError(f"Curve labels must be unique: {labels}")
    if nominal not in labels:
        raise ValueError(f"Nominal {nominal!r} is not among {labels}")

    full = {}
    forward = {}
    backward = {}
    selected_keys = {}
    for index, curve in enumerate(curve_list):
        keys = {
            "cm": _key(curve.cm_key, eta_cut_index),
            "forward": _key(curve.forward_key, eta_cut_index),
            "backward": _key(curve.backward_key, eta_cut_index),
        }
        full[curve.label] = _project_eta(
            path, keys["cm"], ptave_range,
            name=f"dijet_closure_cm_{index}", rebin_eta=rebin_eta,
        )
        forward[curve.label] = _project_eta(
            path, keys["forward"], ptave_range,
            name=f"dijet_closure_forward_{index}", rebin_eta=rebin_eta,
        )
        backward[curve.label] = _project_eta(
            path, keys["backward"], ptave_range,
            name=f"dijet_closure_backward_{index}", rebin_eta=rebin_eta,
        )
        selected_keys[curve.label] = keys

    _validate_common_binning(full, observable="full eta")
    _validate_common_binning(forward, observable="forward eta")
    _validate_common_binning(backward, observable="backward eta")
    for label in labels:
        _validate_common_binning(
            {"forward": forward[label], "backward": backward[label]},
            observable=f"forward/backward eta for {label}",
        )

    eta_shapes = {
        label: normalize_histogram(histogram, normalization)
        for label, histogram in full.items()
    }
    fb_ratios = {
        label: _root_divide(
            forward[label], backward[label],
            name=f"dijet_closure_fb_{index}", option=ratio_option,
        )
        for index, label in enumerate(labels)
    }
    return eta_shapes, fb_ratios, selected_keys
