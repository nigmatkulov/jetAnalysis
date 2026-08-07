"""Fit helpers for systematic-variation histograms."""

from __future__ import annotations

import math
from typing import Mapping, Sequence

import ROOT


def _same_binning(left, right) -> bool:
    if left.GetNbinsX() != right.GetNbinsX():
        return False
    return all(
        abs(
            left.GetXaxis().GetBinLowEdge(index)
            - right.GetXaxis().GetBinLowEdge(index)
        ) < 1e-9
        for index in range(1, left.GetNbinsX() + 2)
    )


def calculate_bin_by_bin_systematic(
    up_ratio,
    down_ratio,
    *,
    name: str,
    up_function=None,
    down_function=None,
    evaluation_range: tuple[float, float] | None = None,
):
    """Reproduce ``calculateSystUncrtBinByBin`` without smoothing.

    The inputs are variation/default ratio histograms and provide the output
    binning. When both functions are supplied, their values at each bin center
    replace the raw bin contents, matching the macro's ``useFit`` path. Each
    output bin is ``(|up - 1| + |down - 1|) / 2``. As in the main-branch macro,
    the result is set to zero when both evaluated values are below ``1e-6``.
    Every bin that overlaps ``evaluation_range`` is included. For a boundary
    bin whose center lies outside the range, the fit is evaluated at the nearest
    range boundary rather than beyond the accepted eta range. The detached output
    has zero bin errors because it stores the systematic estimate itself.
    """

    if not name:
        raise ValueError("name must be non-empty")
    if not _same_binning(up_ratio, down_ratio):
        raise ValueError("Up and Down ratio histograms must have identical binning")
    if (up_function is None) != (down_function is None):
        raise ValueError("up_function and down_function must be supplied together")
    if evaluation_range is not None:
        range_low, range_high = evaluation_range
        if (
            not math.isfinite(range_low)
            or not math.isfinite(range_high)
            or range_low >= range_high
        ):
            raise ValueError(
                "evaluation_range must contain finite values with low < high"
            )

    systematic = up_ratio.Clone(name)
    systematic.SetDirectory(0)
    systematic.Reset()
    for bin_index in range(1, up_ratio.GetNbinsX() + 1):
        bin_center = up_ratio.GetBinCenter(bin_index)
        evaluation_point = bin_center
        if evaluation_range is not None:
            bin_low = up_ratio.GetXaxis().GetBinLowEdge(bin_index)
            bin_high = up_ratio.GetXaxis().GetBinUpEdge(bin_index)
            if bin_high <= range_low or bin_low >= range_high:
                continue
            evaluation_point = min(max(bin_center, range_low), range_high)
        up = (
            up_function.Eval(evaluation_point)
            if up_function is not None else up_ratio.GetBinContent(bin_index)
        )
        down = (
            down_function.Eval(evaluation_point)
            if down_function is not None else down_ratio.GetBinContent(bin_index)
        )
        uncertainty = (abs(up - 1.0) + abs(down - 1.0)) / 2.0
        if up < 1e-6 and down < 1e-6:
            uncertainty = 0.0
        systematic.SetBinContent(bin_index, uncertainty)
        systematic.SetBinError(bin_index, 0.0)
    return systematic


def smooth_systematic_running_max(
    histogram,
    *,
    name: str,
    evaluation_range: tuple[float, float],
    smoothing_origin: float | None = None,
):
    """Return an outward-running-maximum clone within an accepted range.

    With no origin, this reproduces the nonnegative-eta branch used for JER F/B
    uncertainties on ``main``. With an origin, the ROOT bin containing that
    value starts the increasing-x branch, while the preceding bin starts the
    decreasing-x branch. Boundary bins that overlap ``evaluation_range`` are
    included even when their centers lie outside it. Fully disjoint bins are
    explicitly kept at zero so smoothing cannot exceed the analysis acceptance.
    """

    if not name:
        raise ValueError("name must be non-empty")
    range_low, range_high = evaluation_range
    if (
        not math.isfinite(range_low)
        or not math.isfinite(range_high)
        or range_low >= range_high
    ):
        raise ValueError(
            "evaluation_range must contain finite values with low < high"
        )
    if smoothing_origin is not None and (
        not math.isfinite(smoothing_origin)
        or not range_low <= smoothing_origin <= range_high
    ):
        raise ValueError("smoothing_origin must lie within evaluation_range")

    smoothed = histogram.Clone(name)
    smoothed.SetDirectory(0)
    accepted_bins = [
        index
        for index in range(1, smoothed.GetNbinsX() + 1)
        if smoothed.GetXaxis().GetBinUpEdge(index) > range_low
        and smoothed.GetXaxis().GetBinLowEdge(index) < range_high
    ]
    if not accepted_bins:
        raise ValueError("evaluation_range overlaps no histogram bins")

    if smoothing_origin is None:
        branches = (accepted_bins,)
    else:
        start_bin = smoothed.FindBin(smoothing_origin)
        accepted = set(accepted_bins)
        if start_bin not in accepted:
            raise ValueError("smoothing_origin maps outside accepted bins")
        branches = (
            [index for index in accepted_bins if index >= start_bin],
            [index for index in reversed(accepted_bins) if index < start_bin],
        )
    for branch in branches:
        if not branch:
            continue
        running_maximum = smoothed.GetBinContent(branch[0])
        for bin_index in branch[1:]:
            running_maximum = max(
                running_maximum, smoothed.GetBinContent(bin_index)
            )
            smoothed.SetBinContent(bin_index, running_maximum)
    accepted = set(accepted_bins)
    for bin_index in range(1, smoothed.GetNbinsX() + 1):
        if bin_index not in accepted:
            smoothed.SetBinContent(bin_index, 0.0)
            smoothed.SetBinError(bin_index, 0.0)
    return smoothed


def format_fit_summary_lines(
    summaries: Mapping[str, Mapping[str, object]],
    *,
    parameters_per_line: int = 2,
    precision: int = 4,
):
    """Format fitted parameters and chi-square/ndf for ROOT plot text."""

    if parameters_per_line < 1:
        raise ValueError("parameters_per_line must be positive")
    if precision < 1:
        raise ValueError("precision must be positive")

    formatted = {}
    for label, summary in summaries.items():
        parameters = tuple(summary["parameters"])
        lines = []
        for first in range(0, len(parameters), parameters_per_line):
            parameter_text = ", ".join(
                f"p_{{{index}}}={parameters[index]:.{precision}g}"
                for index in range(
                    first, min(first + parameters_per_line, len(parameters))
                )
            )
            prefix = f"{label}: " if first == 0 else ""
            lines.append(prefix + parameter_text)
        lines.append(
            f"#chi^{{2}}/ndf={summary['chi2']:.{precision}g}/{summary['ndf']}"
        )
        formatted[label] = tuple(lines)
    return formatted


def fit_histogram_variations(
    histograms: Mapping[str, object],
    *,
    formula: str,
    fit_range: tuple[float, float],
    name_prefix: str,
    fit_options: str = "RQS0",
    initial_values: Mapping[str, Sequence[float]] | None = None,
):
    """Fit variation histograms and return TF1 objects and numeric summaries.

    ``fit_options`` must contain ROOT's ``S`` option so fit status and parameter
    covariance are available through ``TFitResultPtr``. ROOT selects fit bins by
    their centers, so the numerical fit range is expanded to the centers of the
    first and last bins that overlap ``fit_range``. This includes partially
    accepted rebinned boundary bins. The TF1 drawing range is restored to the
    requested physical range after fitting. The returned summaries contain plain
    Python values suitable for later systematic-uncertainty use.
    """

    if not histograms:
        raise ValueError("No histograms supplied for fitting")
    low, high = fit_range
    if not math.isfinite(low) or not math.isfinite(high) or low >= high:
        raise ValueError("fit_range must contain finite values with low < high")
    if not formula:
        raise ValueError("formula must be non-empty")
    if "S" not in fit_options.upper():
        raise ValueError("fit_options must include 'S' to retain fit results")
    if initial_values is not None:
        unknown_labels = set(initial_values) - set(histograms)
        if unknown_labels:
            raise KeyError(
                f"Initial values have no matching histograms: {sorted(unknown_labels)}"
            )

    functions = {}
    summaries = {}
    for index, (label, histogram) in enumerate(histograms.items()):
        overlapping_bins = [
            bin_index
            for bin_index in range(1, histogram.GetNbinsX() + 1)
            if histogram.GetXaxis().GetBinUpEdge(bin_index) > low
            and histogram.GetXaxis().GetBinLowEdge(bin_index) < high
        ]
        if not overlapping_bins:
            raise ValueError(
                f"fit_range overlaps no histogram bins for {label!r}"
            )
        fit_low = histogram.GetBinCenter(overlapping_bins[0])
        fit_high = histogram.GetBinCenter(overlapping_bins[-1])
        if fit_low >= fit_high:
            raise ValueError(
                f"fit_range must overlap at least two bin centers for {label!r}"
            )
        function = ROOT.TF1(
            f"{name_prefix}_{index}", formula, fit_low, fit_high,
        )
        default_values = (1.0,) + (0.0,) * (function.GetNpar() - 1)
        parameter_values = tuple(
            (initial_values or {}).get(label, default_values)
        )
        if len(parameter_values) != function.GetNpar():
            raise ValueError(
                f"Initial values for {label!r} contain {len(parameter_values)} "
                f"parameters, but {formula!r} requires {function.GetNpar()}"
            )
        if not all(math.isfinite(value) for value in parameter_values):
            raise ValueError(f"Initial values for {label!r} must be finite")
        function.SetParameters(*parameter_values)
        function.SetLineWidth(2)

        fit_result = histogram.Fit(
            function, fit_options, "", fit_low, fit_high,
        )
        status = int(fit_result)
        if status != 0:
            raise RuntimeError(
                f"Fit failed with status {status} for {label!r} using {formula!r}"
            )

        functions[label] = function
        summaries[label] = {
            "formula": formula,
            "range": (low, high),
            "fit_bin_center_range": (fit_low, fit_high),
            "initial_parameters": parameter_values,
            "parameters": tuple(
                function.GetParameter(i) for i in range(function.GetNpar())
            ),
            "parameter_errors": tuple(
                function.GetParError(i) for i in range(function.GetNpar())
            ),
            "covariance": tuple(
                tuple(
                    fit_result.CovMatrix(row, column)
                    for column in range(function.GetNpar())
                )
                for row in range(function.GetNpar())
            ),
            "chi2": function.GetChisquare(),
            "ndf": function.GetNDF(),
            "probability": function.GetProb(),
            "status": status,
        }
        function.SetRange(low, high)
    return functions, summaries
