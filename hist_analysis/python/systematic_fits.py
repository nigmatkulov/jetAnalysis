"""Fit helpers for systematic-variation histograms."""

from __future__ import annotations

import math
from typing import Mapping, Sequence

import ROOT


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
    covariance are available through ``TFitResultPtr``. The returned summaries
    contain plain Python values suitable for later systematic-uncertainty use.
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
        function = ROOT.TF1(
            f"{name_prefix}_{index}", formula, low, high,
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
            function, fit_options, "", low, high,
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
    return functions, summaries
