"""Primary-vertex v_z reweighting: shapes, MC/Data fits, weights and closure.

Production convention (verified in ``processForestSimple.C``):

* ``hVzRaw`` (60 uniform bins in [-15, 15] cm) is filled with unit weight for
  every chain entry, after the beam-direction v_z flip and before
  ``isGoodEvent``.  MC flips v_z when ``isPbGoing``; data flips when
  ``!isPbGoing``.
* ``eventWeight`` multiplies the pPb MC generator weight by
  ``1 / fVzWeight.Eval(vz)``.  ``fVzWeight`` must therefore describe the
  MC/Data shape ratio; the multiplicative MC weight is its reciprocal.
* pPb and pp event selections require ``|vz| <= 15`` cm, so the function is
  never evaluated outside the histogram range in production.

This module contains no notebook state.  Every ROOT object it creates has a
caller-supplied unique name and is detached from ROOT directories.
"""

from __future__ import annotations

from array import array
from dataclasses import dataclass, field
import json
from pathlib import Path
import math
from typing import Iterable, Mapping, Sequence

import ROOT

ROOT.TH1.AddDirectory(False)

# Legacy single production function used for every MC sample before
# 2026-11-27 (old processForestSimple.C:4114-4115). Kept as the fit starting
# point and for old-vs-new comparisons; production now selects a
# generator/direction-specific set (vzWeightParameters in processForestSimple.C).
LEGACY_POL8_PARAMETERS = (
    0.856516, -0.0159813, 0.00436628, -0.00012862, 2.61129e-05,
    -4.16965e-07, 1.73711e-08, -3.11953e-09, 6.24993e-10,
)
PRODUCTION_POL8_PARAMETERS = LEGACY_POL8_PARAMETERS  # backward-compatible alias


# ------------------------------------------------------------------
# Histogram validation and normalization
# ------------------------------------------------------------------
@dataclass(frozen=True)
class BinningInfo:
    n_bins: int
    low: float
    high: float
    uniform: bool
    width: float | None


def binning_info(histogram) -> BinningInfo:
    axis = histogram.GetXaxis()
    widths = [axis.GetBinWidth(i) for i in range(1, histogram.GetNbinsX() + 1)]
    uniform = max(widths) - min(widths) <= 1e-9 * max(abs(w) for w in widths)
    return BinningInfo(histogram.GetNbinsX(), axis.GetXmin(), axis.GetXmax(),
                       uniform, widths[0] if uniform else None)


def same_binning(left, right, *, tolerance: float = 1e-9) -> bool:
    if left.GetNbinsX() != right.GetNbinsX():
        return False
    la, ra = left.GetXaxis(), right.GetXaxis()
    return all(abs(la.GetBinLowEdge(i) - ra.GetBinLowEdge(i)) <= tolerance
               for i in range(1, left.GetNbinsX() + 2))


def validate_histogram(histogram, *, label: str) -> dict:
    """Return a summary; raise if the histogram cannot be used."""
    if histogram is None or not histogram:
        raise KeyError(f"{label}: histogram missing")
    if not histogram.InheritsFrom("TH1") or histogram.GetDimension() != 1:
        raise TypeError(f"{label}: expected a 1D histogram")
    n = histogram.GetNbinsX()
    in_range = sum(histogram.GetBinContent(i) for i in range(1, n + 1))
    if histogram.GetEntries() <= 0 or in_range <= 0:
        raise ValueError(f"{label}: no entries in range")
    negative = [i for i in range(1, n + 1) if histogram.GetBinContent(i) < 0]
    if negative:
        raise ValueError(f"{label}: negative bin contents in bins {negative}")
    underflow, overflow = histogram.GetBinContent(0), histogram.GetBinContent(n + 1)
    total = in_range + underflow + overflow
    info = binning_info(histogram)
    return dict(label=label, entries=float(histogram.GetEntries()),
                in_range=float(in_range), underflow=float(underflow),
                overflow=float(overflow),
                flow_fraction=float((underflow + overflow) / total),
                n_bins=info.n_bins, low=info.low, high=info.high,
                uniform=info.uniform, bin_width=info.width,
                empty_bins=[i for i in range(1, n + 1)
                            if histogram.GetBinContent(i) == 0],
                mean=float(histogram.GetMean()), rms=float(histogram.GetRMS()))


def bin_range(histogram, value_range: tuple[float, float]) -> tuple[int, int]:
    """Inclusive in-range bins whose centres lie in ``value_range``.

    Flow bins are never included.  A range covering the whole axis selects
    bins 1..N.
    """
    low, high = value_range
    if not low < high:
        raise ValueError("range must satisfy low < high")
    axis = histogram.GetXaxis()
    bins = [i for i in range(1, histogram.GetNbinsX() + 1)
            if low <= axis.GetBinCenter(i) <= high]
    if not bins:
        raise ValueError(f"No bin centre lies in {value_range}")
    return bins[0], bins[-1]


def normalize_unit_area(histogram, name: str, *,
                        norm_range: tuple[float, float]):
    """Clone and divide by the sum of bin contents in ``norm_range``.

    The result holds bin fractions p_i = H_i / sum_j H_j (not densities).
    Bin errors are scaled by the same constant: the multinomial correlation
    induced by the shared denominator is neglected (relative effect
    O(1/sqrt(N_total))).  Flow bins are scaled as well but are not part of
    the denominator; bins outside ``norm_range`` keep their scaled values.
    """
    first, last = bin_range(histogram, norm_range)
    total = histogram.Integral(first, last)
    if total <= 0:
        raise ValueError(f"{histogram.GetName()}: non-positive normalization")
    result = histogram.Clone(name)
    result.SetDirectory(0)
    if result.GetSumw2N() == 0:
        result.Sumw2()
    result.Scale(1.0 / total)
    return result, float(total)


def density_from_fraction(fraction, name: str):
    """Convert bin fractions to a probability density (divide by width)."""
    result = fraction.Clone(name)
    result.SetDirectory(0)
    for i in range(0, result.GetNbinsX() + 2):
        width = result.GetXaxis().GetBinWidth(i) if 1 <= i <= result.GetNbinsX() else 1.0
        result.SetBinContent(i, fraction.GetBinContent(i) / width)
        result.SetBinError(i, fraction.GetBinError(i) / width)
    return result


def ratio_histogram(numerator, denominator, name: str,
                    *, bins: tuple[int, int] | None = None):
    """Independent-sample ratio with explicit invalid-bin bookkeeping.

    sigma_R / R = sqrt((sigma_a/a)^2 + (sigma_b/b)^2).  A bin is valid when
    both contents are positive and both errors are positive and finite.
    Invalid bins (including flow and bins outside ``bins``) are set to
    content 0 and error 0, and returned in the ``invalid`` list.
    """
    if not same_binning(numerator, denominator):
        raise ValueError(f"{name}: incompatible binning")
    result = numerator.Clone(name)
    result.SetDirectory(0)
    result.Reset("ICES")
    n = result.GetNbinsX()
    first, last = bins or (1, n)
    invalid = []
    for i in range(0, n + 2):
        a, b = numerator.GetBinContent(i), denominator.GetBinContent(i)
        ea, eb = numerator.GetBinError(i), denominator.GetBinError(i)
        ok = (first <= i <= last and a > 0 and b > 0 and ea > 0 and eb > 0
              and all(map(math.isfinite, (a, b, ea, eb))))
        if not ok:
            result.SetBinContent(i, 0.0)
            result.SetBinError(i, 0.0)
            if 1 <= i <= n:
                invalid.append(i)
            continue
        value = a / b
        result.SetBinContent(i, value)
        result.SetBinError(i, value * math.hypot(ea / a, eb / b))
    result.SetEntries(numerator.GetEntries())
    return result, invalid


# ------------------------------------------------------------------
# Fitting
# ------------------------------------------------------------------
@dataclass
class FitResult:
    name: str
    formula: str
    fit_range: tuple[float, float]
    status: int
    cov_status: int
    parameters: list[float]
    errors: list[float]
    covariance: list[list[float]]
    chi2: float
    ndf: int
    probability: float
    n_points: int
    function: object = field(repr=False)
    linear_parameters: list[float] | None = None
    max_linear_relative_difference: float | None = None

    @property
    def chi2_ndf(self) -> float:
        return self.chi2 / self.ndf if self.ndf > 0 else float("nan")

    def summary(self) -> dict:
        return dict(name=self.name, formula=self.formula,
                    fit_range=list(self.fit_range), status=self.status,
                    cov_status=self.cov_status, parameters=self.parameters,
                    errors=self.errors, covariance=self.covariance,
                    chi2=self.chi2, ndf=self.ndf, chi2_ndf=self.chi2_ndf,
                    probability=self.probability, n_points=self.n_points,
                    linear_parameters=self.linear_parameters,
                    max_linear_relative_difference=self.max_linear_relative_difference)


def fit_points(ratio, fit_range: tuple[float, float]):
    """Valid (centre, value, error) triples whose bin lies inside the range.

    A bin is used when its full extent lies within ``fit_range`` (the
    default [-15.1, 15.1] contains all 60 bins) and its error is positive.
    """
    low, high = fit_range
    axis = ratio.GetXaxis()
    points = []
    for i in range(1, ratio.GetNbinsX() + 1):
        if axis.GetBinLowEdge(i) < low - 1e-9 or axis.GetBinUpEdge(i) > high + 1e-9:
            continue
        error = ratio.GetBinError(i)
        if error > 0 and math.isfinite(error):
            points.append((axis.GetBinCenter(i), ratio.GetBinContent(i), error, i))
    return points


def initial_parameters(order: int, supplied: Sequence[float] | None) -> list[float]:
    """Supplied values (truncated/padded with zeros) for polynomial ``order``."""
    supplied = list(supplied or ())
    return [float(supplied[k]) if k < len(supplied) else 0.0 for k in range(order + 1)]


def _weighted_polynomial_least_squares(points, order: int) -> list[float]:
    """Exact chi2 minimum for a polynomial, solved in a scaled variable.

    Used only as an independent cross-check of the MINUIT fit; it avoids the
    ill conditioning of raw powers of v_z by fitting in u = x / s.
    """
    import numpy as np
    x = np.array([p[0] for p in points])
    y = np.array([p[1] for p in points])
    w = 1.0 / np.array([p[2] for p in points])
    scale = max(1.0, float(np.max(np.abs(x))))
    u = x / scale
    design = np.vander(u, order + 1, increasing=True) * w[:, None]
    coefficients, *_ = np.linalg.lstsq(design, y * w, rcond=None)
    return [float(c / scale ** k) for k, c in enumerate(coefficients)]


def fit_polynomial(ratio, *, name: str, order: int,
                   fit_range: tuple[float, float],
                   start: Sequence[float] | None = None) -> FitResult:
    """Chi2 fit of ``pol<order>`` to the valid ratio points with MINUIT.

    ROOT fits ``polN`` with its linear fitter unless option ``F`` is given;
    ``F`` is used so that the supplied coefficients act as initial values
    (they are never fixed).  The exact linear least-squares solution is
    computed independently and the maximum relative difference of the two
    functions on a 0.05 cm grid is recorded.
    """
    if order < 0:
        raise ValueError("order must be non-negative")
    points = fit_points(ratio, fit_range)
    if len(points) <= order + 1:
        raise ValueError(f"{name}: {len(points)} valid points for pol{order}")
    old = ROOT.gROOT.GetFunction(name)
    if old:
        ROOT.gROOT.GetListOfFunctions().Remove(old)
    function = ROOT.TF1(name, f"pol{order}", *fit_range)
    function.SetNpx(1000)
    linear = _weighted_polynomial_least_squares(points, order)
    initial = initial_parameters(order, start) if start is not None else linear
    for k, value in enumerate(initial):
        function.SetParameter(k, value)
        function.SetParError(k, max(abs(value) * 0.1, 1e-12))
    graph = ROOT.TGraphErrors(len(points))
    graph.SetName(f"g_{name}_fitpoints")
    for index, (x, y, ey, _) in enumerate(points):
        graph.SetPoint(index, x, y)
        graph.SetPointError(index, 0.0, ey)
    fit = graph.Fit(function, "S Q N R F EX0")
    status = int(fit.Status())
    cov_status = int(fit.CovMatrixStatus())
    if status != 0:
        # Retry once starting at the exact solution; record it honestly.
        for k, value in enumerate(linear):
            function.SetParameter(k, value)
        fit = graph.Fit(function, "S Q N R F EX0")
        status = int(fit.Status())
        cov_status = int(fit.CovMatrixStatus())
    n = order + 1
    covariance = [[float(fit.CovMatrix(i, j)) for j in range(n)] for i in range(n)]
    parameters = [float(fit.Parameter(k)) for k in range(n)]
    errors = [float(fit.ParError(k)) for k in range(n)]
    for k in range(n):
        function.SetParameter(k, parameters[k])
        function.SetParError(k, errors[k])
    function.SetChisquare(float(fit.Chi2()))
    function.SetNDF(int(fit.Ndf()))
    reference = ROOT.TF1(f"{name}_linear_check", f"pol{order}", *fit_range)
    for k, value in enumerate(linear):
        reference.SetParameter(k, value)
    grid = _grid(fit_range[0], fit_range[1], 0.05)
    difference = max(abs(function.Eval(x) - reference.Eval(x)) / max(abs(reference.Eval(x)), 1e-12)
                     for x in grid)
    return FitResult(name=name, formula=f"pol{order}", fit_range=tuple(fit_range),
                     status=status, cov_status=cov_status, parameters=parameters,
                     errors=errors, covariance=covariance, chi2=float(fit.Chi2()),
                     ndf=int(fit.Ndf()), probability=float(fit.Prob()),
                     n_points=len(points), function=function,
                     linear_parameters=linear,
                     max_linear_relative_difference=float(difference))


def fit_function(ratio, *, name: str, formula: str, fit_range, start: Sequence[float]):
    """Generic MINUIT chi2 fit, used for diagnostic alternative shapes."""
    points = fit_points(ratio, fit_range)
    function = ROOT.TF1(name, formula, *fit_range)
    function.SetNpx(1000)
    for k, value in enumerate(start):
        function.SetParameter(k, value)
    graph = ROOT.TGraphErrors(len(points))
    for index, (x, y, ey, _) in enumerate(points):
        graph.SetPoint(index, x, y)
        graph.SetPointError(index, 0.0, ey)
    fit = graph.Fit(function, "S Q N R EX0")
    n = function.GetNpar()
    return FitResult(name=name, formula=formula, fit_range=tuple(fit_range),
                     status=int(fit.Status()), cov_status=int(fit.CovMatrixStatus()),
                     parameters=[float(fit.Parameter(k)) for k in range(n)],
                     errors=[float(fit.ParError(k)) for k in range(n)],
                     covariance=[[float(fit.CovMatrix(i, j)) for j in range(n)] for i in range(n)],
                     chi2=float(fit.Chi2()), ndf=int(fit.Ndf()), probability=float(fit.Prob()),
                     n_points=len(points), function=function)


def _grid(low: float, high: float, step: float) -> list[float]:
    count = int(round((high - low) / step))
    return [low + (high - low) * k / count for k in range(count + 1)]


def stability_checks(function, *, check_range: tuple[float, float] = (-15.0, 15.0),
                     step: float = 0.01, ratio=None,
                     fit_range: tuple[float, float] | None = None) -> dict:
    """Positivity, curvature sign changes (oscillations) and edge behaviour.

    ``curvature_sign_changes`` counts sign changes of f'' on the grid,
    ignoring |f''| below 1e-6 of its maximum.  ``extremum_count`` counts local
    extrema of f.  When ``ratio`` is given, pulls (R_i - f(x_i)) / sigma_i over
    the valid fit points are returned as well.
    """
    grid = _grid(*check_range, step)
    values = [function.Eval(x) for x in grid]
    second = [function.Derivative2(x) for x in grid]
    threshold = 1e-6 * max(abs(v) for v in second)
    signs = [math.copysign(1, v) for v in second if abs(v) > threshold]
    curvature_changes = sum(1 for a, b in zip(signs, signs[1:]) if a != b)
    first = [function.Derivative(x) for x in grid]
    # Ignore exact zeros so an extremum that falls on a grid point counts once.
    nonzero = [v for v in first if v != 0.0]
    extrema = sum(1 for a, b in zip(nonzero, nonzero[1:]) if a * b < 0)
    result = dict(
        min_value=float(min(values)), max_value=float(max(values)),
        argmin=float(grid[values.index(min(values))]),
        all_positive=all(v > 0 and math.isfinite(v) for v in values),
        curvature_sign_changes=int(curvature_changes),
        extremum_count=int(extrema),
        edge_values=dict(low=float(values[0]), high=float(values[-1])),
        edge_slopes=dict(low=float(first[0]), high=float(first[-1])),
    )
    if ratio is not None:
        pulls = [(i, x, (y - function.Eval(x)) / ey)
                 for x, y, ey, i in fit_points(ratio, fit_range or check_range)]
        result["pulls"] = [dict(bin=i, x=x, pull=p) for i, x, p in pulls]
        result["max_abs_pull"] = float(max(abs(p) for *_, p in pulls))
        result["edge_pulls"] = [p for i, x, p in pulls[:2] + pulls[-2:]]
    return result


def f_test(chi2_low: float, ndf_low: int, chi2_high: float, ndf_high: int) -> dict:
    """Nested-model F-test for adding (ndf_low - ndf_high) parameters."""
    extra = ndf_low - ndf_high
    if extra <= 0 or ndf_high <= 0 or chi2_high <= 0:
        return dict(F=float("nan"), p_value=float("nan"))
    statistic = ((chi2_low - chi2_high) / extra) / (chi2_high / ndf_high)
    p_value = ROOT.TMath.FDistI(statistic, extra, ndf_high) if statistic > 0 else 0.0
    return dict(F=float(statistic), p_value=float(1.0 - p_value))


def pull_histogram(ratio, function, name: str, *, fit_range):
    """(R_i - f(x_i)) / sigma_i for valid fit bins; other bins are zero."""
    pulls = ratio.Clone(name)
    pulls.SetDirectory(0)
    pulls.Reset("ICES")
    for x, y, ey, i in fit_points(ratio, fit_range):
        pulls.SetBinContent(i, (y - function.Eval(x)) / ey)
        pulls.SetBinError(i, 1.0)
    return pulls


# ------------------------------------------------------------------
# Weights and closure
# ------------------------------------------------------------------
def weight_function(fit_function, name: str, *,
                    valid_range: tuple[float, float] = (-15.0, 15.0)):
    """TF1 for the multiplicative MC weight w(v_z) = 1 / f(v_z).

    The returned function is defined only on ``valid_range`` (the production
    |v_z| <= 15 cm acceptance); it is never extrapolated.
    """
    order = fit_function.GetNpar() - 1
    terms = " + ".join(f"[{k}]*pow(x,{k})" if k else "[0]" for k in range(order + 1))
    weight = ROOT.TF1(name, f"1.0/({terms})", *valid_range)
    weight.SetNpx(1000)
    for k in range(order + 1):
        weight.SetParameter(k, fit_function.GetParameter(k))
    return weight


def bin_averaged_weight(fit_function, low: float, high: float, *, samples: int = 41) -> float:
    """Average of 1/f over [low, high] (Simpson's rule, odd ``samples``)."""
    if samples % 2 == 0:
        samples += 1
    step = (high - low) / (samples - 1)
    total = 0.0
    for k in range(samples):
        coefficient = 1 if k in (0, samples - 1) else (4 if k % 2 else 2)
        total += coefficient / fit_function.Eval(low + k * step)
    return total * step / 3.0 / (high - low)


def reweight_histogram(histogram, fit_function, name: str, *,
                       bins: tuple[int, int] | None = None, mode: str = "average"):
    """Multiply each in-range bin by 1/f (bin average or bin centre).

    This is the histogram-level approximation of the event-level weight; it
    assumes the MC density is approximately uniform within each 0.5 cm bin.
    Bins outside ``bins`` and flow bins are set to zero (no extrapolation).
    Errors scale with the same factor (weight uncertainty not included).
    """
    if mode not in {"average", "centre"}:
        raise ValueError("mode must be 'average' or 'centre'")
    result = histogram.Clone(name)
    result.SetDirectory(0)
    n = result.GetNbinsX()
    first, last = bins or (1, n)
    axis = result.GetXaxis()
    for i in range(0, n + 2):
        if not first <= i <= last:
            result.SetBinContent(i, 0.0)
            result.SetBinError(i, 0.0)
            continue
        factor = (bin_averaged_weight(fit_function, axis.GetBinLowEdge(i), axis.GetBinUpEdge(i))
                  if mode == "average" else 1.0 / fit_function.Eval(axis.GetBinCenter(i)))
        result.SetBinContent(i, histogram.GetBinContent(i) * factor)
        result.SetBinError(i, histogram.GetBinError(i) * factor)
    return result


def agreement(mc, data, *, bins: tuple[int, int]) -> dict:
    """Shape agreement between two unit-normalized histograms.

    chi2 = sum (m_i - d_i)^2 / (s_m^2 + s_d^2) over bins with a positive
    combined error; ndf = number of bins - 1 (one normalization constraint).
    Bins are treated as independent; normalization correlations are ignored.
    """
    chi2 = 0.0
    used = 0
    l1 = 0.0
    max_relative = 0.0
    max_pull = 0.0
    for i in range(bins[0], bins[1] + 1):
        m, d = mc.GetBinContent(i), data.GetBinContent(i)
        variance = mc.GetBinError(i) ** 2 + data.GetBinError(i) ** 2
        l1 += abs(m - d)
        if d > 0:
            max_relative = max(max_relative, abs(m / d - 1.0))
        if variance > 0:
            chi2 += (m - d) ** 2 / variance
            max_pull = max(max_pull, abs(m - d) / math.sqrt(variance))
            used += 1
    ndf = used - 1
    return dict(chi2=float(chi2), ndf=int(ndf),
                chi2_ndf=float(chi2 / ndf) if ndf > 0 else float("nan"),
                sum_abs_difference=float(l1), max_abs_ratio_minus_one=float(max_relative),
                max_abs_pull=float(max_pull),
                mean_mc=float(mc.GetMean()), mean_data=float(data.GetMean()),
                rms_mc=float(mc.GetRMS()), rms_data=float(data.GetRMS()))


# ------------------------------------------------------------------
# Export helpers
# ------------------------------------------------------------------
def cpp_snippet(parameters: Sequence[float], *, order: int, fit_range,
                comment_lines: Iterable[str] = (), name: str = "fVzWeight") -> str:
    """Drop-in replacement for processForestSimple.C:4114-4115."""
    values = ",\n".join(f"    {value:.17g}" for value in parameters)
    header = "\n".join(f"// {line}" for line in comment_lines)
    return (
        f"{header}\n"
        "// fVzWeight is the fitted MC/Data shape ratio f(vz).\n"
        "// The MC event weight DIVIDES by it, as eventWeight() already does:\n"
        "//     double vzWeight = 1. / fVzWeight.Eval(vz);   // w_MC = Data/MC\n"
        "// Valid for |vz| <= 15 cm (the event selection); do not extrapolate.\n"
        f"auto {name} = std::make_unique<TF1>(\n"
        f"    \"{name}\", \"pol{order}\", {fit_range[0]:g}, {fit_range[1]:g}\n"
        ");\n"
        f"{name}->SetParameters(\n{values}\n);\n"
    )


def parse_cpp_parameters(snippet: str) -> list[float]:
    """Read back the SetParameters list from :func:`cpp_snippet` output."""
    body = snippet.split("SetParameters(", 1)[1].split(");", 1)[0]
    return [float(token) for token in body.replace("\n", " ").split(",") if token.strip()]


def to_json(value) -> str:
    return json.dumps(value, sort_keys=True, indent=2, default=float)


def histogram_from_values(name: str, edges: Sequence[float], values: Sequence[float],
                          errors: Sequence[float] | None = None):
    """Small helper for tests and event-level fills."""
    histogram = ROOT.TH1D(name, name, len(edges) - 1, array("d", edges))
    histogram.SetDirectory(0)
    histogram.Sumw2()
    for i, value in enumerate(values, 1):
        histogram.SetBinContent(i, value)
        histogram.SetBinError(i, math.sqrt(value) if errors is None else errors[i - 1])
    histogram.SetEntries(sum(values))
    return histogram


# ------------------------------------------------------------------
# Direction mixtures
# ------------------------------------------------------------------
def mixture(components: Mapping[str, object], fractions: Mapping[str, float], name: str):
    """sum_k f_k * H_k for unit-normalized H_k; fractions must sum to one.

    Used for the luminosity-matched combined MC: each MC direction enters
    with the event fraction of the data direction it is paired with.
    Errors add in quadrature (independent samples).
    """
    total = sum(fractions.values())
    if abs(total - 1.0) > 1e-12:
        raise ValueError(f"fractions sum to {total}, expected 1")
    histograms = [components[k] for k in fractions]
    reference = histograms[0]
    if not all(same_binning(reference, h) for h in histograms[1:]):
        raise ValueError("mixture components have different binning")
    result = reference.Clone(name)
    result.SetDirectory(0)
    result.Reset("ICES")
    for i in range(0, result.GetNbinsX() + 2):
        value = sum(fractions[k] * components[k].GetBinContent(i) for k in fractions)
        error = math.sqrt(sum((fractions[k] * components[k].GetBinError(i)) ** 2 for k in fractions))
        result.SetBinContent(i, value)
        result.SetBinError(i, error)
    result.SetEntries(sum(components[k].GetEntries() for k in fractions))
    return result


# ------------------------------------------------------------------
# Event-level inputs (HiForest), mirroring processForestSimple.C
# ------------------------------------------------------------------
# pT-hat sample ranges, half-open [low, high) (processForestSimple.C:2321-2331).
PTHAT_RANGES = {15: (15., 30.), 30: (30., 50.), 50: (50., 80.), 80: (80., 120.),
                120: (120., 170.), 170: (170., 220.), 220: (220., 280.),
                280: (280., 370.), 370: (370., 460.), 460: (460., 540.),
                540: (540., 8160.)}
# pPb skim filters of isGoodEvent (processForestSimple.C:2379-2385), nominal
# vertex filter pVertexFilterCutdz1p0 (vertexFilterSelection = 0).
PPB_FILTER_EXPRESSION = ("pBeamScrapingFilter == 1 && pPAprimaryVertexFilter == 1 && "
                         "HBHENoiseFilterResultRun2Loose == 1 && phfCoincFilter == 1 && "
                         "pVertexFilterCutdz1p0 == 1")


def pthat_sample_from_name(path) -> int | None:
    import re
    match = re.search(r"pthat(\d+)", str(path), re.I)
    return int(match.group(1)) if match else None


def event_level_histograms(files: Sequence, *, flip: bool, name: str, edges,
                           selection: str, weight_function=None, pp: bool = False) -> dict:
    """Fill v_z histograms from HiForest files with RDataFrame (no IMT).

    ``flip`` reproduces processForestSimple.C:4299-4301 for MC (v_z -> -v_z
    for Pb-going).  ``selection`` is ``'raw'`` (all entries, as hVzRaw) or
    ``'selected'`` (isGoodEvent: |v_z|<=15, pPb filters, pT-hat range of the
    file's sample; pp: |v_z|<=15 only).  The overweight-event rejection that
    precedes hVzUnweighted needs jets and is NOT reproduced, so 'selected'
    is a superset of the hVzUnweighted population.

    Each file yields an even-entry and an odd-entry histogram (entry parity
    within the file, deterministic because implicit MT is not enabled here).
    If ``weight_function`` (the MC/Data fit f) is given, the odd half is also
    filled with the per-event weight 1/f(v_z) for |v_z| <= 15 cm and zero
    outside (no extrapolation; those entries only reach the flow bins).
    """
    if ROOT.IsImplicitMTEnabled():
        raise RuntimeError("Disable ROOT implicit MT: the parity split needs ordered entries")
    if selection not in {"raw", "selected"}:
        raise ValueError("selection must be 'raw' or 'selected'")
    namespace = None
    if weight_function is not None:
        import uuid
        # Cling cannot redefine a namespace: use a fresh one per call.
        namespace = f"vzrw_{uuid.uuid4().hex}"
        order = weight_function.GetNpar() - 1
        coefficients = ", ".join(repr(weight_function.GetParameter(k)) for k in range(order + 1))
        ROOT.gInterpreter.Declare(
            f"namespace {namespace} {{ double f(double x) {{ const double c[] = {{{coefficients}}};"
            f" double r = 0; for (int k = {order}; k >= 0; --k) r = r * x + c[k]; return r; }} }}")
    model = ROOT.RDF.TH1DModel(f"{name}_model", "", len(edges) - 1, array("d", edges))
    keys = ("all", "even", "odd") + (("odd_weighted",) if weight_function is not None else ())
    totals = {}
    entries = 0
    # pp has no per-file pT-hat range: read all files as one chain (faster).
    groups = [list(files)] if pp else [[path] for path in files]
    for group in groups:
        path = group[0]
        chain = ROOT.TChain("hiEvtAnalyzer/HiTree")
        for item in group:
            chain.Add(str(item))
        friend = None
        if selection == "selected" and not pp:
            friend = ROOT.TChain("skimanalysis/HltTree")
            for item in group:
                friend.Add(str(item))
            chain.AddFriend(friend)
        frame = ROOT.RDataFrame(chain)
        frame = frame.Define("vzs", "-vz" if flip else "vz")
        if selection == "selected":
            cut = "std::abs(vzs) <= 15."
            if not pp:
                sample = pthat_sample_from_name(path)
                if sample not in PTHAT_RANGES:
                    raise ValueError(f"Cannot determine pT-hat sample of {path}")
                low, high = PTHAT_RANGES[sample]
                cut += f" && pthat >= {low} && pthat < {high} && {PPB_FILTER_EXPRESSION}"
            frame = frame.Filter(cut)
        frame = frame.Define("parity", "rdfentry_ % 2")
        results = {"all": frame.Histo1D(model, "vzs"),
                   "even": frame.Filter("parity == 0").Histo1D(model, "vzs"),
                   "odd": frame.Filter("parity == 1").Histo1D(model, "vzs")}
        if weight_function is not None:
            results["odd_weighted"] = (frame.Filter("parity == 1")
                                       .Define("w", f"std::abs(vzs) <= 15. ? 1.0 / {namespace}::f(vzs) : 0.0")
                                       .Histo1D(model, "vzs", "w"))
        entries += int(chain.GetEntries())
        per_file = totals.setdefault("per_file", {})
        single = results["all"].GetValue().Clone(f"{name}_file{len(per_file)}")
        single.SetDirectory(0)
        per_file[Path(path).name if len(group) == 1 else f'{len(group)} files'] = single
        for key in keys:
            histogram = results[key].GetValue()
            if key in totals:
                totals[key].Add(histogram)
            else:
                totals[key] = histogram.Clone(f"{name}_{key}")
                totals[key].SetDirectory(0)
        del frame, chain, friend
    totals["n_tree_entries"] = entries
    return totals
