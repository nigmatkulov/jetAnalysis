"""Reusable PyROOT helpers for flattened dijet pTave-eta unfolding.

The analysis uses two coordinates packed into one ROOT bin number.  Each pTave
interval owns one consecutive block of eta bins.  The response matrix therefore
retains migrations in both pTave and eta even though RooUnfold sees a 1D vector.

The factorized correction implemented here is

    inclusive reco --(multiply by purity)--> matched reco
                   --(unfold migrations)--> matched truth
                   --(divide by efficiency)--> inclusive truth.

Purity and efficiency are derived from the same training response marginals.
They are treated as fixed response inputs; their MC uncertainty must be assessed
with response toys or systematic variations.
"""

from __future__ import annotations

from array import array
from contextlib import contextmanager
from dataclasses import dataclass
import json
import hashlib
import math
import os
from pathlib import Path
import random
import struct
from typing import Mapping, Sequence

import ROOT

from hist_analysis.python.histogram_io import load_object, open_root_file
from hist_analysis.python.projections import project_semantic_th2


@contextmanager
def _silence_native_stdout(enabled: bool):
    """Temporarily silence C/C++ writes to stdout when explicitly requested."""

    if not enabled:
        yield
        return
    saved_stdout = os.dup(1)
    null_stdout = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(null_stdout, 1)
        yield
    finally:
        os.dup2(saved_stdout, 1)
        os.close(null_stdout)
        os.close(saved_stdout)


@dataclass(frozen=True)
class FlattenedBinning:
    """Mapping between pT-block/ROOT eta bins and one-based global ROOT bins."""

    pt_intervals: tuple[tuple[float, float], ...]
    n_eta_bins: int

    def __post_init__(self) -> None:
        if not self.pt_intervals:
            raise ValueError("At least one pT interval is required")
        if self.n_eta_bins <= 0:
            raise ValueError("n_eta_bins must be positive")
        previous_high = None
        for low, high in self.pt_intervals:
            if low >= high:
                raise ValueError(f"Invalid half-open pT interval {(low, high)}")
            if previous_high is not None and low < previous_high:
                raise ValueError("pT intervals must be ordered and non-overlapping")
            previous_high = high

    @property
    def n_pt_bins(self) -> int:
        return len(self.pt_intervals)

    @property
    def n_global_bins(self) -> int:
        return self.n_pt_bins * self.n_eta_bins

    def global_bin(self, pt_index: int, eta_bin: int) -> int:
        if not 0 <= pt_index < self.n_pt_bins:
            raise IndexError(f"pT index {pt_index} is outside [0, {self.n_pt_bins})")
        if not 1 <= eta_bin <= self.n_eta_bins:
            raise IndexError(f"ROOT eta bin {eta_bin} is outside [1, {self.n_eta_bins}]")
        return pt_index * self.n_eta_bins + eta_bin

    def indices(self, global_bin: int) -> tuple[int, int]:
        if not 1 <= global_bin <= self.n_global_bins:
            raise IndexError(
                f"Global ROOT bin {global_bin} is outside [1, {self.n_global_bins}]"
            )
        zero_based = global_bin - 1
        return divmod(zero_based, self.n_eta_bins)[0], zero_based % self.n_eta_bins + 1


@dataclass(frozen=True)
class UnfoldingInputKeys:
    truth: str
    measured: str
    response: str
    miss: str
    fake: str
    classification: str | None = None


@dataclass
class UnfoldingInputs:
    truth: object
    measured: object
    response: object
    miss: object
    fake: object
    classification: object | None
    keys: UnfoldingInputKeys


@dataclass
class ResponseDiagnostics:
    matched_truth: object
    matched_measured: object
    effective_miss: object
    effective_fake: object
    boundary_miss: object | None
    boundary_fake: object | None


@dataclass(frozen=True)
class ResponseAccounting:
    """Numerical checks and summary fractions for a flattened response."""

    max_truth_residual: float
    max_measured_residual: float
    max_truth_relative_residual: float
    max_measured_relative_residual: float
    zero_efficiency_truth_bins: tuple[int, ...]
    empty_measured_bins: tuple[int, ...]
    response_sparsity: float
    global_efficiency: float
    global_fake_fraction: float


@dataclass
class RooUnfoldResponseBundle:
    response: object
    scaled_truth: object
    scaled_measured: object
    scaled_matrix: object
    scale: float


@dataclass
class UnfoldingResult:
    histogram: object
    covariance: object
    algorithm: object
    scaled_measured: object


@dataclass
class FactorizedCorrectionInputs:
    """Objects needed by the factorized correction.

    ``measured_signal`` is the input data/reco histogram multiplied by the
    training purity.  ``purity`` lives in measured bins, whereas ``efficiency``
    lives in truth bins; keeping them as separate histograms prevents accidental
    use on the wrong side of the response matrix.
    """

    measured_signal: object
    purity: object
    efficiency: object


@dataclass
class ForwardFoldResult:
    """Forward-folded truth and the fake contribution used to form pseudo-data."""

    histogram: object
    folded_signal: object
    fake: object
    apply_to_truth_includes_fakes: bool


@dataclass(frozen=True)
class ClosureMetrics:
    """Compact closure measures that remain meaningful for sparse spectra."""

    integral_ratio: float
    relative_l1: float
    mean_absolute_relative: float
    compared_bins: int


@dataclass
class ToyMetricResult:
    """Per-bin and aggregate toy metrics for one Bayesian iteration count."""

    iterations: int
    mean: object
    bias: object
    bias_squared: object
    variance: object
    mse: object
    absolute: dict[str, float]
    relative: dict[str, float]
    valid_relative_bins: tuple[int, ...]


@dataclass
class RegularizationScanResult:
    """Complete result of a paired-toy Bayesian iteration scan.

    ``selection_mode`` records whether the recommended iteration minimizes the
    absolute (spectrum-population-weighted) or fractional aggregate MSE.
    """

    metrics: tuple[ToyMetricResult, ...]
    absolute_summary: object
    relative_summary: object
    selected_iterations: int
    selection_mode: str = "absolute"
    absolute_selected_iterations: int | None = None
    relative_selected_iterations: int | None = None


@dataclass
class PtYieldMatchedTruth:
    """Truth scaled so its forward-folded pT-block yields match a target."""

    truth: object
    folded: object
    factors: tuple[float, ...]
    iterations: int
    maximum_relative_residual: float


def as_pt_intervals(bins: Sequence) -> tuple[tuple[float, float], ...]:
    """Accept either interval pairs or a contiguous sequence of bin edges."""

    values = tuple(bins)
    if not values:
        raise ValueError("pT bin configuration is empty")
    if all(isinstance(item, (tuple, list)) and len(item) == 2 for item in values):
        intervals = tuple((float(low), float(high)) for low, high in values)
    else:
        if len(values) < 2:
            raise ValueError("At least two pT bin edges are required")
        intervals = tuple(
            (float(low), float(high)) for low, high in zip(values[:-1], values[1:])
        )
    # Reuse the complete interval validation without imposing an eta layout.
    FlattenedBinning(intervals, 1)
    return intervals


def _same_axis_binning(left, right, *, tolerance: float = 1.0e-12) -> bool:
    if left.GetNbins() != right.GetNbins():
        return False
    for bin_index in range(1, left.GetNbins() + 2):
        if abs(left.GetBinLowEdge(bin_index) - right.GetBinLowEdge(bin_index)) > tolerance:
            return False
    return True


def validate_common_binning(histograms: Sequence, *, label: str = "histograms") -> None:
    if not histograms:
        raise ValueError(f"No {label} supplied")
    reference = histograms[0]
    if not reference.InheritsFrom("TH1"):
        raise TypeError(f"{reference.GetName()} is not a TH1")
    for histogram in histograms[1:]:
        if not histogram.InheritsFrom("TH1"):
            raise TypeError(f"{histogram.GetName()} is not a TH1")
        if not _same_axis_binning(reference.GetXaxis(), histogram.GetXaxis()):
            raise ValueError(
                f"Inconsistent x-axis binning in {label}: "
                f"{reference.GetName()} and {histogram.GetName()}"
            )


def _axis_bin_span(axis, value_range: tuple[float, float]) -> tuple[int, int]:
    """Return the contiguous bins whose centers lie in an inclusive range.

    Center-based inclusion preserves acceptance boundary bins after rebinning.
    For example, a rebinned eta bin spanning 1.8--2.0 has center 1.9 and is
    included for an eta range ending at 1.9.
    """

    low, high = value_range
    if low >= high:
        raise ValueError(f"Invalid axis range {value_range}")
    tolerance = 1.0e-12 * max(1.0, abs(low), abs(high))
    selected = [
        bin_index
        for bin_index in range(1, axis.GetNbins() + 1)
        if low - tolerance <= axis.GetBinCenter(bin_index) <= high + tolerance
    ]
    if not selected:
        raise ValueError(f"Axis range {value_range} contains no bins")
    return selected[0], selected[-1]


def restrict_histogram_range(histogram, value_range, *, name: str):
    """Copy bins whose centers lie in the requested inclusive range."""

    if not histogram.InheritsFrom("TH1") or histogram.InheritsFrom("TH2"):
        raise TypeError("restrict_histogram_range expects a one-dimensional histogram")
    first, last = _axis_bin_span(histogram.GetXaxis(), value_range)
    edges = array("d", [histogram.GetXaxis().GetBinLowEdge(first)])
    edges.extend(
        histogram.GetXaxis().GetBinUpEdge(bin_index)
        for bin_index in range(first, last + 1)
    )
    result = ROOT.TH1D(name, histogram.GetTitle(), last - first + 1, edges)
    result.SetDirectory(0)
    result.Sumw2()
    for output_bin, input_bin in enumerate(range(first, last + 1), 1):
        result.SetBinContent(output_bin, histogram.GetBinContent(input_bin))
        result.SetBinError(output_bin, histogram.GetBinError(input_bin))
    return result


def project_eta_by_pt(
    histogram,
    pt_bins: Sequence,
    *,
    name_prefix: str,
    eta_range: tuple[float, float] | None = None,
) -> list:
    """Project a TH2 eta distribution in configured half-open pT intervals."""

    if not histogram.InheritsFrom("TH2"):
        raise TypeError(f"{histogram.GetName()} is not a TH2")
    intervals = as_pt_intervals(pt_bins)
    projections = [
        project_semantic_th2(
            histogram, "eta", interval, name=f"{name_prefix}_pt{index}",
        )
        for index, interval in enumerate(intervals)
    ]
    if eta_range is None:
        return projections
    return [
        restrict_histogram_range(
            projection, eta_range, name=f"{name_prefix}_pt{index}_restricted",
        )
        for index, projection in enumerate(projections)
    ]


def flatten_pt_eta_projections(
    projections: Sequence,
    *,
    name: str,
    title: str = ";global #eta_{CM} bin;Entries",
    layout: FlattenedBinning | None = None,
    pt_bins: Sequence | None = None,
):
    """Copy one eta histogram per pT interval into a flattened TH1D."""

    projections = tuple(projections)
    validate_common_binning(projections, label="eta projections")
    if layout is None:
        if pt_bins is None:
            raise ValueError("pt_bins are required when layout is not supplied")
        layout = FlattenedBinning(as_pt_intervals(pt_bins), projections[0].GetNbinsX())
    if len(projections) != layout.n_pt_bins:
        raise ValueError("Projection count does not match the flattening layout")
    if projections[0].GetNbinsX() != layout.n_eta_bins:
        raise ValueError("Eta bin count does not match the flattening layout")

    flattened = ROOT.TH1D(
        name, title, layout.n_global_bins, 0.5, layout.n_global_bins + 0.5,
    )
    flattened.SetDirectory(0)
    flattened.Sumw2()
    for pt_index, projection in enumerate(projections):
        for eta_bin in range(1, layout.n_eta_bins + 1):
            global_bin = layout.global_bin(pt_index, eta_bin)
            flattened.SetBinContent(global_bin, projection.GetBinContent(eta_bin))
            flattened.SetBinError(global_bin, projection.GetBinError(eta_bin))
    return flattened, layout


def flatten_pt_eta_histogram(histogram, pt_bins: Sequence, *, name: str):
    projections = project_eta_by_pt(histogram, pt_bins, name_prefix=f"{name}Projection")
    flattened, layout = flatten_pt_eta_projections(
        projections, name=name, pt_bins=pt_bins,
    )
    return flattened, projections, layout


def extract_eta_block(flattened, template, layout: FlattenedBinning,
                      pt_index: int, *, name: str):
    """Extract one pT block from a flattened TH1 using an eta template."""

    if flattened.GetNbinsX() != layout.n_global_bins:
        raise ValueError("Flattened histogram size does not match layout")
    if template.GetNbinsX() != layout.n_eta_bins:
        raise ValueError("Eta template size does not match layout")
    result = template.Clone(name)
    result.Reset("ICES")
    result.SetDirectory(0)
    for eta_bin in range(1, layout.n_eta_bins + 1):
        global_bin = layout.global_bin(pt_index, eta_bin)
        result.SetBinContent(eta_bin, flattened.GetBinContent(global_bin))
        result.SetBinError(eta_bin, flattened.GetBinError(global_bin))
    return result


def unflatten_to_eta_projections(flattened, templates: Sequence,
                                 layout: FlattenedBinning, *, name_prefix: str) -> list:
    if len(templates) != layout.n_pt_bins:
        raise ValueError("Template count does not match layout")
    return [
        extract_eta_block(
            flattened, template, layout, pt_index,
            name=f"{name_prefix}_pt{pt_index}",
        )
        for pt_index, template in enumerate(templates)
    ]


def flatten_sparse_response(
    sparse_response,
    pt_bins: Sequence,
    *,
    name: str = "hResponseEtaCM",
    layout: FlattenedBinning | None = None,
    gen_pt_axis: int = 0,
    gen_eta_axis: int = 1,
    reco_pt_axis: int = 2,
    reco_eta_axis: int = 3,
    eta_range: tuple[float, float] | None = None,
):
    """Flatten a (gen pT, gen eta, reco pT, reco eta) THnSparse response."""

    if not sparse_response.InheritsFrom("THnSparse"):
        raise TypeError(f"{sparse_response.GetName()} is not a THnSparse")
    axis_indices = (gen_pt_axis, gen_eta_axis, reco_pt_axis, reco_eta_axis)
    if len(set(axis_indices)) != 4 or min(axis_indices) < 0:
        raise ValueError("Response axis indices must be four distinct non-negative values")
    if max(axis_indices) >= sparse_response.GetNdimensions():
        raise IndexError("Response axis index exceeds THnSparse dimensionality")
    intervals = as_pt_intervals(pt_bins)
    gen_eta = sparse_response.GetAxis(gen_eta_axis)
    reco_eta = sparse_response.GetAxis(reco_eta_axis)
    if gen_eta.GetNbins() != reco_eta.GetNbins():
        raise ValueError("Gen and reco response eta axes have different bin counts")
    if not _same_axis_binning(gen_eta, reco_eta):
        raise ValueError("Gen and reco response eta axes have different bin edges")
    if eta_range is None:
        gen_eta_first, gen_eta_last = 1, gen_eta.GetNbins()
        reco_eta_first, reco_eta_last = 1, reco_eta.GetNbins()
    else:
        gen_eta_first, gen_eta_last = _axis_bin_span(gen_eta, eta_range)
        reco_eta_first, reco_eta_last = _axis_bin_span(reco_eta, eta_range)
    gen_eta_bins = tuple(range(gen_eta_first, gen_eta_last + 1))
    reco_eta_bins = tuple(range(reco_eta_first, reco_eta_last + 1))
    if len(gen_eta_bins) != len(reco_eta_bins):
        raise ValueError("Selected gen and reco eta ranges have different bin counts")
    n_selected_eta = len(gen_eta_bins)
    layout = layout or FlattenedBinning(intervals, n_selected_eta)
    if layout.pt_intervals != intervals or layout.n_eta_bins != n_selected_eta:
        raise ValueError("Response dimensions do not match flattening layout")

    response = ROOT.TH2D(
        name, ";global reco #eta_{CM} bin;global gen #eta_{CM} bin",
        layout.n_global_bins, 0.5, layout.n_global_bins + 0.5,
        layout.n_global_bins, 0.5, layout.n_global_bins + 0.5,
    )
    response.SetDirectory(0)
    response.Sumw2()
    gen_axis = sparse_response.GetAxis(gen_pt_axis)
    reco_axis = sparse_response.GetAxis(reco_pt_axis)
    try:
        for gen_pt_index, (gen_low, gen_high) in enumerate(intervals):
            gen_axis.SetRange(
                gen_axis.FindBin(gen_low + 0.001), gen_axis.FindBin(gen_high - 0.001)
            )
            for reco_pt_index, (reco_low, reco_high) in enumerate(intervals):
                reco_axis.SetRange(
                    reco_axis.FindBin(reco_low + 0.001),
                    reco_axis.FindBin(reco_high - 0.001),
                )
                block = sparse_response.Projection(gen_eta_axis, reco_eta_axis)
                for local_gen_eta, source_gen_eta in enumerate(gen_eta_bins, 1):
                    global_gen = layout.global_bin(gen_pt_index, local_gen_eta)
                    for local_reco_eta, source_reco_eta in enumerate(reco_eta_bins, 1):
                        global_reco = layout.global_bin(reco_pt_index, local_reco_eta)
                        response.SetBinContent(
                            global_reco, global_gen,
                            block.GetBinContent(source_reco_eta, source_gen_eta),
                        )
                        response.SetBinError(
                            global_reco, global_gen,
                            block.GetBinError(source_reco_eta, source_gen_eta),
                        )
    finally:
        gen_axis.SetRange(0, -1)
        reco_axis.SetRange(0, -1)
    return response, layout


def project_response_eta_blocks(
    sparse_response,
    pt_bins: Sequence,
    *,
    name_prefix: str = "hResponseEtaCM",
    gen_pt_axis: int = 0,
    gen_eta_axis: int = 1,
    reco_pt_axis: int = 2,
    reco_eta_axis: int = 3,
) -> list:
    """Project diagonal eta-response blocks for pT-interval diagnostics."""

    if not sparse_response.InheritsFrom("THnSparse"):
        raise TypeError(f"{sparse_response.GetName()} is not a THnSparse")
    intervals = as_pt_intervals(pt_bins)
    gen_axis = sparse_response.GetAxis(gen_pt_axis)
    reco_axis = sparse_response.GetAxis(reco_pt_axis)
    blocks = []
    try:
        for index, (low, high) in enumerate(intervals):
            gen_axis.SetRange(gen_axis.FindBin(low + 0.001), gen_axis.FindBin(high - 0.001))
            reco_axis.SetRange(reco_axis.FindBin(low + 0.001), reco_axis.FindBin(high - 0.001))
            block = sparse_response.Projection(gen_eta_axis, reco_eta_axis)
            block.SetName(f"{name_prefix}_pt{index}")
            block.SetDirectory(0)
            blocks.append(block)
    finally:
        gen_axis.SetRange(0, -1)
        reco_axis.SetRange(0, -1)
    return blocks


def _assert_no_negative_bins(histogram, label: str, tolerance: float) -> None:
    negative = [
        (index, histogram.GetBinContent(index))
        for index in range(1, histogram.GetNbinsX() + 1)
        if histogram.GetBinContent(index) < -tolerance
    ]
    if negative:
        raise ValueError(f"{label} has negative bins {negative[:10]}")


def calculate_response_diagnostics(
    response,
    training_truth,
    training_measured,
    *,
    explicit_miss=None,
    explicit_fake=None,
    tolerance: float = 1.0e-9,
) -> ResponseDiagnostics:
    """Calculate response marginals and effective/boundary miss/fake spectra.

    With matrix convention ``M[measured, truth]``, ProjectionX is matched reco
    and ProjectionY is matched truth.  Effective components include both the
    explicitly stored selection/matching failures and migrations outside the
    configured flattened acceptance:

    ``effective miss = inclusive truth - matched truth``
    ``effective fake = inclusive reco - matched reco``.
    """

    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative")
    if not response.InheritsFrom("TH2"):
        raise TypeError("Response must be a TH2")
    validate_common_binning(
        [training_truth, response.ProjectionY("_diagnostic_truth_binning")],
        label="training truth and response truth",
    )
    validate_common_binning(
        [training_measured, response.ProjectionX("_diagnostic_reco_binning")],
        label="training measured and response measured",
    )
    matched_measured = response.ProjectionX(f"{response.GetName()}MatchedMeasured")
    matched_truth = response.ProjectionY(f"{response.GetName()}MatchedTruth")
    matched_measured.SetDirectory(0)
    matched_truth.SetDirectory(0)
    effective_miss = training_truth.Clone(f"{response.GetName()}EffectiveMiss")
    effective_miss.Add(matched_truth, -1.0)
    effective_fake = training_measured.Clone(f"{response.GetName()}EffectiveFake")
    effective_fake.Add(matched_measured, -1.0)
    for histogram in (effective_miss, effective_fake):
        histogram.SetDirectory(0)
    _assert_no_negative_bins(effective_miss, "Effective miss", tolerance)
    _assert_no_negative_bins(effective_fake, "Effective fake", tolerance)

    boundary_miss = None
    if explicit_miss is not None:
        validate_common_binning([effective_miss, explicit_miss], label="miss diagnostics")
        boundary_miss = effective_miss.Clone(f"{response.GetName()}BoundaryMiss")
        boundary_miss.Add(explicit_miss, -1.0)
        boundary_miss.SetDirectory(0)
        _assert_no_negative_bins(boundary_miss, "Boundary-migration miss", tolerance)
    boundary_fake = None
    if explicit_fake is not None:
        validate_common_binning([effective_fake, explicit_fake], label="fake diagnostics")
        boundary_fake = effective_fake.Clone(f"{response.GetName()}BoundaryFake")
        boundary_fake.Add(explicit_fake, -1.0)
        boundary_fake.SetDirectory(0)
        _assert_no_negative_bins(boundary_fake, "Boundary-migration fake", tolerance)
    return ResponseDiagnostics(
        matched_truth, matched_measured, effective_miss, effective_fake,
        boundary_miss, boundary_fake,
    )


def validate_response_accounting(
    response,
    training_truth,
    training_measured,
    diagnostics: ResponseDiagnostics,
    *,
    tolerance: float = 1.0e-9,
) -> ResponseAccounting:
    """Validate response marginals and summarize efficiency, fakes, and sparsity."""

    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative")
    validate_common_binning(
        [training_truth, diagnostics.matched_truth, diagnostics.effective_miss],
        label="truth response accounting",
    )
    validate_common_binning(
        [training_measured, diagnostics.matched_measured, diagnostics.effective_fake],
        label="measured response accounting",
    )
    truth_absolute = []
    truth_relative = []
    measured_absolute = []
    measured_relative = []
    zero_efficiency = []
    empty_measured = []
    truth_zero_tolerance = tolerance * max(abs(training_truth.Integral()), 1.0e-300)
    measured_zero_tolerance = tolerance * max(abs(training_measured.Integral()), 1.0e-300)
    for bin_index in range(1, training_truth.GetNbinsX() + 1):
        expected = training_truth.GetBinContent(bin_index)
        accounted = (
            diagnostics.matched_truth.GetBinContent(bin_index)
            + diagnostics.effective_miss.GetBinContent(bin_index)
        )
        residual = abs(expected - accounted)
        truth_absolute.append(residual)
        truth_relative.append(residual / max(abs(expected), 1.0e-300))
        if (expected > truth_zero_tolerance
                and diagnostics.matched_truth.GetBinContent(bin_index) <= truth_zero_tolerance):
            zero_efficiency.append(bin_index)
    for bin_index in range(1, training_measured.GetNbinsX() + 1):
        expected = training_measured.GetBinContent(bin_index)
        accounted = (
            diagnostics.matched_measured.GetBinContent(bin_index)
            + diagnostics.effective_fake.GetBinContent(bin_index)
        )
        residual = abs(expected - accounted)
        measured_absolute.append(residual)
        measured_relative.append(residual / max(abs(expected), 1.0e-300))
        if abs(expected) <= measured_zero_tolerance:
            empty_measured.append(bin_index)
    max_truth = max(truth_absolute, default=0.0)
    max_measured = max(measured_absolute, default=0.0)
    max_truth_relative = max(truth_relative, default=0.0)
    max_measured_relative = max(measured_relative, default=0.0)
    absolute_scale = max(
        1.0,
        abs(training_truth.Integral()),
        abs(training_measured.Integral()),
    )
    if max(max_truth, max_measured) > tolerance * absolute_scale:
        raise ValueError(
            "Response accounting failed: "
            f"truth residual={max_truth:g}, measured residual={max_measured:g}"
        )
    populated = sum(
        response.GetBinContent(reco_bin, truth_bin) != 0.0
        for reco_bin in range(1, response.GetNbinsX() + 1)
        for truth_bin in range(1, response.GetNbinsY() + 1)
    )
    matrix_bins = response.GetNbinsX() * response.GetNbinsY()
    truth_integral = training_truth.Integral()
    measured_integral = training_measured.Integral()
    return ResponseAccounting(
        max_truth, max_measured, max_truth_relative, max_measured_relative,
        tuple(zero_efficiency), tuple(empty_measured),
        1.0 - populated / matrix_bins,
        diagnostics.matched_truth.Integral() / truth_integral if truth_integral else 0.0,
        diagnostics.effective_fake.Integral() / measured_integral if measured_integral else 0.0,
    )


def build_roounfold_response(
    roounfold_module,
    training_truth,
    training_measured,
    response_matrix,
    *,
    diagnostics: ResponseDiagnostics | None = None,
    scale: float = 1.0e12,
    require_fakes: bool = True,
    name: str = "responseEtaCM",
    title: str = "Flattened dijet pTave-etaCM response",
) -> RooUnfoldResponseBundle:
    """Build a scaled RooUnfoldResponse from mutually consistent marginals.

    For the conventional path, pass inclusive truth/reco; RooUnfold then infers
    misses and fakes from their differences with the matrix projections.  For
    the factorized path, pass matched truth/reco and ``require_fakes=False``;
    the response then contains migration probabilities only.

    All inputs are cloned before scaling, so caller-owned analysis histograms
    remain unchanged.  ``roounfold_module`` is retained in the public interface
    to make the required loaded dependency explicit; PyROOT constructs the C++
    object through ``ROOT``.
    """

    if scale <= 0.0:
        raise ValueError("Response scale must be positive")
    validate_common_binning(
        [training_truth, response_matrix.ProjectionY("_response_truth_binning")],
        label="response truth inputs",
    )
    validate_common_binning(
        [training_measured, response_matrix.ProjectionX("_response_reco_binning")],
        label="response measured inputs",
    )
    scaled_truth = training_truth.Clone(f"{name}ScaledTruth")
    scaled_measured = training_measured.Clone(f"{name}ScaledMeasured")
    scaled_matrix = response_matrix.Clone(f"{name}ScaledMatrix")
    for histogram in (scaled_truth, scaled_measured, scaled_matrix):
        histogram.SetDirectory(0)
        histogram.Scale(scale)
    response = ROOT.RooUnfoldResponse(
        scaled_measured, scaled_truth, scaled_matrix, name, title,
    )
    if require_fakes and not response.HasFakes():
        raise RuntimeError("The RooUnfold response has no fake component")
    if diagnostics is not None:
        response_fakes = response.Hfakes()
        for bin_index in range(1, training_measured.GetNbinsX() + 1):
            expected = scale * diagnostics.effective_fake.GetBinContent(bin_index)
            actual = response_fakes.GetBinContent(bin_index)
            tolerance = 1.0e-9 * max(1.0, abs(expected), abs(actual))
            if abs(actual - expected) > tolerance:
                raise ValueError(
                    f"RooUnfold fake mismatch in bin {bin_index}: "
                    f"expected {expected}, got {actual}"
                )
    return RooUnfoldResponseBundle(
        response, scaled_truth, scaled_measured, scaled_matrix, float(scale),
    )


def unfold_bayes(
    roounfold_module,
    response_bundle: RooUnfoldResponseBundle,
    measured,
    *,
    iterations: int,
    handle_fakes: bool = True,
    suppress_native_output: bool = False,
    name: str = "hUnfoldedEtaCM",
    title: str = ";global #eta_{CM} bin;Entries",
) -> UnfoldingResult:
    """Run Bayesian unfolding and restore the original weighted-event scale.

    The response inputs and measured spectrum are multiplied by the common
    response scale before unfolding. RooUnfold therefore returns a spectrum
    scaled by S and a covariance scaled by S^2. This function divides those
    products by S and S^2, respectively. The Bayesian probabilities and the
    regularizing prior shape are unchanged by the common scale.
    """

    if iterations <= 0:
        raise ValueError("Bayesian iteration count must be positive")
    validate_common_binning(
        [measured, response_bundle.scaled_measured], label="unfolding measured input"
    )
    scaled_measured = measured.Clone(f"{name}ScaledMeasured")
    scaled_measured.SetDirectory(0)
    scaled_measured.Scale(response_bundle.scale)
    algorithm = ROOT.RooUnfoldBayes(
        response_bundle.response, scaled_measured, iterations, False, handle_fakes,
    )
    algorithm.SetVerbose(-1)
    with _silence_native_stdout(suppress_native_output):
        unfolded = algorithm.Hunfold().Clone(name)
        covariance = algorithm.Eunfold(ROOT.RooUnfolding.kCovariance)
    unfolded.SetTitle(title)
    unfolded.SetDirectory(0)
    unfolded.Scale(1.0 / response_bundle.scale)
    covariance *= 1.0 / (response_bundle.scale * response_bundle.scale)
    return UnfoldingResult(unfolded, covariance, algorithm, scaled_measured)


def prepare_factorized_corrections(
    measured,
    training_measured,
    matched_measured,
    training_truth,
    matched_truth,
    *,
    name_prefix: str = "hFactorized",
) -> FactorizedCorrectionInputs:
    """Remove fakes with reco-bin purity and calculate truth-bin efficiency.

    In measured bin i, ``purity_i = matched_reco_i / inclusive_reco_i`` and
    ``measured_signal_i = measured_i * purity_i``.  In truth bin j,
    ``efficiency_j = matched_truth_j / inclusive_truth_j``.  A zero training
    reco bin produces zero purity; a zero-efficiency truth bin is rejected later
    because no unfolding can recover a truth bin that is never reconstructed.

    The correction factors are treated as fixed response inputs.  Consequently,
    the corrected measured errors contain the measured-data uncertainty only;
    finite-MC uncertainty in purity and efficiency belongs in response
    statistical/systematic variations rather than being counted independently
    in every multiplication here.
    """

    validate_common_binning(
        [measured, training_measured, matched_measured],
        label="factorized fake correction inputs",
    )
    validate_common_binning(
        [training_truth, matched_truth],
        label="factorized efficiency correction inputs",
    )
    purity = matched_measured.Clone(f"{name_prefix}Purity")
    efficiency = matched_truth.Clone(f"{name_prefix}Efficiency")
    for histogram in (purity, efficiency):
        histogram.SetDirectory(0)
    for bin_index in range(1, measured.GetNbinsX() + 1):
        inclusive = training_measured.GetBinContent(bin_index)
        matched = matched_measured.GetBinContent(bin_index)
        if inclusive < 0.0 or matched < 0.0 or matched > inclusive + 1.0e-9 * max(1.0, inclusive):
            raise ValueError(
                f"Invalid reco purity inputs in bin {bin_index}: "
                f"matched={matched:g}, inclusive={inclusive:g}"
            )
        factor = matched / inclusive if inclusive > 0.0 else 0.0
        # Weighted projections can differ at the last few floating-point bits.
        # Preserve the physical [0, 1] range after the content-level check above.
        factor = min(max(factor, 0.0), 1.0)
        purity.SetBinContent(bin_index, factor)
        purity.SetBinError(bin_index, 0.0)
    for bin_index in range(1, training_truth.GetNbinsX() + 1):
        inclusive = training_truth.GetBinContent(bin_index)
        matched = matched_truth.GetBinContent(bin_index)
        if inclusive < 0.0 or matched < 0.0 or matched > inclusive + 1.0e-9 * max(1.0, inclusive):
            raise ValueError(
                f"Invalid truth efficiency inputs in bin {bin_index}: "
                f"matched={matched:g}, inclusive={inclusive:g}"
            )
        factor = matched / inclusive if inclusive > 0.0 else 0.0
        factor = min(max(factor, 0.0), 1.0)
        efficiency.SetBinContent(bin_index, factor)
        efficiency.SetBinError(bin_index, 0.0)
    measured_signal = apply_purity_correction(
        measured, purity, name=f"{name_prefix}MeasuredSignal",
    )
    return FactorizedCorrectionInputs(measured_signal, purity, efficiency)


def apply_purity_correction(measured, purity, *, name: str):
    """Multiply inclusive reco by a fixed reco-bin purity.

    Both contents and input errors are multiplied by purity.  Purity uncertainty
    is intentionally not added here because purity is part of the response
    model and should be varied coherently with the migration and efficiency.
    """

    validate_common_binning([measured, purity], label="purity correction inputs")
    corrected = measured.Clone(name)
    corrected.SetDirectory(0)
    for bin_index in range(1, measured.GetNbinsX() + 1):
        factor = purity.GetBinContent(bin_index)
        if not math.isfinite(factor) or factor < 0.0 or factor > 1.0 + 1.0e-9:
            raise ValueError(
                f"Invalid purity in measured bin {bin_index}: {factor:g}"
            )
        corrected.SetBinContent(
            bin_index, measured.GetBinContent(bin_index) * factor,
        )
        corrected.SetBinError(
            bin_index, measured.GetBinError(bin_index) * factor,
        )
    return corrected


def apply_efficiency_correction(
    unfolded_matched,
    covariance,
    efficiency,
    *,
    name: str = "hUnfoldedEfficiencyCorrected",
):
    """Convert matched truth to inclusive truth and transform covariance.

    If ``D_jj = 1 / efficiency_j``, the corrected result is ``D t`` and its
    covariance is ``D C D``.  This propagates the unfolding covariance while
    treating the efficiency itself as a fixed response input.
    """

    validate_common_binning(
        [unfolded_matched, efficiency], label="efficiency correction inputs"
    )
    corrected = unfolded_matched.Clone(name)
    corrected.SetDirectory(0)
    corrected_covariance = covariance.Clone()
    for row in range(unfolded_matched.GetNbinsX()):
        row_efficiency = efficiency.GetBinContent(row + 1)
        if row_efficiency <= 0.0:
            raise ValueError(
                f"Cannot efficiency-correct truth bin {row + 1}: "
                f"efficiency={row_efficiency:g}"
            )
        corrected.SetBinContent(
            row + 1, unfolded_matched.GetBinContent(row + 1) / row_efficiency,
        )
        corrected.SetBinError(
            row + 1, unfolded_matched.GetBinError(row + 1) / row_efficiency,
        )
        for column in range(unfolded_matched.GetNbinsX()):
            column_efficiency = efficiency.GetBinContent(column + 1)
            if column_efficiency <= 0.0:
                raise ValueError(
                    f"Cannot efficiency-correct truth bin {column + 1}: "
                    f"efficiency={column_efficiency:g}"
                )
            corrected_covariance[row][column] = (
                covariance[row][column] / (row_efficiency * column_efficiency)
            )
    return corrected, corrected_covariance


def calculate_closure_metrics(
    estimate,
    reference,
    *,
    populated_fraction: float = 1.0e-3,
) -> ClosureMetrics:
    """Summarize closure without letting nearly empty bins dominate.

    ``relative_l1`` is ``sum(|estimate-reference|) / sum(|reference|)`` over all
    bins and is the most stable shape-and-yield summary.  The mean absolute
    relative difference is reported only for reference bins above
    ``populated_fraction * reference.GetMaximum()``.  This threshold affects
    only that diagnostic; no analysis histogram is masked or modified.
    """

    if populated_fraction < 0.0:
        raise ValueError("populated_fraction must be non-negative")
    validate_common_binning([estimate, reference], label="closure metric inputs")
    reference_integral = reference.Integral()
    if reference_integral == 0.0:
        raise ValueError("Cannot calculate closure against zero-integral reference")
    absolute_reference = sum(
        abs(reference.GetBinContent(bin_index))
        for bin_index in range(1, reference.GetNbinsX() + 1)
    )
    if absolute_reference == 0.0:
        raise ValueError("Cannot calculate closure against empty reference")
    relative_l1 = sum(
        abs(estimate.GetBinContent(bin_index) - reference.GetBinContent(bin_index))
        for bin_index in range(1, reference.GetNbinsX() + 1)
    ) / absolute_reference
    threshold = populated_fraction * max(reference.GetMaximum(), 0.0)
    populated = [
        bin_index
        for bin_index in range(1, reference.GetNbinsX() + 1)
        if abs(reference.GetBinContent(bin_index)) > threshold
    ]
    mean_absolute_relative = sum(
        abs(estimate.GetBinContent(bin_index) / reference.GetBinContent(bin_index) - 1.0)
        for bin_index in populated
    ) / len(populated) if populated else math.nan
    return ClosureMetrics(
        estimate.Integral() / reference_integral,
        relative_l1,
        mean_absolute_relative,
        len(populated),
    )


def histogram_fingerprint(histogram) -> str:
    """Return a stable digest of histogram dimensions, contents, and errors.

    Cache metadata uses this digest to reject intermediate results if an input
    ROOT histogram changed without its path or object name changing.
    """

    if not histogram.InheritsFrom("TH1"):
        raise TypeError(f"{histogram.GetName()} is not a histogram")
    dimensions = histogram.GetDimension()
    counts = [histogram.GetNbinsX()]
    if dimensions >= 2:
        counts.append(histogram.GetNbinsY())
    if dimensions >= 3:
        counts.append(histogram.GetNbinsZ())
    digest = hashlib.sha256()
    digest.update(struct.pack("!i", dimensions))
    digest.update(struct.pack(f"!{len(counts)}i", *counts))
    ranges = [range(0, count + 2) for count in counts]
    if dimensions == 1:
        cells = ((x,) for x in ranges[0])
    elif dimensions == 2:
        cells = ((x, y) for x in ranges[0] for y in ranges[1])
    else:
        cells = (
            (x, y, z)
            for x in ranges[0] for y in ranges[1] for z in ranges[2]
        )
    for indices in cells:
        global_bin = histogram.GetBin(*indices)
        digest.update(struct.pack(
            "!dd", histogram.GetBinContent(global_bin),
            histogram.GetBinError(global_bin),
        ))
    return digest.hexdigest()


def _histogram_distance(left, right) -> float:
    validate_common_binning([left, right], label="histogram-distance inputs")
    scale = sum(abs(right.GetBinContent(i)) for i in range(1, right.GetNbinsX() + 1))
    difference = sum(
        abs(left.GetBinContent(i) - right.GetBinContent(i))
        for i in range(1, left.GetNbinsX() + 1)
    )
    return difference / max(scale, 1.0e-300)


def forward_fold_truth(
    response_bundle: RooUnfoldResponseBundle,
    truth,
    *,
    diagnostics: ResponseDiagnostics,
    name: str = "hForwardFoldedTruth",
    add_fakes: bool = True,
    fake_normalization: str = "training_yield",
    detection_tolerance: float = 1.0e-8,
) -> ForwardFoldResult:
    """Apply the inclusive detector response and restore fakes exactly once.

    ``RooUnfoldResponse.ApplyToTruth`` normally returns only the response-matrix
    contribution.  When the bundle was trained with inclusive truth, its matrix
    columns already contain the efficiency: folding therefore maps inclusive
    truth to matched reco.  The training truth is folded first to verify that behavior
    against the matched and inclusive measured spectra.  This protects the
    pseudo-data construction if RooUnfold changes its implementation.
    ``fake_normalization='matched_fraction'`` then preserves the controlled
    fake-to-matched-reco ratio independently in every measured bin, so the fake
    prediction follows the normalization of an independently folded sample.  In
    purity notation this final addition is equivalent to dividing matched reco
    by purity.  Passing a matched-only response here would omit inefficiency and
    is therefore intentionally avoided by the notebooks.
    """

    if fake_normalization not in {"training_yield", "matched_fraction"}:
        raise ValueError(
            "fake_normalization must be 'training_yield' or 'matched_fraction'"
        )
    validate_common_binning([truth, response_bundle.scaled_truth], label="forward-fold truth")
    scaled_truth = truth.Clone(f"{name}ScaledTruth")
    scaled_truth.SetDirectory(0)
    scaled_truth.Scale(response_bundle.scale)
    folded_scaled = response_bundle.response.ApplyToTruth(scaled_truth)
    folded_signal = folded_scaled.Clone(f"{name}Signal")
    folded_signal.SetDirectory(0)
    folded_signal.Scale(1.0 / response_bundle.scale)

    training_folded_scaled = response_bundle.response.ApplyToTruth(
        response_bundle.scaled_truth
    )
    training_folded = training_folded_scaled.Clone(f"{name}TrainingCheck")
    training_folded.SetDirectory(0)
    training_folded.Scale(1.0 / response_bundle.scale)
    matched_distance = _histogram_distance(training_folded, diagnostics.matched_measured)
    inclusive_training = diagnostics.matched_measured.Clone(f"{name}InclusiveCheck")
    inclusive_training.SetDirectory(0)
    inclusive_training.Add(diagnostics.effective_fake)
    inclusive_distance = _histogram_distance(training_folded, inclusive_training)
    includes_fakes = inclusive_distance + detection_tolerance < matched_distance
    if includes_fakes and add_fakes and fake_normalization != "training_yield":
        raise RuntimeError(
            "ApplyToTruth already includes fakes, so their normalization cannot be replaced"
        )

    fake = diagnostics.effective_fake.Clone(f"{name}Fake")
    fake.SetDirectory(0)
    if fake_normalization == "matched_fraction":
        for bin_index in range(1, fake.GetNbinsX() + 1):
            matched = diagnostics.matched_measured.GetBinContent(bin_index)
            matched_error = diagnostics.matched_measured.GetBinError(bin_index)
            training_fake = diagnostics.effective_fake.GetBinContent(bin_index)
            training_fake_error = diagnostics.effective_fake.GetBinError(bin_index)
            folded = folded_signal.GetBinContent(bin_index)
            folded_error = folded_signal.GetBinError(bin_index)
            if matched <= 0.0:
                if training_fake > 0.0:
                    raise ValueError(
                        f"Cannot preserve fake fraction in measured bin {bin_index}: "
                        f"matched={matched:g}, fake={training_fake:g}"
                    )
                fake.SetBinContent(bin_index, 0.0)
                fake.SetBinError(bin_index, 0.0)
                continue
            fake_to_matched = training_fake / matched
            ratio_error = math.hypot(
                training_fake_error / matched,
                training_fake * matched_error / (matched * matched),
            )
            fake.SetBinContent(bin_index, folded * fake_to_matched)
            fake.SetBinError(
                bin_index,
                math.hypot(folded_error * fake_to_matched, folded * ratio_error),
            )
    result = folded_signal.Clone(name)
    result.SetDirectory(0)
    if add_fakes and not includes_fakes:
        result.Add(fake)
    elif not add_fakes or includes_fakes:
        fake.Reset("ICES")
    return ForwardFoldResult(result, folded_signal, fake, includes_fakes)


def match_pt_block_yields(
    source,
    target,
    layout: FlattenedBinning,
    *,
    name: str = "hForwardFoldedPriorDataYield",
) -> tuple[object, tuple[float, ...]]:
    """Scale each flattened pT block to the matching target integral."""

    validate_common_binning([source, target], label="pT-block yield matching")
    if source.GetNbinsX() != layout.n_global_bins:
        raise ValueError("Yield-matching histogram size does not match layout")
    result = source.Clone(name)
    result.SetDirectory(0)
    factors = []
    for pt_index in range(layout.n_pt_bins):
        first = layout.global_bin(pt_index, 1)
        last = layout.global_bin(pt_index, layout.n_eta_bins)
        source_yield = source.Integral(first, last)
        target_yield = target.Integral(first, last)
        if source_yield <= 0.0 or target_yield <= 0.0:
            raise ValueError(
                f"Cannot match pT block {pt_index}: source={source_yield}, "
                f"target={target_yield}"
            )
        factor = target_yield / source_yield
        factors.append(factor)
        for bin_index in range(first, last + 1):
            result.SetBinContent(bin_index, factor * source.GetBinContent(bin_index))
            result.SetBinError(bin_index, target.GetBinError(bin_index))
    return result, tuple(factors)


def scale_pt_blocks(
    histogram,
    layout: FlattenedBinning,
    factors: Sequence[float],
    *,
    name: str,
):
    """Apply one exposure factor to each pT block of a flattened histogram."""

    if histogram.GetNbinsX() != layout.n_global_bins:
        raise ValueError("Histogram size does not match layout")
    if len(factors) != layout.n_pt_bins:
        raise ValueError("Scale-factor count does not match pT blocks")
    result = histogram.Clone(name)
    result.SetDirectory(0)
    for pt_index, factor in enumerate(factors):
        if not math.isfinite(factor) or factor <= 0.0:
            raise ValueError(f"Invalid scale factor {factor} for pT block {pt_index}")
        for eta_bin in range(1, layout.n_eta_bins + 1):
            bin_index = layout.global_bin(pt_index, eta_bin)
            result.SetBinContent(bin_index, factor * histogram.GetBinContent(bin_index))
            result.SetBinError(bin_index, factor * histogram.GetBinError(bin_index))
    return result


def match_truth_to_reco_pt_yields(
    response_bundle: RooUnfoldResponseBundle,
    truth,
    target_reco,
    diagnostics: ResponseDiagnostics,
    layout: FlattenedBinning,
    *,
    name_prefix: str = "hYieldMatchedPrior",
    tolerance: float = 1.0e-8,
    max_iterations: int = 200,
    relaxation: float = 0.5,
) -> PtYieldMatchedTruth:
    """Scale truth pT blocks until its folded reco block yields match data.

    A reco-block scale factor cannot be copied directly to the corresponding
    truth block when the response has off-diagonal pT migration.  This routine
    instead iterates in truth space: scale truth, forward-fold it with the full
    efficiency and fake model, compare reco block yields, and update the truth
    factors.  ``relaxation < 1`` damps migration-induced oscillations.
    """

    if tolerance <= 0.0 or max_iterations <= 0:
        raise ValueError("Positive tolerance and max_iterations are required")
    if not 0.0 < relaxation <= 1.0:
        raise ValueError("relaxation must be in (0, 1]")
    validate_common_binning([target_reco, diagnostics.matched_measured],
                            label="reco yield-matching inputs")
    if truth.GetNbinsX() != layout.n_global_bins:
        raise ValueError("Truth histogram size does not match layout")
    target_yields = []
    for pt_index in range(layout.n_pt_bins):
        first = layout.global_bin(pt_index, 1)
        last = layout.global_bin(pt_index, layout.n_eta_bins)
        target_yield = target_reco.Integral(first, last)
        if target_yield <= 0.0:
            raise ValueError(
                f"Reco target pT block {pt_index} has non-positive yield {target_yield:g}"
            )
        target_yields.append(target_yield)

    factors = [1.0] * layout.n_pt_bins
    for iteration in range(1, max_iterations + 1):
        scaled_truth = scale_pt_blocks(
            truth, layout, factors, name=f"{name_prefix}TruthStep{iteration}",
        )
        fold = forward_fold_truth(
            response_bundle, scaled_truth, diagnostics=diagnostics,
            add_fakes=True, fake_normalization="matched_fraction",
            name=f"{name_prefix}FoldedStep{iteration}",
        )
        ratios = []
        residuals = []
        for pt_index, target_yield in enumerate(target_yields):
            first = layout.global_bin(pt_index, 1)
            last = layout.global_bin(pt_index, layout.n_eta_bins)
            folded_yield = fold.histogram.Integral(first, last)
            if folded_yield <= 0.0:
                raise ValueError(
                    f"Folded pT block {pt_index} has non-positive yield {folded_yield:g}"
                )
            ratio = target_yield / folded_yield
            ratios.append(ratio)
            residuals.append(abs(folded_yield / target_yield - 1.0))
        maximum_residual = max(residuals)
        if maximum_residual <= tolerance:
            scaled_truth.SetName(f"{name_prefix}Truth")
            fold.histogram.SetName(f"{name_prefix}FoldedReco")
            return PtYieldMatchedTruth(
                scaled_truth, fold.histogram, tuple(factors), iteration,
                maximum_residual,
            )
        factors = [
            factor * ratio ** relaxation
            for factor, ratio in zip(factors, ratios)
        ]
    raise RuntimeError(
        f"Truth-to-reco pT yield matching did not converge after {max_iterations} "
        f"iterations; maximum relative residual={maximum_residual:g}"
    )


def copy_histogram_errors(source, error_reference, *, name: str):
    """Clone source contents while taking absolute bin errors from reference."""

    validate_common_binning([source, error_reference], label="error-copy inputs")
    result = source.Clone(name)
    result.SetDirectory(0)
    for bin_index in range(1, source.GetNbinsX() + 1):
        result.SetBinError(bin_index, error_reference.GetBinError(bin_index))
    return result


def make_gaussian_toys(
    expectation,
    *,
    n_toys: int,
    random_seed: int,
    layout: FlattenedBinning | None = None,
    stop_on_negative: bool = True,
    name_prefix: str = "hRecoToy",
) -> list:
    """Generate paired Gaussian toys while preserving input bin errors."""

    if n_toys < 2:
        raise ValueError("At least two toys are required for an unbiased variance")
    if layout is not None and expectation.GetNbinsX() != layout.n_global_bins:
        raise ValueError("Toy expectation size does not match layout")
    rng = random.Random(random_seed)
    toys = []
    for toy_index in range(n_toys):
        toy = expectation.Clone(f"{name_prefix}_{toy_index}")
        toy.SetDirectory(0)
        for bin_index in range(1, expectation.GetNbinsX() + 1):
            mean = expectation.GetBinContent(bin_index)
            sigma = expectation.GetBinError(bin_index)
            if not math.isfinite(mean) or not math.isfinite(sigma) or sigma < 0.0:
                raise ValueError(f"Invalid Gaussian parameters in bin {bin_index}")
            value = rng.gauss(mean, sigma) if sigma > 0.0 else mean
            if value < 0.0 and stop_on_negative:
                pt_text = "unknown"
                if layout is not None:
                    pt_index, eta_bin = layout.indices(bin_index)
                    low, high = layout.pt_intervals[pt_index]
                    pt_text = f"[{low:g}, {high:g}), eta bin {eta_bin}"
                probability = 0.5 * math.erfc(mean / (sigma * math.sqrt(2.0)))
                raise ValueError(
                    f"Negative Gaussian toy content in toy {toy_index}, global bin "
                    f"{bin_index} ({pt_text}): mean={mean:g}, sigma={sigma:g}, "
                    f"P(x<0)={probability:.3g}. Use coarser binning or revise the "
                    "Gaussian approximation."
                )
            toy.SetBinContent(bin_index, value)
            toy.SetBinError(bin_index, sigma)
        toys.append(toy)
    return toys


def calculate_toy_metrics(
    unfolded_toys: Sequence,
    truth,
    *,
    iterations: int,
    relative_truth_threshold: float = 1.0e-12,
    name_prefix: str = "hRegularization",
) -> ToyMetricResult:
    """Calculate per-bin and aggregate bias squared, variance, and MSE."""

    toys = tuple(unfolded_toys)
    if len(toys) < 2:
        raise ValueError("At least two unfolded toys are required")
    validate_common_binning([truth, *toys], label="toy metric inputs")
    n_bins = truth.GetNbinsX()
    sums = [0.0] * (n_bins + 1)
    sums_squared = [0.0] * (n_bins + 1)
    for toy in toys:
        for bin_index in range(1, n_bins + 1):
            value = toy.GetBinContent(bin_index)
            sums[bin_index] += value
            sums_squared[bin_index] += value * value
    return _toy_metrics_from_moments(
        truth, sums, sums_squared, len(toys), iterations=iterations,
        relative_truth_threshold=relative_truth_threshold, name_prefix=name_prefix,
    )


def _toy_metrics_from_moments(
    truth,
    sums: Sequence[float],
    sums_squared: Sequence[float],
    n_toys: int,
    *,
    iterations: int,
    relative_truth_threshold: float = 1.0e-12,
    name_prefix: str = "hRegularization",
) -> ToyMetricResult:
    """Build per-bin metrics and two deliberately different bin averages.

    ``absolute`` averages the metric values themselves, so populated bins have
    the larger influence naturally implied by their squared contents.
    ``relative`` divides each bin by its own truth squared before averaging and
    is retained as a diagnostic of typical fractional performance.
    """

    if n_toys < 2:
        raise ValueError("At least two unfolded toys are required")
    n_bins = truth.GetNbinsX()
    if len(sums) != n_bins + 1 or len(sums_squared) != n_bins + 1:
        raise ValueError("Toy moment arrays do not match truth binning")
    outputs = {}
    for key in ("Mean", "Bias", "BiasSquared", "Variance", "MSE"):
        histogram = truth.Clone(f"{name_prefix}{key}_iter{iterations}")
        histogram.Reset("ICES")
        histogram.SetDirectory(0)
        outputs[key] = histogram
    absolute_sums = {key: 0.0 for key in ("bias_squared", "variance", "mse")}
    relative_sums = {key: 0.0 for key in absolute_sums}
    valid_relative = []
    for bin_index in range(1, n_bins + 1):
        mean = sums[bin_index] / n_toys
        bias = mean - truth.GetBinContent(bin_index)
        bias_squared = bias * bias
        centered_sum = sums_squared[bin_index] - n_toys * mean * mean
        variance = max(0.0, centered_sum / (n_toys - 1))
        mse = bias_squared + variance
        for key, value in (
            ("Mean", mean), ("Bias", bias), ("BiasSquared", bias_squared),
            ("Variance", variance), ("MSE", mse),
        ):
            outputs[key].SetBinContent(bin_index, value)
        absolute_sums["bias_squared"] += bias_squared
        absolute_sums["variance"] += variance
        absolute_sums["mse"] += mse
        truth_value = truth.GetBinContent(bin_index)
        if abs(truth_value) > relative_truth_threshold:
            valid_relative.append(bin_index)
            denominator = truth_value * truth_value
            relative_sums["bias_squared"] += bias_squared / denominator
            relative_sums["variance"] += variance / denominator
            relative_sums["mse"] += mse / denominator
    absolute = {key: value / n_bins for key, value in absolute_sums.items()}
    if not valid_relative:
        raise ValueError("No truth bins are valid for relative toy metrics")
    relative = {
        key: value / len(valid_relative) for key, value in relative_sums.items()
    }
    return ToyMetricResult(
        iterations, outputs["Mean"], outputs["Bias"], outputs["BiasSquared"],
        outputs["Variance"], outputs["MSE"], absolute, relative,
        tuple(valid_relative),
    )


def scan_bayes_regularization(
    roounfold_module,
    response_bundle: RooUnfoldResponseBundle,
    toys: Sequence,
    truth,
    *,
    max_iterations: int,
    selection_mode: str = "absolute",
    handle_fakes: bool = True,
    efficiency=None,
    progress=None,
    unfolded_callback=None,
) -> RegularizationScanResult:
    """Unfold a common toy ensemble at every iteration and summarize metrics.

    Pass purity-corrected toys with a matched-only response and
    ``handle_fakes=False`` for factorized unfolding.  In that mode ``efficiency``
    converts every unfolded matched-truth toy to inclusive truth before bias and
    variance are evaluated, matching the observable represented by ``truth``.
    ``selection_mode='absolute'`` selects the minimum mean MSE in the actual
    bin contents. ``'relative'`` instead gives every non-negligible truth bin
    equal weight after division by its own truth squared. ``'both'`` records
    both minima and uses the absolute minimum as the primary iteration for
    code paths that require one result.
    """

    if max_iterations <= 0:
        raise ValueError("max_iterations must be positive")
    if selection_mode not in {"absolute", "relative", "both"}:
        raise ValueError("selection_mode must be 'absolute', 'relative', or 'both'")
    if efficiency is not None:
        validate_common_binning([truth, efficiency], label="regularization efficiency")
    metrics = []
    for iterations in range(1, max_iterations + 1):
        if progress is not None:
            progress(iterations, max_iterations)
        n_bins = truth.GetNbinsX()
        sums = [0.0] * (n_bins + 1)
        sums_squared = [0.0] * (n_bins + 1)
        for toy_index, toy in enumerate(toys):
            unfolded = unfold_bayes(
                roounfold_module, response_bundle, toy, iterations=iterations,
                handle_fakes=handle_fakes,
                suppress_native_output=True,
                name=f"hUnfoldedToy{toy_index}_iter{iterations}",
            ).histogram
            if efficiency is not None:
                for bin_index in range(1, n_bins + 1):
                    bin_efficiency = efficiency.GetBinContent(bin_index)
                    if bin_efficiency <= 0.0:
                        raise ValueError(
                            f"Cannot efficiency-correct regularization truth bin "
                            f"{bin_index}: efficiency={bin_efficiency:g}"
                        )
                    unfolded.SetBinContent(
                        bin_index,
                        unfolded.GetBinContent(bin_index) / bin_efficiency,
                    )
                    unfolded.SetBinError(
                        bin_index,
                        unfolded.GetBinError(bin_index) / bin_efficiency,
                    )
            if unfolded_callback is not None:
                unfolded_callback(iterations, toy_index, unfolded)
            for bin_index in range(1, n_bins + 1):
                value = unfolded.GetBinContent(bin_index)
                sums[bin_index] += value
                sums_squared[bin_index] += value * value
        metrics.append(_toy_metrics_from_moments(
            truth, sums, sums_squared, len(toys), iterations=iterations,
        ))
    summaries = {}
    for mode in ("absolute", "relative"):
        histogram = ROOT.TH2D(
            f"hRegularization{mode.capitalize()}Summary",
            f";Bayesian iterations;{mode.capitalize()} metric",
            max_iterations, 0.5, max_iterations + 0.5, 3, 0.5, 3.5,
        )
        histogram.SetDirectory(0)
        for metric_index, key in enumerate(("bias_squared", "variance", "mse"), 1):
            histogram.GetYaxis().SetBinLabel(metric_index, key)
            for result in metrics:
                histogram.SetBinContent(
                    result.iterations, metric_index, getattr(result, mode)[key]
                )
        summaries[mode] = histogram
    absolute_selected = min(
        metrics, key=lambda result: result.absolute["mse"]
    ).iterations
    relative_selected = min(
        metrics, key=lambda result: result.relative["mse"]
    ).iterations
    selected = (relative_selected if selection_mode == "relative"
                else absolute_selected)
    return RegularizationScanResult(
        tuple(metrics), summaries["absolute"], summaries["relative"], selected,
        selection_mode, absolute_selected, relative_selected,
    )


def load_unfolding_inputs(filename: str | Path,
                          keys: UnfoldingInputKeys) -> UnfoldingInputs:
    """Load and validate a named response-training input group."""

    objects = {
        field: load_object(str(filename), key)
        for field, key in (
            ("truth", keys.truth), ("measured", keys.measured),
            ("response", keys.response), ("miss", keys.miss), ("fake", keys.fake),
        )
    }
    classification = (
        load_object(str(filename), keys.classification) if keys.classification else None
    )
    for field in ("truth", "measured", "miss", "fake"):
        if not objects[field].InheritsFrom("TH2"):
            raise TypeError(f"{getattr(keys, field)} is not a TH2")
    if not objects["response"].InheritsFrom("THnSparse"):
        raise TypeError(f"{keys.response} is not a THnSparse")
    if classification is not None and not classification.InheritsFrom("TH1"):
        raise TypeError(f"{keys.classification} is not a TH1")
    return UnfoldingInputs(classification=classification, keys=keys, **objects)


def write_unfolding_output(
    filename: str | Path,
    *,
    histograms: Sequence,
    covariance=None,
    covariance_name: str = "hUnfoldedCovariance",
    response=None,
    response_name: str = "responseEtaCM",
    metadata: Mapping | None = None,
    metadata_name: str = "unfoldingConfiguration",
) -> None:
    """Write unfolding artifacts and JSON configuration to a ROOT file."""

    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    output = open_root_file(str(path), "RECREATE")
    try:
        output.cd()
        for histogram in histograms:
            if histogram is not None and histogram.Write() <= 0:
                raise OSError(f"Failed to write {histogram.GetName()} to {path}")
        if covariance is not None and covariance.Write(covariance_name) <= 0:
            raise OSError(f"Failed to write {covariance_name} to {path}")
        if response is not None and response.Write(response_name) <= 0:
            raise OSError(f"Failed to write {response_name} to {path}")
        if metadata is not None:
            payload = json.dumps(metadata, sort_keys=True, separators=(",", ":"))
            if ROOT.TNamed(metadata_name, payload).Write() <= 0:
                raise OSError(f"Failed to write {metadata_name} to {path}")
    finally:
        output.Close()


def _normalized_metadata(metadata: Mapping) -> dict:
    """Normalize tuples and other JSON-compatible containers for comparison."""

    return json.loads(json.dumps(metadata, sort_keys=True))


def _read_cache_metadata(root_file, *, metadata_name: str) -> dict | None:
    metadata_object = root_file.Get(metadata_name)
    if not metadata_object:
        return None
    try:
        return json.loads(metadata_object.GetTitle())
    except (TypeError, json.JSONDecodeError):
        return None


def _cache_matches(actual: Mapping | None, expected: Mapping) -> bool:
    """Require every expected metadata field to match the stored value."""

    if actual is None:
        return False
    expected_normalized = _normalized_metadata(expected)
    return all(actual.get(key) == value for key, value in expected_normalized.items())


class RegularizationDistributionWriter:
    """Stream complete regularization inputs and unfolded toys to one ROOT file.

    The file is first written with a ``.tmp`` suffix and becomes visible at the
    requested path only after :meth:`finalize` records every iteration metric.
    This avoids retaining all unfolded toys in Python memory and prevents a
    partially written file from being mistaken for a reusable result.
    """

    def __init__(
        self,
        filename: str | Path,
        *,
        histogram_groups: Mapping[str, Sequence],
        metadata: Mapping,
    ) -> None:
        self.path = Path(filename)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.temporary = self.path.with_suffix(self.path.suffix + ".tmp")
        self.metadata = _normalized_metadata(metadata)
        self.output = open_root_file(str(self.temporary), "RECREATE")
        self._closed = False
        try:
            for group_name, histograms in histogram_groups.items():
                directory = self.output.mkdir(group_name)
                if not directory:
                    raise OSError(f"Failed to create ROOT directory {group_name}")
                directory.cd()
                for histogram in histograms:
                    if histogram is not None and histogram.Write() <= 0:
                        raise OSError(
                            f"Failed to write {histogram.GetName()} to {self.temporary}"
                        )
            self.output.cd()
            payload = json.dumps(
                {**self.metadata, "complete": False},
                sort_keys=True, separators=(",", ":"),
            )
            ROOT.TNamed("cacheConfiguration", payload).Write()
        except Exception:
            self.abort()
            raise

    def _directory(self, path: str):
        directory = self.output
        for component in path.split("/"):
            child = directory.GetDirectory(component)
            directory = child or directory.mkdir(component)
            if not directory:
                raise OSError(f"Failed to create ROOT directory {path}")
        return directory

    def write_unfolded(self, iterations: int, toy_index: int, histogram) -> None:
        """Write one efficiency-corrected unfolded toy into its iteration folder."""

        if self._closed:
            raise RuntimeError("Cannot write to a closed distribution store")
        directory_name = f"iteration_{iterations:03d}/unfolded_toys"
        directory = self._directory(directory_name)
        directory.cd()
        name = f"hUnfoldedToy_{toy_index:06d}"
        if histogram.Write(name) <= 0:
            raise OSError(f"Failed to write {directory_name}/{name}")

    def finalize(self, scan: RegularizationScanResult) -> None:
        """Write metric distributions, mark the file complete, and install it."""

        if self._closed:
            raise RuntimeError("Distribution store is already closed")
        for result in scan.metrics:
            directory_name = f"iteration_{result.iterations:03d}/metrics"
            directory = self._directory(directory_name)
            directory.cd()
            for histogram in (
                result.mean, result.bias, result.bias_squared,
                result.variance, result.mse,
            ):
                if histogram.Write() <= 0:
                    raise OSError(
                        f"Failed to write {histogram.GetName()} to {directory_name}"
                    )
        summaries = self.output.mkdir("summaries")
        summaries.cd()
        scan.absolute_summary.Write()
        scan.relative_summary.Write()
        self.output.cd()
        payload = json.dumps(
            {
                **self.metadata,
                "complete": True,
                "selected_iterations": scan.selected_iterations,
                "absolute_selected_iterations": scan.absolute_selected_iterations,
                "relative_selected_iterations": scan.relative_selected_iterations,
                "stored_iterations": len(scan.metrics),
            },
            sort_keys=True, separators=(",", ":"),
        )
        ROOT.TNamed("cacheConfiguration", payload).Write(
            "cacheConfiguration", ROOT.TObject.kOverwrite,
        )
        self.output.Close()
        self._closed = True
        os.replace(self.temporary, self.path)

    def abort(self) -> None:
        """Close and remove only this writer's incomplete temporary file."""

        if not self._closed:
            self.output.Close()
            self._closed = True
        self.temporary.unlink(missing_ok=True)


def regularization_distribution_file_matches(
    filename: str | Path,
    *,
    metadata: Mapping,
    max_iterations: int,
    n_toys: int,
) -> bool:
    """Check metadata and terminal objects of a completed distribution file."""

    path = Path(filename)
    if not path.is_file():
        return False
    root_file = open_root_file(str(path), "READ")
    try:
        actual = _read_cache_metadata(root_file, metadata_name="cacheConfiguration")
        expected = {**metadata, "complete": True}
        if not _cache_matches(actual, expected):
            return False
        if actual.get("stored_iterations") != max_iterations:
            return False
        return bool(root_file.Get(
            f"iteration_{max_iterations:03d}/unfolded_toys/"
            f"hUnfoldedToy_{n_toys - 1:06d}"
        ) and root_file.Get(
            f"iteration_{max_iterations:03d}/metrics/"
            f"hRegularizationMSE_iter{max_iterations}"
        ))
    finally:
        root_file.Close()


def write_toy_cache(filename: str | Path, toys: Sequence, metadata: Mapping) -> None:
    """Atomically store generated toy distributions and their provenance."""

    path = Path(filename)
    temporary = path.with_suffix(path.suffix + ".tmp")
    write_unfolding_output(
        temporary, histograms=toys,
        metadata={**metadata, "cache_kind": "regularization_toys"},
        metadata_name="cacheConfiguration",
    )
    os.replace(temporary, path)


def load_toy_cache(
    filename: str | Path,
    *,
    n_toys: int,
    metadata: Mapping,
    name_prefix: str = "hRecoToy",
) -> list | None:
    """Load a compatible toy cache, or return ``None`` on any mismatch."""

    path = Path(filename)
    if not path.is_file():
        return None
    root_file = open_root_file(str(path), "READ")
    try:
        expected = {**metadata, "cache_kind": "regularization_toys"}
        if not _cache_matches(
            _read_cache_metadata(root_file, metadata_name="cacheConfiguration"),
            expected,
        ):
            return None
        toys = []
        for toy_index in range(n_toys):
            source = root_file.Get(f"{name_prefix}_{toy_index}")
            if not source or not source.InheritsFrom("TH1"):
                return None
            toy = source.Clone(f"{name_prefix}_{toy_index}")
            toy.SetDirectory(0)
            toys.append(toy)
        return toys
    finally:
        root_file.Close()


def write_regularization_scan_cache(
    filename: str | Path,
    scan: RegularizationScanResult,
    metadata: Mapping,
) -> None:
    """Atomically store all per-iteration distributions needed to reuse a scan."""

    path = Path(filename)
    temporary = path.with_suffix(path.suffix + ".tmp")
    metric_histograms = tuple(
        histogram
        for result in scan.metrics
        for histogram in (
            result.mean, result.bias, result.bias_squared, result.variance, result.mse,
        )
    )
    cache_metadata = {
        **metadata,
        "cache_kind": "regularization_scan",
        "selected_iterations": scan.selected_iterations,
        "selection_mode": scan.selection_mode,
        "absolute_selected_iterations": scan.absolute_selected_iterations,
        "relative_selected_iterations": scan.relative_selected_iterations,
        "valid_relative_bins": list(scan.metrics[0].valid_relative_bins),
    }
    write_unfolding_output(
        temporary,
        histograms=(
            *metric_histograms, scan.absolute_summary, scan.relative_summary,
        ),
        metadata=cache_metadata,
        metadata_name="cacheConfiguration",
    )
    os.replace(temporary, path)


def load_regularization_scan_cache(
    filename: str | Path,
    *,
    max_iterations: int,
    metadata: Mapping,
) -> RegularizationScanResult | None:
    """Load compatible per-iteration scan distributions without unfolding again."""

    path = Path(filename)
    if not path.is_file():
        return None
    root_file = open_root_file(str(path), "READ")
    try:
        actual = _read_cache_metadata(root_file, metadata_name="cacheConfiguration")
        expected = {**metadata, "cache_kind": "regularization_scan"}
        if not _cache_matches(actual, expected):
            return None
        valid_bins = tuple(int(value) for value in actual["valid_relative_bins"])
        metrics = []
        keys = ("Mean", "Bias", "BiasSquared", "Variance", "MSE")
        summary_keys = ("bias_squared", "variance", "mse")
        absolute_summary_source = root_file.Get("hRegularizationAbsoluteSummary")
        relative_summary_source = root_file.Get("hRegularizationRelativeSummary")
        if not absolute_summary_source or not relative_summary_source:
            return None
        absolute_summary = absolute_summary_source.Clone()
        relative_summary = relative_summary_source.Clone()
        for histogram in (absolute_summary, relative_summary):
            histogram.SetDirectory(0)
        for iterations in range(1, max_iterations + 1):
            histograms = []
            for key in keys:
                source = root_file.Get(f"hRegularization{key}_iter{iterations}")
                if not source:
                    return None
                histogram = source.Clone()
                histogram.SetDirectory(0)
                histograms.append(histogram)
            absolute = {
                key: absolute_summary.GetBinContent(iterations, index)
                for index, key in enumerate(summary_keys, 1)
            }
            relative = {
                key: relative_summary.GetBinContent(iterations, index)
                for index, key in enumerate(summary_keys, 1)
            }
            metrics.append(ToyMetricResult(
                iterations, *histograms, absolute, relative, valid_bins,
            ))
        return RegularizationScanResult(
            tuple(metrics), absolute_summary, relative_summary,
            int(actual["selected_iterations"]), actual.get("selection_mode", "relative"),
            int(actual.get("absolute_selected_iterations")
                or actual["selected_iterations"]),
            int(actual.get("relative_selected_iterations")
                or actual["selected_iterations"]),
        )
    finally:
        root_file.Close()
