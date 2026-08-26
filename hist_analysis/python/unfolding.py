"""Reusable PyROOT helpers for flattened dijet pTave-eta unfolding."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Mapping, Sequence

import ROOT

from hist_analysis.python.histogram_io import load_object, open_root_file
from hist_analysis.python.projections import project_semantic_th2


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


def project_eta_by_pt(histogram, pt_bins: Sequence, *, name_prefix: str) -> list:
    """Project a TH2 eta distribution in configured half-open pT intervals."""

    if not histogram.InheritsFrom("TH2"):
        raise TypeError(f"{histogram.GetName()} is not a TH2")
    intervals = as_pt_intervals(pt_bins)
    return [
        project_semantic_th2(
            histogram, "eta", interval, name=f"{name_prefix}_pt{index}",
        )
        for index, interval in enumerate(intervals)
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
    n_gen_eta = sparse_response.GetAxis(gen_eta_axis).GetNbins()
    n_reco_eta = sparse_response.GetAxis(reco_eta_axis).GetNbins()
    if n_gen_eta != n_reco_eta:
        raise ValueError("Gen and reco response eta axes have different bin counts")
    if not _same_axis_binning(
        sparse_response.GetAxis(gen_eta_axis), sparse_response.GetAxis(reco_eta_axis)
    ):
        raise ValueError("Gen and reco response eta axes have different bin edges")
    layout = layout or FlattenedBinning(intervals, n_gen_eta)
    if layout.pt_intervals != intervals or layout.n_eta_bins != n_gen_eta:
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
                for gen_eta_bin in range(1, layout.n_eta_bins + 1):
                    global_gen = layout.global_bin(gen_pt_index, gen_eta_bin)
                    for reco_eta_bin in range(1, layout.n_eta_bins + 1):
                        global_reco = layout.global_bin(reco_pt_index, reco_eta_bin)
                        response.SetBinContent(
                            global_reco, global_gen,
                            block.GetBinContent(reco_eta_bin, gen_eta_bin),
                        )
                        response.SetBinError(
                            global_reco, global_gen,
                            block.GetBinError(reco_eta_bin, gen_eta_bin),
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
    """Calculate response marginals and effective/boundary miss/fake spectra."""

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
    """Build a scaled RooUnfoldResponse and validate its inferred fakes."""

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
    unfolded = algorithm.Hunfold().Clone(name)
    unfolded.SetTitle(title)
    unfolded.SetDirectory(0)
    unfolded.Scale(1.0 / response_bundle.scale)
    covariance = algorithm.Eunfold(ROOT.RooUnfolding.kCovariance)
    covariance *= 1.0 / (response_bundle.scale * response_bundle.scale)
    return UnfoldingResult(unfolded, covariance, algorithm, scaled_measured)


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
