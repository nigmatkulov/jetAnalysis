"""Eta-dependent hybrid JER inputs and constrained same-event response."""

from array import array
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import ROOT

from hist_analysis.config.histograms import DIJET_PTAVE_BINS
from hist_analysis.python.dijet_unfolding_workflow import (
    BASE_DIR, CUT, ETA_RANGE, check_partition,
)
from hist_analysis.python.histogram_io import (
    load_object, open_root_file, resolve_unfolding_mc_file,
)
from hist_analysis.python.unfolding_validation import (
    _check_pt, edges, load_training, matrix_values, values, variances,
)


PARTITIONS = {"full": "", "controlled": "ControlledSample", "studied": "StudiedSample"}
SIDEBANDS = (
    "fiducial_below_range", "fiducial_other_flow",
    "out_of_selection_below_range", "out_of_selection_in_range",
    "out_of_selection_other_flow",
)


def load_eta_training(filename=None, *, sample="full", intervals=DIJET_PTAVE_BINS):
    """Load compatible geometric and same-event eta-dependent inputs."""
    path = (resolve_unfolding_mc_file(BASE_DIR, "embedding", "combined")
            if filename is None else Path(filename).expanduser().resolve())
    split = PARTITIONS[sample]
    with open_root_file(str(path)) as handle:
        mode, eta = handle.Get("jerSmearingMode"), handle.Get("useEtaDependentJER")
        if not mode or mode.GetTitle() != "hybrid" or not eta or eta.GetTitle() not in ("true", "false"):
            raise ValueError(f"{path}: require documented hybrid JER and eta-dependence metadata")
        eta_all = eta.GetTitle() == "true"
        if split and not eta_all:
            raise ValueError(
                f"{path}: {sample} has base-resolution JER. Supply an eta-enabled "
                "hybrid production with useEtaDependentJER=true for split closure."
            )
        extended_suffix = split if eta_all else "JerDefExtra"
        reco_base = "hRecoDijetPtEtaCMJerDefUnfold" if eta_all else "hRecoDijetPtEtaCMJerDefExtraUnfold"
        matrix_base = ("hGenDijetPtEtaCMVsRecoPtEtaCM" if eta_all
                       else "hGenDijetPtEtaCMVsRecoJerDefExtraPtEtaCM")
        keys = tuple(f"{prefix}{split}_{CUT}" for prefix in
                     (reco_base, "hGenDijetPtEtaCM", matrix_base))
        extended_key = f"hUnfoldingExtendedResponse{extended_suffix}_{CUT}"
        if not handle.GetKey(extended_key):
            raise KeyError(f"{path}: missing eta-dependent same-event object {extended_key}")
    training = load_training(path, keys, intervals, ETA_RANGE, selection="jetId")
    training.metadata["eta_dependent_jer"] = dict(
        mode="hybrid", source_switch=eta_all, sample=sample,
        extended_suffix=extended_suffix, extended_response=extended_key,
    )
    return training


@dataclass(frozen=True)
class Response:
    """Weighted event counts and a forward model for one truth definition."""

    migration: np.ndarray
    migration_sumw2: np.ndarray
    truth: np.ndarray
    truth_sumw2: np.ndarray
    measured: np.ndarray
    measured_sumw2: np.ndarray
    residual: np.ndarray
    residual_sumw2: np.ndarray
    labels: tuple[str, ...]
    n_fiducial: int

    @property
    def kernel(self):
        return np.divide(self.migration, self.truth[None, :],
                         out=np.zeros_like(self.migration), where=self.truth[None, :] > 0)

    @property
    def efficiency(self):
        return self.kernel.sum(axis=0)

    def components(self):
        """Templates at the training exposure, including the residual reco term."""
        return np.column_stack((self.migration[:, self.n_fiducial:], self.residual))

    def background(self, rates=None):
        components = self.components()
        if rates is None:
            return components.sum(axis=1)
        rates = np.asarray(rates, dtype=float)
        if rates.shape != (components.shape[1],) or np.any(rates < 0) or not np.isfinite(rates).all():
            raise ValueError("Provide a finite nonnegative rate for each background component")
        return components @ rates

    def fold(self, fiducial, *, background):
        return self.kernel[:, :self.n_fiducial] @ fiducial + background


def _remainder(total, part, *, label, rtol=2e-7):
    delta = total - part
    if np.any(delta < -rtol * total):
        raise ValueError(f"{label}: response exceeds its inclusive population")
    return np.maximum(delta, 0.)


def geometric_response(training):
    migration = matrix_values(training.matrix)
    sumw2 = matrix_values(training.matrix, errors=True)**2
    reco, reco2 = values(training.reco), variances(training.reco)
    gen, gen2 = values(training.gen), variances(training.gen)
    _remainder(gen, migration.sum(0), label="geometric truth", rtol=1e-5)
    _remainder(gen2, sumw2.sum(0), label="geometric truth Sumw2")
    return Response(migration, sumw2, gen, gen2, reco, reco2,
                    _remainder(reco, migration.sum(1), label="geometric reco", rtol=1e-5),
                    _remainder(reco2, sumw2.sum(1), label="geometric reco Sumw2"),
                    ("geometric_reco_complement",), len(gen))


def _bin_maps(sparse, training):
    reference_eta = edges(training.templates[0].GetXaxis())
    maps = []
    for dim in (0, 2):
        pt, eta = sparse.GetAxis(dim), sparse.GetAxis(dim + 1)
        _check_pt(pt, training.layout.pt_intervals)
        pt_bins = {
            index: i for i, (low, high) in enumerate(training.layout.pt_intervals)
            for index in range(pt.FindFixBin(low + .001), pt.FindFixBin(high - .001) + 1)
        }
        eta_bins = {eta.FindFixBin(x): i for i, x in enumerate(
            (reference_eta[1:] + reference_eta[:-1]) / 2)}
        stored_eta = edges(eta)
        indexed_eta = np.array([stored_eta[b - 1] for b in eta_bins] +
                               [stored_eta[next(reversed(eta_bins))]])
        np.testing.assert_allclose(indexed_eta, reference_eta, rtol=0, atol=1e-9)
        maps.append({(x, y): i * training.layout.n_eta_bins + j
                     for x, i in pt_bins.items() for y, j in eta_bins.items()})
    return maps


def _cells(sparse):
    coordinates = array("i", [0] * sparse.GetNdimensions())
    for index in range(sparse.GetNbins()):
        weight, sumw2 = sparse.GetBinContent(index, coordinates), sparse.GetBinError2(index)
        if not (np.isfinite(weight) and np.isfinite(sumw2) and weight >= 0 and sumw2 >= 0):
            raise ValueError(f"{sparse.GetName()}: invalid weighted response cell")
        yield coordinates, weight, sumw2


def _check_axes(hist, expected):
    if (not hist.InheritsFrom("THnSparse") or hist.GetNdimensions() != len(expected)
            or not hist.GetCalculateErrors()):
        raise ValueError(f"{hist.GetName()}: expected sparse axes and Sumw2")
    for i, axis in enumerate(expected):
        current = hist.GetAxis(i)
        if current.TestBit(ROOT.TAxis.kAxisRange):
            raise ValueError(f"{hist.GetName()}: active source axis range")
        np.testing.assert_allclose(edges(current), edges(axis), rtol=0, atol=1e-9)


def extended_response(training):
    """Aggregate guard populations before conditioning on measured selection."""
    source = training.metadata["identity"]["source"]
    suffix = training.metadata["eta_dependent_jer"]["extended_suffix"]
    with open_root_file(source) as f:
        contract = f.Get("unfoldingExtendedContract")
        if not contract or "truth category 1=outside fiducial,2=fiducial" not in contract.GetTitle():
            raise ValueError(f"{source}: missing same-event category contract")
    sparse, inclusive, misses, measured, fake = (
        load_object(source, f"hUnfoldingExtended{name}{suffix}_{CUT}") for name in
        ("Response", "Truth", "Miss", "Measured", "Fake")
    )
    if not sparse.InheritsFrom("THnSparse") or sparse.GetNdimensions() != 5:
        raise ValueError("Expected a five-dimensional same-event response")
    axes = [sparse.GetAxis(i) for i in range(5)]
    _check_axes(sparse, axes)
    _check_axes(inclusive, (axes[0], axes[1], axes[4]))
    _check_axes(misses, (axes[0], axes[1], axes[4]))
    if axes[4].GetNbins() != 2 or [axes[4].GetBinCenter(i) for i in (1, 2)] != [1., 2.]:
        raise ValueError("Expected truth categories 1 and 2")
    for h in (measured, fake):
        if h.GetDimension() != 2 or not h.GetSumw2N():
            raise ValueError(f"{h.GetName()}: expected a weighted reco TH2")
        for actual, expected in zip((h.GetXaxis(), h.GetYaxis()), axes[2:4]):
            np.testing.assert_allclose(edges(actual), edges(expected), rtol=0, atol=1e-9)
    gen_map, reco_map = _bin_maps(sparse, training)
    lower = axes[0].FindFixBin(training.layout.pt_intervals[0][0] + .001)
    n = training.layout.n_global_bins
    nt = n + len(SIDEBANDS)

    def truth_index(x, y, category):
        if category not in (1, 2):
            raise ValueError("Invalid truth category, including flow")
        active = gen_map.get((x, y))
        if category == 2:
            return active if active is not None else n + (0 if x < lower else 1)
        return n + (2 if x < lower else 3 if active is not None else 4)

    total = np.zeros((2, nt))
    missed = np.zeros_like(total)
    associated = np.zeros_like(total)
    migration = np.zeros((2, n, nt))
    for hist, target in ((inclusive, total), (misses, missed)):
        for (x, y, category), weight, sumw2 in _cells(hist):
            target[:, truth_index(x, y, category)] += (weight, sumw2)
    for (tx, ty, rx, ry, category), weight, sumw2 in _cells(sparse):
        truth_bin = truth_index(tx, ty, category)
        associated[:, truth_bin] += (weight, sumw2)
        reco_bin = reco_map.get((rx, ry))
        if reco_bin is not None:
            migration[:, reco_bin, truth_bin] += (weight, sumw2)
    np.testing.assert_allclose(total, associated + missed, rtol=2e-10, atol=0)
    np.testing.assert_allclose(total[:, :n], [values(training.gen), variances(training.gen)],
                               rtol=2e-7, atol=0)

    def reco_moments(h):
        moments = np.zeros((2, n))
        for (x, y), index in reco_map.items():
            moments[:, index] += (h.GetBinContent(x, y), h.GetBinError(x, y)**2)
        return moments

    reco, residual = reco_moments(measured), reco_moments(fake)
    np.testing.assert_allclose(reco, migration.sum(axis=2) + residual, rtol=2e-10, atol=0)
    np.testing.assert_allclose(reco, [values(training.reco), variances(training.reco)],
                               rtol=2e-7, atol=0)
    _remainder(total, migration.sum(axis=1), label="extended truth")
    return Response(*migration, *total, *reco, *residual, SIDEBANDS + ("no_gen_pair",), n)


def same_event_fiducial(extended):
    """Collapse sidebands to a fixed background for a controlled comparison."""
    end = extended.n_fiducial
    return replace(
        extended,
        migration=extended.migration[:, :end],
        migration_sumw2=extended.migration_sumw2[:, :end],
        truth=extended.truth[:end], truth_sumw2=extended.truth_sumw2[:end],
        residual=extended.background(),
        residual_sumw2=extended.migration_sumw2[:, end:].sum(axis=1) + extended.residual_sumw2,
        labels=("same_event_reco_complement",),
    )


def check_eta_partition(full, controlled, studied):
    """Check both weighted response populations before split closure."""
    check_partition(full[0], controlled[0], studied[0])
    for field in ("migration", "migration_sumw2", "truth", "truth_sumw2",
                  "measured", "measured_sumw2", "residual", "residual_sumw2"):
        np.testing.assert_allclose(getattr(full[1], field),
                                   getattr(controlled[1], field) + getattr(studied[1], field),
                                   rtol=2e-10, atol=0)


def iterate(kernel, measured, background, start, iterations, *, derivative=False,
            background_templates=None):
    """Fixed-background D'Agostini iteration and optional data Jacobian."""
    kernel, measured, background, start = (
        np.asarray(x, dtype=float) for x in (kernel, measured, background, start)
    )
    if (kernel.ndim != 2 or start.shape != (kernel.shape[1],)
            or measured.shape != (kernel.shape[0],) or background.shape != measured.shape
            or not isinstance(iterations, int) or iterations < 1):
        raise ValueError("Incompatible response, spectra or iteration count")
    if any(not np.isfinite(x).all() or np.any(x < 0) for x in
           (kernel, measured, background, start)):
        raise ValueError("D'Agostini iteration requires finite nonnegative inputs")
    efficiency = kernel.sum(axis=0)
    support = efficiency > 0
    if np.any(efficiency > 1 + 1e-5) or np.any(start[support] <= 0):
        raise ValueError("Invalid efficiency or zero start on measured truth support")
    empty = (kernel.sum(axis=1) == 0) & (background == 0)
    if np.any(empty & (measured > 0)):
        raise ValueError("Observed counts lie outside the forward model support")
    estimate = np.where(support, start, 0.)
    inv_eff = np.divide(1., efficiency, out=np.zeros_like(efficiency), where=support)
    if derivative:
        jacobian = np.zeros((len(start), len(measured)))
    if background_templates is not None:
        background_templates = np.asarray(background_templates, dtype=float)
        if (not derivative or background_templates.ndim != 2
                or background_templates.shape[0] != len(measured)
                or not np.isfinite(background_templates).all()):
            raise ValueError("Background templates require derivative=True and measured rows")
        nuisance_jacobian = np.zeros((len(start), background_templates.shape[1]))
    for _ in range(iterations):
        prediction = kernel @ estimate + background
        if np.any((prediction <= 0) & (measured > 0)):
            raise ValueError("An iterative prediction is zero for a populated reco bin")
        inverse = np.divide(1., prediction, out=np.zeros_like(prediction), where=prediction > 0)
        correction = (kernel.T @ (measured * inverse)) * inv_eff
        if derivative:
            step = ((estimate * inv_eff)[:, None] * kernel.T) * inverse
            slope = step * (measured * inverse)
            jacobian = correction[:, None] * jacobian + step - slope @ (kernel @ jacobian)
            if background_templates is not None:
                nuisance_jacobian = (correction[:, None] * nuisance_jacobian
                                     - slope @ (kernel @ nuisance_jacobian + background_templates))
        estimate = estimate * correction
    if background_templates is not None:
        return estimate, jacobian, nuisance_jacobian
    return (estimate, jacobian) if derivative else estimate
