"""Validated geometric response and conditional pPb dijet toy studies."""

from __future__ import annotations

import json
import os
from array import array
from pathlib import Path
import tempfile

import numpy as np

from hist_analysis.config.files import BASE_DIR
from hist_analysis.config.histograms import DIJET_PTAVE_BINS, STANDARD_DIJET_ETA_CUT_INDEX
from hist_analysis.python.data_trigger_stitching import load_trigger_stitching_config
from hist_analysis.python.histogram_io import (
    load_object, open_root_file, resolve_data_file, resolve_unfolding_mc_file,
    resolve_unfolding_truth_file,
)
from hist_analysis.python.notebook_setup import find_project_root, load_root, load_roounfold
from hist_analysis.python.plot_labels import (
    CMS_PRELIMINARY_HEADER, CMS_SIMULATION_HEADER, PPB_MPLHEP_HEADER_RIGHT,
)
from hist_analysis.python.unfolding import (
    flatten_sparse_response, validate_common_binning, write_unfolding_output,
)
from hist_analysis.python.unfolding_validation import (
    _audit_sparse, _bayes_derivative, divide, edges, flatten_spectrum, histogram, load_training,
    matrix_values, ratio, run_unfolding, save_validation, Spectrum, truth_refold_check, validation_views,
    values, variances,
)
from hist_analysis.python.unfolding_validation_plots import ValidationPlots

ETA_RANGE = (-1.9, 1.9)
CUT = STANDARD_DIJET_ETA_CUT_INDEX
ANNOTATIONS = (
    r"$|\eta_{\rm CM}^{\rm jet}|\leq 1.9$, $|\Delta\phi|\geq 2\pi/3$",
    r"$p_{\rm T,1}\geq 50$, $p_{\rm T,2}\geq 40$ GeV",
)


def setup(stage):
    root = find_project_root()
    runtime = load_root()
    runtime.TH1.AddDirectory(False)
    load_roounfold(runtime, project_root=root)
    destination = Path(os.environ.get(
        "DIJET_UNFOLDING_STUDY_OUTPUT", root / "hist_analysis/output/dijet_unfolding_study"
    )) / stage
    destination.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix="run_", dir=destination))


def training_input(sample="full", *, cache=None, reuse=False):
    suffix = {"full": "", "controlled": "ControlledSample", "studied": "StudiedSample"}[sample]
    keys = tuple(f"{key}{suffix}_{CUT}" for key in (
        "hRecoDijetPtEtaCMJerDefUnfold", "hGenDijetPtEtaCM",
        "hGenDijetPtEtaCMVsRecoPtEtaCM",
    ))
    path = resolve_unfolding_mc_file(BASE_DIR, "embedding", "combined")
    return load_training(path, keys, DIJET_PTAVE_BINS, ETA_RANGE,
                         selection="jetId", cache=cache, reuse=reuse)


def reco_gen_covariance(sample, *, suffix=""):
    """Read joint event Sumw2 from the fiducial same-event association."""
    path = sample.metadata["identity"]["source"]
    response_key = f"hUnfoldingExtendedResponse{suffix}_{CUT}"
    f = open_root_file(path)
    try:
        contract = f.Get("unfoldingExtendedContract")
        if not contract or "truth category 1=outside fiducial,2=fiducial" not in contract.GetTitle():
            raise ValueError("Missing documented same-event truth category")
    finally:
        f.Close()
    joint = load_object(path, response_key)
    gen = load_object(path, f"hUnfoldingExtendedTruth{suffix}_{CUT}")
    reco = load_object(path, f"hUnfoldingExtendedMeasured{suffix}_{CUT}")
    if (joint.GetNdimensions() != 5 or gen.GetNdimensions() != 3
            or not joint.GetCalculateErrors() or not gen.GetCalculateErrors()):
        raise ValueError("Unexpected same-event histogram schema")
    joint_axis, truth_axis = joint.GetAxis(4), gen.GetAxis(2)
    for axis in (joint_axis, truth_axis):
        fiducial = axis.FindBin(2.)
        if axis.GetNbins() != 2 or axis.GetBinCenter(fiducial) != 2.:
            raise ValueError("Invalid truth category binning")
        axis.SetRange(fiducial, fiducial)
    for role, h, expected in (("gen", gen.Projection(1, 0, "E"), sample.gen),
                              ("reco", reco, sample.reco)):
        flat, templates, layout = flatten_spectrum(h, sample.layout.pt_intervals, ETA_RANGE,
                                                    name=f"joint_{role}")
        if layout != sample.layout:
            raise ValueError("Same-event marginal differs from geometric spectrum binning")
        validate_common_binning([templates[0], sample.templates[0]], label=role)
        # Legacy inclusive fills may round an event weight to float first.
        np.testing.assert_allclose(values(flat), values(expected), rtol=2e-7, atol=0)
        np.testing.assert_allclose(variances(flat), variances(expected), rtol=2e-7, atol=0)
    fiducial_joint = joint.Projection(4, array("i", [0, 1, 2, 3]), "E")
    for axis, expected in enumerate((gen.GetAxis(0), gen.GetAxis(1), reco.GetXaxis(), reco.GetYaxis())):
        np.testing.assert_allclose(edges(fiducial_joint.GetAxis(axis)), edges(expected),
                                   rtol=0, atol=1e-9)
    matrix, _ = flatten_sparse_response(fiducial_joint, sample.layout.pt_intervals,
                                        name="joint_reco_gen", layout=sample.layout,
                                        eta_range=ETA_RANGE)
    _audit_sparse(fiducial_joint, sample.layout, edges(sample.templates[0].GetXaxis()), matrix)
    covariance = matrix_values(matrix, errors=True)**2
    tolerance = 1 + 2e-7
    if (np.any(covariance.sum(axis=0) > tolerance * variances(sample.gen))
            or np.any(covariance.sum(axis=1) > tolerance * variances(sample.reco))):
        raise ValueError("Same-event Sumw2 exceeds an inclusive variance")
    return covariance, _source_identity(path, response_key)


def plots(output, *, data=False):
    return ValidationPlots(output, ("pdf",),
                           left=CMS_PRELIMINARY_HEADER if data else CMS_SIMULATION_HEADER,
                           right=PPB_MPLHEP_HEADER_RIGHT, annotations=ANNOTATIONS)


def comparisons(plotter, kind, flat, intervals, *, prefix):
    plotter.comparisons(kind, flat, intervals, prefix=prefix, grids=bool(intervals))


def _blocks(training, spectrum):
    width = training.layout.n_eta_bins
    for i, (interval, template) in enumerate(zip(training.layout.pt_intervals, training.templates)):
        block = slice(i * width, (i + 1) * width)
        covariance = spectrum.covariance[block, block]
        h = histogram(template, f"{spectrum.histogram.GetName()}_pt{i}",
                      values(spectrum.histogram)[block], covariance.diagonal())
        yield interval, Spectrum(h, covariance, spectrum.supported[block])


def response_diagnostics(training, output, *, draw=True):
    t = training
    m = matrix_values(t.matrix)
    p = values(t.purity)
    inverse = values(t.inverse_efficiency)
    efficiency = divide(1., inverse)
    # Training already checks the subset relation. Clip only tolerated
    # numerical deviations in the weighted complements.
    fakes = histogram(t.reco, "domain_fake",
                      np.maximum(values(t.reco) - values(t.matched_reco), 0.),
                      np.maximum(variances(t.reco) - variances(t.matched_reco), 0.))
    misses = histogram(t.gen, "domain_miss",
                       np.maximum(values(t.gen) - values(t.matched_gen), 0.),
                       np.maximum(variances(t.gen) - variances(t.matched_gen), 0.))
    eh = histogram(t.gen, "efficiency", efficiency,
                   divide(variances(t.inverse_efficiency), inverse**4))
    diagonal = {"Reco bin purity": divide(m.diagonal(), m.sum(axis=1)),
                "Gen bin stability": divide(m.diagonal(), m.sum(axis=0))}
    weighted_cells = matrix_values(t.matrix, errors=True)**2
    effective = divide(m**2, weighted_cells)
    probabilities = t.probabilities
    singular_values = np.linalg.svd(probabilities, compute_uv=False)
    write_unfolding_output(output / "response.root", histograms=(
        t.matrix, t.reco, t.gen, t.matched_reco, t.matched_gen, fakes, misses, t.purity, eh,
    ), metadata=t.metadata)
    width = t.layout.n_eta_bins
    index = np.arange(t.layout.n_global_bins) // width
    eta_edges = edges(t.templates[0].GetXaxis())
    eta_centers = (eta_edges[1:] + eta_edges[:-1]) / 2
    eta = np.tile(eta_centers, t.layout.n_pt_bins)
    pt = np.repeat(np.mean(t.layout.pt_intervals, axis=1), width)
    # Resolution proxies from discretized matrix coordinates, not per-jet JER.
    d_eta = eta[:, None] - eta
    d_pt = pt[:, None] - pt
    residual_eta = np.histogram(d_eta, bins=76, weights=m)
    residual_pt = np.histogram(d_pt, bins=80, weights=m)
    np.savez_compressed(output / "migration.npz", matrix=m, response_sumw2=weighted_cells,
                        probabilities=probabilities, effective_entries=effective,
                        singular_values=singular_values, eta_edges=eta_edges,
                        pt_intervals=t.layout.pt_intervals, purity=p, efficiency=efficiency,
                        reco_bin_purity=diagonal["Reco bin purity"],
                        gen_bin_stability=diagonal["Gen bin stability"],
                        delta_eta=residual_eta[0], delta_eta_edges=residual_eta[1],
                        delta_pt=residual_pt[0], delta_pt_edges=residual_pt[1])
    summary = dict(matched_weight=float(m.sum()), off_diagonal_pt=float(
        m[index[:, None] != index].sum() / m.sum()),
        zero_efficiency=int((~np.isfinite(inverse)).sum()),
        matrix_rank=int(np.count_nonzero(singular_values > singular_values[0] * len(m)
                                         * np.finfo(float).eps)),
        singular_min=float(singular_values[-1]), singular_max=float(singular_values[0]),
        min_matched_gen_effective=float(np.nanmin(divide(values(t.matched_gen)**2,
                                                        variances(t.matched_gen)))))
    (output / "summary.json").write_text(json.dumps(summary, indent=2))
    if draw:
        plotter = plots(output / "mplhep")
        plotter.response(t)
        groups = {
            "eta": {label: _fixed(h) for label, h in (("Reco", t.reco), ("Gen", t.gen))},
            "purity": {"Purity": _fixed(t.purity)},
            "efficiency": {"Efficiency": _fixed(eh)},
            "bin_migration": {label: _fixed(histogram(t.gen, label.replace(" ", "_"), y,
                                                      np.zeros_like(y)))
                              for label, y in diagonal.items()},
        }
        for kind, curves in groups.items():
            blocks = {interval: {} for interval in t.layout.pt_intervals}
            for label, spectrum in curves.items():
                for interval, block in _blocks(t, spectrum):
                    blocks[interval][label] = block
            comparisons(plotter, kind, curves, blocks, prefix="response")
        for label, (content, bins) in (
            (r"Reco minus Gen dijet $\eta_{\mathrm{CM}}$", residual_eta),
            (r"Reco minus Gen $p_{\mathrm{T}}^{\mathrm{ave}}$ (GeV)", residual_pt),
        ):
            with plotter.plt.style.context(plotter.hep.style.CMS):
                fig, axis = plotter.plt.subplots(figsize=(8, 6))
                fig.subplots_adjust(top=.82, bottom=.16, left=.18)
                plotter._header(fig, left=.18, top=.91, small=True)
                axis.stairs(content, bins)
                axis.set(xlabel=label, ylabel="Matched weighted yield")
                axis.set_title("Residuals from response bin centers", fontsize=13)
                plotter._export(fig, "residual_eta" if content is residual_eta[0] else "residual_pt")
    return summary


def _fixed(h):
    y = values(h)
    return Spectrum(h, np.diag(np.nan_to_num(variances(h))), np.isfinite(y))


def closure(training, measured, truth, output, *, iterations, draw=True, data=False,
            provenance=None, reco_truth_cross=None):
    t = training
    result = run_unfolding(t, measured, iterations=iterations)
    # The estimate and its input share experimental or validation fluctuations.
    p = np.nan_to_num(values(t.purity))
    inverse = np.nan_to_num(values(t.inverse_efficiency))
    _, derivative = _bayes_derivative(t.probabilities, p * values(measured),
                                      values(t.matched_gen), iterations)
    jacobian = inverse[:, None] * derivative * p
    cross = jacobian @ result.reco.covariance
    unfolded_truth_cross = None if reco_truth_cross is None else jacobian @ reco_truth_cross
    flat, slices, metrics = validation_views(t, result, truth, reco_truth_cross=reco_truth_cross,
                                              unfolded_truth_cross=unfolded_truth_cross)
    r = ratio(result.unfolded, result.reco, name="unfolded_over_reco", cross=cross)
    flat["unfolded_over_reco"] = {"Unfolded / Reco": r}
    slices["unfolded_over_reco"] = {interval: {"Unfolded / Reco": spectrum}
                                    for interval, spectrum in _blocks(t, r)}
    save_validation(output / f"iterations_{iterations}.root", t, result, flat, slices, metrics,
                    provenance=dict(provenance or {}, iterations=iterations,
                                    reco_gen_cross=("same-event fiducial Sumw2"
                                                    if reco_truth_cross is not None
                                                    else "not supplied")))
    truth_covariances = {} if reco_truth_cross is None else dict(
        reco_gen_cross=reco_truth_cross, unfolded_gen_cross=unfolded_truth_cross)
    np.savez_compressed(output / f"covariance_iterations_{iterations}.npz",
                        unfolded_reco_cross=cross, unfolded_reco_ratio=r.covariance,
                        **truth_covariances)
    if draw:
        plotter = plots(output / "mplhep", data=data)
        for kind in ("eta", "to_gen", "unfolded_over_reco", "refolded", "refolded_ratio", "fb"):
            if kind in slices:
                comparisons(plotter, kind, flat.get(kind), slices[kind], prefix=f"iter{iterations}")
        closure_residual_plots(plotter, metrics, iterations)
    summary = {key: {k: v for k, v in metrics[key].items() if np.ndim(v) == 0}
               for key in ("inclusive", "closure") if key in metrics}
    for key, info in summary.items():
        pulls = metrics[key]["pulls"]
        finite = np.isfinite(pulls)
        info.update(max_abs_pull=float(np.max(np.abs(pulls[finite]), initial=0.)),
                    pulls_above_3=int((np.abs(pulls[finite]) > 3.).sum()))
    mask = result.refolded.supported
    summary["refold_integral_ratio"] = float(values(result.refolded.histogram)[mask].sum()
                                             / values(measured)[mask].sum())
    if truth is not None:
        mask = result.unfolded.supported
        summary["truth_integral_ratio"] = float(values(result.unfolded.histogram)[mask].sum()
                                                / values(truth)[mask].sum())
    return summary


def closure_residual_plots(plotter, metrics, iterations):
    """Plot conditional signed residuals and pull distributions."""
    labels = {"inclusive": "Refolded − Reco", "closure": "Unfolded − Gen"}
    with plotter.plt.style.context(plotter.hep.style.CMS):
        for kind in ("residual", "pulls"):
            fig, ax = plotter.plt.subplots(figsize=(8, 6))
            fig.subplots_adjust(top=.82, bottom=.16, left=.19)
            plotter._header(fig, left=.19, top=.91, small=True)
            samples = {label: metrics[key][kind] for key, label in labels.items() if key in metrics}
            if kind == "residual":
                for label, y in samples.items():
                    ax.plot(np.arange(1, len(y)+1), y, ".", markersize=3, label=label)
                ax.axhline(0., color=".55", ls="--")
                ax.set(xlabel="Flattened ROOT bin", ylabel="Weighted yield residual")
            else:
                pooled = np.concatenate([y[np.isfinite(y)] for y in samples.values()])
                if pooled.size:
                    bins = np.histogram_bin_edges(pooled, bins=50)
                    for label, y in samples.items():
                        ax.hist(y[np.isfinite(y)], bins=bins, histtype="step", label=label)
                ax.set(xlabel="Conditional pull", ylabel="Bins")
            ax.legend(fontsize=12)
            plotter._export(fig, f"iter{iterations}_{kind}")


def _flat_input(training, path, key, *, name):
    h, templates, layout = flatten_spectrum(load_object(str(path), key),
                                            training.layout.pt_intervals, ETA_RANGE, name=name)
    if layout != training.layout:
        raise ValueError(f"Incompatible flattened layout for {key}")
    validate_common_binning([templates[0], training.templates[0]], label=f"{key} eta bins")
    return h


def _source_identity(path, key):
    f = open_root_file(str(path))
    try:
        return dict(path=str(path), key=key, uuid=str(f.GetUUID().AsString()),
                    cycle=int(f.GetKey(key).GetCycle()))
    finally:
        f.Close()


def alternate_truth(training):
    path = resolve_unfolding_truth_file(BASE_DIR, "pythia_cp5_standalone", "pp")
    key = f"hGenDijetPtEtaCM_{CUT}"
    return _flat_input(training, path, key, name="cp5_pp_truth"), _source_identity(path, key)


def stitching_input():
    override = os.environ.get("TRIGGER_STITCHING_WEIGHTS_FILE")
    if override:
        return Path(override)
    directory = (find_project_root() / "hist_analysis/output/data_trigger_stitching"
                 / "2026-11-11/jetId/root")
    candidates = [p for p in directory.glob("*trigger_stitching_weights_*combined*.json")
                  if not p.name.endswith(".selection.json")]
    if len(candidates) != 1:
        raise ValueError("Set TRIGGER_STITCHING_WEIGHTS_FILE to a combined dijet artifact")
    return candidates[0]


def stitched_data(training, stitching_file):
    config = load_trigger_stitching_config(stitching_file, direction="combined",
                                           selection="jetId", eta_range=ETA_RANGE,
                                           eta_cut_index=CUT)
    files = {}
    provenance = {}
    for trigger, key in config.histogram_keys.items():
        path = resolve_data_file(Path(os.environ.get("PPB_DATA_DIR", BASE_DIR / "exp")), trigger)
        files[trigger] = _flat_input(training, path, key, name=trigger)
        provenance[trigger] = _source_identity(path, key)
        raw = values(files[trigger])
        var = variances(files[trigger])
        if (np.any(raw < 0) or not np.array_equal(raw, np.rint(raw))
                or not np.allclose(var, raw, rtol=1e-10, atol=0)):
            raise ValueError(f"{trigger} Sumw2 is incompatible with unit event counts")
    selected = []
    for lo, hi in training.layout.pt_intervals:
        for upper, trigger in ((80., "MinimumBias"), (100., "Jet60"),
                               (120., "Jet80"), (1000., "Jet100")):
            if lo < upper:
                if hi > upper:
                    raise ValueError(f"Interval {(lo, hi)} crosses a trigger boundary")
                selected.append(trigger)
                break
        else:
            raise ValueError(f"Missing trigger assignment for {(lo, hi)}")
    width = training.layout.n_eta_bins
    counts = np.concatenate([values(files[trigger])[i*width:(i+1)*width]
                             for i, trigger in enumerate(selected)])
    weights = np.repeat([config.weights[trigger] for trigger in selected], width)
    if not np.isfinite(weights).all():
        raise ValueError("Invalid trigger multiplier")
    h = histogram(training.reco, "stitched_reco", counts * weights, counts * weights**2)
    intervals = np.asarray(training.layout.pt_intervals)
    overlap_checks = {}
    for record in config.metadata["scales"]:
        low, high = record["pt_range"]
        if low not in intervals[:, 0] or high not in intervals[:, 1]:
            raise ValueError("Overlap boundaries must align with response intervals")
        block_mask = (intervals[:, 0] >= low) & (intervals[:, 1] <= high)
        mask = np.repeat(block_mask, width)
        reference = "MinimumBias" if record["reference"] == "MB" else record["reference"]
        candidate = record["candidate"]
        earlier, later = (values(files[stream]) for stream in (reference, candidate))
        for raw, field in ((earlier, "reference_integral"), (later, "candidate_integral")):
            if float(raw[mask].sum()) != float(record[field]):
                raise ValueError(f"Current inputs differ from the {candidate} overlap artifact")
        if candidate == "Jet60":
            overlap_mask, mb_yield, mb_range = mask, float(earlier[mask].sum()), [low, high]
        ratio_pt = divide(later.reshape(-1, width).sum(axis=1) * config.weights[candidate],
                          earlier.reshape(-1, width).sum(axis=1) * config.weights[reference])
        overlap_checks[candidate] = dict(reference=reference, normalization_range=[low, high],
                                         scaled_ratio_by_pt=ratio_pt.tolist())
    return h, weights, overlap_mask, mb_yield, dict(
        sources=provenance, artifact=str(config.source), trigger_intervals=selected,
        weights=config.weights, weight_uncertainties=config.uncertainties,
        mb_overlap=mb_range, mb_overlap_yield=mb_yield, overlap_checks=overlap_checks,
        covariance="conditional on fixed stream scales; diagonal counting covariance",
    )


def check_partition(full, controlled, studied):
    """Check stored entry complements in both weighted moments."""
    for field in ("reco", "gen", "matrix"):
        for errors in (False, True):
            arrays = [getattr(t, field) for t in (full, controlled, studied)]
            if field == "matrix":
                a, b, c = [matrix_values(h, errors=errors)**(2 if errors else 1)
                           for h in arrays]
            else:
                a, b, c = [(variances if errors else values)(h) for h in arrays]
            # THnSparseF sums weights at float precision in different orders.
            tolerance = 1e-5 if field == "matrix" and not errors else 1e-10
            np.testing.assert_allclose(a, b + c, rtol=tolerance, atol=0.)


def fold_expected(training, truth):
    """Apply inclusive response probabilities and the fixed Reco fake model."""
    t = training
    if (truth.shape != values(t.gen).shape or not np.isfinite(truth).all()
            or np.any(truth < 0) or truth.sum() <= 0):
        raise ValueError("Toy truth must have finite nonnegative yields and positive integral")
    m = matrix_values(t.matrix)
    gen, p, inv = values(t.gen), values(t.purity), values(t.inverse_efficiency)
    valid = np.isfinite(inv) & (inv > 0)
    if np.any(truth[~valid] > 0):
        raise ValueError("Alternative truth has weight in unmeasured Gen bins")
    inclusive_probability = np.divide(m, gen, out=np.zeros_like(m), where=gen > 0)
    signal = inclusive_probability @ truth
    supported = np.isfinite(p) & (p > 0)
    if np.any((~supported) & (values(t.reco) > 0)):
        raise ValueError("Training contains populated Reco bins with undefined fake restoration")
    efficiency = np.divide(1., inv, out=np.zeros_like(inv), where=valid)
    np.testing.assert_allclose(signal, t.probabilities @ (efficiency * truth),
                               rtol=1e-5, atol=1e-15 * signal.max())
    return np.divide(signal, p, out=np.zeros_like(signal), where=supported)


def poisson_toys(expected, weights, count, rng):
    """Independent Poisson counts transported onto the fixed MB exposure."""
    if (count < 2 or expected.shape != weights.shape or not np.isfinite(expected).all()
            or not np.isfinite(weights).all() or np.any(weights <= 0) or np.any(expected < 0)):
        raise ValueError("Invalid scaled Poisson toy configuration")
    return rng.poisson(expected / weights, size=(count, expected.size)) * weights


def bayes_toy_trajectory(training, measured, max_iterations):
    """Advance all toys together through the fixed-response Bayes recursion."""
    a = training.probabilities
    p = np.nan_to_num(values(training.purity))
    inverse = np.nan_to_num(values(training.inverse_efficiency))
    truth_prior = values(training.matched_gen)
    prior = np.broadcast_to(truth_prior / truth_prior.sum(), measured.shape)
    signal = measured * p
    for _ in range(max_iterations):
        predicted = prior @ a.T
        if np.any((predicted <= 0) & (signal > 0)):
            raise ValueError("Pseudodata falls outside the iterative response support")
        transport = np.divide(signal, predicted, out=np.zeros_like(signal), where=predicted > 0)
        prior = prior * (transport @ a)
        yield prior * inverse


def _fb_array(y, n_pt, width):
    a = y.reshape(*y.shape[:-1], n_pt, width)
    half = width // 2
    return divide(a[..., half:], a[..., :half][..., ::-1])


def risk_metrics(ensemble, truth):
    mean = ensemble.mean(axis=0)
    variance = ensemble.var(axis=0, ddof=0)
    bias = mean - truth
    mse = np.mean((ensemble - truth)**2, axis=0)
    np.testing.assert_allclose(mse, bias**2 + variance, rtol=1e-11, atol=0.)
    return mean, bias, variance, mse


def regularization(training, truth_hist, data_weights, overlap, mb_yield, output, *,
                   iterations, count, seed, precision=1., draw=True, provenance=None):
    iterations = tuple(sorted(set(iterations)))
    if (not iterations or any(not isinstance(i, int) or i < 1 for i in iterations)
            or not np.isfinite(precision) or precision <= 0):
        raise ValueError("Use positive iteration counts and precision")
    t = training
    raw_truth = values(truth_hist)
    raw_expected = fold_expected(t, raw_truth)
    if raw_expected[overlap].sum() <= 0 or mb_yield <= 0:
        raise ValueError("MB overlap must contain positive data and prediction weight")
    exposure = mb_yield / raw_expected[overlap].sum()
    weights = data_weights * precision
    truth, expected = raw_truth * exposure, raw_expected * exposure
    rng = np.random.default_rng(seed)
    toys = poisson_toys(expected, weights, count, rng)
    width, blocks = t.layout.n_eta_bins, t.layout.n_pt_bins
    eta_axis = edges(t.templates[0].GetXaxis())
    if width % 2 or not np.allclose(eta_axis, -eta_axis[::-1], rtol=0, atol=1e-9):
        raise ValueError("F/B requires mirrored eta edges with a boundary at zero")
    fb_truth = _fb_array(truth, blocks, width)
    # Within each pT block normalize risk by the squared mean truth bin yield.
    # Give blocks equal influence and retain even bins whose alternative truth is zero.
    reference = np.repeat(truth.reshape(blocks, width).mean(axis=1), width)
    support = np.isfinite(values(t.inverse_efficiency))
    if np.any(reference[support] <= 0):
        raise ValueError("Alternative truth must populate each scored pT interval")
    sample = histogram(t.reco, "toy_expectation", expected, expected * weights)
    # Verify the batched propagation on both the exact expectation and a
    # fluctuated toy. RooUnfold also provides nominal linearized covariance.
    trajectories = bayes_toy_trajectory(t, np.vstack([expected, toys]), max(iterations))
    arrays, risks, diagnostics = {}, [], []
    for iteration, all_estimates in enumerate(trajectories, 1):
        if iteration not in iterations:
            continue
        ensemble = all_estimates[1:]
        result = run_unfolding(t, sample, iterations=iteration)
        toy_hist = histogram(t.reco, "fluctuated_toy", toys[0], toys[0] * weights)
        toy_check = run_unfolding(t, toy_hist, iterations=iteration)
        for estimate, roo in ((all_estimates[0], result), (ensemble[0], toy_check)):
            np.testing.assert_allclose(estimate, values(roo.unfolded.histogram),
                                       rtol=1e-9, atol=1e-13 * expected.max())
        mean, bias, variance, mse = risk_metrics(ensemble, truth)
        sigma = np.sqrt(result.unfolded.covariance.diagonal())
        pulls = divide(ensemble - truth, sigma)
        coverage = np.where(support & (sigma > 0), np.mean(np.abs(pulls) <= 1, axis=0), np.nan)
        fb = _fb_array(ensemble, blocks, width)
        scaled_errors = (ensemble[:, support] - truth[support]) / reference[support]
        # Compute equal block weights even when response support differs by block.
        counts = support.reshape(blocks, width).sum(axis=1)
        if np.any(counts == 0):
            raise ValueError("Missing truth support for an entire pT interval")
        risk_weights = np.repeat(1. / (counts * blocks), width)[support]
        risks.append(scaled_errors**2 @ risk_weights)
        diagnostics.append(np.stack([a[support] @ risk_weights for a in (
            (bias / reference)**2, variance / reference**2, mse / reference**2,
        )]))
        arrays.update({f"{key}_{iteration}": value for key, value in dict(
            unfolded=ensemble, mean=mean, bias=bias, variance=variance, mse=mse,
            coverage=coverage, pulls=pulls, fb=fb, conditional_covariance=result.unfolded.covariance,
        ).items()})
    risks = np.stack(risks, axis=1)
    loss = risks.mean(axis=0)
    best = int(np.argmin(loss))
    # Use the same bootstrap draw at all iterations to preserve pairing.
    draws = rng.multinomial(count, np.full(count, 1. / count), size=1000)
    bootstrap_mean = draws @ risks / count
    winners = np.argmin(bootstrap_mean, axis=1)
    differences = risks - risks[:, best, None]
    difference_se = differences.std(axis=0, ddof=1) / np.sqrt(count)
    summary = dict(iterations=iterations, toys=count, seed=seed, exposure=exposure,
                   precision_variance_multiplier=precision, best_iteration=iterations[best],
                   global_mse=loss.tolist(), paired_difference_se=difference_se.tolist(),
                   bootstrap_winner_fraction=(np.bincount(winners, minlength=len(iterations))
                                              / len(winners)).tolist(),
                   score="equal pT block means of squared errors divided by squared block mean truth",
                   coverage="fixed linearized covariance at expected reco; nominal 1 sigma vs truth",
                   provenance=provenance or {})
    (output / "summary.json").write_text(json.dumps(summary, indent=2))
    np.savez_compressed(output / "toys.npz", truth=truth, expected=expected, weights=weights,
                        measured_toys=toys, support=support, pt_intervals=t.layout.pt_intervals,
                        eta_edges=eta_axis, fb_truth=fb_truth, iterations=iterations,
                        fb_measured_toys=_fb_array(toys, blocks, width),
                        diagnostics=np.asarray(diagnostics), toy_risk=risks, **arrays)
    if draw:
        _plot_scan(output, iterations, summary, np.asarray(diagnostics), arrays,
                   truth.reshape(blocks, width), fb_truth, eta_axis, t.layout.pt_intervals,
                   support.reshape(blocks, width))
    return summary


def redraw_regularization(directory):
    """Render figures from the archived ensemble without running unfolding."""
    summary = json.loads((directory / "summary.json").read_text())
    iterations = summary["iterations"]
    selected = (iterations[0], summary["best_iteration"], iterations[-1])
    with np.load(directory / "toys.npz") as data:
        intervals, eta = data["pt_intervals"], data["eta_edges"]
        shape = (len(intervals), len(eta)-1)
        arrays = {f"{field}_{k}": data[f"{field}_{k}"]
                  for field in ("mean", "bias", "variance", "mse") for k in iterations}
        arrays.update({f"{field}_{k}": data[f"{field}_{k}"]
                       for field in ("fb", "pulls", "coverage") for k in set(selected)})
        _plot_scan(directory, iterations, summary, data["diagnostics"], arrays,
                   data["truth"].reshape(shape), data["fb_truth"], eta, intervals,
                   data["support"].reshape(shape))


def _plot_scan(output, iterations, summary, diagnostics, arrays, truth, fb_truth,
               eta, intervals, support):
    plotter = plots(output / "mplhep")
    best = summary["best_iteration"]
    with plotter.plt.style.context(plotter.hep.style.CMS):

        def canvas():
            fig, ax = plotter.plt.subplots(figsize=(8, 6))
            fig.subplots_adjust(top=.82, bottom=.16, left=.19)
            plotter._header(fig, left=.19, top=.91, small=True)
            return fig, ax

        def iteration_label(k):
            return "1 iteration" if k == 1 else f"{k} iterations"

        def interval_title(ax, title):
            ax.set_title(title, fontsize=16)

        fig, ax = canvas()
        for row, label in enumerate(("Bias squared", "Variance", "MSE")):
            ax.plot(iterations, diagnostics[:, row], "o-", label=label)
        ax.axvline(best, color="0.4", ls="--")
        ax.set(xlabel="Bayesian iterations", ylabel="Global normalized risk")
        ax.set_xticks(iterations)
        ax.legend(fontsize=12)
        plotter._export(fig, "global_mse")
        center = (eta[:-1] + eta[1:]) / 2
        width = len(center)
        blocks = len(intervals)
        for i, interval in enumerate(intervals):
            block = slice(i*width, (i+1)*width)
            pt_label = rf"${interval[0]:g}\leq p_{{\rm T}}^{{\rm ave}}<{interval[1]:g}$ GeV"
            normalizer = truth[i].mean()
            mask = support[i]
            stem = f"pt{interval[0]:g}_{interval[1]:g}"
            fig, ax = canvas()
            for key, label in (("bias", "Bias squared"), ("variance", "Variance"), ("mse", "MSE")):
                rows = [arrays[f"{key}_{k}"][block][mask] for k in iterations]
                a = np.asarray(rows) / normalizer**(1 if key == "bias" else 2)
                ax.plot(iterations, (a**2 if key == "bias" else a).mean(axis=1), "o-", label=label)
            ax.axvline(best, color="0.4", ls="--")
            ax.set(xlabel="Bayesian iterations", ylabel="Normalized risk")
            ax.set_xticks(iterations)
            interval_title(ax, pt_label)
            ax.legend(fontsize=12)
            plotter._export(fig, f"{stem}_mse")
            fig, ax = canvas()
            bias = np.array([arrays[f"bias_{k}"][block] / normalizer for k in iterations])
            mesh = ax.pcolormesh(eta, np.arange(len(iterations)+1)-.5,
                                 np.where(mask, bias, np.nan), cmap="coolwarm")
            fig.colorbar(mesh, ax=ax, label="Bias / block mean truth")
            ax.set(xlabel=r"Dijet $\eta_{\mathrm{CM}}$", ylabel="Bayesian iterations")
            ax.set_yticks(np.arange(len(iterations)), iterations)
            interval_title(ax, pt_label)
            plotter._export(fig, f"{stem}_bias")
            fig, ax = canvas()
            ax.stairs(truth[i], eta, label="Input truth")
            mean = arrays[f"mean_{best}"][block]
            sigma = np.sqrt(arrays[f"variance_{best}"][block])
            ax.errorbar(center[mask], mean[mask], yerr=sigma[mask], xerr=(np.diff(eta)/2)[mask],
                        fmt="o", label=f"Toy mean and SD, {iteration_label(best)}")
            ax.set(xlabel=r"Dijet $\eta_{\mathrm{CM}}$", ylabel="Yield / bin")
            interval_title(ax, pt_label)
            ax.legend(fontsize=12)
            plotter._export(fig, f"{stem}_spectra")
            fig, ax = canvas()
            ax.stairs(fb_truth[i], eta[width//2:], baseline=None, label="Input truth")
            for iteration in dict.fromkeys((iterations[0], best, iterations[-1])):
                ensemble = arrays[f"fb_{iteration}"][:, i]
                valid = np.isfinite(ensemble).any(axis=0)
                low, median, high = np.nanquantile(ensemble[:, valid], [.16, .5, .84], axis=0)
                ax.errorbar(center[width//2:][valid], median, yerr=[median-low, high-median],
                            xerr=(np.diff(eta)[width//2:]/2)[valid],
                            fmt="o", label=iteration_label(iteration))
            ax.set(xlabel=r"Dijet $|\eta_{\mathrm{CM}}|$",
                   ylabel=r"Forward ($\eta_{\mathrm{CM}}>0$) / Backward ($\eta_{\mathrm{CM}}<0$)")
            ax.yaxis.label.set_size(14)
            interval_title(ax, pt_label)
            ax.legend(fontsize=11)
            plotter._export(fig, f"{stem}_fb")

        fig, ax = canvas()
        for i in dict.fromkeys((iterations[0], best, iterations[-1])):
            pulls = arrays[f"pulls_{i}"].ravel()
            ax.hist(pulls[np.isfinite(pulls)], bins=100, density=True, histtype="step",
                    label=iteration_label(i))
        ax.set(xlabel="Pull relative to input truth", ylabel="Density")
        ax.legend(fontsize=12)
        plotter._export(fig, "pulls")
        fig, ax = canvas()
        for i in dict.fromkeys((iterations[0], best, iterations[-1])):
            coverage = arrays[f"coverage_{i}"].reshape(blocks, width)
            ax.plot(np.arange(1, blocks+1), np.nanmean(coverage, axis=1), "o-",
                    label=iteration_label(i))
        ax.axhline(.682689, color=".5", ls="--")
        ax.set(xlabel="pTave interval index", ylabel="Fraction within one linearized SD")
        ax.legend(fontsize=12)
        plotter._export(fig, "coverage")


def save_stitched_input(training, measured, overlap, mb_yield, provenance, output, *, draw=True):
    """Archive measured yields and the corresponding CM forward/backward ratios."""
    from hist_analysis.python.unfolding_validation import forward_backward

    spectrum = _fixed(measured)
    blocks = dict(_blocks(training, spectrum))
    fb = {interval: forward_backward(block, name=f"fb_stitched_reco_pt{i}")
          for i, (interval, block) in enumerate(blocks.items())}
    exposure = mb_yield / values(training.reco)[overlap].sum()
    mc = {
        label: _fixed(histogram(h, name, exposure * values(h), exposure**2 * variances(h)))
        for label, name, h in (("Embedded Reco", "embedding_reco_mb", training.reco),
                               ("Embedded Gen", "embedding_gen_mb", training.gen))
    }
    provenance = dict(provenance, embedded_exposure_scale=exposure,
                      embedded_normalization="one scale from Embedded Reco in direct MB overlap")
    write_unfolding_output(output / "stitched_input.root",
                           histograms=[measured, *(b.histogram for b in blocks.values()),
                                       *(b.histogram for b in fb.values()),
                                       *(s.histogram for s in mc.values())],
                           metadata=provenance)
    (output / "stitching.json").write_text(json.dumps(provenance, indent=2))
    if draw:
        plotter = plots(output / "mplhep", data=True)
        for kind, spectra in (("eta", blocks), ("fb", fb)):
            comparisons(plotter, kind, None, {i: {"Reco": s} for i, s in spectra.items()},
                        prefix="data_input")
        for label, s in mc.items():
            mc[label] = dict(_blocks(training, s))
        comparison = {interval: {"Reco": b, **{label: slices[interval]
                                               for label, slices in mc.items()}}
                      for interval, b in blocks.items()}
        comparisons(plotter, "eta", None, comparison, prefix="data_embedded_mb_exposure")
