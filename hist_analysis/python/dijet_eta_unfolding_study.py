"""Response diagnostics, closure and regularization for eta-dependent dijets."""

import json

import numpy as np
import ROOT

from hist_analysis.python.dijet_eta_unfolding import (
    extended_response, geometric_response, iterate, load_eta_training,
    same_event_fiducial,
)
from hist_analysis.python.dijet_unfolding_workflow import (
    ANNOTATIONS, _blocks, _fb_array, poisson_toys, risk_metrics,
)
from hist_analysis.python.histogram_io import open_root_file
from hist_analysis.python.plot_labels import (
    CMS_PRELIMINARY_HEADER, CMS_SIMULATION_HEADER, PPB_MPLHEP_HEADER_RIGHT,
)
from hist_analysis.python.unfolding import write_unfolding_output
from hist_analysis.python.unfolding_validation import (
    Spectrum, _bayes_derivative, agreement, forward_backward, histogram, ratio,
    root_matrix, run_unfolding, values, variances,
)
from hist_analysis.python.unfolding_validation_plots import ValidationPlots


def inputs(filename=None, *, sample="full", intervals=None):
    options = {} if intervals is None else {"intervals": intervals}
    training = load_eta_training(filename, sample=sample, **options)
    extended = extended_response(training)
    return training, {
        "A": geometric_response(training),
        "bridge": same_event_fiducial(extended),
        "B": extended,
    }


def plotter(output, *, data=False, formats=("pdf",)):
    return ValidationPlots(output / "mplhep", formats,
                           left=CMS_PRELIMINARY_HEADER if data else CMS_SIMULATION_HEADER,
                           right=PPB_MPLHEP_HEADER_RIGHT,
                           annotations=(*ANNOTATIONS, "Eta-dependent hybrid JER"), figsize=(6, 6))


def response_diagnostics(training, responses, output, *, draw=False, formats=("pdf",)):
    summary = {}
    f = open_root_file(str(output / "eta_responses.root"), "RECREATE")
    try:
        for name, response in responses.items():
            directory = f.mkdir(name)
            directory.cd()
            kernel = response.kernel
            singular = np.linalg.svd(kernel, compute_uv=False)
            threshold = singular[0] * max(kernel.shape) * np.finfo(float).eps
            rank = int(np.count_nonzero(singular > threshold))
            summary[name] = dict(shape=kernel.shape, rank=rank, nullity=kernel.shape[1] - rank,
                                 condition=float(singular[0]/singular[rank - 1]) if rank else None,
                                 zero_efficiency=int(np.count_nonzero(response.efficiency == 0)),
                                 background_components=response.labels)
            for key, matrix in (
                ("migration", response.migration), ("migration_sumw2", response.migration_sumw2),
                ("kernel", kernel), ("background_templates", response.components()),
                ("truth_and_sumw2", np.stack((response.truth, response.truth_sumw2))),
                ("measured_and_sumw2", np.stack((response.measured, response.measured_sumw2))),
                ("efficiency", response.efficiency[None, :]), ("singular_values", singular[None, :]),
            ):
                root_matrix(matrix).Write(key)
        f.cd()
        ROOT.TNamed("provenance", json.dumps(dict(training=training.metadata, methods=summary))).Write()
    finally:
        f.Close()
    (output / "response_summary.json").write_text(json.dumps(summary, indent=2))
    if draw:
        charts = plotter(output, formats=formats)
        charts.response(training)
        response = responses["B"]
        with charts.plt.style.context(charts.hep.style.CMS):
            fig, ax = charts.plt.subplots(figsize=(6, 6))
            fig.subplots_adjust(left=.17, right=.84, top=.84, bottom=.15)
            charts._header(fig, left=.17, top=.91, small=True)
            n_reco, n_truth = response.migration.shape
            image = ax.imshow(np.ma.masked_equal(response.kernel.T, 0), origin="lower", aspect="auto",
                              extent=(.5, n_reco+.5, .5, n_truth+.5), interpolation="nearest")
            colorbar = fig.colorbar(image, ax=ax)
            colorbar.set_label("Migration probability", fontsize=13)
            colorbar.ax.tick_params(labelsize=11)
            ax.set_xlabel("Flattened reco bin", fontsize=15)
            ax.set_ylabel("Flattened extended truth bin", fontsize=15)
            ax.tick_params(labelsize=12)
            ax.axhline(response.n_fiducial + .5, color="white", linewidth=1.)
            charts._export(fig, "method_b_extended_response")

            fig, ax = charts.plt.subplots(figsize=(6, 6))
            fig.subplots_adjust(left=.31, right=.83, top=.85, bottom=.20)
            charts._header(fig, left=.17, top=.91, small=True)
            components = response.components().reshape(training.layout.n_pt_bins,
                                                       training.layout.n_eta_bins, -1).sum(1)
            denominators = response.measured.reshape(training.layout.n_pt_bins, -1).sum(1)
            fraction = np.divide(components, denominators[:, None],
                                 out=np.zeros_like(components), where=denominators[:, None] > 0)
            x = np.arange(training.layout.n_pt_bins)
            ax.imshow(np.ma.masked_equal(fraction.T, 0), origin="upper", aspect="auto",
                      interpolation="nearest")
            ax.set_xticks(x, [f"{lo:g}–{hi:g}" for lo, hi in training.layout.pt_intervals],
                          rotation=65, ha="right", fontsize=8)
            ax.set_yticks(np.arange(len(response.labels)), (
                "Passes cuts, below range", "Passes cuts, other flow", "Fails cuts, below range",
                "Fails cuts, in range", "Fails cuts, other flow", "No gen pair",
            ), fontsize=8)
            ax.set_xlabel(r"Reco $p_{\mathrm{T}}^{\mathrm{ave}}$ (GeV)", fontsize=14)
            colorbar = fig.colorbar(ax.images[0], ax=ax, fraction=.05)
            colorbar.set_label("Fraction of selected reco", fontsize=12)
            colorbar.ax.tick_params(labelsize=10)
            charts._export(fig, "method_b_sideband_fractions")
    return summary


def estimate(training, response, measured, *, iterations, exposure=1., rates=None,
             rate_covariance=None, start=None, measured_covariance=None, factorized=False):
    """Estimate fiducial yields and propagate conditional data and external-rate covariance."""
    if not np.isfinite(exposure) or exposure <= 0:
        raise ValueError("Training-to-measurement exposure must be positive")
    data = values(measured)
    vy = np.diag(variances(measured)) if measured_covariance is None else measured_covariance
    if factorized:
        if any(x is not None for x in (rates, rate_covariance, start)):
            raise ValueError("Use the explicit model for nuisance or prior variations")
        r = run_unfolding(training, measured, iterations=iterations, measured_covariance=vy)
        p = np.nan_to_num(values(training.purity))
        inv = np.nan_to_num(values(training.inverse_efficiency))
        _, jac = _bayes_derivative(training.probabilities, p*data,
                                   values(training.matched_gen), iterations)
        return r.unfolded, r.refolded, inv[:, None]*jac*p, r.refolded_cross
    n = response.n_fiducial
    kernel = response.kernel[:, :n]
    background = exposure * response.background(rates)
    templates = exposure * response.components()
    prior = exposure * response.truth[:n] if start is None else start
    x, jac, rate_jac = iterate(kernel, data, background, prior, iterations,
                               derivative=True, background_templates=templates)
    fold_jac = kernel @ jac
    vx = jac @ vy @ jac.T
    vf = fold_jac @ vy @ fold_jac.T
    if rate_covariance is not None:
        covariance = np.asarray(rate_covariance, dtype=float)
        if (covariance.shape != (templates.shape[1],)*2 or not np.isfinite(covariance).all()
                or not np.allclose(covariance, covariance.T)
                or np.linalg.eigvalsh(covariance).min() < -1e-12):
            raise ValueError("Supply a positive semidefinite external covariance for component rates")
        vx += rate_jac @ covariance @ rate_jac.T
        fold_rate_jac = kernel @ rate_jac + templates
        vf += fold_rate_jac @ covariance @ fold_rate_jac.T
    refold = kernel @ x + background
    return (
        Spectrum(histogram(training.gen, "unfolded", x, vx.diagonal()), vx, kernel.sum(0) > 0),
        Spectrum(histogram(training.reco, "refolded", refold, vf.diagonal()),
                 vf, (kernel.sum(1) > 0) | (background > 0)),
        jac, fold_jac @ vy,
    )


def _fb_jacobian(spectrum):
    y = values(spectrum.histogram)
    half = len(y)//2
    neg = np.arange(half - 1, -1, -1)
    reciprocal = np.divide(1., y[neg], out=np.zeros(half), where=y[neg] > 0)
    j = np.zeros((half, len(y)))
    j[np.arange(half), np.arange(half, len(y))] = reciprocal
    j[np.arange(half), neg] = -y[half:]*reciprocal**2
    return j


def comparisons(training, responses, measured, output, *, iterations, exposure=1.,
                truth=None, reco_truth_cross=None, rate_covariances=None, draw=False,
                data=False, formats=("pdf",), provenance=None):
    raw = values(measured)
    vy = np.diag(variances(measured))
    reco = Spectrum(measured, vy, np.isfinite(raw))
    gen = None if truth is None else Spectrum(
        truth, np.diag(variances(truth)), np.isfinite(values(truth))
    )
    results = {}
    for name, model in responses.items():
        k = iterations[name] if isinstance(iterations, dict) else iterations
        unfolded, refolded, jac, refold_cross = estimate(
            training, model, measured, iterations=k, exposure=exposure,
            rate_covariance=(rate_covariances or {}).get(name), factorized=name == "A_factorized",
        )
        unfolded_reco_cross = jac @ vy
        unfolded_gen_cross = None if reco_truth_cross is None else jac @ reco_truth_cross
        methods = {"Reco": reco, "Unfolded": unfolded}
        spectra = {"eta": methods, "refolded": {"Reco": reco, "Refolded": refolded},
                   "unfolded_over_reco": {"Unfolded / Reco": ratio(
                       unfolded, reco, name="unfolded_over_reco", cross=unfolded_reco_cross)},
                   "refolded_ratio": {"Refolded / Reco": ratio(
                       refolded, reco, name="refolded_over_reco", cross=refold_cross)}}
        if gen is not None:
            methods["Gen"] = gen
            spectra["to_gen"] = {
                "Reco / Gen": ratio(reco, gen, name="reco_over_gen", cross=reco_truth_cross),
                "Unfolded / Gen": ratio(unfolded, gen, name="unfolded_over_gen",
                                       cross=unfolded_gen_cross),
            }
        sliced = {kind: {p: {} for p in training.layout.pt_intervals} for kind in spectra}
        for kind, overlays in spectra.items():
            for label, item in overlays.items():
                for pt, block in _blocks(training, item):
                    sliced[kind][pt][label] = block
        sliced["fb"] = {}
        sliced["fb_double"] = {}
        fb_transforms = {label: np.zeros((len(raw)//2, len(raw))) for label in methods}
        for i, (pt, block) in enumerate(sliced["eta"].items()):
            fb = {label: forward_backward(item, name=f"fb_{label}_pt{i}")
                  for label, item in block.items()}
            sliced["fb"][pt] = fb
            width = training.layout.n_eta_bins
            local = slice(i*width, (i+1)*width)
            folded = slice(i*width//2, (i+1)*width//2)
            for label, item in block.items():
                fb_transforms[label][folded, local] = _fb_jacobian(item) * fb[label].supported[:, None]
            cross = (_fb_jacobian(block["Unfolded"]) @ unfolded_reco_cross[local, local]
                     @ _fb_jacobian(block["Reco"]).T)
            sliced["fb_double"][pt] = {"Unfolded / Reco": ratio(
                fb["Unfolded"], fb["Reco"], cross=cross, name=f"fb_double_pt{i}")}
            if gen is not None:
                cross_rg = np.zeros_like(vy) if reco_truth_cross is None else reco_truth_cross
                cross_ug = np.zeros_like(vy) if unfolded_gen_cross is None else unfolded_gen_cross
                for label, covariance in (("Reco", cross_rg), ("Unfolded", cross_ug)):
                    cross = _fb_jacobian(block[label]) @ covariance[local, local] @ _fb_jacobian(block["Gen"]).T
                    sliced["fb_double"][pt][f"{label} / Gen"] = ratio(
                        fb[label], fb["Gen"], cross=cross, name=f"fb_{label}_to_gen_pt{i}")
        checks = {"refolded": agreement(refolded, reco, refold_cross)}
        if gen is not None:
            cross = np.zeros_like(vy) if unfolded_gen_cross is None else unfolded_gen_cross
            checks["closure"] = agreement(unfolded, gen, cross)
            gen_yield = values(gen.histogram)
            checks["closure"]["truth_fraction_without_response"] = (
                float(gen_yield[~unfolded.supported].sum()/gen_yield.sum())
                if gen_yield.sum() else None
            )
        archive = output / f"{name}_iterations_{k}.root"
        quantities = {item.histogram.GetName(): item
                      for collection in (*spectra.values(), *(g for pts in sliced.values() for g in pts.values()))
                      for item in collection.values()}
        write_unfolding_output(archive, histograms=[item.histogram for item in quantities.values()],
                               metadata=dict(training=training.metadata, method=name, iterations=k,
                                             exposure=exposure, provenance=provenance or {},
                                             uncertainty="conditional data and supplied external rates",
                                             unsupported={key: np.flatnonzero(~item.supported).tolist()
                                                          for key, item in quantities.items()}))
        with open_root_file(str(archive), "UPDATE") as f:
            f.cd()
            for key, item in quantities.items():
                root_matrix(item.covariance).Write(f"{key}_covariance")
            root_matrix(unfolded_reco_cross).Write("unfolded_reco_cross")
            root_matrix(refold_cross).Write("refolded_reco_cross")
            if unfolded_gen_cross is not None:
                root_matrix(unfolded_gen_cross).Write("unfolded_gen_cross")
            for label, transform in fb_transforms.items():
                root_matrix(transform @ methods[label].covariance @ transform.T).Write(
                    f"fb_{label.lower()}_full_covariance"
                )
        if draw:
            charts = plotter(output, data=data, formats=formats)
            for kind, pts in sliced.items():
                charts.comparisons(kind, spectra.get(kind), pts, prefix=name)
            with charts.plt.style.context(charts.hep.style.CMS):
                variance = np.sqrt(np.maximum(unfolded.covariance.diagonal(), 0.))
                product = variance[:, None]*variance
                correlation = np.divide(unfolded.covariance, product,
                                        out=np.full_like(unfolded.covariance, np.nan), where=product > 0)
                for label, matrix in (("covariance", unfolded.covariance), ("correlation", correlation)):
                    fig, ax = charts.plt.subplots(figsize=(6, 6))
                    fig.subplots_adjust(left=.17, right=.84, bottom=.15, top=.84)
                    charts._header(fig, left=.17, top=.91, small=True)
                    limits = {"vmin": -1, "vmax": 1} if label == "correlation" else {}
                    visible = np.ma.masked_equal(np.ma.masked_invalid(matrix), 0)
                    im = ax.imshow(visible, origin="lower", cmap="coolwarm", **limits)
                    ax.set(xlabel="Flattened truth bin", ylabel="Flattened truth bin")
                    fig.colorbar(im, ax=ax, label=label.capitalize())
                    charts._export(fig, f"{name}_unfolded_{label}")
        results[name] = dict(spectrum=unfolded, jacobian=jac, checks=checks)
    if {"A", "B"} <= results.keys():
        a, b = results["A"], results["B"]
        cross = b["jacobian"] @ vy @ a["jacobian"].T
        change = ratio(b["spectrum"], a["spectrum"], name="method_b_over_a", cross=cross)
        archive = output / "method_b_over_a.root"
        write_unfolding_output(archive, histograms=(change.histogram,),
                               metadata=dict(training=training.metadata,
                                             cross_covariance="shared measured data only"))
        with open_root_file(str(archive), "UPDATE") as f:
            f.cd()
            root_matrix(change.covariance).Write("method_b_over_a_covariance")
        if draw:
            charts = plotter(output, data=data, formats=formats)
            charts.comparisons("method_comparison", {"B / A": change},
                               {pt: {"B / A": block} for pt, block in _blocks(training, change)},
                               prefix="method_b_over_a")
    (output / "closure_summary.json").write_text(json.dumps({
        method: {test: {key: value for key, value in info.items() if np.ndim(value) == 0}
                 for test, info in result["checks"].items()}
        for method, result in results.items()
    }, indent=2))
    return results


def regularization(training, responses, target, data_weights, anchor, anchor_yield, output,
                   *, iterations=12, count=1000, seed=1, sideband_rates=None,
                   draw=False, formats=("pdf",), store_trajectories=False, provenance=None):
    """Paired conditional toys from a common same-event generating model."""
    b = responses["B"]
    n = b.n_fiducial
    target = np.asarray(target, dtype=float)
    if (target.shape != (n,) or not np.isfinite(target).all() or np.any(target < 0)
            or not isinstance(iterations, int) or iterations < 1):
        raise ValueError("Use nonnegative finite fiducial truth and a positive scan length")
    background = b.background(sideband_rates)
    unscaled = b.fold(target, background=background)
    if anchor_yield <= 0 or unscaled[anchor].sum() <= 0:
        raise ValueError("The common MB anchor needs positive observed and predicted counts")
    exposure = anchor_yield / unscaled[anchor].sum()
    truth = target*exposure
    expected = unscaled*exposure
    rng = np.random.default_rng(seed)
    toys = poisson_toys(expected, data_weights, count, rng)
    width = training.layout.n_eta_bins
    block_scale = np.repeat((exposure*b.truth[:n]).reshape(-1, width).mean(1), width)
    common = np.logical_and.reduce([model.efficiency[:n] > 0 for model in responses.values()])
    per_block = common.reshape(-1, width).sum(1)
    if np.any(per_block == 0) or np.any(block_scale[common] <= 0) or np.any(truth[~common] > 0):
        raise ValueError("The toy truth or scoring blocks extend outside common response support")
    weights = np.repeat(1/(len(per_block)*per_block), width)[common]
    arrays, summary = dict(expected=expected, truth=truth, expected_variance=data_weights*expected), {}
    if store_trajectories:
        arrays["measured_toys"] = toys
    truth_fb = _fb_array(truth, len(per_block), width)
    for name, model in responses.items():
        kernel = model.kernel[:, :n]
        estimate_values = np.broadcast_to(exposure*model.truth[:n], (count, n))
        inverse = np.divide(1., kernel.sum(0), out=np.zeros(n), where=kernel.sum(0) > 0)
        baseline = exposure*model.background()
        if np.any((kernel.sum(1) == 0) & (baseline == 0) & (expected > 0)):
            raise ValueError(f"{name}: pseudodata outside model support")
        risks = []
        for k in range(1, iterations + 1):
            mean = estimate_values @ kernel.T + baseline
            transport = np.divide(toys, mean, out=np.zeros_like(mean), where=mean > 0)
            estimate_values = estimate_values * (transport @ kernel) * inverse
            ensemble_mean, bias, variance, mse = risk_metrics(estimate_values, truth)
            _, dy = iterate(kernel, expected, baseline, exposure*model.truth[:n], k, derivative=True)
            sigma = np.sqrt((dy**2) @ (data_weights*expected))
            risk = ((estimate_values[:, common] - truth[common])/block_scale[common])**2 @ weights
            risks.append(risk)
            fb = _fb_array(estimate_values, len(per_block), width)
            fb_available = np.isfinite(fb) & np.isfinite(truth_fb)
            fb_count = fb_available.sum(axis=0)
            fb_squared_error = np.where(fb_available, (fb - truth_fb)**2, 0.).sum(axis=0)
            fb_mse = np.divide(fb_squared_error, fb_count,
                               out=np.full_like(truth_fb, np.nan), where=fb_count > 0)
            for key, item in dict(mean=ensemble_mean, bias=bias, variance=variance, mse=mse,
                                  sigma=sigma, fb_mse=fb_mse, fb_valid_fraction=fb_count/count).items():
                arrays[f"{name}_{k}_{key}"] = item
            if store_trajectories:
                arrays[f"{name}_{k}_unfolded"] = estimate_values
                arrays[f"{name}_{k}_fb"] = fb
            for level, z in ((68, 1.), (95, 1.959964)):
                coverage = np.mean(abs(estimate_values-truth) <= z*sigma, axis=0)
                arrays[f"{name}_{k}_coverage{level}"] = np.where(common & (sigma > 0), coverage, np.nan)
        risk = np.asarray(risks).T
        loss = risk.mean(0)
        best = int(np.argmin(loss))
        se = (risk - risk[:, best, None]).std(0, ddof=1)/np.sqrt(count)
        summary[name] = dict(loss=loss.tolist(), paired_se=se.tolist(), minimum=best+1,
                             candidate=int(np.flatnonzero(loss <= loss[best]+se)[0])+1)
        candidate = summary[name]["candidate"]
        candidate_coverage = {}
        for level in (68, 95):
            coverage = arrays[f"{name}_{candidate}_coverage{level}"]
            finite = np.isfinite(coverage)
            candidate_coverage[str(level)] = float(coverage[finite].mean()) if finite.any() else None
        summary[name]["mean_coverage_at_candidate"] = candidate_coverage
        arrays[f"{name}_risk"] = risk
    np.savez_compressed(output / "paired_toys.npz", **arrays)
    (output / "regularization_summary.json").write_text(json.dumps(
        dict(seed=seed, count=count, iterations=iterations, exposure=exposure,
             response_training=training.metadata,
             provenance=provenance or {},
             generator_sideband_rates=None if sideband_rates is None else list(sideband_rates),
             uncertainty="measurement only at fixed response",
             methods=summary), indent=2))
    if draw:
        charts = plotter(output, formats=formats)
        with charts.plt.style.context(charts.hep.style.CMS):
            for name in responses:
                fig, ax = charts.plt.subplots(figsize=(6, 6))
                fig.subplots_adjust(left=.19, right=.94, top=.84, bottom=.14)
                charts._header(fig, left=.19, top=.91, small=True)
                x = np.arange(1, iterations + 1)
                for label, metric in (("Squared bias", "bias"), ("Variance", "variance"), ("MSE", "mse")):
                    y = np.asarray([arrays[f"{name}_{k}_{metric}"][common]/block_scale[common]**(1 if metric == "bias" else 2)
                                    for k in x])
                    ax.plot(x, (y*y if metric == "bias" else y) @ weights,
                            marker="o", markersize=4, label=label)
                ax.set(xlabel="D'Agostini iterations", ylabel="Mean scaled error squared")
                ax.legend(frameon=False)
                charts._export(fig, f"method_{name}_risk")
    return summary


def bottom_line(model, hypotheses, weights, output, *, iterations, count=1000, seed=1,
                exposure=1., observed=None, provenance=None):
    """Calibrate conditional reco and unfolded hypothesis tests on paired toys."""
    from scipy.special import xlogy

    target = np.asarray(hypotheses, dtype=float)
    n = model.n_fiducial
    if (target.shape != (2, n) or not np.isfinite(target).all() or np.any(target < 0)
            or exposure <= 0 or not np.isfinite(exposure)):
        raise ValueError("Supply two finite nonnegative fiducial predictions at the same exposure")
    weights = np.asarray(weights, dtype=float)
    if (weights.shape != (len(model.measured),) or np.any(weights <= 0)
            or not np.isfinite(weights).all()):
        raise ValueError("Supply positive stream weights for every reco bin")
    kernel = model.kernel[:, :n]
    prior = exposure*model.truth[:n]
    background = exposure*model.background()
    target = target*exposure
    predicted = target @ kernel.T + background
    rng = np.random.default_rng(seed)
    # Center estimation, threshold calibration and test toys are disjoint.
    ensembles = np.stack([poisson_toys(mu, weights, 3*count, rng) for mu in predicted])
    observed = None if observed is None else np.asarray(observed, dtype=float)
    if observed is not None and (observed.shape != weights.shape or np.any(observed < 0)
                                 or not np.isfinite(observed).all()):
        raise ValueError("Observed spectrum must contain finite nonnegative measured counts")
    rows = ensembles.reshape(-1, len(weights))
    if observed is not None:
        rows = np.vstack((rows, observed))
    unfold = np.stack([iterate(kernel, y, background, prior, iterations) for y in rows])
    reco_scores, unfolded_scores, calibrated_scores, ranks = [], [], [], []
    for i, (truth, mu) in enumerate(zip(target, predicted)):
        expected_n = mu / weights
        observed_n = rows / weights
        # xlogy implements the empty-count limit and preserves infinity for impossible counts.
        deviance = 2*np.sum(xlogy(observed_n, observed_n) - xlogy(observed_n, expected_n)
                           - observed_n + expected_n, axis=1)
        reco_scores.append(deviance)
        _, dy = iterate(kernel, mu, background, prior, iterations, derivative=True)
        center = unfold[3*i*count:(3*i+1)*count].mean(axis=0)
        covariance = (dy*(mu*weights)) @ dy.T
        eigenvalues, modes = np.linalg.eigh(covariance)
        keep = eigenvalues > np.max(eigenvalues)*max(covariance.shape)*np.finfo(float).eps
        ranks.append(int(keep.sum()))
        whiten = modes[:, keep]/np.sqrt(eigenvalues[keep])
        unfolded_scores.append(np.sum(((unfold - truth) @ whiten)**2, axis=1))
        calibrated_scores.append(np.sum(((unfold - center) @ whiten)**2, axis=1))
    statistics = dict(reco=np.stack(reco_scores), unfolded_truth=np.stack(unfolded_scores),
                      unfolded_expected=np.stack(calibrated_scores))
    summary = dict(iterations=iterations, count=count, seed=seed, unfolded_covariance_ranks=ranks,
                   exposure=exposure, covariance="conditional at each hypothesis", statistics={},
                   provenance=provenance or {})
    for name, score in statistics.items():
        per_test = []
        for i in range(2):
            start = 3*i*count
            other_start = 3*(1-i)*count
            calibration = score[i, start+count:start+2*count]
            own = score[i, start+2*count:start+3*count]
            other = score[i, other_start+2*count:other_start+3*count]
            cut = np.quantile(calibration, .95, method="higher")
            result = dict(null_rejection_at_empirical_95=float(np.mean(own > cut)),
                          alternative_rejection=float(np.mean(other > cut)))
            if observed is not None:
                result["observed_pvalue"] = float((1 + np.count_nonzero(calibration >= score[i, -1]))/(count+1))
            per_test.append(result)
        summary["statistics"][name] = per_test
    np.savez_compressed(output / "bottom_line.npz", reco_predictions=predicted,
                        truth_hypotheses=target, **statistics)
    (output / "bottom_line.json").write_text(json.dumps(summary, indent=2))
    return summary
