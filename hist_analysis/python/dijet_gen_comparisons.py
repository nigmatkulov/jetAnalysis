"""Compare unfolding stages to Gen and audit sub-threshold response migration."""

from __future__ import annotations

from array import array
import argparse
import json
from pathlib import Path

import numpy as np

from hist_analysis.python.dijet_unfolding_workflow import (
    CUT, ETA_RANGE, _blocks, _fixed, comparisons, plots, setup, training_input,
)
from hist_analysis.python.histogram_io import load_object, open_root_file
from hist_analysis.python.unfolding_validation import (
    Spectrum, edges, histogram, load_training, make_training, numpy_matrix, ratio, run_unfolding,
    values, variances,
)


def _spectrum(template, name, content, covariance, support=None):
    support = np.isfinite(content) if support is None else support
    return Spectrum(histogram(template, name, content, covariance.diagonal()),
                    covariance, support)


def _fold_operator(training):
    efficiency = values(training.inverse_efficiency)
    purity = values(training.purity)
    inverse_purity = np.divide(1., purity, out=np.zeros_like(purity), where=purity > 0)
    forward_efficiency = np.divide(
        1., efficiency, out=np.zeros_like(efficiency), where=efficiency > 0,
    )
    return inverse_purity[:, None] * training.probabilities * forward_efficiency


def _archive_curves(training, output, reference, curves, crosses, *, data=False, notes=()):
    output.mkdir()
    plotter = plots(output / "mplhep", data=data)
    plotter.annotations += notes
    ratios = {f"{label} / Gen": ratio(s, reference, name=f"ratio_{i}", cross=crosses.get(label))
              for i, (label, s) in enumerate(curves.items())}
    overlay = curves | {"Gen": reference}
    for kind, group in (("eta", overlay), ("to_gen", ratios)):
        blocks = {interval: {} for interval in training.layout.pt_intervals}
        for label, s in group.items():
            for interval, block in _blocks(training, s):
                blocks[interval][label] = block
        comparisons(plotter, kind, None, blocks, prefix=output.name)
    archive = {"eta_edges": edges(training.templates[0].GetXaxis()),
               "pt_intervals": training.layout.pt_intervals}
    for group in (overlay, ratios):
        for label, s in group.items():
            archive[f"{label}::yield"] = values(s.histogram)
            archive[f"{label}::covariance"] = s.covariance
            archive[f"{label}::supported"] = s.supported
    np.savez_compressed(output / "curves.npz", **archive)


def mc_comparisons(training, saved, output, iterations):
    operator = _fold_operator(training)
    for iteration in iterations:
        filename = saved / f"iterations_{iteration}.root"
        root = open_root_file(str(filename))
        try:
            gen = _fixed(load_object(str(filename), "training_gen"))
            curves = {}
            for label, key, covariance in (
                ("Reco", "training_reco", "reco_covariance"),
                ("Unfolded", "corrected_unfolded", "unfolded_covariance"),
                ("Refolded", "refolded", "refolded_covariance"),
            ):
                h = load_object(str(filename), key)
                curves[label] = _spectrum(h, key, values(h), numpy_matrix(root.Get(covariance)))
        finally:
            root.Close()
        np.testing.assert_allclose(operator @ values(curves["Unfolded"].histogram),
                                   values(curves["Refolded"].histogram), rtol=1e-9, atol=1e-18)
        with np.load(saved / f"covariance_iterations_{iteration}.npz") as cov:
            unfolded_cross = cov["unfolded_gen_cross"]
            crosses = {"Reco": cov["reco_gen_cross"], "Unfolded": unfolded_cross,
                       "Refolded": operator @ unfolded_cross}
        _archive_curves(training, output / f"iterations_{iteration}", gen, curves, crosses,
                        notes=(f"Bayes iterations: {iteration}",))


def toy_comparisons(training, saved, output, iteration):
    with np.load(saved) as archive:
        expected, truth = archive["expected"], archive["truth"]
        ensemble = archive[f"unfolded_{iteration}"]
        folded = ensemble @ _fold_operator(training).T
        samples = {
            "Reco": (expected, expected * archive["weights"]),
            "Unfolded": (ensemble.mean(axis=0), ensemble.var(axis=0)),
            "Refolded": (folded.mean(axis=0), folded.var(axis=0)),
        }
    gen = _spectrum(training.gen, "toy_gen", truth, np.zeros((truth.size, truth.size)))
    curves = {label: _spectrum(training.gen, label, mean, np.diag(var))
              for label, (mean, var) in samples.items()}
    _archive_curves(training, output, gen, curves, {},
                    notes=(f"CP5 toys, {iteration} iteration", "Bars show per-bin toy spread"))


def data_comparisons(training, saved, output, iteration):
    reco, gen = (load_object(str(saved), key) for key in ("stitched_reco", "embedding_gen_mb"))
    result = run_unfolding(training, reco, iterations=iteration)
    _archive_curves(training, output, _fixed(gen),
                    {"Data": result.reco, "Unfolded data": result.unfolded,
                     "Refolded data": result.refolded}, {}, data=True,
                    notes=(f"Exploratory, {iteration} Bayes iteration", "Gen scaled in MB overlap"))


def _feed_in_file(training, filename, sample):
    """Count disjoint sub-threshold Gen categories in each selected Reco bin."""
    width = training.layout.n_eta_bins
    sources = ("fiducial_below", "outside_fiducial_below", "fiducial_in_analysis",
               "outside_fiducial_above", "fiducial_elsewhere")
    case_yields = np.zeros((len(sources), 2, training.layout.n_global_bins))
    joint = load_object(str(filename), f"hUnfoldingExtendedResponse{sample}_{CUT}")
    if joint.GetNdimensions() != 5 or not joint.GetCalculateErrors():
        raise ValueError("Unexpected same-event extended response")
    ref = load_object(str(filename), f"hGenDijetPtEtaCMVsRecoPtEtaCM{sample}_{CUT}")
    geometric_below = np.zeros((2, training.layout.n_global_bins))
    eligible_all = np.zeros_like(geometric_below)
    eta = edges(training.templates[0].GetXaxis())
    lower_pt = training.layout.pt_intervals[0][0]
    cutoff = joint.GetAxis(0).FindFixBin(lower_pt)
    if not np.isclose(joint.GetAxis(0).GetBinLowEdge(cutoff), lower_pt):
        raise ValueError("The analysis threshold must be a stored Gen bin edge")

    for sparse, extended in ((joint, True), (ref, False)):
        coords = array("i", [0] * sparse.GetNdimensions())
        axes = [sparse.GetAxis(i) for i in range(4)]
        for axis, expected in zip(axes, (joint.GetAxis(i) for i in range(4))):
            np.testing.assert_allclose(edges(axis), edges(expected), rtol=0, atol=1e-9)
        maps = []
        for pt, y in ((axes[0], axes[1]), (axes[2], axes[3])):
            pt_map = {pt.FindFixBin(x): i for i, (lo, hi) in enumerate(training.layout.pt_intervals)
                      for x in edges(pt)[:-1] + np.diff(edges(pt)) / 2 if lo <= x < hi}
            eta_map = {y.FindFixBin(x): i for i, x in enumerate((eta[:-1] + eta[1:]) / 2)}
            maps.append((pt_map, eta_map))
        gen_map, reco_map = maps
        for i in range(sparse.GetNbins()):
            weight = sparse.GetBinContent(i, coords)
            reco_pt, reco_eta = reco_map[0].get(coords[2]), reco_map[1].get(coords[3])
            if reco_pt is None or reco_eta is None:
                continue
            j = reco_pt * width + reco_eta
            below = coords[0] < cutoff
            moments = (weight, sparse.GetBinError2(i))
            if extended:
                if coords[4] not in (1, 2):
                    raise ValueError("Invalid truth category")
                fiducial = coords[4] == 2
                inside = coords[0] in gen_map[0] and coords[1] in gen_map[1]
                category = (0 if fiducial else 1) if below else (
                    (2 if inside else 4) if fiducial else 3)
                case_yields[category, :, j] += moments
            else:
                eligible_all[:, j] += moments
                if below:
                    geometric_below[:, j] += moments

    reco_h = training.reco
    total = np.stack((values(reco_h), variances(reco_h)))
    joint_reco = case_yields.sum(axis=0)
    no_pair = total - joint_reco
    if np.any(no_pair < -2e-7 * total):
        raise ValueError("Extended joint yields exceed the inclusive Reco")
    matched = np.stack((values(training.matched_reco), variances(training.matched_reco)))
    np.testing.assert_array_less(matched, eligible_all + 2e-7 * total)
    np.testing.assert_array_less(geometric_below, case_yields[0] + 2e-7 * total)
    return dict(zip(sources, case_yields)) | {
        "geometric_below": geometric_below, "geometric_domain_fake": np.maximum(total-matched, 0.),
        "no_gen_pair": np.maximum(no_pair, 0.), "all_reco": total,
    }


def _without_measured_guard(sample):
    """Keep the guard truth while making its Reco rows unobserved."""
    matrix = sample.matrix.Clone("truth_guard_matrix")
    matrix.SetDirectory(0)
    width = sample.layout.n_eta_bins
    for i in range(1, width+1):
        for j in range(1, matrix.GetNbinsY()+1):
            matrix.SetBinContent(i, j, 0.)
            matrix.SetBinError(i, j, 0.)
    y, variance = values(sample.reco), variances(sample.reco)
    y[:width] = variance[:width] = 0.
    reco = histogram(sample.reco, "reco_without_guard", y, variance)
    return make_training(matrix, reco, sample.gen, sample.layout, sample.templates, sample.metadata)


def _guard_closure(full, control, studied):
    """Test the stored [40,50) Gen bin with and without measured guard counts."""
    intervals = ((40., 50.), *full.layout.pt_intervals)
    extended = []
    for sample in (control, studied):
        identity = sample.metadata["identity"]
        extended.append(load_training(identity["source"], identity["roles"], intervals, ETA_RANGE,
                                      selection="jetId"))
    truth_only = [_without_measured_guard(sample) for sample in extended]
    summaries = {}
    for iteration in (1, 3, 4):
        outputs = {}
        for name, train, validation, offset in (
            ("nominal", control, studied, 0),
            ("guard_truth_only", *truth_only, full.layout.n_eta_bins),
            ("guard_both_axes", *extended, full.layout.n_eta_bins),
        ):
            result = run_unfolding(train, validation.reco, iterations=iteration)
            truth = values(validation.gen)[offset:]
            unfolded = values(result.unfolded.histogram)[offset:]
            error = (unfolded-truth).reshape(-1, full.layout.n_eta_bins)
            counts = truth.reshape(error.shape).sum(axis=1)
            outputs[name] = dict(l1=float(np.abs(error).sum()/counts.sum()),
                                 l1_by_pt=(np.abs(error).sum(axis=1)/counts).tolist(),
                                 total_ratio=float(unfolded.sum()/truth.sum()),
                                 max_first_bin_relative=float(
                                     np.max(np.abs(error[0]/truth[:full.layout.n_eta_bins]))))
        summaries[iteration] = outputs
    return summaries


def response_audit(full, controlled, studied, output):
    output.mkdir()
    exposure = {}
    arrays = {}
    width = full.layout.n_eta_bins
    for label, training, suffix in (("full", full, ""), ("studied", studied, "StudiedSample")):
        source = training.metadata["identity"]["source"]
        counts = _feed_in_file(training, source, suffix)
        arrays.update({f"{label}::{key}": value for key, value in counts.items()})
        total_pt = counts["all_reco"][0].reshape(-1, width).sum(axis=1)
        exposure[label] = {key: (v[0].reshape(-1, width).sum(axis=1) / total_pt).tolist()
                           for key, v in counts.items()}
    guard = _guard_closure(full, controlled, studied)
    (output / "summary.json").write_text(json.dumps(dict(
        pt_intervals=full.layout.pt_intervals, fractions_of_all_reco=exposure,
        controlled_to_studied_guard_closure=guard,
    ), indent=2))
    np.savez_compressed(output / "feed_in.npz", eta_edges=edges(full.templates[0].GetXaxis()),
                        pt_intervals=full.layout.pt_intervals, **arrays)
    plotter = plots(output / "mplhep")
    with plotter.plt.style.context(plotter.hep.style.CMS):
        fig, ax = plotter.plt.subplots(figsize=(8, 6))
        fig.subplots_adjust(top=.83, bottom=.18, left=.18)
        plotter._header(fig, left=.18, top=.93, small=True)
        fractions = exposure["full"]
        x = np.arange(1, len(full.layout.pt_intervals)+1)
        for key, label in (("geometric_below", "Matched fiducial Gen < 50 GeV"),
                           ("fiducial_below", "Same-event fiducial Gen < 50 GeV"),
                           ("outside_fiducial_below", "Same-event nonfiducial Gen < 50 GeV")):
            ax.plot(x, 100*np.asarray(fractions[key]), "o-", label=label)
        ax.set(xlabel=r"Reco $p_{\mathrm{T}}^{\mathrm{ave}}$ interval index (lowest = 50–60 GeV)",
               ylabel="Fraction of selected Reco (%)")
        ax.xaxis.label.set_size(15)
        ax.legend(fontsize=11)
        plotter._export(fig, "subthreshold_gen_feed_in")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("full-mc", "split-mc", "cp5-toys", "stitched-data"):
        parser.add_argument(f"--{arg}", type=Path, required=True)
    args = parser.parse_args()
    output = setup("gen_comparisons")
    print(output, flush=True)
    full = training_input()
    controlled = training_input("controlled")
    studied = training_input("studied")
    response_audit(full, controlled, studied, output / "response_audit")
    for name, response, archive, iterations in (
        ("self_closure", full, args.full_mc, (3, 4)),
        ("held_out_closure", controlled, args.split_mc, (1, 3, 4)),
    ):
        directory = output / name
        directory.mkdir()
        mc_comparisons(response, archive, directory, iterations)
    toy_comparisons(full, args.cp5_toys, output / "cp5_toys", 1)
    data_comparisons(full, args.stitched_data, output / "exploratory_data", 1)


if __name__ == "__main__":
    main()
