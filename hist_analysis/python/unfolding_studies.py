"""Controlled correction-chain diagnostics for the joint dijet pTave–eta response."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from hist_analysis.python.data_distributions import forward_backward_from_full
from hist_analysis.python.mplhep_plotting import save_comparison
from hist_analysis.python.unfolding import (
    apply_efficiency_correction,
    extract_eta_block,
    unfold_bayes,
    write_unfolding_output,
)


CHAIN_LABELS = {
    "purity_only": "A · Purity → unfold",
    "efficiency_only": "B · Unfold → efficiency",
    "uncorrected": "C · Unfold",
    "fully_corrected": "Full · Purity → unfold → efficiency",
}


def _values(histogram):
    return np.array([histogram.GetBinContent(i) for i in range(1, histogram.GetNbinsX() + 1)])


def _variances(histogram):
    return np.array([histogram.GetBinError(i)**2 for i in range(1, histogram.GetNbinsX() + 1)])


def _fraction(numerator, numerator_variance, total, total_variance, *, tolerance=1e-5):
    """Weighted subset ratio with disjoint numerator and complement Sumw2."""
    for part, whole in ((numerator, total), (numerator_variance, total_variance)):
        limit = tolerance * np.maximum(np.abs(whole), np.abs(part)) + 1e-30
        if np.any((part < -limit) | (whole < -limit) | (part > whole + limit)):
            raise ValueError("Response subset exceeds inclusive yield or Sumw2")
    defined = total > 0
    fraction = np.divide(numerator, total, out=np.full_like(total, np.nan), where=defined)
    fraction = np.clip(fraction, 0, 1)
    complement_variance = np.maximum(total_variance - numerator_variance, 0)
    variance = np.divide(
        (1 - fraction)**2 * numerator_variance + fraction**2 * complement_variance,
        total**2, out=np.full_like(total, np.nan), where=defined,
    )
    return fraction, variance


def _response_quantities(bundle, diagnostics, truth, reco, explicit_miss, explicit_fake):
    response = bundle.response
    if response.HasFakes() or response.UseOverflowStatus():
        raise ValueError("Expected a response without fakes or flow bins")
    matrix = bundle.scaled_matrix
    size = truth.GetNbinsX()
    # Bayes normalizes the raw RooUnfold matrix by Vtruth and skips zero
    # truth priors. RooUnfold may pad isolated bins with a unit diagonal.
    raw_matrix = _matrix_array(response.Mresponse(False))
    truth_prior = np.array([response.Vtruth()[i] for i in range(size)])
    matched_truth = _values(diagnostics.matched_truth)
    if np.any(raw_matrix < 0) or not np.array_equal(truth_prior > 0, matched_truth > 0):
        raise ValueError("Invalid migration probabilities or truth support")
    populated = truth_prior > 0
    if not np.allclose(raw_matrix[:, populated].sum(axis=0), truth_prior[populated],
                       rtol=1e-10, atol=1e-12):
        raise ValueError("RooUnfold columns must sum to one on matched truth support")
    diagonal = np.array([matrix.GetBinContent(i, i) / bundle.scale for i in range(1, size + 1)])
    diagonal_variance = np.array([
        (matrix.GetBinError(i, i) / bundle.scale)**2 for i in range(1, size + 1)
    ])
    for inclusive, matched, explicit in (
        (truth, diagnostics.matched_truth, explicit_miss),
        (reco, diagnostics.matched_measured, explicit_fake),
    ):
        if explicit is not None:
            _fraction(_values(matched) + _values(explicit),
                      _variances(matched) + _variances(explicit),
                      _values(inclusive), _variances(inclusive))

    truth_values, truth_variances = _values(truth), _variances(truth)
    reco_values, reco_variances = _values(reco), _variances(reco)
    matched_truth_variance = _variances(diagnostics.matched_truth)
    matched_reco = _values(diagnostics.matched_measured)
    matched_reco_variance = _variances(diagnostics.matched_measured)
    efficiency = _fraction(matched_truth, matched_truth_variance, truth_values, truth_variances)
    purity = _fraction(matched_reco, matched_reco_variance, reco_values, reco_variances)
    quantities = {
        "efficiency": efficiency,
        "matched_fraction": purity,
        "fake_fraction": (1 - purity[0], purity[1]),
        "diagonal_purity": _fraction(diagonal, diagonal_variance, reco_values, reco_variances),
        "stability": _fraction(diagonal, diagonal_variance, matched_truth, matched_truth_variance),
    }
    diagonal_fractions = {
        "reco": _fraction(diagonal, diagonal_variance, matched_reco, matched_reco_variance)[0],
        "truth": quantities["stability"][0],
    }
    effective_entries = {
        side: np.divide(values**2, variance, out=np.zeros_like(values), where=variance > 0)
        for side, values, variance in (
            ("reco", matched_reco, matched_reco_variance),
            ("truth", matched_truth, matched_truth_variance),
        )
    }
    return quantities, effective_entries, diagonal_fractions


def _matrix_array(matrix):
    return np.array([[matrix[i][j] for j in range(matrix.GetNcols())]
                     for i in range(matrix.GetNrows())])


def _projection(histogram, covariance, template, layout, index, mask, name):
    # A selected interval is one block of the flattened 2D result. Its
    # covariance is the corresponding principal submatrix of the full one.
    projected = extract_eta_block(histogram, template, layout, index, name=name)
    start = index * layout.n_eta_bins
    for local in range(layout.n_eta_bins):
        global_index = start + local
        if not mask[global_index]:
            projected.SetBinContent(local + 1, 0)
            projected.SetBinError(local + 1, 0)
        elif covariance is not None:
            projected.SetBinError(local + 1, math.sqrt(max(covariance[global_index, global_index], 0)))
    return projected


def _neighbor_rows(stage, iterations, histogram, covariance, support, layout, templates,
                   quantities, effective_entries, diagonal_fractions, side):
    values = _values(histogram)
    variances = _variances(histogram) if covariance is None else covariance.diagonal()
    factor = quantities["matched_fraction" if side == "reco" else "efficiency"][0]
    diagonal_fraction = diagonal_fractions[side]
    for pt_index, (low, high) in enumerate(layout.pt_intervals):
        axis = templates[pt_index].GetXaxis()
        for local in range(1, layout.n_eta_bins):
            left = layout.global_bin(pt_index, local) - 1
            right = left + 1
            if not (support[left] and support[right]):
                continue
            a, b = values[left], values[right]
            cross = 0 if covariance is None else covariance[left, right]
            variance = variances[left] + variances[right] - 2 * cross
            yield {
                "stage": stage, "iterations": iterations, "pt_low": low, "pt_high": high,
                "eta_boundary": axis.GetBinUpEdge(local),
                "left_yield": a, "right_yield": b, "difference": b - a,
                "relative_difference": 2 * (b - a) / (a + b) if a + b > 0 else "",
                "difference_over_sigma": (b - a) / math.sqrt(variance) if variance > 0 else "",
                "left_relative_error": math.sqrt(max(variances[left], 0)) / a if a > 0 else "",
                "right_relative_error": math.sqrt(max(variances[right], 0)) / b if b > 0 else "",
                "left_factor": factor[left], "right_factor": factor[right],
                "left_diagonal_fraction": diagonal_fraction[left],
                "right_diagonal_fraction": diagonal_fraction[right],
                "left_matched_effective_entries": effective_entries[side][left],
                "right_matched_effective_entries": effective_entries[side][right],
            }


def run_correction_study(
    roounfold, bundle, diagnostics, factorized, matched, full, full_covariance,
    *, measured, training_reco, training_truth, reference_truth, layout, eta_templates,
    iterations, chains, scan_iterations=(), explicit_miss=None, explicit_fake=None,
    output_dir, output_tag, annotations, cms_header_left, cms_header_right,
    data_reference=False, combined_ratios=False, save_png=False,
):
    """Run the raw and purity-corrected inputs through the same matched response.

    The nominal purity-only and full results are supplied by the primary
    notebook calculation. Iteration scans repeat both inputs at each setting.
    Covariance is conditional on the trained response and fixed corrections.
    """
    if not chains or len(set(chains)) != len(chains) or set(chains) - CHAIN_LABELS.keys():
        raise ValueError("Select unique chains from CHAIN_LABELS")
    quantities, effective_entries, diagonal_fractions = _response_quantities(
        bundle, diagnostics, training_truth, training_reco, explicit_miss, explicit_fake,
    )
    size = layout.n_global_bins
    reco_training_support = _values(training_reco) > 0
    reco_support = _values(diagnostics.matched_measured) > 0
    truth_support = _values(diagnostics.matched_truth) > 0
    observed = (_values(measured) != 0) | (_variances(measured) != 0)
    unsupported_reco = (np.flatnonzero(~reco_support & observed) + 1).tolist()
    unsupported_truth = (np.flatnonzero(~truth_support & (_values(reference_truth) > 0)) + 1).tolist()
    if unsupported_reco or unsupported_truth:
        print(f"Response support (global ROOT bins): reco input without migration {unsupported_reco}; "
              f"Gen reference without migration {unsupported_truth}")
    if np.any(~reco_training_support & observed):
        print("Purity is undefined for observed bins with no training reco; "
              "the matched response cannot unfold those entries.")

    gen = reference_truth
    scale = 1.0
    if data_reference:
        scale = measured.Integral() / training_reco.Integral()
        gen = reference_truth.Clone(f"{output_tag}_gen_model")
        gen.SetDirectory(0)
        gen.Scale(scale)

    path = Path(output_dir)
    path.mkdir(parents=True, exist_ok=True)
    all_support = np.ones(size, dtype=bool)
    rows = []
    for stage, hist, mask, side in (
        ("reco", measured, all_support, "reco"),
        ("purity_corrected_reco", factorized.measured_signal, reco_training_support, "reco"),
        ("gen", gen, all_support, "truth"),
    ):
        rows.extend(_neighbor_rows(stage, "", hist, None, mask, layout, eta_templates,
                                   quantities, effective_entries, diagonal_fractions, side))

    nominal_results = {}
    eta_axis = eta_templates[0].GetXaxis()
    eta_edges = [eta_axis.GetBinLowEdge(i) for i in range(1, layout.n_eta_bins + 2)]
    factor_histograms = {}
    undefined_factors = {}
    for key, (value, variance) in quantities.items():
        histogram = training_truth.Clone(key)
        histogram.SetDirectory(0)
        histogram.Reset()
        undefined_factors[key] = (np.flatnonzero(~np.isfinite(value)) + 1).tolist()
        for index in np.flatnonzero(np.isfinite(value)):
            histogram.SetBinContent(int(index) + 1, value[index])
            histogram.SetBinError(int(index) + 1, math.sqrt(variance[index]))
        factor_histograms[key] = histogram
    write_unfolding_output(
        path / f"{output_tag}_response_factors.root",
        histograms=(*factor_histograms.values(), factorized.measured_signal),
        metadata={"pt_intervals": layout.pt_intervals, "eta_edges": eta_edges,
                  "undefined_factor_bins": undefined_factors},
    )

    for count in dict.fromkeys((iterations, *scan_iterations)):
        if count == iterations:
            purity_only = matched
            corrected, corrected_covariance = full, full_covariance
        else:
            purity_only = unfold_bayes(
                roounfold, bundle, factorized.measured_signal, iterations=count,
                handle_fakes=False, name=f"{output_tag}_purity_only_iter{count}",
            )
            corrected, corrected_covariance = apply_efficiency_correction(
                purity_only.histogram, purity_only.covariance, factorized.efficiency,
                name=f"{output_tag}_fully_corrected_iter{count}",
            )
        raw = unfold_bayes(
            roounfold, bundle, measured, iterations=count, handle_fakes=False,
            name=f"{output_tag}_uncorrected_iter{count}",
        )
        efficiency_only, efficiency_covariance = apply_efficiency_correction(
            raw.histogram, raw.covariance, factorized.efficiency,
            name=f"{output_tag}_efficiency_only_iter{count}",
        )
        results = {
            "purity_only": (purity_only.histogram, purity_only.covariance),
            "efficiency_only": (efficiency_only, efficiency_covariance),
            "uncorrected": (raw.histogram, raw.covariance),
            "fully_corrected": (corrected, corrected_covariance),
        }
        for chain in chains:
            hist, cov = results[chain]
            array = _matrix_array(cov)
            rows.extend(_neighbor_rows(chain, count, hist, array, truth_support, layout,
                                       eta_templates, quantities, effective_entries,
                                       diagonal_fractions, "truth"))
            if count == iterations:
                nominal_results[chain] = (hist, array)
                write_unfolding_output(
                    path / f"{output_tag}_{chain}.root",
                    histograms=(hist, measured, gen), covariance=cov,
                    metadata={
                        "chain": chain, "iterations": count,
                        "pt_intervals": layout.pt_intervals, "eta_edges": eta_edges,
                        "gen_scale": scale,
                        "unsupported_truth_bins": (np.flatnonzero(~truth_support) + 1).tolist(),
                        "unsupported_observed_reco_bins": unsupported_reco,
                        "undefined_purity_bins": (np.flatnonzero(~reco_training_support) + 1).tolist(),
                    },
                )

    with (path / f"{output_tag}_neighbor_differences.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]) if rows else [])
        writer.writeheader()
        writer.writerows(rows)

    _plot_study(path, output_tag, measured, gen, nominal_results, factor_histograms,
                truth_support, layout, eta_templates, annotations, cms_header_left,
                cms_header_right, data_reference, combined_ratios, save_png)
    return nominal_results


def _plot_study(path, tag, measured, gen, results, factors, support, layout,
                templates, annotations, left, right, data_reference, combined_ratios, save_png):
    diagnostic_titles = {
        "efficiency": ("Efficiency", "Fraction of Gen"),
        "matched_fraction": ("Matched fraction", "Fraction of Reco"),
        "fake_fraction": ("Fake fraction", "Fraction of Reco"),
        "diagonal_purity": ("Diagonal purity", "Fraction of Reco"),
        "stability": ("Stability", "Fraction of matched Gen"),
    }
    gen_label = "MC Gen reference" if data_reference else "Gen"
    folded_histograms = []
    undefined_bins = {}
    for index, ((low, high), template) in enumerate(zip(layout.pt_intervals, templates, strict=True)):
        interval = f"{low:g}_{high:g}"
        pt_line = rf"${low:g} \leq p_{{\mathrm{{T}}}}^{{\mathrm{{ave}}}} < {high:g}$ GeV"
        lines = (*annotations, pt_line)

        def save(histograms, suffix, y_label, *, extra=(), reference_y=None, y_range=None,
                 x_label=r"$\eta_{\mathrm{CM}}^{\mathrm{dijet}}$"):
            for histogram in histograms.values():
                histogram.SetMinimum()
                histogram.SetMaximum()
            save_comparison(
                histograms, output=path / f"{tag}_{suffix}_pt_{interval}.pdf",
                x_label=x_label, y_label=y_label,
                annotations=(*lines, *extra), save_png=save_png,
                cms_header_left=left, cms_header_right=right,
                headroom=1.65, reference_y=reference_y, y_range=y_range,
            )

        reco_eta = extract_eta_block(measured, template, layout, index, name=f"reco_pt{index}")
        gen_eta = extract_eta_block(gen, template, layout, index, name=f"gen_pt{index}")
        def fold(histogram, label, supported_pairs=None):
            fb = forward_backward_from_full(histogram, name=f"{label}_fb_pt{index}")
            forward, backward = fb._forward_backward_inputs
            undefined = []
            for bin_index in range(1, fb.GetNbinsX() + 1):
                supported = supported_pairs is None or supported_pairs[bin_index - 1]
                if not supported or backward.GetBinContent(bin_index) == 0:
                    undefined.append(bin_index)
                    fb.SetBinContent(bin_index, 0)
                    fb.SetBinError(bin_index, 0)
            undefined_bins[fb.GetName()] = undefined
            folded_histograms.extend((fb, forward, backward))
            return fb

        reco_fb = fold(reco_eta, "reco")
        gen_fb = fold(gen_eta, "gen")
        start = index * layout.n_eta_bins
        eta_support = support[start:start + layout.n_eta_bins]
        half = reco_fb.GetNbinsX()
        pair_support = eta_support[half:] & eta_support[:half][::-1]
        reco_ratio = reco_eta.Clone(f"reco_ratio_pt{index}")
        reco_ratio.SetDirectory(0)
        reco_ratio.Divide(gen_eta)
        unfolded_ratios = {}
        for chain, (hist, covariance) in results.items():
            eta = _projection(hist, covariance, template, layout, index, support,
                              f"{chain}_pt{index}")
            ratio = eta.Clone(f"{chain}_ratio_pt{index}")
            ratio.SetDirectory(0)
            ratio.Divide(gen_eta)
            extra = (CHAIN_LABELS[chain],)
            if data_reference:
                extra += ("MC Gen scaled globally using Reco yields",)
            save({"Reco": reco_eta, "Unfolded": eta, gen_label: gen_eta},
                 f"{chain}_eta", "Weighted yield per bin", extra=extra)
            save({f"Reco / {gen_label}": reco_ratio, f"Unfolded / {gen_label}": ratio},
                 f"{chain}_ratio_to_gen", "Ratio to Gen", extra=extra, reference_y=1)
            unfolded_fb = fold(eta, chain, pair_support)
            save({"Reco": reco_fb, "Unfolded": unfolded_fb, gen_label: gen_fb},
                 f"{chain}_forward_backward",
                 r"Forward ($\eta_{\mathrm{CM}}^{\mathrm{dijet}}>0$) / "
                 r"Backward ($\eta_{\mathrm{CM}}^{\mathrm{dijet}}<0$)",
                 x_label=r"$|\eta_{\mathrm{CM}}^{\mathrm{dijet}}|$",
                 extra=(CHAIN_LABELS[chain],), reference_y=1)
            unfolded_ratios[CHAIN_LABELS[chain]] = ratio
        if combined_ratios:
            save(unfolded_ratios, "all_chains_ratio_to_gen", "Unfolded / Gen",
                 reference_y=1)

        for key, histogram in factors.items():
            factor = extract_eta_block(histogram, template, layout, index, name=f"{key}_pt{index}")
            label, axis_label = diagnostic_titles[key]
            save({label: factor}, key, axis_label, y_range=(-0.05, 1.65))

    write_unfolding_output(
        path / f"{tag}_forward_backward.root",
        histograms=folded_histograms,
        metadata={"pt_intervals": layout.pt_intervals, "undefined_fb_bins": undefined_bins,
                  "errors": "ROOT independent forward/backward statistical errors"},
    )
