"""Optional systematic checks, isolated from nominal histograms and CSVs."""
from pathlib import Path

from hist_analysis.python.histogram_ops import ratio_to_nominal
from hist_analysis.python.systematic_methods import construct_systematic_variations
from hist_analysis.python.systematic_fits import (
    calculate_bin_by_bin_systematic, calculate_one_sided_systematic,
    fit_histogram_variations, smooth_systematic_running_max,
)


def estimate_ratio_error_variant(
    numerators, denominator, *, ratio_option, method, extraction,
    formula, fit_options, initial_values, evaluation_range,
    apply_smoothing, smoothing_origin, name_prefix, combination=None,
):
    """Recompute one diagnostic variant from the existing projected inputs.

    Uses TH1::Divide(...,1,1,option) through ratio_to_nominal, not a Python
    approximation of option B. Two insertion-ordered numerators mean Up/Down;
    a single numerator means the source's one-sided estimate. Only fresh ratio
    histograms are fitted. The inputs and all nominal objects remain untouched.
    """
    if extraction not in ('fit', 'bin_by_bin'):
        raise ValueError('extraction must be fit or bin_by_bin')
    if len(numerators) not in (1, 2):
        raise ValueError('Supply one one-sided ratio or two ordered Up/Down ratios')
    if len(numerators) == 2 and combination not in ('average', 'maximum'):
        raise ValueError('Two-sided diagnostics require average or maximum')
    raw_ratios = {
        label: ratio_to_nominal(histogram, denominator,
            name=f'{name_prefix}_ratio_{index}', option=ratio_option)
        for index, (label, histogram) in enumerate(numerators.items())
    }
    ratios = construct_systematic_variations(raw_ratios, method=method)
    functions, summaries = ({}, {})
    if extraction == 'fit':
        functions, summaries = fit_histogram_variations(
            ratios, formula=formula, fit_options=fit_options,
            fit_range=evaluation_range, initial_values=initial_values,
            name_prefix=f'{name_prefix}_fit',
        )
    labels = tuple(ratios)
    if len(labels) == 2:
        up, down = labels
        unsmoothed = calculate_bin_by_bin_systematic(
            ratios[up], ratios[down], name=f'{name_prefix}_unsmoothed',
            up_function=functions.get(up), down_function=functions.get(down),
            evaluation_range=evaluation_range, combination=combination,
        )
    else:
        label = labels[0]
        unsmoothed = calculate_one_sided_systematic(
            ratios[label], name=f'{name_prefix}_unsmoothed',
            variation_function=functions.get(label), evaluation_range=evaluation_range,
        )
    systematic = (
        smooth_systematic_running_max(unsmoothed,
            name=f'{name_prefix}_smoothed', evaluation_range=evaluation_range,
            smoothing_origin=smoothing_origin)
        if apply_smoothing else unsmoothed
    )
    return dict(raw_ratios=raw_ratios, ratios=ratios, functions=functions,
        fit_summaries=summaries, unsmoothed=unsmoothed, systematic=systematic,
        ratio_option=ratio_option)


def make_systematic_diagnostics(
    *, unsmoothed, nominal, numerators, denominator,
    plot_smoothed_vs_unsmoothed, plot_binomial_vs_standard_errors,
    method, extraction, formula, fit_options, initial_values,
    evaluation_range, apply_smoothing, smoothing_origin, combination=None,
    renderer, nominal_csv, source, observable, x_title, x_range, y_range,
    annotations, style, grid=False, save_png=False,
):
    """Plot fractional estimates on percent clones; never write a CSV.

    Smoothing comparison reuses the retained pre-smoothing estimate. If
    smoothing is nominal, its retained result is reused too. Otherwise only
    its would-be smoothed clone is computed. Error-treatment variants keep the
    selected construction, fit/raw mode, combination and nominal smoothing.
    """
    if not (plot_smoothed_vs_unsmoothed or plot_binomial_vs_standard_errors):
        return {}
    prefix = Path(nominal_csv).with_suffix('')
    name = prefix.name
    method_label = 'Folded deviations' if method == 'folded_deviation' else 'Direct ratios'
    extraction_label = 'Fit-based' if extraction == 'fit' else 'Bin-by-bin'
    estimator_label = f'; Up/Down {combination}' if combination else ''
    context = (*annotations, f'{source}: {observable}',
        f'{method_label}; {extraction_label}{estimator_label}')
    diagnostics = {}

    def draw(curves, kind):
        percent = {}
        for index, (label, histogram) in enumerate(curves.items()):
            clone = histogram.Clone(f'{name}_diagnostic_{kind}_percent_{index}')
            clone.SetDirectory(0)
            clone.Scale(100.0)
            percent[label] = clone
        output = prefix.with_name(f'{name}_diagnostic_{kind}.pdf')
        canvas = renderer(percent, title='', x_title=x_title,
            y_title='Relative systematic uncertainty (%)',
            x_range=x_range, y_range=y_range,
            annotations=context + ((('Outward running maximum' if apply_smoothing
                else 'No post-estimator smoothing'),) if kind == 'binomial_vs_standard_errors' else ()),
            grid=grid,
            style_indices={label: index for index, label in enumerate(curves)},
            style=style, output=output, save_png=save_png,
            mplhep=True, mplhep_png=save_png,
            canvas_name=f'{name}_diagnostic_{kind}')
        return dict(canvas=canvas, percent_curves=percent, output=output)

    if plot_smoothed_vs_unsmoothed:
        smoothed = nominal if apply_smoothing else smooth_systematic_running_max(
            unsmoothed, name=f'{name}_diagnostic_smoothed',
            evaluation_range=evaluation_range, smoothing_origin=smoothing_origin,
        )
        curves = {'Unsmoothed': unsmoothed, 'Smoothed': smoothed}
        diagnostics['smoothed_vs_unsmoothed'] = dict(
            curves=curves, **draw(curves, 'smoothed_vs_unsmoothed'))

    if plot_binomial_vs_standard_errors:
        variants = {}
        for label, option, tag in (
            ('Standard ratio errors', '', 'standard'),
            ('Binomial ratio errors ("B")', 'B', 'binomial'),
        ):
            variants[label] = estimate_ratio_error_variant(
                numerators, denominator, ratio_option=option,
                method=method, extraction=extraction, formula=formula,
                fit_options=fit_options, initial_values=initial_values,
                evaluation_range=evaluation_range, apply_smoothing=apply_smoothing,
                smoothing_origin=smoothing_origin, combination=combination,
                name_prefix=f'{name}_diagnostic_{tag}',
            )
        curves = {label: variant['systematic'] for label, variant in variants.items()}
        diagnostics['binomial_vs_standard_errors'] = dict(
            variants=variants, curves=curves,
            **draw(curves, 'binomial_vs_standard_errors'))
    return diagnostics
