"""Construction choices for systematic fit and bin-by-bin inputs."""
import math

SYSTEMATIC_METHODS = ('direct_ratio', 'folded_deviation')


def validate_systematic_method(method):
    if method not in SYSTEMATIC_METHODS:
        raise ValueError(f'SYSTEMATIC_METHOD must be one of {SYSTEMATIC_METHODS}, got {method!r}')
    return method


def construct_systematic_variations(ratios, *, method):
    """Keep direct inputs or reflect finite TH1 ratios about unity.

    ROOT division stores missing/empty bins as content=error=0. Preserve those
    sentinels (including flows), so they remain excluded from weighted fits and
    retain the existing raw estimator's empty-bin behavior. Nonfinite inputs
    fail explicitly in folded mode. The current fit interface uses symmetric
    TH1 errors; no asymmetric interval or folded-normal likelihood is invented.
    """
    validate_systematic_method(method)
    if method == 'direct_ratio':
        return ratios
    folded = {}
    for label, ratio in ratios.items():
        if ratio.GetBinErrorOption() != 0:
            raise ValueError('Folded systematic inputs require symmetric TH1 bin errors')
        histogram = ratio.Clone(f'{ratio.GetName()}_folded_deviation')
        histogram.SetDirectory(0)
        for index in range(ratio.GetNbinsX() + 2):
            value, error = ratio.GetBinContent(index), ratio.GetBinError(index)
            if not math.isfinite(value) or not math.isfinite(error):
                raise ValueError(f'Nonfinite systematic ratio/error in {label}, bin {index}')
            if value == 0.0 and error == 0.0:
                continue
            histogram.SetBinContent(index, 1.0 + abs(1.0 - value))
            histogram.SetBinError(index, error)
        folded[label] = histogram
    return folded


def systematic_variation_labels(values, *, method):
    """Translate display keys only; numerical fit/initial-value keys stay stable."""
    validate_systematic_method(method)
    if method == 'direct_ratio':
        return values
    return {f'1 + |1 - {label}|': value for label, value in values.items()}
