"""Stored nPDF comparisons using shared prepared Data/PYTHIA histograms.

The input schema is the existing EPPS21 directory/TH1/TGraphAsymmErrors schema.
Register another calculation only after verifying its observable definitions,
normalization, eta frame, central keys and independently provided total bands.
"""
from pathlib import Path
import warnings

from hist_analysis.config.histograms import NPDF_PTAVE_BINS
from hist_analysis.config.theory import theory_file
from hist_analysis.python.histogram_io import load_histogram
from hist_analysis.python.data_mc_plotting import (
    accepted_eta_range, draw_comparison_figure,
)
from hist_analysis.python.theory_plotting import validate_total_band


def prediction_registry(project_root):
    return {
        'epps21_ct18': dict(label='EPPS21 x CT18', color=3,
            path=theory_file('epps21_ct18', project_root),
            intervals=(*NPDF_PTAVE_BINS, (50, 70)),
            keys={'full': ('ratio_central', 'ratio_total_band'),
                  'fb': ('forward_backward_pPb_over_pp', 'forward_backward_total_band')}),
        'ncteq15hq': dict(label='nCTEQ15HQ', color=6,
            path=theory_file('ncteq15hq', project_root),
            intervals=(*NPDF_PTAVE_BINS, (50, 70)),
            keys={'full': ('ratio_central', 'ratio_total_band'),
                  'fb': ('forward_backward_pPb_over_pp', 'forward_backward_total_band')}),
    }


def select_predictions(registry, selected):
    if len(set(selected)) != len(selected):
        raise ValueError('Select unique nPDF calculations')
    available = {}
    for model in selected:
        if model not in registry:
            raise ValueError(f'Unknown nPDF calculation: {model}')
        definition = registry[model]
        if definition is None or not Path(definition['path']).is_file():
            warnings.warn(f'Unavailable optional nPDF calculation: {model}')
            continue
        available[model] = definition
    return available


def load_prediction(definition, interval, observable, eta_cut):
    import ROOT
    directory = f'pt{interval[0]:03d}_{interval[1]:03d}_eta{eta_cut:.1f}'.replace('.', 'p')
    if interval not in definition['intervals']:
        raise ValueError(f'No exact {definition["label"]} prediction for {interval}')
    central_key, band_key = definition['keys'][observable]
    central = load_histogram(str(definition['path']), f'{directory}/{central_key}')
    file = ROOT.TFile.Open(str(definition['path']), 'READ')
    if not file or file.IsZombie():
        raise OSError(f'Cannot open {definition["path"]}')
    try:
        confidence = file.Get('pdf_confidence_level')
        if not confidence or confidence.GetVal() != 68:
            raise ValueError(f'{definition["path"]}: expected stored 68% PDF confidence level')
        source = file.Get(f'{directory}/{band_key}')
        if not source or not source.InheritsFrom('TGraphAsymmErrors'):
            raise KeyError(f'{definition["path"]}:{directory}/{band_key}')
        # Keep the historical prefix recognized by the ROOT/MPLHEP renderer.
        band = source.Clone(f'g_epps21_total_{definition["label"]}_{observable}_{directory}')
    finally:
        file.Close()
    validate_total_band(central, band, context=f'{definition["label"]}/{directory}/{band_key}')
    band.SetFillColorAlpha(definition['color'], .32)
    band.SetLineColor(definition['color']); band.SetFillStyle(1001); band.SetMarkerSize(0)
    return central, band


def statistical_band(histogram, *, name, color, hatch=1001):
    """Independent bin rectangles centered on the unchanged MC F/B values."""
    import ROOT
    band = ROOT.TGraphAsymmErrors(histogram)
    band.SetName('g_mc_stat_total_' + name)
    band.SetFillColorAlpha(color, .26); band.SetLineColor(color)
    band.SetFillStyle(hatch); band.SetMarkerSize(0)
    return band


def draw_npdf_comparisons(entry, predictions, config):
    """Draw Data/PYTHIA + stored pPb/pp, and Data/PYTHIA F/B + stored theory F/B."""
    interval, result = entry['ptave_range'], entry['result']
    data = config['DATA_LABEL']; eta_cut = config['ETA_CUT_VALUE']
    models = config['COMPARISON_MODELS']; registry = config['MODEL_REGISTRY']
    labels = tuple(registry[m]['label'] for m in models)
    canvases = {}
    for observable in ('full', 'fb'):
        source = result['ratios'] if observable == 'full' else result['forward_backward']
        missing = set(labels) - source.keys()
        if missing:
            raise ValueError(f'Missing selected PYTHIA inputs for {interval}/{observable}: {sorted(missing)}')
        reference = result['histograms'][data] if observable == 'full' else result['forward_backward'][data]
        x_range = accepted_eta_range(reference, eta_cut, observable=observable)
        bands = []
        for model, definition in predictions.items():
            if interval not in definition['intervals']:
                warnings.warn(f'{model}: no exact prediction for {interval}'); continue
            central, band = load_prediction(definition, interval, observable, eta_cut)
            for i in reversed(range(band.GetN())):
                x = band.GetPointX(i)
                if not x_range[0] - 1e-9 <= x <= x_range[1] + 1e-9:
                    band.RemovePoint(i); continue
                band.SetPointError(i, min(band.GetErrorXlow(i), max(0., x-x_range[0])),
                    min(band.GetErrorXhigh(i), max(0., x_range[1]-x)),
                    band.GetErrorYlow(i), band.GetErrorYhigh(i))
            suffix = ' pPb/pp' if observable == 'full' else ' (pPb/pp) F/B'
            bands.append((band, definition['label'] + suffix))
        total = config['data_totals'](entry['sample'], interval, observable)
        if observable == 'full':
            if not labels:
                continue
            plotted = {f'Data / {label}': source[label] for label in labels}
            styles = {f'Data / {label}': config['CURVE_STYLES'][label] for label in labels}
            systematic_keys = tuple(plotted)
        else:
            # MC histograms are represented ONLY by graphs, so neither backend
            # can accidentally redraw markers or conventional error bars.
            plotted = {data: reference}; styles = config['CURVE_STYLES']
            systematic_keys = (data,)
            for index, label in enumerate(labels):
                band = statistical_band(source[label], name=f'{index}_{interval}',
                    color=styles[label][0], hatch=1001 if index == 0 else 3004)
                bands.append((band, config['display_label'](label) + ' stat.'))
        model_tag = '_'.join(predictions) or 'no_theory'
        selection = config['DATA_SELECTION']
        eta_tag = f'{eta_cut:g}'.replace('.', 'p')
        stem = (f'npdf_{observable}_{model_tag}_{selection}_etaCM{eta_tag}'
                f'_ptave_{interval[0]:g}_{interval[1]:g}')
        canvas = draw_comparison_figure(plotted, name='c_'+stem,
            x_title='#eta_{CM}^{dijet}' if observable == 'full' else '|#eta_{CM}^{dijet}|',
            y_title='Data / MC' if observable == 'full' else 'Forward / Backward',
            x_range=x_range, styles=styles,
            annotations=config['selection_lines'](config['PROJECT_ROOT'], interval,
                eta_cut=eta_cut, direction=config['DATA_DIRECTION'], selection=config['DATA_SELECTION']),
            totals=tuple((key, total) for key in systematic_keys) if total is not None else (),
            y_range=config['NPDF_Y_RANGES'][observable], reference=True,
            bands=bands, grid=config['DRAW_GRID'], legend_text_size=.027,
            output=config['OUTPUT_DIR']/f'{stem}.pdf' if config['SAVE_PDF'] else None,
            save_png=config['SAVE_PNG'], labels=config['display_label'])
        canvases[observable] = canvas
    return canvases
