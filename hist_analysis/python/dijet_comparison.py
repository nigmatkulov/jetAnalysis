"""Shared Data/MC preparation, normalization and rendering workflow.

Each call creates independent closures and a detached TH2 cache. Configuration
is supplied explicitly by the notebook; no analysis state is shared across runs.
"""
from IPython.display import display
from pathlib import Path
from hist_analysis.python.histogram_production import ETA_KEY_PATTERN, load_cached_eta_projection
from hist_analysis.config.histograms import SYSTEMATIC_PTAVE_BINS


def build_comparison_workflow(configuration):
    COMPARISONS = configuration['COMPARISONS']
    CP5_GEN_FILE = configuration['CP5_GEN_FILE']
    CP5_GEN_HISTOGRAM = configuration['CP5_GEN_HISTOGRAM']
    CP5_GEN_LABEL = configuration['CP5_GEN_LABEL']
    CUEP8M1_GEN_FILE = configuration['CUEP8M1_GEN_FILE']
    CUEP8M1_GEN_HISTOGRAM = configuration['CUEP8M1_GEN_HISTOGRAM']
    CUEP8M1_GEN_LABEL = configuration['CUEP8M1_GEN_LABEL']
    CURVE_STYLES = configuration['CURVE_STYLES']
    DATA_DIRECTION = configuration['DATA_DIRECTION']
    DATA_FILES = configuration['DATA_FILES']
    DATA_LABEL = configuration['DATA_LABEL']
    DATA_SELECTION = configuration['DATA_SELECTION']
    DRAW_GRID = configuration['DRAW_GRID']
    DRAW_MC_FIGURES = configuration['DRAW_MC_FIGURES']
    EMBEDDING_GEN_LABEL = configuration['EMBEDDING_GEN_LABEL']
    ETA_CUTS = configuration['ETA_CUTS']
    ETA_CUT_INDICES = configuration['ETA_CUT_INDICES']
    ETA_NORMALIZATION = configuration['ETA_NORMALIZATION']
    ETA_Y_TITLE = configuration['ETA_Y_TITLE']
    FB_MODELS = configuration['FB_MODELS']
    FULL_MODELS = configuration['FULL_MODELS']
    LOG_Y = configuration['LOG_Y']
    MC_FILE = configuration['MC_FILE']
    MC_RECO_LABEL = configuration['MC_RECO_LABEL']
    MODEL_REGISTRY = configuration['MODEL_REGISTRY']
    NORMALIZATION_RANGE = configuration['NORMALIZATION_RANGE']
    OUTPUT_DIR = configuration['OUTPUT_DIR']
    PLOT_Y_LIMITS = configuration['PLOT_Y_LIMITS']
    PROJECT_ROOT = configuration['PROJECT_ROOT']
    PTAVE_BINS = configuration['PTAVE_BINS']
    PTAVE_DISPLAY_RANGE = configuration['PTAVE_DISPLAY_RANGE']
    RATIO_COMBINATIONS = configuration['RATIO_COMBINATIONS']
    RATIO_Y_RANGES = configuration['RATIO_Y_RANGES']
    REBIN_ETA = configuration['REBIN_ETA']
    REBIN_PTAVE = configuration['REBIN_PTAVE']
    SAVE_PDF = configuration['SAVE_PDF']
    SAVE_PNG = configuration['SAVE_PNG']
    accepted_eta_range = configuration['accepted_eta_range']
    apply_curve_style = configuration['apply_curve_style']
    check_distinct_generator_inputs = configuration['check_distinct_generator_inputs']
    check_distinct_generator_predictions = configuration['check_distinct_generator_predictions']
    data_over_mc = configuration['data_over_mc']
    data_totals = configuration['data_totals']
    display_label = configuration['display_label']
    draw_comparison_figure = configuration['draw_comparison_figure']
    forward_backward_from_full = configuration['forward_backward_from_full']
    load_comparison_source = configuration['load_comparison_source']
    normalize_histogram = configuration['normalize_histogram']
    project_semantic_th2 = configuration['project_semantic_th2']
    resolve_catalog_file = configuration['resolve_catalog_file']
    selection_lines = configuration['selection_lines']
    # Inputs are listed once. COMPARISONS alone selects the active MC labels.
    SAMPLE_INPUTS = {
        MC_RECO_LABEL: (MC_FILE, 'hRecoDijetPtEtaCMJerDefExtra_{eta_cut_index}'),
        EMBEDDING_GEN_LABEL: (MC_FILE, 'hGenDijetPtEtaCM_{eta_cut_index}'),
        CUEP8M1_GEN_LABEL: (CUEP8M1_GEN_FILE, CUEP8M1_GEN_HISTOGRAM),
        CP5_GEN_LABEL: (CP5_GEN_FILE, CP5_GEN_HISTOGRAM),
        DATA_LABEL: (None, 'hRecoDijetPtEtaCM_{eta_cut_index}'),
    }
    for model, definition in MODEL_REGISTRY.items():
        if 'source_id' in definition:
            SAMPLE_INPUTS[definition['label']] = (
                resolve_catalog_file(definition['source_id'], selection=DATA_SELECTION, require_exists=False),
                'hGenDijetPtEtaCM_{eta_cut_index}')
    requested_models = set(FULL_MODELS) | set(FB_MODELS) | {m for pair in RATIO_COMBINATIONS for m in pair if m != 'data'}
    ACTIVE_MC_LABELS = {label for labels in COMPARISONS.values() for label in labels}
    ACTIVE_MC_LABELS.update(MODEL_REGISTRY[m]['label'] for m in requested_models)
    _input_histograms = {}
    import warnings
    for label in tuple(ACTIVE_MC_LABELS):
        filename, template = SAMPLE_INPUTS[label]
        try:
            for index in ETA_CUT_INDICES:
                load_comparison_source(filename, template.format(eta_cut_index=index), _input_histograms)
        except (FileNotFoundError, OSError, KeyError, RuntimeError) as error:
            warnings.warn(f'Disabling unavailable optional model {label}: {error}')
            ACTIVE_MC_LABELS.remove(label)
    COMPARISONS = {tag: tuple(label for label in labels if label in ACTIVE_MC_LABELS)
                   for tag, labels in COMPARISONS.items()}
    COMPARISONS = {tag: labels for tag, labels in COMPARISONS.items() if labels}
    CURVES = tuple((label, *source) for label, source in SAMPLE_INPUTS.items()
                   if label == DATA_LABEL or label in ACTIVE_MC_LABELS)
    check_distinct_generator_inputs(CURVES, EMBEDDING_GEN_LABEL)

    # Exact-key preflight above also populates the shared detached-source cache.

    def _project_ptave(filename, key, name):
        source = load_comparison_source(filename, key, _input_histograms)
        projection = source.ProjectionX(
            name, 1, source.GetYaxis().GetNbins(), 'e',
        )
        projection.SetDirectory(0)
        if REBIN_PTAVE > 1:
            projection.Rebin(REBIN_PTAVE)
        return projection

    def _project_eta(filename, key, ptave_range, name):
        projection = None
        if (Path(filename).parent.name == 'produced_histograms'
                and (ETA_KEY_PATTERN.match(key) or key == 'hGenDijetPtEtaCM')
                and tuple(ptave_range) in SYSTEMATIC_PTAVE_BINS):
            try:
                projection = load_cached_eta_projection(filename, key, 'standard', tuple(ptave_range))
                projection.SetName(name)
            except KeyError:
                pass  # A compatible legacy cache retains its raw TH2.
        if projection is None:
            source = load_comparison_source(filename, key, _input_histograms)
            projection = project_semantic_th2(
                source, 'eta', ptave_range, name=name,
            )
        if REBIN_ETA > 1:
            projection.Rebin(REBIN_ETA)
        return projection

    def _normalize_in_value_range(histogram, low, high):
        axis = histogram.GetXaxis()
        bin_low = axis.FindBin(low + 0.001)
        bin_high = axis.FindBin(high - 0.001)
        integral = histogram.Integral(bin_low, bin_high)
        if integral <= 0.0:
            raise ValueError(
                f'Cannot normalize {histogram.GetName()}: integral over bins '
                f'{bin_low}--{bin_high} ({low:g}--{high:g} GeV) is {integral:g}'
            )
        histogram.Scale(1.0 / integral)
        return bin_low, bin_high, integral

    def _assert_compatible_binning(reference, candidate):
        if reference.GetNbinsX() != candidate.GetNbinsX():
            raise ValueError('Cannot overlay histograms with different bin counts')
        for edge_index in range(1, reference.GetNbinsX() + 2):
            reference_edge = reference.GetXaxis().GetBinLowEdge(edge_index)
            candidate_edge = candidate.GetXaxis().GetBinLowEdge(edge_index)
            if abs(reference_edge - candidate_edge) > 1.0e-9:
                raise ValueError('Cannot overlay histograms with different bin edges')

    def _prepared_ratios(histograms, tag):
        # One independent-error calculation per denominator, shared with all renderers.
        return {label: data_over_mc(histograms[DATA_LABEL], h,
                    name=f'{h.GetName()}_data_over_mc_{tag}', style=CURVE_STYLES[label])
                for label, h in histograms.items() if label != DATA_LABEL}


    def _draw_comparisons(histograms, ratios, *, sample, interval, observable, eta_cut_index):
        eta_cut = ETA_CUTS[eta_cut_index]
        ptave = observable == 'ptave'
        x_range = PTAVE_DISPLAY_RANGE if ptave else accepted_eta_range(
            histograms[DATA_LABEL], eta_cut, observable=observable)
        x_title = 'p_{T}^{ave} (GeV)' if ptave else (
            '#eta_{CM}^{dijet}' if observable == 'full' else '|#eta_{CM}^{dijet}|')
        y_title = ('1/N dN/dp_{T}^{ave}' if ptave else
                   ETA_Y_TITLE if observable == 'full' else 'Forward / Backward')
        lines = selection_lines(PROJECT_ROOT, interval, eta_cut=eta_cut,
                                direction=DATA_DIRECTION, selection=DATA_SELECTION)
        if ptave:
            lines = (*lines[:3], f'normalized in {interval[0]:g} #leq p_{{T}}^{{ave}} < {interval[1]:g} GeV')
        total = None if ptave else data_totals(sample, interval, observable)
        interval_tag = '' if ptave else f'_ptave_{interval[0]:g}_{interval[1]:g}'
        stem = f'{sample}_dijet_{observable}_etaCM_{int(round(10*eta_cut))}{interval_tag}'
        range_key = {'ptave':'ptave_data_over_mc', 'full':'eta_data_over_mc', 'fb':'fb_data_over_mc'}[observable]
        results = {}
        families = COMPARISONS or {'data_only': ()}
        for tag, mc_labels in families.items():
            selected = {label: histograms[label] for label in (*mc_labels, DATA_LABEL)}
            selected_ratios = {f'Data / {label}': ratios[label] for label in mc_labels}
            ratio_styles = {f'Data / {label}': CURVE_STYLES[label] for label in mc_labels}
            overlay_path = OUTPUT_DIR / f'{stem}_data_vs_{tag}.pdf' if SAVE_PDF else None
            ratio_path = OUTPUT_DIR / f'{stem}_data_over_{tag}.pdf' if SAVE_PDF else None
            canvas = draw_comparison_figure(selected, name=f'c_{stem}_data_vs_{tag}',
                x_title=x_title, y_title=y_title, x_range=x_range, styles=CURVE_STYLES,
                annotations=lines, totals=((DATA_LABEL,total),) if total is not None else (),
                log_y=LOG_Y if ptave else False, grid=DRAW_GRID,
                y_range=PLOT_Y_LIMITS[observable] if not ptave else None,
                reference=observable == 'fb',
                output=overlay_path, save_png=SAVE_PNG, labels=display_label)
            display(canvas)
            ratio_canvas = None
            if selected_ratios:
                ratio_canvas = draw_comparison_figure(selected_ratios, name=f'c_{stem}_data_over_{tag}',
                    x_title=x_title, y_title='Data / MC', x_range=x_range, styles=ratio_styles,
                    annotations=lines, totals=tuple((key,total) for key in selected_ratios) if total is not None else (),
                    y_range=RATIO_Y_RANGES[range_key] if ptave else PLOT_Y_LIMITS['ratio'],
                    reference=True, grid=DRAW_GRID,
                    output=ratio_path, save_png=SAVE_PNG, labels=display_label)
                display(ratio_canvas)
            results[tag] = dict(histograms=selected, ratios=selected_ratios,
                                canvas=canvas, ratio_canvas=ratio_canvas)
        return results


    def _draw_eta_cut(sample, eta_cut_index):
        histograms, normalization_details = {}, {}
        for label, configured_file, key_template in CURVES:
            filename = DATA_FILES[sample] if configured_file is None else configured_file
            h = _project_ptave(filename, key_template.format(eta_cut_index=eta_cut_index),
                               f'h_{sample}_{label.lower().replace(" ", "_")}_eta_{eta_cut_index}')
            if histograms:
                _assert_compatible_binning(next(iter(histograms.values())), h)
            normalization_details[label] = _normalize_in_value_range(h, *NORMALIZATION_RANGE)
            apply_curve_style(h, *CURVE_STYLES[label])
            histograms[label] = h
        check_distinct_generator_predictions(histograms, EMBEDDING_GEN_LABEL,
                                             (CP5_GEN_LABEL,CUEP8M1_GEN_LABEL))
        ratios = _prepared_ratios(histograms, 'ptave')
        comparisons = _draw_comparisons(histograms, ratios, sample=sample,
            interval=NORMALIZATION_RANGE, observable='ptave', eta_cut_index=eta_cut_index)
        return dict(histograms=histograms, ratios_to_mc=ratios,
                    normalization=normalization_details, comparisons=comparisons)


    def analyze_sample(sample):
        if not PTAVE_BINS[sample]:
            return {}
        return {index: _draw_eta_cut(sample, index) for index in ETA_CUT_INDICES}


    def _draw_combined(histograms, *, sample, interval, observable, eta_cut_index):
        from hist_analysis.python.data_mc_plotting import prepare_configured_ratios
        models = FULL_MODELS if observable == 'full' else FB_MODELS
        selected = {MODEL_REGISTRY[m]['label']: histograms[MODEL_REGISTRY[m]['label']]
                    for m in models if MODEL_REGISTRY[m]['label'] in histograms}
        selected[DATA_LABEL] = histograms[DATA_LABEL]
        x_range = accepted_eta_range(histograms[DATA_LABEL], ETA_CUTS[eta_cut_index], observable=observable)
        lines = selection_lines(PROJECT_ROOT, interval, eta_cut=ETA_CUTS[eta_cut_index],
                                direction=DATA_DIRECTION, selection=DATA_SELECTION)
        total = data_totals(sample, interval, observable)
        stem = f'{sample}_{observable}_etaCM_{int(round(10*ETA_CUTS[eta_cut_index]))}_ptave_{interval[0]:g}_{interval[1]:g}'
        common = dict(x_title='#eta_{CM}^{dijet}' if observable == 'full' else '|#eta_{CM}^{dijet}|',
                      x_range=x_range, annotations=lines, grid=DRAW_GRID, labels=display_label)
        overlay = draw_comparison_figure(selected, name='c_combined_'+stem,
            y_title=ETA_Y_TITLE if observable == 'full' else 'Forward / Backward',
            styles=CURVE_STYLES, totals=((DATA_LABEL,total),) if total is not None else (),
            y_range=PLOT_Y_LIMITS[observable], reference=observable == 'fb',
            output=OUTPUT_DIR / f'{stem}_data_vs_combined.pdf' if SAVE_PDF else None,
            save_png=SAVE_PNG, **common)
        display(overlay)
        operands = {'data': histograms[DATA_LABEL]}
        operands.update({m: histograms[v['label']] for m,v in MODEL_REGISTRY.items() if v['label'] in histograms})
        prepared = prepare_configured_ratios(operands, RATIO_COMBINATIONS,
            styles={m: CURVE_STYLES[v['label']] for m,v in MODEL_REGISTRY.items()}, name=stem)
        ratios = {f"{DATA_LABEL if a == 'data' else MODEL_REGISTRY[a]['label']} / {MODEL_REGISTRY[b]['label']}": h
                  for (a,b),h in prepared.items()}
        systematic_keys = [f"Data / {MODEL_REGISTRY[b]['label']}" for a,b in prepared if a == 'data']
        ratio_canvas = None
        if ratios:
            ratio_canvas = draw_comparison_figure(ratios, name='c_combined_ratio_'+stem,
                y_title='Ratio', styles={k: (h.GetMarkerColor(), h.GetMarkerStyle()) for k,h in ratios.items()}, totals=tuple((k,total) for k in systematic_keys) if total is not None else (),
                y_range=PLOT_Y_LIMITS['ratio'], reference=True,
                output=OUTPUT_DIR / f'{stem}_configured_ratios.pdf' if SAVE_PDF else None,
                save_png=SAVE_PNG, **common)
            display(ratio_canvas)
        return dict(canvas=overlay, ratio_canvas=ratio_canvas, ratios=prepared)


    def _draw_eta_distribution(sample, eta_cut_index, ptave_range):
        low, high = ptave_range
        raw_histograms, histograms = {}, {}
        for label, configured_file, key_template in CURVES:
            filename = DATA_FILES[sample] if configured_file is None else configured_file
            h = _project_eta(filename, key_template.format(eta_cut_index=eta_cut_index), ptave_range,
                f'h_{sample}_{label.lower().replace(" ", "_")}_etaCM_{eta_cut_index}_ptave_{low:g}_{high:g}')
            if raw_histograms:
                _assert_compatible_binning(next(iter(raw_histograms.values())), h)
            raw_histograms[label] = h
            histograms[label] = normalize_histogram(h, ETA_NORMALIZATION)
        check_distinct_generator_predictions(histograms, EMBEDDING_GEN_LABEL,
                                             (CP5_GEN_LABEL,CUEP8M1_GEN_LABEL))
        # F/B uses the unnormalized, disjoint hemispheres, never normalized shapes.
        fb = {label: forward_backward_from_full(h,
            name=f'h_{sample}_{label.lower().replace(" ", "_")}_fb_{eta_cut_index}_ptave_{low:g}_{high:g}')
            for label,h in raw_histograms.items()}
        ratios, fb_ratios = _prepared_ratios(histograms,'full'), _prepared_ratios(fb,'fb')
        full_figures = {} if not DRAW_MC_FIGURES else _draw_comparisons(histograms,ratios,sample=sample,interval=ptave_range,
                                        observable='full',eta_cut_index=eta_cut_index)
        fb_figures = {} if not DRAW_MC_FIGURES else _draw_comparisons(fb,fb_ratios,sample=sample,interval=ptave_range,
                                      observable='fb',eta_cut_index=eta_cut_index)
        combined = {obs: _draw_combined(h, sample=sample, interval=ptave_range, observable=obs, eta_cut_index=eta_cut_index)
                    for obs,h in (('full',histograms),('fb',fb))} if DRAW_MC_FIGURES else {}
        return dict(combined=combined, raw_histograms=raw_histograms,histograms=histograms,forward_backward=fb,
                    ratios=ratios,fb_ratios=fb_ratios,comparisons=full_figures,fb_comparisons=fb_figures)


    def analyze_eta_sample(sample):
        return {(index,interval): _draw_eta_distribution(sample,index,interval)
                for index in ETA_CUT_INDICES for interval in PTAVE_BINS[sample]}

    return {name: value for name, value in locals().items() if name.startswith('_') or name in ('CURVES', 'analyze_sample', 'analyze_eta_sample')}
