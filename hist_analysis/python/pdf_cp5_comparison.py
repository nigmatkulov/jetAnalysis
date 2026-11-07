"""pp member-0 theory versus standalone CP5 using unchanged shared operations."""
import math

from hist_analysis.config.theory import theory_file
from hist_analysis.config.histograms import NPDF_PTAVE_BINS
from hist_analysis.python.histogram_io import load_histogram, open_root_file
from hist_analysis.python.histogram_ops import (
    align_common_binning, normalize_histogram, ratio_to_nominal,
)
from hist_analysis.python.projections import project_semantic_th2
from hist_analysis.python.forward_backward import fold_eta_histogram, forward_backward_ratio
from hist_analysis.python.output_naming import verified_reference_eta_cuts

PDF_MODELS = {'CT18': 'epps21_ct18', 'nCTEQ15HQ': 'ncteq15hq'}
PDF_COLORS = {'CT18': 418, 'nCTEQ15HQ': 6}  # existing green+2 / magenta
PDF_INTERVALS = (*NPDF_PTAVE_BINS, (50, 70))
PDF_CENTRAL_KEYS = ('pdf_central', 'pdf_central_mirrored')


def load_pdf_central(pdf_family, interval, *, eta_cut=1.9, project_root=None,
                     pdf_key='pdf_central_mirrored'):
    """Load the explicitly selected stored pp central; never substitute keys."""
    if pdf_key not in PDF_CENTRAL_KEYS:
        raise ValueError(f'Unknown PDF central key: {pdf_key}')
    if pdf_family not in PDF_MODELS:
        raise ValueError(f'Unknown PDF family: {pdf_family}')
    if tuple(interval) not in PDF_INTERVALS or eta_cut != 1.9:
        raise ValueError(f'No exact stored theory prediction for {interval}, eta={eta_cut}')
    path = theory_file(PDF_MODELS[pdf_family], project_root)
    directory = f'pt{interval[0]:03d}_{interval[1]:03d}_eta1p9'
    file = open_root_file(str(path))
    try:
        definition = file.Get('central_distribution_definition')
        if not definition or 'pdf_central=pp member 0' not in definition.GetTitle():
            raise ValueError(f'{path}: pp central definition is unverified')
        if pdf_key == 'pdf_central_mirrored':
            mirrored = file.Get('central_distribution_covariance_definition')
            if not mirrored or 'pdf_central_mirrored uses M=(I+P)/2' not in mirrored.GetTitle():
                raise ValueError(f'{path}: mirrored pp central definition is unverified')
    finally:
        file.Close()
    histogram = load_histogram(str(path), directory+'/'+pdf_key)
    # The stored 50--70 combination retains its first 50--60 component title.
    # Its exact directory and the file's width-weighted combination metadata
    # identify the combined interval; no adjacent histogram is substituted.
    title_directory = 'pt050_060_eta1p9' if tuple(interval) == (50, 70) else directory
    expected = f'hDijetEta_pp_cm_{title_directory}_pdf000_mu1p0'
    if histogram.GetDimension() != 1 or histogram.GetTitle() != expected:
        raise ValueError(f'{path}:{directory}: unexpected central spectrum definition')
    return histogram


def load_cp5_source(filename, *, eta_cut=1.9):
    key = 'hGenDijetPtEtaCM'
    cuts = verified_reference_eta_cuts([('CP5', filename, key)])
    if 'CP5' not in cuts or not math.isclose(float(cuts['CP5']), eta_cut, abs_tol=1e-9):
        raise ValueError(f'{filename}: CP5 acceptance {cuts} differs from eta={eta_cut}')
    source = load_histogram(str(filename), key)
    if source.GetDimension() != 2 or 'CM' not in source.GetYaxis().GetTitle():
        raise ValueError(f'{filename}:{key}: expected pTave versus etaCM TH2')
    return source


def prepare_pdf_cp5(pdf_family, interval, cp5_source, *, eta_cut=1.9,
                    normalization='integral', project_root=None,
                    pdf_key='pdf_central_mirrored'):
    """Sum whole common eta bins; normalize shapes only after raw F/B folding.

    Source pp theory eta bins are uniform (0.2). Thus an overall source density
    factor cancels from both normalized shapes and F/B. CP5 weighted yields
    have 0.1 bins, combined by summation with independent variances. No partial
    source bins or interpolated values are used, including at the acceptance.
    """
    if normalization not in {'integral', 'bin_width'}:
        raise ValueError('Compare shapes using integral or bin_width normalization')
    axis = cp5_source.GetXaxis()
    edges = [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins()+2)]
    if not all(any(math.isclose(value, edge, abs_tol=1e-8) for edge in edges)
               for value in interval):
        raise ValueError(f'CP5 pTave bins do not exactly cover {interval}')
    theory = load_pdf_central(pdf_family, interval, eta_cut=eta_cut,
                             project_root=project_root, pdf_key=pdf_key)
    widths = [theory.GetBinWidth(i) for i in range(1, theory.GetNbinsX()+1)]
    if not all(math.isclose(w, widths[0], abs_tol=1e-9) for w in widths):
        raise ValueError('Nonuniform theory bins require an explicit units convention')
    cp5 = project_semantic_th2(cp5_source, 'eta', interval,
                              name=f'cp5_pt{interval[0]}_{interval[1]}')
    original = {pdf_family: theory, 'PYTHIA CP5': cp5}
    for label, h in original.items():
        if any(not math.isfinite(h.GetBinContent(i)) or not math.isfinite(h.GetBinError(i))
               for i in range(1, h.GetNbinsX()+1)):
            raise ValueError(f'{label}/{interval}: nonfinite spectrum')
    raw = align_common_binning(original)
    summary = {'ptave': tuple(interval), 'normalization': normalization,
               'pdf_key': pdf_key, 'samples': {}}
    fb, shapes = {}, {}
    for label, h in raw.items():
        before = original[label].Integral()
        if not math.isclose(before, h.Integral(), rel_tol=1e-10, abs_tol=1e-9):
            raise ValueError(f'{label}: common eta range would discard nonzero yield')
        f, b = fold_eta_histogram(h, name=h.GetName()+'_fb')
        # Project F/B convention uses independent marginal errors even for a
        # mirrored theory central. M=(I+P)/2 correlates reflected bins; that
        # covariance is intentionally not imported here, so unity F/B can have
        # nonzero displayed errors. These are not a full theory uncertainty.
        ratio = forward_backward_ratio(f, b, name=h.GetName()+'_RFB')
        ratio._valid_ratio_bins = tuple(i for i in range(1, ratio.GetNbinsX()+1)
                                       if b.GetBinContent(i) != 0)
        # Check ROOT's weighted independent error formula on every usable bin.
        for i in ratio._valid_ratio_bins:
            F, B, sF, sB = f.GetBinContent(i), b.GetBinContent(i), f.GetBinError(i), b.GetBinError(i)
            expected = math.hypot(sF/B, F*sB/(B*B))
            if not math.isclose(ratio.GetBinError(i), expected, rel_tol=1e-10, abs_tol=1e-12):
                raise AssertionError('ROOT F/B error propagation mismatch')
        shapes[label] = normalize_histogram(h, normalization)
        after = shapes[label].Integral('width' if normalization == 'bin_width' else '')
        if not math.isclose(after, 1., abs_tol=1e-10):
            raise AssertionError('Shape normalization failed')
        fb[label] = ratio
        summary['samples'][label] = dict(raw_integral=before, aligned_integral=h.Integral(),
            normalized_integral=after, bins=h.GetNbinsX(),
            zero_denominator_bins=[i for i in range(1, ratio.GetNbinsX()+1)
                                   if i not in ratio._valid_ratio_bins])
    comparison_ratios = {}
    for observable, histograms in (('shapes', shapes), ('fb', fb)):
        numerator, denominator = histograms['PYTHIA CP5'], histograms[pdf_family]
        ratio = ratio_to_nominal(numerator, denominator,
            name=f'{numerator.GetName()}_over_{pdf_family}_{observable}', option='')
        valid = set(range(1, ratio.GetNbinsX()+1))
        for source in (numerator, denominator):
            valid.intersection_update(getattr(source, '_valid_ratio_bins', valid))
        ratio._valid_ratio_bins = tuple(i for i in sorted(valid)
                                       if denominator.GetBinContent(i) != 0)
        comparison_ratios[observable] = ratio
    summary['comparison_ratio_invalid_bins'] = {
        observable: [i for i in range(1, ratio.GetNbinsX()+1)
                     if i not in ratio._valid_ratio_bins]
        for observable, ratio in comparison_ratios.items()
    }
    return dict(raw=raw, shapes=shapes, fb=fb,
                comparison_ratios=comparison_ratios, validation=summary)
