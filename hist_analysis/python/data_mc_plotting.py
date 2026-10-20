"""Display stored experimental totals without recomputing systematic combinations."""
from pathlib import Path
import json
import math
import re
import numpy as np
from hist_analysis.python.systematic_plotting import build_systematic_band
from hist_analysis.python.root_style import draw_text_block


def load_data_total(directory, interval, *, observable, direction, selection, binning,
                    eta_cut, sample, configuration=None):
    """Resolve a unique producer manifest and validate its fractional CSV and audit."""
    directory = Path(directory) / observable
    low, high = interval
    matches = []
    for path in directory.glob(f'*_ptave_{low:g}_{high:g}_provenance.json'):
        records = json.loads(path.read_text())
        expected = dict(observable=observable, frame='CM', direction=direction,
                        selection=selection, binning=binning, eta_cut=eta_cut,
                        units='fraction', convention='symmetric_relative',
                        target_sample='pPb 8.16 TeV dijets', pt_low=low, pt_high=high)
        if not records or not all(all(r.get(k) == v for k, v in expected.items()) for r in records):
            continue
        if configuration and configuration not in path.name:
            continue
        data_records = [r for r in records if r['production_sample'] == 'data']
        trigger = {'mb': 'minimumbias', 'jet60': 'jet60', 'jet80': 'jet80', 'jet100': 'jet100'}[sample]
        if not data_records or not all(Path(r['csv']).name.startswith(trigger + '_') for r in data_records):
            continue
        matches.append(path)
    if len(matches) != 1:
        raise ValueError(f'{directory}: expected one total for {interval}, {sample}, '
                         f'{direction}/{selection}/{binning}, found {len(matches)}. '
                         'Run systematics_total with matching inputs or select DATA_SYSTEMATICS_CONFIGURATION.')
    manifest = matches[0]
    stem = str(manifest).removesuffix('_provenance.json')
    csv = Path(stem + '_total.csv')
    values = np.loadtxt(csv, delimiter=',', ndmin=2)
    audit = np.genfromtxt(stem + '.csv', delimiter=',', names=True)
    if (values.shape[1] != 4 or not np.isfinite(values).all() or
            not np.allclose(values[:, 1], 1, rtol=0, atol=1e-12) or
            np.any(values[:, 2] <= 0) or np.any(values[:, 3] < 0)):
        raise ValueError(f'{csv}: expected finite center, unity, half-width, fractional total')
    edges = np.column_stack((values[:, 0]-values[:, 2], values[:, 0]+values[:, 2]))
    audit_values = np.column_stack((audit['x_low'], audit['x_high'], audit['total']))
    if (len(audit_values) != len(values) or
            not np.allclose(np.column_stack((edges, values[:, 3])), audit_values, rtol=1e-10, atol=1e-10) or
            np.any(np.diff(edges[:, 0]) <= 0) or np.any(edges[1:, 0] < edges[:-1, 1]-1e-9)):
        raise ValueError(f'{csv}: audit/total CSV alignment mismatch')
    return dict(edges=edges, fraction=values[:, 3], csv=csv, provenance=manifest,
                eta_cut=eta_cut, observable=observable)


def accepted_eta_range(histogram, eta_cut, *, observable='full'):
    """Whole bins intersecting the jet-selected dijet eta support.

    The producer cuts each jet value and fills their mean, not a bin-center
    selection. Keep the partially populated boundary bin; do not include the
    next bin when the cut coincides with an edge (up to edge roundoff).
    """
    if observable not in {'full', 'fb'} or not math.isfinite(eta_cut) or eta_cut <= 0:
        raise ValueError('Expected full/fb and a finite positive eta cut')
    low, high = (-eta_cut, eta_cut) if observable == 'full' else (0., eta_cut)
    axis = histogram.GetXaxis()
    bins = [i for i in range(1, histogram.GetNbinsX()+1)
            if min(axis.GetBinUpEdge(i), high)-max(axis.GetBinLowEdge(i), low) > 1e-9]
    if not bins:
        raise ValueError(f'{histogram.GetName()}: no bins intersect eta acceptance')
    return float(axis.GetBinLowEdge(bins[0])), float(axis.GetBinUpEdge(bins[-1]))


def check_distinct_generator_inputs(curves, official_gen_label):
    """Reject aliases of the same ROOT input before showing duplicate curves."""
    if not any(label == official_gen_label for label, _, _ in curves):
        return
    official = next((Path(path).resolve(), key) for label, path, key in curves
                    if label == official_gen_label)
    for label, path, key in curves:
        if label == official_gen_label or path is None or label == 'Data':
            continue
        if (Path(path).resolve(), key) == official:
            raise ValueError(f'{label} aliases the official generator input {official}; '
                             'disable the optional duplicate or configure a distinct source')


def data_total_band(histogram, total, *, name):
    """Scale stored fractions by plotted y (also valid for independent Data/MC ratios)."""
    delta = histogram.Clone(name + '_fraction')
    delta.SetDirectory(0)
    delta.Reset()
    axis = histogram.GetXaxis()
    bins = []
    acceptance = accepted_eta_range(histogram, total['eta_cut'], observable=total['observable'])
    for i in range(1, histogram.GetNbinsX()+1):
        left, right = axis.GetBinLowEdge(i), axis.GetBinUpEdge(i)
        y = histogram.GetBinContent(i)
        if min(right, acceptance[1])-max(left, acceptance[0]) <= 1e-9 or not math.isfinite(y) or y == 0:
            continue
        found = np.flatnonzero(np.all(np.isclose(total['edges'], (left, right), rtol=0, atol=1e-9), axis=1))
        if len(found) != 1:
            raise ValueError(f'{total["csv"]}: no unique systematic bin for Data edges {(left, right)}')
        j = int(found[0]); bins.append(i)
        delta.SetBinContent(i, float(total['fraction'][j]))
    graph = build_systematic_band(delta, evaluation_range=acceptance, name=name)
    # The shared helper clips edge rectangles to the fiducial interval.
    # Remove absent Data bins rather than turning them into artificial zeros.
    for p in reversed(range(graph.GetN())):
        i = axis.FindBin(graph.GetPointX(p))
        if i not in bins:
            graph.RemovePoint(p)
            continue
        y = histogram.GetBinContent(i)
        fraction = delta.GetBinContent(i)
        graph.SetPoint(p, graph.GetPointX(p), y)
        graph.SetPointError(p, graph.GetErrorXlow(p), graph.GetErrorXhigh(p),
                            abs(y)*fraction, abs(y)*fraction)
    graph.SetFillStyle(0)
    graph.SetMarkerSize(0.)
    return graph


def selection_lines(project_root, interval, *, eta_cut, direction, selection):
    """Read the producer contract instead of duplicating its momentum thresholds."""
    contract = (Path(project_root) / 'AnalysisContract.h').read_text()
    def constant(name):
        match = re.search(r'constexpr float ' + name + r'\s*=\s*([0-9.]+)f?;', contract)
        if not match:
            raise ValueError(f'Missing producer selection constant {name}')
        return float(match[1])
    lead, sub = constant('kLeadingJetPtMin'), constant('kSubleadingJetPtMin')
    angle = constant('kTwoPiOverThree')
    angular = '2#pi/3' if math.isclose(angle, 2*math.pi/3, abs_tol=1e-6) else f'{angle:g}'
    # The producer explicitly reads this anti-kT R=0.4 PF-jet collection.
    producer = (Path(project_root) / 'processForestSimple.C').read_text()
    if 'ak4PFJetAnalyzer/t' not in producer:
        raise ValueError('Cannot verify anti-kT R=0.4 collection')
    low, high = interval
    return ('#it{anti}-k_{T} R = 0.4',
            f'p_{{T,1}} #geq {lead:g}, p_{{T,2}} #geq {sub:g} GeV',
            f'|#eta_{{CM}}^{{jet}}| #leq {eta_cut:g}; |#Delta#phi_{{1,2}}| #geq {angular}',
            f'{low:g} #leq p_{{T}}^{{ave}} < {high:g} GeV')


def decorate_pad(pad, histograms, *, totals=(), annotations=(), labels=None,
                 annotation_size=.030):
    """Add stored systematic caps and selection text to prepared sample objects."""
    import ROOT
    pad.cd()
    retained = []
    primitives = pad.GetListOfPrimitives()
    # Remove old annotation blocks; the complete selection block replaces them.
    if annotations:
        for item in list(primitives):
            if item.InheritsFrom('TLatex') and 'CMS' not in item.GetTitle() and '8.16' not in item.GetTitle():
                primitives.Remove(item)
        retained.extend(draw_text_block(pad, annotations, x=pad.GetLeftMargin()+.025,
                                         y=1-pad.GetTopMargin()-.035,
                                         text_size=annotation_size, line_spacing=.050))
    legends = [p for p in primitives if p.InheritsFrom('TLegend')]
    for leg in legends:
        if labels:
            for entry in leg.GetListOfPrimitives():
                entry.SetLabel(labels(entry.GetLabel()))
        leg.SetTextSize(.028)
    bands = []
    for key, total in totals:
        graph = data_total_band(histograms[key], total, name='g_data_total_' + histograms[key].GetName())
        color = histograms[key].GetMarkerColor()
        graph.SetLineColor(color)
        graph.SetMarkerColor(color)
        graph.SetLineWidth(2)
        graph.Draw('PX SAME')  # retain numeric graph for MPLHEP; ROOT markers/errors are invisible
        retained.extend(draw_root_systematic_brackets(pad, graph))
        bands.append(graph)
    if bands:
        for key, h in histograms.items():
            h.Draw('E1 SAME')
    retained.extend(bands)
    pad.Modified(); pad.Update(); pad.RedrawAxis()
    return retained


def set_acceptance_frame(pad, histogram, x_range):
    """Use exact fiducial edges without changing histogram binning or contents."""
    pad.cd()
    pad.Update()
    primitives = pad.GetListOfPrimitives()
    objects = [(obj, '') for obj in primitives
               if not obj.InheritsFrom('TFrame')]
    for obj, _ in objects:
        primitives.Remove(obj)
    import ROOT
    frame = ROOT.TH1F('hframe_' + pad.GetName(), '', 1, *x_range)
    frame.SetDirectory(0)
    frame.SetMinimum(pad.GetUymin())
    frame.SetMaximum(pad.GetUymax())
    frame.Draw('AXIS')
    for source, target in ((histogram.GetXaxis(), frame.GetXaxis()),
                           (histogram.GetYaxis(), frame.GetYaxis())):
        target.SetTitle(source.GetTitle())
        for attr in ('TitleSize', 'TitleOffset', 'LabelSize', 'LabelOffset',
                     'TitleFont', 'LabelFont', 'Ndivisions', 'TickLength'):
            getattr(target, 'Set' + attr)(getattr(source, 'Get' + attr)())
        target.CenterTitle(source.GetCenterTitle())
    for obj, option in objects:
        if obj.InheritsFrom('TH1'):
            option = 'E1 SAME'
        primitives.Add(obj, option)
    pad.Modified()
    pad.Update()
    return frame


def check_distinct_generator_predictions(histograms, official_gen_label, optional_labels):
    """Require a decision rather than draw an exact duplicate under a new label."""
    if official_gen_label not in histograms:
        return
    official = histograms[official_gen_label]
    def signature(histogram):
        return (tuple(histogram.GetXaxis().GetBinLowEdge(i)
                      for i in range(1, histogram.GetNbinsX()+2)),
                tuple((histogram.GetBinContent(i), histogram.GetBinError(i))
                      for i in range(histogram.GetNbinsX()+2)))
    reference = signature(official)
    for label in optional_labels:
        if label in histograms and signature(histograms[label]) == reference:
            raise ValueError(f'{label} is numerically identical to official MC gen; '
                             'choose which representation to retain before plotting')


def load_comparison_source(filename, key, cache):
    """Load each detached TH2 once per notebook configuration/kernel run."""
    from hist_analysis.python.histogram_io import load_histogram
    identity = (str(Path(filename).resolve()), key)
    if identity not in cache:
        source = load_histogram(*identity)
        if not source.InheritsFrom('TH2'):
            raise TypeError(f'{key!r} in {filename} is not a TH2')
        cache[identity] = source
    return cache[identity]


BRACKET_WIDTH_POINTS = 7.0
BRACKET_END_FRACTION = 0.55


def draw_root_systematic_brackets(pad, graph):
    """Two open brackets at y±syst; lengths are fixed in display points."""
    from array import array
    import ROOT
    from hist_analysis.python.mplhep_plotting import DEFAULT_MPLHEP_STYLE
    pad.cd(); pad.Update()
    # ROOT canvas pixels map to the shared 8-inch figure, independent of y scale.
    width = BRACKET_WIDTH_POINTS / (72 * DEFAULT_MPLHEP_STYLE.figure_size)
    dx = width * (pad.GetUxmax()-pad.GetUxmin()) / (1-pad.GetLeftMargin()-pad.GetRightMargin()) / 2
    dy = width * BRACKET_END_FRACTION * (pad.GetUymax()-pad.GetUymin()) / (1-pad.GetTopMargin()-pad.GetBottomMargin()) / 2
    retained = []
    for i in range(graph.GetN()):
        x, y = graph.GetPointX(i), graph.GetPointY(i)
        for level, sign in ((y-graph.GetErrorYlow(i), 1), (y+graph.GetErrorYhigh(i), -1)):
            if pad.GetLogy() and level <= 0:
                continue
            v = math.log10(level) if pad.GetLogy() else level
            inner = v + sign*dy
            inner = 10**inner if pad.GetLogy() else inner
            line = ROOT.TPolyLine(4, array('d', [x-dx,x-dx,x+dx,x+dx]),
                                 array('d', [inner,level,level,inner]))
            line.SetLineColor(graph.GetLineColor()); line.SetLineWidth(2)
            line.Draw('SAME'); retained.append(line)
    return retained


def draw_mpl_systematic_brackets(axis, arrays, color):
    """Fixed-point open marker paths follow the data transforms during layout."""
    from matplotlib.path import Path as MplPath
    for level, sign in (('ylow', 1), ('yhigh', -1)):
        path = MplPath([(-1, sign*BRACKET_END_FRACTION), (-1,0), (1,0),
                        (1,sign*BRACKET_END_FRACTION)])
        y = np.asarray(arrays['y']) + (-1 if sign == 1 else 1)*np.asarray(arrays[level])
        axis.plot(arrays['x'], y, linestyle='none', marker=path,
                  markersize=BRACKET_WIDTH_POINTS, markerfacecolor='none',
                  markeredgecolor=color, markeredgewidth=1.3, zorder=3, label='_nolegend_')


def data_over_mc(data, mc, *, name, style):
    """ROOT's existing independent-error division, with the MC visual identity."""
    from hist_analysis.python.root_style import apply_curve_style
    result = data.Clone(name); result.SetDirectory(0)
    result.Divide(mc)
    apply_curve_style(result, *style)
    # Keep ROOT numerical zero-denominator behavior in saved products. Display
    # copies omit undefined bins instead of presenting those zeros as measurements.
    result._valid_ratio_bins = tuple(i for i in range(1, result.GetNbinsX()+1)
                                    if all(math.isfinite(v) for v in (mc.GetBinContent(i), mc.GetBinError(i),
                                        data.GetBinContent(i), data.GetBinError(i)))
                                    and mc.GetBinContent(i) != 0)
    return result


def draw_comparison_figure(histograms, *, name, x_title, y_title, x_range,
                           styles, annotations, totals=(), y_range=None, log_y=False,
                           reference=False, bands=(), grid=True, output=None, save_png=False,
                           labels=None):
    """One single-panel ROOT description, optionally exported by shared MPLHEP."""
    from hist_analysis.python.plotting import draw_overlay
    from hist_analysis.python.root_style import save_canvas
    import ROOT
    copies = {}
    for label, h in histograms.items():
        copy = h.Clone(name + '_' + h.GetName()); copy.SetDirectory(0)
        valid = getattr(h, '_valid_ratio_bins', None)
        if valid is not None:
            for i in range(1, copy.GetNbinsX()+1):
                if i not in valid:
                    copy.SetBinContent(i, float('nan')); copy.SetBinError(i, 0.)
        copies[label] = copy
    # Derive visible headroom from actual points/errors, including the stored
    # Data systematics. The reference window is a starting range, never a clip.
    values = [h.GetBinContent(i) for h in copies.values() for i in range(1,h.GetNbinsX()+1)
              if x_range[0] <= h.GetBinCenter(i) <= x_range[1] and math.isfinite(h.GetBinContent(i))]
    if reference and not bands and y_range is None:
        intervals = [(h.GetBinContent(i)-h.GetBinError(i), h.GetBinContent(i)+h.GetBinError(i))
                     for h in copies.values() for i in range(1,h.GetNbinsX()+1)
                     if x_range[0] <= h.GetBinCenter(i) <= x_range[1] and math.isfinite(h.GetBinContent(i))]
        for key, total in totals:
            g = data_total_band(copies[key], total, name='range_only')
            intervals += [(g.GetPointY(i)-g.GetErrorYlow(i),g.GetPointY(i)+g.GetErrorYhigh(i)) for i in range(g.GetN())]
        lo = min(1., *(v[0] for v in intervals)) if intervals else .75
        hi = max(1., *(v[1] for v in intervals)) if intervals else 1.25
        y_range = (lo, hi + 1.6*max(hi-lo, .1))
    elif y_range is None and not log_y:
        y_range = (0., max(values, default=1.) * 2.8)
    elif y_range is None:
        positive = [y for y in values if y > 0]
        if not positive:
            raise ValueError('Logarithmic overlay has no positive points')
        y_range = (min(positive)*.3, max(positive)*5000.)
    # ROOT cannot establish axes from a TH1 containing NaNs, even with explicit
    # limits. Establish the finite frame first, then omit invalid ratio bins on
    # display copies. Numerical exports keep ROOT's original division result.
    missing = [(h, i) for h in copies.values() for i in range(1, h.GetNbinsX()+1)
               if not math.isfinite(h.GetBinContent(i))]
    for h, i in missing:
        h.SetBinContent(i, 0.)
    canvas = draw_overlay(copies, title='', x_title=x_title, y_title=y_title,
        grid=grid, log_y=log_y, output=None, headroom=3000 if log_y else 2.8,
        canvas_name=name, annotations=(), x_range=x_range, y_range=y_range,
        legend_bounds=(.43,.48,.93,.68), curve_styles=styles,
        reference_y=None)
    canvas._data_mc_plot = True
    canvas._marker_only_legend = True  # applies to every notebook plot family
    canvas._comparison_reference = reference
    frame = set_acceptance_frame(canvas, next(iter(copies.values())), x_range) if not log_y else None
    for h, i in missing:
        h.SetBinContent(i, float('nan'))
    extras = []
    for graph, label in bands:
        graph.Draw('2 SAME')
        legend = next(p for p in canvas.GetListOfPrimitives() if p.InheritsFrom('TLegend'))
        legend.AddEntry(graph, label, 'f'); extras.append(graph)
    if reference:
        line = ROOT.TLine(x_range[0],1.,x_range[1],1.)
        line.SetLineStyle(2); line.SetLineColor(ROOT.kGray+1)
        line.Draw('SAME'); extras.append(line)
    extras += decorate_pad(canvas, copies, totals=totals, annotations=annotations, labels=labels)
    for h in copies.values():
        h.Draw('E1 SAME')
    legend = next(p for p in canvas.GetListOfPrimitives() if p.InheritsFrom('TLegend'))
    legend.SetFillStyle(0)
    for entry in legend.GetListOfPrimitives():
        obj = entry.GetObject()
        if obj and obj.InheritsFrom('TH1'):
            entry.SetOption('p')
    # Align the first legend row with the top-left selection's first line.
    # Compact rows avoid moving the lower entries into the distribution.
    legend_top = 1 - canvas.GetTopMargin() - .025
    legend.SetX1NDC(.56); legend.SetX2NDC(.94)
    legend.SetMargin(.16)
    legend.SetY2NDC(legend_top)
    legend.SetY1NDC(legend_top - .05 * legend.GetNRows())
    canvas._comparison_figure_objects = [frame, *copies.values(), *extras]
    canvas.Modified(); canvas.Update()
    save_canvas(canvas, output, save_png=save_png)
    return canvas


def prepare_configured_ratios(histograms, combinations, *, styles, name):
    """Divide compatible prepared shapes with existing independent ROOT errors.

    Missing optional operands are diagnosed, never replaced. Reco/gen uses the
    same independent-error convention as the existing closure comparisons;
    no unavailable event-level covariance is inferred.
    """
    import warnings
    results = {}
    for numerator, denominator in combinations:
        if 'epps21_ct18' in (numerator, denominator):
            raise ValueError('EPPS21 cannot be a Data/MC ratio operand')
        missing = {numerator, denominator} - histograms.keys()
        if missing:
            warnings.warn(f'Skipping {numerator}/{denominator}: unavailable {sorted(missing)}')
            continue
        a, b = histograms[numerator], histograms[denominator]
        edges = lambda h: np.array([h.GetXaxis().GetBinLowEdge(i)
                                   for i in range(1, h.GetNbinsX()+2)])
        ae, be = edges(a), edges(b)
        if ae.shape != be.shape or not np.allclose(ae, be, rtol=0, atol=1e-9):
            raise ValueError(f'Incompatible binning: {numerator}/{denominator}')
        results[numerator, denominator] = data_over_mc(a, b,
            name=f'{name}_{numerator}_over_{denominator}',
            style=styles[denominator if numerator == 'data' else numerator])
    return results
