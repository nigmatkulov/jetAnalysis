"""Strict adapters for the existing headerless systematic TGraphErrors CSVs.

Legacy files do not contain complete production metadata. The caller supplies
an explicit manifest for those facts; this module cannot certify their origin.
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from hist_analysis.config.histograms import pt_bins


COMMON_FIELDS = ('observable', 'frame', 'direction', 'selection', 'eta_cut',
                 'units', 'convention', 'target_sample')
REQUIRED_FIELDS = COMMON_FIELDS + ('binning', 'production_sample', 'extraction')


def archive_systematic_inputs(provenance, output_directory, *, filename_prefix=""):
    """Copy the exact consumed CSVs and record their relative archive paths."""
    output_directory = Path(output_directory)
    archived = provenance.copy()
    paths = []
    for record in provenance.to_dict('records'):
        source = Path(record['csv'])
        relative = Path('individual_uncertainties') / record['source'] / (filename_prefix + source.name)
        destination = output_directory / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        if hashlib.sha256(source.read_bytes()).hexdigest() != record['sha256']:
            raise ValueError(f'Input changed after combination: {source}')
        shutil.copyfile(source, destination)
        if hashlib.sha256(destination.read_bytes()).hexdigest() != record['sha256']:
            raise ValueError(f'Archived input hash mismatch: {destination}')
        paths.append(relative.as_posix())
    archived['copied_csv'] = paths
    return archived


def write_total_systematic_csv(table, output):
    """Write headerless center, unity, half-width, fractional total uncertainty."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    values = np.column_stack(((table.x_low + table.x_high) / 2,
                              np.ones(len(table)),
                              (table.x_high - table.x_low) / 2, table.total))
    np.savetxt(output, values, delimiter=',', fmt='%.17g')
    return output


def load_source(source, config, bins, *, binning):
    """Load exactly one requested variant per interval; never infer missing facts.

    ``pattern`` is relative to ``directory`` and accepts low/high fields. All
    other variant choices should be literal filename tags or explicit globs;
    ambiguous globs fail. ``metadata`` records facts verified by the caller
    against the source run, including facts not encoded in the legacy filename.
    """
    directory = Path(config['directory'])
    metadata = config['metadata']
    missing = [key for key in REQUIRED_FIELDS if metadata.get(key) is None]
    if missing:
        raise ValueError(f'{source}: declare verified metadata: {missing}')
    if metadata['binning'] != binning:
        raise ValueError(f'{source}: manifest binning {metadata["binning"]!r} != {binning!r}')
    canonical = set(map(tuple, pt_bins('dijet', binning)))
    if not bins or len(set(map(tuple, bins))) != len(bins):
        raise ValueError('Select nonempty, unique ptAve bins')
    if not set(map(tuple, bins)) <= canonical:
        raise ValueError(f'{source}: requested bins are incompatible with {binning} binning')
    if (metadata['units'], metadata['convention']) != ('fraction', 'symmetric_relative'):
        raise ValueError(f'{source}: adapter requires symmetric fractional uncertainties')
    tables, provenance = {}, []
    for low, high in bins:
        pattern = config['pattern'].format(low=f'{low:g}'.replace('.', 'p'),
                                          high=f'{high:g}'.replace('.', 'p'))
        candidates = sorted(directory.glob(pattern))
        context = f'{source}, ptAve [{low:g}, {high:g}), input {directory / pattern}'
        if not candidates:
            raise FileNotFoundError(f'Missing requested systematic/variant: {context}')
        if len(candidates) != 1:
            raise ValueError(f'Ambiguous requested variant: {context}: {candidates}')
        path = candidates[0].resolve()
        # Filename tags and path facts must agree with the explicit manifest.
        pt_tag = f'_ptave_{low:g}_{high:g}'.replace('.', 'p') + '.csv'
        if not path.name.endswith(pt_tag):
            raise ValueError(f'Filename ptAve mismatch: {context}')
        frame_tag = f'_eta{metadata["frame"]}_{round(10 * metadata["eta_cut"]):g}_'
        if frame_tag not in path.name or f'_{metadata["observable"]}_relative_' not in path.name:
            raise ValueError(f'Observable/frame/acceptance mismatch: {context}')
        if f'_relative_{metadata["selection"]}_' not in path.name:
            raise ValueError(f'Selection mismatch: {context}')
        if any(part in ('normal', 'test') and part != binning for part in path.parts):
            raise ValueError(f'Input directory binning mismatch: {context}')
        if source in ('JER', 'pointing_resolution'):
            prefix = f'{metadata["production_sample"]}_{metadata["direction"]}_'
            if not path.name.startswith(prefix):
                raise ValueError(f'MC sample/direction mismatch: {context}')
        try:
            values = np.loadtxt(path, delimiter=',', ndmin=2)
        except ValueError as error:
            raise ValueError(f'Invalid four-column CSV: {context}: {error}') from error
        if values.shape[1] != 4 or not len(values) or not np.isfinite(values).all():
            raise ValueError(f'Expected finite x, unity, half-width, fraction: {context}')
        x, unity, half_width, uncertainty = values.T
        if (not np.allclose(unity, 1, rtol=0, atol=1e-12)
                or (half_width <= 0).any() or (uncertainty < 0).any()):
            raise ValueError(f'Invalid unity/width/nonnegative symmetric uncertainty: {context}')
        left, right = x - half_width, x + half_width
        if (np.diff(x) <= 0).any() or not np.allclose(left[1:], right[:-1], rtol=0, atol=1e-8):
            raise ValueError(f'Unordered, duplicate, overlapping or gapped x bins: {context}')
        acceptance = metadata['eta_cut']
        expected = (-acceptance, acceptance) if metadata['observable'] == 'full' else (0, acceptance)
        if left[0] > expected[0] + 1e-8 or right[-1] < expected[1] - 1e-8:
            raise ValueError(f'Incomplete eta acceptance: {context}')
        if (right <= expected[0]).any() or (left >= expected[1]).any():
            raise ValueError(f'Bins outside eta acceptance: {context}')
        tables[(low, high)] = pd.DataFrame({'x_low': left, 'x_high': right, source: uncertainty})
        provenance.append(dict(source=source, pt_low=low, pt_high=high, csv=str(path),
                               sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                               variant=config['pattern'], **metadata))
    return tables, provenance


def combine_systematics(config, bins, *, binning):
    """Return per-bin audit tables and manifest, assuming independent sources."""
    enabled = {name: entry for name, entry in config.items() if entry.get('enabled', True)}
    if not enabled:
        raise ValueError('Enable at least one systematic source')
    if {'total', 'x_low', 'x_high'} & enabled.keys():
        raise ValueError('Systematic name collides with an audit column')
    loaded, provenance = {}, []
    for source, entry in enabled.items():
        loaded[source], records = load_source(source, entry, bins, binning=binning)
        provenance.extend(records)
    manifest = pd.DataFrame(provenance)
    if manifest['csv'].duplicated().any():
        raise ValueError('A CSV was counted more than once (source alias or overlapping bins)')
    for field in COMMON_FIELDS:
        if manifest[field].nunique(dropna=False) != 1:
            raise ValueError(f'Incompatible {field}: {manifest[["source", "csv", field]].to_dict("records")}')
    results = {}
    for interval in bins:
        first = next(iter(loaded))
        table = loaded[first][tuple(interval)].copy()
        for source in list(loaded)[1:]:
            other = loaded[source][tuple(interval)]
            if len(other) != len(table) or not np.allclose(
                    table[['x_low', 'x_high']], other[['x_low', 'x_high']], rtol=0, atol=1e-8):
                raise ValueError(f'{source}, ptAve {interval}: incompatible x bins with {first}; no rebinning')
            table[source] = other[source].to_numpy()
        table['total'] = np.sqrt(np.sum(np.square(table[list(enabled)].to_numpy()), axis=1))
        if not np.isfinite(table['total']).all():
            raise ValueError(f'Nonfinite total in ptAve {interval}')
        results[tuple(interval)] = table
    return results, manifest


def total_systematic_configuration_tag(manifest):
    """Stable identity of declared variants/production roots, independent of run/bin.

    Content hashes are deliberately excluded: rerunning the same configuration
    with refreshed inputs updates its outputs. Full metadata and consumed hashes
    remain in each interval's provenance sidecar.
    """
    ignored = {'csv', 'sha256', 'copied_csv', 'pt_low', 'pt_high'}
    records = set()
    for record in manifest.to_dict('records'):
        identity = {k: v for k, v in record.items() if k not in ignored}
        identity['input_directory'] = str(Path(record['csv']).parent)
        records.add(json.dumps(identity, sort_keys=True, default=str))
    payload = '\n'.join(sorted(records)).encode('utf-8')
    return 'cfg_' + hashlib.sha256(payload).hexdigest()[:12]


def total_systematic_output_stem(manifest, interval):
    """Use existing source filename fields plus a deterministic variant tag."""
    def token(value):
        return re.sub(r'[^A-Za-z0-9]+', '_', str(value)).strip('_')

    def shared(field):
        values = manifest[field].drop_duplicates().tolist() if field in manifest else []
        return token(values[0]) if len(values) == 1 and values[0] is not None else 'mixed'

    sources = '_'.join(token(s) for s in sorted(manifest.source.unique()))
    metadata = manifest.iloc[0]
    low, high = (f'{v:g}'.replace('.', 'p') for v in interval)
    return (
        f"{token(metadata.direction)}_total_{sources}_{token(metadata.selection)}"
        f"_{token(metadata.binning)}_{token(metadata.observable)}"
        f"_fullFit_{shared('full_fit')}_fbFit_{shared('fb_fit')}"
        f"_method_{shared('method')}_extraction_{shared('extraction')}"
        f"_eta{token(metadata.frame)}_{round(10 * metadata.eta_cut):g}"
        f"_{total_systematic_configuration_tag(manifest)}_ptave_{low}_{high}"
    )


def _curve_intersections(ax, bbox):
    """Count displayed curve primitives entering a label rectangle."""
    conflicts = 0
    for line in ax.lines:
        if line.get_linestyle() in ('None', 'none', ' ', ''):
            xy = line.get_transform().transform(line.get_xydata())
            conflicts += sum(bbox.contains(x, y) for x, y in xy)
        elif line.get_path().transformed(line.get_transform()).intersects_bbox(bbox, filled=False):
            conflicts += 1
    for artist in (*ax.collections, *ax.patches):
        paths = artist.get_paths() if hasattr(artist, 'get_paths') else [artist.get_path()]
        for path in paths:
            if path.transformed(artist.get_transform()).intersects_bbox(bbox, filled=False):
                conflicts += 1
    return conflicts


def _place_total_annotations(ax, artists):
    """Keep the shared top-left block unless visible curves occupy that space."""
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    best = None
    for top in (.955, .865, .775, .685, .595):
        for left in (.04, .12, .20, .28, .36, .44):
            for index, artist in enumerate(artists):
                artist.set_position((left, top-index*.045))
            boxes = [artist.get_window_extent(renderer).padded(8) for artist in artists]
            if any(not ax.bbox.contains(box.x0, box.y0) or not ax.bbox.contains(box.x1, box.y1) for box in boxes):
                continue
            score = sum(_curve_intersections(ax, box) for box in boxes)
            if best is None or score < best[0]:
                best = (score, left, top)
            if score == 0:
                break
        if best is not None and best[0] == 0:
            break
    if best is not None:
        for index, artist in enumerate(artists):
            artist.set_position((best[1], best[2]-index*.045))


def _place_total_legend(ax):
    """Prefer the shared upper-right legend; move within the axes if crowded."""
    legend = ax.legend(loc='upper right', fontsize=12, frameon=False)
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    text_boxes = [text.get_window_extent(renderer) for text in ax.texts]
    candidates = [(x, y) for y in np.arange(.98, .25, -.04)
                  for x in np.arange(.98, .41, -.04)]
    best = None
    for anchor in candidates:
        legend.set_bbox_to_anchor(anchor)
        bbox = legend.get_window_extent(renderer).padded(8)
        if not ax.bbox.contains(bbox.x0, bbox.y0) or not ax.bbox.contains(bbox.x1, bbox.y1):
            continue
        score = _curve_intersections(ax, bbox)
        score += 100 * sum(bbox.overlaps(text_box) for text_box in text_boxes)
        if best is None or score < best[0]:
            best = (score, anchor)
        if score == 0:
            break
    if best is not None:
        legend.set_bbox_to_anchor(best[1])
    return legend


def plot_systematics(table, manifest, interval, *, y_range=None,
                     cms_header_left=None, cms_header_right='pPb (8.16 TeV)'):
    """Plot percent display values using the individual-systematic CMS scheme."""
    from hist_analysis.python.mplhep_plotting import (
        _mplhep, _curve_style, _draw_cms_header, DEFAULT_MPLHEP_STYLE,
    )
    if y_range is not None and (
        len(y_range) != 2 or not np.isfinite(y_range).all() or y_range[0] >= y_range[1]
    ):
        raise ValueError('y_range must contain finite values with low < high')
    plt, hep = _mplhep()
    hep.style.use('CMS')
    metadata = manifest.iloc[0]
    samples = list(manifest.production_sample.unique())
    if cms_header_left is None:
        cms_header_left = 'CMS' if 'data' in samples else 'CMS Simulation'
    labels = {'JER': 'JER', 'JEU': 'JES', 'pileup': 'Pileup',
              'pointing_resolution': 'Pointing resolution'}
    source_indices = {'JER': 0, 'JEU': 1, 'pileup': 2, 'pointing_resolution': 3}
    sources = list(manifest.source.unique())
    style = DEFAULT_MPLHEP_STYLE
    fig, ax = plt.subplots(figsize=(style.figure_size, style.figure_size), layout='constrained')
    edges = np.r_[table.x_low.to_numpy(), table.x_high.iloc[-1]]
    for index, source in enumerate(sources):
        label = labels.get(source, source)
        _, color, marker, face = _curve_style(label, source_indices.get(source, index))
        hep.histplot(100 * table[source].to_numpy(), bins=edges, yerr=False, xerr=True,
            histtype='errorbar', ax=ax, color=color, marker=marker, markerfacecolor=face,
            markeredgecolor=color, markeredgewidth=1.5, markersize=style.marker_size,
            capsize=style.capsize, zorder=3, label=label)
    # Total is drawn last and remains visible when components coincide.
    ax.stairs(100 * table.total, edges, color='black', linewidth=2.6, zorder=4)
    ax.plot((table.x_low + table.x_high)/2, 100 * table.total, color='black',
        marker='s', markersize=6, linestyle='none', zorder=5, label='Total')
    x_title = (r'$\eta_{\mathrm{CM}}^{\mathrm{dijet}}$' if metadata.observable == 'full'
               else r'$|\eta_{\mathrm{CM}}^{\mathrm{dijet}}|$')
    ax.set_xlabel(x_title, fontsize=20)
    ax.set_ylabel(r'Relative systematic uncertainty $\delta_{\mathrm{syst}}$ (%)', fontsize=18)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(*(y_range if y_range is not None else
                  (0, max(float(table.total.max()) * 100 * 1.65, 1.0))))
    ax.tick_params(labelsize=16)
    ax.grid(True, alpha=0.2)
    sample_names = {'data': 'Data', 'embedding': 'Embedding MC', 'pythia': 'PYTHIA 8'}
    direction_names = {'combined': 'Combined', 'pgoing': 'p-going', 'Pbgoing': 'Pb-going'}
    selection_names = {'jetId': 'Jet ID', 'trkMax': 'Track maximum selection', 'noSel': 'No jet selection'}
    observable = 'Unit-normalized CM dijet shape' if metadata.observable == 'full' else 'Dijet Forward / Backward'
    annotations = (
        ' + '.join(sample_names.get(s, s) for s in samples),
        direction_names.get(metadata.direction, metadata.direction) + '; ' +
            selection_names.get(metadata.selection, metadata.selection),
        rf'${interval[0]:g} \leq p_{{\mathrm{{T}}}}^{{\mathrm{{ave}}}} < {interval[1]:g}$ GeV',
        rf'$|\eta_{{\mathrm{{CM}}}}^{{\mathrm{{jet}}}}| \leq {metadata.eta_cut:g}$',
        observable,
        'Symmetric; quadrature total',
    )
    annotation_artists = [ax.text(.04, .955-index*.045, text, transform=ax.transAxes,
        ha='left', va='top', fontsize=12) for index, text in enumerate(annotations)]
    _draw_cms_header(hep, ax, left=cms_header_left, right=cms_header_right,
                     style=DEFAULT_MPLHEP_STYLE)
    _place_total_annotations(ax, annotation_artists)
    _place_total_legend(ax)
    return fig
