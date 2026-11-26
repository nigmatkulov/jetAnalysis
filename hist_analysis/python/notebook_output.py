"""One deterministic ROOT archive per notebook, without changing live objects.

Canvas snapshots preserve helper-local display copies at export time. The final
namespace pass also saves unplotted results in ordinary containers/dataclasses.
Only explicitly registered output directories participate; ordinary plotting
calls outside a notebook archive retain their existing behavior.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import fields, is_dataclass
import hashlib
import json
from pathlib import Path
import re
import weakref

import ROOT

from hist_analysis.python.histogram_io import open_root_file
from hist_analysis.python.output_naming import legacy_output_name


_ARCHIVES = weakref.WeakValueDictionary()


def _token(value):
    text = str(value)
    clean = re.sub(r'[^A-Za-z0-9_]+', '_', text).strip('_') or 'object'
    # Distinguish e.g. tuple intervals, labels containing punctuation, and paths.
    if clean != text or len(clean) > 100:
        clean = clean[:100] + '_' + hashlib.sha256(text.encode()).hexdigest()[:12]
    return clean


def _objects(value, path=(), ancestors=frozenset()):
    """Walk supported result containers; never introspect arbitrary Python/ROOT objects."""
    if id(value) in ancestors:
        return
    ancestors = ancestors | {id(value)}
    if isinstance(value, ROOT.TObject):
        if not value:
            return
        if any(value.InheritsFrom(cls) for cls in ('TH1', 'THnBase', 'TGraph')):
            yield path, value
        elif value.InheritsFrom('TPad'):
            # Ignore canvas/pad names: some helpers use Python's randomized hash.
            for index, primitive in enumerate(value.GetListOfPrimitives()):
                yield from _objects(primitive, (*path, f'primitive_{index}'), ancestors)
        elif value.InheritsFrom('THStack'):
            for index, histogram in enumerate(value.GetHists() or ()):
                yield from _objects(histogram, (*path, f'histogram_{index}'), ancestors)
    elif isinstance(value, Mapping):
        for key, item in value.items():
            if not str(key).startswith('_'):
                yield from _objects(item, (*path, _token(key)), ancestors)
    elif isinstance(value, (tuple, list)):
        for index, item in enumerate(value):
            yield from _objects(item, (*path, str(index)), ancestors)
    elif is_dataclass(value) and not isinstance(value, type):
        for field in fields(value):
            yield from _objects(getattr(value, field.name), (*path, field.name), ancestors)


class NotebookRootOutput:
    """Stream native ROOT objects to a fresh file and replace keys on cell reruns.

    Keys use stable result-container paths rather than potentially random object
    names. ``plots`` contains the exact displayed objects; ``results`` contains
    retained numerical products. Legacy explicit writers retain their object
    keys under ``exports``. No input histogram is renamed or re-owned.
    """

    def __init__(self, directory, notebook):
        self.directory = Path(directory).resolve()
        self.path = directory / f'{notebook}.root'
        self.notebook = notebook
        self.exports = {}
        self.manifest = {}
        self.selection = getattr(directory, 'selection', None)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with ROOT.TDirectory.TContext():
            output = open_root_file(str(self.path), 'RECREATE')
            output.Close()
        _ARCHIVES[str(self.directory)] = self

    def _write(self, objects):
        with ROOT.TDirectory.TContext():
            output = open_root_file(str(self.path), 'UPDATE')
            try:
                for parts, obj in objects:
                    target = output
                    for part in parts[:-1]:
                        target = target.GetDirectory(part) or target.mkdir(part)
                        if not target:
                            raise OSError(f'Cannot create ROOT directory: {parts}')
                    name = parts[-1]
                    target.cd()
                    # Clone only one object at a time; preserve Sumw2, flows,
                    # variable bins and profile state through ROOT streaming.
                    # Keep explicit legacy export names; newly collected
                    # histograms also carry their full semantic container path
                    # as a name, rather than an ambiguous sequence index.
                    clone_name = name if parts[0] in {'exports', 'production'} else '__'.join(parts)
                    clone = obj.Clone(clone_name)
                    if hasattr(clone, 'SetDirectory'):
                        clone.SetDirectory(0)
                    if clone.Write(name, ROOT.TObject.kOverwrite) <= 0:
                        raise OSError(f'Cannot write {parts} to {self.path}')
                    self.manifest['/'.join(parts)] = dict(
                        original_name=obj.GetName(), class_name=obj.ClassName())
                if output.TestBit(ROOT.TFile.kWriteError):
                    raise OSError(f'ROOT write failed: {self.path}')
            finally:
                output.Close()

    def collect(self, value, *, prefix='results'):
        self._write(_objects(value, (_token(prefix),)))

    def save_histograms(self, filename, histograms, mode='RECREATE'):
        """Notebook replacement for the existing explicit histogram writer."""
        if mode != 'RECREATE':
            raise ValueError('Notebook archives require RECREATE semantics')
        self.write_unfolding_output(filename, histograms=histograms)

    def write_unfolding_output(self, filename, *, histograms, covariance=None,
                               covariance_name='hUnfoldedCovariance', response=None,
                               response_name='responseEtaCM', metadata=None,
                               metadata_name='unfoldingConfiguration'):
        """Keep legacy keys and specialized unfolding objects in this archive."""
        scope = _token(legacy_output_name(filename).removesuffix('.root'))
        objects = {}
        for histogram in histograms:
            if histogram is not None:
                name = histogram.GetName()
                if name in objects:
                    if objects[name] is histogram:
                        continue
                    suffix = 2
                    while f'{name}_{suffix}' in objects:
                        suffix += 1
                    name = f'{name}_{suffix}'
                objects[name] = histogram
        if covariance is not None:
            objects[covariance_name] = covariance
        if response is not None:
            objects[response_name] = response
        if metadata is not None:
            objects[metadata_name] = ROOT.TNamed(metadata_name, json.dumps(metadata, sort_keys=True))
        self.exports[scope] = str(filename)
        self._write(((('exports', scope, name), obj) for name, obj in objects.items()))

    def import_file(self, filename, *, prefix):
        """Copy producer products without replacing its downstream input catalog."""
        with ROOT.TDirectory.TContext():
            source = open_root_file(str(filename))
            try:
                def walk(directory, parts):
                    for key in directory.GetListOfKeys():
                        obj = key.ReadObj()
                        if obj.InheritsFrom('TDirectory'):
                            yield from walk(obj, (*parts, key.GetName()))
                        elif any(obj.InheritsFrom(cls) for cls in ('TH1', 'THnBase', 'TGraph')):
                            yield (*parts, key.GetName()), obj
                self._write(walk(source, ('production', _token(prefix))))
            finally:
                source.Close()

    def finish(self, namespace):
        # Exclude IPython history and unrelated canvases from previous executions.
        self.collect({key: value for key, value in namespace.items()
                      if not key.startswith('_') and key not in {'In', 'Out'} })
        def configuration(value):
            if value is None or isinstance(value, (str, bool, int, float)):
                return value
            if isinstance(value, Path):
                return str(value)
            if isinstance(value, Mapping):
                return {str(key): configuration(item) for key, item in value.items()}
            if isinstance(value, (tuple, list)):
                return [configuration(item) for item in value]
            return f'<{type(value).__name__}>'
        payload = dict(notebook=self.notebook, schema=1,
                       selection=self.selection.identity() if self.selection else None,
                       configuration={key: configuration(value) for key, value in namespace.items()
                                      if key.isupper() and not key.startswith('_')},
                       exports=self.exports, objects=self.manifest)
        metadata = ROOT.TNamed('notebookConfiguration', json.dumps(payload, sort_keys=True))
        self._write([(('notebookConfiguration',), metadata)])
        print(f'ROOT archive: {self.path} ({len(self.manifest) - 1} objects)')
        return self.path


def capture_canvas(canvas, output):
    """Persist helper-local plot objects only for a registered notebook output."""
    if output is None:
        return
    path = Path(output).resolve()
    for parent in path.parents:
        archive = _ARCHIVES.get(str(parent))
        if archive is not None:
            logical = path.relative_to(parent).with_name(legacy_output_name(path)).with_suffix('')
            archive.collect(canvas, prefix='plots_' + logical.as_posix())
            return


def systematic_table_histograms(table):
    """Represent the total-systematic plot's table as native variable-bin TH1Ds.

    Fractional curves and their percent display counterparts are both explicit.
    The values themselves are uncertainty magnitudes, not statistical errors;
    no Poisson errors are assigned to these deterministic table entries.
    """
    from array import array
    edges = [*map(float, table.x_low), float(table.x_high.iloc[-1])]
    histograms = {}
    for column in table.columns:
        if column in {'x_low', 'x_high'} or not table[column].dtype.kind in 'fiu':
            continue
        for unit, scale in (('fraction', 1.0), ('percent', 100.0)):
            name = f'{column}_{unit}'
            histogram = ROOT.TH1D(name, f'{column};eta;Uncertainty ({unit})',
                                  len(edges) - 1, array('d', edges))
            histogram.SetDirectory(0)
            histogram.Sumw2()
            for index, value in enumerate(table[column], 1):
                histogram.SetBinContent(index, scale * float(value))
                histogram.SetBinError(index, 0.0)
            histograms[name] = histogram
    return histograms
