"""Exact selection identities and non-destructive notebook output layouts.

New runs use dated renderer directories and unique filenames; existing artifacts are read-only.
An output sidecar records the complete bin set and the original descriptive
stem. Readers match that declaration rather than guessing selections from old
rounded filenames. No analysis calculation or histogram key is changed here.
"""
from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from fnmatch import fnmatch
import hashlib
import json
import math
import os
from pathlib import Path
import re
import time
import uuid
from datetime import datetime
from zoneinfo import ZoneInfo


JET_SELECTIONS = frozenset(('noSel', 'trkMax', 'jetId'))


def number_token(value):
    """Lossless decimal spelling of a configured number, safe in filenames."""
    if isinstance(value, bool):
        raise ValueError('Boolean is not a physics boundary')
    number = Decimal(str(value))
    if not number.is_finite():
        raise ValueError(f'Non-finite physics boundary: {value}')
    text = format(number, 'f')
    if '.' in text:
        text = text.rstrip('0').rstrip('.')
    if number == 0:
        text = '0'
    return text.replace('-', 'm').replace('.', 'p')


def _numbers(values):
    return [number_token(value) for value in values]


def _intervals(values):
    if values is None:
        return []
    if isinstance(values, dict):
        return {str(key): _intervals(value) for key, value in values.items()}
    result = []
    for interval in values:
        # Unfolding PtInterval objects are explicit low/high records.
        low, high = (interval.low, interval.high) if hasattr(interval, 'low') else interval
        if Decimal(str(low)) >= Decimal(str(high)):
            raise ValueError(f'Invalid half-open interval: {interval}')
        result.append(_numbers((low, high)))
    return result


def _flat_intervals(values):
    if isinstance(values, dict):
        return [interval for group in values.values() for interval in _flat_intervals(group)]
    return values


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


@dataclass(frozen=True)
class OutputSelection:
    selection: str | tuple[str, ...] | None
    eta_cuts: object = None
    pt_bins: object = None
    pt_kind: str = 'ptAve'

    def identity(self):
        roles = self.selection if isinstance(self.selection, dict) else {}
        selections = tuple(dict.fromkeys(roles.values())) if roles else ((self.selection,) if isinstance(self.selection, str) else tuple(self.selection or ()))
        if any(selection not in JET_SELECTIONS for selection in selections):
            raise ValueError(f'Unsupported jet selection: {selections}')
        cuts = self.eta_cuts
        if cuts is None:
            eta = {}
        elif isinstance(cuts, dict):
            eta = {str(frame): _numbers(values if isinstance(values, (list,tuple)) else (values,))
                   for frame, values in cuts.items()}
        else:
            eta = {'CM': _numbers(cuts if isinstance(cuts, (list,tuple)) else (cuts,))}
        identity = dict(selection=list(selections), eta_cuts=eta,
                        pt_bins=_intervals(self.pt_bins), pt_kind=self.pt_kind)
        if roles:
            identity['selection_roles'] = roles
        return identity

    def token(self):
        identity = self.identity()
        selection = '_'.join(identity['selection']) or 'jetSelectionNA'
        if len(identity['selection']) > 1 and 'selection_roles' in identity:
            selection = '_'.join(role+'-'+method for role,method in identity['selection_roles'].items())
        eta = '_'.join('eta'+frame+'_'.join(values) for frame, values in identity['eta_cuts'].items()) or 'etaCutNA'
        intervals = _flat_intervals(identity['pt_bins'])
        if not intervals:
            pt = self.pt_kind+'NA'
        elif len(intervals) == 1:
            pt = self.pt_kind+intervals[0][0]+'to'+intervals[0][1]
        else:
            # The manifest distinguishes gaps, overlaps, ordering and trigger maps,
            # not just endpoints. Full precision participates in the identifier.
            pt = self.pt_kind+'Bins'+_digest(identity['pt_bins'])[:16]
        return f'{selection}_{eta}_{pt}'


def output_filename(observable, comparison, *, selection, eta_cut, pt_interval=None,
                    pt_bins=None, extension='pdf'):
    """Public single-interval / complete-bin-set filename builder."""
    if pt_interval is not None and pt_bins is not None:
        raise ValueError('Specify an interval or a bin set, not both')
    context = OutputSelection(selection, eta_cut, [pt_interval] if pt_interval is not None else pt_bins)
    return _name(f'{observable}_{comparison}.{extension.lstrip(".")}', context)


def _name(legacy, selection):
    path = Path(legacy)
    stem = selection.token()+'_'+re.sub(r'_+', '_', path.stem)
    # Long systematic method stems already approach NAME_MAX. Preserve their
    # exact original names in the sidecar and include a checked digest in the
    # shortened filename; never truncate the physics identity.
    if len((stem+path.suffix).encode()) > 220:
        marker = '_name'+_digest(path.stem)[:20]
        prefix = selection.token()+'_'
        budget = 220-len((prefix+marker+path.suffix).encode())
        if budget < 0:
            raise ValueError('Selection identity exceeds safe filename length; split the output bin set')
        readable = path.stem.encode()[:min(55,budget)].decode('utf-8','ignore')
        stem = prefix+readable+marker
    if len((stem+path.suffix).encode()) > 220:
        raise ValueError('Selection identity exceeds safe filename length; split the output bin set')
    if '/' in stem or '\\' in stem:
        raise ValueError('Unsafe filename component in output selection')
    return re.sub(r'_+', '_', stem+path.suffix)


def output_date():
    return datetime.now(ZoneInfo('America/Chicago')).date().isoformat()


def notebook_output_base(path):
    """Workflow root, accepting legacy selection/direction suffixes."""
    path = Path(path)
    parts = path.parts
    for index in range(len(parts)-2):
        if re.fullmatch(r'\d{4}-\d{2}-\d{2}', parts[index]) and parts[index+1] in JET_SELECTIONS and parts[index+2] in {'root','mplhep'}:
            return Path(*parts[:index])
    for index in range(len(parts)-2):
        if parts[index:index+2] == ('hist_analysis', 'output'):
            return Path(*parts[:index+3])
    return path


def mplhep_output_path(path):
    path = Path(path)
    parts = list(path.parts)
    for index in range(2, len(parts)):
        if parts[index] == 'root' and parts[index-1] in JET_SELECTIONS and re.fullmatch(r'\d{4}-\d{2}-\d{2}', parts[index-2]):
            parts[index] = 'mplhep'
            return Path(*parts)
    return path.parent / 'mplhep' / path.name


class OutputDirectory(os.PathLike):
    """A directory whose file joins include an explicitly declared selection.

    Pass it through ``as_output_directory`` in helpers accepting output dirs.
    Files returned by a join are ordinary Paths; no global naming state exists.
    """
    def __init__(self, root, selection, *, run_id=None, new_run=True, directory_selection=None, workflow_root=None):
        self.selection = selection
        identity = selection.identity()
        self.run_id = run_id or f'{time.time_ns()}_{uuid.uuid4().hex[:12]}'
        root = Path(root)
        self.namespace = ''
        if new_run:
            methods = identity['selection']
            method = directory_selection or identity.get('selection_roles', {}).get('data')
            if method is None and len(methods) == 1:
                method = methods[0]
            if method not in JET_SELECTIONS:
                raise ValueError('Declare directory_selection for generator or selection-comparison outputs')
            base = Path(workflow_root) if workflow_root is not None else notebook_output_base(root)
            self.namespace = '_'.join(root.relative_to(base).parts)
            self.path = base / output_date() / method / 'root'
            (self.path.parent / 'mplhep').mkdir(parents=True, exist_ok=True)
        else:
            self.path = root
        self.path.mkdir(parents=True, exist_ok=True)

    def __fspath__(self):
        return str(self.path)

    def __str__(self):
        return str(self.path)

    def __truediv__(self, component):
        component = Path(component)
        if component.is_absolute() or '..' in component.parts:
            raise ValueError('Output joins must stay beneath the declared directory')
        if component.suffix.lower() not in {'.pdf','.png','.root','.csv','.json','.npz','.npy','.txt','.html','.svg'}:
            target = mplhep_output_path(self.path / 'placeholder.pdf').parent if str(component) == 'mplhep' else self.path/component
            child = OutputDirectory(target, self.selection, run_id=self.run_id, new_run=False)
            child.namespace = self.namespace
            return child
        return self.file(component)

    def file(self, legacy, *, selection=None):
        legacy = Path(legacy)
        if legacy.is_absolute() or '..' in legacy.parts:
            raise ValueError('Output joins must stay beneath the declared directory')
        context = selection or self._file_selection(legacy)
        descriptive = '_'.join(part for part in (self.namespace, legacy.stem, 'run'+self.run_id) if part)
        path = self.path/legacy.parent/_name(descriptive+legacy.suffix, context)
        path.parent.mkdir(parents=True, exist_ok=True)
        record = dict(schema=1, identity=context.identity(), legacy_name=legacy.name,
                      run_id=self.run_id, filename=path.name)
        sidecar = selection_sidecar(path)
        if sidecar.exists():
            if json.loads(sidecar.read_text()) != record:
                raise ValueError(f'Output identity collision: {path}')
        else:
            with sidecar.open('x') as stream:
                json.dump(record, stream, sort_keys=True, indent=2)
                stream.write('\n')
        return path

    def _file_selection(self, legacy):
        # Only current explicitly declared configuration is consulted. This does
        # not recover or guess selections of any pre-existing file.
        if legacy.suffix.lower() == '.root':
            return self.selection  # Bundles may contain many pT intervals.
        selection = self.selection
        if isinstance(selection.pt_bins, dict) and {'ptAve','jetPt'} <= selection.pt_bins.keys():
            inclusive = bool(re.search(r'inclusive|single', legacy.name, re.I))
            kind = 'jetPt' if inclusive else 'ptAve'
            eta = ({k:v for k,v in selection.eta_cuts.items() if k != 'CM'} if inclusive and isinstance(selection.eta_cuts, dict) else
                   None if inclusive else selection.eta_cuts)
            selection = OutputSelection(selection.selection, eta, selection.pt_bins[kind], kind)
        identity = selection.identity()
        eta_cuts = selection.eta_cuts
        frame = re.search(r'_eta(CM|Lab)_', legacy.name)
        if frame and isinstance(eta_cuts, dict) and frame.group(1) in eta_cuts:
            eta_cuts = {k:v for k,v in eta_cuts.items() if k not in {'CM','Lab'} or k == frame.group(1)}
        cut_match = re.search(r'_eta(CM|Lab)_([m-]?[0-9]+(?:[p.][0-9]+)?)(?=_|\.[A-Za-z]|$)',legacy.name)
        if cut_match:
            frame_name, spelling = cut_match.groups()
            declared_cuts = identity['eta_cuts'].get(frame_name, [])
            matches = []
            for token in declared_cuts:
                value = Decimal(token.replace('m','-').replace('p','.'))
                if spelling in {token, str(value), str(round(10*float(value)))}:
                    matches.append(value)
            if len(matches)>1:
                raise ValueError(f'Ambiguous rounded eta cut in {legacy}; use explicit .file selection')
            if len(matches)==1:
                eta_cuts = dict(eta_cuts) if isinstance(eta_cuts,dict) else {}
                eta_cuts[frame_name] = matches
        candidates = []
        for low, high in _flat_intervals(identity['pt_bins']):
            values = [token.replace('m','-').replace('p','.') for token in (low,high)]
            spellings = [(a,b) for a in (values[0], format(float(values[0]),'g'))
                               for b in (values[1], format(float(values[1]),'g'))]
            if any(re.search(r'(?:ptave_|pt_|ptAve)'+re.escape(a)+r'(?:_|to)'+re.escape(b)+r'(?:_|\.|$)',legacy.name,re.I)
                   for a,b in spellings):
                candidates.append((low,high))
        candidates = list(dict.fromkeys(candidates))
        if len(candidates)>1:
            raise ValueError(f'Ambiguous rounded interval in {legacy}; use .file(..., selection=OutputSelection(...))')
        if not candidates:
            return OutputSelection(selection.selection, eta_cuts, selection.pt_bins, selection.pt_kind)
        low, high = (Decimal(v.replace('m','-').replace('p','.')) for v in candidates[0])
        return OutputSelection(selection.selection, eta_cuts, [(low, high)], selection.pt_kind)

    def joinpath(self, *parts):
        result = self
        for part in parts:
            result = result/part
        return result

    def __getattr__(self, attribute):
        return getattr(self.path, attribute)


def as_output_directory(directory):
    return directory if isinstance(directory, OutputDirectory) else Path(directory)


def selection_sidecar(path):
    path = Path(path)
    return path.with_name(path.name+'.selection.json')


def output_record(path):
    sidecar = selection_sidecar(path)
    return json.loads(sidecar.read_text()) if sidecar.is_file() else None


def legacy_output_name(path):
    record = output_record(path)
    return record['legacy_name'] if record else Path(path).name


def find_outputs(directory, legacy_pattern, *, selection=None, eta_cut=None, interval=None, frame='CM'):
    """Prefer declared outputs; choose the newest matching run per logical name.

    Distinct legacy names remain ambiguous rather than being merged. Existing
    legacy readers perform their original strict provenance/format checks when
    no declared output is available.
    """
    directory = Path(directory)
    search_root = notebook_output_base(directory)
    # Legacy configured subdirectories may live beneath a custom workflow root.
    # Discover only explicitly dated layouts; physics is still checked below.
    if search_root == directory:
        for parent in (directory, *directory.parents):
            if parent == Path(parent.anchor):
                break
            if parent.is_dir() and any(child.is_dir() and re.fullmatch(r'\d{4}-\d{2}-\d{2}', child.name) for child in parent.iterdir()):
                search_root = parent
                break
    extension = Path(legacy_pattern).suffix
    found = {}
    for sidecar in search_root.rglob('*'+extension+'.selection.json'):
        record = json.loads(sidecar.read_text())
        path = sidecar.with_name(record['filename'])
        if not path.is_file() or not fnmatch(record['legacy_name'], legacy_pattern):
            continue
        identity = record['identity']
        if selection is not None and identity['selection'] != [selection]:
            continue
        if eta_cut is not None and number_token(eta_cut) not in identity['eta_cuts'].get(frame, []):
            continue
        if interval is not None and _numbers(interval) not in _flat_intervals(identity['pt_bins']):
            continue
        old = found.get(record['legacy_name'])
        if old is None or record['run_id'] > old[0]:
            found[record['legacy_name']] = (record['run_id'], path)
    if found:
        return sorted(value[1] for value in found.values())
    return sorted(path for path in directory.glob(legacy_pattern)
                  if path.is_file() and output_record(path) is None)


def sibling_output(path, old_suffix, new_suffix):
    """Locate a related total/audit artifact in the same declared run."""
    record = output_record(path)
    if record is None:
        return Path(str(path).removesuffix(old_suffix)+new_suffix)
    legacy = record['legacy_name'].removesuffix(old_suffix)+new_suffix
    matches = [p for p in find_outputs(Path(path).parent, legacy)
               if output_record(p)['run_id'] == record['run_id']]
    if len(matches) != 1:
        raise ValueError(f'Missing/ambiguous sibling {legacy} for {path}')
    return matches[0]


def nominal_jet_eta_cut(project_root):
    """Read the nominal producer contract for unindexed dijet/run diagnostics."""
    contract = (Path(project_root)/'AnalysisContract.h').read_text()
    match = re.search(r'kNominalJetEtaMax\s*=\s*([0-9.]+)f', contract)
    if not match:
        raise ValueError('Producer contract has no nominal jet eta acceptance')
    return float(match.group(1))


def related_output(path, suffix, extension='pdf'):
    """Build a diagnostic sibling while preserving its exact selection record."""
    record = output_record(path)
    if record is None:
        path = Path(path)
        return path.with_name(path.stem+suffix+'.'+extension)
    identity = record['identity']
    def value(token):
        return Decimal(token.replace('m','-').replace('p','.'))
    cuts = {frame: [value(v) for v in values] for frame, values in identity['eta_cuts'].items()}
    def bins(values):
        return {k:bins(v) for k,v in values.items()} if isinstance(values,dict) else [[value(v) for v in pair] for pair in values]
    context = OutputSelection(identity.get('selection_roles',tuple(identity['selection'])), cuts, bins(identity['pt_bins']), identity['pt_kind'])
    directory = OutputDirectory(Path(path).parent, context, run_id=record['run_id'], new_run=False)
    return directory.file(Path(record['legacy_name']).stem+suffix+'.'+extension)


def bundle_selection(selection, *, standalone=False):
    from hist_analysis.config.histograms import DIJET_ETA_CUTS, SYSTEMATIC_PTAVE_BINS, TEST_DIJET_PTAVE_BINS
    # Bundles retain the full indexed cut catalog, including standalone gen
    # arrays. The standalone collection's historical jetId stem is a storage
    # namespace; reconstructed jet ID is not applied to pure generator jets.
    return OutputSelection(None if standalone else selection, DIJET_ETA_CUTS,
        {'standard': SYSTEMATIC_PTAVE_BINS, 'test': TEST_DIJET_PTAVE_BINS})


def bundle_filename(legacy, selection, *, standalone=False):
    path = Path(legacy)
    return path.stem+'_'+bundle_selection(selection, standalone=standalone).token()+path.suffix


def unused_version(path):
    """Allocate a destination without renaming/removing an existing artifact."""
    path = Path(path)
    if not path.exists():
        return path
    index = 2
    while path.with_name(f'{path.stem}_v{index}{path.suffix}').exists():
        index += 1
    return path.with_name(f'{path.stem}_v{index}{path.suffix}')


def verified_reference_eta_cuts(inputs):
    """Identify fixed pp-reference acceptances from declared indexed histograms.

    The producer's unsuffixed pp reference is compared with its indexed cut
    catalog, including contents, errors and axes. No acceptance is assumed from
    a filename or from the consumer's Data cut. Unselected inputs are not read.
    """
    import ROOT
    import numpy as np
    result = {}
    for label, filename, key in inputs:
        if not Path(filename).is_file():
            continue  # Optional-source preflight reports unavailable inputs.
        file = ROOT.TFile.Open(str(filename))
        if not file or file.IsZombie():
            raise OSError(f'Cannot inspect reference acceptance: {filename}')
        try:
            source = file.Get(key)
            if not source:
                continue  # The existing exact-key preflight handles this case.
            metadata = file.Get('histogramProductionMetadata')
            table = json.loads(metadata.GetTitle()).get('eta_cut_index_to_value',{}) if metadata else {}
            if not table:
                # Original producer titles explicitly state each indexed cut.
                for item in file.GetListOfKeys():
                    match = re.fullmatch(r'hGenDijetPtEtaCM_(\d+)',item.GetName())
                    if match:
                        title = file.Get(item.GetName()).GetTitle()
                        cut = re.search(r'\|#eta[^|]*\|\s*(?:<|#leq|<=)\s*([0-9.]+)',title)
                        if cut:
                            table[match.group(1)] = cut.group(1)
            if not table:
                raise ValueError(f'{filename}: reference jet acceptance is undeclared')
            indexed = re.fullmatch(r'hGenDijetPtEtaCM_(\d+)',key)
            if indexed and indexed.group(1) in table:
                result[label] = table[indexed.group(1)]
                continue
            if key != 'hGenDijetPtEtaCM' or not source.InheritsFrom('TH2D'):
                raise ValueError(f'{filename}:{key}: cannot verify the reference acceptance')
            values = np.frombuffer(source.GetArray(),dtype=np.float64,count=source.GetNcells())
            matches = []
            for index,cut in table.items():
                other = file.Get('hGenDijetPtEtaCM_'+str(index))
                if not other or not other.InheritsFrom('TH2D') or other.GetNcells()!=source.GetNcells():
                    continue
                candidate = np.frombuffer(other.GetArray(),dtype=np.float64,count=other.GetNcells())
                if not np.array_equal(values,candidate):
                    continue
                same_axes = all(source.GetXaxis().GetBinLowEdge(i)==other.GetXaxis().GetBinLowEdge(i)
                    for i in range(1,source.GetNbinsX()+2)) and all(
                    source.GetYaxis().GetBinLowEdge(i)==other.GetYaxis().GetBinLowEdge(i)
                    for i in range(1,source.GetNbinsY()+2))
                if same_axes and all(source.GetBinError(i)==other.GetBinError(i) for i in range(source.GetNcells())):
                    matches.append(cut)
            if len(matches)!=1:
                raise ValueError(f'{filename}:{key}: no unique verified jet acceptance ({matches})')
            result[label] = matches[0]
        finally:
            file.Close()
    return result
