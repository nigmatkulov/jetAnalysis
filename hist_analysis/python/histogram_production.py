"""Explicit production and read-only resolution of reusable ROOT histograms."""

from __future__ import annotations

import json
from pathlib import Path
import re
from typing import Iterable

import ROOT

from hist_analysis.config.files import BASE_DIR
from hist_analysis.config.histograms import DIJET_PTAVE_BINS, TEST_DIJET_PTAVE_BINS, DIJET_ETA_CUTS
from hist_analysis.python.histogram_io import (
    open_root_file, resolve_combined_file, resolve_data_file,
    produced_histogram_path, find_produced_histogram_file,
    resolve_direction_file,
)
from hist_analysis.python.projections import project_semantic_th2


ROOT_FILES_DIR = Path(__file__).resolve().parents[1] / "root_files"
PRODUCTION_DIR = ROOT_FILES_DIR / "produced_histograms"
ETA_KEY_PATTERN = re.compile(
    r"^(?!.*(?:XPb|Xp|Run_))(h(?:Gen|Ref|RefSel|Reco)Dijet[A-Za-z0-9]*PtEta(?:CM|LabUnflipped|Lab|Forward|Backward)[A-Za-z0-9]*)_(\d+)$"
)
UNFOLDING_MARKERS = (
    "Unfold", "VsReco", "Miss", "Fake", "PairClassification",
    "ControlledSample", "StudiedSample", "MixedPrior",
)


def _interval_tag(interval: tuple[float, float]) -> str:
    def number(value: float) -> str:
        return f"{value:g}".replace("-", "m").replace(".", "p")
    return f"{number(interval[0])}to{number(interval[1])}"


def projection_key(source_key: str, binning: str,
                   interval: tuple[float, float], *,
                   eta_cut_values=DIJET_ETA_CUTS) -> str:
    """Return the stable flat ROOT key for a persisted dijet eta projection."""

    match = ETA_KEY_PATTERN.match(source_key)
    if not match:
        raise ValueError(f"Unsupported dijet eta source key {source_key!r}")
    base, eta_index = match.groups()
    index = int(eta_index)
    if not 0 <= index < len(eta_cut_values):
        raise ValueError(f"Unknown eta-cut index {index} in {source_key!r}")
    cut = f"{eta_cut_values[index]:.1f}".replace(".", "p")
    return f"{base.replace('PtEta', 'Eta')}_etaCut{cut}_{binning}_ptAve{_interval_tag(interval)}"


def _source_cut_values(source, keys):
    """Cross-check titles before assigning physical labels to indexed inputs."""
    observed = {}
    for key in keys:
        index = int(ETA_KEY_PATTERN.match(key).group(2))
        title = source.Get(key).GetTitle()
        match = re.search(r"\|#eta[^|]*\|\s*<\s*(\d+(?:\.\d+)?)", title)
        if match:
            value = float(match.group(1))
            if index in observed and observed[index] != value:
                raise ValueError(f"Conflicting eta-cut titles for index {index}")
            observed[index] = value
    values = DIJET_ETA_CUTS
    for index, value in observed.items():
        if index >= len(values) or values[index] != value:
            raise ValueError(
                f"Source eta-cut index {index} is {value} according to its title; "
                "it differs from the current production table. Regenerate the source "
                "with the current cut table before producing reusable histograms."
            )
    return values


def resolve_eta_source_key(filename: str | Path, source_key: str) -> str:
    """Resolve a current cut selection against a produced file's cut mapping."""
    match = ETA_KEY_PATTERN.match(source_key)
    if not match:
        return source_key
    requested = DIJET_ETA_CUTS[int(match.group(2))]
    file = open_root_file(str(filename))
    try:
        metadata = file.Get("histogramProductionMetadata")
        if metadata:
            table = json.loads(metadata.GetTitle()).get("eta_cut_index_to_value")
            if table:
                indices = [int(i) for i, value in table.items() if float(value) == requested]
                if len(indices) != 1:
                    raise ValueError(f"{filename} has no unique source cut {requested}; regenerate the source")
                return f"{match.group(1)}_{indices[0]}"
        histogram = file.Get(source_key)
        if histogram:
            title = re.search(r"\|#eta[^|]*\|\s*<\s*(\d+(?:\.\d+)?)", histogram.GetTitle())
            if title and float(title.group(1)) != requested:
                raise ValueError(f"{source_key} source title disagrees with requested eta cut {requested}")
        return source_key
    finally:
        file.Close()


def load_eta_projection(filename: str | Path, source_key: str, binning: str,
                        interval: tuple[float, float]):
    """Read a physical-cut projection, accepting the old index-named key."""
    resolved = resolve_eta_source_key(filename, source_key)
    key = projection_key(source_key, binning, interval)
    try:
        return load_produced_histogram(filename, key)
    except KeyError:
        base, index = ETA_KEY_PATTERN.match(resolved).groups()
        legacy = f"{base.replace('PtEta', 'Eta')}_etaCut{index}_{binning}_ptAve{_interval_tag(interval)}"
        return load_produced_histogram(filename, legacy)


def produced_file(*, category: str, dataset: str, direction: str,
                  trigger: str | None = None) -> Path:
    """Resolve one dedicated production output beneath ``root_files``."""

    parts = [category, dataset, direction]
    if trigger is not None:
        parts.append(trigger.lower())
    return produced_histogram_path(*parts)


def source_catalog(base_dir: str | Path = BASE_DIR) -> dict[str, dict]:
    """Return every explicitly supported production source and output path."""

    base = Path(base_dir)
    catalog: dict[str, dict] = {}
    for direction in ("Pbgoing", "pgoing", "combined"):
        for trigger in ("MinimumBias", "Jet60", "Jet80", "Jet100"):
            identifier = f"data_{direction}_{trigger.lower()}"
            catalog[identifier] = {
                "category": "data", "dataset": "data", "direction": direction,
                "trigger": trigger,
                "source": resolve_data_file(
                    base / "exp", trigger, direction, prefer_produced=False,
                ),
                "output": produced_file(
                    category="data", dataset="data", direction=direction,
                    trigger=trigger,
                ),
                "copy_unfolding": False,
            }
    for dataset in ("embedding", "pythia"):
        for direction in ("Pbgoing", "pgoing", "combined"):
            source = (
                resolve_combined_file(
                    base, dataset, "jetId", prefer_produced=False,
                )
                if direction == "combined"
                else resolve_direction_file(
                    base, dataset, direction, "jetId", prefer_produced=False,
                )
            )
            identifier = f"official_mc_{dataset}_{direction}"
            catalog[identifier] = {
                "category": "official_mc", "dataset": dataset,
                "direction": direction, "trigger": None, "source": source,
                "output": produced_file(
                    category="official_mc", dataset=dataset,
                    direction=direction,
                ),
                "copy_unfolding": True,
            }
    catalog["official_mc_pythia_cp5_pp"] = {
        "category": "official_mc", "dataset": "pythia_cp5",
        "direction": "pp", "trigger": None,
        "source": base / "pythia" / "cp5tune" / "HiForestMiniAOD_hist.root",
        "output": produced_file(
            category="official_mc", dataset="pythia_cp5", direction="pp",
        ),
        "copy_unfolding": True,
    }
    standalone = {
        "cuep8m1": base / "pythia" / "cuep8m1tune" / "cuep8m1_gen_dijet_histograms.root",
        "cp5": base / "pythia" / "cp5tune" / "cp5_gen_dijet_histograms.root",
    }
    for dataset, source in standalone.items():
        identifier = f"standalone_{dataset}_pp"
        catalog[identifier] = {
            "category": "standalone", "dataset": dataset,
            "direction": "pp", "trigger": None, "source": source,
            "output": produced_file(
                category="standalone", dataset=dataset, direction="pp",
            ),
            "copy_unfolding": False,
        }
    return catalog


def resolve_catalog_file(identifier: str, *, base_dir: str | Path = BASE_DIR) -> Path:
    """Prefer the dedicated catalog output, retaining its original-source fallback."""

    catalog = source_catalog(base_dir)
    identifier = identifier.replace("__", "_")
    if identifier not in catalog:
        raise ValueError(f"Unsupported production source {identifier!r}")
    entry = catalog[identifier]
    if entry["output"].is_file():
        return entry["output"]
    parts = [entry["category"], entry["dataset"], entry["direction"]]
    if entry["trigger"] is not None:
        parts.append(entry["trigger"].lower())
    legacy = find_produced_histogram_file(*parts)
    return legacy if legacy is not None else entry["source"]


def _clone_detached(obj, name: str | None = None):
    clone = obj.Clone(name or obj.GetName())
    if hasattr(clone, "SetDirectory"):
        clone.SetDirectory(0)
    return clone


def _write_object(output, obj, name: str) -> None:
    output.cd()
    obj.Write(name, ROOT.TObject.kOverwrite)


def produce_source(entry: dict, *, overwrite: bool = True) -> dict:
    """Produce one dedicated ROOT file from a catalog entry."""

    source_path = Path(entry["source"])
    output_path = Path(entry["output"])
    if not source_path.is_file():
        raise FileNotFoundError(source_path)
    if source_path.resolve() == output_path.resolve():
        raise ValueError("Production output must not overwrite its source")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    source = open_root_file(str(source_path))
    output = None
    projected_keys: list[str] = []
    copied_keys: list[str] = []
    try:
        source_keys = [key.GetName() for key in source.GetListOfKeys()]
        eta_keys = sorted(
            name for name in source_keys
            if ETA_KEY_PATTERN.match(name)
            and source.Get(name).InheritsFrom("TH2")
        )
        cut_values = _source_cut_values(
            source, eta_keys)
        # Validate the complete source-cut mapping before opening the output.
        for source_key in eta_keys:
            projection_key(source_key, "standard", DIJET_PTAVE_BINS[0], eta_cut_values=cut_values)
        output = open_root_file(str(output_path), "RECREATE" if overwrite else "UPDATE")
        for source_key in eta_keys:
            source_histogram = source.Get(source_key)
            if not source_histogram or not source_histogram.InheritsFrom("TH2"):
                raise TypeError(f"{source_path}:{source_key} is not a TH2")
            for binning, intervals in (
                ("standard", DIJET_PTAVE_BINS),
                ("test", TEST_DIJET_PTAVE_BINS),
            ):
                for interval in intervals:
                    key = projection_key(source_key, binning, interval, eta_cut_values=cut_values)
                    projection = project_semantic_th2(
                        source_histogram, "eta", interval, name=key,
                    )
                    _write_object(output, projection, key)
                    projected_keys.append(key)

        for key in sorted(source_keys):
            obj = source.Get(key)
            if not obj or (
                not obj.InheritsFrom("TH1")
                and not obj.InheritsFrom("THnSparse")
                and not obj.InheritsFrom("TNamed")
            ):
                continue
            _write_object(output, _clone_detached(obj, key), key)
            copied_keys.append(key)

        source_stat = source_path.stat()
        metadata = {
            "schema_version": 2,
            "eta_cut_index_to_value": dict(enumerate(cut_values)),
            "eta_cut_mapping_source": "current_production_table",
            "eta_cut_convention": "upstream cuts on both jets; no additional dijet-eta cut",
            "projection_key_convention": "etaCut<physical value with decimal p>",
            "category": entry["category"], "dataset": entry["dataset"],
            "direction": entry["direction"], "trigger": entry["trigger"],
            "source": {
                "path": str(source_path.resolve()), "size": source_stat.st_size,
                "mtime_ns": source_stat.st_mtime_ns,
            },
            "binning": {
                "standard": DIJET_PTAVE_BINS, "test": TEST_DIJET_PTAVE_BINS,
                "interval_convention": "half-open [low, high)",
            },
            "projected_keys": projected_keys,
            "copied_source_keys": copied_keys,
        }
        ROOT.TNamed(
            "histogramProductionMetadata",
            json.dumps(metadata, sort_keys=True, separators=(",", ":")),
        ).Write("histogramProductionMetadata", ROOT.TObject.kOverwrite)
        output.Write("", ROOT.TObject.kOverwrite)
    finally:
        if output is not None:
            output.Close()
        source.Close()
    return {
        "source": source_path, "output": output_path,
        "projected": len(projected_keys), "copied": len(copied_keys),
    }


def audit_source_cuts(entry: dict) -> dict:
    """Preflight an input cut catalog without creating or modifying output files."""
    source = open_root_file(str(entry["source"]))
    try:
        keys = [key.GetName() for key in source.GetListOfKeys()
                if ETA_KEY_PATTERN.match(key.GetName())
                and source.Get(key.GetName()).InheritsFrom("TH2")]
        values = _source_cut_values(source, keys)
        for key in keys:
            projection_key(key, "standard", DIJET_PTAVE_BINS[0], eta_cut_values=values)
        return {"source": Path(entry["source"]), "eta_histograms": len(keys),
                "eta_cut_index_to_value": dict(enumerate(values))}
    finally:
        source.Close()


def produce_all(entries: Iterable[dict] | None = None) -> list[dict]:
    """Produce the complete configured catalog in deterministic order."""

    selected = list(source_catalog().values()) if entries is None else list(entries)
    for entry in selected:
        audit_source_cuts(entry)
    return [produce_source(entry) for entry in selected]


def validate_epps21(path: str | Path = ROOT_FILES_DIR / "epps21_pPb_pp_pdf_scale_uncertainties_symmetric_fitresidual.root") -> dict:
    """Validate the existing theory artifact without changing its bin catalog."""

    path = Path(path)
    root_file = open_root_file(str(path))
    try:
        keys = [key.GetName() for key in root_file.GetListOfKeys()]
    finally:
        root_file.Close()
    if not keys:
        raise ValueError(f"EPPS21 ROOT file is empty: {path}")
    return {"path": path, "top_level_keys": keys}


def load_produced_histogram(filename: str | Path, key: str):
    """Read one already-produced histogram; never access or project a source."""

    root_file = open_root_file(str(filename))
    try:
        obj = root_file.Get(key)
        if not obj:
            raise KeyError(f"{key!r} not found in produced file {filename}")
        if not obj.InheritsFrom("TH1"):
            raise TypeError(f"{filename}:{key} is not a histogram")
        return _clone_detached(obj)
    finally:
        root_file.Close()
