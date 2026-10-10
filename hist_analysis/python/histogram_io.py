"""
histogram_io.py

Utilities for reading and writing ROOT histograms.

Author:
    Grigory Nigmatkulov

Requirements:
    ROOT >= 6.28
"""

import os
from pathlib import Path
from typing import TYPE_CHECKING, Any

import ROOT

# Prevent histograms from being automatically attached
# to the currently opened ROOT file.
ROOT.TH1.AddDirectory(False)

if TYPE_CHECKING:
    ROOTObject = Any
else:
    ROOTObject = Any

_SUPPORTED_HISTOGRAM_CLASSES = ("TH1", "TH2", "TH3", "THnSparse")
_SUPPORTED_GENERIC_CLASSES = _SUPPORTED_HISTOGRAM_CLASSES + ("TGraph",)
_DATA_TRIGGER_STEMS = {
    "MinimumBias": "mb",
    "MB": "mb",
    "mb": "mb",
    "Jet60": "jet60",
    "jet60": "jet60",
    "Jet80": "jet80",
    "jet80": "jet80",
    "Jet100": "jet100",
    "jet100": "jet100",
}
_DATA_DIRECTIONS = ("combined", "Pbgoing", "pgoing")
_DATA_RECO_SELECTIONS = ("jetId", "trkMax", "noSel")
_DATA_VERTEX_FILTERS = ("dz1p0", "Gplus", "Vtx1")


# ------------------------------------------------------------
# Path handling
# ------------------------------------------------------------
def list_root_files(base_dir: Path) -> list[Path]:
    """
    Return all ROOT files under a directory.
    """

    if not base_dir.exists():
        return []
    return sorted(p for p in base_dir.rglob("*.root") if p.is_file())

#------------------------------------------------------------
# Object handling
#------------------------------------------------------------
def _is_supported_object(obj, class_names: tuple[str, ...]) -> bool:
    return any(obj.InheritsFrom(class_name) for class_name in class_names)

#------------------------------------------------------------
# Object cloning
#------------------------------------------------------------
def _clone_object(obj):
    # A borrowed TFile::Get object dies with its file. Clone while the file is
    # open, then detach TH1-derived objects before the caller closes the file.
    clone = obj.Clone()
    if hasattr(clone, "SetDirectory"):
        clone.SetDirectory(0)
    return clone

#------------------------------------------------------------
# Path resolution
#------------------------------------------------------------
def produced_histogram_path(*parts: str) -> Path:
    """Canonical single-underscore filename for a reusable histogram bundle."""
    return (Path(__file__).resolve().parents[1] / "root_files"
            / "produced_histograms" / ("_".join(parts) + ".root"))


def find_produced_histogram_file(*parts: str) -> Path | None:
    """Prefer canonical bundles, accepting legacy double-underscore filenames."""
    canonical = produced_histogram_path(*parts)
    legacy = canonical.parent / ("__".join(parts) + ".root")
    return next((path for path in (canonical, legacy) if path.is_file()), None)


def resolve_direction_file(base_dir: Path,
                           generator: str,
                           direction: str,
                           stem: str,
                           *,
                           prefer_produced: bool = True) -> Path:
    """
    Build a direction-specific output path.
    """
    if prefer_produced and stem == "jetId":
        produced = find_produced_histogram_file("official_mc", generator, direction)
        if produced is not None:
            return produced
    return base_dir / generator / direction / f"{generator}_{direction}_{stem}.root"


#------------------------------------------------------------
# Path resolution
#------------------------------------------------------------
def resolve_combined_file(base_dir: Path,
                          generator: str,
                          stem: str,
                          *,
                          prefer_produced: bool = True) -> Path:
    """
    Build a combined-orientation output path.
    """

    if prefer_produced and stem == "jetId":
        produced = find_produced_histogram_file("official_mc", generator, "combined")
        if produced is not None:
            return produced
    return base_dir / generator / f"{generator}_{stem}.root"


def normalize_unfolding_direction(dataset: str, direction: str) -> str:
    """Return the physical direction label for an unfolding dataset.

    The pp CP5 production has no beam orientation.  Treat any configured
    direction as ``pp`` so switching the dataset does not also require editing
    a stale pPb direction setting.
    """

    if dataset in {"pythia_cp5", "pythia_cuep8m1", "pythia_cp5_standalone"}:
        return "pp"
    if dataset in {"embedding", "pythia"}:
        if direction not in {"combined", "Pbgoing", "pgoing"}:
            raise ValueError(
                f"{dataset} direction must be 'combined', 'Pbgoing', or "
                f"'pgoing', got {direction!r}"
            )
        return direction
    raise ValueError(f"Unsupported unfolding dataset {dataset!r}")


def resolve_unfolding_mc_file(base_dir: Path,
                              dataset: str,
                              direction: str,
                              stem: str = "jetId") -> Path:
    """Resolve an MC file that contains a complete unfolding response.

    ``embedding`` and ``pythia`` use the standard direction-aware production
    layout. ``pythia_cp5`` is the pp MiniAOD analysis output and therefore has
    no beam direction; callers must select ``pp`` for that dataset.
    """

    direction = normalize_unfolding_direction(dataset, direction)
    produced = find_produced_histogram_file("official_mc", dataset, direction)
    if dataset == "pythia_cp5" and "CP5_UNFOLDING_FILE" in os.environ:
        return Path(os.environ["CP5_UNFOLDING_FILE"]).expanduser()
    if stem == "jetId" and produced is not None:
        return produced
    if dataset in {"embedding", "pythia"}:
        if direction == "combined":
            return resolve_combined_file(base_dir, dataset, stem)
        return resolve_direction_file(base_dir, dataset, direction, stem)
    if dataset == "pythia_cp5":
        return Path(os.environ.get(
            "CP5_UNFOLDING_FILE",
            base_dir / "pythia" / "cp5tune" / "HiForestMiniAOD_hist.root",
        )).expanduser()
    raise ValueError(
        f"Unsupported unfolding MC dataset {dataset!r}; expected "
        "'embedding', 'pythia', or 'pythia_cp5'"
    )


def resolve_unfolding_truth_file(base_dir: Path,
                                 dataset: str,
                                 direction: str,
                                 stem: str = "jetId") -> Path:
    """Resolve a generator-truth source used as an unfolding/closure target."""

    direction = normalize_unfolding_direction(dataset, direction)
    standalone_names = {
        "pythia_cuep8m1": "cuep8m1",
        "pythia_cp5_standalone": "cp5",
    }
    overrides = {
        "pythia_cuep8m1": "CUEP8M1_GEN_FILE",
        "pythia_cp5_standalone": "CP5_GEN_FILE",
    }
    if dataset in overrides and overrides[dataset] in os.environ:
        return Path(os.environ[overrides[dataset]]).expanduser()
    if dataset in standalone_names:
        produced = find_produced_histogram_file("standalone", standalone_names[dataset], "pp")
        if produced is not None:
            return produced
    if dataset in {"embedding", "pythia", "pythia_cp5"}:
        return resolve_unfolding_mc_file(base_dir, dataset, direction, stem)
    standalone_defaults = {
        "pythia_cuep8m1": (
            "CUEP8M1_GEN_FILE",
            base_dir / "pythia" / "cuep8m1tune"
            / "cuep8m1_gen_dijet_histograms.root",
        ),
        "pythia_cp5_standalone": (
            "CP5_GEN_FILE",
            base_dir / "pythia" / "cp5tune" / "cp5_gen_dijet_histograms.root",
        ),
    }
    if dataset in standalone_defaults:
        environment_name, default = standalone_defaults[dataset]
        return Path(os.environ.get(environment_name, default)).expanduser()
    raise ValueError(
        f"Unsupported unfolding truth dataset {dataset!r}; expected embedding, "
        "pythia, pythia_cp5, pythia_cuep8m1, or pythia_cp5_standalone"
    )


def resolve_data_file(data_dir: Path,
                      trigger: str,
                      direction: str = "combined",
                      selection: str = "jetId",
                      vertex_filter: str = "dz1p0",
                      *,
                      prefer_produced: bool = True) -> Path:
    """Build a merged experimental-data output path.

    Trigger aliases used by the data notebooks map to the lowercase stems
    written by ``processing/merge_data_outputs.py``. Combined files live in
    ``data_dir``; direction-specific files live in their matching subdirectory
    and include the direction in the filename. Non-nominal pileup-filter
    productions append ``_Gplus`` or ``_Vtx1``; nominal ``dz1p0`` remains
    unsuffixed.
    """

    try:
        trigger_stem = _DATA_TRIGGER_STEMS[trigger]
    except KeyError as error:
        supported = ", ".join(sorted(_DATA_TRIGGER_STEMS))
        raise ValueError(
            f"Unsupported data trigger {trigger!r}; expected one of: {supported}"
        ) from error
    if direction not in _DATA_DIRECTIONS:
        supported = ", ".join(_DATA_DIRECTIONS)
        raise ValueError(
            f"Unsupported data direction {direction!r}; expected one of: {supported}"
        )
    if selection not in _DATA_RECO_SELECTIONS:
        supported = ", ".join(_DATA_RECO_SELECTIONS)
        raise ValueError(
            f"Unsupported reco-jet selection {selection!r}; expected one of: {supported}"
        )
    if vertex_filter not in _DATA_VERTEX_FILTERS:
        supported = ", ".join(_DATA_VERTEX_FILTERS)
        raise ValueError(
            f"Unsupported vertex filter {vertex_filter!r}; expected one of: "
            f"{supported}"
        )

    data_dir = Path(data_dir)
    if prefer_produced and selection == "jetId" and vertex_filter == "dz1p0":
        produced_trigger = "minimumbias" if trigger_stem == "mb" else trigger_stem
        produced = find_produced_histogram_file("data", "data", direction, produced_trigger)
        if produced is not None:
            return produced
    vertex_suffix = "" if vertex_filter == "dz1p0" else f"_{vertex_filter}"
    if direction == "combined":
        return data_dir / f"{trigger_stem}_ak4_{selection}{vertex_suffix}.root"
    return (
        data_dir / direction
        / f"{trigger_stem}_{direction}_ak4_{selection}{vertex_suffix}.root"
    )

# ------------------------------------------------------------
# File handling
# ------------------------------------------------------------
def open_root_file(filename: str,
                   mode: str = "READ") -> ROOT.TFile:
    """
    Open a ROOT file safely.

    Parameters
    ----------
    filename
        Path to ROOT file.

    mode
        READ, RECREATE, UPDATE, ...

    Returns
    -------
    ROOT.TFile
    """

    path = Path(filename)
    if mode == "READ" and not path.exists():
        raise FileNotFoundError(path)

    root_file = ROOT.TFile.Open(str(path), mode)
    if root_file is None or root_file.IsZombie():
        raise OSError(f"Cannot open ROOT file: {filename}")

    return root_file


# ------------------------------------------------------------
# Object listing
# ------------------------------------------------------------
def list_keys(filename: str) -> List[str]:
    """
    Return names of all top-level objects.
    """

    root_file = open_root_file(filename)
    try:
        names = [key.GetName() for key in root_file.GetListOfKeys()]
        return sorted(names)
    finally:
        root_file.Close()


#------------------------------------------------------------
# Object printing
#------------------------------------------------------------
def print_keys(filename: str) -> None:
    """
    Print ROOT file content.
    """

    keys = list_keys(filename)
    if len(keys) == 0:
        print(f"{filename} is empty.")
    else:
        print(f"Keys in {filename}:")
        for name in list_keys(filename):
            print(name)


# ------------------------------------------------------------
# Histogram loading
# ------------------------------------------------------------
def load_histogram(filename: str,
                   histogram_name: str):
    """
    Load one histogram-like ROOT object.

    Supported classes:

        TH1
        TH2
        TH3
        THnSparse

    Returns a detached clone when the object supports directory ownership.
    """

    root_file = open_root_file(filename)
    try:
        obj = root_file.Get(histogram_name)
        if not obj:
            raise KeyError(
                f"Histogram '{histogram_name}' not found."
            )

        if not _is_supported_object(obj, _SUPPORTED_HISTOGRAM_CLASSES):
            raise TypeError(
                f"{histogram_name} is not a supported histogram type."
            )

        return _clone_object(obj)
    finally:
        root_file.Close()


#------------------------------------------------------------
# Histogram loading
#------------------------------------------------------------
def load_histograms(filename: str,
                    histogram_names: List[str]
                    ) -> dict[str, ROOTObject]:
    """
    Load multiple histogram-like ROOT objects.

    Missing or unsupported keys are reported and omitted. Supported classes:

        TH1
        TH2
        TH3
        THnSparse
    """

    histograms: dict[str, ROOTObject] = {}
    root_file = open_root_file(filename)
    try:
        for name in histogram_names:
            obj = root_file.Get(name)
            if not obj:
                print(f"WARNING: {name} not found")
                continue

            if not _is_supported_object(obj, _SUPPORTED_HISTOGRAM_CLASSES):
                print(f"WARNING: {name} is not a supported histogram type")
                continue

            histograms[name] = _clone_object(obj)
        return histograms
    finally:
        root_file.Close()


# ------------------------------------------------------------
# Generic object loading
# ------------------------------------------------------------
def load_object(filename: str,
                object_name: str):
    """
    Load a supported histogram or graph as a detached clone.

    Supported classes:

        TH1
        TH2
        TH3
        THnSparse
        TGraph
    """

    root_file = open_root_file(filename)
    try:
        obj = root_file.Get(object_name)
        if not obj:
            raise KeyError(object_name)

        if not _is_supported_object(obj, _SUPPORTED_GENERIC_CLASSES):
            raise TypeError(
                f"{object_name} is not a supported ROOT object type."
            )

        return _clone_object(obj)
    finally:
        root_file.Close()


# ------------------------------------------------------------
# Histogram saving
# ------------------------------------------------------------
def save_histograms(filename: str,
                    histograms,
                    mode: str = "RECREATE") -> None:
    """
    Write non-null histogram-like ROOT objects using the requested file mode.
    """

    root_file = open_root_file(filename, mode)
    try:
        for hist in histograms:
            if hist is None:
                continue
            hist.Write()
    finally:
        root_file.Close()


# ------------------------------------------------------------
# Histogram information
# ------------------------------------------------------------
def histogram_info(hist) -> None:
    """
    Print a ROOT object summary.
    """

    print("--------------------------------")
    print("Name      :", hist.GetName())
    print("Title     :", hist.GetTitle())
    print("Class     :", hist.ClassName())
    print("Entries   :", hist.GetEntries())
    print("Integral  :", hist.Integral())
    print("--------------------------------")


# ------------------------------------------------------------
# Histogram existence
# ------------------------------------------------------------
def histogram_exists(filename: str,
                     histogram_name: str) -> bool:
    """
    Check whether a named object exists in the ROOT file.
    """

    root_file = open_root_file(filename)
    try:
        return root_file.Get(histogram_name) is not None
    finally:
        root_file.Close()
