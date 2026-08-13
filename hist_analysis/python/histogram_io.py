"""
histogram_io.py

Utilities for reading and writing ROOT histograms.

Author:
    Grigory Nigmatkulov

Requirements:
    ROOT >= 6.28
"""

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
    clone = obj.Clone()
    if hasattr(clone, "SetDirectory"):
        clone.SetDirectory(0)
    return clone

#------------------------------------------------------------
# Path resolution
#------------------------------------------------------------
def resolve_direction_file(base_dir: Path,
                           generator: str,
                           direction: str,
                           stem: str) -> Path:
    """
    Build a direction-specific output path.
    """
    return base_dir / generator / direction / f"{generator}_{direction}_{stem}.root"


#------------------------------------------------------------
# Path resolution
#------------------------------------------------------------
def resolve_combined_file(base_dir: Path,
                          generator: str,
                          stem: str) -> Path:
    """
    Build a combined-orientation output path.
    """

    return base_dir / generator / f"{generator}_{stem}.root"


def resolve_data_file(data_dir: Path,
                      trigger: str,
                      direction: str = "combined",
                      selection: str = "jetId",
                      vertex_filter: str = "dz1p0") -> Path:
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
