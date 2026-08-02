"""Shared environment setup for the histogram-analysis notebooks."""

from __future__ import annotations

import importlib
import os
from pathlib import Path
import shutil
import subprocess
import sys


def find_project_root(start: str | Path | None = None) -> Path:
    """Return the nearest parent containing this project's stable sentinels."""

    current = Path.cwd() if start is None else Path(start)
    current = current.expanduser().resolve()
    candidates = (current, *current.parents)
    for candidate in candidates:
        if ((candidate / "CMakeLists.txt").is_file()
                and (candidate / "hist_analysis").is_dir()):
            return candidate
    raise RuntimeError(
        "Cannot locate the jetAnalysis repository. Start Jupyter from the "
        "repository root with py-env/bin/python -m jupyter notebook."
    )


def load_root(*, batch: bool = True):
    """Import PyROOT from the active environment and configure notebook mode."""

    try:
        root = importlib.import_module("ROOT")
    except ModuleNotFoundError as error:
        if error.name != "ROOT":
            raise

        root_config = shutil.which("root-config")
        if root_config is not None:
            root_libdir = Path(
                subprocess.check_output(
                    [root_config, "--libdir"], text=True,
                ).strip()
            )
            if root_libdir.is_dir() and str(root_libdir) not in sys.path:
                sys.path.insert(0, str(root_libdir))
            try:
                root = importlib.import_module("ROOT")
            except ModuleNotFoundError as retry_error:
                if retry_error.name != "ROOT":
                    raise
                error = retry_error
            else:
                error = None

        if error is not None:
            raise ModuleNotFoundError(
                "PyROOT is unavailable in the active interpreter "
                f"({sys.executable}) and could not be resolved with root-config. "
                "Start Jupyter from the repository root with "
                "py-env/bin/python -m jupyter notebook."
            ) from error

    root.gROOT.SetBatch(batch)
    return root


def load_roounfold(root, *, project_root: str | Path,
                   checkout: str | Path | None = None) -> tuple[Path, Path]:
    """Load RooUnfold and return its checkout and shared-library paths."""

    project_root = Path(project_root).expanduser().resolve()
    configured_root = checkout or os.environ.get("ROOUNFOLD_ROOT")
    roounfold_root = (
        Path(configured_root).expanduser().resolve()
        if configured_root
        else project_root.parent / "RooUnfold"
    )

    extension = str(root.gSystem.GetSoExt()).lstrip(".")
    library_name = f"libRooUnfold.{extension}"
    candidates = (
        roounfold_root / "build" / library_name,
        roounfold_root / "build" / "lib" / library_name,
    )
    library = next((path for path in candidates if path.is_file()), None)
    if library is None:
        searched = ", ".join(str(path) for path in candidates)
        raise FileNotFoundError(
            f"RooUnfold library not found; searched: {searched}. Set "
            "ROOUNFOLD_ROOT to the RooUnfold checkout."
        )

    for path in (
        roounfold_root / "src",
        roounfold_root / "build",
        roounfold_root / "build" / "lib",
    ):
        if path.is_dir() and str(path) not in sys.path:
            sys.path.insert(0, str(path))

    load_status = root.gSystem.Load(str(library))
    if load_status < 0:
        raise RuntimeError(f"ROOT failed to load {library}")

    return roounfold_root, library
