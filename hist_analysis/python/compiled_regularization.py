"""ctypes bridge to the standalone compiled Bayesian regularization scan."""

from __future__ import annotations

from array import array
import ctypes
from dataclasses import dataclass
import os
from pathlib import Path
from typing import Sequence


@dataclass(frozen=True)
class CompiledScanMoments:
    n_truth_bins: int
    n_toys: int
    max_iterations: int
    sums: tuple[tuple[float, ...], ...]
    sums_squared: tuple[tuple[float, ...], ...]


def resolve_compiled_scan_library(project_root: str | Path) -> Path:
    """Resolve an explicitly configured or standard CMake-built scan library."""

    configured = os.environ.get("JETANALYSIS_REGULARIZATION_SCAN_LIBRARY")
    if configured:
        candidates = (Path(configured).expanduser(),)
    else:
        suffixes = (".dylib", ".so", ".dll")
        build = Path(project_root) / "build" / "regularization-scan"
        candidates = tuple(
            build / f"libjetanalysis_regularization_scan{suffix}"
            for suffix in suffixes
        ) + tuple(
            build / f"jetanalysis_regularization_scan{suffix}"
            for suffix in suffixes
        )
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    searched = ", ".join(str(candidate) for candidate in candidates)
    raise FileNotFoundError(
        "Compiled regularization scan library not found. Build it with "
        "`cmake -S hist_analysis/cpp -B build/regularization-scan "
        "-DCMAKE_BUILD_TYPE=Release` followed by "
        "`cmake --build build/regularization-scan --parallel`, or set "
        f"JETANALYSIS_REGULARIZATION_SCAN_LIBRARY. Searched: {searched}"
    )


def _double_pointer(values: array):
    return (ctypes.c_double * len(values)).from_buffer(values)


def scan_compiled_bayes_moments(
    response_matrix,
    training_prior,
    truth_efficiency,
    measured_toys: Sequence,
    *,
    max_iterations: int,
    project_root: str | Path,
    worker_count: int = 0,
) -> CompiledScanMoments:
    """Run the numeric toy scan and return iteration-wise ensemble moments."""

    toys = tuple(measured_toys)
    if len(toys) < 2 or max_iterations < 1 or worker_count < 0:
        raise ValueError("At least two toys, one iteration, and workers >= 0 are required")
    n_measured = response_matrix.GetNbinsX()
    n_truth = response_matrix.GetNbinsY()
    if training_prior.GetNbinsX() != n_truth:
        raise ValueError("Compiled scan prior does not match response truth bins")
    if truth_efficiency.GetNbinsX() != n_truth:
        raise ValueError("Compiled scan efficiency does not match response truth bins")
    if any(toy.GetNbinsX() != n_measured for toy in toys):
        raise ValueError("Compiled scan toy does not match response measured bins")

    # ABI is row-major [measured, truth], matching ROOT x/y semantics above.
    # Keep all array buffers alive through the ctypes call; the C++ code reads
    # their addresses synchronously and returns moments, not a full covariance.
    response = array("d", (
        response_matrix.GetBinContent(measured, truth)
        for measured in range(1, n_measured + 1)
        for truth in range(1, n_truth + 1)
    ))
    prior = array("d", (
        training_prior.GetBinContent(truth) for truth in range(1, n_truth + 1)
    ))
    efficiency = array("d", (
        truth_efficiency.GetBinContent(truth) for truth in range(1, n_truth + 1)
    ))
    toy_values = array("d", (
        toy.GetBinContent(measured)
        for toy in toys
        for measured in range(1, n_measured + 1)
    ))
    output_size = max_iterations * n_truth
    sums = array("d", [0.0]) * output_size
    sums_squared = array("d", [0.0]) * output_size

    library_path = resolve_compiled_scan_library(project_root)
    library = ctypes.CDLL(str(library_path))
    function = library.jetanalysis_regularization_scan_moments
    double_pointer = ctypes.POINTER(ctypes.c_double)
    function.argtypes = [
        ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t,
        double_pointer, double_pointer, double_pointer, double_pointer,
        ctypes.c_size_t, double_pointer, double_pointer,
        ctypes.POINTER(ctypes.c_char), ctypes.c_size_t,
    ]
    function.restype = ctypes.c_int
    error = ctypes.create_string_buffer(2048)
    status = function(
        n_measured, n_truth, len(toys), max_iterations,
        _double_pointer(response), _double_pointer(prior),
        _double_pointer(efficiency), _double_pointer(toy_values), worker_count,
        _double_pointer(sums), _double_pointer(sums_squared), error, len(error),
    )
    if status != 0:
        message = error.value.decode("utf-8", errors="replace")
        raise RuntimeError(f"Compiled regularization scan failed: {message}")

    return CompiledScanMoments(
        n_truth, len(toys), max_iterations,
        tuple(tuple(sums[iteration * n_truth:(iteration + 1) * n_truth])
              for iteration in range(max_iterations)),
        tuple(tuple(sums_squared[iteration * n_truth:(iteration + 1) * n_truth])
              for iteration in range(max_iterations)),
    )
