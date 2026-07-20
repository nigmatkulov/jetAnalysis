"""Existence checks for the configured histogram-analysis directory tree."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path


def _check_path_tree(node, prefix: str = "") -> dict[str, bool]:
    results: dict[str, bool] = {}

    if isinstance(node, Path):
        results[prefix] = node.exists()
        return results

    if isinstance(node, Mapping):
        for key, value in node.items():
            child_prefix = f"{prefix}.{key}" if prefix else str(key)
            results.update(_check_path_tree(value, child_prefix))
        return results

    raise TypeError(f"Unsupported node type for inventory check: {type(node)!r}")


def config_path_status(base_dir: Path,
                       data_dirs: Mapping[str, Path],
                       mc_dirs: Mapping[str, Mapping[str, Path]]) -> dict[str, bool]:
    """
    Check that the configured inventory paths exist.
    """

    results = {"BASE_DIR": base_dir.exists()}
    results.update(_check_path_tree(data_dirs, "DATA_DIRS"))
    results.update(_check_path_tree(mc_dirs, "MC_DIRS"))
    return results
