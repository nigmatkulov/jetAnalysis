#!/usr/bin/env python3
"""Check availability of every ROOT file listed for pPb8160 data.

The expected production coverage is fixed by the analysis workflow:

* MB Pb-going: PD1 through PD20, inclusive;
* MB p-going: PD1 through PD8, inclusive;
* PAEGJet: one list for each beam direction.

Blank lines and comments are ignored.  If a list entry has trailing metadata,
only its first whitespace-delimited field is interpreted as the ROOT path,
matching ``processForestSimple`` list handling.
"""

from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


PROCESSING_DIR = Path(__file__).resolve().parent
DEFAULT_FILELIST_ROOT = PROCESSING_DIR / "filelists" / "pPb8160"


@dataclass(frozen=True)
class FileListSpec:
    """One expected production file list and its readable display label."""

    label: str
    path: Path


@dataclass(frozen=True)
class FileListResult:
    """Availability results for one file list."""

    spec: FileListSpec
    available: tuple[Path, ...]
    missing: tuple[Path, ...]


class ProgressBar:
    """Render overall file-check progress without external dependencies."""

    def __init__(self, total: int, width: int = 36) -> None:
        self.total = total
        self.width = width
        self.completed = 0
        self.interactive = sys.stderr.isatty()
        self.last_render_time = 0.0
        if self.interactive:
            self._render()
        else:
            print(f"Checking {total} ROOT files...", file=sys.stderr)

    def advance(self) -> None:
        """Advance the bar by one completed availability check."""

        self.completed += 1
        now = time.monotonic()
        if self.interactive and (
            self.completed == self.total or now - self.last_render_time >= 0.1
        ):
            self._render()

    def finish(self) -> None:
        """Finish the current progress line."""

        if self.interactive:
            if self.completed < self.total:
                self.completed = self.total
                self._render()
            print(file=sys.stderr)
        else:
            print(f"Completed {self.completed}/{self.total} checks.", file=sys.stderr)

    def _render(self) -> None:
        """Draw the current interactive progress state."""

        fraction = self.completed / self.total if self.total else 1.0
        filled = min(self.width, int(self.width * fraction))
        bar = "#" * filled + "-" * (self.width - filled)
        percent = 100.0 * fraction
        print(
            f"\rChecking ROOT files [{bar}] {percent:6.2f}% "
            f"({self.completed}/{self.total})",
            end="",
            file=sys.stderr,
            flush=True,
        )
        self.last_render_time = time.monotonic()


def positive_int(value: str) -> int:
    """Parse a strictly positive integer for ``--workers``."""

    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("value must be at least 1")
    return parsed


def parse_args() -> argparse.Namespace:
    """Parse checker configuration."""

    parser = argparse.ArgumentParser(
        description=(
            "Check every ROOT path in the pPb8160 MB and PAEGJet data file lists."
        )
    )
    parser.add_argument(
        "--filelist-root",
        type=Path,
        default=DEFAULT_FILELIST_ROOT,
        help=f"pPb8160 file-list directory (default: {DEFAULT_FILELIST_ROOT})",
    )
    parser.add_argument(
        "--workers",
        type=positive_int,
        default=16,
        help="Concurrent filesystem checks (default: 16)",
    )
    parser.add_argument(
        "--show-available",
        action="store_true",
        help="Print available paths in addition to missing paths",
    )
    return parser.parse_args()


def expected_filelists(root: Path) -> list[FileListSpec]:
    """Return all expected MB and PAEGJet lists in production order."""

    specs: list[FileListSpec] = []
    for direction, maximum_pd in (("Pbgoing", 20), ("pgoing", 8)):
        for pd_number in range(1, maximum_pd + 1):
            specs.append(
                FileListSpec(
                    label=f"MB {direction} PD{pd_number}",
                    path=(
                        root
                        / "DATA_MB"
                        / direction
                        / f"MB_PD{pd_number}_{direction}.txt"
                    ),
                )
            )

    for direction in ("Pbgoing", "pgoing"):
        specs.append(
            FileListSpec(
                label=f"PAEGJet {direction}",
                path=(
                    root
                    / "DATA_PAEGJet"
                    / direction
                    / f"PAEGJet_{direction}.txt"
                ),
            )
        )
    return specs


def read_root_paths(filelist: Path) -> list[Path]:
    """Read ROOT paths from one list, ignoring blank and comment lines."""

    paths: list[Path] = []
    for raw_line in filelist.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        paths.append(Path(line.split(maxsplit=1)[0]).expanduser())
    return paths


def path_is_available(path: Path) -> bool:
    """Return whether a path is a regular file, treating I/O errors as missing."""

    try:
        return path.is_file()
    except OSError:
        return False


def check_paths(
    paths: Iterable[Path], workers: int, progress: ProgressBar | None = None
) -> tuple[tuple[Path, ...], tuple[Path, ...]]:
    """Check paths concurrently and preserve their original list order."""

    ordered_paths = list(paths)
    availability = [False] * len(ordered_paths)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures: dict[Future[bool], int] = {
            executor.submit(path_is_available, path): index
            for index, path in enumerate(ordered_paths)
        }
        for future in as_completed(futures):
            availability[futures[future]] = future.result()
            if progress is not None:
                progress.advance()

    available = tuple(
        path for path, exists in zip(ordered_paths, availability) if exists
    )
    missing = tuple(
        path for path, exists in zip(ordered_paths, availability) if not exists
    )
    return available, missing


def main() -> int:
    """Check all expected lists and print a grouped availability report."""

    args = parse_args()
    filelist_root = args.filelist_root.expanduser().resolve()
    specs = expected_filelists(filelist_root)

    absent_lists = [spec for spec in specs if not spec.path.is_file()]
    if absent_lists:
        print("Missing required file lists:")
        for spec in absent_lists:
            print(f"  [{spec.label}] {spec.path}")
        print(f"\nCannot check data files: {len(absent_lists)} required list(s) missing.")
        return 2

    paths_by_spec = [(spec, read_root_paths(spec.path)) for spec in specs]
    progress = ProgressBar(sum(len(paths) for _, paths in paths_by_spec))

    results: list[FileListResult] = []
    for spec, paths in paths_by_spec:
        available, missing = check_paths(paths, args.workers, progress)
        results.append(FileListResult(spec, available, missing))
    progress.finish()

    total_available = 0
    total_missing = 0
    for result in results:
        total_available += len(result.available)
        total_missing += len(result.missing)
        total = len(result.available) + len(result.missing)
        print(
            f"[{result.spec.label}] total={total} "
            f"available={len(result.available)} missing={len(result.missing)}"
        )
        if args.show_available:
            for path in result.available:
                print(f"  AVAILABLE: {path}")
        for path in result.missing:
            print(f"  MISSING: {path}")

    print(
        f"\nChecked {len(results)} file lists: "
        f"available={total_available} missing={total_missing}"
    )
    if total_missing:
        print("Some listed ROOT files are unavailable.")
        return 1

    print("All listed ROOT files are available.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
