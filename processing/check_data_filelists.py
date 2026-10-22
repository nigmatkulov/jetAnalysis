#!/usr/bin/env python3
"""Check availability and optional ROOT integrity of pPb8160 data or specific files.

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
import json
import shutil
import subprocess
import tempfile
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
class FileProblem:
    """A specific input that failed availability or integrity validation."""

    path: Path | str
    category: str
    detail: str


@dataclass(frozen=True)
class FileListResult:
    """Availability/integrity results for one file list."""

    spec: FileListSpec
    available: tuple[Path | str, ...]
    missing: tuple[Path | str, ...]
    problems: tuple[FileProblem, ...] = ()


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
    parser.add_argument("root_files", nargs="*", help="Check specific files/ROOT URLs instead of the production lists; implies --check-root")
    parser.add_argument("--check-root", action="store_true", help="Open files with ROOT and check integrity, including stored objects")
    parser.add_argument("--deep", action="store_true", help="Also read every TTree entry/branch; implies --check-root (potentially expensive)")
    parser.add_argument("--timeout", type=positive_int, default=120, help="Maximum seconds per ROOT check (default: 120; increase for --deep)")
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
        help="Concurrent file checks (default: 16)",
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


def read_root_paths(filelist: Path) -> list[Path | str]:
    """Read ROOT paths from one list, ignoring blank and comment lines."""

    paths: list[Path | str] = []
    for raw_line in filelist.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        token = line.split(maxsplit=1)[0]
        paths.append(token if "://" in token else Path(token).expanduser())
    return paths


def path_is_available(path: Path) -> bool:
    """Return whether a path is a regular file, treating I/O errors as missing."""

    try:
        return Path(path).is_file()
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


def check_root_file(path: Path | str, deep: bool = False, timeout: int = 120) -> FileProblem | None:
    """Isolate ROOT in a subprocess so corrupt files cannot crash the scanner.

    Ordinary checks deserialize objects recursively; deep checks also read all
    TTree entries. Neither check establishes physics/processing completeness.
    """
    if "://" not in str(path):
        if not path_is_available(Path(path)):
            return FileProblem(path, "MISSING", "not a readable regular file")
        try:
            if Path(path).stat().st_size == 0:
                return FileProblem(path, "EMPTY", "zero-byte file")
        except OSError as error:
            return FileProblem(path, "UNREADABLE", str(error))
    with tempfile.TemporaryDirectory(prefix="root-input-check-") as directory:
        macro = Path(directory) / "check.C"
        macro.write_text(r"""#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TTree.h>
#include <TSystem.h>
#include <TError.h>
#include <memory>
#include <iostream>
#include <string>
bool readError = false;
void captureError(int level, Bool_t abort, const char *where, const char *message) {
 if (level >= kError) readError = true;
 DefaultErrorHandler(level, abort, where, message);
}
bool inspect(TDirectory *directory, bool deep) {
 TIter next(directory->GetListOfKeys());
 while (TKey *key = (TKey*)next()) {
  std::unique_ptr<TObject> object(key->ReadObj());
  if (!object || readError) return false;
  if (auto child = dynamic_cast<TDirectory*>(object.get())) {
   if (!inspect(child, deep)) return false;
  }
  if (deep) {
   if (auto tree = dynamic_cast<TTree*>(object.get())) {
    tree->SetBranchStatus("*", 1);
    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
     if ((tree->GetEntry(entry) <= 0 && tree->GetListOfBranches()->GetEntries() > 0) || readError) return false;
    }
   }
  }
 }
 return !readError;
}
void finish(int code) {
 std::cout << "ROOT_CHECK_RESULT=" << code << std::endl;
 gSystem->Exit(code);
}
void check() {
 SetErrorHandler(captureError);
 const std::string path = INPUT_PATH;
 std::unique_ptr<TFile> file(path.find("://") == std::string::npos
                           ? new TFile(path.c_str(), "READ") : TFile::Open(path.c_str(), "READ"));
 if (!file) finish(13);
 if (file->IsZombie()) finish(10);
 if (file->TestBit(TFile::kRecovered)) finish(11);
 if (readError || !inspect(file.get(), DEEP_CHECK)) finish(12);
 file->Close();
 finish(readError ? 12 : 0);
}
""".replace("INPUT_PATH", json.dumps(str(path))).replace("DEEP_CHECK", "true" if deep else "false"))
        try:
            result = subprocess.run(["root", "-l", "-b", "-q", str(macro)],
                                    capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return FileProblem(path, "TIMEOUT", f"ROOT check exceeded {timeout}s; integrity unknown")
        except OSError as error:
            return FileProblem(path, "CHECK_ERROR", str(error))
        codes = {10: ("ZOMBIE", "ROOT zombie file"),
                 11: ("RECOVERED", "ROOT recovered an improperly closed file"),
                 12: ("CORRUPT", "ROOT reported an error reading stored objects or tree data"),
                 13: ("UNREADABLE", "ROOT could not open the file (corruption, access or I/O failure)")}
        if f"ROOT_CHECK_RESULT={result.returncode}" not in result.stdout:
            return FileProblem(path, "CHECK_ERROR", f"ROOT crashed or failed to run (exit {result.returncode}); integrity unknown")
        if result.returncode == 0:
            return None
        category, detail = codes.get(result.returncode, ("CHECK_ERROR", "unexpected ROOT check result"))
        return FileProblem(path, category, detail)


def check_integrity_paths(paths, workers, progress, deep=False, timeout=120):
    """Collect deterministic per-file results despite concurrent completion."""
    ordered = list(paths)
    issues = [None] * len(ordered)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(check_root_file, path, deep, timeout): index
                   for index, path in enumerate(ordered)}
        for future in as_completed(futures):
            issues[futures[future]] = future.result()
            progress.advance()
    good = tuple(path for path, issue in zip(ordered, issues) if issue is None)
    problems = tuple(issue for issue in issues if issue is not None)
    missing = tuple(issue.path for issue in problems if issue.category == "MISSING")
    return good, missing, problems


def main() -> int:
    """Check production lists or explicit files and summarize every problem."""
    args = parse_args()
    integrity = args.check_root or args.deep or bool(args.root_files)
    if args.root_files:
        paths = [value if "://" in value else Path(value).expanduser() for value in args.root_files]
        paths_by_spec = [(FileListSpec("Specified files", Path(".")), paths)]
    else:
        specs = expected_filelists(args.filelist_root.expanduser().resolve())
        absent_lists = [spec for spec in specs if not spec.path.is_file()]
        if absent_lists:
            print("Missing required file lists:")
            for spec in absent_lists:
                print(f"  [{spec.label}] {spec.path}")
            print(f"\nCannot check data files: {len(absent_lists)} required list(s) missing.")
            return 2
        paths_by_spec = [(spec, read_root_paths(spec.path)) for spec in specs]
    integrity = integrity or any("://" in str(path) for _, paths in paths_by_spec for path in paths)
    if integrity and not shutil.which("root"):
        print("ROOT integrity checks require the root executable in PATH.", file=sys.stderr)
        return 2
    progress = ProgressBar(sum(len(paths) for _, paths in paths_by_spec))
    results = []
    for spec, paths in paths_by_spec:
        if integrity:
            available, missing, problems = check_integrity_paths(paths, args.workers, progress, args.deep, args.timeout)
        else:
            available, missing = check_paths(paths, args.workers, progress)
            problems = tuple(FileProblem(path, "MISSING", "not a regular file") for path in missing)
        results.append(FileListResult(spec, available, missing, problems))
    progress.finish()
    for result in results:
        total = len(result.available) + len(result.problems)
        print(f"[{result.spec.label}] total={total} available={len(result.available)} "
              f"missing={len(result.missing)} problems={len(result.problems)}")
        if args.show_available:
            for path in result.available:
                print(f"  {'VALID' if integrity else 'AVAILABLE'}: {path}")
    total_available = sum(len(result.available) for result in results)
    total_missing = sum(len(result.missing) for result in results)
    total_problems = sum(len(result.problems) for result in results)
    print(f"\nChecked {len(results)} group(s): available={total_available} missing={total_missing} problems={total_problems}")
    print("Integrity mode: " + ("full TTree entry reads" if args.deep else "ROOT open and stored-object checks" if integrity else "availability only; ROOT integrity not checked"))
    if total_problems:
        print("Problem files:")
        for result in results:
            for issue in result.problems:
                print(f"  [{result.spec.label}] {issue.category}: {issue.path} — {issue.detail}")
        return 1
    print("All listed ROOT files passed the requested checks.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
