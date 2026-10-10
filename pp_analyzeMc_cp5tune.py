#!/usr/bin/env python3
"""Process local pp CP5 forests, one output and log per input file."""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Sequence

ROOT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_DIR = Path.home() / "cernbox" / "ana" / "pPb8160" / "pythia" / "cp5tune"
DEFAULT_INPUT_DIR = DEFAULT_OUTPUT_DIR / "forest"
MC_TYPE = 2
IS_PB_GOING_DIR = 1  # Retains the existing correction configuration; pp eta is unflipped.
TRIGGER_ID = 0
RECO_JET_SEL_METHOD = 2
VERTEX_FILTER_SELECTION = 0


def build_command(executable: Path, input_file: Path, output_file: Path,
                  jer_smearing_mode: int, use_eta_dependent_jer: bool) -> list[str]:
    return [
        str(executable), str(input_file), str(output_file), str(MC_TYPE),
        str(IS_PB_GOING_DIR), "0",  # ptHatSample is unused in pp mode.
        str(TRIGGER_ID), str(RECO_JET_SEL_METHOD), str(VERTEX_FILTER_SELECTION), "pp",
        str(jer_smearing_mode), str(int(use_eta_dependent_jer)),
    ]


def run_file(command: list[str]) -> tuple[Path, int]:
    output_file = Path(command[2])
    log_file = output_file.with_suffix(".log")
    print(f"Processing {Path(command[1]).name} -> {output_file.name}", flush=True)
    with log_file.open("w") as log:
        log.write(shlex.join(command) + "\n")
        log.flush()
        result = subprocess.run(command, cwd=ROOT_DIR, stdout=log, stderr=subprocess.STDOUT)
    returncode = result.returncode
    if returncode == 0 and (not output_file.is_file() or output_file.stat().st_size == 0):
        print(f"Missing or empty output: {output_file}", flush=True)
        returncode = 1
    return output_file, returncode


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=1,
                        help="Concurrent executable processes (default: 1).")
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--executable", type=Path, default=ROOT_DIR / "build" / "processForestSimple")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without writing files.")
    parser.add_argument("--overwrite", action="store_true", help="Allow replacement of existing outputs/logs.")
    parser.add_argument("--jer-smearing-mode", type=int, choices=(0, 1, 2), default=0,
                        help="0=CMS hybrid, 1=forced stochastic, 2=truth-to-pseudo-reco.")
    parser.add_argument("--eta-dependent-jer", action="store_true",
                        help="Multiply the central JER by the analysis-derived eta shape.")
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    input_dir = args.input_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    executable = args.executable.expanduser().resolve()
    if not input_dir.is_dir():
        parser.error(f"Input directory does not exist: {input_dir}")
    inputs = sorted(path for path in input_dir.glob("*.root") if path.is_file())
    if not inputs:
        parser.error(f"No ROOT files found in {input_dir}")
    if input_dir == output_dir:
        parser.error("Input and output directories must differ")
    commands = [build_command(executable, path, output_dir / f"{path.stem}_hist.root",
                              args.jer_smearing_mode, args.eta_dependent_jer)
                for path in inputs]
    workers = min(args.workers, os.cpu_count() or 1, len(commands))
    print(f"Found {len(commands)} files; using {workers} worker(s).")
    if args.dry_run:
        for command in commands:
            print(shlex.join(command))
        return 0
    if not executable.is_file() or not os.access(executable, os.X_OK):
        parser.error(f"Executable is missing or not executable: {executable}; build it first")
    if not args.overwrite:
        for command in commands:
            output = Path(command[2])
            for path in (output, output.with_suffix(".log")):
                if path.exists():
                    parser.error(f"Output/log already exists: {path}; use --overwrite to replace it")
    output_dir.mkdir(parents=True, exist_ok=True)
    if workers == 1:
        results = [run_file(command) for command in commands]
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(run_file, commands))
    failed = [(output, code) for output, code in results if code != 0]
    print(f"Completed: {len(results) - len(failed)} succeeded, {len(failed)} failed.")
    for output, code in failed:
        print(f"  Failed (exit {code}): {output}; log: {output.with_suffix('.log')}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
