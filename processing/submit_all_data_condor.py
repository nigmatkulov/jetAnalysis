#!/usr/bin/env python3
"""Submit all pPb data trigger samples in both beam orientations."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


PROCESSING_DIR = Path(__file__).resolve().parent
SUBMITTER = PROCESSING_DIR / "submit_process_forest_condor.py"
FILELISTS = PROCESSING_DIR / "filelists" / "pPb8160"
DEFAULT_OUTPUT_BASE = Path("/eos/user/g/gnigmatk/ana/pPb8160/exp")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Submit MB, Jet60, Jet80, and Jet100 data for both beam orientations."
    )
    parser.add_argument(
        "--output-base", type=Path, default=DEFAULT_OUTPUT_BASE,
        help=f"Parent of Pbgoing/ and pgoing/ output directories (default: {DEFAULT_OUTPUT_BASE})",
    )
    parser.add_argument("--files-per-job", type=int, default=10)
    parser.add_argument("--reco-jet-selection", type=int, default=2)
    parser.add_argument("--job-flavour", help="Override the default espresso flavour")
    parser.add_argument("--proxy", type=Path, help="Forward a non-default proxy path")
    parser.add_argument("--work-dir", type=Path, help="Forward a non-default Condor work directory")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def submission_command(args: argparse.Namespace, beam: str, trigger_id: int) -> list[str]:
    if trigger_id == 0:
        input_list = FILELISTS / "DATA_MB" / beam / f"MB_PD{{pd}}_{beam}.txt"
    else:
        input_list = FILELISTS / "DATA_PAEGJet" / beam / f"PAEGJet_{beam}.txt"

    command = [
        sys.executable,
        str(SUBMITTER),
        str(input_list),
        str(args.output_base / beam),
        "--mc-type", "0",
        "--beam", beam,
        "--trigger-id", str(trigger_id),
        "--files-per-job", str(args.files_per_job),
        "--reco-jet-selection", str(args.reco_jet_selection),
    ]
    if trigger_id == 0:
        command.extend(("--pd-number", "all"))
    if args.job_flavour:
        command.extend(("--job-flavour", args.job_flavour))
    if args.proxy:
        command.extend(("--proxy", str(args.proxy)))
    if args.work_dir:
        command.extend(("--work-dir", str(args.work_dir)))
    if args.dry_run:
        command.append("--dry-run")
    return command


def main() -> int:
    args = parse_args()
    if args.files_per_job < 1:
        raise SystemExit("--files-per-job must be at least 1")
    if not 0 <= args.reco_jet_selection <= 3:
        raise SystemExit("--reco-jet-selection must be between 0 and 3")

    for beam in ("Pbgoing", "pgoing"):
        for trigger_id in range(4):
            completed = subprocess.run(
                submission_command(args, beam, trigger_id), check=False
            )
            if completed.returncode != 0:
                return completed.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
