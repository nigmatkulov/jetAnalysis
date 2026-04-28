#!/usr/bin/env python3
"""Sequential batch runner for the compiled processForestSimple executable.

This script is intentionally simple:
- no command-line options
- fixed pPb/Pythia, p-going settings
- runs all ptHat samples one after another
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parent
EXECUTABLE = ROOT_DIR / "build" / "processForestSimple"

# Hardcoded analysis settings
MC_TYPE = 2  # 0=data, 1=embedding, 2=pythia
IS_PB_GOING_DIR = 1
JEU_SYST = 0
JER_SYST = -99
TRIGGER_ID = 0
RECO_JET_SEL_METHOD = 0

PT_HAT_SAMPLES = [15, 30, 50, 80, 120, 170, 220, 280, 370, 460, 540]
# PT_HAT_SAMPLES = [30]

INPUT_BASE = Path.home() / "cernbox" / "ana" / "pPb8160"
OUTPUT_DIR = ROOT_DIR / "macro" / "eta_shift"
GENERATOR = "pythia"
DIRECTION = "Pbgoing"
TAG = "unembedded"


def build_input_file(pt_hat_sample: int) -> Path:
    return (
        INPUT_BASE
        / GENERATOR
        / DIRECTION
        / "forest"
        / f"HiForestSkim_pPb_MC_pthat{pt_hat_sample}_{DIRECTION}_{TAG}.root"
    )


def build_output_file(pt_hat_sample: int) -> Path:
    return OUTPUT_DIR / f"{GENERATOR}_{DIRECTION}_ptHat{pt_hat_sample}.root"


def run_sample(pt_hat_sample: int) -> int:
    input_file = build_input_file(pt_hat_sample)
    output_file = build_output_file(pt_hat_sample)

    if not input_file.exists():
        print(f"Skipping ptHat {pt_hat_sample}: input file not found")
        print(f"  {input_file}")
        return 0

    output_file.parent.mkdir(parents=True, exist_ok=True)

    command = [
        str(EXECUTABLE),
        str(input_file),
        str(output_file),
        str(MC_TYPE),
        str(IS_PB_GOING_DIR),
        str(pt_hat_sample),
        str(JEU_SYST),
        str(JER_SYST),
        str(TRIGGER_ID),
        str(RECO_JET_SEL_METHOD),
    ]

    print("\nRunning:")
    print(" ".join(command))
    completed = subprocess.run(command, cwd=ROOT_DIR)
    return completed.returncode


def main() -> int:
    if not EXECUTABLE.exists():
        print(f"Error: executable not found at {EXECUTABLE}", file=sys.stderr)
        print("Build it first with: cmake --build build -j", file=sys.stderr)
        return 1

    for pt_hat_sample in PT_HAT_SAMPLES:
        return_code = run_sample(pt_hat_sample)
        if return_code != 0:
            print(f"Stopped at ptHat {pt_hat_sample} because the executable returned {return_code}")
            return return_code

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
