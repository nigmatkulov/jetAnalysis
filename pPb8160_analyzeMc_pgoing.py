#!/usr/bin/env python3
"""Batch runner for the compiled processForestSimple executable."""

from __future__ import annotations

import sys

from pPb8160_runner import RunnerConfig, main as run_batch


MC_TYPE = 1  # 0=data, 1=embedding, 2=pythia
IS_PB_GOING_DIR = 0
TRIGGER_ID = 0  
RECO_JET_SEL_METHOD = 2
JER_SMEARING_MODE = 0  # 0=CMS hybrid, 1=forced stochastic, 2=truth-to-pseudo-reco
USE_ETA_DEPENDENT_JER = 0  # 0=base standard + eta-dependent Extra; 1=eta-dependent in both
PT_HAT_SAMPLES = [15, 30, 50, 80, 120, 170, 220, 280, 370, 460, 540]

CONFIG = RunnerConfig(
    mc_type=MC_TYPE,
    is_pb_going_dir=IS_PB_GOING_DIR,
    trigger_id=TRIGGER_ID,
    reco_jet_sel_method=RECO_JET_SEL_METHOD,
    jer_smearing_mode=JER_SMEARING_MODE,
    use_eta_dependent_jer=USE_ETA_DEPENDENT_JER,
)


def main() -> int:
    return run_batch(CONFIG, PT_HAT_SAMPLES, sys.argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
