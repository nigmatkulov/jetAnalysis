#!/usr/bin/env python3
"""Shared helpers for the pPb8160 batch runner scripts."""

from __future__ import annotations

import argparse
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT_BASE = Path.home() / "cernbox" / "ana" / "pPb8160"
DEFAULT_OUTPUT_DIR = Path.home() / "cernbox" / "ana" / "pPb8160"


@dataclass(frozen=True)
class RunnerConfig:
    mc_type: int
    is_pb_going_dir: int
    trigger_id: int
    reco_jet_sel_method: int
    input_base: Path = DEFAULT_INPUT_BASE
    output_dir: Path = DEFAULT_OUTPUT_DIR


@dataclass(frozen=True)
class SampleResult:
    pt_hat_sample: int
    status: str
    returncode: int
    input_file: Path
    output_file: Path


def _generator_name(mc_type: int) -> str:
    return "pythia" if mc_type == 2 else "embedding"


def _tag_name(mc_type: int) -> str:
    return "unembedded" if mc_type == 2 else "embedded"


def _direction_name(is_pb_going_dir: int) -> str:
    return "Pbgoing" if is_pb_going_dir == 1 else "pgoing"


def build_input_file(config: RunnerConfig, pt_hat_sample: int) -> Path:
    generator = _generator_name(config.mc_type)
    direction = _direction_name(config.is_pb_going_dir)
    tag = _tag_name(config.mc_type)
    return (
        config.input_base
        / generator
        / direction
        / "forest"
        / f"HiForestSkim_pPb_MC_pthat{pt_hat_sample}_{direction}_{tag}.root"
    )


def build_output_file(config: RunnerConfig, pt_hat_sample: int) -> Path:
    generator = _generator_name(config.mc_type)
    direction = _direction_name(config.is_pb_going_dir)
    return config.output_dir / f"{generator}" / f"{direction}" / f"{generator}_{direction}_ptHat{pt_hat_sample}.root"


def build_command(config: RunnerConfig, pt_hat_sample: int) -> list[str]:
    input_file = build_input_file(config, pt_hat_sample)
    output_file = build_output_file(config, pt_hat_sample)
    return [
        str(ROOT_DIR / "build" / "processForestSimple"),
        str(input_file),
        str(output_file),
        str(config.mc_type),
        str(config.is_pb_going_dir),
        str(pt_hat_sample),
        str(config.trigger_id),
        str(config.reco_jet_sel_method),
    ]


def run_sample_task(task: tuple[RunnerConfig, int]) -> SampleResult:
    config, pt_hat_sample = task
    input_file = build_input_file(config, pt_hat_sample)
    output_file = build_output_file(config, pt_hat_sample)

    if not input_file.exists():
        return SampleResult(
            pt_hat_sample=pt_hat_sample,
            status="skipped",
            returncode=0,
            input_file=input_file,
            output_file=output_file,
        )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    command = build_command(config, pt_hat_sample)

    print(f"\nRunning ptHat {pt_hat_sample}:")
    print(" ".join(command))
    completed = subprocess.run(command, cwd=ROOT_DIR)
    status = "succeeded" if completed.returncode == 0 else "failed"

    return SampleResult(
        pt_hat_sample=pt_hat_sample,
        status=status,
        returncode=completed.returncode,
        input_file=input_file,
        output_file=output_file,
    )


def _clamp_workers(requested_workers: int) -> int:
    cpu_count = os.cpu_count() or 1
    return max(1, min(requested_workers, cpu_count))


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run processForestSimple over pPb8160 pt-hat samples."
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of pt-hat samples to run concurrently (default: 1).",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    return args


def run_samples(
    config: RunnerConfig,
    pt_hat_samples: list[int],
    workers: int,
) -> list[SampleResult]:
    if workers == 1:
        return [run_sample_task((config, pt_hat_sample)) for pt_hat_sample in pt_hat_samples]

    tasks = [(config, pt_hat_sample) for pt_hat_sample in pt_hat_samples]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(run_sample_task, tasks))


def print_summary(results: list[SampleResult]) -> int:
    skipped = [result for result in results if result.status == "skipped"]
    failed = [result for result in results if result.status == "failed"]

    if skipped:
        print("\nSkipped samples:")
        for result in skipped:
            print(f"  ptHat {result.pt_hat_sample}: missing input file")
            print(f"    {result.input_file}")

    if failed:
        print("\nFailed samples:")
        for result in failed:
            print(f"  ptHat {result.pt_hat_sample}: return code {result.returncode}")
            print(f"    input:  {result.input_file}")
            print(f"    output: {result.output_file}")
        return 1

    print("\nAll requested ptHat samples completed successfully.")
    return 0


def main(config: RunnerConfig, pt_hat_samples: list[int], argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    workers = _clamp_workers(args.workers)
    if workers != args.workers:
        print(f"Using {workers} worker(s) based on available CPU cores.")

    results = run_samples(config, pt_hat_samples, workers)
    return print_summary(results)
