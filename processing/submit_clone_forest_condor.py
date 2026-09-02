#!/usr/bin/env python3
"""Prepare one HTCondor ``cloneForest.C`` job per input ROOT file."""

from __future__ import annotations

import argparse
import datetime as dt
import re
import shutil
import subprocess
import uuid
from pathlib import Path, PurePosixPath


REPOSITORY = Path(__file__).resolve().parents[1]
PROCESSING_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT_LIST = (
    PROCESSING_DIR / "filelists/pPb8160/PYTHIA_CP5_Tune/PYTHIA8160.txt"
)
DEFAULT_OUTPUT_DIR = Path("/eos/user/g/gnigmatk/ana/pPb8160/pythia/cp5tune")
DEFAULT_MACRO = REPOSITORY / "macro/cloneForest.C"
WORKER = PROCESSING_DIR / "run_clone_forest_condor.sh"
PT_HAT_SAMPLES = (15, 30, 50, 80, 120, 170, 220, 280, 370, 460, 540)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Submit one cloneForest.C job per ROOT file in a text list."
    )
    parser.add_argument("input_list", nargs="?", type=Path)
    parser.add_argument("output_dir", nargs="?", type=Path)
    parser.add_argument(
        "--sample", choices=("cp5", "unembedded", "embedded"), default="cp5",
        help="Built-in file-list/output preset (default: cp5)",
    )
    parser.add_argument(
        "--beam", choices=("Pbgoing", "pgoing", "all"), default="Pbgoing",
        help="Beam direction for embedded/unembedded presets (default: Pbgoing)",
    )
    parser.add_argument(
        "--pt-hat", choices=("all", *(str(value) for value in PT_HAT_SAMPLES)),
        default="all", help="pT-hat list for embedded/unembedded presets (default: all)",
    )
    parser.add_argument("--macro", type=Path, default=DEFAULT_MACRO)
    parser.add_argument("--name")
    parser.add_argument("--job-flavour", default="microcentury")
    parser.add_argument("--request-memory", default="2GB")
    parser.add_argument(
        "--requirements",
        default='(OpSysAndVer =?= "AlmaLinux9") && (CERNEnvironment =?= "qa")',
    )
    parser.add_argument("--proxy", type=Path, default=PROCESSING_DIR / "voms_proxy.txt")
    parser.add_argument("--work-dir", type=Path, default=PROCESSING_DIR / "condor")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def condor_token(value: str | Path) -> str:
    text = str(value)
    if any(character.isspace() for character in text):
        raise SystemExit(f"whitespace is unsupported in Condor paths: {text!r}")
    if any(character in text for character in ('"', "'", "\\")):
        raise SystemExit(f"quotes and backslashes are unsupported in paths: {text!r}")
    return text


def read_inputs(path: Path) -> list[str]:
    inputs = []
    for line_number, raw_line in enumerate(path.read_text().splitlines(), start=1):
        entry = raw_line.strip()
        if not entry or entry.startswith("#"):
            continue
        if not entry.endswith(".root"):
            raise SystemExit(f"{path}:{line_number}: expected a .root path: {entry}")
        condor_token(entry)
        if entry.startswith("/store/"):
            entry = f"root://cms-xrd-global.cern.ch/{entry}"
        inputs.append(entry)
    if not inputs:
        raise SystemExit(f"input list has no ROOT files: {path}")
    return inputs


def preset_jobs(args: argparse.Namespace) -> list[tuple[str, Path]]:
    """Return ``(input file, output directory)`` pairs for a preset."""

    if args.sample == "cp5":
        return [(item, DEFAULT_OUTPUT_DIR) for item in read_inputs(DEFAULT_INPUT_LIST)]

    beams = ("Pbgoing", "pgoing") if args.beam == "all" else (args.beam,)
    pt_hats = PT_HAT_SAMPLES if args.pt_hat == "all" else (int(args.pt_hat),)
    list_kind = "MC_Unembedded" if args.sample == "unembedded" else "MC_embedded"
    output_kind = "pythia" if args.sample == "unembedded" else "embedding"
    jobs = []
    for beam in beams:
        for pt_hat in pt_hats:
            input_list = (
                PROCESSING_DIR / "filelists/pPb8160" / list_kind / beam
                / f"MC_pthat{pt_hat}_{beam}_{args.sample}.txt"
            )
            if not input_list.is_file():
                raise SystemExit(f"input list does not exist: {input_list}")
            output_dir = Path(
                f"/eos/user/g/gnigmatk/ana/pPb8160/{output_kind}/{beam}/forest"
            )
            jobs.extend((item, output_dir) for item in read_inputs(input_list))
    return jobs


def main() -> int:
    args = parse_args()
    macro = args.macro.expanduser().resolve()
    work_dir = args.work_dir.expanduser().resolve()
    proxy = args.proxy.expanduser().resolve()

    if (args.input_list is None) != (args.output_dir is None):
        raise SystemExit("input_list and output_dir must be provided together")
    if not macro.is_file():
        raise SystemExit(f"macro does not exist: {macro}")
    if not WORKER.is_file():
        raise SystemExit(f"worker script does not exist: {WORKER}")
    if not proxy.is_file():
        raise SystemExit(f"proxy does not exist: {proxy}")
    name = args.name or (
        "PYTHIA8160_CP5" if args.input_list is None and args.sample == "cp5"
        else f"MC_{args.sample}_{args.beam}"
        if args.input_list is None else "cloneForest"
    )
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
        raise SystemExit("--name may contain only letters, numbers, '.', '_' and '-'")
    if "\n" in args.requirements:
        raise SystemExit("--requirements must be a single line")

    if args.input_list is not None:
        input_list = args.input_list.expanduser().resolve()
        if not input_list.is_file():
            raise SystemExit(f"input list does not exist: {input_list}")
        jobs = [
            (item, args.output_dir.expanduser()) for item in read_inputs(input_list)
        ]
    else:
        jobs = preset_jobs(args)

    outputs = [str(output / PurePosixPath(item).name) for item, output in jobs]
    duplicate_outputs = sorted({item for item in outputs if outputs.count(item) > 1})
    if duplicate_outputs:
        raise SystemExit("duplicate output paths: " + ", ".join(duplicate_outputs))

    campaign_id = (
        dt.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        + f"_{uuid.uuid4().hex[:8]}"
    )
    campaign_dir = work_dir / campaign_id / name
    logs_dir = campaign_dir / "log"
    logs_dir.mkdir(parents=True)
    shutil.copy2(proxy, campaign_dir / "voms_proxy.txt")

    rows = []
    for job_number, (input_file, output_dir) in enumerate(jobs, start=1):
        job_name = f"{name}_job{job_number:04d}"
        rows.append(
            f"{job_name} {condor_token(input_file)} {condor_token(output_dir)}"
        )

    submit_file = campaign_dir / f"{name}_{campaign_id}.sub"
    submit_lines = [
        "universe = vanilla",
        "initialdir = .",
        f"executable = {condor_token(WORKER)}",
        f'+JobFlavour = "{args.job_flavour}"',
        "getenv = True",
        "request_cpus = 1",
        f"request_memory = {args.request_memory}",
        f"requirements = {args.requirements}",
        "transfer_input_files = voms_proxy.txt",
        'environment = "X509_USER_PROXY=voms_proxy.txt"',
        "arguments = " + " ".join(
            (condor_token(macro), "$(input_file)", "$(output_dir)")
        ),
        f"output = {condor_token(logs_dir)}/$(job_name).out",
        f"error = {condor_token(logs_dir)}/$(job_name).err",
        f"log = {condor_token(logs_dir)}/$(job_name).log",
        "queue job_name, input_file, output_dir from (",
        *rows,
        ")",
        "",
    ]
    submit_file.write_text("\n".join(submit_lines))
    print(f"Prepared {len(jobs)} job(s): {submit_file}")

    if args.dry_run:
        print("Dry run: jobs were not submitted.")
        return 0
    if shutil.which("condor_submit") is None:
        raise SystemExit("condor_submit was not found; rerun with --dry-run")
    return subprocess.run(
        ["condor_submit", submit_file.name], cwd=campaign_dir, check=False
    ).returncode


if __name__ == "__main__":
    raise SystemExit(main())
