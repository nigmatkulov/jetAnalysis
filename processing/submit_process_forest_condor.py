#!/usr/bin/env python3
"""Prepare and submit ``processForestSimple`` jobs to HTCondor.

The analysis executable accepts one ROOT file *or* one text file containing a
list of ROOT files.  Processing a large production list in a single invocation
is inconvenient for batch execution, so this script performs the orchestration
needed by HTCondor:

1. Read a user-supplied master file list, ignoring blank lines and comments.
2. Split the entries into deterministic, numbered sublists.
3. Assign one output ROOT file to every sublist.
4. Write one Condor submit description whose queue table contains all jobs.
5. Call ``condor_submit`` unless ``--dry-run`` was requested.

The physics-facing command-line options deliberately mirror the current
``processForestSimple.C`` executable interface.  In particular, this script
does not expose the JEU, JER, pT-hat-upper-bound, or event-limit arguments used
by older analysis executables.

Generated files are grouped into a uniquely named *campaign directory*::

    <work-dir>/<YYYYMMDD_HHMMSSffffff_uuid>/<campaign-label>/
        <campaign-label>_<campaign-id>.sub
        input/<campaign-label>_job0001.list
        log/<campaign-label>_job0001.{out,err,log}

Submit files and campaign directories include a microsecond timestamp plus a
random UUID suffix, so simultaneous submissions cannot collide.  ROOT outputs
use stable physics labels and job numbers; rerunning an identical configuration
therefore replaces its earlier output, while different triggers, directions,
PDs, pT-hat samples, and selections remain distinct.

The repository, analysis executable, generated sublists, and output directory
are passed to jobs as absolute paths.  They therefore need to be visible from
the worker nodes, for example through CERN AFS or EOS.  An X.509 proxy is copied
to each sample directory and transferred by Condor; the analysis installation
itself is not packaged or transferred.
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
import shutil
import subprocess
import uuid
from pathlib import Path


# Derive installation paths from this file instead of the caller's working
# directory.  This permits invocation from either the repository root or an
# arbitrary directory and ensures the worker starts in the correct directory
# for relative resources such as aux_files/pPb_8160.
REPOSITORY = Path(__file__).resolve().parents[1]
PROCESSING_DIR = Path(__file__).resolve().parent
DEFAULT_EXECUTABLE = REPOSITORY / "build" / "processForestSimple"
WORKER = PROCESSING_DIR / "run_process_forest_condor.sh"
PT_HAT_SAMPLES = (15, 30, 50, 80, 120, 170, 220, 280, 370, 460, 540)


def bounded_int(name: str, minimum: int, maximum: int):
    """Build an argparse converter for an integer in an inclusive range.

    ``choices=range(...)`` would enforce the same constraint, but this helper
    gives users an error that names the offending option and its valid bounds.
    It is shared by the enumerated integer arguments defined by the C++
    executable: MC type, trigger ID, and reco-jet selection method.
    """

    def convert(value: str) -> int:
        """Parse and range-check one command-line value."""

        parsed = int(value)
        if not minimum <= parsed <= maximum:
            raise argparse.ArgumentTypeError(
                f"{name} must be between {minimum} and {maximum}"
            )
        return parsed

    return convert


def positive_int(value: str) -> int:
    """Parse an integer that must be positive, for ``--files-per-job``."""

    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("value must be at least 1")
    return parsed


def positive_int_or_all(value: str) -> int | str:
    """Parse a positive integer or the case-insensitive keyword ``all``."""

    if value.lower() == "all":
        return "all"
    return positive_int(value)


def parse_args() -> argparse.Namespace:
    """Define and parse all submission and analysis command-line options.

    The two positional paths identify what to process and where ROOT results
    should be written.  Options through ``--reco-jet-selection`` map directly
    to ``processForestSimple``.  The remaining options control job splitting,
    naming, Condor resource requests, credentials, and submission behavior.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Split a ROOT-file list and create HTCondor jobs for "
            "processForestSimple."
        )
    )
    parser.add_argument(
        "input_list",
        help=(
            "Input-list path; use {pd} with --pd-number all or {pt_hat} with "
            "--pt-hat all"
        ),
    )
    parser.add_argument("output_dir", type=Path, help="Directory for output ROOT files")
    parser.add_argument(
        "--mc-type",
        type=bounded_int("--mc-type", 0, 2),
        default=0,
        help="0=data, 1=embedding, 2=pythia (default: 0)",
    )
    parser.add_argument(
        "--beam",
        choices=("pgoing", "Pbgoing"),
        default="Pbgoing",
        help="Beam orientation (default: Pbgoing)",
    )
    parser.add_argument(
        "--pt-hat",
        type=positive_int_or_all,
        default=0,
        help="MC pT-hat sample or 'all'; ignored for data (default: 0)",
    )
    parser.add_argument(
        "--trigger-id",
        type=bounded_int("--trigger-id", 0, 3),
        default=0,
        help="0=MB, 1=Jet60, 2=Jet80, 3=Jet100 (default: 0)",
    )
    parser.add_argument(
        "--pd-number",
        type=positive_int_or_all,
        help="MB primary-dataset number or 'all'; ignored for other triggers",
    )
    parser.add_argument(
        "--reco-jet-selection",
        type=bounded_int("--reco-jet-selection", 0, 3),
        default=2,
        help="0=none, 1=trkMaxPt/RawPt, 2=tight ID+lepton veto, 3=loose ID",
    )
    parser.add_argument(
        "--files-per-job", type=positive_int, default=50, help="Default: 50"
    )
    parser.add_argument(
        "--name",
        help="Optional prefix placed before the mandatory readable physics labels",
    )
    parser.add_argument(
        "--executable", type=Path, default=DEFAULT_EXECUTABLE,
        help=f"Analysis executable (default: {DEFAULT_EXECUTABLE})",
    )
    parser.add_argument(
        "--job-flavour", default="longlunch", help="CERN JobFlavour (default: longlunch)"
    )
    parser.add_argument(
        "--request-memory", default="2GB", help="Condor memory request (default: 2GB)"
    )
    parser.add_argument(
        "--requirements",
        default='(OpSysAndVer =?= "AlmaLinux9") && (CERNEnvironment =?= "qa")',
        help="Condor requirements expression (default: CERN AlmaLinux 9 worker)",
    )
    parser.add_argument(
        "--proxy",
        type=Path,
        default=PROCESSING_DIR / "voms_proxy.txt",
        help="X.509 proxy (default: processing/voms_proxy.txt)",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=PROCESSING_DIR / "condor",
        help="Directory for generated inputs and Condor files",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Generate files but do not call condor_submit"
    )
    return parser.parse_args()


def condor_quote(value: Path | str) -> str:
    """Quote a literal for use in a Condor submit description.

    Condor submit syntax is not shell syntax.  Double quotes preserve spaces in
    paths placed in the inline queue table.  Embedded quotes and newlines would
    require Condor-specific escaping and could alter the generated description,
    so reject them explicitly instead of producing an ambiguous submit file.

    Condor macros such as ``$(input_list)`` are intentionally *not* passed
    through this function where the ``arguments`` line is assembled; they must
    remain available for expansion for each queued job.
    """

    text = str(value)
    if '"' in text or "\n" in text:
        raise ValueError(f"unsupported character in Condor argument: {text!r}")
    return f'"{text}"'


def make_name(
    args: argparse.Namespace, pd_number: int | None, pt_hat: int
) -> str:
    """Return a readable, filesystem-safe physics label.

    Trigger and direction are always present.  MB data additionally carries the
    required primary-dataset number, while MC carries its type and pT-hat bin.
    A custom ``--name`` is only a prefix and therefore cannot hide the mandatory
    identifiers.  Restricting it to a small safe character set also prevents
    whitespace from complicating Condor queue and log paths.
    """

    if args.name and not re.fullmatch(r"[A-Za-z0-9_.-]+", args.name):
        raise ValueError("--name may contain only letters, numbers, '.', '_' and '-'")

    kind = ("data", "embedding", "pythia")[args.mc_type]
    trigger = ("MB", "Jet60", "Jet80", "Jet100")[args.trigger_id]
    reco_selection = ("noSel", "trkMax", "jetId", "jetId")[
        args.reco_jet_selection
    ]
    if args.mc_type != 0 and args.trigger_id == 0:
        trigger = "NoTrigger"
    labels = [trigger]
    if args.mc_type == 0 and args.trigger_id == 0:
        labels.append(f"PD{pd_number}")
    labels.extend((args.beam, kind))
    if args.mc_type != 0:
        labels.append(f"ptHat{pt_hat}")
    labels.append(reco_selection)
    if args.name:
        labels.insert(0, args.name)
    return "_".join(labels)


def selected_pd_numbers(args: argparse.Namespace) -> list[int | None]:
    """Expand the MB PD selection and enforce direction-dependent boundaries.

    PD numbering is an MB-only data convention.  The option is deliberately
    ignored for Jet60/80/100 and MC so it cannot multiply those submissions.
    """

    if args.mc_type != 0 or args.trigger_id != 0:
        return [None]
    if args.pd_number is None:
        raise SystemExit("--pd-number is required when submitting MB data")
    maximum = 20 if args.beam == "Pbgoing" else 8
    if args.pd_number == "all":
        return list(range(1, maximum + 1))
    if not 1 <= args.pd_number <= maximum:
        raise SystemExit(
            f"MB {args.beam} PD number must be between 1 and {maximum}, inclusive"
        )
    return [args.pd_number]


def selected_pt_hat_samples(args: argparse.Namespace) -> list[int]:
    """Expand ``--pt-hat all`` for MC while using a harmless zero for data."""

    if args.mc_type == 0:
        return [0]
    if args.pt_hat == "all":
        return list(PT_HAT_SAMPLES)
    if args.pt_hat not in PT_HAT_SAMPLES:
        supported = ", ".join(str(value) for value in PT_HAT_SAMPLES)
        raise SystemExit(f"MC --pt-hat must be 'all' or one of: {supported}")
    return [args.pt_hat]


def resolve_input_list(
    template: str, beam: str, pd_number: int | None, pt_hat: int
) -> Path:
    """Resolve one input-list template for a concrete PD or pT-hat sample."""

    try:
        rendered = template.format(pd=pd_number, pt_hat=pt_hat, beam=beam)
    except KeyError as error:
        raise SystemExit(f"Unknown input-list placeholder: {error.args[0]}") from error
    return Path(rendered).expanduser().resolve()


def make_campaign_id() -> str:
    """Create a readable ID that remains unique across parallel invocations.

    Microseconds make nearby submissions easy to order by eye.  The UUID suffix
    supplies collision resistance even if separate processes or hosts obtain
    the same clock value.
    """

    timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    return f"{timestamp}_{uuid.uuid4().hex[:8]}"


def read_inputs(path: Path) -> list[str]:
    """Read usable ROOT-file entries while preserving their original order.

    Blank lines and lines whose first non-whitespace character is ``#`` are
    ignored.  Inline comments are not stripped because URL-like input paths may
    legitimately contain characters that should be passed to ROOT unchanged.
    """

    lines = [line.strip() for line in path.read_text().splitlines()]
    return [line for line in lines if line and not line.startswith("#")]


def main() -> int:
    """Prepare every requested PD/pT-hat sample and submit them sequentially."""

    args = parse_args()
    executable = args.executable.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    work_dir = args.work_dir.expanduser().resolve()
    proxy = args.proxy.expanduser().resolve()

    if not executable.is_file():
        raise SystemExit(f"Executable does not exist: {executable}")
    if not proxy.is_file():
        raise SystemExit(f"Proxy does not exist: {proxy}")
    if "\n" in args.requirements:
        raise SystemExit("--requirements must be a single line")

    pd_numbers = selected_pd_numbers(args)
    pt_hat_samples = selected_pt_hat_samples(args)
    if len(pd_numbers) > 1 and "{pd}" not in args.input_list:
        raise SystemExit("MB --pd-number all requires {pd} in input_list")
    if len(pt_hat_samples) > 1 and "{pt_hat}" not in args.input_list:
        raise SystemExit("MC --pt-hat all requires {pt_hat} in input_list")

    campaign_id = make_campaign_id()
    output_dir.mkdir(parents=True, exist_ok=True)
    submit_files: list[Path] = []

    # Each concrete PD or pT-hat sample gets its own submit description.  The
    # outer loop later invokes condor_submit once per description, in order.
    for pd_number in pd_numbers:
        for pt_hat in pt_hat_samples:
            input_list = resolve_input_list(
                args.input_list, args.beam, pd_number, pt_hat
            )
            if not input_list.is_file():
                raise SystemExit(f"Input list does not exist: {input_list}")
            inputs = read_inputs(input_list)
            if not inputs:
                raise SystemExit(f"Input list has no usable entries: {input_list}")

            label = make_name(args, pd_number, pt_hat)
            sample_dir = work_dir / campaign_id / label
            lists_dir = sample_dir / "input"
            logs_dir = sample_dir / "log"
            lists_dir.mkdir(parents=True)
            logs_dir.mkdir(parents=True)

            # Condor must see these exact credential directives.  The submit
            # command runs with cwd=sample_dir, so initialdir="." and the
            # relative voms_proxy.txt both resolve inside this sample campaign.
            shutil.copy2(proxy, sample_dir / "voms_proxy.txt")

            jobs: list[tuple[str, Path, Path]] = []
            for offset in range(0, len(inputs), args.files_per_job):
                job_id = len(jobs) + 1
                job_name = f"{label}_job{job_id:04d}"
                sublist = lists_dir / f"{job_name}.list"
                sublist.write_text(
                    "\n".join(inputs[offset : offset + args.files_per_job]) + "\n"
                )
                # Rerunning the same configuration intentionally reuses this
                # readable name.  Direction and trigger labels prevent the
                # cross-sample collisions the naming contract is designed for.
                output = output_dir / f"{job_name}.root"
                jobs.append((job_name, sublist, output))

            rows = [
                " ".join(
                    (
                        str(job_id),
                        job_name,
                        condor_quote(sublist),
                        condor_quote(output),
                    )
                )
                for job_id, (job_name, sublist, output) in enumerate(jobs, start=1)
            ]
            is_pb_going = 1 if args.beam == "Pbgoing" else 0
            submit_file = sample_dir / f"{label}_{campaign_id}.sub"
            submit_lines = [
                "universe = vanilla",
                "initialdir = .",
                f"executable = {WORKER}",
                f'+JobFlavour = "{args.job_flavour}"',
                "getenv = True",
                "request_cpus = 1",
                f"request_memory = {args.request_memory}",
                f"requirements = {args.requirements}",
                "transfer_input_files  = voms_proxy.txt",
                'environment = "X509_USER_PROXY=voms_proxy.txt"',
                "arguments = " + " ".join(
                    (
                        condor_quote(REPOSITORY),
                        condor_quote(executable),
                        "$(input_list)",
                        "$(output_file)",
                        str(args.mc_type),
                        str(is_pb_going),
                        str(pt_hat),
                        str(args.trigger_id),
                        str(args.reco_jet_selection),
                    )
                ),
                f"output = {logs_dir}/$(job_name).out",
                f"error = {logs_dir}/$(job_name).err",
                f"log = {logs_dir}/$(job_name).log",
                "queue job_id, job_name, input_list, output_file from (",
                *rows,
                ")",
                "",
            ]
            submit_file.write_text("\n".join(submit_lines))
            submit_files.append(submit_file)
            print(f"Prepared {len(jobs)} job(s) for {label}: {submit_file}")

    if args.dry_run:
        print(f"Dry run: {len(submit_files)} sample(s) were not submitted.")
        return 0
    if shutil.which("condor_submit") is None:
        raise SystemExit("condor_submit was not found; rerun with --dry-run to generate only")

    # Sequential calls make each MB PD or MC pT-hat sample an independent
    # Condor submission.  Stop immediately if any submission is rejected.
    for submit_file in submit_files:
        completed = subprocess.run(
            ["condor_submit", submit_file.name], cwd=submit_file.parent, check=False
        )
        if completed.returncode != 0:
            return completed.returncode
    return 0


if __name__ == "__main__":
    # Converting the return value to SystemExit makes success/failure visible to
    # callers while leaving ``main`` directly testable.
    raise SystemExit(main())
