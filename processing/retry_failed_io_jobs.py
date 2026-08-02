#!/usr/bin/env python3
"""Find completed Condor jobs with ROOT read errors and retry only those jobs."""

from __future__ import annotations

import argparse
import datetime as dt
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


PROCESSING_DIR = Path(__file__).resolve().parent
DEFAULT_CONDOR_DIR = PROCESSING_DIR / "condor"
QUEUE_HEADER = "queue job_name, input_list, output_file from ("
READ_ERROR = re.compile(
    r"(?:"
    r"SysError in <TFile::ReadBuffer>:.*Input/output error"
    r"|Error in <TBranch::GetBasket>:.*badread="
    r"|Failed to read chain entry .* from tree .* file "
    r")"
)


@dataclass(frozen=True)
class QueueRow:
    """One generated Condor queue row."""

    job_name: str
    input_list: Path
    output_file: Path
    source_line: str


@dataclass(frozen=True)
class FailedJob:
    """One completed job selected for retry."""

    row: QueueRow
    stdout: Path
    stderr: Path
    event_log: Path


def parse_args() -> argparse.Namespace:
    """Parse recovery options."""

    parser = argparse.ArgumentParser(
        description=(
            "Find completed processForestSimple jobs with ROOT input read "
            "errors and generate retries containing only those jobs."
        )
    )
    parser.add_argument(
        "condor_dir",
        nargs="?",
        type=Path,
        default=DEFAULT_CONDOR_DIR,
        help=f"Condor campaign root (default: {DEFAULT_CONDOR_DIR})",
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help=(
            "Archive failed-job logs, quarantine stale outputs, and call "
            "condor_submit. Without this flag, report a dry run only."
        ),
    )
    return parser.parse_args()


def original_submit_files(condor_dir: Path) -> list[Path]:
    """Return original generated submit files, excluding prior retries."""

    return sorted(
        path
        for path in condor_dir.rglob("*.sub")
        if "_retry_" not in path.stem
    )


def parse_queue(submit_file: Path) -> tuple[list[str], int, int, list[QueueRow]]:
    """Parse the queue block and return its lines, bounds, and rows."""

    lines = submit_file.read_text().splitlines()
    try:
        queue_start = next(
            index
            for index, line in enumerate(lines)
            if line.strip().lower() == QUEUE_HEADER
        )
    except StopIteration as error:
        raise ValueError(f"Queue block not found in {submit_file}") from error

    try:
        queue_end = next(
            index
            for index in range(queue_start + 1, len(lines))
            if lines[index].strip() == ")"
        )
    except StopIteration as error:
        raise ValueError(f"Queue block is not closed in {submit_file}") from error

    rows: list[QueueRow] = []
    for line in lines[queue_start + 1 : queue_end]:
        fields = line.split()
        if not fields:
            continue
        if len(fields) != 3:
            raise ValueError(f"Malformed queue row in {submit_file}: {line!r}")
        rows.append(
            QueueRow(
                job_name=fields[0],
                input_list=Path(fields[1]),
                output_file=Path(fields[2]),
                source_line=line,
            )
        )
    return lines, queue_start, queue_end, rows


def contains_read_error(stderr: Path) -> bool:
    """Return whether stderr contains one targeted ROOT input-read failure."""

    with stderr.open(errors="replace") as stream:
        return any(READ_ERROR.search(line) for line in stream)


def job_has_terminated(event_log: Path) -> bool:
    """Return whether the latest relevant Condor lifecycle event is termination."""

    terminated = False
    with event_log.open(errors="replace") as stream:
        for line in stream:
            if "Job executing on host" in line:
                terminated = False
            elif "Job terminated." in line:
                terminated = True
    return terminated


def failed_jobs(submit_file: Path, rows: list[QueueRow]) -> list[FailedJob]:
    """Select completed jobs whose stderr contains a targeted read error."""

    logs_dir = submit_file.parent / "log"
    selected: list[FailedJob] = []
    for row in rows:
        stdout = logs_dir / f"{row.job_name}.out"
        stderr = logs_dir / f"{row.job_name}.err"
        event_log = logs_dir / f"{row.job_name}.log"
        if not stderr.is_file() or not event_log.is_file():
            continue
        if not contains_read_error(stderr):
            continue
        if not job_has_terminated(event_log):
            print(f"Skipping active or unterminated job: {row.job_name}")
            continue
        selected.append(
            FailedJob(
                row=row,
                stdout=stdout,
                stderr=stderr,
                event_log=event_log,
            )
        )
    return selected


def write_retry_submit(
    submit_file: Path,
    lines: list[str],
    queue_start: int,
    queue_end: int,
    jobs: list[FailedJob],
    retry_tag: str,
) -> Path:
    """Write a submit description containing only selected failed jobs."""

    retry_file = submit_file.with_name(f"{submit_file.stem}_retry_{retry_tag}.sub")
    retry_lines = [
        *lines[: queue_start + 1],
        *(job.row.source_line for job in jobs),
        *lines[queue_end:],
    ]
    retry_file.write_text("\n".join(retry_lines) + "\n")
    return retry_file


def clear_failed_artifacts(
    submit_file: Path, jobs: list[FailedJob], retry_tag: str,
) -> None:
    """Archive failed logs and quarantine any stale output files."""

    archive_dir = submit_file.parent / "retry_archive" / retry_tag
    archive_dir.mkdir(parents=True, exist_ok=False)
    for job in jobs:
        job_archive = archive_dir / job.row.job_name
        job_archive.mkdir()
        for log_file in (job.stdout, job.stderr, job.event_log):
            if log_file.exists():
                shutil.move(str(log_file), job_archive / log_file.name)

        output = job.row.output_file
        if output.exists():
            quarantined = output.with_name(f"{output.name}.invalid-{retry_tag}")
            output.replace(quarantined)
            print(f"Quarantined stale output: {quarantined}")


def main() -> int:
    """Scan campaigns, prepare retry files, and optionally submit them."""

    args = parse_args()
    condor_dir = args.condor_dir.expanduser().resolve()
    if not condor_dir.is_dir():
        raise SystemExit(f"Condor directory does not exist: {condor_dir}")
    if args.submit and shutil.which("condor_submit") is None:
        raise SystemExit("condor_submit was not found; omit --submit for a dry run")

    submit_files = original_submit_files(condor_dir)
    if not submit_files:
        raise SystemExit(f"No original submit files found below {condor_dir}")

    plans: list[tuple[Path, list[str], int, int, list[FailedJob]]] = []
    total_failed = 0
    for submit_file in submit_files:
        lines, queue_start, queue_end, rows = parse_queue(submit_file)
        jobs = failed_jobs(submit_file, rows)
        if not jobs:
            continue
        plans.append((submit_file, lines, queue_start, queue_end, jobs))
        total_failed += len(jobs)
        print(f"{submit_file}: {len(jobs)} failed job(s)")
        for job in jobs:
            print(f"  {job.row.job_name}: {job.row.output_file}")

    if total_failed == 0:
        print("No completed jobs with matching ROOT input read errors were found.")
        return 0
    if not args.submit:
        print(f"Dry run: {total_failed} failed job(s) would be retried.")
        return 0

    retry_tag = dt.datetime.now().strftime("%Y%m%d_%H%M%S%f")
    for submit_file, lines, queue_start, queue_end, jobs in plans:
        retry_file = write_retry_submit(
            submit_file, lines, queue_start, queue_end, jobs, retry_tag,
        )
        clear_failed_artifacts(submit_file, jobs, retry_tag)
        subprocess.run(
            ["condor_submit", retry_file.name],
            cwd=submit_file.parent,
            check=True,
        )
        print(f"Submitted {len(jobs)} retry job(s): {retry_file}")

    print(f"Submitted {total_failed} failed job(s) in total.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
