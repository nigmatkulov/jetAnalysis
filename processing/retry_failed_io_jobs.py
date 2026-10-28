#!/usr/bin/env python3
"""Audit expected processing outputs and prepare bounded, explicit recovery."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


PROCESSING_DIR = Path(__file__).resolve().parent
DEFAULT_CONDOR_DIR = PROCESSING_DIR / "condor"
QUEUE_HEADER = "queue job_name, input_list, output_file from ("
READ_ERROR = re.compile(
    r"(?:SKIPPING INPUT FILE|TLS error:|"
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


# Ordered CERN runtime classes; recovery must strictly increase the allowance.
FLAVOURS = ("espresso", "microcentury", "longlunch", "workday", "tomorrow", "testmatch", "nextweek")


def lifecycle(text: str) -> tuple[str, str | None]:
    """Use the latest lifecycle event, including removal/hold/eviction."""
    events = list(re.finditer(r"(?m)^(000|001|004|005|009|012|013) \((\d+\.\d+)\.\d+\)", text))
    if not events:
        return "unknown", None
    event = events[-1]
    body = text[event.start():]
    code = event.group(1)
    job_id = event.group(2)
    if code in ("000", "001", "004", "013"):
        return "active", job_id
    if code == "012":
        return "held", job_id
    if code == "009":
        return ("walltime" if re.search(r"wall.?time|allowed max|MaxRuntime", body, re.I) else "removed"), job_id
    if re.search(r"Normal termination.*return value 0\)", body, re.S):
        return "success", job_id
    return "execution_failure", job_id


def classify(row, event_text, error_text, active_ids, validator):
    state, job_id = lifecycle(event_text)
    # Condor can append a new attempt to the same log. Diagnostics from an
    # earlier submission must not turn its successful replacement into failure.
    submissions = list(re.finditer(r"(?m)^000 \(\d+\.\d+\.\d+\)", event_text))
    if submissions:
        event_text = event_text[submissions[-1].start():]
    transfer_error = bool(re.search(r"failed to transfer|transfer.*failed|transfer.*error", event_text + error_text, re.I))
    if state == "held" and transfer_error:
        return "held_transfer_failure", "scheduler still owns job; inspect HoldReason before release"
    if job_id in active_ids or state in ("active", "held"):
        return state if state == "held" else "active", "scheduler still owns job"
    ok, reason = validator(row.output_file)
    if reason.startswith("validation unavailable"):
        return "manual", reason
    if state == "walltime":
        return "walltime", reason
    if state == "unknown":
        return "manual", "no terminal scheduler evidence; " + reason
    if transfer_error:
        return "transfer_failure", reason
    if READ_ERROR.search(error_text):
        return "transient", reason
    if state == "success" and ok:
        return "complete", reason
    if state == "success":
        return "invalid_output", reason
    if state == "execution_failure" and "status 75" in error_text:
        return "invalid_output", reason
    return "manual", state + "; " + reason


def query_active_ids() -> set[str]:
    """Require a valid live queue response before authorizing recovery."""
    try:
        result = subprocess.run(["condor_q", "-json"], capture_output=True, text=True)
    except OSError as error:
        raise ValueError(f"condor_q could not run: {error}") from error
    if result.returncode:
        raise ValueError(f"condor_q exited {result.returncode}: {result.stderr.strip()}")
    try:
        ads = json.loads(result.stdout)
    except json.JSONDecodeError as error:
        raise ValueError("condor_q -json returned empty or invalid JSON; "
                         "inspect its stdout and stderr directly") from error
    if not isinstance(ads, list) or any(
        not isinstance(ad, dict)
        or type(ad.get("ClusterId")) is not int
        or type(ad.get("ProcId")) is not int
        for ad in ads
    ):
        raise ValueError("condor_q -json returned an unexpected job-record format")
    return {f"{ad['ClusterId']}.{ad['ProcId']}" for ad in ads}


def main() -> int:
    """Inspect every expected row; explicit submission is bounded and locked."""
    import collections
    import fcntl
    import hashlib
    from validate_job_output import validate

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("condor_dir", nargs="?", type=Path, default=DEFAULT_CONDOR_DIR)
    parser.add_argument("--submit", action="store_true")
    parser.add_argument("--max-retries", type=int, default=3)
    parser.add_argument("--walltime-flavour", choices=FLAVOURS,
                        help="Explicit longer runtime for wall-time failures; no event splitting")
    parser.add_argument("--report", type=Path, help="Write a JSON recovery plan")
    args = parser.parse_args()
    if args.max_retries < 0:
        parser.error("--max-retries must be nonnegative")
    root = args.condor_dir.expanduser().resolve()
    if not root.is_dir():
        parser.error("campaign directory does not exist")
    # A failed queue query must never authorize a submission. Dry runs remain
    # usable off-cluster, with event-log evidence and explicit diagnostics.
    active_ids = set()
    queue_available = False
    if shutil.which("condor_q"):
        try:
            active_ids = query_active_ids()
            queue_available = True
        except ValueError as error:
            print(f"Queue check failed: {error}")
    if args.submit and (not queue_available or not shutil.which("condor_submit")):
        parser.error("submission requires successful condor_q -json and condor_submit")
    if not queue_available:
        print("Queue unavailable: dry-run classifications use retained event logs only.")
    report = []
    # Serialize submitters across scans of the same campaign/sample tree.
    # Each original submit file has its own lock and durable attempt journal.
    submits = original_submit_files(root)
    parsed = {}
    outputs = collections.Counter()
    for sub in submits:
        if QUEUE_HEADER not in sub.read_text():
            report.append(dict(submit_file=str(sub), category="manual", eligible=False,
                               reason="unsupported queue schema"))
            continue
        parsed[sub] = parse_queue(sub)
        outputs.update(str(row.output_file.resolve()) for row in parsed[sub][3])
    for sub, (lines, start, end, rows) in parsed.items():
        if not any("run_process_forest_condor.sh" in line for line in lines):
            report.append(dict(submit_file=str(sub), category="manual", eligible=False,
                               reason="unsupported worker"))
            continue
        manifest_path = sub.with_suffix(".manifest.json")
        manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else None
        lock = sub.with_suffix(".recovery.lock")
        with lock.open("a") as stream:
            fcntl.flock(stream, fcntl.LOCK_EX)
            journal_path = sub.with_suffix(".recovery.json")
            journal = json.loads(journal_path.read_text()) if journal_path.exists() else {}
            def save():
                temporary = journal_path.with_suffix(".json.tmp")
                temporary.write_text(json.dumps(journal, indent=2) + "\n")
                temporary.replace(journal_path)
            for row in rows:
                logs = sub.parent / "log"
                event = logs / f"{row.job_name}.log"
                err = logs / f"{row.job_name}.err"
                event_text = event.read_text(errors="replace") if event.exists() else ""
                error_text = err.read_text(errors="replace") if err.exists() else ""
                state, reason = classify(row, event_text, error_text, active_ids, validate)
                if not row.input_list.is_file():
                    state, reason = "manual", "original input list is missing"
                elif manifest:
                    expected = manifest["jobs"].get(row.job_name, {})
                    if (manifest["submit_sha256"] != hashlib.sha256(sub.read_bytes()).hexdigest()
                            or expected.get("input_sha256") != hashlib.sha256(row.input_list.read_bytes()).hexdigest()):
                        state, reason = "manual", "workload changed since submission; reconcile manifest"
                if outputs[str(row.output_file.resolve())] > 1:
                    state, reason = "manual", "multiple original queue rows target this output; select one campaign"
                attempts = journal.get(row.job_name, [])
                eligible = state in ("transient", "invalid_output", "transfer_failure", "walltime")
                flavour = attempts[-1].get("flavour") if attempts else None
                if attempts and attempts[-1]["state"] == "uncertain":
                    eligible = False
                    state = "submission_uncertain"
                    reason += "; submission outcome uncertain: reconcile journal with scheduler"
                if attempts and attempts[-1]["state"] == "submitted":
                    last_id = attempts[-1]["job_id"]
                    log_state, log_id = lifecycle(event_text)
                    if last_id in active_ids or log_id != last_id or log_state in ("active", "unknown", "held"):
                        eligible = False
                        state = "recovery_pending"
                        reason += "; recovery pending scheduler/log confirmation"
                attempt_count = max(len(attempts), len(list((sub.parent / "retry_archive").glob(f"*/{row.job_name}"))))
                if attempt_count >= args.max_retries:
                    eligible = False
                    reason += "; retry limit reached"
                if state == "walltime":
                    match = re.search(r'(?im)^\+JobFlavour\s*=\s*"([^"\n]+)"', "\n".join(lines))
                    previous = attempts[-1].get("flavour") if attempts else None
                    previous = previous or (match.group(1) if match else None)
                    flavour = args.walltime_flavour
                    explicit_runtime = any(re.match(r"(?i)^\+?MaxRuntime\s*=", line.strip()) for line in lines)
                    if explicit_runtime or previous not in FLAVOURS or flavour is None or FLAVOURS.index(flavour) <= FLAVOURS.index(previous):
                        eligible = False
                        reason += "; requires explicitly selected strictly longer --walltime-flavour"
                item = dict(submit_file=str(sub), job=row.job_name, input=str(row.input_list),
                            output=str(row.output_file), category=state, reason=reason,
                            attempts=attempt_count, eligible=eligible)
                report.append(item)
                print(f"{row.job_name}: {state}; eligible={eligible}; {reason}")
                if not args.submit or not eligible:
                    continue
                # Check live scheduler again under the lock before mutation.
                try:
                    active_ids = query_active_ids()
                except ValueError as error:
                    parser.error(str(error))
                if lifecycle(event_text)[1] in active_ids:
                    raise SystemExit("Job is still in queue; rerun to inspect current state")
                tag = dt.datetime.now().strftime("%Y%m%d_%H%M%S%f")
                job = FailedJob(row, logs / f"{row.job_name}.out", err, event)
                retry = write_retry_submit(sub, lines, start, end, [job], tag)
                if flavour:
                    text = retry.read_text()
                    text = re.sub(r'(?im)^\+JobFlavour\s*=.*$', f'+JobFlavour = "{flavour}"', text)
                    retry.write_text(text)
                attempt = dict(state="uncertain", submit_file=str(retry), flavour=flavour,
                               original_input=str(row.input_list), output=str(row.output_file))
                attempts.append(attempt)
                journal[row.job_name] = attempts
                save()  # persist before any external side effect; crash fails closed
                clear_failed_artifacts(sub, [job], tag)
                result = subprocess.run(["condor_submit", "-terse", retry.name], cwd=sub.parent,
                                        capture_output=True, text=True)
                match = re.fullmatch(r"\s*(\d+\.\d+)(?:\s*-\s*\d+\.\d+)?\s*", result.stdout)
                if result.returncode or not match:
                    raise SystemExit("Submission outcome uncertain; inspect journal and scheduler: " + result.stderr)
                attempt.update(state="submitted", job_id=match.group(1))
                save()
                print(f"Submitted {row.job_name}: {match.group(1)}")
                # Refresh baseline for the next independent submission.
                try:
                    active_ids = query_active_ids()
                except ValueError as error:
                    parser.error(str(error))
    counts = collections.Counter(item["category"] for item in report)
    print("Summary:", dict(counts), "eligible:", sum(item["eligible"] for item in report))
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    return 0 if report and all(item["category"] == "complete" for item in report) else 1


if __name__ == "__main__":
    raise SystemExit(main())
