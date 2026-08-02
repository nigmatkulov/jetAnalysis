#!/usr/bin/env python3
"""Safely merge one data trigger/direction with bounded ``hadd`` fan-in.

Job outputs are discovered recursively from the naming convention produced by
``submit_process_forest_condor.py``. Large input sets are reduced through
validated intermediate ROOT files, and the final candidate is atomically moved
into place. Source job outputs are deleted only after the installed final file
passes the same validation.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path


# These labels deliberately match the submitter's case-sensitive output names.
TRIGGERS = ("MB", "Jet60", "Jet80", "Jet100")
# The submitter maps numeric reco-selection methods 2 (tight+lepton veto) and 3
# (loose) to the same jetId label. This merger therefore selects by the stored
# label and cannot separate methods 2 and 3 if their outputs share a directory.
RECO_JET_SELECTIONS = ("jetId", "trkMax", "noSel")


def positive_batch_size(value: str) -> int:
    """Parse a fan-in that always stays below ROOT's problematic 1,000 files."""

    parsed = int(value)
    if not 2 <= parsed <= 999:
        raise argparse.ArgumentTypeError("batch size must be between 2 and 999")
    return parsed


def parse_args() -> argparse.Namespace:
    """Define discovery, output, fan-in, and source-retention controls."""

    parser = argparse.ArgumentParser(
        description=(
            "Recursively hadd one trigger and beam orientation in batches, "
            "verify the final ROOT file, and remove the merged job outputs."
        )
    )
    parser.add_argument("input_dir", type=Path, help="Directory searched recursively")
    parser.add_argument("--trigger", choices=TRIGGERS, required=True)
    parser.add_argument("--beam", choices=("Pbgoing", "pgoing"), required=True)
    parser.add_argument(
        "--reco-jet-selection",
        choices=RECO_JET_SELECTIONS,
        default="jetId",
        help="Selection label in source and output names (default: jetId)",
    )
    parser.add_argument(
        "--output", type=Path,
        help=(
            "Final ROOT file (default: "
            "INPUT_DIR/<trigger>_<beam>_ak4_<reco-selection>.root)"
        ),
    )
    parser.add_argument("--batch-size", type=positive_batch_size, default=500)
    parser.add_argument(
        "--keep-inputs", action="store_true",
        help="Keep source job files after a successful, verified merge",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print the selected files and output without running hadd or deleting files",
    )
    return parser.parse_args()


def input_pattern(trigger: str, beam: str, reco_jet_selection: str) -> str:
    """Return the submitter pattern for one trigger/direction/selection tuple.

    MB includes a primary-dataset label between the trigger and direction;
    PAEGJet-derived trigger samples do not. The selection label is included so
    neighboring campaigns produced with another selection are never deleted.
    """

    if trigger == "MB":
        return f"MB_PD*_{beam}_data_{reco_jet_selection}_job*.root"
    return f"{trigger}_{beam}_data_{reco_jet_selection}_job*.root"


def validate_root_file(path: Path, rootls: str) -> None:
    """Require a nonempty ROOT file that exposes at least one top-level key.

    A successful ``hadd`` exit alone is not used as permission to delete input
    files. ``rootls`` independently reopens the product through ROOT and fails
    for unreadable files; requiring output also rejects structurally empty
    products.
    """

    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"missing or empty merged ROOT file: {path}")
    completed = subprocess.run(
        [rootls, "-1", str(path)], capture_output=True, text=True, check=False
    )
    if completed.returncode != 0 or not completed.stdout.strip():
        detail = completed.stderr.strip() or "ROOT file contains no top-level keys"
        raise RuntimeError(f"failed to validate {path}: {detail}")


def run_hadd(output: Path, inputs: list[Path], hadd: str, rootls: str) -> None:
    """Merge and immediately validate one bounded group of ROOT files."""

    print(f"Merging {len(inputs)} file(s) -> {output}")
    # -f208 both permits replacement of a temporary target and requests the
    # user-specified LZMA compression algorithm (2), level 8.
    completed = subprocess.run(
        [hadd, "-f208", str(output), *(str(path) for path in inputs)],
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"hadd failed with exit code {completed.returncode}: {output}")
    validate_root_file(output, rootls)


def merge_in_stages(
    sources: list[Path], output: Path, batch_size: int, hadd: str, rootls: str
) -> None:
    """Reduce sources in bounded stages and atomically install the final file.

    Intermediates are created in the output directory so the final ``replace``
    remains on one filesystem. TemporaryDirectory removes all partial products
    on either success or failure; it never contains the original job outputs.
    """

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".hadd-data-", dir=output.parent) as raw_temp:
        temp_dir = Path(raw_temp)
        current = sources
        stage = 0
        while len(current) > batch_size:
            next_stage: list[Path] = []
            # Every child invocation receives at most batch_size inputs. The
            # resulting partials become the complete input set for the next
            # stage, repeating until one final bounded merge is possible.
            for batch_number, offset in enumerate(range(0, len(current), batch_size), 1):
                partial = temp_dir / f"stage{stage:02d}_batch{batch_number:04d}.root"
                run_hadd(partial, current[offset : offset + batch_size], hadd, rootls)
                next_stage.append(partial)
            current = next_stage
            stage += 1
        final_candidate = temp_dir / "final.root"
        run_hadd(final_candidate, current, hadd, rootls)
        # Do not let a failed hadd damage an older final product. Only a fully
        # validated candidate replaces it, atomically on the same filesystem.
        final_candidate.replace(output)


def main() -> int:
    """Discover, merge, validate, and optionally remove one source set."""

    args = parse_args()
    input_dir = args.input_dir.expanduser().resolve()
    if not input_dir.is_dir():
        raise SystemExit(f"input directory does not exist: {input_dir}")

    default_output = input_dir / (
        f"{args.trigger.lower()}_{args.beam}_ak4_{args.reco_jet_selection}.root"
    )
    output = args.output.expanduser().resolve() if args.output else default_output
    pattern = input_pattern(args.trigger, args.beam, args.reco_jet_selection)
    # Resolve and deduplicate paths before merging. Excluding output is a
    # second guard against ever feeding the target back into hadd if a custom
    # output name happens to match the discovery pattern.
    sources = sorted({
        path.resolve()
        for path in input_dir.rglob(pattern)
        if path.resolve() != output
    })
    if not sources:
        raise SystemExit(f"no files matched {pattern!r} beneath {input_dir}")

    print(f"Selected {len(sources)} source file(s)")
    print(f"Final output: {output}")
    if args.dry_run:
        for source in sources:
            print(source)
        print("Dry run: no files were merged or removed.")
        return 0

    hadd = shutil.which("hadd")
    rootls = shutil.which("rootls")
    if hadd is None or rootls is None:
        raise SystemExit("both hadd and rootls must be available in PATH")

    try:
        merge_in_stages(sources, output, args.batch_size, hadd, rootls)
        # Reopen the installed path rather than relying only on validation of
        # its temporary predecessor.
        validate_root_file(output, rootls)
    except RuntimeError as error:
        raise SystemExit(str(error)) from error

    if args.keep_inputs:
        print(f"Verified {output}; keeping {len(sources)} source file(s).")
        return 0

    # This is intentionally the last state-changing phase: neither partial
    # merges nor a failed final validation can remove original job outputs.
    for source in sources:
        source.unlink()
    print(f"Verified {output}; removed {len(sources)} source file(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
