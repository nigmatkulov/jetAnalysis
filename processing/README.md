# HTCondor processing

`submit_process_forest_condor.py` splits a text file of ROOT input paths and
submits one `processForestSimple` job per sublist. Its analysis options match the
current executable interface; the obsolete JEU, JER, and upper-pT-hat arguments
from the historical `submit_dijetAna_pPb8160.sh` workflow are intentionally not
included.

The commands in this document are intended to be run from the `processing/`
directory, which is also the conventional submission directory:

```bash
cd /afs/cern.ch/user/g/gnigmatk/soft/jetAnalysis/processing
```

Build the executable from this directory before submission:

```bash
cmake -S .. -B ../build -DCMAKE_BUILD_TYPE=Release
cmake --build ../build -j
```

Create or copy a valid proxy to the default protected location:

```bash
voms-proxy-init --voms cms --valid 24:00
cp "$(voms-proxy-info -path)" voms_proxy.txt
voms-proxy-info -timeleft
```

The proxy is ignored by Git. Generate a data campaign without submitting it:

```bash
../py-env/bin/python submit_process_forest_condor.py \
  filelists/pPb8160/DATA_MB/Pbgoing/MB_PD12_Pbgoing.txt \
  /eos/user/g/gnigmatk/ana/pPb8160/exp/Pbgoing \
  --mc-type 0 --beam Pbgoing --trigger-id 0 --pd-number 12 \
  --reco-jet-selection 2 --files-per-job 10 --dry-run
```

Remove `--dry-run` to call `condor_submit`.

The general `submit_process_forest_condor.py` submitter uses 10 input files per
job by default. This limits the amount of otherwise valid work invalidated by
one transient input read failure and makes failed jobs cheaper to retry.
Override `--files-per-job` explicitly when a different split has been validated
for the campaign. The data-only bulk wrapper documented below deliberately
overrides this with 50 files per job.

The submitter requests CERN's `espresso` flavour (20 minutes) by default for
experimental data and `microcentury` (1 hour) for embedding and Pythia jobs.
Use `--job-flavour FLAVOUR` to override the sample-dependent default.

If `../py-env/bin/python` is not installed in the CERN checkout, use the
environment's `python3` command. Do not prefix the script or input paths with
`processing/` when the current directory is already `processing/`.

### Three-file Jet60 submission test

From the `processing/` directory, create a temporary list from the first three
Pb-going PAEGJet entries:

```bash
mkdir -p input
head -n 3 filelists/pPb8160/DATA_PAEGJet/Pbgoing/PAEGJet_Pbgoing.txt \
  > input/Jet60_Pbgoing_3files.list
```

Generate and inspect one Condor job without submitting it:

```bash
../py-env/bin/python submit_process_forest_condor.py \
  input/Jet60_Pbgoing_3files.list \
  /eos/user/g/gnigmatk/ana/pPb8160/test/Jet60/Pbgoing \
  --mc-type 0 --beam Pbgoing --trigger-id 1 \
  --reco-jet-selection 2 --files-per-job 3 \
  --name test3files --dry-run
```

After inspecting the generated `.sub` file under `condor/`, repeat the command
without `--dry-run` to submit it.

`--files-per-job 3` creates one Condor job containing all three test files. Use
`--files-per-job 1` to create three independent jobs, one for each input file.

Before the Condor test, the same list can be processed directly from the
repository root to verify list parsing and event access independently of the
scheduler:

```bash
cd ..
./build/processForestSimple \
  processing/input/Jet60_Pbgoing_3files.list \
  /tmp/Jet60_Pbgoing_3files.root 0 1 0 1 2 0
cd processing
```

Confirm that the reported accepted-plus-skipped count matches the sublist
length and that the event count is nonzero. Current HiForest list entries need
to contain `Forest` and `.root`; `AOD` is not required in the filename. Missing,
unreadable, or incomplete ROOT files print an explicit `SKIPPING INPUT FILE`
message on the error stream and processing continues with all remaining files.
The job fails if no usable files remain, a list entry is malformed, or the
aggregate tree entry counts are inconsistent.

Submit from a shell where the CMS/ROOT runtime environment is already loaded;
the generated description uses `getenv = True`. By default jobs require the
CERN AlmaLinux 9 QA environment, matching the previous submission workflow.
Override that expression with `--requirements` when submitting elsewhere.

Options passed to `processForestSimple` are:

- `--mc-type`: `0` data, `1` embedding, or `2` Pythia.
- `--beam`: `pgoing` (`isPbGoingDir=0`) or `Pbgoing` (`isPbGoingDir=1`).
- `--pt-hat`: one supported MC pT-hat sample or `all`; ignored for data.
- `--trigger-id`: `0` MB/no trigger, `1` Jet60, `2` Jet80, or `3` Jet100.
- `--pd-number`: required for MB data; ignored for Jet60/80/100 and MC.
- `--reco-jet-selection`: `0` none, `1` `trkMaxPt/RawPt`, `2` tight jet ID
  with lepton veto, or `3` loose jet ID.
- `--vertex-filter-selection`: `0` `pVertexFilterCutdz1p0` (nominal), `1`
  `pVertexFilterCutGplus`, or `2` `pVertexFilterCutVtx1`.

The default proxy is `voms_proxy.txt` in the `processing/` directory; override it with
`--proxy /path/to/x509_proxy`. The proxy is copied into every sample directory
under that fixed basename and transferred to each job. The
generated sublists, submit files, and logs are stored below
`condor/<campaign-id>/<sample-name>/`; outputs go to the requested
output directory. Repository, executable, generated-list, and output paths must
be visible from the worker nodes (as they are on CERN shared filesystems).
Whitespace is not supported in repository, executable, work, input-list, or
output paths because those values are emitted as unquoted Condor tokens.

## Unique and readable names

Every submission receives a campaign ID consisting of a microsecond timestamp
and an eight-character random UUID suffix. This keeps campaign directories and
submission files unique. ROOT output and log basenames contain the trigger, MB
primary-dataset number when applicable, beam orientation, sample type, MC
pT-hat value when applicable, jet selection (`noSel`, `trkMax`, or `jetId`), and
file-list job number. `--name` adds an optional prefix.

For example, MB PD12 job 3 may produce:

```text
MB_PD12_Pbgoing_data_jetId_job0003.root
```

The generated input list and Condor log files have the same readable basename.
Rerunning exactly the same configuration intentionally replaces matching ROOT
outputs. Different beam directions, triggers, PDs, pT-hat samples, selection
families, and job numbers cannot overwrite one another because those values are
part of the filename. Reco selection methods `2` and `3` both use the requested
short label `jetId`; do not write tight- and loose-ID campaigns to the same
output directory unless replacement is intended.

`voms_proxy.txt` is used by default. Every generated submit file
contains:

```text
transfer_input_files  = voms_proxy.txt
environment = "X509_USER_PROXY=voms_proxy.txt"
```

## Standard submissions

The `{pd}` and `{pt_hat}` text below is replaced by the submitter. Quote paths
containing these placeholders so the shell passes them unchanged. The commands
use the standard production output directories beneath
`/eos/user/g/gnigmatk/ana/pPb8160/`:

- experimental data: `exp/Pbgoing` and `exp/pgoing`;
- embedding: `embedding/Pbgoing` and `embedding/pgoing`;
- Pythia: `pythia/Pbgoing` and `pythia/pgoing`.

### Data

Submit all four trigger configurations (MB, Jet60, Jet80, and Jet100) for both
beam orientations with one command:

```bash
../py-env/bin/python submit_all_data_condor.py
```

This uses the standard output directories, 50 files per job, tight jet ID with
the lepton veto, every direction-appropriate MB PD, and the shared PAEGJet list
for each jet trigger. Add `--dry-run` to prepare all eight configurations
without calling `condor_submit`. Use `--help` to see the common output, split,
selection, job-flavour, proxy, and work-directory overrides.

The equivalent individual commands are listed below.

```bash
# MB Pb-going: PD1 through PD20, submitted sequentially
../py-env/bin/python submit_process_forest_condor.py \
  'filelists/pPb8160/DATA_MB/Pbgoing/MB_PD{pd}_Pbgoing.txt' \
  /eos/user/g/gnigmatk/ana/pPb8160/exp/Pbgoing \
  --mc-type 0 --beam Pbgoing --trigger-id 0 \
  --pd-number all --files-per-job 50 --reco-jet-selection 2

# MB p-going: PD1 through PD8, submitted sequentially
../py-env/bin/python submit_process_forest_condor.py \
  'filelists/pPb8160/DATA_MB/pgoing/MB_PD{pd}_pgoing.txt' \
  /eos/user/g/gnigmatk/ana/pPb8160/exp/pgoing \
  --mc-type 0 --beam pgoing --trigger-id 0 \
  --pd-number all --files-per-job 50 --reco-jet-selection 2

# Jet60 Pb-going
../py-env/bin/python submit_process_forest_condor.py \
  filelists/pPb8160/DATA_PAEGJet/Pbgoing/PAEGJet_Pbgoing.txt \
  /eos/user/g/gnigmatk/ana/pPb8160/exp/Pbgoing \
  --mc-type 0 --beam Pbgoing --trigger-id 1 \
  --files-per-job 50 --reco-jet-selection 2

# Jet60 p-going
../py-env/bin/python submit_process_forest_condor.py \
  filelists/pPb8160/DATA_PAEGJet/pgoing/PAEGJet_pgoing.txt \
  /eos/user/g/gnigmatk/ana/pPb8160/exp/pgoing \
  --mc-type 0 --beam pgoing --trigger-id 1 \
  --files-per-job 50 --reco-jet-selection 2

# Jet80 Pb-going
../py-env/bin/python submit_process_forest_condor.py \
  filelists/pPb8160/DATA_PAEGJet/Pbgoing/PAEGJet_Pbgoing.txt \
  /eos/user/g/gnigmatk/ana/pPb8160/exp/Pbgoing \
  --mc-type 0 --beam Pbgoing --trigger-id 2 \
  --files-per-job 50 --reco-jet-selection 2

# Jet80 p-going
../py-env/bin/python submit_process_forest_condor.py \
  filelists/pPb8160/DATA_PAEGJet/pgoing/PAEGJet_pgoing.txt \
  /eos/user/g/gnigmatk/ana/pPb8160/exp/pgoing \
  --mc-type 0 --beam pgoing --trigger-id 2 \
  --files-per-job 50 --reco-jet-selection 2

# Jet100 Pb-going
../py-env/bin/python submit_process_forest_condor.py \
  filelists/pPb8160/DATA_PAEGJet/Pbgoing/PAEGJet_Pbgoing.txt \
  /eos/user/g/gnigmatk/ana/pPb8160/exp/Pbgoing \
  --mc-type 0 --beam Pbgoing --trigger-id 3 \
  --files-per-job 50 --reco-jet-selection 2

# Jet100 p-going
../py-env/bin/python submit_process_forest_condor.py \
  filelists/pPb8160/DATA_PAEGJet/pgoing/PAEGJet_pgoing.txt \
  /eos/user/g/gnigmatk/ana/pPb8160/exp/pgoing \
  --mc-type 0 --beam pgoing --trigger-id 3 \
  --files-per-job 50 --reco-jet-selection 2
```

`--pd-number` has no effect for Jet60, Jet80, or Jet100; each uses the single
direction-specific PAEGJet list.

### Check data-file availability

Before submitting data, check every ROOT path in all MB lists (PD1–20 for
Pb-going and PD1–8 for p-going) and both direction-specific PAEGJet lists:

```bash
../py-env/bin/python check_data_filelists.py
```

The report is grouped by direction and PD and prints every unavailable path as
`MISSING`. It exits with status `0` when all files are available, `1` when any
listed ROOT file is unavailable, and `2` when a required file list itself is
missing. Availability checks use 16 concurrent workers by default; change this
with `--workers`. An overall progress bar shows the number and percentage of
completed checks. Add `--show-available` to print successful paths as well.

### Monte Carlo

Each command below creates and sequentially submits one independent Condor
description for every supported pT-hat sample: 15, 30, 50, 80, 120, 170, 220,
280, 370, 460, 540.

```bash
# Embedded Pb-going
../py-env/bin/python submit_process_forest_condor.py \
  'filelists/pPb8160/MC_embedded/Pbgoing/MC_pthat{pt_hat}_Pbgoing_embedded.txt' \
  /eos/user/g/gnigmatk/ana/pPb8160/embedding/Pbgoing \
  --mc-type 1 --beam Pbgoing --pt-hat all \
  --files-per-job 10 --reco-jet-selection 2

# Embedded p-going
../py-env/bin/python submit_process_forest_condor.py \
  'filelists/pPb8160/MC_embedded/pgoing/MC_pthat{pt_hat}_pgoing_embedded.txt' \
  /eos/user/g/gnigmatk/ana/pPb8160/embedding/pgoing \
  --mc-type 1 --beam pgoing --pt-hat all \
  --files-per-job 10 --reco-jet-selection 2

# Pythia Pb-going
../py-env/bin/python submit_process_forest_condor.py \
  'filelists/pPb8160/MC_Unembedded/Pbgoing/MC_pthat{pt_hat}_Pbgoing_unembedded.txt' \
  /eos/user/g/gnigmatk/ana/pPb8160/pythia/Pbgoing \
  --mc-type 2 --beam Pbgoing --pt-hat all \
  --files-per-job 10 --reco-jet-selection 2

# Pythia p-going
../py-env/bin/python submit_process_forest_condor.py \
  'filelists/pPb8160/MC_Unembedded/pgoing/MC_pthat{pt_hat}_pgoing_unembedded.txt' \
  /eos/user/g/gnigmatk/ana/pPb8160/pythia/pgoing \
  --mc-type 2 --beam pgoing --pt-hat all \
  --files-per-job 10 --reco-jet-selection 2
```

## Monitoring and validation

Monitor all of your jobs or one returned cluster ID:

```bash
condor_q "$USER"
condor_q CLUSTER_ID
```

If jobs stay idle unexpectedly, inspect matching constraints with:

```bash
condor_q -better-analyze CLUSTER_ID
```

After completion, list nonempty error files and inspect the corresponding
standard output for the file and event counts:

```bash
find condor -type f -name '*.err' -size +0 -print
find condor -type f -name '*.out' | sort | tail
```

A successful worker output reports the input sublist, a nonzero
`Total number of files in chain`, a nonzero entry count for nonempty inputs, and
the final processed/selected-event counts. ROOT outputs are written directly to
the requested shared output directory.

### Merge completed data outputs

First preview the recursively discovered MB Pb-going job outputs and standard
final filename:

```bash
../py-env/bin/python merge_data_outputs.py \
  /eos/user/g/gnigmatk/ana/pPb8160/exp/Pbgoing \
  --trigger MB --beam Pbgoing --vertex-filter-selection dz1p0 --dry-run
```

Review the selected paths, then remove `--dry-run` to merge with at most 500
inputs in each `hadd -f208` call. The successful merge verifies and installs
`mb_Pbgoing_ak4_jetId.root`, then removes only the matched per-job inputs.

The reconstruction-selection option controls both source discovery and the
default final filename:

| Merger option | Submitter method | Matched filename label | Final suffix |
| --- | --- | --- | --- |
| `jetId` (default) | `2` or `3` | `_data_jetId_job*.root` | `_ak4_jetId.root` |
| `trkMax` | `1` | `_data_trkMax_job*.root` | `_ak4_trkMax.root` |
| `noSel` | `0` | `_data_noSel_job*.root` | `_ak4_noSel.root` |

Reco-selection methods `2` (tight ID with lepton veto) and `3` (loose ID) both
use the `jetId` filename label. The merger cannot distinguish them, so do not
place both campaigns beneath the same merge directory. This is the same naming
constraint described in the submission section above.

Use `--vertex-filter-selection Gplus` or `Vtx1` to merge a pileup systematic
campaign. These labels are added to both the matched job filenames and the
default merged filename, for example
`mb_Pbgoing_ak4_jetId_Gplus.root`. The nominal `dz1p0` selection retains the
existing filenames.

The corresponding p-going command uses the `pgoing` directory and
`--beam pgoing`, producing `mb_pgoing_ak4_jetId.root`. Jet-trigger samples use
the same interface with `--trigger Jet60`, `Jet80`, or `Jet100`; their default
outputs are `jet60_<beam>_ak4_jetId.root`, `jet80_<beam>_ak4_jetId.root`, and
`jet100_<beam>_ak4_jetId.root`.

Pass `--output /path/to/name.root` to replace the default final path. This does
not change source discovery: `--trigger`, `--beam`, and
`--reco-jet-selection` and `--vertex-filter-selection` still determine exactly which job outputs are merged
and, after verification, removed.

The script requires every intermediate file and the final file to be nonempty
and reopenable by `rootls` with at least one top-level key. It builds the final
candidate beside the destination and replaces the destination atomically only
after validation. Only after reopening that installed final file successfully
does it remove the matched per-job ROOT files. A failed merge leaves all
original job outputs and any older final output intact. Use
`--dry-run` to inspect the recursive file selection without merging or deleting,
or `--keep-inputs` to merge and verify while retaining the source files. Change
the fan-in with `--batch-size`; values are restricted to 2–999.

### Retry completed jobs with ROOT input read errors

Wait for the campaigns being checked to finish, or expect active and
unterminated jobs to be reported and skipped. A nonempty `.err` file alone does
not mark a job for retry. From `processing/`, scan all retained campaigns
without changing files or submitting jobs:

```bash
../py-env/bin/python retry_failed_io_jobs.py
```

The scanner selects only jobs that satisfy both conditions:

- the Condor event log records that the job terminated;
- stderr contains a ROOT `TFile::ReadBuffer` `Input/output error`, a
  `TBranch::GetBasket` `badread` error, or the analysis `Failed to read chain
  entry` error.

Normal jobs, active jobs, and jobs with unrelated errors are not touched. The
final `Dry run: N failed job(s) would be retried.` line is the total number of
selected queue rows. Review the listed job names and output paths before
enabling submission. If `condor/` contains unrelated or historical campaigns,
pass the intended campaign subtree explicitly, for example:

```bash
../py-env/bin/python retry_failed_io_jobs.py \
  condor/20260803_211646_314197_f0ca9cb3
```

The retry submit descriptions reuse the proxy copied into each sample campaign
at its original submission. Check one campaign copy before resubmitting:

```bash
voms-proxy-info \
  -file condor/20260803_211646_314197_f0ca9cb3/MB_PD1_Pbgoing_data_jetId/voms_proxy.txt \
  -timeleft
```

Renewing only `processing/voms_proxy.txt` after the original submission does
not replace those campaign copies. If a campaign proxy is expired, prepare a
fresh campaign instead of submitting its generated retry description unchanged.

After reviewing the dry run and confirming the campaign proxy is valid, archive
the selected jobs' old logs, quarantine any stale ROOT outputs with an
`.invalid-<timestamp>` suffix, generate retry submit files containing only the
affected queue rows, and call `condor_submit` with:

```bash
../py-env/bin/python retry_failed_io_jobs.py --submit
```

The default `condor` argument is appropriate only when every retained campaign
is in scope. For a controlled recovery, add `--submit` to the same explicit
campaign path used for its dry run:

```bash
../py-env/bin/python retry_failed_io_jobs.py \
  condor/20260803_211646_314197_f0ca9cb3 --submit
```

Archived logs are stored below each sample campaign in
`retry_archive/<timestamp>/<job-name>/`. Any pre-existing output for a selected
job is preserved beside the original output path as
`<name>.root.invalid-<timestamp>` before resubmission. Successful jobs and their
outputs are not moved. The generated `_retry_<timestamp>.sub` files are not
treated as original campaigns during later scans, preventing duplicate queue
definitions from being scanned.

Record the cluster IDs printed by `condor_submit`, monitor them with
`condor_q`, and rerun the same dry-run command after they terminate. A clean
campaign reports `No completed jobs with matching ROOT input read errors were
found.` If some jobs repeatedly fail on the same EOS input, test that file
directly and retry it in a smaller job rather than repeatedly processing the
same multi-file queue row. Do not merge a trigger/direction until this check is
clean. Run `retry_failed_io_jobs.py --help` for its campaign-path and `--submit`
interface.
