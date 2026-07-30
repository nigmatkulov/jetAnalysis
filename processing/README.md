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
  --reco-jet-selection 2 --files-per-job 50 --dry-run
```

Remove `--dry-run` to call `condor_submit`.

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
  /tmp/Jet60_Pbgoing_3files.root 0 1 0 1 2
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
  --files-per-job 50 --reco-jet-selection 2

# Embedded p-going
../py-env/bin/python submit_process_forest_condor.py \
  'filelists/pPb8160/MC_embedded/pgoing/MC_pthat{pt_hat}_pgoing_embedded.txt' \
  /eos/user/g/gnigmatk/ana/pPb8160/embedding/pgoing \
  --mc-type 1 --beam pgoing --pt-hat all \
  --files-per-job 50 --reco-jet-selection 2

# Pythia Pb-going
../py-env/bin/python submit_process_forest_condor.py \
  'filelists/pPb8160/MC_Unembedded/Pbgoing/MC_pthat{pt_hat}_Pbgoing_unembedded.txt' \
  /eos/user/g/gnigmatk/ana/pPb8160/pythia/Pbgoing \
  --mc-type 2 --beam Pbgoing --pt-hat all \
  --files-per-job 50 --reco-jet-selection 2

# Pythia p-going
../py-env/bin/python submit_process_forest_condor.py \
  'filelists/pPb8160/MC_Unembedded/pgoing/MC_pthat{pt_hat}_pgoing_unembedded.txt' \
  /eos/user/g/gnigmatk/ana/pPb8160/pythia/pgoing \
  --mc-type 2 --beam pgoing --pt-hat all \
  --files-per-job 50 --reco-jet-selection 2
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

Run `--help` for naming, memory, job-flavour, executable, and work-directory
options.
