# HTCondor processing

`submit_process_forest_condor.py` splits a text file of ROOT input paths and
submits one `processForestSimple` job per sublist. Its analysis options match the
current executable interface; the obsolete JEU, JER, and upper-pT-hat arguments
from the historical `submit_dijetAna_pPb8160.sh` workflow are intentionally not
included.

For pPb MC, the producer normalizes by the entry count of each invocation.
Splitting one pT-hat sample across sublists and summing outputs therefore needs
an explicit common-denominator or merge-rescaling convention; it is not
partition invariant. Data has unit weights and pp MC uses stored weights, so
this specific limitation concerns pPb MC. Neither the submitter nor the data merger
applies an MC normalization correction. Keep the split fixed for comparisons;
absolute yields and combinations of differently partitioned samples require a
validated sample-wide normalization. The data merger below is not an MC merger.

The commands in this document are intended to be run from the `processing/`
directory, which is also the conventional submission directory:

```bash
cd /path/to/jetAnalysis/processing
```

Load the same CMS/ROOT runtime used for the executable in the submission shell.
For the current EL8 CMSSW setup:

```bash
source ~/setup_cmsenv.sh
echo "$CMSSW_BASE"
echo "$SCRAM_ARCH"       # expected: el8_amd64_gcc12
root-config --version   # current setup: 6.26/11
```

Source the setup script rather than running it as a child process so its
exports remain in the submission shell. Build from this directory:


```bash
cmake -S .. -B ../build -DCMAKE_BUILD_TYPE=Release
cmake --build ../build -j
```

Load the ROOT/CMS runtime on the submission host before building, and use a
Linux build compatible with the selected worker environment. The repository,
executable, calibration files under `aux_files/`, and generated sublists remain
on shared storage; they are not bundled into the job. Do not rebuild or edit
that shared analysis installation while a campaign is running.

## Native XRootD inputs

`submit_process_forest_condor.py` (also invoked by `submit_all_data_condor.py`)
converts EOS paths in generated job sublists automatically:

- `/eos/cms/...` becomes `root://eoscms.cern.ch//eos/cms/...`.
- `/eos/user/...` becomes `root://eosuser.cern.ch//eos/user/...`.
- `/eos/home-g/USER/...` becomes `root://eosuser.cern.ch//eos/user/g/USER/...`
  (and similarly for other home letters).

Existing URLs and ordinary local paths remain unchanged. Unknown EOS mount
names require an explicit URL. Original input lists are not modified; ordering,
duplicates, job grouping, and output names are preserved. ROOT streams input
through XRootD without copying forests to scratch. Output paths and the shared
executable/runtime setup remain unchanged. No disk request, runtime bundle,
or new automatic whole-job retries are added. Existing I/O failure reporting
and `retry_failed_io_jobs.py` remain available.

Use the usual submission command with `--dry-run` first and inspect the
resulting `input/*.list` files. Remove `--dry-run` to submit. No staging option
is required. A valid transferred proxy and the matching ROOT/XRootD runtime
are still required; verify one CERN pilot job before a full campaign.

Create or copy a valid proxy to the default ignored location:

```bash
voms-proxy-init --voms cms --valid 24:00
cp "$(voms-proxy-info -path)" voms_proxy.txt
chmod 600 voms_proxy.txt
voms-proxy-info -file voms_proxy.txt -timeleft
```

A dry run still requires an existing executable (analysis submitter), macro
(cloning submitter), input list, and proxy file. It creates campaign directories,
sublists, submit files and proxy copies; the analysis submitter also creates the
output directory. It does not validate proxy lifetime, open the ROOT inputs or
run the analysis. Check proxy lifetime against queue wait plus runtime. Inspect
all generated descriptions before repeating a command without `--dry-run`.
The personal EOS paths below are examples and script defaults; replace them
with shared writable paths for your account (`--output-base` for the bulk data
wrapper, positional output directory for the general submitters).

### Clone the CP5 Pythia forests

`submit_clone_forest_condor.py` submits one `macro/cloneForest.C` invocation
per ROOT file. With no positional arguments it uses
`filelists/pPb8160/PYTHIA_CP5_Tune/PYTHIA8160.txt` and writes the cloned files
to `/eos/user/g/gnigmatk/ana/pPb8160/pythia/cp5tune`. Entries beginning with
`/store/` are converted to global CMS XRootD URLs on submission.

After the cloned files are available locally under
`$HOME/cernbox/ana/pPb8160/pythia/cp5tune/forest/`, process them from the
repository root with the dedicated pp runner:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
py-env/bin/python pp_analyzeMc_cp5tune.py --workers 4 --dry-run
py-env/bin/python pp_analyzeMc_cp5tune.py --workers 4
```

The runner writes one ROOT output and one log per forest to the parent
`cp5tune/` directory. The main `README.md` documents the pp physics behavior,
file naming, parallel execution, and overwrite policy.

Prepare and inspect the campaign without submitting it (one job per usable
input-file entry; the printed count reflects the current list):

```bash
../py-env/bin/python submit_clone_forest_condor.py --dry-run
```

Submit after checking the generated description beneath `condor/`:

```bash
../py-env/bin/python submit_clone_forest_condor.py
```

An alternative input list and output directory can be passed positionally.
Use `--help` for macro, campaign-name, proxy, resource, and work-directory
overrides. Input basenames must be unique because the macro preserves them as
output filenames. The worker also fails the Condor job if ROOT exits without
creating a nonempty output file.

The same submitter has presets for every pT-hat list in the unembedded and
embedded MC samples. Their default output directories are respectively
`pythia/<direction>/forest` and `embedding/<direction>/forest`:

```bash
# Prepare all 11 pT-hat samples for one direction
../py-env/bin/python submit_clone_forest_condor.py \
  --sample unembedded --beam Pbgoing --dry-run
../py-env/bin/python submit_clone_forest_condor.py \
  --sample embedded --beam pgoing --dry-run

# Prepare all 22 pT-hat lists in one campaign (one job per input file)
../py-env/bin/python submit_clone_forest_condor.py \
  --sample unembedded --beam all --dry-run

# Prepare only one pT-hat sample
../py-env/bin/python submit_clone_forest_condor.py \
  --sample embedded --beam Pbgoing --pt-hat 120 --dry-run
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

Before the Condor test, a generated job sublist can be processed directly from
the repository root to verify list parsing and event access independently of the
scheduler. Replace the placeholder below with the generated sublist path; it
contains the XRootD URLs, whereas the original source list is unchanged:

```bash
cd ..
./build/processForestSimple \
  processing/condor/CAMPAIGN/SAMPLE/input/JOB.list \
  /tmp/Jet60_Pbgoing_3files.root 0 1 0 1 2 0
cd processing
```

Confirm that the reported accepted-plus-skipped count matches the sublist
length and that the event count is nonzero. Current HiForest list entries need
to contain `Forest` and `.root`; `AOD` is not required in the filename. Missing,
unreadable, or incomplete ROOT files print an explicit `SKIPPING INPUT FILE`
message on the error stream and processing continues with all remaining files.
The job fails if no usable files remain, a list entry is malformed, or the
per-file or aggregate tree entry counts are inconsistent. The Condor wrapper
also converts recognized ROOT read errors, TLS errors, and `SKIPPING INPUT FILE`
messages into exit status 74, even if the executable returns success. Require
`skipped=0` when checking a pilot and before accepting production outputs. The
worker resolves the transferred proxy to an absolute path before changing to
the repository directory; an explicitly configured unreadable proxy fails early.

Submit from a shell where the CMS/ROOT runtime environment is already loaded;
the generated description uses `getenv = True`. By default jobs require the
CERN AlmaLinux 9 QA environment, matching the previous submission workflow.
This default is retained from the existing workflow; it does not select an
EL8 container. For the current `el8_amd64_gcc12` build, select a compatible EL8
execution environment before submission. `--requirements` overrides worker
matching in the general submitter; the bulk wrapper does not expose that
option. Runtime compatibility is independent of the input URL conversion.

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
- `--residual-mode`: data only, `2` new residual (default), `1` old residual,
  `0` no residual. Ignored for MC; the beam-dependent L2 correction is always
  retained. Mode `1` uses `Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PF.txt`;
  mode `2` uses `Summer16_07Aug2017GH_V11_DATA_L2L3Residual_AK4PF.txt`.
- `--jer-smearing-mode`: `0` CMS hybrid (default), `1` forced CMS stochastic,
  or `2` truth-to-pseudo-reco diagnostic. It affects MC only.
- `--eta-dependent-jer`: optional override making standard JER branches use
  eta-dependent resolution too. Omit it for the default production of both
  distinct resolution families in one pass. It affects MC only.

This submitter supports pPb processing only: it has no collision-system option.
For the pp CP5 sample use `pp_analyzeMc_cp5tune.py`, as described above.
Embedding selects direction-specific Embedded JEC; standalone pPb PYTHIA selects
direction-specific Unembedded JEC. The pp runner selects Pb-going pPb Unembedded
JEC by the agreed analysis convention. No separate JEC payload option is exposed.

Both `submit_process_forest_condor.py` and `submit_all_data_condor.py` accept
`--residual-mode`. New-residual names stay unchanged. Old/no-residual data
campaigns add `_oldResiduals`/`_noResiduals` before `_jobNNNN.root`, e.g.
`Jet60_Pbgoing_data_jetId_oldResiduals_job0001.root`. The same suffix separates
campaign directories and logs. MC names remain unchanged.

Use the same mode with `merge_data_outputs.py`: its discovery pattern excludes
the other residual variants, and its default merged filename adds the suffix,
e.g. `jet60_Pbgoing_ak4_jetId_oldResiduals.root`. Mode `2` retains existing
merged names. For example, from the repository root:

```bash
py-env/bin/python processing/submit_all_data_condor.py --residual-mode 1 --dry-run
py-env/bin/python processing/merge_data_outputs.py OUTPUT_DIR --trigger Jet60 --beam Pbgoing --residual-mode 1 --dry-run
```

The generated worker command passes the residual mode as the final execution
argument. Existing 13-argument worker calls remain supported with mode `2`.

The generated worker command always passes `pPb` and both JER settings
explicitly. The default MC production uses CMS Hybrid mode with
`useEtaDependentJER=0`: one pass writes both base-resolution
`JerDef/Up/Down` histograms and eta-dependent
`JerDefExtra/UpExtra/DownExtra` histograms, together with both unfolding
response/fake/miss families. No second run is needed.

The MC commands below deliberately omit `--eta-dependent-jer`. Adding that
flag makes corresponding standard and `Extra` branches use the same
eta-dependent resolution and cached momenta; it does not enable an otherwise
missing `Extra` family. In Hybrid mode the resolution changes matching
eligibility and stochastic fallback; matched jets retain residual scaling.
Use distinct output directories if you intentionally produce the override too.

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
families, and job numbers have distinct labels. However, campaign IDs, input
list contents, split size and JER settings are NOT encoded in output basenames.
A fresh campaign directory does not protect an existing ROOT output. Changing
`--files-per-job` or the input list changes which files contribute to the same
job number; use a separate output directory or a distinct `--name` for reruns
that must coexist. Reco selection methods `2` and `3` both use the requested
short label `jetId`; do not write tight- and loose-ID campaigns to the same
output directory unless replacement is intended.

### AFS-safe logging

Generated Condor descriptions send successful job stdout to `/dev/null`, retain
one small scheduler lifecycle `.log` per job for status and retry checks, and
retain a per-job `.err` only when the worker reports a failure. Workers capture
the application output in node-local temporary storage and emit at most the
last 64 KiB to `.err` on failure. Email notification and live output streaming
are disabled. This bounds persistent text output and avoids exhausting the
2 GB AFS quota during large campaigns.

For `processForestSimple`, the worker also recognizes the ROOT read-error
signatures used by `retry_failed_io_jobs.py`. If one appears even when the
executable otherwise returns success, the worker fails the Condor job and
retains the bounded diagnostic. Normal ROOT progress and warning output is not
copied to AFS. Physics ROOT outputs continue to be written to the requested EOS
or shared-filesystem directory and are not subject to this text-log cap.

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
  --files-per-job 10 --reco-jet-selection 2 \
  --jer-smearing-mode 0

# Embedded p-going
../py-env/bin/python submit_process_forest_condor.py \
  'filelists/pPb8160/MC_embedded/pgoing/MC_pthat{pt_hat}_pgoing_embedded.txt' \
  /eos/user/g/gnigmatk/ana/pPb8160/embedding/pgoing \
  --mc-type 1 --beam pgoing --pt-hat all \
  --files-per-job 10 --reco-jet-selection 2 \
  --jer-smearing-mode 0

# Pythia Pb-going
../py-env/bin/python submit_process_forest_condor.py \
  'filelists/pPb8160/MC_Unembedded/Pbgoing/MC_pthat{pt_hat}_Pbgoing_unembedded.txt' \
  /eos/user/g/gnigmatk/ana/pPb8160/pythia/Pbgoing \
  --mc-type 2 --beam Pbgoing --pt-hat all \
  --files-per-job 10 --reco-jet-selection 2 \
  --jer-smearing-mode 0

# Pythia p-going
../py-env/bin/python submit_process_forest_condor.py \
  'filelists/pPb8160/MC_Unembedded/pgoing/MC_pthat{pt_hat}_pgoing_unembedded.txt' \
  /eos/user/g/gnigmatk/ana/pPb8160/pythia/pgoing \
  --mc-type 2 --beam pgoing --pt-hat all \
  --files-per-job 10 --reco-jet-selection 2 \
  --jer-smearing-mode 0
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

After completion, list nonempty error files and inspect scheduler termination
states. Successful jobs intentionally have no `.out` file; their `.err` may be empty:

```bash
find condor -type f -name '*.err' -size +0 -print
find condor -type f -name '*.log' | sort | tail
```

Use the recovery audit below before merging. A filename or an empty stderr file
alone does not establish completeness. Outputs are written directly to shared
storage, rather than transferred back by Condor; output write/access failures
are therefore distinct from Condor sandbox/proxy transfer failures.

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

The merger does not consult Condor logs, verify campaign completeness or deduplicate
events. Its recursive search can combine matching outputs from multiple runs.
A `--name` prefix changes the source basename and is not matched by this merger's
unprefixed patterns. Merge a directory containing exactly one intended campaign
configuration; use a separate merge procedure for prefixed outputs.

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

### Audit and recover incomplete processing jobs

The existing retry entry point now audits **every queue row**, including after
jobs leave the queue. Run from `processing/` with the intended campaign subtree:

```bash
../py-env/bin/python retry_failed_io_jobs.py condor/CAMPAIGN --report /tmp/recovery.json
```

This prints each expected unit and a category/eligibility summary. Exit status
is 0 only when every supported expected unit is complete, and 1 when any remains
incomplete (including after submission). Invalid arguments/unsafe submission
prerequisites return 2. The JSON report includes original input/output paths,
category, reason, attempts, and retry eligibility. Dry runs do not submit,
archive, or move outputs; a small advisory lock file may be created.

Classification combines current `condor_q -json` records, the latest retained
scheduler lifecycle event, executable status/read errors and output validation:

- `complete`: normal scheduler termination with exit 0, no known read errors,
  and a valid ROOT output.
- `active`, `held`, `held_transfer_failure`: scheduler-owned jobs; never duplicated
  or automatically released/cancelled.
- `transient`: terminal job with known EOS/ROOT read-error diagnostics.
- `walltime`: scheduler removal explicitly identifying the wall-time limit.
- `invalid_output`: missing, empty or invalid output after successful termination,
  or worker output-validation failure (75).
- `transfer_failure`: terminal transfer errors; held transfer failures require intervention.
- `recovery_pending`, `submission_uncertain`: pending or ambiguous recovery attempts;
  never duplicated until scheduler records are reconciled.
- `manual`: missing lifecycle evidence, other executable/removal failures,
  unavailable validator, missing/changed input lists, or other uncertainty.

ROOT validation uses PyROOT, or the `root` CLI if PyROOT is unavailable. It checks
nonempty storage, successful open, zombie/recovered flags, the existing `hVzRaw`
key, and readability of all top-level objects. This can be expensive for large
outputs. It does not prove event-level completeness alone: normal execution and
absence of known skipped-input/read errors are also required. The worker rejects
existing outputs (73), retains bounded failure diagnostics, returns 74 for known
read errors, returns 75 for invalid outputs, and records `JETANALYSIS_COMPLETED`
after exit 0 and validation. Set `JETANALYSIS_PYTHON` to a worker Python with
PyROOT, or ensure `python3` and the ROOT CLI are available in its runtime.

New submission descriptions have a `.manifest.json` recording submit/input-list
SHA-256 digests and expected output paths. Recovery refuses changed workloads.
Older campaigns use their retained queue rows/lists; their original content
cannot be cryptographically verified. Keep those files and logs immutable.

After inspecting the report and checking each campaign's copied
`voms_proxy.txt` with `voms-proxy-info -file ... -timeleft`, submit eligible
retries explicitly:

```bash
../py-env/bin/python retry_failed_io_jobs.py condor/CAMPAIGN --submit --max-retries 3
```

Retries retain the original input list, arguments, output path, and resources.
Logs are archived under `retry_archive/<timestamp>/<job>/`; incomplete outputs
are preserved as `.root.invalid-<timestamp>`. Completed jobs are untouched.
The archive count also bounds legacy retries. The per-submit `.recovery.json`
journal records attempts, generated submit descriptions and returned scheduler
IDs. Do not delete this journal or the archives while recovering a campaign.
A per-submit advisory lock serializes concurrent recovery scans, including
scans invoked from different parent directories. Submit only through this tool
while recovering: manually submitting generated descriptions bypasses its guards.

A crash or ambiguous `condor_submit` result leaves an `uncertain` journal entry
and blocks further retries. Reconcile its generated description and scheduler
records manually before marking it `submitted` with its `job_id`, or removing
that entry only after proving no job was submitted. A confirmed submission stays
blocked until its own terminal event arrives; old logs cannot trigger duplicates.
Submission requires a successful live queue query. Off-cluster dry runs use logs
and explicitly warn that live status could not be checked.

### Wall-time recovery without changing MC partitions

MC defaults to `microcentury`; data defaults to `espresso`. `--job-flavour`
overrides this in the original submitter. File-count splitting (default 10) is
configurable, but does not balance events, remote-I/O latency or runtime.

Do not automatically subdivide existing MC jobs: pPb MC weights depend on the
invocation entry count, and the controlled-sample RNG restarts per invocation.
Both affect merged physics results when partitions change. Event-range processing
is not exposed by this executable. No physics code or partition semantics have
been changed by recovery.

Instead, explicitly choose a **strictly longer** flavour for failed wall-time
units after reviewing the plan, for example:

```bash
../py-env/bin/python retry_failed_io_jobs.py condor/CAMPAIGN \
  --walltime-flavour longlunch --report /tmp/walltime-plan.json
../py-env/bin/python retry_failed_io_jobs.py condor/CAMPAIGN \
  --walltime-flavour longlunch --max-retries 3 --submit
```

Only wall-time retries get the selected runtime. Other parameters and input
coverage remain identical. A second wall-time failure cannot use the same
flavour again. Unknown original flavours, explicit MaxRuntime overrides, no longer available allowance, and
exhausted retry limits require intervention. This is a bounded recovery path;
it does not guarantee that a slow single file will finish. Changing partitions
or introducing event ranges requires a separate physics review/approval.

Rerun the same audit after retries finish. Merge only after the audit returns 0
and the summary accounts for all expected rows. Retain campaign records until
after verification; deleting per-job files during merging makes later audits
report them missing. Clone campaigns use a different queue schema and are not
supported by this processing-output audit.

Run the controlled recovery and worker tests from the repository root:

```bash
py-env/bin/python -m unittest discover -s tests -p 'test_condor*.py' -v
```

Tests simulate scheduler submissions and exercise real ROOT validation when the
ROOT CLI is installed. They do not submit cluster work.
