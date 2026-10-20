# jetAnalysis

This repository contains the `processForestSimple` analysis used to process CMS jet forest ROOT files. The code reads input trees through ROOT `TChain`s, applies event selection and jet corrections, fills histograms, and writes the output to a ROOT file.

## Physics conventions and review

The optional local `PHYSICS_REVIEW.md` records the annotated execution map,
confirmed calibration/response/smearing issues, and statistical limitations;
it is not distributed with the repository.
The follow-up fixes stored truth misses for events with fewer than two selected
reco jets and makes pPb pT-hat selection/weight intervals consistently half-open.
JER variations now cache one shared draw per physical jet/replica before sorting.
MC JEC selection uses sample type: embedding uses Embedded and standalone
PYTHIA uses Unembedded payloads. pp PYTHIA uses the Pb-going pPb Unembedded
payload by convention. Job-dependent pPb normalization remains unresolved. Histogram titles alone do not imply bin-width or luminosity normalization.

## Project layout

- `processForestSimple.C`: main analysis macro and compiled program source.
- `CMakeLists.txt`: CMake build file that finds ROOT and builds the executable.
- `JetCorrector.cc` and `JetUncertainty.cc`: jet-correction and uncertainty helpers.
- `LinkDef.h`: ROOT dictionary declarations used by the compiled executable.
- `aux_files/`: correction files and supporting inputs.
- `build/`: out-of-source build directory created by CMake.
- `pPb8160_analyzeMc_Pbgoing.py` and `pPb8160_analyzeMc_pgoing.py`: preset batch runners for the compiled executable.
- `pPb8160_runner.py`: shared implementation used by the preset batch runners.
- `pp_analyzeMc_cp5tune.py`: local per-file runner for the pp CP5 sample, with optional parallel workers.
- `processing/`: pPb 8.16 TeV file lists, data-availability checks, and the
  HTCondor submission workflow.
- `hist_analysis/`: Python, PyROOT, and Jupyter tools for ROOT-file inventory and histogram comparisons.
- `AGENTS.md`: repository-specific working, physics, implementation, validation, and handoff instructions for coding agents.

## Build

Make sure ROOT is available in your environment and that `ROOTSYS` is set:

```bash
echo "$ROOTSYS"
root-config --version
```

### Release build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Startup rejects equal input/output paths, including relative paths with `.` or
`..`, and checks each filename in an input list. This simple path check does not
resolve symlinks or hard links. Input preflight checks equal entry counts across
the four required trees separately in each file; identical event ordering is
assumed.

The jet energy corrector builds a candidate-row grid for one- and two-coordinate
payloads at load time. Lookup retains inclusive bounds, original row precedence,
dependency clamping, formulas, and sequential correction levels. Unsupported or
oversized layouts and non-finite query coordinates use the original scan.
`ctest --test-dir build --output-on-failure` includes exact indexed/scan
comparisons at payload boundaries, neighboring floating-point values, and
synthetic overlapping bins.

Generator smearing caches each deterministic Gaussian draw per sorted generator
jet and replica, reusing it across the four scale families. Replica counts,
weights, histogram fill order, and smearing formulas are unchanged.

Input preflight requires a consistent rho layout across all accepted files:
`rho` takes precedence over `fixedGridRhoFastjetAll`, and the selected leaf must
be scalar `Float_t`. A change of name, type, shape, or presence stops the run.
When neither branch exists in any file, the existing zero-rho fallback remains.

The local output directory must already exist and be writable. The executable
reports output creation, object-write, flush, and close failures on stderr and
returns a nonzero status. ROOT handles remote output URLs at file creation.

Use Release for normal validation and production execution. If compilation or
execution fails, reconfigure in Debug and reproduce the problem:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
```

The executable is built as:

```bash
build/processForestSimple
```

## Optional ROOT MCP server for Codex

[`root-mcp`](https://github.com/mohamedelashri/root-mcp) exposes ROOT-file
inspection, histogramming, fitting, plotting, and optional native PyROOT tools
through the Model Context Protocol. Codex supports project-scoped MCP settings
in `.codex/config.toml` for trusted projects; see the
[Codex MCP documentation](https://developers.openai.com/codex/extend/mcp).

Keep this optional tool separate from `py-env`, which is the analysis
environment. From the repository root, create a dedicated environment and
install the server:

```bash
python3 -m venv --system-site-packages .root-mcp-env
.root-mcp-env/bin/python -m pip install --upgrade pip
.root-mcp-env/bin/python -m pip install root-mcp
```

For XRootD URLs, install the optional remote-I/O support instead:

```bash
.root-mcp-env/bin/python -m pip install 'root-mcp[xrootd]'
```

`--system-site-packages` lets the environment reuse PyROOT from the selected
ROOT-compatible Python installation. Native ROOT tools such as
`run_root_code`, `run_rdataframe`, and `run_root_macro` require PyROOT to be
importable by this same interpreter.
Verify that before enabling them:

```bash
.root-mcp-env/bin/python -c 'import ROOT; print(ROOT.gROOT.GetVersion())'
```

If that command fails, either create `.root-mcp-env` with the Python
interpreter supplied by the active ROOT installation, or remove
`--enable-root` from `.codex/run-root-mcp.sh`. The server can still inspect ROOT
files through its non-native backend.

The tracked `.codex/config.toml` registers the project-local `root_cern` STDIO
server:

```toml
[mcp_servers.root_cern]
command = "./.codex/run-root-mcp.sh"
cwd = "."
startup_timeout_sec = 60
tool_timeout_sec = 120
enabled = true
default_tools_approval_mode = "prompt"
```

The tracked `.codex/run-root-mcp.sh` launcher resolves the checkout path at
runtime, keeps the machine-specific environment out of Git, and passes that
path as both `--data-path` and `--allowed-root`. This limits local file access
to this checkout. Native execution can run Python or ROOT macro code, so review
tool calls before approval. Do not expose this local STDIO configuration as an
unauthenticated HTTP service.

Restart Codex after adding the configuration, trust the project if prompted,
and verify the connection:

```bash
codex mcp list
```

In the Codex terminal UI, `/mcp` should list `root_cern`. Ask it to call
`get_server_info` and confirm that `root_native_available` and
`root_native_enabled` are both `true` when native ROOT support is expected.
For a direct server smoke test outside Codex, run:

```bash
./.codex/run-root-mcp.sh
```

The command waits for MCP JSON-RPC input; press `Ctrl-C` after confirming that
it starts without an import or configuration error. If Codex cannot start the
server, check the absolute paths, run the PyROOT import test above, and raise
`startup_timeout_sec` for slow ROOT initialization.

## Run with ROOT

You can run the macro directly from ROOT and pass the same arguments used by the executable:

```bash
root -l -b -q 'processForestSimple.C("input.root", "output.root", 2, 1, 30, 0, 2, 0)'
```

Argument order:

1. input file or file list
2. output file
3. `mcType` (`0` data, `1` embedding, `2` pythia)
4. `isPbGoingDir` (`0` p-going, `1` Pb-going)
5. `ptHatSample` (selects the pPb MC range; accepted but unused for pp)
6. `triggerId` (`0` no trigger/MB, `1` jet60, `2` jet80, `3` jet100)
7. `recoJetSelMethod` (`0` none, `1` `trkMaxPt/RawPt`, `2` tight jet ID with lepton veto, `3` loose jet ID)
8. `vertexFilterSelection` (`0` `pVertexFilterCutdz1p0`, `1` `pVertexFilterCutGplus`, `2` `pVertexFilterCutVtx1`)
9. optional `collisionSystem` (`"pPb"` by default, or `"pp"`); selects frame transformations, MC branches, weighting, and event cuts
10. optional `jerSmearingMode` (`0` by default): `0` CMS hybrid, `1` forced CMS stochastic, or `2` truth-to-pseudo-reco diagnostic
11. optional `useEtaDependentJER` (`0` by default): `0` keeps base-resolution standard and eta-dependent `Extra` families in one run; `1` uses eta-dependent resolution for both
12. optional `residualMode` (`2` by default): data only, `2` new residual,
    `1` old residual, `0` no residual; ignored for MC.

Data retains the beam-dependent L2 correction in all three modes. Mode `1`
then applies `Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PF.txt`; mode `2`
applies `Summer16_07Aug2017GH_V11_DATA_L2L3Residual_AK4PF.txt` instead.
Mode `0` applies neither residual. For example, select old residuals with:

```bash
./build/processForestSimple input.root output_oldResiduals.root 0 1 0 1 2 0 pPb 0 0 1
```

The ROOT-call signature has the same final optional argument. Explicit output
paths are honored as supplied; automatically generated data output names add
`_oldResiduals` or `_noResiduals` for modes `1` and `0`. Default names and ROOT
object keys remain unchanged.


Existing eight-argument calls retain the pPb boost, weighting, selection, and
orientation conventions.
For pp, append `pp` to the executable command, or `"pp"` to the ROOT call:

```bash
./build/processForestSimple input.root output.root 2 1 0 0 2 0 pp
root -l -b -q 'processForestSimple.C("input.root", "output.root", 2, 1, 0, 0, 2, 0, "pp")'
```

### Jet-energy-resolution smearing

The default (`jerSmearingMode=0`) follows the CMS hybrid prescription. A match
must satisfy `deltaR < 0.2` and
`|pT_reco-pT_ref| < 3*resolution*pT_reco`; accepted matches use
`pT_smear = pT_ref + SF * (pT_reco - pT_ref)`, while other jets fall back to
the stochastic branch. Mode `1` forces the CMS
stochastic branch on the reconstructed jet,
`pT_smear = pT_reco * (1 + N(0,resolution) * sqrt(max(SF^2-1,0)))`. Mode `2`
retains the analysis-specific truth-to-pseudo-reco construction for generator-level
tests, `pT_pseudo = pT_ref * (1 + N(0,SF*resolution))`.

Stochastic normal deviates are generated deterministically from a hash of jet
pT, jet eta, event rho when available, run, lumi, event ID, and replica index.
The JER variation is not part of the hash. The initial corrected-pT jet index
is frozen before response sorting and is combined with the replica index.
All nine existing momentum branches are cached on each jet before selecting
any varied pair; response sorting therefore cannot change its random draw.
`RecoJet::refPtSmeared` stores unit-SF smearing with base resolution,
`refPtSmearedEta` stores unit-SF smearing with eta-dependent resolution, and
`refPtSmearedJerDef/Up/Down` store the scale-factor variations using the
configured resolution. The internal `jerMomenta` cache also retains x2, x0.9,
and eta-dependent Down/Default/Up branches. Hybrid/Stochastic unit-SF smearing
adds no stochastic width; mode 2 retains full-resolution smearing around truth.
The named fields describe the current replica (the last replica after mode 2).
With `useEtaDependentJER=0`, standard templates use the base resolution and
`Extra` templates use the eta-dependent resolution. With the switch set to 1,
both families use the eta-dependent resolution and corresponding momenta agree.
For good Hybrid matches, resolution affects match eligibility only: the matched
residual is scaled by SF, without an additional eta multiplier.
If a resulting correction is non-finite or non-positive, the correction factor
is reset to one. Forests without a supported rho branch use `rho=0` with a
runtime warning; EventID is required.

Official hybrid and stochastic production modes use one deterministic smear
realization per jet. The retained truth-to-pseudo-reco diagnostic keeps ten
replicas for compatibility with the existing generator-smearing studies.

The data/MC scale-factor arrays `fJerDef`, `fJerLow`, and `fJerHi` come from
[`Summer16_25nsV1b_MC_SF_AK4PF.txt`](https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Summer16_25nsV1b_MC/Summer16_25nsV1b_MC_SF_AK4PF.txt).
The base resolution and optional eta-dependent multiplier are extracted from
this analysis's MC response. The eta-independent mode uses only the central
pT-dependent resolution; the eta-dependent mode multiplies it by the fitted,
centrally normalized eta shape.

For every selected MC event, `processRecoDijets` constructs the Down, Default,
and Up jet momenta with the corresponding eta-binned scale factor. It reorders
the complete reconstructed-jet collection independently for each variation,
then reapplies the leading/subleading pT and delta-phi selections. The
`hRecoDijetPtEtaCMJer{Down,Def,Up}_<cut>` histograms store the varied dijet
pTave on x and the average reconstructed jet etaCM on y. Forward and Backward are derived from their full etaCM bins during analysis;
no separate hemisphere objects are stored. Matched-reference companions use the same varied reco pTave and selected
reco pair but fill reference dijet eta. The first Default realization also
fills the nominal JER unfolding response, miss, fake, and pair-classification
objects. MC output files record `jerSmearingMode` and `useEtaDependentJER` as
string-valued `TNamed` metadata so their meaning remains stable after `hadd`.

In pp mode, Lab and CM eta both equal stored detector eta: neither a boost
nor a beam-direction sign flip is applied, for data or MC. This covers Gen,
Reco, reference jets, dijets, systematics, and unfolding. The shift-scan
histograms keep their existing names/titles for output compatibility, but every
entry uses stored eta in pp mode, regardless of the shift printed in its name.
CM eta acceptance and forward/backward assignments use these unshifted values.

For pp MC, `hiEvtAnalyzer/HiTree/pthat` and `weight` must be scalar `Float_t`
branches. The stored `weight` is used directly for each event, including zero
and negative values; it is not divided by the file's event count or multiplied
by a vertex weight. Non-finite weights fail explicitly. `ptHatSample` is unused:
there is no pT-hat range cut. The pp event selection uses only `|vz| <= 15 cm`
(the pPb skim filters and vertex-filter selector are ignored). Existing jet
selection and pT-hat-dependent overweight-event rejection remain in effect.

Both pp and pPb MC read `hiEvtAnalyzer/HiTree/pthat` (lowercase). The pPb
pT-hat range cut and weight-table lookup both use `[low, high)`: an exact shared
boundary belongs to the next interval. Generator normalization, vertex weighting,
and event filters otherwise retain their existing behavior. The pPb weight still
divides by the input entries in the current invocation, so splitting a sample
into jobs and directly summing their outputs changes its normalization. Use a
complete sample per invocation until a shared normalization convention is defined.
Required MC branch names and scalar `Float_t` types are
checked in every input file.

Corrections, vertex orientation, default macro paths, and other analysis settings
remain the existing pPb settings; `isPbGoingDir` still controls those settings.
Unknown collision-system values are rejected before input processing.

## Run the compiled binary

After building, run the executable directly:

```bash
./build/processForestSimple input.root output.root 2 1 30 0 2 0
```

Example for a file list:

```bash
./build/processForestSimple filelist.txt output.root 2 1 30 0 2 0
```

A file list contains one ROOT path per line. Blank lines and comment lines are
ignored; candidate entries must contain both `Forest` and `.root`. Every listed
file should be readable and contain the four required analysis trees. Missing,
unreadable, or incomplete ROOT files are skipped with an explicit message on
the error stream, while the remaining files continue processing. Malformed
list entries, no usable remaining files, or inconsistent aggregate tree entry
counts terminate the job with a nonzero status. The filename does not need to
contain `AOD`.

Before submitting a new dataset to Condor, run one small list directly and
confirm that the executable reports a nonzero file count and the expected event
count:

```bash
./build/processForestSimple test.list /tmp/processForestSimple-test.root \
  0 1 0 1 2 0
```

## Preset batch scripts

The repository also includes two sequential batch runners:

- `pPb8160_analyzeMc_Pbgoing.py`
- `pPb8160_analyzeMc_pgoing.py`

Each script loops over its explicit `PT_HAT_SAMPLES` list and runs the compiled
executable. By default the jobs run one at a time, but you can use multiple local
cores with `--workers`:

```bash
py-env/bin/python pPb8160_analyzeMc_Pbgoing.py --workers 4
py-env/bin/python pPb8160_analyzeMc_pgoing.py --workers 4
```

Both runners pass the JER configuration explicitly. The default is CMS
Hybrid with `USE_ETA_DEPENDENT_JER=0`: a single run writes both base-resolution
standard JER distributions and eta-dependent `Extra` distributions, including
both unfolding response families. No second MC production is needed.
Leave `--eta-dependent-jer` off to keep these families distinct.
`--no-eta-dependent-jer` explicitly selects this default for the pPb runners.
The optional `--eta-dependent-jer` override makes both families eta-dependent;
corresponding cached momenta then agree. `--jer-smearing-mode {0,1,2}` selects
the smearing algorithm independently. Direction-specific scripts expose the
presets as `JER_SMEARING_MODE` and `USE_ETA_DEPENDENT_JER`; command-line options
take precedence when the scripts are run.

The runner keeps going after individual sample failures and prints a final summary of skipped and failed pt-hat jobs.
Inputs are read from
`$HOME/cernbox/ana/pPb8160/<generator>/<direction>/forest/`, and outputs are
written beneath `$HOME/cernbox/ana/pPb8160/<generator>/<direction>/`. The
pt-hat samples to run are the explicit `PT_HAT_SAMPLES` list in each preset
script; review that list before starting a production run.

## Local pp CP5 processing

`pp_analyzeMc_cp5tune.py` discovers every `*.root` file directly under
`$HOME/cernbox/ana/pPb8160/pythia/cp5tune/forest/` and processes one file per job.
It writes `<input-stem>_hist.root` and `<input-stem>_hist.log` under
`$HOME/cernbox/ana/pPb8160/pythia/cp5tune/`. On this workstation those directories
are `/Users/gnigmat/cernbox/ana/pPb8160/pythia/cp5tune/forest/` and its parent.
The runner invokes `processForestSimple` as Pythia MC in `pp` mode. It passes
direction 1 to select the Pb-going pPb Unembedded MC JEC payload, as intended
for pp PYTHIA. The selector tests `mcType == 1` for embedding; `mcType == 2`
uses Unembedded. pp does not have a separate calibration payload.
The direction does not alter pp eta. It uses tight jet ID and trigger ID 0.
These analysis choices are constants near the top of the script. JER mode and
eta dependence are selected with the same `--jer-smearing-mode` and
`--eta-dependent-jer` options as the pPb runners.

Run the complete local pp sample as follows:

```bash
cd /Users/gnigmat/work/cms/jetAnalysis
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# Confirm the discovered inputs and generated commands without processing.
py-env/bin/python pp_analyzeMc_cp5tune.py --workers 4 --dry-run

# Process all discovered forests with four concurrent jobs.
py-env/bin/python pp_analyzeMc_cp5tune.py --workers 4
```

For example, `HiForestMiniAOD_1.root` produces:

```text
/Users/gnigmat/cernbox/ana/pPb8160/pythia/cp5tune/HiForestMiniAOD_1_hist.root
/Users/gnigmat/cernbox/ana/pPb8160/pythia/cp5tune/HiForestMiniAOD_1_hist.log
```

Inspect failures in the matching `.log` files. A successful campaign ends with
`Completed: N succeeded, 0 failed.`, where `N` is the number of discovered inputs. The script returns a nonzero exit status
if any file fails.

The default is one worker. Choose `--workers` according to available memory and
I/O bandwidth; each worker is a separate ROOT analysis process. `--input-dir`,
`--output-dir`, and `--executable` override the paths; `--dry-run` only prints
commands. Existing outputs or logs
are rejected unless `--overwrite` is provided. Each subprocess runs from the
repository root to resolve correction files. Individual processing failures are
reported with their log paths; remaining files continue, and any failed job
makes the runner exit nonzero. Outputs are per input file and are not merged.

To intentionally rerun the complete campaign and replace existing per-file
outputs and logs, use:

```bash
py-env/bin/python pp_analyzeMc_cp5tune.py --workers 4 --overwrite
```

To process one pp forest manually, bypassing the batch runner, use the compiled
executable directly. The fifth numeric argument is `0` because the pT-hat sample
selector is unused in pp mode; the final `pp` argument enables the pp branches,
stored event weights, event selection, and eta convention. Use an unused output
path when calling the executable directly because the runner's overwrite guard
does not apply:

```bash
./build/processForestSimple \
  /Users/gnigmat/cernbox/ana/pPb8160/pythia/cp5tune/forest/HiForestMiniAOD_1.root \
  /Users/gnigmat/cernbox/ana/pPb8160/pythia/cp5tune/HiForestMiniAOD_1_hist.root \
  2 1 0 0 2 0 pp
```

## HTCondor processing

The `processing/submit_process_forest_condor.py` workflow splits the tracked
pPb 8.16 TeV file lists, creates readable output and log names, and invokes
`processForestSimple` with the complete current interface. It passes `pPb`, the
JER mode, and the eta-dependence switch explicitly. The submitter exposes these
as `--jer-smearing-mode` and `--eta-dependent-jer`. By default, one MC pass
keeps base standard and eta-dependent `Extra` families; omit
`--eta-dependent-jer` to retain both distinct families. It supports
individual or all MB primary datasets and individual or all MC pT-hat samples.
Prepare
`processing/voms_proxy.txt`, then follow the data and MC production examples in
`processing/README.md`. The general submitter uses 10 input files per job by
default so one transient input failure invalidates a smaller unit of work.
Both Condor submitters accept `--residual-mode {0,1,2}` (default `2`). Data
job labels add `_oldResiduals` or `_noResiduals` before `_jobNNNN`; new-residual
names remain unchanged. Use the same `--residual-mode` with
`processing/merge_data_outputs.py` to select and name the matching merged output.

`processing/submit_all_data_condor.py` is the convenience wrapper for submitting
all four data triggers in both beam orientations; it uses 50 input files per
job by default.
The default CERN job flavour is `espresso` (20 minutes) for experimental data
and `microcentury` (1 hour) for MC; `--job-flavour` overrides these defaults.
To protect the 2 GB AFS quota, successful Condor stdout is discarded. Worker
output is captured in node-local temporary storage, and only the final 64 KiB
is written to a per-job `.err` when a job fails. Small per-job scheduler `.log`
files are retained because the targeted ROOT I/O retry workflow uses their
termination records; live output streaming and email notifications are disabled.

Standard production ROOT outputs are written beneath
`/eos/user/g/gnigmatk/ana/pPb8160/`: experimental data under `exp/<direction>`,
embedding under `embedding/<direction>`, and Pythia under
`pythia/<direction>`, where `<direction>` is `Pbgoing` or `pgoing`.

Use `--dry-run` to generate and inspect all sublists and submit descriptions
without calling `condor_submit`. Run submission commands from `processing/`;
that directory convention and the complete production/test commands are
documented in [processing/README.md](processing/README.md).
A dry run writes preparation files and proxy copies; it does not run the
analysis or validate input availability. Output names omit the campaign ID,
input partition and JER settings, so use separate output directories for
configurations that must coexist.

After jobs finish, `processing/retry_failed_io_jobs.py` can scan generated
campaigns for ROOT input read failures and resubmit only the affected queue
rows. It reports a dry run by default; confirm its selected total and the copied
campaign proxy before using the explicit `--submit` mode. The archive,
quarantine, campaign-scoping, and post-retry checks are documented in
`processing/README.md`.

`processing/merge_data_outputs.py` recursively merges one completed data
trigger/direction in batches using `hadd -f208`, verifies the final ROOT file,
atomically installs it, and then removes the merged per-job ROOT files. See
`processing/README.md` for MB and jet-trigger examples, validation details, and
the `jetId`/`trkMax`/`noSel` selection mapping, custom output naming, and
non-destructive options.

Before a data campaign, run `processing/check_data_filelists.py` on a CERN host
with EOS mounted. It checks every MB PD list in both directions and both
direction-specific PAEGJet lists, prints missing ROOT paths, and returns a
nonzero status when files or required lists are unavailable.

`hRecoInclusiveJetRawPtEtaLabUnflipped` stores selected raw jet pT versus
detector eta using the same configured selection as the nominal corrected-pT
inclusive-jet histogram. It is written for both data and MC so the raw and
corrected inclusive-jet spectra can be compared directly. Data outputs also
retain the run-integrated and run-specific `hRecoDijetPtEtaCMRunArr` histograms
for the configured run-ID comparison. MC outputs also contain
`hRefSelInclusiveJet{PtEtaLabUnflipped,PtEtaLab,PtEtaCM}` and the corresponding
`TrkMax` and `NoSel` families, filled from reference jets matched to the
selected reconstructed-jet population. For every configured dijet eta cut, MC
outputs additionally contain nominal Lab-frame Reco, matched-Ref, and Gen
histograms named `h{Reco,Ref,Gen}DijetPtEtaLab_<index>`. The pointing-resolution
workflow uses the matched-pair `h{Reco,Ref}DijetPtEta{CM,Lab}PointingRes_<index>`
pairs. Both histograms retain the reconstructed dijet pTave and require matched
jets. The Reco histogram applies its eta acceptance to reconstructed jets and
the Ref histogram applies it to the matched reference jets.
The configured per-jet acceptances are 1.5, 1.6, 1.7, 1.8, 1.9, 2.4, 2.5,
and 3.0 (indices 0–7); the standard 1.9 acceptance is index 4. Older
nine-cut outputs beginning at 1.4 must be regenerated before use with the
current shared cut table.

The optional mixed-prior histograms are disabled by default. Configure with
`-DDO_MIXED_PRIOR=ON` to compile their constants, allocation, filling, and ROOT
output back into `processForestSimple`.

## Histogram analysis

The `hist_analysis/` package provides a lightweight PyROOT workflow for
inspecting and comparing analysis outputs:

- `notebooks/01_inventory.ipynb`: checks configured directories, lists ROOT files, and inventories object keys.
- `notebooks/02_mc_event_histograms.ipynb`: plots event-level pt-hat and vertex-z distributions and reproduces the gen/reco overweight-protection upper-tail diagnostic.
- `notebooks/02_mc_jet_efficiency_fakes.ipynb`: calculates Lab-frame inclusive-jet efficiency using matched reference pT and eta (`hRefSelInclusiveJetPtEtaLab / hGenInclusiveJetPtEtaLab`), matched rate, and fake rate without pre-division normalization. It draws the matched-reference 2D yield distribution, two-dimensional rate maps, eta projections in configured jet-pT intervals, and pT projections in configurable eta intervals; matched and fake 1D rates use style indices 0 and 1 on shared log-y canvases.
- `notebooks/02_mc_jet_selection.ipynb`: overlays Reco, matched-Ref, RefSel, and Gen inclusive-jet eta distributions in the Lab and CM frames for jet-ID, track-maximum, and no-selection stages in configured jet-pT intervals after displaying the underlying two-dimensional maps. Full-distribution overlays support integral, eta-range, stored-yield, and bin-width normalization, with the `LOG_Y` option controlling the y-axis scale. Reco/Gen, Ref/Gen, and RefSel/Gen ratios are drawn together on individual canvases using independent-error propagation by default; CM canvases mark `|eta_CM|=1.9` and `2.4` with dashed and solid black lines.
- `notebooks/02_mc_jet_JES_JER.ipynb`: extracts inclusive-jet JES/JER, compares systematic variations, and draws response maps.
- `notebooks/02_mc_jet_JEC_check.ipynb`: compares raw and corrected MC inclusive-jet eta projections in configurable pT intervals for Pb-going and p-going samples. It writes square four-distribution overlays and two-panel corrected/raw and beam-direction double-ratio figures, with configurable eta and panel ranges and ROOT ratio-error options.
- `notebooks/03_mc_dijet_embedding_vs_pythia.ipynb`: compares unit-normalized Reco, matched-Ref, and Gen dijet CM-pseudorapidity shapes between Embedding and PYTHIA using the selected standard or test pTave binning, with PYTHIA/Embedding ratio panels. Optional CP5 input supports `cernbox`, `local`, and `auto`; `local` reads `$HOME/cernbox/ana/pPb8160/pythia/cp5tune/cp5_gen_dijet_histograms.root` and contributes to the Gen comparisons because the standalone file is generator-only. It also compares the corresponding unnormalized forward/backward ratios and their PYTHIA/Embedding ratios using independent-error propagation. Each level runs from its own cell, and each overlay family uses a common y-axis range across all levels and pTave intervals.
- `notebooks/03_pythia_standalone_check.ipynb`: compares standalone pp 8.16 TeV CP5 and CUEP8M1 generator outputs, preferring notebook 00 products with original-source fallback. It draws inclusive-jet and dijet eta shapes and F/B ratios built from unnormalized yields with independent errors.
- `notebooks/03_mc_dijet_reco_to_gen_closures.ipynb`: compares configurable Reco and matched-Ref dijet eta distributions with nominal Gen for the intervals defined by `hist_analysis.config.histograms.DIJET_PTAVE_BINS`. It draws raw full-distribution overlays and ratios before normalization, normalized full shapes, and unnormalized forward/backward ratios. Forward/backward construction always uses independent error propagation.
- `notebooks/03_mc_dijet_reco_smeared_to_gen_closures.ipynb`: compares eta-dependent JER-smeared Reco with nominal Gen using normalized full distributions and unnormalized forward/backward ratios, with separately configurable full-distribution and F/B-to-F/B ratio errors.
- `notebooks/03_mc_dijet_smearing_effect.ipynb`: always compares nominal Gen, eta-dependent smeared Gen, and Reco dijet pseudorapidity distributions; other curves can be enabled from a catalog. It writes raw full-distribution overlays and ratios after correcting accumulated smear trials, normalized full-shape comparisons, and forward/backward comparisons for each configured pTave interval. Constructing F/B always uses independent errors; subsequent comparison-ratio options are configurable.
- `notebooks/04_systematics_beam_orientation.ipynb`: compares p-going, Pb-going, and combined beam-orientation eta projections for single jets and dijets, with direction-comparison ratio panels and same-direction frame overlays.
- `notebooks/04_systematics_JER.ipynb`: estimates the reconstructed-dijet JER shape systematic in embedding or Pythia MC from unit-integral JER Up, Down, and default CM eta projections and unnormalized forward/backward ratios. It draws dedicated Up/Def and Down/Def comparisons with configurable polynomial fits, calculates bin-by-bin relative uncertainties from the fitted values, and optionally applies outward running-maximum smoothing to the full CM and F/B uncertainties; F/B construction always uses independent errors.
- `notebooks/04_systematics_JEU.ipynb`: estimates the reconstructed-dijet JEU shape systematic for the MB and jet-triggered data samples from unit-integral JEU Up, Down, and default CM eta projections and unnormalized forward/backward ratios. Configurable polynomial fits provide the bin-center values used to calculate the relative uncertainty, with optional outward smoothing for full-CM and F/B results.
- `notebooks/04_systematics_pileup.ipynb`: estimates the reconstructed-dijet pileup-filter shape systematic for the MB and jet-triggered data samples by comparing unit-integral `Gplus` and optional diagnostic `Vtx1` CM eta projections with nominal `dz1p0`, together with their unnormalized forward/backward ratios. The symmetric uncertainty is `|Gplus/dz1p0 - 1|`, using fitted or direct bin-by-bin values and optional outward smoothing.
- `notebooks/04_systematics_pointing_resolution.ipynb`: estimates the MC dijet pointing-resolution shape systematic in the CM and Lab frames from unit-normalized Reco and pointing-reference projections. It extracts `|Reco/Ref - 1|` using fitted or direct bin-by-bin values and optional outward smoothing; CM additionally includes F/B ratios folded from unnormalized full projections, while Lab uses full distributions only.

Unfolding population definitions and current statistical limitations are
documented in [the histogram workflow](hist_analysis/README.md#direct-2d-unfolding-populations).

- `notebooks/05_unfold2D_basics.ipynb`: builds a flattened dijet pTave-eta response for the eta-dependent JER-default measured distribution, runs a Bayesian RooUnfold closure test with explicit miss/fake handling, and produces Gen/Reco/unfolded eta overlays and closure ratios for every configured pTave interval. Set `COMPARISON_TARGET` to `gen` or `ref`, and use `PTAVE_BIN_SET` to select the shared test or standard pTave intervals. ROOT results and PDF diagnostics are written to `output/unfold2D/`; set `SAVE_PNG = True` for matching PNGs.
- `notebooks/05_unfold2D_data.ipynb`: unfolds trigger-specific reconstructed data with a coherent MC response, combines the trigger intervals using the trigger-stitching weights, and compares the unfolded shapes with configured generator predictions. See the histogram workflow documentation for response and stitching configuration.
- `notebooks/05_unfold2D_mc_direction.ipynb`: independently selects the coherent response production and the MC Gen/Reco production to unfold, including direction-aware Embedding/PYTHIA and the full pp CP5 production. It preserves the response scaling, miss/fake handling, diagnostics, and output conventions of `05_unfold2D_basics.ipynb`.
- `notebooks/05_unfold2D_studied_to_controlled.ipynb`: builds the response from the reproducible 60% controlled MC subsample, unfolds the independent 40% studied JER-default Reco spectrum, and compares it with studied Gen truth. Reco/Gen, unfolded/Gen, and forward-folded/reco ratio-error propagation can be configured independently as standard or binomial.
- `notebooks/05_unfold2D_gen_and_prior.ipynb`: compares nominal Gen with the mixed two-Gaussian prior in either the shared test or standard pTave intervals, with optional raw-yield canvases and normalized shape ratios.
- `notebooks/05_unfold2D_regularization.ipynb`: scans Bayesian iterations with independently selectable response-training MC (Embedding, PYTHIA, or full CP5) and truth (response Gen, mixed prior, direction-aware Embedding/PYTHIA, standalone CUEP8M1 or CP5, full CP5, an arbitrary compatible ROOT histogram, or a generated Gaussian). File and histogram-key overrides make the truth choice independent of the built-in catalog. One global scale factor places the complete toy truth on the MinimumBias exposure using the direct MB--Jet60 overlap recorded in the trigger-stitching artifact; individual pTave blocks are not tuned. Trigger-stitched MB, Jet60, Jet80, and Jet100 distributions then provide per-bin relative statistical uncertainties. If MC is populated in an empty data bin, the nearest usable relative uncertainty from the same pTave block is used. Populated toy bins are Gaussian-smeared with the resulting absolute errors, while originally empty bins remain zero. Its default C++17 backend is built independently with `cmake -S hist_analysis/cpp -B build/regularization-scan -DCMAKE_BUILD_TYPE=Release` and `cmake --build build/regularization-scan --parallel`; it advances every toy through the iteration range once and parallelizes toy chunks. Output names and plot annotations identify both roles.
- `notebooks/05_unfold2D_regularization_tutorial.ipynb`: a student-facing, step-by-step version of the regularization study. It plots the response-training prior and toy truth, purity and efficiency, MB normalization, forward-folded expectation, a selected fluctuated toy before and after purity correction, that toy after every Bayesian iteration, the sequential ensemble bias/variance/MSE scan, and the mean unfolded ensemble at every iteration. It deliberately uses repeated PyROOT/RooUnfold calls and never invokes the compiled parallel scan.
- `notebooks/06_data_jet_distributions.ipynb`: compares Pb-going, p-going, and combined inclusive-jet data for configurable jet-ID, track-maximum, or no-selection outputs. It draws common-scale orientation 2D maps and pT/eta overlays with ratios to combined in the unflipped-Lab, flipped-Lab, and CM frames, plus selected raw-jet-pT versus unflipped detector-eta diagnostics.
- `notebooks/06_data_dijet_distributions.ipynb`: provides the corresponding selection-configurable dijet orientation maps and pTave/eta comparisons using the full unflipped-Lab, flipped-Lab, and CM dijet histograms, plus CM forward/backward comparisons and ratios to the combined orientation. Eta projections use unit-area bin-width normalization for true `1/N dN/deta` densities. It also collects the 16 selected combined-orientation intervals into separate square 4x4 CM-shape and F/B summary figures.
- `notebooks/06_data_trigger_stitching.ipynb`: derives trigger-overlap weights and writes the ROOT/JSON/CSV artifacts used by data unfolding and regularization studies.
- `notebooks/06_data_runId_dependence.ipynb`: compares the configured run-specific data histograms with the run-integrated distributions.
- `notebooks/07_dijet_data_vs_mc.ipynb`: prepares configured Data and optional official gen/reco, CP5 and CUETP8M1 references once, then writes separate single-panel overlay and Data/MC figures for each enabled comparison. Independent `DRAW_OFFICIAL_MC`, `DRAW_CP5`, `DRAW_CUETP8M1`, `DRAW_EPPS21`, and `DRAW_DATA_SYSTEMATIC` switches control samples and stored-total systematic brackets. Full eta, raw F/B and normalized pTave comparisons use Data/MC throughout (the pTave direction change is explicitly approved). ROOT/MPLHEP share styles, clean physics-only selection text and accepted eta-bin edges. EPPS21 uses stored dimensionless pPb/pp bands with explicit observable labels. Combined model selections, configurable ratio operands, and independent full/ratio/F/B y limits are documented in `hist_analysis/README.md`.
- `config/`: input locations, common pseudorapidity range, configured single-jet pT and dijet pTave projection bins, and the standard dijet eta-cut index.
- `python/`: reusable ROOT object I/O, projection, histogram-operation, inventory, plotting, and flattened RooUnfold helpers. `unfolding.py` owns common projection, flattening, response, diagnostic, Bayesian-unfolding, and ROOT-output logic; `unfolding_plots.py` owns the corresponding diagnostics and closure figures.
- `output/`: generated figures and ROOT artifacts; this directory is ignored by Git.

Run `hist_analysis/notebooks/00_produce_reusable_histograms.ipynb` to build the
shared ROOT inputs before plotting or unfolding. Consumers prefer its files,
retain source fallback, and can project copied multidimensional inputs in
memory. See the [histogram workflow documentation](hist_analysis/README.md)
and [cut and naming review](hist_analysis/REUSABLE_HISTOGRAM_REVIEW.md).

Reusable projection keys encode the physical per-jet acceptance, such as
`etaCut1p5` and `etaCut1p9`, instead of an opaque cut index. Readers retain
compatibility with older reusable files; copied upstream ROOT keys remain
unchanged. Regenerate reusable bundles with notebook 00 to adopt the new names.

Dijet Forward and Backward histograms are no longer stored by the producer,
including nominal, systematic, smearing, and eta-shift families. Analysis folds
raw full-eta bins through `hist_analysis/python/forward_backward.py`, copies
contents and errors, and uses ROOT independent-error division for F/B. Older
ROOT files remain usable through their full histograms; stored hemisphere keys
are ignored. Exact-edge entries can differ from the former event-by-event
`abs(eta)` filling because ROOT's reflected bin intervals have the opposite
edge convention. This small, explicitly accepted difference does not change
full distributions, selections, weights, or responses. See the
[folding and validation details](hist_analysis/README.md#full-eta-forwardbackward-folding).

MC production also writes selected generator-level CM correlations of jet-based
`xPb`, `xp`, and `xPb/xp` with dijet eta and pTave. Existing selections and
histogram keys are preserved; see the definitions and new sparse keys in the
histogram workflow documentation.

Start Jupyter with the project Python environment from the repository root:

```bash
py-env/bin/python -m jupyter notebook
```

The default input locations are configured in `hist_analysis/config/files.py`.
Review them before running the notebooks on another machine or with a different
output layout.

## Validation

After C++ changes, configure and build in Release mode first:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/processForestSimple
ctest --test-dir build --output-on-failure
```

Running the executable without arguments should print its usage and exit with a
nonzero status. If compilation or execution fails, rebuild in Debug mode and
reproduce the failure before diagnosing it.

After changes to the batch runners, check their syntax with:

```bash
py-env/bin/python -m py_compile \
  pPb8160_analyzeMc_Pbgoing.py \
  pPb8160_analyzeMc_pgoing.py \
  pPb8160_runner.py \
  pp_analyzeMc_cp5tune.py
```

After changes to pp processing, run the runner checks and ROOT weight/selection
checks (these do not need forest inputs):

```bash
py-env/bin/python -m py_compile pp_analyzeMc_cp5tune.py
py-env/bin/python -m unittest discover -s tests -p test_pp_runner.py -v
root -l -b -q -e 'gROOT->ProcessLine(".L JetCorrector.cc"); gROOT->ProcessLine(".L JetUncertainty.cc"); gROOT->ProcessLine(".L processForestSimple.C"); gROOT->ProcessLine(".L tests/test_pp_processing.C"); test_pp_processing();'
```

After changes to the HTCondor workflow, check its syntax and generate a small
campaign without submitting it:

```bash
PYTHONDONTWRITEBYTECODE=1 py-env/bin/python -m py_compile \
  processing/submit_process_forest_condor.py \
  processing/submit_all_data_condor.py \
  processing/submit_clone_forest_condor.py \
  processing/merge_data_outputs.py \
  processing/retry_failed_io_jobs.py \
  processing/check_data_filelists.py
bash -n processing/run_process_forest_condor.sh
bash -n processing/run_clone_forest_condor.sh
py-env/bin/python processing/check_data_filelists.py --help
py-env/bin/python processing/submit_process_forest_condor.py \
  processing/filelists/pPb8160/DATA_PAEGJet/Pbgoing/PAEGJet_Pbgoing.txt \
  /tmp/jetAnalysis-condor-check --mc-type 0 --beam Pbgoing --trigger-id 1 \
  --proxy processing/voms_proxy.txt --dry-run
```

After changes to the histogram-analysis Python modules, check their syntax with:

```bash
py-env/bin/python -m compileall -q hist_analysis/config hist_analysis/python
py-env/bin/python -m unittest tests.test_unfolding tests.test_histogram_production -v
```

Generated build products, ROOT outputs, notebook checkpoints, and files beneath
`hist_analysis/output/` should not normally be committed.

The complete Python regression suite, including the compiled regularization
bridge and bounded Condor logging, can be run without production forests:

```bash
py-env/bin/python -m unittest discover -s tests -p 'test_*.py' -v
```

For response-pair classification and generator jet-x correlations, run the
synthetic ROOT checks from the repository root (no forest inputs required):

```bash
root -l -b -q -e 'gROOT->ProcessLine(".L JetCorrector.cc"); gROOT->ProcessLine(".L JetUncertainty.cc"); gROOT->ProcessLine(".L processForestSimple.C"); gROOT->ProcessLine(".L tests/test_unfolding_processing.C"); test_unfolding_processing(); gROOT->ProcessLine(".L tests/test_jer_cached_momenta.C"); test_jer_cached_momenta(); gROOT->ProcessLine(".L tests/test_unfolding_missing_reco.C"); test_unfolding_missing_reco(); gROOT->ProcessLine(".L tests/test_x_correlations.C"); test_x_correlations();'
```

These checks complement a small real-forest run; they do not validate remote
input availability or complete notebook execution.

To compare producer versions on the same forest input, keep the old and new
ROOT outputs and run this read-only check:

```bash
PYTHONPATH=. py-env/bin/python tests/compare_full_eta_outputs.py before.root after.root
```

It requires exact agreement of every retained object and that only hemisphere
keys were removed. It reports content/error differences in derived hemispheres
and F/B projections separately, including the accepted bin-edge differences.

To verify that both real Pb-going JEC payloads load and evaluate distinctly:

```bash
root -l -b -q -e 'gROOT->ProcessLine(".L JetCorrector.cc"); gROOT->ProcessLine(".L tests/test_jec_payloads.C"); test_jec_payloads();'
```

This payload check does not certify the main program's selector. Inspect its
`JEC files` log and accepted input paths separately; the current standalone
PYTHIA selector now chooses Unembedded JEC (P01 fixed); embedding retains
Embedded JEC.

Calibration input checks (missing files, oversized lines, truncated uncertainty
triples, and a valid payload):

```bash
root -l -b -q -e 'gROOT->ProcessLine(".L JetCorrector.cc"); gROOT->ProcessLine(".L JetUncertainty.cc"); gROOT->ProcessLine(".L tests/test_calibration_input_validation.C"); test_calibration_input_validation();'
```

Nonempty missing calibration paths and oversized lines fail explicitly. Empty
JEU paths remain the intentional unused-MC sentinel. These checks do not certify
all payload rows or the physics validity of extrapolation.

### Direct unfolding pair diagnostics

MC output includes six direct-fill dijet pTave/etaCM population TH2s, separating
matched selected pairs and found/not-found subdivisions of misses and fakes.
See [definitions and ROOT keys](hist_analysis/README.md#direct-2d-unfolding-populations).
Validate them without forest inputs:

```bash
root -l -b -q -e 'gROOT->ProcessLine(".L JetCorrector.cc"); gROOT->ProcessLine(".L JetUncertainty.cc"); gROOT->ProcessLine(".L processForestSimple.C"); gROOT->ProcessLine(".L tests/test_unfolding_pair_populations.C"); test_unfolding_pair_populations();'
```

### Extended unfolding

MC production additionally writes the same-event global-leading-pair extended
response, with explicit fiducial/outside-fiducial truth categories and complete
truth, measured, miss, fake and matched marginals. This includes rank/composition
migrations without changing the original geometrically matched response.
See [the extended producer contract and keys](hist_analysis/README.md#extended-leading-pair-response).
Validate the additional family with:

```bash
root -l -b -q -e 'gROOT->ProcessLine(".L JetCorrector.cc"); gROOT->ProcessLine(".L JetUncertainty.cc"); gROOT->ProcessLine(".L processForestSimple.C"); gROOT->ProcessLine(".L tests/test_unfolding_extended.C"); test_unfolding_extended();'
```
Existing notebooks require a new category-aware consumer to use these objects.
