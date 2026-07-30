# jetAnalysis

This repository contains the `processForestSimple` analysis used to process CMS jet forest ROOT files. The code reads input trees through ROOT `TChain`s, applies event selection and jet corrections, fills histograms, and writes the output to a ROOT file.

## Project layout

- `processForestSimple.C`: main analysis macro and compiled program source.
- `CMakeLists.txt`: CMake build file that finds ROOT and builds the executable.
- `JetCorrector.cc` and `JetUncertainty.cc`: jet-correction and uncertainty helpers.
- `LinkDef.h`: ROOT dictionary declarations used by the compiled executable.
- `aux_files/`: correction files and supporting inputs.
- `build/`: out-of-source build directory created by CMake.
- `pPb8160_analyzeMc_Pbgoing.py` and `pPb8160_analyzeMc_pgoing.py`: preset batch runners for the compiled executable.
- `pPb8160_runner.py`: shared implementation used by the preset batch runners.
- `processing/`: pPb 8.16 TeV file lists and the HTCondor submission workflow.
- `hist_analysis/`: Python, PyROOT, and Jupyter tools for ROOT-file inventory and histogram comparisons.
- `steps/`: approval-gated dijet-measurement roadmap and validation records.
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

## Run with ROOT

You can run the macro directly from ROOT and pass the same arguments used by the executable:

```bash
root -l -b -q 'processForestSimple.C("input.root", "output.root", 2, 1, 30, 0, 2)'
```

Argument order:

1. input file or file list
2. output file
3. `mcType` (`0` data, `1` embedding, `2` pythia)
4. `isPbGoingDir` (`0` p-going, `1` Pb-going)
5. `ptHatSample`
6. `triggerId` (`0` no trigger/MB, `1` jet60, `2` jet80, `3` jet100)
7. `recoJetSelMethod` (`0` none, `1` `trkMaxPt/RawPt`, `2` tight jet ID with lepton veto, `3` loose jet ID)

## Run the compiled binary

After building, run the executable directly:

```bash
./build/processForestSimple input.root output.root 2 1 30 0 2
```

Example for a file list:

```bash
./build/processForestSimple filelist.txt output.root 2 1 30 0 2
```

A file list contains one ROOT path per line. Blank lines are not meaningful to
the C++ loader; candidate entries must contain both `Forest` and `.root`. Each
candidate is opened and accepted only when it is non-zombie and contains ROOT
keys. The filename does not need to contain `AOD`.

Before submitting a new dataset to Condor, run one small list directly and
confirm that the executable reports a nonzero file count and the expected event
count:

```bash
./build/processForestSimple test.list /tmp/processForestSimple-test.root \
  0 1 0 1 2
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

The runner keeps going after individual sample failures and prints a final summary of skipped and failed pt-hat jobs.
Inputs are read from
`$HOME/cernbox/ana/pPb8160/<generator>/<direction>/forest/`, and outputs are
written beneath `$HOME/cernbox/ana/pPb8160/<generator>/<direction>/`. The
pt-hat samples to run are the explicit `PT_HAT_SAMPLES` list in each preset
script; review that list before starting a production run.

## HTCondor processing

The `processing/submit_process_forest_condor.py` workflow splits the tracked
pPb 8.16 TeV file lists, creates readable output and log names, and submits the
current seven-argument `processForestSimple` executable. It supports individual
or all MB primary datasets and individual or all MC pT-hat samples. Prepare
`processing/voms_proxy.txt`, then follow the data and MC production examples in
`processing/README.md`.

Standard production ROOT outputs are written beneath
`/eos/user/g/gnigmatk/ana/pPb8160/`: experimental data under `exp/<direction>`,
embedding under `embedding/<direction>`, and Pythia under
`pythia/<direction>`, where `<direction>` is `Pbgoing` or `pgoing`.

Use `--dry-run` to generate and inspect all sublists and submit descriptions
without calling `condor_submit`. Run submission commands from `processing/`;
that directory convention and the complete production/test commands are
documented in `processing/README.md`.

## Histogram analysis

The `hist_analysis/` package provides a lightweight PyROOT workflow for
inspecting and comparing analysis outputs:

- `notebooks/01_inventory.ipynb`: checks configured directories, lists ROOT files, and inventories object keys.
- `notebooks/02_event_histograms.ipynb`: plots event-level pt-hat and vertex-z distributions and reproduces the gen/reco overweight-protection upper-tail diagnostic.
- `notebooks/02_jet_efficiency_fakes.ipynb`: calculates Lab-frame inclusive-jet efficiency, matched rate, and fake rate without pre-division normalization. It draws two-dimensional maps, eta projections in configured jet-pT intervals, and pT projections in configurable eta intervals; matched and fake 1D rates use style indices 0 and 1 on shared log-y canvases.
- `notebooks/02_jet_selection.ipynb`: overlays unit-normalized Reco, matched-Ref, and Gen inclusive-jet eta distributions in the Lab and CM frames, with Reco/Gen and Ref/Gen ratios, for jet-ID, track-maximum, and no-selection stages in configured jet-pT intervals after displaying the underlying two-dimensional maps. Full-range or configurable eta-range integral normalization is supported; stored-yield and bin-width modes and regular or binomial ratio errors remain configurable.
- `notebooks/03_dijet_reco_to_gen_closures.ipynb`: compares configurable Reco and matched-Ref dijet eta distributions with nominal Gen for the intervals defined by `hist_analysis.config.histograms.DIJET_PTAVE_BINS`. It draws raw full-distribution overlays and ratios before normalization, normalized full shapes, and unnormalized forward/backward ratios. Forward/backward construction always uses independent error propagation.
- `notebooks/03_dijet_reco_smeared_to_gen_closures.ipynb`: compares eta-dependent JER-smeared Reco with nominal Gen using normalized full distributions and unnormalized forward/backward ratios, with separately configurable full-distribution and F/B-to-F/B ratio errors.
- `notebooks/03_dijet_smearing_effect.ipynb`: always compares nominal Gen, eta-dependent smeared Gen, and Reco dijet pseudorapidity distributions; other curves can be enabled from a catalog. It writes raw full-distribution overlays and ratios after correcting accumulated smear trials, normalized full-shape comparisons, and forward/backward comparisons for each configured pTave interval. Constructing F/B always uses independent errors; subsequent comparison-ratio options are configurable.
- `notebooks/03_dijet_data_vs_mc.ipynb`: compares unit-normalized reconstructed-data, reconstructed-MC, and generator-level-MC dijet pTave spectra for each configured jet-eta acceptance and trigger sample.
- `notebooks/03_jet_JES_JER.ipynb`: extracts inclusive-jet JES/JER, compares systematic variations, and draws response maps.
- `notebooks/04_systematics_beam_orientation.ipynb`: compares p-going, Pb-going, and combined beam-orientation eta projections for single jets and dijets, with direction-comparison ratio panels and same-direction frame overlays.
- `notebooks/04_systematics_JER.ipynb`: estimates the reconstructed-dijet JER shape systematic in embedding or Pythia MC from unit-integral JER Up, Down, and default CM eta projections and unnormalized forward/backward ratios. It also draws dedicated Up/Def and Down/Def comparisons; F/B construction always uses independent errors.
- `notebooks/04_systematics_JEU.ipynb`: estimates the reconstructed-dijet JEU shape systematic for the MB and jet-triggered data samples from unit-integral JEU Up, Down, and default CM eta projections and unnormalized forward/backward ratios.
- `notebooks/05_unfold2D.ipynb`: builds a flattened dijet pTave-eta response for the eta-dependent JER-default measured distribution, runs a Bayesian RooUnfold closure test with explicit miss/fake handling, and produces Gen/Reco/unfolded eta overlays and closure ratios for every configured pTave interval. ROOT results and PDF diagnostics are written to `output/unfold2D/`; set `SAVE_PNG = True` for matching PNGs.
- `notebooks/05_unfold2D_mc_direction.ipynb`: trains the same flattened response on Pb-going embedding and applies it to the independent p-going embedding reco distribution, then compares the unfolded result with p-going generator truth. It preserves the response scaling, miss/fake handling, diagnostics, and output conventions of `05_unfold2D.ipynb`.
- `notebooks/06_data_check.ipynb`: plots reconstructed-data dijet pTave spectra, normalized CM eta projections, and unnormalized forward/backward ratios for configurable trigger-specific pTave intervals and jet-eta acceptances.
- `config/`: input locations, common pseudorapidity range, configured single-jet pT and dijet pTave projection bins, and the standard dijet eta-cut index.
- `python/`: reusable ROOT object I/O, projection, histogram-operation, inventory, and plotting helpers.
- `output/`: generated figures and ROOT artifacts; this directory is ignored by Git.

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
```

Running the executable without arguments should print its usage and exit with a
nonzero status. If compilation or execution fails, rebuild in Debug mode and
reproduce the failure before diagnosing it.

After changes to the batch runners, check their syntax with:

```bash
py-env/bin/python -m py_compile \
  pPb8160_analyzeMc_Pbgoing.py \
  pPb8160_analyzeMc_pgoing.py \
  pPb8160_runner.py
```

After changes to the HTCondor workflow, check its syntax and generate a small
campaign without submitting it:

```bash
PYTHONDONTWRITEBYTECODE=1 py-env/bin/python -m py_compile \
  processing/submit_process_forest_condor.py
bash -n processing/run_process_forest_condor.sh
py-env/bin/python processing/submit_process_forest_condor.py \
  processing/filelists/pPb8160/DATA_PAEGJet/Pbgoing/PAEGJet_Pbgoing.txt \
  /tmp/jetAnalysis-condor-check --mc-type 0 --beam Pbgoing --trigger-id 1 \
  --proxy processing/voms_proxy.txt --dry-run
```

After changes to the histogram-analysis Python modules, check their syntax with:

```bash
py-env/bin/python -m compileall -q hist_analysis/config hist_analysis/python
```

Generated build products, ROOT outputs, notebook checkpoints, and files beneath
`hist_analysis/output/` should not normally be committed.
