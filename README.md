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

## Histogram analysis

The `hist_analysis/` package provides a lightweight PyROOT workflow for
inspecting and comparing analysis outputs:

- `notebooks/01_inventory.ipynb`: checks configured directories, lists ROOT files, and inventories object keys.
- `notebooks/02_event_histograms.ipynb`: plots event-level pt-hat and vertex-z distributions and reproduces the gen/reco overweight-protection upper-tail diagnostic.
- `notebooks/02_jet_efficiency_fakes.ipynb`: calculates Lab-frame inclusive-jet efficiency, matched rate, and fake rate without pre-division normalization. It draws two-dimensional maps, eta projections in configured jet-pT intervals, and pT projections in configurable eta intervals; matched and fake 1D rates use style indices 0 and 1 on shared log-y canvases.
- `notebooks/02_jet_selection.ipynb`: overlays unit-normalized Reco, matched-Ref, and Gen inclusive-jet eta distributions in the Lab and CM frames, with Reco/Gen and Ref/Gen ratios, for jet-ID, track-maximum, and no-selection stages in configured jet-pT intervals after displaying the underlying two-dimensional maps. Full-range or configurable eta-range integral normalization is supported; stored-yield and bin-width modes and regular or binomial ratio errors remain configurable.
- `notebooks/03_dijet_reco_to_gen_closures.ipynb`: compares configurable reconstructed dijet eta distributions with nominal Gen for the intervals defined by `hist_analysis.config.histograms.DIJET_PTAVE_BINS`. It draws full eta shapes and forward/backward ratios with separate upper- and lower-panel axis ranges, ROOT normalization/division, explicit Reco-red and Gen-blue styles, and the standard dijet-selection annotation. Outputs follow `<generator>_<direction>_full_etaCM_<eta*10>_ptave_<low>_<high>.pdf` and the corresponding `..._fb_etaCM_...` pattern.
- `notebooks/03_dijet_smearing_effect.ipynb`: compares configurable nominal Gen, smeared Gen, eta-dependent smeared Gen, and Reco dijet pseudorapidity distributions for one eta acceptance and every configured pTave interval. Density, forward/backward, ratio-to-Gen, and ratio-to-eta-dependent-smeared-Gen overlays are saved on separate single-panel canvases.
- `notebooks/03_jet_JES_JER.ipynb`: extracts inclusive-jet JES/JER, compares systematic variations, and draws response maps.
- `notebooks/04_systematics_beam_orientation.ipynb`: compares p-going, Pb-going, and combined beam-orientation eta projections for single jets and dijets, with direction-comparison ratio panels and same-direction frame overlays.
- `notebooks/05_unfold2D.ipynb`: builds a flattened dijet pTave-eta response, runs a Bayesian RooUnfold closure test with explicit fake handling, and produces gen/reco/unfolded eta overlays with reco/gen and unfolded/gen ratios for every configured pTave interval. ROOT results and PDF diagnostics are written to `output/unfold2D/`; set `SAVE_PNG = True` for matching PNGs.
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

After changes to the histogram-analysis Python modules, check their syntax with:

```bash
py-env/bin/python -m compileall -q hist_analysis/config hist_analysis/python
```

Generated build products, ROOT outputs, notebook checkpoints, and files beneath
`hist_analysis/output/` should not normally be committed.
