# AGENTS.md

Guidance for agents working in this repository.

## Project Stack

- This is a CERN ROOT-based CMS jet analysis project.
- The main analysis source is `processForestSimple.C`. It is both a ROOT macro, for CLING/CINT execution, and the compiled executable source.
- `CMakeLists.txt` builds one executable, `processForestSimple`, and discovers ROOT through `ROOTSYS` and `find_package(ROOT CONFIG REQUIRED)`.
- `JetCorrector.h` and `JetUncertainty.h` are header-only ROOT/C++ helpers for JEC and uncertainty text files.
- `aux_files/` contains calibration/correction text inputs for pp, pPb, and PbPb workflows. Treat these as data files; do not reformat them.
- `pPb8160_analyzeMc_Pbgoing.py` and `pPb8160_analyzeMc_pgoing.py` are sequential Python 3 batch runners around the compiled executable.
- Generated/local outputs currently live under `build/`, `macro/`, and ROOT files such as `*.root`; avoid committing or rewriting generated analysis outputs unless explicitly requested.

## Build Commands

Make sure ROOT is available in the shell environment before configuring:

```bash
echo "$ROOTSYS"
root-config --version
```

Configure and build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
```

For an optimized build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The executable path is:

```bash
./build/processForestSimple
```

## Run Commands

Compiled executable:

```bash
./build/processForestSimple input.root output.root 2 1 30 0 -99 0 2
```

File-list input is also accepted:

```bash
./build/processForestSimple filelist.txt output.root 2 1 30 0 -99 0 2
```

ROOT macro mode:

```bash
root -l -b -q 'processForestSimple.C("input.root", "output.root", 2, 1, 30, 0, -99, 0, 2)'
```

Preset pPb 8.16 TeV MC batch runners:

```bash
python3 pPb8160_analyzeMc_Pbgoing.py
python3 pPb8160_analyzeMc_pgoing.py
```

The Python runners expect the executable at `build/processForestSimple` and input forests under `~/cernbox/ana/pPb8160/...`. They create outputs in `macro/eta_shift/`.

## Runtime Arguments

`processForestSimple` requires exactly 9 runtime arguments after the program name:

1. input ROOT file or text file list
2. output ROOT file
3. `mcType`: `0` data, `1` embedding, `2` pythia
4. `isPbGoingDir`: `0` p-going, `1` Pb-going
5. `ptHatSample`
6. `jeuSyst`: `-1`, `0`, or `1`
7. `jerSyst`: `-99`, `-1`, `0`, or `1`
8. `triggerId`: `0` no trigger/MB, `1` jet60, `2` jet80, `3` jet100
9. `recoJetSelMethod`: `0` none, `1` `trkMaxPt/RawPt`, `2` `jetIdTightLeptVeto`, `3` `jetIdLoose`

## Coding Rules

- Keep changes narrowly scoped. This repository is analysis code with hardcoded physics selections, binning, weights, and correction filenames.
- Preserve the dual-use structure of `processForestSimple.C`: ROOT macro code inside `#if defined(__CINT__) || defined(__CLING__)`, compiled executable code in the `#else` branch.
- Use ROOT types and local style consistently: `TString`, `TChain`, `TH1D`/`TH2D`/`TH3D`, `Form`, and the existing helper structs.
- Avoid broad refactors of histogram definitions, branch setup, event selection, correction paths, or physics constants unless the task explicitly asks for them.
- Do not reformat `aux_files/**/*.txt`; their spacing and content are external calibration inputs.
- Python scripts should remain simple sequential runners using `pathlib`, explicit constants, and `subprocess.run`.
- Generated data such as `build/`, `macro/eta_shift/*.root`, `extraJEC_*.root`, caches, and local output directories should generally remain untracked.
- Be careful with the current working tree: there may be local edits and generated ROOT files. Do not revert user changes.

## Self-Check Instructions

After C++ changes:

```bash
cmake --build build -j
```

If CMake has not been configured yet, run:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
```

Check executable argument validation without needing data:

```bash
./build/processForestSimple
```

This should print usage and exit nonzero.

For analysis-behavior changes, run a small known input if available:

```bash
./build/processForestSimple input.root /tmp/processForestSimple_test.root 2 1 30 0 -99 0 2
```

Then inspect the output with ROOT:

```bash
root -l -b -q '/tmp/processForestSimple_test.root'
```

After Python runner changes:

```bash
python3 -m py_compile pPb8160_analyzeMc_Pbgoing.py pPb8160_analyzeMc_pgoing.py
```

Before handing off, also run:

```bash
git status --short
```

Confirm only intended source/documentation files changed, and call out any checks skipped because ROOT or input ROOT files are unavailable.
