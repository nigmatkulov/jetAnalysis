# jetAnalysis

This repository contains the `processForestSimple` analysis used to process CMS jet forest ROOT files. The code reads input trees through ROOT `TChain`s, applies event selection and jet corrections, fills histograms, and writes the output to a ROOT file.

## Project layout

- `processForestSimple.C`: main analysis macro and compiled program source.
- `CMakeLists.txt`: CMake build file that finds ROOT and builds the executable.
- `aux_files/`: correction files and supporting inputs.
- `build/`: out-of-source build directory created by CMake.
- `pPb8160_analyzeMc_Pbgoing.py` and `pPb8160_analyzeMc_pgoing.py`: preset batch runners for the compiled executable.

## Build

Make sure ROOT is available in your environment and that `ROOTSYS` is set.

### Debug build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
```

### Release build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The executable is built as:

```bash
build/processForestSimple
```

## Run with ROOT

You can run the macro directly from ROOT and pass the same arguments used by the executable:

```bash
root -l -b -q 'processForestSimple.C("input.root", "output.root", 2, 1, 30, 0, -99, 0, 2)'
```

Argument order:

1. input file or file list
2. output file
3. `mcType` (`0` data, `1` embedding, `2` pythia)
4. `isPbGoingDir` (`0` p-going, `1` Pb-going)
5. `ptHatSample`
6. `jeuSyst`
7. `jerSyst`
8. `triggerId`
9. `recoJetSelMethod`

## Run the compiled binary

After building, run the executable directly:

```bash
./build/processForestSimple input.root output.root 2 1 30 0 -99 0 2
```

Example for a file list:

```bash
./build/processForestSimple filelist.txt output.root 2 1 30 0 -99 0 2
```

## Preset batch scripts

The repository also includes two sequential batch runners:

- `pPb8160_analyzeMc_Pbgoing.py`
- `pPb8160_analyzeMc_pgoing.py`

Each script loops over all available `ptHat` samples and runs the compiled executable one sample at a time.
