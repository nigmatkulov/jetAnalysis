# Compiled regularization scan kernel

`regularization_scan.cpp` is the C++17 numerical kernel used by
`05_unfold2D_regularization.ipynb`. It remains a standalone CMake project and is
not part of the repository's top-level build.

The kernel accepts plain row-major arrays, advances each toy from iteration 1
through the configured maximum exactly once, and accumulates only the per-bin
sums and sums of squares needed for ensemble bias and variance. Toy chunks run
in worker threads, but no ROOT histogram or RooUnfold object is shared between
threads. Complete unfolded-toy persistence intentionally remains outside this
kernel so ROOT I/O cannot serialize the numerical hot path.

`tests/test_compiled_regularization.py` validates every iteration of a seeded
ensemble bin by bin against `RooUnfoldBayes`, including the matched-only
response convention and external efficiency correction.

## Standalone build

The adjacent CMake project builds a shared library without modifying the main
analysis executable or requiring ROOT:

```sh
cmake -S hist_analysis/cpp -B build/regularization-scan \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build/regularization-scan --parallel
```

The resulting library is named `libjetanalysis_regularization_scan` with the
platform-specific shared-library suffix. The Python bridge searches this build
directory by default. Set `JETANALYSIS_REGULARIZATION_SCAN_LIBRARY` to load a
library from another location.
