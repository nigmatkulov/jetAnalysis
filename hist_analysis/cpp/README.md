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

For response probability `R(i,j)` and prior `p(j)`, each step forms
`d(i) = sum_j R(i,j)p(j)` and applies the inverse probability
`B(j,i) = R(i,j)p(j)/d(i)` to the measured toy. The unfolded contents become
the normalized prior for the next step. After the external efficiency
correction, the kernel adds each bin to the iteration's sum and squared sum.
Every worker owns a disjoint toy chunk and returns partial raw moments. The
kernel merges those moments across all workers in a fixed order; only the
global totals and total toy count are passed to Python. The Python layer then
calculates the ensemble mean, unbiased variance, bias squared, and MSE. It does
not average independently calculated per-worker statistics. A scan minimum at
the configured final iteration is boundary-limited; extend the maximum until
the minimum is bracketed by evaluated higher-iteration points.

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
