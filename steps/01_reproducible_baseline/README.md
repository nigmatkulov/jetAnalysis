# Step 01: Reproducible software baseline

**Status:** In progress. An AppleClang Release build and one scoped Pb-going
embedded-MC run passed; the remaining reproducibility and four-channel exit
criteria are not complete.

## Goal and problem

Obtain a clean, documented build and deterministic small runs before interpreting
physics plots. The earlier CMake cache referenced compilers that were no longer
installed; the current Release cache uses `/usr/bin/cc` and `/usr/bin/c++`, and a
scoped Pb-going embedded-MC run has passed. RooUnfold is still resolved from a
configurable external checkout in the prototype notebook and is not part of the
CMake build.

## Recommended solution

1. Record OS, compiler, CMake, ROOT, Python, and RooUnfold versions.
2. Configure a fresh build directory without altering generated outputs or the
   analyst's current build until the replacement is proven.
3. Add focused automated checks for argument validation, convention helpers,
   bin boundaries, delta-phi wrapping, and histogram existence.
4. Validate the executable on a small subset or representative file from each
   existing embedded-MC direction and save logs plus a machine-readable
   configuration record. This is a software check, not a new production campaign.
5. Decide whether RooUnfold belongs in CMake, a dedicated Python environment, or
   both; remove absolute user paths from the reproducible workflow.
6. Verify that the executable/runner configuration can represent all four channels:
   Pb-going MC, p-going MC, Pb-going data, and p-going data. Small data runs may
   wait for approved inputs, but their configurations and dry-run tests must exist.

## Validation and exit criteria

- Clean configure and build succeed with documented commands.
- Python runners compile and a dry-run/configuration test passes.
- Repeated small runs agree exactly, or stochastic JER uses recorded seeds and
  agrees within a predefined statistical test.
- Required ROOT objects exist, have expected axes, and contain finite errors.
- The four channel configurations cannot overwrite or silently merge one another.
- Compiler warnings are reviewed and no warning relevant to correctness remains.

## Required figures

- One baseline event/jet cut-flow figure and one nominal gen/reco
  \((p_T^{\mathrm{ave}},\eta_{\mathrm{dijet}})\) figure, each PDF+PNG.

## Results

In progress. The stale compiler-cache failure was superseded by the AppleClang
Release configuration and successful Pb-going ptHat-120 execution recorded in
Step 00. A clean-checkout reproduction, repeatability test, warning review, and
the remaining channel configurations/runs are still required before this step is
complete.
