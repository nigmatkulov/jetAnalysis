# Agent Guide

## Repository context

- CERN ROOT-based CMS jet analysis.
- `processForestSimple.C` is both a ROOT macro and the source of the compiled executable.
- `README.md` is the source of truth for layout, commands, runtime arguments, workflows, and validation.
- `aux_files/` contains external calibration and correction inputs; treat it as data.

## Working agreement

- Keep the change limited to the user’s request. Do not fold in cleanup or refactors.
- Inspect the relevant code, `README.md`, and `git status --short` before editing.
- Preserve unrelated working-tree changes. Never revert or overwrite them.
- Preserve existing behavior and public interfaces unless the request explicitly changes them.
- Do not edit generated ROOT files, build products, caches, notebook checkpoints, or `hist_analysis/output/` unless requested.
- Avoid new dependencies unless they are necessary and justified.
- If a decision depends on an undocumented analysis convention, ask one focused question instead of guessing.

## Physics guardrails

Treat these as protected analysis behavior unless the user explicitly requests a change:

- event and jet selections, matching, and trigger logic;
- correction constants, correction-file paths, weights, thresholds, and systematic switches;
- histogram names, keys, dimensions, binning, axes, and output paths;
- ROOT output structure;
- response-matrix, fake, and miss definitions.

For an authorized physics change:

1. State the assumption and expected physics effect.
2. Make the smallest implementation change.
3. Validate the directly affected counts or objects.
4. Check that unrelated selections and yields do not change unexpectedly.

## Implementation invariants

### C++ and ROOT

- Preserve both execution modes in `processForestSimple.C`: CLING/CINT code under `#if defined(__CINT__) || defined(__CLING__)` and compiled code in the corresponding `#else` branch.
- Keep the executable linked with `JetCorrector.cc`, `JetUncertainty.cc`, and the ROOT dictionary generated from `LinkDef.h`.
- Preserve ROOT ownership and lifetime semantics; avoid copying large ROOT objects without need.
- Do not rename public functions, TTree branches, histogram keys, or output objects without explicit approval.
- Avoid formatting-only changes to `processForestSimple.C`; they hide physics-relevant diffs.

### Python and notebooks

- Use `py-env/bin/python` and run notebooks from the repository root.
- Keep sequential execution as the default; local parallelism remains opt-in.
- Treat uppercase notebook configuration cells as user-facing workflow parameters. Preserve their defaults unless the requested change requires otherwise.
- Prefer `pathlib.Path`, the file resolvers in `hist_analysis/python/histogram_io.py`, and the shared plotting conventions in `hist_analysis/python/root_style.py`.
- Preserve PyROOT ownership and lifetime semantics: clone objects that must outlive an input file, detach histograms with `SetDirectory(0)`, and close files reliably.
- Put newly introduced reusable analysis, projection, fitting, normalization, and plotting logic in `hist_analysis/python/`. Existing notebooks contain legacy inline implementations; do not refactor them unless required by the task.
- Keep machine-specific data and external-tool paths in configuration or environment variables. Do not introduce additional absolute paths in notebooks.
- Keep generated plots and derived ROOT files beneath `hist_analysis/output/` unless the user explicitly requests another location.

## Validation

Run the smallest applicable checks documented in `README.md`.

- C++ changes: configure and build Release first. Check executable argument handling when relevant.
- If C++ compilation or execution fails, reproduce it in Debug before diagnosing it.
- Python changes: run the documented syntax/compile check and the smallest affected workflow when practical.
- Notebook changes: verify that the notebook JSON is valid and execute the smallest affected workflow when its ROOT inputs and external libraries are available. Confirm that cells execute in order from a clean kernel and that the expected output files and ROOT keys are produced.
- Analysis changes: when a small known input is available, verify output creation, affected keys and dimensions, relevant event or jet counts, and unexpected yield changes.

Do not claim a check passed unless it was run. If data, ROOT, or another prerequisite is unavailable, report the skipped check and the exact reason.

## Handoff

Before finishing:

- run `git status --short` and confirm only intended files were touched by the agent;
- summarize modified files and behavior changes;
- report checks run and their results;
- report skipped checks, remaining uncertainty, and expected physics/output differences (or explicitly state that none are expected).
