# Agent Guide

## Repository context

- CERN ROOT-based CMS jet analysis.
- `processForestSimple.C` is both a ROOT macro and the source of the compiled executable.
- `README.md` is the source of truth for layout, commands, runtime arguments, workflows, and validation.
- `hist_analysis/README.md` documents the notebook workflows in detail.
- `aux_files/` contains external calibration and correction inputs; treat it as data.

## Working agreement

- Keep the change limited to the user’s request. Do not fold in cleanup or refactors.
- Inspect the relevant code, `README.md`, and `git status --short` before editing.
- Preserve unrelated working-tree changes. Never revert or overwrite them.
- Preserve existing behavior and public interfaces unless the request explicitly changes them.
- Do not edit generated ROOT files, build products, caches, notebook checkpoints, or `hist_analysis/output/` unless requested.
- Avoid new dependencies unless they are necessary and justified.
- If a decision depends on an undocumented analysis convention, ask one focused question instead of guessing.

## Project review reports

- Save Markdown findings from project, code, architecture, physics analysis, documentation, performance, or similar reviews under the repository's `review/` directory, creating it if needed.
- Do not create review-result Markdown files in the repository root or other project directories. This applies to all future reviews of any project component.
- Review reports are temporary development artifacts and must not be committed to Git.

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

## Multi-agent collaboration

- The main agent defaults to GPT-6.1 Sol Low. Keep routine, low-impact tasks with one agent; use specialized subagents when their expertise or independent evidence justifies the cost.
- For scientifically significant changes, delegate implementation, independent physics review, and independent numerical validation to separate agents. Configure their models and effort independently; do not let an implementer certify its own work as independent validation.
- Available roles: `root-cpp-engineer` (Sol Low), `scientific-python-engineer` (Sol Low), `notebook-engineer` (Luna High), `physics-reviewer` (Sol Low), `validation-engineer` (Sol Low), and `batch-workflow-engineer` (Luna High). Their files under `.codex/agents/` define scope; the main agent retains integration and conflict resolution.
- Before scientific implementation, state the observable, provenance, units/frame, selections/weights, bin and flow treatment, normalization, uncertainty/correlation prescription, expected changes, and acceptance criteria. Ask about undocumented analysis choices before dependent work.
- Assign explicit file ownership: one writer per file, including notebooks and the producer. Serialize shared-file integration and ROOT output writes. Use isolated temporary outputs and avoid concurrent builds in the same directory. Do not assume a clean worktree contains the user's uncommitted changes.
- Give independent reviewers the requirements and source evidence before implementation explanations. Derive expectations using analytical examples, independently specified events, direct weighted fills or separate reference methods. Compilation and notebook execution alone are not scientific validation.
- Validate a frozen candidate identified by revision or relevant file hashes. Revalidate affected results after subsequent edits. Require exact discrete/schema comparisons where appropriate and justify floating-point tolerances and statistical closure criteria.
- Resolve conflicting findings using definitions and independent evidence; unresolved scientific assumptions require the user's decision. Report every verification as passed, failed or skipped, with exact reasons and material limitations.
- Keep geometric and extended response definitions separate. Preserve the F/B skill's independent-error prescription. Distinguish supplied theory-band validation from calculating PDF error members. Do not silently add covariance or change uncertainty conventions.
- Preserve campaign partitions: pPb normalization depends on invocation entry counts. Batch work defaults to dry runs; prepare merges with `--keep-inputs`. Production submission, job cancellation, deletion and replacement of validated outputs require explicit authorization for those operations.
- Keep at most three subagents active. Start with the smallest useful set and do not delegate recursively unless explicitly assigned. Reuse concise task-local code/test maps; refresh them when relevant files change. Return findings and evidence summaries instead of large logs.
- Escalate ambiguous or high-risk scientific implementation from Sol Low to Medium; use Luna XHigh only for a specific unresolved complex workflow problem. Explain the escalation and retain independent review. Custom role files take precedence over spawn effort overrides, so select a fresh explicitly configured agent when a different effort is needed rather than assuming an override changed a fixed role.
- Role instructions do not create filesystem isolation or grant permissions: subagents inherit the parent sandbox and tool permissions. A read-only reviewer is read-only by instruction, not a separate enforced sandbox.
