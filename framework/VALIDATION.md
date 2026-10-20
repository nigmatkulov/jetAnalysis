# Numerical validation — 2026-10-20

## Scope and preserved reference

Development is on `macro2framework`, based on `macro` commit
`d934635819cc6bbcfa29eb830bcfb70d63d28c21`. `main` remains at
`6905d0b72076aa5126152f00377686bf301d16be`. The original macro is unchanged;
there are no migration commits or remote pushes. Existing untracked notebook
review documents and `hist_analysis/root_files/` were preserved.

Environment: macOS ARM64, AppleClang 21.0.0, ROOT 6.40.04. The framework and
reference use the same correction classes, payloads, ROOT installation, compile
configuration, inputs, argument order, and event ordering. No physics optimization
was introduced. A one-time textual check confirmed that all 17 migrated physics
function bodies (selection, weights, matching, nominal/systematic processing) and
the histogram creation body are unchanged. Input conversion, jet loading, loop
orchestration, and conditional output writing are the deliberately changed parts.

## Full real-embedding comparison

Input:
`~/cernbox/ana/pPb8160/embedding/pgoing/forest/HiForestSkim_pPb_MC_pthat540_pgoing_embedded.root`.
Arguments: `1 0 540 0 2 0`, default pPb/hybrid/base-resolution configuration.
The unmodified reference baseline was produced before physics migration.

| Quantity | Reference | Framework |
|---|---:|---:|
| Processed events | 502,400 | 502,400 |
| Selected events | 490,102 | 490,102 |
| Input reconstructed jets | 31,451,653 | 31,451,653 |
| Input generator jets | 7,686,551 | 7,686,551 |
| ROOT objects compared | 997 | 997 |
| Histogram cells / sparse bins compared | 18,457,544 | 18,457,544 |

`validate_forest_reader` independently used the macro's buffers as an oracle for
every event and every consumed reco/gen jet field. Metadata, filters, reference
associations, and kinematic inputs agreed exactly; floating jet inputs were
compared by representation. MC HLT is intentionally not read because the original
analysis bypasses it. The file has no supported rho branch; both implementations
use the unchanged zero fallback with a warning.

`compare_root_outputs` found exact agreement, without tolerances, in every retained
object: keys/cycles, concrete types, titles, axes and labels, fixed/variable
binning, entries, statistics, contents, errors, stored Sumw2, all flow cells,
response objects and THnSparse coordinates/content/error squares. Nominal,
JER, JEU where applicable, pointing, sample splits, and unfolding families are
included whenever produced by the selected configuration. ROOT serialization
metadata (file UUIDs, timestamps, compression layout) is deliberately not compared.

Initial artifacts:
`/tmp/jetAnalysis-macro2framework-reference/embedding-pgoing-540.{root,log}` and
`/tmp/jetAnalysis-framework-release/embedding-pgoing-540.{root,log}`.
The final analysis build also repeated this full-input comparison exactly, with
resources and comparison logs in `/tmp/jetAnalysis-framework-final-validation/`.
The resource benchmark runs also rerun exact output comparisons; their reports
are described in PERFORMANCE.md.

## Additional real-embedding and profiling check

A 2,000-entry aligned fixture was copied from the Pb-going embedding pT-hat 540
forest with `make_validation_fixture.py`. Repeating this same fixture ten times
produced 20,000 input events for the profiling run. This is not an independent
20,000-event sample and its normalization differs from full production.

Reference and instrumented framework both selected 19,480 events. They matched
exactly across 997 objects and 15,586,362 cells/sparse bins. Native profiling also
completed successfully, and its output was compared against the same reference.
Artifacts are in `/tmp/jetAnalysis-framework-profile-run/` and
`/tmp/jetAnalysis-framework-profile-native/`.

## Automated integration coverage

`framework/tests/test_workflow.py` creates two small forests and runs seven full
reference/framework comparisons:

- embedding, both orientations, pPb vertex selectors;
- hybrid, stochastic, truth-to-pseudo-reco modes and both eta-resolution settings;
- pp PYTHIA with positive, zero, and negative stored weights;
- pPb and pp data, three trigger selections, and all reco jet-ID choices across cases;
- missing reco/gen pairs, failed jet ID, reference association changes, and flow bins.

The same test checks event-level reader equivalence across file transitions,
unsigned event IDs above 2^63, both supported rho branch names in separate
invocations, optional tracks and vector particles, jets disabled in reader-only
mode, AK5PF/AK4CSPF aliases, rejection of non-AK4PF physics analysis, a malformed
second file, mixed rho-layout rejection, output-prefix suppression, and rejection
of an intentionally altered histogram. Profiling instrumentation is also checked
against an uninstrumented reference output.

The branch-contract unit test rejects missing leaves, wrong signedness/precision,
scalar/array mismatches, wrong array count leaves, and incorrect static shape.

Release and Debug CTest suites passed all five tests. A Release build with
`DO_MIXED_PRIOR=ON` also passed all five tests. The suites include the existing
analysis-contract test. Python syntax checks, pp runner checks (3), and Condor
logging checks (3) passed. The Condor dry-run generated a 610-job campaign and
forwarded framework arguments correctly; it submitted nothing. A clearly named
empty placeholder proxy was used only for that dry-run because this checkout
has no production proxy.

## Deliberate constraints and remaining validation limits

- By user decision, DiJetAnalysis supports AK4PF only until other collections'
  JEC/JER calibrations are supplied and approved. Other jet trees are reader-only.
- Full production equivalence was established for the stated p-going embedding
  sample/configuration, not every physics sample or systematic configuration.
  Other configurations have synthetic coverage and the stated Pb-going fixture.
- Production individual-track/particle trees were unavailable; those paths were
  checked with synthetic schemas matching the existing `main` reader.
- Mixed rho layouts now fail explicitly instead of inheriting ambiguous legacy
  friend-chain behavior. All schemas are validated per physical input file.
- Equal tree lengths do not prove consistent event ordering; input producers must
  preserve that forest invariant.
- No cluster jobs or notebook workflows were executed. Default ROOT compatibility
  is supported by exact output comparison, not by rerunning every downstream
  notebook. Linux/cluster compilation remains to be checked on that platform.
- The existing pPb invocation-dependent normalization and signed vertex diagnostic
  conventions are preserved, not repaired by this architectural refactoring.

There are no unexplained numerical discrepancies in the comparisons reported here.

## File inventory and version-control handoff

Modified tracked files: `CMakeLists.txt`, root `README.md`, `pPb8160_runner.py`,
`pp_analyzeMc_cp5tune.py`, `processing/submit_process_forest_condor.py`, and
`processing/run_process_forest_condor.sh`.

Added under `framework/`:

| Group | Files |
|---|---|
| Public headers (`include/jetanalysis/`) | `Config.h`, `Event.h`, `EventReader.h`, `CMSForestReader.h`, `DiJetAnalysis.h`, `RootBranchContract.h`, `Output.h` |
| Implementation (`src/`) | `main.cpp`, `CMSForestReader.cpp`, `DiJetAnalysis.cpp`, `Histograms.h`, `Histograms.cpp`, `Constants.h`, `Output.cpp`, `StageTimings.h` |
| C++ tools (`tools/`) | `check_forest_schema.cpp`, `validate_reader.cpp`, `compare_root_outputs.cpp` |
| Tests (`tests/`) | `test_root_branch_contract.cpp`, `test_workflow.py` |
| Workflow tools (`python/`) | `make_validation_fixture.py`, `validate_and_benchmark.py`, `profile_analysis.py` |
| Documentation | `README.md`, `VALIDATION.md`, `PERFORMANCE.md` |
| Recorded metrics (`reports/`) | `embedding_benchmark.json`, `final_resources.json`, `profile_summary.json` |

No files were removed. Changes remain uncommitted. No commits, merges, rebases,
resets, or pushes were performed. `main` and the original `macro` branch pointers
are unchanged. The original macro, correction helpers, calibration files, and
LinkDef are unchanged. Existing unrelated untracked files remain untouched.
Generated builds, ROOT files, logs, and the dry-run campaign are not source changes.

## Shared-library/task revision — 2026-10-20

Reviewed `main`'s CMake shared-library target, `BaseAnalysis.h`, `Manager.cc`,
and `dijetAna_pPb8160.cxx`. The current framework now uses that separation with
ordinary C++ objects: shared core (`JetAnalysis::Framework`), a separate shared
DiJet module (`JetAnalysis::DiJet`), and task applications linked to those libraries.
`AnalysisManager` owns tasks and dispatches one immutable raw event to all of them.
The complete preselection entry count is passed in `RunContext`. DiJet retains
its own selection, correction and sample-split state and output ownership.

The shared-library executable repeated the full real-embedding comparison:
502,400 events, 490,102 selected, 31,451,653 raw reco jets, 7,686,551 gen jets,
and exact agreement of all 997 ROOT objects / 18,457,544 cells or sparse bins.
Artifacts: `/tmp/jetAnalysis-library-full-validation/`.

Task tests verify callback ordering, exact 64-bit event IDs, multiple tasks,
single-use ownership, and failure propagation without normal finish callbacks.
The integration workflow runs two independently configured DiJet task instances
and an observer through one reader: both DiJet outputs match the reference and
the observer sees every raw event, including those rejected by DiJet. It checks
that DiJet's private vertex transformations do not affect the shared input.

An installed-consumer test installs into a temporary prefix, builds a separate
CMake project using `find_package(JetAnalysisFramework CONFIG REQUIRED)`, and runs
from outside the checkout. Its track-inventory example links the core only;
its DiJet example links the DiJet module and passes exact ROOT comparison.
The test also checks the installed executable, runtime library lookup, and the
correction PCM. ROOT's `TreePlayer` dependency is now explicit; the initial link
failure was reproduced in Debug before correcting it. The initially missing
installed PCM was fixed and the clean installed consumer rerun passed.

The normal Release and mixed-prior Release builds passed all seven CTest tests.
The Debug build passed all seven tests covering the same task and integration paths. Current build
artifacts are in `build/`, `/tmp/jetAnalysis-library-release/`,
`/tmp/jetAnalysis-library-debug/`, and `/tmp/jetAnalysis-library-mixed/`.
No femtoscopy physics was added; that module's scientific definitions remain
outside this revision. No branch switches, commits, merges, or pushes occurred.

Additional files in this revision:

- `include/jetanalysis/AnalysisTask.h`, `AnalysisManager.h`;
- `src/AnalysisManager.cpp`;
- `cmake/JetAnalysisFrameworkConfig.cmake.in`;
- `examples/CMakeLists.txt`, `examples/dijet_task.cpp`, `examples/track_inventory.cpp`;
- `tests/test_analysis_manager.cpp`, `tests/test_dijet_composition.cpp`,
  `tests/test_library_consumer.py`;
- `reports/library_task_resources.json`.

The existing DiJet implementation, public header, driver, CMake, workflow test,
and documentation were updated for the task lifecycle. The original macro and
calibration helpers remain unchanged. The reader's ROOT type checks are unchanged.
