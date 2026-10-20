# Modular C++ analysis framework

The `jetAnalysis` executable runs the migrated dijet analysis. The untouched
`processForestSimple.C` and its compiled `processForestSimple` executable remain
available as independent physics references. All migration work is on
`macro2framework`, forked from `macro` at
`d934635819cc6bbcfa29eb830bcfb70d63d28c21`.

## Build

Use the existing ROOT installation and correction payloads; no additional C++
libraries are required. C++17 or the newer standard required by ROOT is used.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

For debugging use a separate directory and `-DCMAKE_BUILD_TYPE=Debug`. For
profiling use `-DCMAKE_BUILD_TYPE=RelWithDebInfo`. Optional mixed-prior output
retains the existing build switch `-DDO_MIXED_PRIOR=ON` for both executables.
The correction classes and their ROOT dictionary are built once in
`jetanalysis_corrections`, shared by the reference and framework.

CTest runs task lifecycle/dispatch and installed-consumer checks, branch-contract checks, argument checks, existing analysis-contract
check, and synthetic end-to-end physics comparisons. The latter uses Python and
PyROOT; `JETANALYSIS_PYTHON` can select the interpreter when configuring CMake.
This checkout defaults to its `py-env/bin/python` when available, using the
existing ROOT loader. A runtime ROOT/PyROOT failure is a failed test, not a skip.

## Shared libraries and analysis tasks

The shared-library organization follows `main`'s `BaseAnalysis` / `Manager`
pattern, using ordinary C++ interfaces and explicit ownership:

| CMake target | Shared library | Responsibility |
|---|---|---|
| `JetAnalysis::Framework` | `libjetanalysis_framework` | Event model, typed CMS reader, task interface, shared event loop |
| `JetAnalysis::DiJet` | `libjetanalysis_dijet` | DiJetAnalysis, its corrections, histogram definitions and output |

The filename suffix is `.dylib` on macOS and `.so` on Linux. The `jetAnalysis`
executable is an application linked to the DiJet library, which in turn links the
core. A track-based task can link only the core; it has no dependency on DiJet.

`AnalysisTask` has three callbacks: `initialize(const RunContext&)`,
`process(const Event&)`, and `finish()`. `RunContext::entries` is the complete
invocation count before any task cuts, preserving the pPb normalization convention.
Each task owns its selections, corrections, derived collections, output, and any
mixing history. The manager borrows its reader and exclusively owns tasks through
`std::unique_ptr`. All tasks receive the same immutable raw event in registration
order; DiJet's private vertex transformations, sorted jets, and rejection do not
modify the input or stop other tasks from receiving it. Borrowed event data must
not be retained past the next reader call; a mixing task must own the selected
particle values it retains. Processing is sequential.

```cpp
#include <jetanalysis/AnalysisManager.h>
#include <jetanalysis/CMSForestReader.h>
#include <jetanalysis/DiJetAnalysis.h>

jetanalysis::AnalysisConfig config;
// Set input, output, sample, beam direction, pT-hat, and other task options.
jetanalysis::CMSForestReader reader(config.reader);
jetanalysis::AnalysisManager manager(reader);
manager.addTask(std::make_unique<jetanalysis::DiJetAnalysis>(config));
// manager.addTask(std::make_unique<MyTrackTask>(trackConfig));
manager.run();
```

The DiJet task and manager are single-use. On reader/process exceptions the manager
propagates the error and does not call normal finish/output callbacks. Finish
callbacks run sequentially; outputs already written by an earlier finish callback
are not rolled back if a later one fails. Each task must use its own output path.
The existing `DiJetAnalysis(config).run(reader)` remains a single-task convenience
using this same manager lifecycle.

Install the library, headers, ROOT correction PCM, executable, and CMake package:

```sh
cmake --install build --prefix /tmp/jetanalysis-install
cmake -S framework/examples -B /tmp/jetanalysis-consumer \
  -DCMAKE_PREFIX_PATH=/tmp/jetanalysis-install
cmake --build /tmp/jetanalysis-consumer -j
```

A separate application's CMake file can use:

```cmake
find_package(JetAnalysisFramework CONFIG REQUIRED)
add_executable(my_task my_task.cpp)
target_link_libraries(my_task PRIVATE JetAnalysis::Framework)
# For a DiJet application, link JetAnalysis::DiJet instead; core is transitive.
```

The same target names work with `add_subdirectory()` on the source checkout.
The installed package requires the compatible ROOT installation used to build it.
Calibration data are not bundled into the library installation; set the task's
`auxDirectory` (or CLI `--aux-directory`) when running outside the checkout.

`examples/track_inventory.cpp` is a compiled core-only task example and counts
tracks without adding physics definitions. `examples/dijet_task.cpp` is a compiled
DiJet application, accepting input, output and an optional auxiliary-data directory.
Its concrete preset is p-going embedding with pT-hat 540; adjust its task
configuration for other samples.
A future `TwoParticleFemtoscopy` module would implement this same interface and
link the core. Femtoscopy pair cuts, particle identification, mixing pools, weights,
and correlation normalization have not been specified or implemented here.

## Run

The existing positional arguments and their defaults are retained:

```sh
build/jetAnalysis input.root output.root 1 0 540 0 2 0
build/jetAnalysis inputs.list output.root 1 0 540 0 2 0 pPb 0 0
build/processForestSimple inputs.list reference.root 1 0 540 0 2 0 pPb 0 0
build/compare_root_outputs reference.root output.root
```

Arguments are input, output, MC type, beam direction, pT-hat sample, trigger,
reco jet selection, vertex filter, optional collision system, JER mode, and
eta-dependent JER switch. Their meaning is unchanged; see the root README.
Run `build/jetAnalysis --help` for extra options. Empty input/output arguments
retain the reference's default path construction. Run from the repository root,
or set `--aux-directory` to the correction directory. `--jec FILE` may be
repeated to override the JEC chain; `--jeu FILE` overrides the uncertainty file.
Default correction paths and selection by sample/direction remain unchanged.

Output ROOT names, classes, titles, axes, statistics, flow bins, uncertainty
storage, response matrices, and systematic families are unchanged by default.
Metadata remains mandatory. Existing output groups can be suppressed by their
ROOT key prefix, without changing selection, weights, or histogram filling:

```sh
build/jetAnalysis input.root output.root 1 0 540 0 2 0 \
  --exclude-output hRef --exclude-output hGen
```

Exclusion controls writing only; it does not promise lower event-processing time
or histogram memory usage. Enabling/disabling output does not affect physics.

## Collections and readers

`CMSForestReader` implements `EventReader`, with exclusive ownership of each
physical ROOT file, tree readers, and typed branch handles. It preflights every
accepted file, sums the full input entry count, and visits files in list order.
It checks tree alignment per file. Unreadable/incomplete forests are skipped with
explicit diagnostics as in the reference; a missing explicitly requested optional
collection or a malformed schema fails the invocation. No partial analysis output
is written after such a failure. ROOT file-list comments and entry validation are
retained. Configured trees must have the same entry ordering; alignment checks
cannot prove that two equal-length trees describe identical event IDs.

Default collections are reco jets and MC gen jets from `ak4PFJetAnalyzer/t`.
Individual tracks and particles are off by default. Their absence is represented
by availability flags, distinct from an available but empty collection.

```sh
# Inspect an available alternative jet collection without running dijet physics.
build/jetAnalysis input.root unused.root 1 0 540 0 2 0 \
  --reader-only --jet-collection AK4CSPF

# Read tracks/particles without jets, when the corresponding trees exist.
build/jetAnalysis input.root unused.root 1 0 540 0 2 0 \
  --reader-only --no-reco-jets --no-gen-jets --read-tracks --read-particles
```

`AK4PF`, `AK4CSPF`, and `AK5PF` map to `ak4PFJetAnalyzer/t`,
`akCs4PFJetAnalyzer/t`, and `ak5PFJetAnalyzer/t`. Other PF radii use the same
naming convention. `--jet-tree PATH` selects an explicit path. Track and particle
tree paths can be changed with `--track-tree` and `--particle-tree`.

**Approved scope:** DiJetAnalysis is restricted to AK4PF until collection-specific
JEC/JER calibrations are provided and approved. Alternative trees work only in
reader-only mode. DiJetAnalysis requires reco jets and, for MC, gen jets.
Enabling tracks/particles does not add track cuts, weights, or histograms.

The initial track schema follows the existing `main` reader: `ppTrack/trackTree`
with `nTrk`, `trkPt/Eta/Phi`, `trkCharge`, and `highPurity`. Generator particles
use `HiGenParticleAna/hi` vectors `pt/eta/phi/chg/pdg/sube`. These are kinematics,
identity, and quality fields; this migration does not claim a complete track
analysis data model. Formats with other branch types are rejected explicitly.
Production embedding skims used here contain no individual-track/particle trees;
these optional paths were validated with synthetic forests, not production data.

## Input types and memory

Transient `Event`, `RecoJet`, `GenJet`, `Track`, and `GenParticle` are ordinary
C++ value objects; no event dictionary or ROOT inheritance is required. Event
IDs use fixed-width integers, and physics floating-point precision is preserved.
Input handles use exact ROOT types (`Float_t`, `Int_t`, `UInt_t`, `ULong64_t`,
`Bool_t`) and typed `TTreeReaderValue`/`TTreeReaderArray` access. The contract checks
leaf types, signedness, scalar/array shape, count-leaf identity, and vector element
types before reading. Required representations are checked at compile time.
No branch buffer is reinterpreted or bound through a differently typed pointer.
Jet counts retain the original 1500-jet limit; dynamic ROOT arrays avoid a fixed
buffer overflow even when rejecting an oversized event. Optional collection
sizes are checked across their component branches.

ROOT aliases are not universal fixed-width guarantees: `Long_t` follows native
`long`, and `Float16_t`/`Double32_t` have special storage semantics. Neither is
silently accepted in place of the schema's declared types.

The reference's rho-absent fallback remains zero with a warning. Mixed rho branch
names/presence within an input list are rejected: the legacy friend chain can
retain stale rho values at such transitions. The new reader does not silently
change those smearing seeds. Homogeneous absence or either supported homogeneous
rho branch works. This is an explicit stricter input check.

## Architecture and physics ordering

- `include/jetanalysis/Event.h`: ROOT-independent transient objects and availability.
- `Config.h`, `EventReader.h`, `CMSForestReader.h`: configuration and input interface.
- `src/CMSForestReader.cpp`: per-file schema checks, typed access, conversion.
- `AnalysisTask.h`, `AnalysisManager.h`, `src/AnalysisManager.cpp`: task callbacks, ownership, and the shared event loop.
- `src/DiJetAnalysis.cpp`: per-task physics state, corrections, selections, matching,
  nominal and systematic filling. No mutable physics state is shared between runs.
- `src/Histograms.*`: original histogram definitions and ownership.
- `src/Output.cpp`: original output ordering and optional key-prefix exclusions.
- `src/StageTimings.h`: opt-in stage measurements (`--profile-stages`).

The architecture in `main` informed the reader/event/analysis split. Its
ROOT-derived event hierarchy and pointer-owned object collections were not copied.
The physics implementation remains deliberately close to the macro; separating
it does not imply redesigning its formulas or calibration model.

The manager reads each event once, then invokes tasks. DiJet takes its
controlled/studied random draw once per input event; its sequence is unchanged.
The remaining physics order is vertex orientation;
raw vertex fill; event cuts; MC weight; reco jets with JEC; gen jets and overweight
rejection; selected-event diagnostics; gen jets; gen dijets; reference dijets;
reco jets with configured ID and JEU; reco dijets, JER variants, and unfolding.
The original gen-before-reco ordering, deterministic smear seeds, reference
associations, sorting, misses/fakes, and sample split are retained.

## Validation and performance tools

```sh
build/check_forest_schema input.root 1 pPb
build/validate_forest_reader input.root 1 pPb
build/compare_root_outputs reference.root candidate.root
py-env/bin/python framework/python/validate_and_benchmark.py \
  --build build --input input.root --output-dir /tmp/jetanalysis-benchmark \
  --repeat 2 --analysis-args 1 0 540 0 2 0
py-env/bin/python framework/python/profile_analysis.py \
  --executable build/jetAnalysis --input input.root --output-dir /tmp/jetanalysis-profile \
  --analysis-args 1 0 540 0 2 0
```

Add `--sample` to the profiling command on macOS to collect native CPU stacks.
macOS may require permission to inspect the child process. Stage timings work
without that profiler. Benchmarking uses per-child `wait4` resources, records
elapsed/user/system time, peak RSS with platform unit conversion, event counts,
throughput, logs, and exact numerical comparisons in JSON. It alternates run
order on successive repetitions. It does not control the host load or caches.

`make_validation_fixture.py INPUT OUTPUT --entries 2000` creates a bounded,
aligned forest for development. It refuses an existing output. A pPb fixture's
normalization differs from the full source: compare only reference/framework
runs on identical fixtures, not fixture yields against full-sample yields.
The existing pPb per-invocation normalization remains unchanged and still makes
arbitrary job splitting unsuitable for direct yield summation without an approved
shared normalization convention.

`validate_forest_reader` includes the untouched macro as an oracle and compares
every consumed input field and reco/gen object. `compare_root_outputs` recursively
checks key names/cycles, concrete types, titles, axes/labels, stored errors,
statistics, entries, all flow and regular cells, and sparse coordinates/content/
errors. It uses exact equality, refuses unsupported object classes, and excludes
file serialization metadata such as timestamps and compression buffers.

See [VALIDATION.md](VALIDATION.md) and [PERFORMANCE.md](PERFORMANCE.md) for results.

## Submission

The pPb and pp runners and `processing/submit_process_forest_condor.py` default to
`build/jetAnalysis`. `--executable` selects the reference or another build.
Repeat `--analysis-arg` to forward framework options without shell evaluation:

```sh
py-env/bin/python pp_analyzeMc_cp5tune.py --dry-run \
  --analysis-arg=--exclude-output --analysis-arg=hRef
```

The same forwarding is available in the pPb runners and Condor submitter. The
worker retains the original 13 positional fields and accepts trailing framework
arguments. Existing splitting, output labels, sequential defaults, and bounded
failure logging remain unchanged. No jobs were submitted during migration.
