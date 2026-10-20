# Performance measurements — 2026-10-20

Physics equivalence was established before profiling. No physics or calibration
optimization was applied during this migration.

## Full embedding benchmark

Input: complete p-going embedding pT-hat 540 forest; 502,400 events, 490,102
selected. Same Release build environment and analysis arguments `1 0 540 0 2 0`.
The benchmark runs sequentially, alternating reference/framework order on the
second repetition. Every pair passes the exact ROOT-output comparison.

| Run | Executable | Elapsed (s) | User CPU (s) | System CPU (s) | Events/s | Peak RSS (MiB) |
|---|---|---:|---:|---:|---:|---:|
| 1 | Reference | 174.320 | 172.652 | 0.824 | 2,882.1 | 932.44 |
| 1 | Framework | 175.412 | 173.881 | 0.785 | 2,864.1 | 962.08 |
| 2 | Framework | 175.209 | 174.086 | 0.776 | 2,867.4 | 964.80 |
| 2 | Reference | 175.123 | 174.022 | 0.794 | 2,868.8 | 922.16 |

Mean elapsed: reference 174.721 s, framework 175.311 s, a measured difference
of +0.34%. Mean peak RSS: reference 927.30 MiB, framework 963.44 MiB (+3.90%).
This shows comparable throughput on this workload, not a statistically established
speed change. Two repetitions do not characterize all production workloads.

Timings include process startup, input validation, event processing, and output
writing. Peak RSS is obtained per child with `wait4`, with macOS byte units
converted explicitly to MiB. OS caches and host load were not controlled; some
development checks ran concurrently. These measurements must not be interpreted
as a controlled optimization result. Raw commands, per-child measurements,
counts, and comparison results are retained in
[reports/embedding_benchmark.json](reports/embedding_benchmark.json).

This two-pair baseline was built before the final schema/CLI guards and opt-in
stage instrumentation. Those refinements do not change supported-input physics;
a final-build full-input confirmation is recorded separately below and in the
validation report.

The final analysis build was then rerun on the complete same input: 174.605 s
elapsed, 173.272 s user CPU, 0.797 s system CPU, 2,877.4 events/s, and 954.09 MiB
peak RSS. Its output again matched all 997 objects and 18,457,544 cells/sparse bins
exactly. Resources are in [reports/final_resources.json](reports/final_resources.json).
This confirmation is not an additional independently controlled reference pair.

## Stage timings and native stack profile

The RelWithDebInfo framework processed a real Pb-going embedding fixture: the
same leading 2,000 aligned entries were repeated ten times (20,000 total,
19,480 selected). This bounds profiling time but is not an independent sample
or a full-production normalization. Its output matches the reference exactly:
997 objects and 15,586,362 cells/sparse bins. The sampled and unsampled runs both
pass comparison.

| Timed stage | Seconds |
|---|---:|
| Read and convert | 0.268 |
| Reco loading, JEC, and sorting | 3.388 |
| Gen loading and sorting | 0.001 |
| Gen jet processing | 0.073 |
| Gen dijet processing | 1.369 |
| Reference dijets | 0.010 |
| Reco jets | 0.238 |
| Reco dijets and response/systematic processing | 0.741 |
| ROOT output writing | 5.146 |

These are elapsed stage times from one instrumented run, not CPU-time fractions.
They exclude initialization, histogram creation, event cuts, weights, and some
other loop bookkeeping. On this short fixture, output serialization is comparable
to all event processing combined. Full production amortizes that fixed cost.

macOS `sample` collected native stacks for three seconds at a one-millisecond
interval during event processing, with symbols from the profiling build:

- 1,316 of 2,480 main-thread sampled stacks included `JetCorrector::GetCorrectedPT`
  (about 53%); `SingleJetCorrector::GetCorrection` payload lookup dominates there.
- The gen-dijet path accounted for roughly 24% of the sampled stacks, including
  deterministic normal-deviate generation and histogram work.
- The stage timers attribute a small fraction of this fixture's processing time
  to input conversion and vector sorting relative to calibration and physics.

The three-second window does not include final output writing and must not be
used to estimate whole-process CPU percentages. It identifies calibration lookup
as a candidate for a later isolated optimization. Any indexing, caching, random
number, or sorting change needs another exact numerical validation and benchmark;
none was introduced here.

The compact stage/sample record is retained in
[reports/profile_summary.json](reports/profile_summary.json).
Raw native profile, stage log, status and comparison are in
`/tmp/jetAnalysis-framework-profile-native/`. The initial sandbox attempt could
not inspect the process; rerunning with approved access succeeded. Stage timing
itself requires no native profiler permissions. Linux sampling was not tested.

## Reproduce

See `framework/README.md` for `validate_and_benchmark.py` and
`profile_analysis.py` commands. Use separate output directories for independent
runs. Retain the full same input list for each pair: pPb normalization depends on
the invocation's entry count. Compare outputs before interpreting any speed result.

## Shared-library/task revision confirmation

After separating the shared core, DiJet module, and generic manager callbacks,
the complete same embedding sample was processed again: 182.106 s elapsed,
180.072 s user CPU, 0.793 s system CPU, 2,758.8 events/s, and 913.47 MiB peak RSS.
All 997 objects / 18,457,544 cells or sparse bins matched the original exactly.
The measurement is retained in
[reports/library_task_resources.json](reports/library_task_resources.json).
Other builds/tests were active during this run. It is a regression confirmation,
not a controlled comparison of static versus shared libraries or callback costs.
The earlier paired benchmarks and native profile describe the preceding migration
revision; they should not be treated as fresh paired measurements of the task API.
Read timing now belongs to the shared manager, while physics/output stage timing
belongs to each analysis task.
