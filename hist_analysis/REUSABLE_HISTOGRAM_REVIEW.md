# Reusable histogram cut and naming review

Reviewed October 2, 2026 against the final user-updated C++ cut table:
`(1.5, 1.6, 1.7, 1.8, 1.9, 2.3, 2.4, 3.0)`.

## Findings and corrections

- Indices 0–7 now map to those eight physical per-jet acceptances. Persisted
  projections use `etaCut1p5` through `etaCut3p0`. The shared Python
  `DIJET_ETA_CUTS` table is checked against C++ by a regression test.
- All 19 notebook-local cut tables reference the shared table. Defaults that
  selected 1.9 by index 5 now use `STANDARD_DIJET_ETA_CUT_INDEX`, derived from
  the physical value (index 4). Trigger stitching uses the same index. The
  pointing-resolution Lab default preserves 2.3 by deriving index 5 from its
  value. This avoids changing the intended physical acceptance when indices shift.
- Cuts apply upstream to both jets using the source family's frame and jet
  collection. Production must not reapply them to projected dijet eta, since
  reference jets associated with selected reco jets can have different eta.
  Projection only adds the configured half-open pTave interval. Nominal C++
  selection uses leading/subleading pT >= 50/40 GeV, |delta phi| >= 2pi/3,
  and per-jet |eta| <= the cut. Reusable intervals start at 50 GeV, excluding
  the otherwise allowed pTave interval 45–50 GeV. Flow bins are excluded.
- The broad projection-key pattern mistook run numbers such as `_285480` for
  cut indices. Run diagnostics are now copied unchanged without generating
  misleading eta-cut projection keys. Sparse responses remain copied objects.
- Producer filenames used single underscores while some readers used doubles.
  Producers and readers now share canonical construction. Legacy filenames and
  projection keys remain readable where the requested physical selection is
  available and consistent. Metadata remaps archived indices by physical cut;
  a missing requested cut raises an error rather than substituting another cut.
- Schema 2 metadata records the current cut table. Copied upstream ROOT keys
  retain their indices. Batch production preflights all sources before writing;
  conflicting source titles, unknown indices and source/output identity are
  rejected. Files are closed reliably if opening an output fails.

## Validation and source migration

The audit covered 6,516 eta TH2s across 21 source files; all configured pTave
boundaries aligned with source-bin edges. Existing data files use an older
8-cut table starting at 1.4 without 2.3, while existing official/standalone MC
files use the earlier 9-cut table starting at 1.4. All must be rebuilt for the
final 8-cut table starting at 1.5. Some titles omit a cut; titles are a
consistency check, not proof of the event selection.

The standalone CP5 and CUEP8M1 analyzer macros were inspected read-only and
already contain the new eight-cut table. Rerun histogram analysis on the stored
generator inputs; new PYTHIA event generation is not needed for this change.
After rebuilding official/data and standalone source ROOT outputs, run notebook
00 to create the new physical-name bundles, then rerun downstream notebooks.

Python tests cover every physical label, source-order mismatch rejection, both
consumer paths, archived metadata remapping, weighted interval boundaries,
reference-eta preservation, run-ID handling, filenames, and file overrides.
All shared helpers passed compilation, all 31 notebooks passed JSON/Python
syntax checks, and the complete Python unit suite passed (44 tests). The Release executable
build and CTest passed; the compiler reported two unused-parameter warnings. Representative
production was exercised under `/tmp` while source and shared cut tables still
matched the earlier nine-cut order. Final source preflight correctly rejects
the existing outputs after the user's final eight-cut update. The final eight-cut
production is exercised with synthetic ROOT inputs; full production/notebook
execution awaits rebuilt source files. Persistent ROOT outputs were not changed.

Only the requested C++ cut-table update is included in this commit; other
C++ physics edits and unrelated working-tree changes remain uncommitted. The commit contains reusable-production integration,
shared cut selections, reader compatibility, tests and focused documentation.

## Possible improvements

Store the actual cut table and selection provenance in upstream ROOT metadata
so older cut orders can be verified without relying on titles. Optionally check
source freshness before selecting a reusable bundle, and publish completed
bundles through atomic replacement to prevent partial files after an I/O failure.
