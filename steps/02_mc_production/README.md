# Step 02: Existing embedded-MC inventory and bookkeeping

**Status:** Proposed; all required products exist, but their provenance and
validation are not yet documented.

## Goal and problem

Establish that the existing embedded Pythia+EPOS products for p-going and
Pb-going directions are complete, statistically correct, and usable, preserving
direction-dependent corrections and enough metadata to reproduce their merge
across pt-hat samples. Regenerating the existing production is not part of this
step.

This is an existing-product validation step: Pb-going MC and p-going MC are
mandatory separate inputs. Existing Pb-going and p-going experimental-data
products are inventoried before Step 07 but are not corrected or combined here.

## Recommended solution

1. Inventory every pt-hat input with production tag, generator settings, event
   count, direction/orientation, and checksum or immutable dataset identifier.
2. Validate pt-hat intervals, cross-section/event weights, vertex weights,
   overweight rejection, and boundary handling.
3. Locate the existing direction-specific outputs and recover or reconstruct
   their exact configurations and immutable logs where available.
4. Validate the existing pt-hat merge by comparing overlaps and normalization.
5. Keep direction outputs separate through closure and systematic evaluation.
6. Create any combined-MC diagnostic only after both existing directions pass
   their individual checks, and label it as a derived combined result.

The current Python runners document how the products were or can be generated.
They should be changed only if required for provenance or a separately approved
targeted rerun. A full rerun is not an exit criterion.

## Validation and exit criteria

- All expected existing samples and products are present or explicitly excluded
  with justification.
- Weighted pt-hat spectra are continuous within a predefined statistical metric.
- Existing merge contents are reproducible from per-sample products or validated
  against recorded production metadata.
- No missing, duplicate, corrupt, or silently skipped files remain.
- Direction-specific JEC and orientation metadata are recorded in output metadata.

## Required figures

- Weighted and unweighted pt-hat spectra, sample-boundary ratios, vertex-z
  closure, and overweight diagnostics, each PDF+PNG.

## Unresolved inputs

- Dataset identifiers and the cross sections/effective event counts used by the
  current hard-coded weighting function must be confirmed.
- Training/test split strategy is deferred to Step 04 but must be compatible with
  available statistics.

## Results

Not started. All data and MC productions are reported to exist. Runners cover
pt-hat values 15 through 540 for both directions and document embedded-MC, JER,
and nominal jet-ID settings; no full production rerun is planned.
