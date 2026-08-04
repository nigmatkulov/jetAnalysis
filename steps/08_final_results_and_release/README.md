# Step 08: Final results and reproducible release

**Status:** Proposed.

## Goal and problem

Produce the final self-normalized distributions, complete uncertainty record, and
a reproducible GitHub release suitable for internal review and later publication.

The release must retain the four-channel provenance chain: Pb-going MC and
p-going MC validation products, separate corrected Pb-going and p-going data
results, and the derived combined data result.

## Recommended solution

1. Freeze final pT intervals and eta reporting bins from the approved Step 04 rule.
2. Produce final plots and numerical tables with statistical and systematic
   covariance matrices.
3. Verify normalization, units, labels, bin-width treatment, and common CM sign.
4. Document inputs, environment, configuration, commands, random seeds, and Git
   commit; tag the reviewed release on GitHub.
5. Archive lightweight configuration and final tables in Git; record stable
   locations/checksums for large ROOT products.
6. Conduct an independent reproduction from a clean checkout before release.

## Validation and exit criteria

- Each reported pT interval satisfies its documented normalization exactly:
  either the per-bin fractions sum to one (`integral`) or the density has a
  bin-width-weighted integral of one (`bin_width`). Final results must use one
  convention consistently and label the y axis with the corresponding units.
- Tables reproduce plotted points and uncertainties exactly.
- Statistical, systematic, and total covariance matrices accompany the result.
- Every main figure exists as both PDF and PNG.
- Separate direction results and their covariance inputs reproduce the released
  combined result using the documented combination procedure.
- Clean-checkout reproduction succeeds and the review checklist is signed off.

## Required figures

- Final distributions for every pT interval, uncertainty breakdowns, direction
  comparison, and key validation summary plots; PDF+PNG.

## Results

Not started.
