# Step 07: Data correction and direction combination

**Status:** Proposed; blocked until Step 05 is complete.

## Goal and problem

Apply frozen direction-specific corrections to p-going and Pb-going data, verify
compatibility, and combine them into the common physics orientation without
masking an orientation or calibration error.

The required execution order is explicit:

1. inventory and validate the existing Pb-going data product independently;
2. inventory and validate the existing p-going data product independently;
3. correct Pb-going data with the validated Pb-going MC response;
4. correct p-going data with the validated p-going MC response;
5. transform both corrected results to the common convention;
6. test compatibility and only then combine them.

## Recommended solution

1. Freeze all selections, binning, response, regularization, and systematic rules
   before examining unfolded distributions beyond necessary technical checks.
2. Validate trigger efficiency/combination and data-quality selections in the
   existing products separately by direction. Cross-section normalization is not
   required, but event and trigger weighting can still affect the self-normalized
   shape.
3. Unfold each direction independently with its corresponding corrections.
4. Refold and compare with measured distributions.
5. Test p-going/Pb-going compatibility using the full covariance after transforming
   both to the agreed common sign convention.
6. Combine using an explicit covariance-aware estimator; document which systematic
   sources are correlated between directions.

## Validation and exit criteria

- Direction-specific refolding tests pass.
- Separate Pb-going and p-going data products, correction outputs, and covariances
  exist before the combined product is created.
- Control distributions reveal no unresolved pileup, trigger, or orientation bias.
- Compatibility statistic and degrees of freedom are documented.
- The combined estimator and covariance pass toy tests.
- No binning or method choice is tuned to improve the observed final data shape.

## Required figures

- Separate direction results, refolding comparisons, direction ratio/difference,
  compatibility pulls, and combined self-normalized spectra; PDF+PNG.

## Results

Not started; deliberately gated by MC closure approval. The required Pb-going and
p-going data productions already exist and will be consumed rather than regenerated.
