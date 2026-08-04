# Step 05: MC unfolding and closure

**Status:** Proposed; a RooUnfoldBayes notebook prototype with explicit
miss/fake handling exists, but no independent closure result is approved.

## Goal and problem

Select an unfolding method and regularization that recover independent embedded-MC
truth without unacceptable bias or variance in the full two-dimensional space.

The complete closure suite is mandatory separately for Pb-going MC and p-going
MC. The corresponding data direction remains blocked unless its own MC direction
passes; a successful combined-MC closure is insufficient.

## Recommended solution

1. Build a reusable, non-interactive unfolding driver from the validated response.
2. Scan Bayesian iterations as the first candidate and compare with at least one
   meaningful alternative or limiting case, such as SVD or bin-by-bin only as a
   diagnostic.
3. Perform same-generator statistically independent closure, reweighted-shape
   stress tests, and direction-specific closures.
4. Propagate the full covariance, including weighted-MC statistical uncertainty;
   use toys or bootstrap where RooUnfold's analytic treatment is insufficient.
5. Select regularization with a rule defined before data unfolding.
6. Unflatten to physical 2D histograms, then project/merge into final pT intervals
   and self-normalize each eta distribution with covariance propagation.

## Validation and exit criteria

- Closure residuals and pulls satisfy approved thresholds in all reported bins.
- Pull means are compatible with zero and widths with one within MC precision.
- Results are stable under reasonable regularization and prior variations.
- Covariance matrices are finite, symmetric, and positive semidefinite within
  numerical tolerance; normalization-induced correlations are retained.
- Both directions close independently before any combined closure is attempted.
- The analyst approves the algorithm and regularization rule.

## Required figures

- Truth/measured/unfolded spectra, ratio and pull panels, iteration scan, covariance
  and correlation matrices, refolding test, and stress-test closures; PDF+PNG.

## Results

Not started. `hist_analysis/notebooks/05_unfold2D.ipynb` currently loads a local
RooUnfold build and defaults to a five-iteration Bayesian prototype. It includes
response inefficiency for misses, explicitly enables fake handling, rescales
small weighted response entries before RooUnfold matrix sanitization, and writes
per-pTave projection and closure diagnostics. The current configuration uses the
same sample for response construction and closure; no statistically independent
closure criterion, regularization-selection rule, or approved covariance study
is documented.
