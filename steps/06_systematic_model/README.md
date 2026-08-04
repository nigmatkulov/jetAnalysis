# Step 06: Systematic-uncertainty model

**Status:** Proposed; several detector variations are implemented but unvalidated.

## Goal and problem

Define, propagate, and combine uncertainties for the self-normalized unfolded
distribution without double counting detector, model, and statistical effects.

Each applicable variation is produced separately for Pb-going MC, p-going MC,
Pb-going data, and p-going data. A source may be MC-only or data-only, but that
scope must be stated explicitly. Correlation assumptions are applied only during
the later direction combination, after the four relevant channel results exist.

## Candidate sources from the current design

- JER nominal/up/down;
- eta-dependent JER model variation;
- CM-boost variation around 0.465;
- jet pointing/angular resolution, after definition;
- underlying-event/model dependence from Pythia versus Pythia embedded in EPOS;
- nominal JEC/JEU components where applicable;
- unfolding regularization and prior/model dependence;
- finite response statistics;
- jet ID and selection stability;
- data pileup as a control study unless evidence requires a systematic assignment.

## Recommended solution

For each source, specify its physical meaning, variation, affected response and/or
measured distribution, propagation method, symmetrization rule, and correlations
between directions and bins. Rebuild the response when the variation changes
migration or efficiency; varying only the final histogram is insufficient.

Self-normalize every varied result using the same definition as nominal and
propagate normalization correlations. Compare signed shifts before reducing them
to an uncertainty envelope or covariance.

## Validation and exit criteria

- Every included source has evidence, a reproducible implementation, and no known
  overlap with another source.
- Direction and bin correlations are explicitly documented.
- Availability and validation status are tabulated for all four channels for every
  applicable source.
- MC statistical fluctuations are separated from systematic shifts.
- Pileup is either excluded with quantitative evidence or included with a defined
  variation.
- Total covariance is numerically valid and dominant sources are understood.

## Required figures

- Per-source signed shifts, fractional uncertainties, correlation matrices,
  direction comparisons, and total uncertainty breakdown; PDF+PNG.

## Results

Not started. JER and boost variants already exist in processing/plotting code,
but their exact interpretation and propagation have not passed this design.
