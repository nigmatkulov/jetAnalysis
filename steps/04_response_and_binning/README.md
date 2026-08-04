# Step 04: Two-dimensional response and analysis binning

**Status:** Proposed; 4D sparse response storage is implemented but unvalidated.

## Goal and problem

Construct a statistically stable detector response mapping truth to reconstructed
\((p_T^{\mathrm{ave}},\eta_{\mathrm{dijet}})\), including cross-bin migration,
misses, and fakes. Choose fine internal bins and a defensible rule for merging
them into final pT ranges.

Separate nominal responses are constructed and validated for Pb-going embedded
MC and p-going embedded MC. These responses correct the corresponding Pb-going
and p-going data in Step 07. A combined response may be studied only as a closure
cross-check and is not the default correction.

## Recommended solution

1. Freeze truth and reconstructed phase spaces, including underflow/overflow.
2. Define event-level dijet matching and distinguish matched, missed, fake, and
   swapped-leading/subleading cases.
3. Use statistically independent training and testing samples, preferably by a
   deterministic event split if independent productions are unavailable.
4. Validate every THnSparse axis order and projection against direct event-loop
   histograms using integrals and selected cells.
5. Flatten the 2D bins only as a RooUnfold interface; retain an explicit reversible
   mapping between global and physical bins.
6. Choose fine bins using resolution and occupancy. Merge final pT intervals using
   documented thresholds on response statistics, stability, purity, and closure—not
   by visual preference after examining data.

## Validation and exit criteria

- Response bookkeeping satisfies weighted integral identities for matched,
  missed, fake, truth, and measured populations.
- Off-diagonal pT and eta migration blocks are present and numerically verified.
- All retained bins satisfy approved occupancy/stability/purity requirements, or
  are merged/masked by a recorded algorithm.
- Projection and flatten/unflatten round-trip tests pass exactly within floating
  precision.
- Training and closure samples are statistically independent.
- Pb-going and p-going response bookkeeping and closure tests pass independently.

## Required figures

- Full flattened response, representative physical response blocks, efficiency,
  fake rate, purity, stability, condition diagnostics, and bin-merging rationale;
  each PDF+PNG.

## Unresolved decisions

- Exact fine eta and pT bin edges.
- Quantitative bin-merging thresholds.
- Generator-jet versus reference-match nominal response, pending Step 00.

## Results

Not started. `hist_analysis/notebooks/05_unfold2D.ipynb` demonstrates pT-block
flattening and includes off-diagonal migration blocks, but uses provisional
`PT_AVE_BINS = [0, 40, 80, 180, 250, 300, 500, 1000]` and is not an approved
binning or response implementation. The closure-plot notebook independently
uses an edge array plus adjacent-edge iteration, so its diagnostic intervals do
not freeze the final response or reporting bins.
