---
name: root-forward-backward-ratio
description: Create, modify, review, or validate forward/backward (F/B) histogram workflows in this project's histogram analysis notebooks and PyROOT helpers. Use whenever a request mentions a forward/backward ratio, F/B distribution, reflected eta ratio, positive-to-negative eta folding, or an F/B comparison. Require CERN ROOT/PyROOT and standard independent-error propagation for constructing F/B; never use ROOT's binomial division option for the F/B construction.
---

# ROOT Forward/Backward Ratio

Implement forward/backward histograms consistently with this CMS jet analysis. Treat the convention as physics-sensitive behavior and make the smallest task-scoped change.

## Non-negotiable rules

- Use CERN ROOT through PyROOT and produce a ROOT `TH1` histogram.
- Construct the ratio from unnormalized forward and backward yields. Do not normalize either input independently before division.
- Treat forward and backward as separate weighted populations, not as a subset efficiency.
- Use ROOT's standard independent-error propagation:

  ```python
  ratio = forward.Clone(name)
  ratio.SetDirectory(0)
  ratio.Divide(forward, backward, 1.0, 1.0, "")
  ```

- Never pass `"B"`, `"b"`, or another binomial option when constructing F/B. Do not expose a user configuration that permits it.
- Preserve input uncertainties. Ensure newly constructed histograms have `Sumw2()` before filling, and copy both bin contents and bin errors when folding.
- Preserve ROOT ownership and lifetimes: detach returned histograms with `SetDirectory(0)` and retain any temporary input histograms needed by PyROOT.
- Preserve existing histogram keys, axes, binning, selections, output paths, and public interfaces unless the user explicitly authorizes a change.

For one bin with values `F`, `B` and independent uncertainties `sigma_F`, `sigma_B`, ROOT's standard division corresponds to

```text
sigma(F/B)^2 = (sigma_F/B)^2 + (F*sigma_B/B^2)^2
```

This expression is a review check, not a replacement for `TH1::Divide` in the implementation.

## Select the input pattern

Inspect the target notebook, `hist_analysis/README.md`, and the relevant helper before editing. Reuse existing helpers instead of duplicating notebook code.

### Stored forward and backward histograms

Use `hist_analysis/python/dijet_closures.py` when the ROOT file already contains separate CM-forward and CM-backward TH2 histograms.

1. Resolve the configured ROOT keys without renaming them.
2. Project both TH2 inputs over the same half-open pT or pTave interval using `project_semantic_th2` or the existing closure helper.
3. Apply the same eta rebinning to both projections.
4. Check identical ROOT bin edges.
5. Divide the unnormalized projections with option `""`.

Call `build_dijet_gen_comparisons(..., ratio_option="")` where it fits. Keep any later comparison of an already-constructed F/B curve to another curve separate from F/B construction.

### Fold a full CM eta histogram

Use `hist_analysis.python.data_distributions.forward_backward_from_full` where practical. Follow its project convention when implementing an equivalent extension:

1. Find the first positive-side bin with `axis.FindBin(0.0 + 0.001)`.
2. Map positive eta to Forward and the reflected negative eta bin to Backward.
3. Use `abs(eta)` bin edges on the result.
4. Limit the output to paired bins present on both sides of zero; do not invent or interpolate bins.
5. Copy each source bin's content and error into detached Forward and Backward TH1 histograms with `Sumw2()` enabled.
6. Divide Forward by Backward using option `""`.
7. Retain the Forward and Backward objects on the returned ratio when necessary to protect their PyROOT lifetimes.

Do not silently decide how to handle an eta bin that straddles zero or asymmetric binning if the established helper cannot represent it. Ask one focused question when the requested convention is not documented.

## Distinguish F/B from later ratios

An F/B histogram is `Forward / Backward`; its construction is always independent-error ROOT division.

A later ratio such as `(F/B)_variation / (F/B)_nominal` is a different observable. Keep its error option separate in code and configuration. Default it to standard independent errors. Only preserve an existing binomial setting for that later ratio when it is already documented and statistically justified by a numerator-subset relationship; never let that setting flow into F/B construction.

## Plotting and output

- Use shared styles and plotting helpers from `hist_analysis/python/root_style.py` and `hist_analysis/python/plotting.py`.
- Label the y axis `Forward / Backward` unless the surrounding workflow specifies a more precise established label.
- Use the existing configured x and y ranges and annotations.
- Write generated plots and derived ROOT files beneath `hist_analysis/output/` unless the user requests another location.
- Keep reusable projection, ratio, and plotting logic in `hist_analysis/python/`; do not add another legacy inline notebook implementation.

## Validation

Perform the smallest applicable checks from the repository documentation:

1. Check notebook JSON after notebook edits.
2. Run the documented Python syntax/compile check for changed Python helpers.
3. Execute the smallest affected workflow when its ROOT inputs and external libraries are available.
4. Confirm Forward and Backward have identical binning and remain unnormalized at division time.
5. Confirm the implementation passes `""` to the F/B `TH1::Divide` call and offers no binomial path.
6. On a small known input, compare representative bin contents and errors with ROOT's independent-propagation result, including a weighted bin when available.
7. Inspect zero-denominator behavior rather than replacing it with an undocumented convention.
8. Confirm expected output files, ROOT keys, dimensions, and relevant yields; check that unrelated results did not change.

If ROOT data or another prerequisite is unavailable, report exactly which runtime validation was skipped. Before handoff, run `git status --short`, identify only the files changed for the request, and state the expected physics/output effect.
