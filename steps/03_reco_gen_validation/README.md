# Step 03: Generator and reconstruction validation

**Status:** In progress; many diagnostic histograms and plots are implemented,
including configurable dijet closures and inclusive-jet efficiency, matched-rate,
fake-rate, and selection notebooks.

## Goal and problem

Demonstrate that selections, corrections, resolutions, orientation transforms,
and detector acceptance are understood before building the unfolding response.
The observed eta-dependent distortions must be separated into calibration,
resolution, acceptance, and plausible physics components.

## Recommended solution

Validate, separately by direction and pt-hat/sample range:

1. event selection and vertex reweighting;
2. jet-ID efficiency and its eta/pT dependence;
3. JEC closure in eta and pT;
4. nominal, up, down, and eta-dependent JER models;
5. leading/subleading pT, delta-phi, pT balance, and dijet eta;
6. generator/reconstructed efficiencies, migrations, misses, and fakes;
7. acceptance scans, treating \(|\eta_{\mathrm{CM}}|<1.9\) as nominal;
8. CM-boost scan around 0.465 and a defined pointing/angular-resolution study.

Plots already available in `macro/plotMcClosures.C` should be catalogued and
tested against their Python counterparts before either implementation is used
for a validation decision.

All generator/reconstruction checks are performed independently for Pb-going MC
and p-going MC. Data control checks that do not require unfolding are later
repeated independently for Pb-going data and p-going data; data and MC directions
must never be pooled to conceal a direction-specific discrepancy.

## Validation and exit criteria

- Every correction has a documented closure metric and uncertainty.
- No unexplained direction-dependent sign or normalization discrepancy remains.
- Both MC directions pass independently; a combined result cannot substitute for
  a failed or missing direction-level check.
- Efficiency/purity and migration behavior support the proposed response range.
- JER variations are physically defined and reproducible, not merely alternative
  arbitrary scale factors.
- Effects near acceptance boundaries are quantified before binning is frozen.

## Required figures

At minimum: JEC/JER closure maps, eta and pT resolutions, selection efficiencies,
miss/fake rates, direction comparisons, acceptance scans, boost scan, and dijet
kinematic closures. Every main plot is PDF+PNG.

## Unresolved question

“Pointing resolution” must be defined: jet-axis \(\eta/\phi\) resolution from
reco-to-gen matching, primary-vertex pointing, or another detector effect?

## Results

`hist_analysis/notebooks/03_mc_dijet_reco_to_gen_closures.ipynb` provides a configurable
counterpart to `plotDiJetClosures` for full eta shapes, forward/backward ratios,
and both observables divided by nominal Gen. It uses the explicit intervals in
`config.histograms.DIJET_PTAVE_BINS`, canonical eta-cut indices, ROOT
scaling/division, and explicit normalization choices. The integral-normalized
projections and unnormalized F/B
ratios were checked bin-by-bin against direct ROOT operations for the combined
embedding file. This is an implemented diagnostic, not a direction-specific
physics validation: p-going and Pb-going acceptance criteria and pass/fail
results remain outstanding.

`hist_analysis/notebooks/02_mc_jet_efficiency_fakes.ipynb` now provides direct 2D
and projected 1D Lab-frame efficiency, matched-rate, and fake-rate diagnostics.
For the configured combined embedding file, clean-kernel execution produced the
expected plots; matched plus fake closed to unity within floating-point precision,
and matched plus unmatched reconstructed yields reproduced the inclusive-reco
map. `02_jet_selection.ipynb` separately compares Reco, matched Ref, and Gen
shapes across jet-selection stages. These are implemented combined-sample
diagnostics, not substitutes for the required p-going and Pb-going validations.
