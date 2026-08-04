# Dijet pseudorapidity measurement plan

**Plan status:** Approved for incremental execution on 2026-07-14. Individual
steps still require their documented completion review before the next dependent
step begins.

## Measurement objective

Measure the self-normalized dijet pseudorapidity distribution

\[
\frac{1}{N_{\mathrm{dijet}}}\frac{dN_{\mathrm{dijet}}}{d\eta_{\mathrm{dijet}}}
\]

in several final average-dijet transverse-momentum intervals, where

\[
p_T^{\mathrm{ave}}=\frac{p_{T,1}+p_{T,2}}{2},\qquad
\eta_{\mathrm{dijet}}=\frac{\eta_{1,\mathrm{CM}}+\eta_{2,\mathrm{CM}}}{2}.
\]

The measurement uses pPb collisions at 8.16 TeV, AK4PF jets with radius 0.4,
and combines Pb-going and p-going results only after direction-specific
corrections. The nominal selection is:

- `jetIdTightLeptVeto`;
- leading-jet \(p_T\geq 50\) GeV;
- subleading-jet \(p_T\geq 40\) GeV;
- \(|\Delta\phi|\geq 2\pi/3\);
- \(|\eta_{\mathrm{CM}}|\leq 1.9\) for each selected jet.

Other jet-acceptance limits already present in the processing macro are retained
as diagnostic variations. The nominal correction is intended to be two
dimensional in \((p_T^{\mathrm{ave}},\eta_{\mathrm{dijet}})\). Bayesian
unfolding is the current candidate but is not yet approved as the nominal method.

## Scientific workflow

Work proceeds in numbered, approval-gated steps. A step may begin only after its
design is accepted. Completion requires code review, reproducible commands,
validation results, and updated documentation. “Implemented” below means that
code exists; it does not imply that the physics result has been validated.

| Step | Purpose | Current assessment |
|---|---|---|
| [00](00_analysis_contract/README.md) | Freeze conventions and audit current behavior | In progress; Pb-going embedded MC validated |
| [01](01_reproducible_baseline/README.md) | Establish a reproducible build and small MC run | In progress; AppleClang Release build and one Pb-going embedded-MC run passed |
| [02](02_mc_production/README.md) | Inventory and validate existing direction-specific MC products | Products exist; provenance and validation pending |
| [03](03_reco_gen_validation/README.md) | Validate selections, weights, JEC, JER, and orientations | Dijet closure and inclusive-jet efficiency/fake diagnostics implemented; direction-level acceptance criteria still missing |
| [04](04_response_and_binning/README.md) | Validate the 2D response and choose analysis binning | Flattened-response prototype implemented with provisional edge bins; design incomplete |
| [05](05_mc_unfolding_closure/README.md) | Select and validate unfolding using MC closure | `hist_analysis/notebooks/05_unfold2D.ipynb` prototype exists; no approved independent closure result |
| [06](06_systematic_model/README.md) | Quantify systematic variations | Several variations implemented; model incomplete |
| [07](07_data_correction_and_combination/README.md) | Correct data and combine beam directions | Not started; gated by MC closure |
| [08](08_final_results_and_release/README.md) | Final plots, documentation, and reproducible release | Not started |

## Common requirements for every step

Each step document contains facts, recommendations, unresolved questions, an
implementation proposal, validation criteria, and a results section. Results
must record the Git commit, input sample identifiers, configuration, commands,
ROOT and RooUnfold versions, random seeds where applicable, and numerical tests.

Main figures must be saved in both PDF and PNG. ROOT canvases or source
histograms must also be retained where practical. Generated ROOT files and large
plot collections remain outside Git; the Markdown document records their stable
location and provenance. Small tables, configuration files, and final figures
may be versioned after review.

No experimental data unfolding begins before Step 05 is approved as complete.
Final \(p_T^{\mathrm{ave}}\) reporting intervals are derived from validated fine
2D/3D histograms; they are not fixed by the present prototype intervals.

## Mandatory orientation workflow

Every applicable analysis stage must be executed and documented independently
for these four channels:

1. Pb-going embedded MC;
2. p-going embedded MC;
3. Pb-going data;
4. p-going data.

Direction-specific orientation transformations, JEC/JER treatment, response
matrices, corrections, validation plots, systematic variations, and covariance
matrices must be produced before combination. MC directions are combined only for
an explicitly defined combined-MC check. Data directions are combined only after
each direction has passed MC closure and direction-specific data validation.
No histogram may be combined merely by adding files whose names appear to have
compatible orientations; the transformations implemented by `etaLab(...)` and
`etaCM(...)` must be applied first.

All required data and MC productions already exist. They are inputs to inventory
and validate, not jobs to regenerate. A targeted reprocessing is proposed only if
a documented correctness defect cannot be addressed downstream, and it requires
analyst approval before execution.

## Status vocabulary

- **Proposed:** design awaits approval.
- **In progress:** design approved and work started.
- **Implemented, unvalidated:** code exists but has not passed documented tests.
- **Complete:** all exit criteria passed and results documented.
- **Blocked:** a named external input or decision is required.

## Current unresolved project-level decisions

1. Inventory of existing input datasets, production versions, event counts,
   checksums or immutable identifiers, and file locations.
2. Fine internal bin edges and the quantitative rule for merging them into final
   \(p_T^{\mathrm{ave}}\) intervals.
3. MC split or independent samples used for training, test closure, and systematic
   comparisons.
4. RooUnfold algorithm, regularization-selection rule, and covariance convention.
5. Pointing-resolution definition and available inputs.
6. Final systematic-source list and correlation model between directions and bins.

The authoritative starting convention for stored orientation, physical beam
direction, and sign is the current implementation of `etaLab(...)` and
`etaCM(...)` in `processForestSimple.C`. Step 00 must convert this implementation
into an explicit truth table and validate all four data/MC × beam-direction cases.
