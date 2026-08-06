# hist_analysis

Workspace for ROOT-file inventory, event diagnostics, inclusive-jet efficiency and selection studies, closure comparisons, JES/JER studies, beam-orientation checks, and unfolding prototypes built around the direction-specific and combined `pPb8160` outputs.

## Layout

- `notebooks/`: Jupyter workflows for inventory, event diagnostics, inclusive-jet validation, closures, JES/JER, beam orientation, and unfolding.
- `python/`: reusable PyROOT helpers for histogram I/O, projections, normalization, fitting, and plotting. `python/root_style.py` defines the shared `PlotStyle`, canvas geometry, axes, annotations, legends, and palette styling.
- `config/`: input locations, closure bins, the common eta range, and the standard dijet eta-cut index.
- `output/`: generated PDFs, optional PNGs, and derived ROOT files. This directory is ignored by Git.

Notebook filenames use `_mc_` for MC-only validation workflows and `_data_`
for data-only workflows. The numeric prefix groups the analysis stage; the
data-versus-MC comparison is stage 07.

## Data locations

The direction-dependent outputs are stored under:

- `/Users/gnigmat/cernbox/ana/pPb8160/exp/`
- `/Users/gnigmat/cernbox/ana/pPb8160/embedding/`
- `/Users/gnigmat/cernbox/ana/pPb8160/pythia/`

Within each generator directory, the direction-specific subdirectories are used for `Pbgoing` and `pgoing` outputs. Combined-orientation files live directly in the generator directory.

## Python environment

Use the project environment from the repository root:

```bash
py-env/bin/python
```

Example:

```bash
py-env/bin/python -m jupyter notebook
```

Notebook setup is shared through `python/notebook_setup.py`. It locates PyROOT
from the active environment or dynamically through `root-config`; notebooks do
not add version- or machine-specific ROOT paths. The unfolding notebooks load
RooUnfold separately and use the adjacent `../RooUnfold` checkout by default;
set `ROOUNFOLD_ROOT` to override that location.

## Starting point

Run `notebooks/01_inventory.ipynb` first to check the configured directory tree,
find ROOT files, and inspect top-level object keys.

Use `notebooks/02_mc_event_histograms.ipynb` for event-level `pThat` and vertex-z
plots and for the gen/reco dijet overweight-protection diagnostic. The latter
uses the same discrete per-`pThat`-bin upper-tail threshold definition as
`macro/plotMcClosures.C::plotOverweightProtection`, followed by the configured
two-exponential threshold fits over 15–950 GeV. Its two-dimensional maps use
square canvases, the Bird palette, logarithmic z axes, shared palette geometry,
and optional threshold-fit overlays.

Use `notebooks/02_mc_jet_efficiency_fakes.ipynb` for Lab-frame inclusive-jet
efficiency (`Ref matched / Gen inclusive`), matched rate
(`Reco matched / Reco inclusive`), and fake rate
(`Reco unmatched / Reco inclusive`). It draws direct two-dimensional rate maps,
eta-dependent rates in `config.histograms.SINGLE_JET_PT_BINS`, and pT-dependent
rates in configurable eta intervals. Numerator and denominator are projected
before each 1D division and are not normalized. All 1D curves use
`set_1d_style`; matched and fake rates use indices 0 and 1 on shared log-y
canvases. The default ROOT binomial errors describe weighted-MC efficiencies and
must not be interpreted as unweighted counting uncertainties.

Use `notebooks/02_mc_jet_selection.ipynb` to compare Reco, matched Ref, and Gen
inclusive-jet eta shapes for the standard jet ID, track-maximum-only, and
no-selection stages in Lab and CM frames. It displays the underlying 2D maps and
then produces configured pT-interval projections with Reco/Gen and Ref/Gen ratio
panels. Integral, bin-width, and stored-yield modes and regular or binomial ratio
errors are configurable.

Use `notebooks/02_mc_jet_JES_JER.ipynb` to extract JES and JER from Gaussian fits to
inclusive-jet response slices, compare the configured reco/gen smearing cases,
scan pT and eta dependence, and draw the corresponding two-dimensional response
maps. The notebook also normalizes JER-versus-eta curves to unit central-bin mean

within `-0.8 < eta < 0.8` and writes those derived histograms to a ROOT file in
`output/jes_jer/`.

Use `notebooks/02_mc_jet_JEC_check.ipynb` to compare raw and JEC-corrected MC
inclusive-jet eta projections in configurable pT intervals. It overlays the
Pb-going and p-going yields, plots corrected/raw ratios, and compares the two
beam directions with independent-error propagation.

Use `notebooks/03_mc_dijet_embedding_vs_pythia.ipynb` to compare Embedding and
PYTHIA Reco, matched-Ref, and Gen dijet distributions for a configurable stored
jet-eta acceptance and every half-open interval in
`config.histograms.TEST_DIJET_PTAVE_BINS`. Reco, Ref, and Gen comparisons are
launched from separate cells. Each full CM-pseudorapidity distribution is
normalized independently to unit bin-content sum, then overlaid with a
PYTHIA/Embedding ratio panel. The notebook derives one common full-distribution
y-axis range from all levels, samples, and pTave intervals before drawing. It
also overlays the unnormalized Forward/Backward ratios with a common configured
`FB_RANGE` and draws `(F/B)_PYTHIA / (F/B)_Embedding`. Both F/B construction and
the later sample-comparison ratios use ROOT independent-error propagation; no
binomial option is exposed. Combined and direction-specific inputs use the
standard MC filename resolvers. Set `DIJET_EMBEDDING_VS_PYTHIA_OUTPUT_DIR` to
redirect output from `output/dijet_embedding_vs_pythia/`, and `SAVE_PNG = True`
to write PNG files alongside the default PDFs.

Use `notebooks/03_mc_dijet_smearing_effect.ipynb` to compare nominal Gen, Gen with
the default pT smearing, Gen with eta-dependent smearing, nominal Reco, and Reco
with eta-dependent default JER smearing for one configured eta-cut index and
every interval in `config.histograms.DIJET_PTAVE_BINS`. `CURVE_CATALOG` holds
the available histogram keys, styles, and raw scale factors. Nominal Gen,
eta-dependent smeared Gen, and nominal Reco are required; add or remove other
catalog entries with `OPTIONAL_CURVE_LABELS`. Before normalization, raw
full-distribution overlays and ratios are written after scaling accumulated generator-smeared curves by
`1/N_SMEAR_RUNS` (20 by default) and accumulated Reco-smeared curves by
`1/N_RECO_SMEAR_RUNS` (10 by default). Normalized full-shape overlays, ratios to
Gen and eta-dependent smeared Gen, unnormalized forward/backward overlays, and
their ratios to both references are also written, for nine separate single-panel
canvases per pTave interval using the shared ROOT styles. Generator-level and
Reco-like full-distribution ratios have separate configurable ROOT division
options. Forward/Backward construction always uses independent propagation.
The later ratio of one F/B histogram to another has a separate user-configurable
error option, `''` or `'B'`.

Use `notebooks/07_dijet_data_vs_mc.ipynb` to compare nominal reconstructed-data,
eta-dependent JER-default reconstructed-MC, and nominal generator-level-MC dijet
pTave spectra for MinimumBias, Jet60, Jet80, and Jet100. The MC-reco curve uses
`hRecoDijetPtEtaCMJerDefExtra`; data reco and MC gen use
`hRecoDijetPtEtaCM` and `hGenDijetPtEtaCM`, respectively. Each curve is
normalized independently over the configured half-open pTave interval, and a
separate comparison is produced for every configured jet-eta acceptance. The
notebook also projects and unit-normalizes pseudorapidity shapes in every
half-open interval from `config.histograms.TEST_DIJET_PTAVE_BINS`, producing
data/MC-gen and MC-reco/MC-gen ratios for each trigger and eta acceptance. The
default ratio ranges are 0.75--1.25 for pTave and 0.65--1.15 for eta, configured
independently with `PTAVE_RATIO_RANGE` and `ETA_RATIO_RANGE`. The input data
directory and output directory can be overridden with `PPB_DATA_DIR` and
`DIJET_DATA_VS_MC_OUTPUT_DIR`. Set `DATA_DIRECTION` to `combined`,
`Pbgoing`, or `pgoing`, and `DATA_SELECTION` to `jetId`, `trkMax`, or `noSel`;
the same trigger/direction/selection filename resolver used by the other data
notebooks selects all four data inputs.

Use `notebooks/03_mc_dijet_reco_smeared_to_gen_closures.ipynb` for the focused
closure comparison between nominal Gen and Reco smeared with the default
eta-dependent JER. It writes normalized full-distribution and unnormalized
Forward/Backward comparisons for every configured pTave interval. As in the
other closure notebooks, F/B construction always uses independent errors;
full-distribution and F/B-to-F/B comparison error options are configured
separately.

Use `notebooks/04_systematics_beam_orientation.ipynb` to compare p-going and
Pb-going single-jet and dijet eta projections. Gen and Reco levels are enabled
by default; Ref is supported through the commented `LEVELS` configuration. Direction
comparisons show ratios to the configured nominal, while the same-direction
frame overlays intentionally do not form ratios. The flipped-lab and CM
direction comparisons also include the combined-orientation file and use it as
the ratio denominator.

Use `notebooks/04_systematics_JER.ipynb` to estimate the reconstructed-dijet
JER shape uncertainty with embedding or Pythia MC. For every configured pTave
interval it overlays the JER Up, Down, and default CM eta projections after
scaling each to unit integral, and overlays their unnormalized Forward/Backward
ratios. Dedicated Up/Def and Down/Def plots are produced for both observables.
Forward/Backward construction always uses independent error propagation;
full-shape and subsequent F/B-to-F/B comparison ratios independently support
ROOT's standard or `B` division option.

Use `notebooks/04_systematics_JEU.ipynb` to estimate the reconstructed-dijet
JEU shape uncertainty in MinimumBias, Jet60, Jet80, and Jet100 data. Each
trigger has its own configured pTave intervals. The notebook compares
unit-integral JEU Up, Down, and default CM eta projections and their
unnormalized Forward/Backward ratios, then draws dedicated Up/Def and Down/Def
comparisons. Forward/Backward construction always uses independent errors; the
later comparison-ratio error options are configured separately. Set
`DATA_DIRECTION` and `DATA_SELECTION` to choose the shared four-trigger data
production, and `DIJET_JEU_SYSTEMATICS_OUTPUT_DIR` to redirect its generated
plots.

Use `notebooks/05_unfold2D.ipynb` to construct a flattened dijet
`pTave`-eta response and run a Bayesian RooUnfold closure test for the configured
MC sample using the eta-dependent JER-default measured distribution and its
corresponding response, miss, fake, and pair-classification objects. Misses are
represented by the response inefficiency, while fake handling is enabled
explicitly in `RooUnfoldBayes`. The notebook writes the
flattened input spectra, response and miss/fake diagnostics, unfolded result,
closure ratios, covariance matrix, serialized response, and configuration
metadata to
`output/unfold2D/<generator>_<direction>_unfold2D_jerDefExtra_eta_<cut>_iter_<iterations>.root`.
Set `DIJET_UNFOLD2D_OUTPUT_DIR` to redirect these generated artifacts without
editing the notebook configuration.
The projection/response diagnostics are produced for every configured pTave
interval; the flattened-response and final closure canvases are also saved as
PDFs in the same directory. For every pTave interval, the notebook extracts the
unfolded eta distribution, overlays it with the original gen and reco
projections, and plots reco/gen and unfolded/gen ratios.
Those component and ratio histograms are stored in the ROOT output. Unfolding
plots use the semantic styles from `python/root_style.py`: gen blue, reco red,
miss orange, fake green, and unfolded black. Set `SAVE_PNG = True` in the
notebook configuration to also write matching PNG files.

Before constructing `RooUnfoldResponse`, the notebook scales cloned truth,
measured, and response histograms by the common `RESPONSE_SCALE`. This keeps
small weighted response entries above RooUnfold's absolute matrix-sanitization
threshold. The unfolded histogram is scaled back by `1/RESPONSE_SCALE`, and its
covariance by `1/RESPONSE_SCALE^2`, before results are plotted or written.

Use `notebooks/05_unfold2D_mc_direction.ipynb` for the beam-direction closure
test. It constructs the response and its training truth/reco marginals, misses,
fakes, and pair classification from Pb-going embedding, then unfolds the
independent p-going embedding eta-dependent JER-default reco distribution and
compares it with p-going generator truth. The notebook uses the same flattened
bin mapping, RooUnfold response scaling, explicit fake handling, covariance,
plotting styles, and ROOT output objects as `05_unfold2D.ipynb`. Its output tag
records both the response-training and test directions. Configure
`FLATTENED_RATIO_TO_GEN_Y_RANGE` and `ETA_RATIO_TO_GEN_Y_RANGE` to set the
ratio-to-gen y-axis ranges for the flattened closure and per-pTave eta panels,
respectively. Set `DRAW_GRID` to enable or disable grid lines on all diagnostic
and closure canvases.

Use `notebooks/06_data_jet_distributions.ipynb` and
`notebooks/06_data_dijet_distributions.ipynb` for trigger-by-trigger
reconstructed-data beam-orientation checks. Set `SELECTION` to `jetId`,
`trkMax`, or `noSel`; the notebooks resolve the matching combined, Pb-going,
and p-going ROOT files for MinimumBias, Jet60, Jet80, and Jet100. Combined
filenames are `<trigger>_ak4_<selection>.root`; direction-specific filenames
are `<direction>/<trigger>_<direction>_ak4_<selection>.root`, with lowercase
trigger stems (`mb`, `jet60`, `jet80`, and `jet100`). For each trigger they draw raw two-dimensional maps in
the unflipped-Lab, flipped-Lab, and CM frames, using one common z-axis range for
all three orientations. The inclusive-jet notebook also draws selected raw-pT
versus unflipped detector-eta diagnostics from
`hRecoInclusiveJetRawPtEtaLabUnflipped`. They overlay pT or pTave spectra and
normalized eta projections with Pb-going/combined and p-going/combined ratio
panels. The pT
spectra are normalized in the configured 110--130 GeV interval. Inclusive-jet
and dijet projection intervals are independently configurable for MinimumBias,
Jet60, Jet80, and Jet100; inclusive-jet pT spectra also have a configurable eta
selection. Set `PPB_DATA_DIR`,
`DATA_JET_OUTPUT_DIR`, or `DATA_DIJET_OUTPUT_DIR` to override input and output
locations.

For `DIRECTION = 'combined'`, the unfolding notebook reads
`<generator>/<generator>_<stem>.root`. For `pgoing` or `Pbgoing`, it reads
`<generator>/<direction>/<generator>_<direction>_<stem>.root`, matching the
production output naming convention.

Use `notebooks/03_mc_dijet_reco_to_gen_closures.ipynb` for the configurable
counterpart of `macro/plotMcClosures.C::plotDiJetClosures`. It projects the selected dijet
intervals from `config.histograms.DIJET_PTAVE_BINS` for every selected eta-cut index. It
compares full eta shapes to nominal Gen and compares forward/backward ratios to
the Gen forward/backward ratio. Before normalizing the shapes, it also writes a
raw-yield overlay of Reco, matched Ref, and Gen with Reco/Gen and Ref/Gen in the
lower panel. `NORMALIZATION='integral'` produces fractions per bin whose
bin-content sum is one; `bin_width` produces a density whose bin-width-weighted
integral is one. Projection normalization uses ROOT `TH1::Scale`, while
forward/backward and closure ratios use ROOT `TH1::Divide`. Forward and backward
inputs are not normalized before division.

Reco is enabled in red (`set_1d_style` index 0) and Gen in blue (index 1). Ref,
smeared Gen, or smeared Reco curves can be added by supplying their CM, Forward,
and Backward key templates. Individual full-distribution and F/B x-axis ranges
extend 0.1 beyond the selected eta cut; eta-cut overlays extend 0.1 beyond the
largest selected cut. Full-shape ratio limits, F/B limits, and F/B double-ratio
limits are configured separately.
The annotation records the generator, levels, CM frame, pTave interval, jet eta
acceptance, leading/subleading thresholds, and delta-phi requirement. The ROOT
Forward/Backward ratios always use independent error propagation; ROOT's
binomial division option is invalid because Forward and Backward are separate
weighted populations. Options for subsequent closure ratios are configured
separately and may use binomial propagation only when their numerator is a
statistical subset of the denominator.

PDF names are
`<generator>_<direction>_full_etaCM_<eta*10>_ptave_<low>_<high>.pdf` and
`<generator>_<direction>_fb_etaCM_<eta*10>_ptave_<low>_<high>.pdf`; for example,
the nominal 1.9 acceptance and configured 60--80 GeV interval use
`etaCM_19_ptave_60_80`.

The closure projection helpers identify pT and eta axes from their physical
ranges because some existing reference-jet files contain misleading axis
titles. Configured selection intervals are interpreted as half-open
`[low, high)` ranges.
