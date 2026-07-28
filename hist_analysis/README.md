# hist_analysis

Workspace for ROOT-file inventory, event diagnostics, inclusive-jet efficiency and selection studies, closure comparisons, JES/JER studies, beam-orientation checks, and unfolding prototypes built around the direction-specific and combined `pPb8160` outputs.

## Layout

- `notebooks/`: Jupyter workflows for inventory, event diagnostics, inclusive-jet validation, closures, JES/JER, beam orientation, and unfolding.
- `python/`: reusable PyROOT helpers for histogram I/O, projections, normalization, fitting, and plotting. `python/root_style.py` defines the shared `PlotStyle`, canvas geometry, axes, annotations, legends, and palette styling.
- `config/`: input locations, closure bins, the common eta range, and the standard dijet eta-cut index.
- `output/`: generated PDFs, optional PNGs, and derived ROOT files. This directory is ignored by Git.

## Data locations

The direction-dependent outputs are stored under:

- `/Users/gnigmat/cernbox/ana/pPb8160/exp/`
- `/Users/gnigmat/cernbox/ana/pPb8160/embedding/`
- `/Users/gnigmat/cernbox/ana/pPb8160/pythia/`

Within each generator directory, the direction-specific subdirectories are used for `Pbgoing` and `pgoing` outputs. Combined-orientation files live directly in the generator directory.

## Python environment

Use the project environment at:

```bash
/Users/gnigmat/work/cms/jetAnalysis/py-env
```

Example:

```bash
/Users/gnigmat/work/cms/jetAnalysis/py-env/bin/python -m jupyter notebook
```

## Starting point

Run `notebooks/01_inventory.ipynb` first to check the configured directory tree,
find ROOT files, and inspect top-level object keys.

Use `notebooks/02_event_histograms.ipynb` for event-level `pThat` and vertex-z
plots and for the gen/reco dijet overweight-protection diagnostic. The latter
uses the same discrete per-`pThat`-bin upper-tail threshold definition as
`macro/plotMcClosures.C::plotOverweightProtection`, followed by the configured
two-exponential threshold fits over 15–950 GeV. Its two-dimensional maps use
square canvases, the Bird palette, logarithmic z axes, shared palette geometry,
and optional threshold-fit overlays.

Use `notebooks/02_jet_efficiency_fakes.ipynb` for Lab-frame inclusive-jet
efficiency (`Ref matched / Gen inclusive`), matched rate
(`Reco matched / Reco inclusive`), and fake rate
(`Reco unmatched / Reco inclusive`). It draws direct two-dimensional rate maps,
eta-dependent rates in `config.histograms.SINGLE_JET_PT_BINS`, and pT-dependent
rates in configurable eta intervals. Numerator and denominator are projected
before each 1D division and are not normalized. All 1D curves use
`set_1d_style`; matched and fake rates use indices 0 and 1 on shared log-y
canvases. The default ROOT binomial errors describe weighted-MC efficiencies and
must not be interpreted as unweighted counting uncertainties.

Use `notebooks/02_jet_selection.ipynb` to compare Reco, matched Ref, and Gen
inclusive-jet eta shapes for the standard jet ID, track-maximum-only, and
no-selection stages in Lab and CM frames. It displays the underlying 2D maps and
then produces configured pT-interval projections with Reco/Gen and Ref/Gen ratio
panels. Integral, bin-width, and stored-yield modes and regular or binomial ratio
errors are configurable.

Use `notebooks/03_jet_JES_JER.ipynb` to extract JES and JER from Gaussian fits to
inclusive-jet response slices, compare the configured reco/gen smearing cases,
scan pT and eta dependence, and draw the corresponding two-dimensional response
maps. The notebook also normalizes JER-versus-eta curves to unit central-bin mean
within `-0.8 < eta < 0.8` and writes those derived histograms to a ROOT file in
`output/jes_jer/`.

Use `notebooks/04_systematics_beam_orientation.ipynb` to compare p-going and
Pb-going single-jet and dijet eta projections. Gen and Reco levels are enabled
by default; Ref is supported through the commented `LEVELS` configuration. Direction
comparisons show ratios to the configured nominal, while the same-direction
frame overlays intentionally do not form ratios. The flipped-lab and CM
direction comparisons also include the combined-orientation file and use it as
the ratio denominator.

Use `notebooks/05_unfold2D.ipynb` to construct a flattened dijet
`pTave`-eta response and run a Bayesian RooUnfold closure test for the configured
MC sample. Misses are represented by the response inefficiency, while fake
handling is enabled explicitly in `RooUnfoldBayes`. The notebook writes the
flattened input spectra, response and miss/fake diagnostics, unfolded result,
closure ratios, covariance matrix, serialized response, and configuration
metadata to
`output/unfold2D/<generator>_<direction>_unfold2D_eta_<cut>_iter_<iterations>.root`.
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

For `DIRECTION = 'combined'`, the unfolding notebook reads
`<generator>/<generator>_<stem>.root`. For `pgoing` or `Pbgoing`, it reads
`<generator>/<direction>/<generator>_<direction>_<stem>.root`, matching the
production output naming convention.

Use `notebooks/03_dijet_reco_to_gen_closures.ipynb` for the configurable
counterpart of `macro/plotMcClosures.C::plotDiJetClosures`. It projects the selected dijet
intervals from `config.histograms.DIJET_PTAVE_BINS` for every selected eta-cut index. It
compares full eta shapes to nominal Gen and compares forward/backward ratios to
the Gen forward/backward ratio. `NORMALIZATION='integral'` produces fractions per
bin whose bin-content sum is one; `bin_width` produces a density whose
bin-width-weighted integral is one. Projection normalization uses ROOT
`TH1::Scale`, while forward/backward and closure ratios use ROOT `TH1::Divide`.
Forward and backward inputs are not normalized before division.

Reco is enabled in red (`set_1d_style` index 0) and Gen in blue (index 1). Ref,
smeared Gen, or smeared Reco curves can be added by supplying their CM, Forward,
and Backward key templates. Individual full-distribution and F/B x-axis ranges
extend 0.1 beyond the selected eta cut; eta-cut overlays extend 0.1 beyond the
largest selected cut. Full-shape ratio limits, F/B limits, and F/B double-ratio
limits are configured separately.
The annotation records the generator, levels, CM frame, pTave interval, jet eta
acceptance, leading/subleading thresholds, and delta-phi requirement. The ROOT
binomial division option reproduces the macro when requested, but the notebook
can use independent error propagation for weighted forward/backward histograms.

PDF names are
`<generator>_<direction>_full_etaCM_<eta*10>_ptave_<low>_<high>.pdf` and
`<generator>_<direction>_fb_etaCM_<eta*10>_ptave_<low>_<high>.pdf`; for example,
the nominal 1.9 acceptance and configured 60--80 GeV interval use
`etaCM_19_ptave_60_80`.

The closure projection helpers identify pT and eta axes from their physical
ranges because some existing reference-jet files contain misleading axis
titles. Configured selection intervals are interpreted as half-open
`[low, high)` ranges.
