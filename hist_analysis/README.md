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

Every notebook begins with a detailed workflow guide defining its observable,
normalization, ratios or fits, statistical assumptions, and interpretation.
Each code cell also starts with a `Cell role` and `Interpretation` comment so
the executable step can be read together with the preceding physics Markdown.
Across the project, a histogram bin content denotes the stored weighted yield
`sum(w)`, while its usual variance is `sum(w^2)`. Unit-area normalization,
bin-width density normalization, and unnormalized yields are distinct and are
identified explicitly. Forward/backward ratios always use independent-error
propagation because their positive- and reflected-negative-eta bins are not a
subset relationship.

## Data locations

The direction-dependent outputs are stored under:

- `/Users/gnigmat/cernbox/ana/pPb8160/exp/`
- `/Users/gnigmat/cernbox/ana/pPb8160/embedding/`
- `/Users/gnigmat/cernbox/ana/pPb8160/pythia/`

Within each generator directory, the direction-specific subdirectories are used for `Pbgoing` and `pgoing` outputs. Combined-orientation files live directly in the generator directory.

Per-forest pp CP5 outputs from `pp_analyzeMc_cp5tune.py` live directly under
`/Users/gnigmat/cernbox/ana/pPb8160/pythia/cp5tune/` and use the suffix
`_hist.root`; their matching execution logs use `_hist.log`. These files retain
the same histogram keys as pPb outputs, but Lab and CM eta are both stored eta,
and weighted histograms contain sums of the input `hiEvtAnalyzer/HiTree/weight`
values. They are not part of the direction-combination workflows unless a
notebook input is explicitly configured to use them.

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

ROOT figures produced by MC, systematic, and simulation-only unfolding
notebooks default to a `CMS Simulation` left header and
`pPb 174.6 nb^{-1} (8.16 TeV)` on the right. Notebooks whose filenames contain
`data` instead default to `CMS Preliminary` with the same pPb luminosity and
energy; their configurable left label can later be changed to `CMS`. Notebooks expose
`CMS_HEADER_LEFT` and `CMS_HEADER_RIGHT` in their setup configuration and enable
the header through the shared `root_style.configure_cms_header` save path, so
both helper-generated and notebook-local canvases receive the same header.
The standalone pp tune workflow, `03_pythia_standalone_check.ipynb`, instead
defaults to `pp (8.16 TeV)` for both ROOT and mplhep output. Shared literal
labels live in `python/plot_labels.py`; mplhep helpers otherwise default to the
pPb MC header.

Every ROOT canvas saved through `python/root_style.py` can also be re-rendered
with mplhep while retaining the original ROOT output. Set the environment
variable `HIST_ANALYSIS_SAVE_MPLHEP=1` before starting Jupyter, or call
`configure_mplhep_export(True)` in a notebook. Python PDFs are written to an
`mplhep/` directory beside the ROOT PDFs. The canvas exporter handles TH1
plots, overlays, ratio pads, TH2 maps, and recursively nested multi-panel
canvases; set `HIST_ANALYSIS_SAVE_MPLHEP_PNG=1` or pass `save_png=True` to its
configuration function for PNG copies.
Every notebook exposes `SAVE_MPLHEP` and `SAVE_MPLHEP_PNG` in its setup cell;
the default is `False` except in the standalone PYTHIA tune workflow, whose
existing dedicated mplhep export remains enabled by default.
The generic exporter copies the final displayed x- and y-axis limits from each
ROOT pad, including log-scale limits and independent limits in multi-panel
canvases. A specific call can override them with `mplhep_x_range` or
`mplhep_y_range` in `save_canvas`; the standalone histogram-aware mplhep
helpers retain their existing explicit `x_range`, `y_range`, and `ratio_range`
arguments.
The shared mplhep curve cycle supports eight simultaneous curves and alternates
filled and open circles, triangles, squares, and diamonds. CP5 remains a filled
red circle and CUETP8M1 an open blue circle regardless of their overlay order.

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

Use `notebooks/02_mc_jet_selection.ipynb` to compare Reco, matched Ref, RefSel,
and Gen inclusive-jet eta shapes for the standard jet ID, track-maximum-only,
and no-selection stages in Lab and CM frames. It displays the underlying 2D maps,
then produces configured pT-interval full-distribution overlays and separate
ratio canvases containing Reco/Gen, Ref/Gen, and RefSel/Gen. Ratios always use
ROOT option `B` binomial errors. Integral, eta-range, bin-width, and stored-yield
modes are supported; `LOG_Y` controls the full-distribution y-axis scale. CM
1D canvases include solid black boundaries at `|eta_CM|=2.4` and dashed black
boundaries at `|eta_CM|=1.9`, with both boundary types listed in the legend.
The MC outputs also provide RefSel inclusive-jet histograms for these three
selection stages in the unflipped-Lab, Lab, and CM frames.

Use `notebooks/02_mc_jet_JES_JER.ipynb` to extract JES and JER from Gaussian fits to
inclusive-jet response slices, compare the configured reco/gen smearing cases,
scan pT and eta dependence, and draw the corresponding two-dimensional response
maps. The notebook also normalizes JER-versus-eta curves to unit central-bin mean

within `-0.8 < eta < 0.8` and writes those derived histograms to a ROOT file in
`output/jes_jer/`.

Use `notebooks/02_mc_jet_JEC_check.ipynb` to compare raw and JEC-corrected MC
inclusive-jet eta projections in configurable pT intervals. Each raw or
corrected histogram is projected in its own pT coordinate, exposing the bin
migration introduced by JEC. For each interval, the notebook writes a square
four-distribution overlay and a two-panel figure containing the Pb-going and
p-going corrected/raw ratios and their Pb-going/p-going double ratio. Select
the production and jet selection with `GENERATOR` and `FILE_STEM`, limit every
eta axis with `X_AXIS_RANGE`, and configure the distribution, corrected/raw,
and double-ratio y ranges independently. `CORRECTED_RAW_RATIO_OPTION` and
`CORRECTED_RAW_DIRECTION_RATIO_OPTION` accept ROOT `TH1::Divide` options `''`
(independent-error propagation) or `'B'` (binomial errors). Outputs use
`<generator>_<jetSelection>_<overlay|ratio>_pt_<low>_<high>.pdf` beneath
`output/mc_jet_JEC_check/`; set `SAVE_PNG = True` to also write PNG files.

Use `notebooks/03_mc_dijet_embedding_vs_pythia.ipynb` to compare Embedding and
standard PYTHIA Reco, matched-Ref, and Gen dijet distributions for a
configurable stored jet-eta acceptance. Set `INCLUDE_CP5 = True` to also load
the configurable pp `CP5_FILE` and produce standard-PYTHIA/CP5 and
Embedding/CP5 comparisons. `CP5_INPUT_SOURCE` accepts `cernbox`, `local`, or
`auto`; `local` uses
`$HOME/cernbox/ana/pPb8160/pythia/cp5tune/cp5_gen_dijet_histograms.root`, while
`CP5_DIJET_COMPARISON_FILE` overrides every default. Because the local
standalone file contains generator histograms only, its CP5 curves are added
only to the Gen section; Reco and matched-Ref retain the Embedding/PYTHIA
comparison. `PTAVE_BIN_SET` selects every half-open interval in
either `config.histograms.DIJET_PTAVE_BINS` or `TEST_DIJET_PTAVE_BINS`. Reco,
Ref, and Gen comparisons are launched from separate cells. For each pair, the
full CM-pseudorapidity distributions are normalized independently to unit
bin-content sum and overlaid above a ratio of the second sample to the first.
The notebook derives one common full-distribution y-axis range from all levels,
enabled samples, and pTave intervals before drawing. It also overlays the
unnormalized Forward/Backward ratios with a common configured `FB_RANGE` and
draws the corresponding pairwise F/B ratio. Both F/B construction and later
sample comparisons use ROOT independent-error propagation; no binomial option
is exposed. Combined and direction-specific pPb inputs use the standard MC
filename resolvers. Set `DIJET_EMBEDDING_VS_PYTHIA_OUTPUT_DIR` to redirect
output from `output/dijet_embedding_vs_pythia/`, and `SAVE_PNG = True` to write
PNG files alongside the default PDFs.

Use `notebooks/03_mc_fb_generator_check.ipynb` to inspect one Embedding,
standard PYTHIA pPb, or PYTHIA CP5 production at generator level. For either
the standard or compact test pTave bins, it overlays the independently
normalized CM- and Lab-frame dijet pseudorapidity shapes with a Lab/CM ratio
panel (default range 0.5--1.5), and draws the CM Forward/Backward ratio. The
latter is formed from the unnormalized stored
Forward and Backward projections with independent-error propagation; no
binomial option is exposed. `DIRECTION` selects Pb-going, p-going, or combined
pPb input, while `CP5_FILE` configures the pp input. Outputs go beneath
`output/mc_fb_generator_check/` by default.

Use `notebooks/03_pythia_standalone_check.ipynb` to compare the standalone pp
8.16 TeV PYTHIA CP5 and CUETP8M1 generator samples. The defaults resolve to
`$HOME/cernbox/ana/pPb8160/pythia/cp5tune/cp5_gen_dijet_histograms.root` and
`$HOME/cernbox/ana/pPb8160/pythia/cuep8m1tune/cuep8m1_gen_dijet_histograms.root`;
environment variables can override either file. `TUNES` enables both samples or
selects one for a single-tune check. The notebook overlays the tunes and draws
CUETP8M1/CP5 ratios for inclusive-jet and dijet full-eta and Forward/Backward
distributions. It also compares inclusive-jet pT spectra in five eta windows
(`[-3.2,-2.4)`, `[-2.4,-1.4)`, `[-1.4,1.4)`, `[1.4,2.4)`, and
`[2.4,3.2)`), event-weight distributions, and separate unweighted and weighted
pThat distributions. Full eta shapes use the selected display normalization. F/B
always divides unnormalized forward and backward weighted yields, and both F/B
and tune-to-tune ratios use independent-error propagation. Outputs go beneath
`output/pythia_standalone_check/` by default. Set `SAVE_MPLHEP = True` to retain
all ROOT figures and additionally export the final 1D comparisons, ratios, and
2D response maps as 8 x 8 inch mplhep CMS-style PDFs beneath the `mplhep/`
subdirectory. `MPLHEP_STYLE` configures the figure, axis, tick, annotation,
legend, CMS-label, marker, and colorbar font sizes; `annotation_y` controls the
vertical position of in-frame descriptive text. CP5 uses filled red circles,
CUETP8M1 uses open blue circles, and six additional deterministic color/marker
combinations support overlays with up to eight curves. `SAVE_MPLHEP_PNG` adds
200 DPI PNG copies. Dijet x-axis limits extend one projected bin width beyond
the configured absolute-eta selection so the final accepted bin is fully shown.
`MPLHEP_CMS_HEADER_LEFT` accepts `CMS Simulation`, `CMS Preliminary`, or `CMS`.
`MPLHEP_CMS_HEADER_RIGHT` accepts a math-text-capable label such as
`pp (8.16 TeV)` or `pPb 174.6 nb$^{-1}$ (8.16 TeV)`.
The notebook applies the same configurable header and tune semantics to its
ROOT canvases through `ROOT_PLOT_STYLE` and `ROOT_CURVE_STYLES`: CP5 is a filled
red circle and CUETP8M1 is an open blue circle. These settings are passed as
optional shared-plotting arguments, so other notebooks retain their defaults.
The default `DIJET_Y_RANGE = (0.00, 0.14)` fixes a common y-axis range for all
full dijet eta distributions. `MPLHEP_STYLE.capsize` controls the end-cap size
for both vertical statistical errors and horizontal bin-width bars.
Set `FULL_ETA_Y_RANGE` to an ordered `(low, high)` pair to enforce one common
y-axis range across every inclusive-jet and dijet full-eta figure; when unset,
the existing `INCLUSIVE_Y_RANGE` and `DIJET_Y_RANGE` controls remain active.
The mplhep renderer omits markers and error bars for bins whose content is zero
without modifying the source ROOT histograms. Its one-dimensional markers carry
horizontal error bars spanning half the bin width in each direction, in
addition to the stored vertical statistical uncertainties.
For each enabled tune it also draws the weighted
`hGenDijetPtAveOverPTHatVsPTHat` map and applies `TH2::FitSlicesY` to extract
the Gaussian mean of `pTave/pTHat` in every pTHat bin. The configurable
`FIT_SLICES_MIN_ENTRIES` controls the minimum slice population; the response
map and fitted-mean display ranges are configured independently. These
diagnostics are filled before the final dijet pT and delta-phi cuts and are
intended to expose rare upward migrations, not to impose an event veto.

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

Use `notebooks/07_dijet_data_vs_mc.ipynb` to compare nominal reconstructed data,
eta-dependent JER-default reconstructed embedding MC, generator-level embedding
MC, and generator-level PYTHIA CP5 dijet pTave spectra for MinimumBias, Jet60,
Jet80, and Jet100. The MC-reco curve uses
`hRecoDijetPtEtaCMJerDefExtra`; data reco and embedding MC gen use
`hRecoDijetPtEtaCM` and `hGenDijetPtEtaCM_{eta_cut_index}`, respectively. Each
curve is normalized independently over the configured half-open pTave interval.
The CP5 generator distribution is `hGenDijetPtEtaCM`. Set
`CP5_INPUT_SOURCE` to `cernbox`, `local`, or `auto`; `local` defaults to
`$HOME/cernbox/ana/pPb8160/pythia/cp5tune/cp5_gen_dijet_histograms.root`.
The `CP5_GEN_FILE` environment variable overrides every source default.
Combined outputs receive an `_embedding_cp5`
suffix. At
the standard eta-cut index 5, the notebook projects pseudorapidity shapes in the
same 16 trigger-specific half-open pTave intervals used by
`06_data_dijet_distributions.ipynb`. `ETA_NORMALIZATION` selects unit
bin-content-sum (`integral`, the default) or true unit-area density
(`bin_width`) normalization, with matching y-axis labels and simultaneous
generator-referenced lower panels. Full-eta plots are produced in three groups:
MC Reco (MC Emb.) with MC Gen (MC Emb.) and Data; MC Gen (MC Emb.) with MC Gen
(PYTHIA CP5) and Data; and MC Gen (PYTHIA CP5) with Data. The same groups are
used for Forward/Backward ratios built from the unnormalized projections using
independent-error propagation. The three-curve generator comparison shows both
Data/MC Gen (MC Emb.) and Data/MC Gen (PYTHIA CP5) in its lower panel. Every
lower panel is titled `Ratio to Gen`. Separate 1200 x 1200 summary
canvases collect all 16 intervals for every comparison group. Their annotations
and legends are placed across the top inside each plot frame, with vertical
headroom reserved so they do not overlap markers.
Their 1200-pixel PNG outputs provide at least 200 DPI up to a 6-inch figure,
and trigger names are omitted from the summary legends.
F/B construction and the subsequent comparisons use independent errors. The
default ratio ranges are 0.75--1.25 for pTave, 0.75--1.15 for the full-eta
ratios to generator MC, 0.75--1.25 for Forward/Backward, and 0.8--1.2 for the
subsequent Forward/Backward ratios to generator MC. All four y-axis ranges are
configured independently in `RATIO_Y_RANGES`. With `SAVE_RATIO_HISTOGRAMS`
enabled, the notebook writes every displayed ratio as an individually named
`TH1` in `dijet_data_vs_mc_ratios<comparison_suffix>.root`; configure the stem
with `RATIO_ROOT_BASENAME`. The input data
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

Use `notebooks/05_unfold2D_data.ipynb` to unfold the reconstructed data in the
same 16 trigger-specific pTave intervals used by
`06_data_dijet_distributions.ipynb`. MinimumBias supplies 50--80 GeV, Jet60
80--100 GeV, Jet80 100--120 GeV, and Jet100 120--500 GeV. The combined-direction
embedding sample supplies the flattened two-dimensional response,
reconstructed-bin purity, and truth-bin efficiency, and all unfolded inputs are
the corresponding combined-direction data files.
The factorized calculation applies purity, unfolds matched migrations with
RooUnfoldBayes, and then corrects for efficiency; 12 iterations are the
default. Before any pT projection or unfolding, every trigger TH2 is multiplied
by its cumulative MinimumBias-scale weight from the provenance-compatible JSON
written by `06_data_trigger_stitching.ipynb`; ROOT scales the stored `Sumw2`
errors at the same time. The notebook requires exactly one JSON matching the
configured data direction, jet selection, dijet CM observable, eta projection,
and dijet eta-cut index. Set `TRIGGER_STITCHING_WEIGHTS_FILE` to select an
explicit artifact when unfolding one direction. Stitching-weight uncertainties
are correlated normalization uncertainties: they are retained in output
metadata and are not added independently to each eta-bin statistical variance.
Output filenames include `data_stitched` so they cannot be confused with older
unweighted results. The workflow writes one ROOT result per direction with the
unfolded yield, transformed covariance, response diagnostics, trigger
provenance, and forward-folded reconstructed-level validation. Per-interval plots compare data
before and after unfolding in both unnormalized yield and independently
normalized eta shape; the normalized distributions are presentation-only
supplements to the unnormalized unfolded results.

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
ROOT's standard or `B` division option. Both ratio families have independently
configurable polynomial fits, currently `pol2` for full distributions and
`pol1` for F/B; the notebook retains the
coefficients, parameter errors, covariance matrix, chi-square, degrees of
freedom, and fit probability, and prints the scalar results for later
systematic-uncertainty extraction. Initial polynomial coefficients are
configured separately for Up/Def and Down/Def in both the full-distribution
and F/B fit families. `FIT_WEIGHT_OPTION = 'W'` uses ROOT's unit-weight fit for
every non-empty bin; set it to `''` to use the histogram bin uncertainties.
Set `SHOW_FIT_RESULTS = False` to hide the compact parameter and
chi-square/ndf box while retaining the fitted curves and numeric summaries.
Following `macro/systematics.C::calculateSystUncrtBinByBin` on `main`, the
notebook also calculates the unsmoothed relative JER uncertainty in each bin as
`(|Up/Def - 1| + |Down/Def - 1|) / 2`, evaluating the fitted Up/Def and Down/Def
functions at the bin center as in the macro's `useFit` path. Set
`SYSTEMATIC_COMBINATION = 'average'` to use the mean instead. The default is
`'maximum'`, namely `max(|Up/Def - 1|, |Down/Def - 1|)`.
Before the optional
smoothing step, the retained systematic histograms remain fractional; display clones
are scaled by 100 and plotted with the y-axis title
`JER Rel. Syst. Uncrt. (%)`, matching the macro's percent-plot convention.
Systematic bins that overlap the selected eta acceptance are included, even
when rebinning places their centers beyond the cut. For such boundary bins the
fit is evaluated at `ETA_CUT` (or `-ETA_CUT`) so it is never evaluated beyond
the accepted range; only bins fully outside the acceptance remain zero. Fits
also include every overlapping boundary bin by expanding ROOT's numerical fit
range to that bin's center, while the displayed fit curve remains restricted to
the physical eta acceptance.
Set `APPLY_SYSTEMATIC_SMOOTHING = True` to smooth both JER uncertainty families.
The F/B relative uncertainty becomes a running maximum from low to high
absolute eta. Full-CM smoothing runs outward independently in both directions
from the ROOT bin selected by `FULL_SMOOTHING_ORIGIN`, whose default is
`-0.465 + 0.00001` for the CM boost. Original unsmoothed histograms are retained.
Systematic-plot filenames contain `systematic_smoothed_` or
`systematic_nonsmoothed_` according to this option. All generated PDF and CSV
filenames also record the full/F/B fit functions, fit weighting, and systematic
combination before the final eta and pTave fields, for example
`fullFit_pol4_fbFit_pol2_weightsStd_systCombMax_etaCM_19_ptave_60_100`.

Use `notebooks/04_systematics_JEU.ipynb` to estimate the reconstructed-dijet
JEU shape uncertainty in MinimumBias, Jet60, Jet80, and Jet100 data. Each
trigger has its own configured pTave intervals. The notebook compares
unit-integral JEU Up, Down, and default CM eta projections and their
unnormalized Forward/Backward ratios, then draws dedicated Up/Def and Down/Def
comparisons. Forward/Backward construction always uses independent errors; the
later comparison-ratio error options are configured separately. The Up/Def and
Down/Def curves use independently configurable `pol2` full-CM and `pol1` F/B
fits, including configurable initial parameters and optional compact fit-result
text. `FIT_WEIGHT_OPTION = 'W'` gives every non-empty bin unit weight in the
ROOT fits; set it to `''` to use the histogram bin uncertainties. The relative
JEU systematic is evaluated from the fitted curves as
`(|Up/Def - 1| + |Down/Def - 1|) / 2` by default; set
`SYSTEMATIC_COMBINATION = 'maximum'` to use the larger absolute deviation.
The result is retained fractionally and plotted in
percent as `JEU Rel. Syst. Uncrt. (%)`. Bins overlapping the eta acceptance are
included; boundary-bin fit evaluation is clamped to `ETA_CUT`, while fully
outside bins remain zero. The fits themselves also include these overlapping
boundary bins, while displayed fit curves remain restricted to the acceptance.
`APPLY_SYSTEMATIC_SMOOTHING` optionally applies a low-to-high running maximum
for F/B and independent outward smoothing for full CM from the bin selected by
`FULL_SMOOTHING_ORIGIN` (default `-0.465 + 0.00001`). Both smoothed and original
unsmoothed objects are retained, and systematic filenames identify the selected
mode. Set
`DATA_DIRECTION` and `DATA_SELECTION` to choose the shared four-trigger data
production, and `DIJET_JEU_SYSTEMATICS_OUTPUT_DIR` to redirect its generated
plots.

The JEU output filenames use the same configuration tag convention as JER.
Pileup and pointing-resolution outputs use `systCombMax` because their
uncertainties are constructed from one variation.

Use `notebooks/04_systematics_pileup.ipynb` to estimate the reconstructed-dijet
pileup-filter shape uncertainty in MinimumBias, Jet60, Jet80, and Jet100 data.
The notebook overlays unit-integral `dz1p0`, `Gplus`, and optionally `Vtx1` full
CM eta projections and their unnormalized Forward/Backward ratios. It draws
`Gplus/dz1p0` and optional diagnostic `Vtx1/dz1p0` comparisons; only
`|Gplus/dz1p0 - 1|` defines the symmetric, double-sided uncertainty. Set
`SYSTEMATIC_EXTRACTION` to `fit` or `bin_by_bin`; outward running-maximum
smoothing is enabled by default. For fitted extraction,
`FIT_WEIGHT_OPTION = 'W'` gives every non-empty bin unit weight; set it to `''`
to use histogram bin uncertainties. Full-distribution and later F/B-to-F/B ratios
independently support standard or binomial errors, while F/B construction is
always performed with standard independent-error propagation.

Use `notebooks/04_systematics_pointing_resolution.ipynb` to estimate the dijet
pointing-resolution uncertainty in embedding or Pythia MC. It projects the Y
axis of nominal `hRecoDijetPtEtaCM_<eta cut>` and pointing-reference
`hRefDijetPtEtaCMPointingRes_<eta cut>` in common half-open reconstructed-pTave
intervals. Full shapes are independently normalized to unit integral.
Each Forward/Backward ratio is folded from its unnormalized full projection,
starting Forward at `FindBin(0 + 0.001)` and reflecting the bins below zero,
with independent errors. The one-sided uncertainty is
`|Reco/Ref - 1|`, evaluated from `pol2` full-shape and `pol1` F/B fits by
default or directly bin by bin. Outward running-maximum smoothing is optional
and enabled by default. For fitted extraction, `FIT_WEIGHT_OPTION = 'W'` gives
every non-empty bin unit weight; set it to `''` to use histogram bin
uncertainties. Full-shape Reco/Ref and the later F/B-to-F/B comparison
independently allow standard or ROOT option `B` errors.

The same notebook also produces a full-distribution-only Lab-frame estimator
from `hRecoDijetPtEtaLab_<lab eta cut>` and
`hRefDijetPtEtaLabPointingRes_<lab eta cut>`. `LAB_ETA_CUT_INDEX` is configured
separately from the CM selection and defaults to index 6, corresponding to
`|eta_Lab^jet| < 2.3`. It reuses the configured pTave intervals, Reco/Ref fit,
one-sided uncertainty, ratio-band, CSV, and outward smoothing workflows, with
the Lab smoothing origin fixed at zero. No Lab Forward/Backward ratio is made.

JER, JEU, pileup, and pointing-resolution write the evaluated full-CM and F/B
relative uncertainties as headerless CSV files. The columns are eta-bin center,
unity, eta half-width, and fractional systematic uncertainty, respectively, so
they can be loaded directly with
`ROOT.TGraphErrors(path, "%lg,%lg,%lg,%lg")`. Bins overlapping the configured
eta acceptance are written (including a boundary bin cut through by rebinning),
while bins entirely beyond the acceptance are omitted. Beam orientation only
produces comparison plots and therefore has no extracted uncertainty to export.

Use `notebooks/05_unfold2D_gen_and_prior.ipynb` to compare nominal Gen dijet
`etaCM` projections with the mixed-Gaussian prior
`hGenDijetPtEtaCMMixedPrior_<eta cut>` in every configured half-open interval
from the selected shared pTave bin set. The prior is an equal mixture of a
broad left Gaussian (`mean=-0.5`, `sigma=0.96`) and a narrow right Gaussian
(`mean=0.7`, `sigma=0.64`), formed by reweighting the nominal
`mean=0`, `sigma=0.8` Gen shape. It writes raw weighted-yield and
independently normalized shape overlays with `Mixed prior / Gen` lower panels,
using the shared file resolvers, annotations, ROOT styles, and output helpers.
Set `UNFOLD2D_GEN_PRIOR_OUTPUT_DIR` to redirect PDFs from
`output/unfold2D_gen_and_prior/`. Set `PLOT_RAW_DISTRIBUTIONS = False` to skip
the unnormalized yield canvases while retaining the normalized shape plots, and
set `SAVE_PNG = True` for PNG copies.

All `05_unfold2D*` notebooks expose `PTAVE_BIN_SET`. Select `test` to use
`config.histograms.TEST_DIJET_PTAVE_BINS` or `standard` to use
`config.histograms.DIJET_PTAVE_BINS`; both are interpreted as ordered half-open
intervals. Shared projection, flattening/unflattening, sparse-response assembly,
miss/fake diagnostics, RooUnfold construction, Bayesian execution, ROOT output,
and unfolding plots live in `python/unfolding.py` and
`python/unfolding_plots.py`.

Use `notebooks/05_unfold2D_basics.ipynb` to construct a flattened dijet
`pTave`-eta response and run a Bayesian RooUnfold closure test for the configured
MC sample using the eta-dependent JER-default measured distribution and its
corresponding response, miss, fake, and pair-classification objects. Misses are
handled with an explicit factorization: the measured distribution is multiplied
bin by bin by the MC matched-reco/inclusive-reco purity, the signal-only result
is unfolded with the matched migration response, and the unfolded matched-truth
distribution is divided bin by bin by the MC matched-truth/inclusive-truth
efficiency. The notebook retains the conventional RooUnfold fake/miss treatment
as a validation reference and writes a direct comparison of the two methods.
This factorization keeps all pTave-eta migrations and is valid when simulation
describes the reco-bin purity and truth-bin efficiency. Zero-efficiency truth
bins stop the workflow because they cannot be recovered. Purity and efficiency
are treated as response inputs; their finite-MC uncertainty must be assessed
with response variations or toys rather than added again during each bin-wise
operation. The notebook writes the
flattened input spectra, response and miss/fake diagnostics, unfolded result,
closure ratios, covariance matrix, serialized response, and configuration
metadata to
`output/unfold2D/<generator>_<direction>_unfold2D_jerDefExtra_eta_<cut>_iter_<iterations>.root`.
Set `DIJET_UNFOLD2D_OUTPUT_DIR` to redirect these generated artifacts without
editing the notebook configuration.
Set `COMPARISON_TARGET` to `gen` or `ref` to choose the denominator and target
shown in the flattened and per-pTave closure plots. This is a comparison-only
choice: the RooUnfold response, training truth, inefficiency, and Bayesian prior
remain Gen-based. The Ref option loads and projects `hRefDijetPtEtaCM_<eta cut>`.
Set `PLOT_MISS_AND_FAKES` to include or suppress Miss and Fake curves and their
legend entries in the projection and flattened diagnostic plots; it does not
disable their use in constructing or validating the correction factors and the
reference unfolding response.
Because the basics notebook uses its response-training sample as pseudo-data,
the factorized reco input, response prior, and Gen target are algebraically
consistent. Its unfolded/Gen and refolded/reco ratios should therefore be unity
to numerical precision. This is a technical identity check, not an independent
closure test.
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

Use `notebooks/05_unfold2D_studied_to_controlled.ipynb` for the independent
subsample closure test. It constructs the response, inefficiency, and fake
contribution exclusively from the reproducible 60% controlled MC subsample,
then unfolds the JER-default reconstructed distribution from the complementary
40% studied subsample and compares it with studied Gen truth. It applies the
controlled matched-reco/inclusive-reco purity to studied reco, unfolds with the
controlled matched-event migration response, and divides by the controlled
matched-truth/inclusive-truth efficiency. The conventional RooUnfold fake/miss
treatment is retained as a validation reference. The producer keys
`hRecoDijetPtEtaCMJerDefUnfoldControlledSample_<eta cut>` and
`hRecoDijetPtEtaCMJerDefUnfoldStudiedSample_<eta cut>` use the same data-driven
JER-smeared reconstructed coordinate as the controlled response matrix. Configure
`PLOT_MISS_AND_FAKES` to show or hide controlled Miss/Fake diagnostics and use
`DIJET_UNFOLD2D_STUDIED_CONTROLLED_OUTPUT_DIR` to redirect output from
`output/unfold2D_studied_to_controlled/`.
The reconstructed AK4 jets use the data-driven default JER smearing. The
flattened space is restricted to the
configured eta acceptance on both response axes and all spectra. Standard
contiguous pTave bins are the default; the gapped test bins are only a quick
execution check. Before unfolding, a separate accounting check verifies that
each controlled truth/reco marginal equals its matched response projection plus
the effective miss/fake component. A separate diagnostic canvas distinguishes
explicit selection/matching failures, pTave-boundary migrations, and the
effective response components. The response, purity, efficiency, misses, fakes,
and their diagnostics remain controlled-sample-only; studied objects supply
only the independent measured input and closure truth. The configured Bayesian
iteration count remains fixed for this workflow. The flattened controlled response matrix
is also saved by itself as `<output tag>_response_matrix.pdf` with a logarithmic
color scale; the matrix histogram remains stored in the ROOT output.
Eta-range restriction uses inclusive bin centers. Consequently, after eta
rebinning a bin spanning 1.8--2.0 (center 1.9) is retained for an acceptance
ending at 1.9, consistently for projections and both response axes.
Before unfolding, the notebook also saves one two-panel comparison per pTave
interval: controlled versus studied Gen shapes and controlled versus studied
reco shapes. Each displayed histogram is cloned and independently normalized
as `1/N dN/detaCM` using its bin widths, leaving all response and unfolding
inputs unchanged. Legends explicitly identify controlled/studied Gen or Reco. Lower
panels on the same canvas show the studied/controlled ratios with independently
propagated histogram uncertainties. The normalized comparison histograms and
ratios are stored in the notebook ROOT output.
The independent studied closure is not expected to be exactly unity in every
bin. Controlled and studied samples contain different weighted-event
fluctuations. Iterative Bayesian regularization also produces a bias--variance
tradeoff: increasing iterations generally improves refolded/reco agreement
while amplifying fluctuations and potentially worsening unfolded/studied-Gen
agreement. Both notebooks print an integral ratio, a yield-weighted L1
difference over all bins, and a mean relative difference over populated bins.
These summaries prevent nearly empty bins with unstable ratios from dominating
the closure interpretation.
The notebook also forward-folds the unfolded studied spectrum through the
controlled response and compares it with studied reco in flattened global bins
and in every pTave interval. Miss fractions persist through the controlled
response efficiency. Fake fractions are transferred bin by bin by multiplying
the folded matched-reco prediction by the controlled fake/matched-reco ratio;
no controlled/studied sample-size normalization is applied. The folded signal,
transferred fake component, total forward-folded spectrum, and
forward-folded/reco ratios are stored in the ROOT output; the fake-handling
decision is recorded in its metadata. `RECO_TO_GEN_RATIO_OPTION`,
`UNFOLDED_TO_GEN_RATIO_OPTION`, and `FORWARD_FOLDED_TO_RECO_RATIO_OPTION`
independently select ROOT standard (`''`) or binomial (`'B'`) displayed ratio
errors for flattened and per-pTave comparisons; the selected options are saved
in the ROOT metadata.

The recommended unfolding validation order is:

1. Run `05_unfold2D_basics.ipynb` as a same-sample technical closure check.
2. Run `05_unfold2D_studied_to_controlled.ipynb` as the independent 60/40 MC
   closure test.
3. Inspect the prior deformation with `05_unfold2D_gen_and_prior.ipynb`.
4. Run `05_unfold2D_regularization.ipynb` to study the Bayesian iteration count.

Use `05_unfold2D_regularization_tutorial.ipynb` when teaching or auditing that
last step. It exposes every important intermediate distribution and follows one
configured toy through purity correction, sequential RooUnfoldBayes iterations,
and efficiency correction before running a deliberately small sequential toy
ensemble. It never loads the compiled parallel regularization backend.

The regularization notebook uses the same factorized unfolding convention as
the closure notebooks and separates two truth roles explicitly.
`RESPONSE_GENERATOR`, `RESPONSE_DIRECTION`, and `RESPONSE_FILE_STEM` select the
MC that constructs the response matrix, purity, efficiency, and RooUnfold
Bayesian prior. `TOY_TRUTH_MODEL` independently chooses the known truth that is
forward-folded to generate pseudo-data: `mixed_prior`, `cp5`, or `gaussian`.
This choice does not replace RooUnfold's response-training Bayesian prior.
The CP5 option reads `CP5_TOY_TRUTH_HISTOGRAM` from `CP5_TOY_TRUTH_FILE`. The
Gaussian option generates only the eta shape; its pTave marginal is copied from
the configured response generator using `GAUSSIAN_PT_SOURCE_ETA_CUT_INDEX`.
Output and cache names explicitly contain both `response-...` and
`toytruth-...` tags, and figures repeat both roles in their annotations. The
notebook forward-folds the selected toy truth with the inclusive response and
restores fakes through the binwise response-MC fake/matched-reco ratio. It does
one global exposure normalization before constructing the toys. The direct
MinimumBias--Jet60 overlap range is read from the validated trigger-stitching
artifact. The ratio of the MinimumBias data yield to the unscaled forward-folded
toy yield in that range multiplies the complete toy truth once. The notebook
then forward-folds the scaled truth; it does not tune individual pTave blocks.
The original truth, scaled truth, reference projection, two normalization
integrals, and global scale factor are stored in the ROOT output and metadata.
The compact `test` binning cannot resolve the 90--140 GeV overlap boundary; for
that code-check mode only, normalization falls back to its complete
MinimumBias-assigned interval and records that provenance explicitly.

### Statistical-uncertainty propagation

Trigger-stitched experimental distributions contribute their **relative per-bin
statistical uncertainties**. For reconstructed bin `i` assigned to trigger `T`,
let the unstitched data content and statistical error be `d_i` and `delta_i`,
and let the cumulative stitching weight be `w_T`. ROOT scaling gives stitched
content `D_i = w_T d_i` and error `Delta_i = w_T delta_i`, so the fractional
statistical uncertainty is

```text
r_i = Delta_i / abs(D_i) = delta_i / abs(d_i).
```

The common stitching weight cancels in the ratio. The fractional uncertainty is
then transferred to the independently normalized forward-folded toy mean:

```text
content_i = mu_i = forward_fold(toy_truth)_i
error_i   = sigma_i = abs(mu_i) r_i.
```

The stitched content `D_i` supplies the single global MinimumBias exposure
factor described above but does not reshape individual pTave blocks. With the configured
`EMPTY_DATA_UNCERTAINTY_POLICY = 'nearest_in_pt_block'`, a populated MC bin with
no usable data uncertainty copies `r_i` from the nearest populated data bin in
the same pTave block. Ties select the lower global eta bin, and the search never
crosses a pTave boundary. The workflow raises if an entire pTave block has no
data bin with positive content and positive statistical error. A source content
no larger than `TOY_EXPECTATION_ZERO_RELATIVE_TOLERANCE` times the largest
source bin is treated as numerical zero and set exactly to zero together with
its error.

A seeded common ensemble of `N_TOYS` inclusive-reco histograms is generated
with `x_i ~ Normal(mu_i, sigma_i)` for bins with nonzero original expectation.
Originally empty bins remain exactly zero. For populated bins, non-positive
Gaussian draws are rejected and redrawn until a positive dijet yield is
obtained. Only the bin contents fluctuate, while every toy retains the data
error `sigma_i`. Each toy is multiplied by full-MC purity,
unfolded with the matched-only response and its response-training Bayesian
prior, and divided by response-MC efficiency before its bias and variance are
evaluated against the known toy truth.
Before those errors are copied, every experimental trigger TH2 is multiplied
by its cumulative MinimumBias-scale weight from the provenance-compatible JSON
written by `06_data_trigger_stitching.ipynb`; ROOT scales its stored `Sumw2`
errors at the same time. The notebook requires exactly one JSON matching the
configured data direction, jet selection, dijet CM observable, eta projection,
and dijet eta-cut index. Set `TRIGGER_STITCHING_WEIGHTS_FILE` to select an
explicit artifact when needed. Stitching-weight uncertainties are correlated
normalization uncertainties and are recorded in metadata rather than added as
independent eta-bin statistical errors. The artifact path, tag, cumulative
weights, and histogram keys participate in cache validation, so an unstitched
or differently stitched toy cache is not reused.
The default `RESAMPLE_NEGATIVE_TOY_BIN = True` rejects non-positive Gaussian
draws in populated bins and redraws until each value is positive. This is a
lower-truncated Gaussian: for bins with small `mu_i/sigma_i`, the realized toy
mean and variance need not equal the nominal `mu_i` and `sigma_i^2`. Setting it
to `False` restores the strict bin-level stop controlled by
`STOP_ON_NEGATIVE_TOY_BIN`. `PTAVE_BIN_SET = 'test'` is a
fast workflow check, while `standard` is the intended analysis binning.
All response-MC, toy-truth, sparse-response, and data projections use the same configured
eta acceptance, expanded to complete histogram-bin edges when necessary. The
broad test intervals use MinimumBias for 50--100 GeV, Jet80 for 100--180 GeV,
and Jet100 above 180 GeV; the standard intervals retain the finer established
trigger transitions.

Purity, response, and efficiency are fixed in this statistical ensemble. Their
own uncertainties, the uncertainty of the stitching weights, and other
systematic uncertainties are not sampled. Nevertheless, the data statistical
fluctuation in each toy is numerically propagated through the fixed purity
correction, Bayesian migration, and efficiency correction before the ensemble
moments are calculated.

For every iteration, bin `i` uses
`bias_i = mean(unfolded_i) - toy_truth_i`, unbiased toy variance with the `N_TOYS-1`
denominator, and `MSE_i = bias_i^2 + variance_i`. The notebook plots averages
of these per-bin quantities in both absolute and truth-relative form; it never
**Mean metric** is `(1/N_bins) sum_i q_i` and retains yield-squared units. The
**Mean fractional metric** is `(1/N_valid) sum_i q_i/toy_truth_i^2`; it is
dimensionless and gives each valid bin equal fractional weight rather than
applying one global normalization to the absolute mean. Equal eta-bin widths
need no additional correction. `REGULARIZATION_SELECTION_MODE = 'absolute'`
is the default, so the recommended iteration preserves the intended larger
contribution from bins with larger absolute discrepancies. Set it to `relative`
when typical fractional per-bin performance is the desired objective, or to
`both` to report and store both minima. Because later closure and output code
requires one iteration, `both` uses the absolute minimum as that primary result.
Bins with zero or negligible toy truth
are reported and excluded from relative summaries. The first minimum of the global
configured MSE is reported as the recommended iteration count.

For the matched-only migration matrix, let `R_ij = P(reco i | truth j)` and let
`pi_j^(0)` be the normalized response-training truth prior. For toy `k`, Bayesian
step `n` calculates

```text
d_i^(n)   = sum_l R_il pi_l^(n-1)
B_ji^(n)  = R_ij pi_j^(n-1) / d_i^(n)
u_j^(k,n) = sum_i B_ji^(n) m_i^(k)
pi_j^(n)  = u_j^(k,n) / sum_l u_l^(k,n).
```

Zero-denominator reconstructed bins contribute zero. The unfolded matched-truth
content `u_j` is divided by the separately derived inclusive-truth efficiency
before entering the toy moments. The compiled backend performs this recurrence
from step 1 through `MAX_ITERATIONS` once per toy. At every step it accumulates
`S_jn = sum_k u_j^(k,n)` and `Q_jn = sum_k (u_j^(k,n))^2`; the mean and unbiased
variance follow from `S/K` and `(Q - S^2/K)/(K-1)`. This is algebraically the
same content calculation as the validated RooUnfold path, but avoids restarting
the first `n-1` updates when requesting result `n` and never constructs unused
per-toy covariance matrices.
Each worker thread accumulates partial `S` and `Q` arrays for its disjoint toy
chunk. These raw partial moments are summed over all workers first, in a fixed
worker order. Only then are the mean, unbiased variance, bias squared, and MSE
calculated using the total `K = N_TOYS`. The code does not calculate per-worker
metrics or average worker variances, which would give an incorrect result when
chunk sizes differ and would omit between-chunk mean fluctuations.

A minimum at iteration 1 or `MAX_ITERATIONS` is a boundary result, not evidence
that the physical optimum has been located. In particular, if the selected
minimum equals `MAX_ITERATIONS`, the configured scan cannot distinguish a true
minimum there from an MSE curve that is still falling or will turn later. Keep
the same seed, toys, response, truth, binning, and selection metric, increase
`MAX_ITERATIONS`, and rerun until the preferred minimum is interior with scan
points on its higher-iteration side. If it remains at the new upper boundary,
extend again. Changing any other configuration changes the statistical study
and is not a pure extension. This recommendation is conditional on the selected
response, toy truth, binning, and data uncertainty model.
`LOG_REGULARIZATION_METRICS_Y = True` uses logarithmic y axes
for the global and optional per-pTave Bias-squared, Variance, and MSE comparison
canvases; set it to `False` for linear axes. The legacy backend performs
`N_TOYS * MAX_ITERATIONS` separate RooUnfold calls. The default compiled C++17
backend instead advances every toy through the full iteration range once and
parallelizes toy chunks over `REGULARIZATION_SCAN_WORKERS`; zero selects the
available hardware threads. Build it with `cmake -S hist_analysis/cpp -B
build/regularization-scan -DCMAKE_BUILD_TYPE=Release` followed by `cmake --build
build/regularization-scan --parallel`. The bridge also accepts
`JETANALYSIS_REGULARIZATION_SCAN_LIBRARY` for a nonstandard library location.
Both backends skip per-toy covariance construction by default because those
matrices do not enter ensemble bias or variance. The selected final unfolding
still calculates and stores its full covariance. Set
`CALCULATE_TOY_COVARIANCE = True` in the notebook only when per-toy covariance
construction is explicitly needed by the legacy backend; it substantially
increases runtime. The compiled backend returns metric moments rather than all
unfolded histograms, so `STORE_ITERATION_DISTRIBUTIONS` must be `False`; select
the legacy backend for the optional complete unfolded-toy archive. Outputs are written
below `output/unfold2D_regularization/`; set
`DIJET_UNFOLD2D_REGULARIZATION_OUTPUT_DIR` to redirect them.
Every generated filename includes `_toyseed_<TOY_ENSEMBLE_RANDOM_SEED>`. With
`USE_INTERMEDIATE_CACHE = True`, the inclusive-reco toy distributions and all
per-iteration mean, bias, bias-squared, variance, and MSE distributions are
stored as compact ROOT caches beneath `output/unfold2D_regularization/intermediate/`.
With the legacy backend and `STORE_ITERATION_DISTRIBUTIONS = True`, a separate,
potentially large ROOT file also records the eta projections, flattened training inputs,
response matrix and accounting diagnostics, purity and efficiency, forward-fold
components, inclusive and purity-corrected toys, every efficiency-corrected
unfolded toy grouped by iteration, and every metric distribution. This is the
complete reproducibility input for the metric calculation; eta slices of a
flattened result can be recovered from its documented pT-block layout.
Cache metadata includes input-histogram fingerprints, binning, correction
factors, toy settings, iteration range, and seed. Compatible caches are loaded
on rerun, skipping toy generation and completed RooUnfold scans. An atomic scan
checkpoint is replaced after every completed Bayesian-iteration scan point. If
execution stops, the next run resumes at the following iteration. When complete
iteration distributions are enabled, their temporary ROOT archive is flushed
at the same boundary and must agree with the metric checkpoint before that
prefix is reused. Mismatched caches are regenerated, and a partial iteration is
never advertised as complete.
Because `MAX_ITERATIONS` participates in cache metadata and filenames, extending
the range creates a distinct scan cache rather than silently treating a shorter
boundary-limited result as complete.
The scan always unfolds and evaluates the complete flattened distribution.
Optional per-pT metric figures are diagnostic slices of the global-bin results
and never select independent iteration counts.
The inspection cell below the scan accepts a zero-based toy model number, a
tuple of Bayesian iterations, and a stored random seed. It reads the complete
distribution ROOT file without unfolding again, reconstructs each pTave eta
distribution, and writes separate target-truth/unfolded overlays and
unfolded/target-truth ratio canvases.
A second inspection cell accepts a seed and a tuple of zero-based toy numbers.
It compares the stored inclusive-reco toy fluctuations with their common
forward-folded expectation and draws the toy/expectation ratio for every pTave
eta interval, again without generating or unfolding toys.

Use `notebooks/05_unfold2D_mc_direction.ipynb` for the beam-direction closure
test. It constructs the response and its training truth/reco marginals, misses,
fakes, and pair classification from Pb-going embedding, then unfolds the
independent p-going embedding eta-dependent JER-default reco distribution and
compares it with p-going generator truth. The notebook uses the same flattened
bin mapping, RooUnfold response scaling, explicit fake handling, covariance,
plotting styles, and ROOT output objects as `05_unfold2D_basics.ipynb`. Its output tag
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
spectra are normalized in the configured 110--130 GeV interval. The dijet
notebook uses bin-width normalization for true unit-area `1/N dN/deta`
densities and also folds each raw CM eta projection into an unnormalized
Forward/Backward ratio, overlays the three orientations, and plots the
Pb-going/combined and p-going/combined F/B ratios with independent-error
propagation. Its eta projections are displayed through the configured eta cut
plus 0.1, and a final section collects the 16 selected intervals into two
legend-free square 4x4 figures: one for combined-orientation normalized CM
shapes and one for combined-orientation F/B curves. Inclusive-jet
and dijet projection intervals are independently configurable for MinimumBias,
Jet60, Jet80, and Jet100; inclusive-jet pT spectra also have a configurable eta
selection. In the inclusive-jet notebook, set `ETA_DISPLAY_RANGE` to limit the
eta axis and `ETA_Y_RANGE` to set a common y-axis range for the normalized eta
orientation overlays; leave either as `None` for the full axis or automatic
scaling. For every trigger and pT interval, the notebook also overlays the
Pb-going and p-going unflipped-Lab corrected/raw eta-yield ratios. These use
independent-error propagation and are displayed within
`CORRECTED_RAW_RATIO_RANGE`; a lower panel compares Pb-going to p-going within
`CORRECTED_RAW_DIRECTION_RATIO_RANGE`. Set `PPB_DATA_DIR`,
`DATA_JET_OUTPUT_DIR`, or `DATA_DIJET_OUTPUT_DIR` to override input and output
locations.

Use `notebooks/06_data_trigger_stitching.ipynb` to extract multiplicative
MinimumBias-scale weights for the MinimumBias, Jet60, Jet80, and Jet100 data
streams. The default projects `hRecoInclusiveJetPtEtaLab` over the half-open
range `-1.6 <= eta_Lab < 1.6`; `OBSERVABLE` and `FRAME` select inclusive jets
or dijets in Lab or CM coordinates. For dijets, `DIJET_ETA_CUT_INDEX` selects
the production-level single-jet acceptance, while `PROJECTION_ETA_RANGE`
selects the dijet-axis interval used for the pTave projection. The notebook
matches Jet60 to MB in 90--140 GeV, Jet80 to Jet60 in 110--180 GeV, and Jet100
to Jet80 in 130--200 GeV using ratios of unnormalized overlap integrals with
independent statistical-error propagation. It writes raw and rescaled adjacent
trigger comparisons, a final four-trigger overlay with ratios to MB, and the
local and cumulative weights in JSON and CSV beneath
`output/data_trigger_stitching/`. Every artifact filename records the data
direction, jet selection, observable, frame, projected eta range, applicable
dijet jet-eta cut, and all three trigger-matching ranges. The CSV repeats that
provenance in flat columns alongside the histogram keys, overlap integrals,
local factors, and cumulative MB-scale weights. Set
`DATA_TRIGGER_STITCHING_OUTPUT_DIR` to redirect those outputs.

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
