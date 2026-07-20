# hist_analysis

Workspace for ROOT-file inventory, event diagnostics, closure comparisons, and inclusive-jet JES/JER studies built around the direction-specific and combined `pPb8160` outputs.

## Layout

- `notebooks/`: Jupyter workflows for inventory, event diagnostics, closures, and JES/JER studies.
- `python/`: reusable PyROOT helpers for histogram I/O, projections, normalization, fitting, and plotting.
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
two-exponential threshold fits over 15–950 GeV.

Use `notebooks/03_basic_closures.ipynb` for single-jet and dijet overlays and
ratios to an explicitly selected nominal. It supports p-going, Pb-going, and
existing combined files; gen/ref/reco comparisons within one MC generator; and
embedding/Pythia comparisons at a selected reconstruction level. Explicit data
curves require a chosen input file, histogram key, trigger combination, and
normalization convention.

Use `notebooks/03_jet_JES_JER.ipynb` to extract JES and JER from Gaussian fits to
inclusive-jet response slices, compare the configured reco/gen smearing cases,
scan pT and eta dependence, and draw the corresponding two-dimensional response
maps. The notebook also normalizes JER-versus-eta curves to unit central-bin mean
within `-0.8 < eta < 0.8` and writes those derived histograms to a ROOT file in
`output/jes_jer/`.

The closure projection helpers identify pT and eta axes from their physical
ranges because some existing reference-jet files contain misleading axis
titles. Configured selection intervals are interpreted as half-open
`[low, high)` ranges.
