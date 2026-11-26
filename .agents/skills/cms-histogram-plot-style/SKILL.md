---
name: cms-histogram-plot-style
description: Create, modify, review, or validate CMS histogram figures in this project's notebooks and Python plotting helpers using ROOT or mplhep. Use for requests involving CMS headers, typography, annotations, legends, markers, error bars, canvas geometry, shared-axis grids, or consistent styling across histogram workflows.
---

# CMS Histogram Plot Style

Apply a consistent presentation to ROOT and mplhep histogram figures without changing histogram contents, selections, normalization, binning, or statistical calculations.

## Project defaults

- Use the shared labels in `hist_analysis/python/plot_labels.py`; do not duplicate luminosity or collision-energy strings.
- MC, systematic, and simulation-only unfolding notebooks default to `CMS Simulation` on the left and `pPb 174.6 nb^{-1} (8.16 TeV)` on the right.
- Any notebook whose filename contains `data` defaults to `CMS Preliminary` on the left and the same pPb luminosity and energy on the right. Keep the left label configurable so it can later switch to `CMS` without changing plotting code.
- `03_pythia_standalone_check.ipynb` defaults to `CMS Simulation` and `pp (8.16 TeV)`.
- Keep headers configurable as `CMS Simulation`, `CMS Preliminary`, or `CMS` on the left and as a math-capable collision-system label on the right.
- Draw the header above the frame. Keep descriptive annotations inside the frame without overlap.

## Shared visual conventions

- Reuse `hist_analysis/python/root_style.py`, `hist_analysis/python/plotting.py`, and `hist_analysis/python/mplhep_plotting.py`; do not add notebook-local styling when a shared option suffices.
- Keep inclusive-jet and dijet figures on the same typography, marker, legend, annotation, margin, and header configuration.
- Use proper ROOT `TLatex` or matplotlib math text for subscripts, superscripts, Greek symbols, and units.
- Make figure size, font sizes, marker sizes, error-bar cap size, annotation placement, legend placement, and axis ranges configurable.
- For tune comparisons, use a filled red circle for CP5 and an open blue circle for CUETP8M1. Retain deterministic distinct color/marker choices for overlays with additional curves.
- Use the shared eight-curve mplhep style cycle. It alternates filled and open circle, triangle, square, and diamond markers while retaining distinct colors; open markers use `markerfacecolor="none"` with a colored edge.
- On mplhep point plots, show vertical statistical uncertainties and horizontal half-bin-width uncertainties. Suppress points and both error bars where bin content is zero without modifying the source ROOT histogram.
- Preserve original ROOT figures when an mplhep export is added.
- Every canvas saved through `root_style.save_canvas` must support an optional parallel mplhep PDF, including 1D overlays, ratio pads, TH2 maps, and recursively nested multi-panel canvases. Write mirrors beneath a sibling `mplhep/` output directory.
- By default, copy each ROOT pad's final displayed x- and y-axis ranges to its mplhep counterpart, including conversion from ROOT log-pad coordinates. Preserve an explicitly configured mplhep range as an override.

## Shared-axis histogram grids

- For dijet closure plots, use square standalone canvases and square grid panels. Reuse the standalone curves, colors, marker fills, uncertainties, and reference lines in the grid. For closure comparisons, style each ratio by its numerator.
- For the dijet closure 4×4 summaries, fill panels in row-major order using the first 16 configured intervals. Make adjacent frames meet, share the x and y scales within each grid, and show tick labels only on the left column and bottom row. Use one x-axis and one y-axis title for the grid.
- Preserve configured limits, including ratio limits. When standalone panels use automatic limits, use a common range covering their displayed ranges. Check text against points and error bars at the configured limits.
- Use one CMS header above the grid with its left edge aligned to the leftmost plotting axis. Enlarge outer tick labels, axis titles, and in-panel text for legibility. Center the interval label at the top of each panel and keep the other annotations below it.
- In the dijet closure grid, put the generator and orientation in row 1, column 1 and the marker-only legend in row 1, column 2. Append `p-going` or `Pb-going` to the generator with a comma for directional samples. For combined samples, show the generator alone.
- Put the jet eta and delta-phi cuts on two lines in row 2, column 1. Put the leading and subleading jet pT cuts on two lines in row 2, column 2. Keep the physical cuts visible and omit the internal Reco selection identifier from the grid.

## Dijet forward/backward axis title

- For the dijet F/B observable, use `Forward (#eta_{CM}^{dijet} > 0) / Backward (#eta_{CM}^{dijet} < 0)` in ROOT TLatex and the equivalent math text in mplhep. Put the sign convention in that y-axis title and use a separate observable label for ratios of F/B curves.

## Scope and validation

- Keep physics-specific ranges, selections, and annotations in notebook configuration rather than freezing them in this skill.
- For Forward/Backward plots, also follow `root-forward-backward-ratio`; styling must not alter its unnormalized construction or independent-error propagation.
- Validate changed Python modules with the project interpreter, parse changed notebook JSON and code cells, and execute the smallest affected notebook workflows when inputs are available. Visually inspect representative ROOT and mplhep outputs when geometry changes.
