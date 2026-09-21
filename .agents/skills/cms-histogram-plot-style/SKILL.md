---
name: cms-histogram-plot-style
description: Create, modify, review, or validate CMS histogram figures in this project's notebooks and Python plotting helpers using ROOT or mplhep. Use for requests involving CMS headers, typography, annotations, legends, markers, error bars, canvas geometry, or consistent styling across histogram workflows.
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

## Scope and validation

- Keep physics-specific ranges, selections, and annotations in notebook configuration rather than freezing them in this skill.
- For Forward/Backward plots, also follow `root-forward-backward-ratio`; styling must not alter its unnormalized construction or independent-error propagation.
- Validate changed Python modules with the project interpreter, parse changed notebook JSON and code cells, and execute the smallest affected notebook workflows when inputs are available. Visually inspect representative ROOT and mplhep outputs when geometry changes.
