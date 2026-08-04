"""Run-by-run reconstructed-dijet CM pseudorapidity comparisons."""

from __future__ import annotations

from pathlib import Path

from IPython.display import display

from hist_analysis.python.data_distributions import project_data_histogram
from hist_analysis.python.histogram_ops import normalize_histogram, ratio_to_nominal
from hist_analysis.python.plotting import draw_overlay


# Keep this catalog synchronized with diagnosticRunIds and
# diagnosticRunPileup in processForestSimple.C, where the histograms are filled.
DIAGNOSTIC_RUN_IDS = (285480, 285505, 285517, 285832, 285993)
DIAGNOSTIC_RUN_PILEUP = {
    285480: 0.04,
    285505: 0.25,
    285517: 0.1,
    285832: 0.004,
    285993: 0.2,
}
RUN_INCLUSIVE_KEY = "hRecoDijetPtEtaCMRunIntegrated"


def _tag(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def draw_run_dependence(
    trigger: str,
    filename: str | Path,
    ptave_bins,
    *,
    output_dir: str | Path,
    rebin_eta: int,
    eta_range: tuple[float, float] | None,
    ratio_range: tuple[float, float],
    ratio_option: str,
    save_png: bool,
    grid: bool,
    style,
):
    """Draw normalized run overlays and separate ratios to the inclusive run set.

    ``ratio_option`` is used only for the later run/inclusive comparison.
    ROOT's independent (``''``) and binomial (``'B'``) modes are supported.
    """

    path = Path(filename)
    output = Path(output_dir)
    if ratio_option not in {"", "B"}:
        raise ValueError("ratio_option must be '' or 'B'")
    sample_tag = trigger.lower().replace(" ", "_")
    keys = {
        **{
            f"{run_id} (PU {DIAGNOSTIC_RUN_PILEUP[run_id]:.3g})":
                f"hRecoDijetPtEtaCMRun_{run_id}"
            for run_id in DIAGNOSTIC_RUN_IDS
        },
        "Run inclusive": RUN_INCLUSIVE_KEY,
    }
    style_indices = {label: index for index, label in enumerate(keys)}
    results = {}

    for low, high in ptave_bins:
        pt_tag = f"{_tag(low)}_{_tag(high)}"
        histograms = {}
        for curve_index, (label, key) in enumerate(keys.items()):
            projection = project_data_histogram(
                path, key, "eta", (low, high),
                name=f"h_{sample_tag}_run_curve_{curve_index}_{pt_tag}",
                rebin=rebin_eta,
            )
            # Compare eta shapes rather than the different run event totals.
            histograms[label] = normalize_histogram(projection, "integral")

        annotations = (
            "pPb 8.16 TeV data", trigger,
            f"{low:g} #leq p_{{T}}^{{ave}} < {high:g} GeV",
        )
        overlay_canvas = draw_overlay(
            histograms, title="", x_title="#eta_{CM}^{dijet}",
            y_title="1/N dN/d#eta_{CM}^{dijet}", x_range=eta_range,
            annotations=annotations, style_indices=style_indices,
            grid=grid, style=style, headroom=1.55,
            output=output / f"{sample_tag}_dijet_etaCM_run_overlay_ptave_{pt_tag}.pdf",
            save_png=save_png,
            canvas_name=f"{sample_tag}_dijet_etaCM_run_overlay_{pt_tag}",
        )
        display(overlay_canvas)

        # The integrated control is the common denominator for every run.
        inclusive = histograms["Run inclusive"]
        ratios = {
            label: ratio_to_nominal(
                histogram, inclusive,
                name=f"h_{sample_tag}_run_curve_{curve_index}_over_inclusive_{pt_tag}",
                option=ratio_option,
            )
            for curve_index, (label, histogram) in enumerate(histograms.items())
            if label != "Run inclusive"
        }
        ratio_canvas = draw_overlay(
            ratios, title="", x_title="#eta_{CM}^{dijet}",
            y_title="Run / run inclusive", x_range=eta_range,
            y_range=ratio_range, annotations=annotations,
            style_indices=style_indices, grid=grid, style=style,
            reference_y=1.0,
            output=output / f"{sample_tag}_dijet_etaCM_run_ratio_ptave_{pt_tag}.pdf",
            save_png=save_png,
            canvas_name=f"{sample_tag}_dijet_etaCM_run_ratio_{pt_tag}",
        )
        display(ratio_canvas)
        results[(low, high)] = {
            "histograms": histograms,
            "ratios_to_inclusive": ratios,
            "overlay_canvas": overlay_canvas,
            "ratio_canvas": ratio_canvas,
        }

    return results
