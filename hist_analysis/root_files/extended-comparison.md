# EPPS21 × CT18 campaign comparison

`extended-comparison.root` holds central predictions at 8.16 TeV for Nov. 09 (15,001,049,074 saved samples), Nov. 18 partial (15,000,050,152), and their native pool (30,001,099,226). The Nov. 18 merge omits two jobs.

There are 24 ROOT directories, e.g. `pt050_060_eta1p9`, spanning seven production pTave bins plus 50–70 GeV, each at three jet acceptances (1.5, 1.9, 2.4). Observable prefixes are `ratio_central` (RpPb), `forward_backward_pPb_over_pp` (F/B of RpPb), `pdf_central` (pp), `npdf_central` (pPb), `pp_forward_backward`, and `pPb_forward_backward`.

For each prefix, `_<nov09|nov18|combined>_propagated` stores the central histogram with pre-fit statistical errors. `_<nov09|nov18>_over_combined` stores the campaign-to-pool ratio with shared-sample statistical covariance included. Invalid bins are NaN. The top-level `statistical_definition` describes the covariance approximation. PDF and scale bands, fit residuals, and individual PDF members live in separate analysis outputs.
