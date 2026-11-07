"""Exact supplied theory artifacts shared by production and comparisons."""
from pathlib import Path

THEORY_FILENAMES = {
    'epps21_ct18': 'epps21_pPb_pp_pdf_scale_uncertainties_eta1p9_symmetric_cl68_fitresidual.root',
    'ncteq15hq': 'ncteq15hq_pPb_pp_pdf_scale_uncertainties_eta1p9_symmetric_cl68_fitresidual.root',
}


def theory_file(model, project_root=None):
    root = Path(project_root) if project_root is not None else Path(__file__).resolve().parents[2]
    return root / 'hist_analysis/root_files' / THEORY_FILENAMES[model]
