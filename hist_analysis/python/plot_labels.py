"""Shared experiment-header labels for histogram figures."""

CMS_SIMULATION_HEADER = "CMS Simulation"
CMS_PRELIMINARY_HEADER = "CMS Preliminary"
CMS_HEADER = "CMS"
CMS_HEADER_CHOICES = (
    CMS_SIMULATION_HEADER,
    CMS_PRELIMINARY_HEADER,
    CMS_HEADER,
)

PP_HEADER_RIGHT = "pp (8.16 TeV)"
PPB_HEADER_RIGHT = "pPb 174.6 nb^{-1} (8.16 TeV)"
PPB_MPLHEP_HEADER_RIGHT = r"pPb 174.6 nb$^{-1}$ (8.16 TeV)"


_BEAM_ORIENTATION_LABELS = {
    'pgoing': 'p-going',
    'Pbgoing': 'Pb-going',
    'combined': None,
    # Some comparison helpers already use formatted direction keys.
    'p-going': 'p-going',
    'Pb-going': 'Pb-going',
    # pp is a collision system, retained in workflows that distinguish it.
    'pp': 'pp',
}


def beam_orientation_label(direction: str | None) -> str | None:
    """Format a beam orientation for display only."""
    return None if direction is None else _BEAM_ORIENTATION_LABELS[direction]


def beam_sample_label(sample: str, direction: str | None, *, sep: str = ', ',
                      orientation_first: bool = False) -> str:
    """Keep the sample text and append (or prepend) any orientation."""
    orientation = beam_orientation_label(direction)
    if not orientation:
        return sample
    return sep.join((orientation, sample) if orientation_first else (sample, orientation))
