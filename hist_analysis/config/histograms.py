"""Shared closure selections and the standard dijet eta-cut index."""

from __future__ import annotations

SINGLE_JET_PT_BINS = [(50, 80), (180, 200), (400, 500)]
DIJET_PTAVE_BINS = [
    (50, 60), (60, 70), (70, 80), (80, 90), (90, 100),
    (100, 110), (110, 120), (120, 130), (130, 140), (140, 150),
    (150, 160), (160, 180), (180, 200), (200, 250), (250, 300),
    (300, 500), (500, 1000),
]
COMMON_ETA_CM_RANGE = (-1.9, 1.9)
DIJET_DELTA_PHI_SELECTION_LABEL = "|#Delta#phi| > 2#pi/3"

# Compact bins for fast notebook checks; canonical bins remain user-tweakable.
TEST_SINGLE_JET_PT_BINS = [(50, 80), (180, 200), (400, 500)]
TEST_DIJET_PTAVE_BINS = [(50, 100), (100, 180), (180, 300), (300, 500), (500, 1000)]

DIJET_ETA_CUTS = (1.5, 1.6, 1.7, 1.8, 1.9, 2.4, 2.5, 3.0)

# Producer acceptance is inclusive (<=1.9) on BOTH selected jets; this index
# does not impose an additional cut on the average dijet eta. The legacy plot
# selection label above uses '>' although the implemented delta-phi cut is >=.
# processForestSimple.C etaCuts index for |eta_CM^jet| <= 1.9.
STANDARD_DIJET_ETA_CUT_INDEX = DIJET_ETA_CUTS.index(1.9)


def pt_bins(observable: str, mode: str = "normal") -> list[tuple[float, float]]:
    """Select existing single-jet or dijet intervals without changing edges."""
    if mode not in {"normal", "standard", "test"}:
        raise ValueError(f"Unknown pT binning {mode!r}; use normal (standard) or test")
    tables = {
        "single": (SINGLE_JET_PT_BINS, TEST_SINGLE_JET_PT_BINS),
        "dijet": (DIJET_PTAVE_BINS, TEST_DIJET_PTAVE_BINS),
    }
    if observable not in tables:
        raise ValueError(f"Unknown observable {observable!r}; use single or dijet")
    return list(tables[observable][mode == "test"])


def eta_cut_index(value: float = 1.9) -> int:
    """Translate physical per-jet acceptance to the current producer index."""
    try:
        return DIJET_ETA_CUTS.index(value)
    except ValueError as error:
        raise ValueError(f"Unsupported eta cut {value!r}; choose {DIJET_ETA_CUTS}") from error

# Existing regularization notebook test-trigger assignment, confirmed for the
# other downstream data notebooks. This is an input assignment, not stitching.
TEST_DIJET_TRIGGER_BY_INTERVAL = {
    (50, 100): 'MinimumBias', (100, 180): 'Jet80',
    (180, 300): 'Jet100', (300, 500): 'Jet100', (500, 1000): 'Jet100',
}


def trigger_dijet_bins(normal_bins: dict, mode: str = 'normal') -> dict:
    """Preserve notebook-specific normal bins or use the confirmed test map."""
    if mode in {'normal', 'standard'}:
        return {trigger: list(intervals) for trigger, intervals in normal_bins.items()}
    if mode != 'test':
        raise ValueError(f'Unknown dijet binning {mode!r}; use normal or test')
    aliases = {'mb': 'MinimumBias', 'jet60': 'Jet60', 'jet80': 'Jet80', 'jet100': 'Jet100'}
    return {trigger: [interval for interval in TEST_DIJET_PTAVE_BINS
                      if TEST_DIJET_TRIGGER_BY_INTERVAL[interval] == aliases.get(trigger, trigger)]
            for trigger in normal_bins}


def trigger_single_jet_bins(normal_bins: dict, mode: str = 'normal') -> dict:
    """Use existing single-jet test intervals within confirmed trigger ranges."""
    if mode in {'normal', 'standard'}:
        return {trigger: list(intervals) for trigger, intervals in normal_bins.items()}
    if mode != 'test':
        raise ValueError(f'Unknown single-jet binning {mode!r}; use normal or test')
    assigned = {}
    for interval in TEST_SINGLE_JET_PT_BINS:
        triggers = {trigger for (low, high), trigger in TEST_DIJET_TRIGGER_BY_INTERVAL.items()
                    if low <= interval[0] and interval[1] <= high}
        if len(triggers) != 1:
            raise ValueError(f'No unique confirmed test trigger for single-jet bin {interval}')
        assigned[interval] = triggers.pop()
    aliases = {'mb': 'MinimumBias', 'jet60': 'Jet60', 'jet80': 'Jet80', 'jet100': 'Jet100'}
    return {trigger: [interval for interval in TEST_SINGLE_JET_PT_BINS
                      if assigned[interval] == aliases.get(trigger, trigger)]
            for trigger in normal_bins}

# Exact stored EPPS21 intervals. Wider intervals overlap the standard analysis
# bins and must be projected from raw yields, never averaged from normalized bins.
NPDF_TRIGGER_BINS = {
    'mb': ((50, 60), (60, 70), (70, 80)),
    'jet60': ((80, 90), (90, 120)),
    'jet100': ((120, 160), (200, 300)),
}
NPDF_PTAVE_BINS = tuple(interval for bins in NPDF_TRIGGER_BINS.values() for interval in bins)
SYSTEMATIC_PTAVE_BINS = tuple(dict.fromkeys((*DIJET_PTAVE_BINS, *NPDF_PTAVE_BINS)))


def systematic_pt_bins(mode='normal'):
    """Canonical plus exact nPDF bins; retain existing compact test mode."""
    return list(SYSTEMATIC_PTAVE_BINS) if mode in {'normal', 'standard'} else pt_bins('dijet', mode)


def with_npdf_trigger_bins(normal_bins):
    """Extend existing assignments without changing any other interval."""
    aliases = {'MinimumBias': 'mb', 'Jet60': 'jet60', 'Jet80': 'jet80', 'Jet100': 'jet100'}
    assigned = {interval: trigger for trigger, intervals in NPDF_TRIGGER_BINS.items() for interval in intervals}
    for trigger, intervals in normal_bins.items():
        for interval in intervals:
            owner = assigned.get(tuple(interval))
            if owner is not None and aliases.get(trigger, trigger) != owner:
                raise ValueError(f'nPDF interval {interval} requires {owner}, not {trigger}')
    return {trigger: list(dict.fromkeys((*bins, *NPDF_TRIGGER_BINS.get(aliases.get(trigger, trigger), ()))))
            for trigger, bins in normal_bins.items()}
