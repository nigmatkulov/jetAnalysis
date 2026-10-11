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
