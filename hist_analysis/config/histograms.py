"""Shared closure selections and the standard dijet eta-cut index."""

from __future__ import annotations

SINGLE_JET_PT_BINS = [(50, 80), (180, 200), (400, 500)]
DIJET_PTAVE_BINS = [(50, 100), (150, 200), (300, 500)]
COMMON_ETA_CM_RANGE = (-1.9, 1.9)

# Backward-compatible names exported by hist_analysis.config.
TEST_SINGLE_JET_PT_BINS = SINGLE_JET_PT_BINS
TEST_DIJET_PTAVE_BINS = DIJET_PTAVE_BINS

# processForestSimple.C etaCuts index for |eta_CM^jet| < 1.9.
STANDARD_DIJET_ETA_CUT_INDEX = 5
