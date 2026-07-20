"""Configuration for hist_analysis."""

from .files import BASE_DIR, DATA_DIRS, MC_DIRS
from .histograms import (
    COMMON_ETA_CM_RANGE,
    TEST_DIJET_PTAVE_BINS,
    TEST_SINGLE_JET_PT_BINS,
)

__all__ = [
    "BASE_DIR",
    "DATA_DIRS",
    "MC_DIRS",
    "COMMON_ETA_CM_RANGE",
    "TEST_DIJET_PTAVE_BINS",
    "TEST_SINGLE_JET_PT_BINS",
]
