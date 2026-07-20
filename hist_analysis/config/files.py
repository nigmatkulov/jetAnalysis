"""Machine-local pPb8160 data, embedding, and Pythia output locations."""

from __future__ import annotations

from pathlib import Path

BASE_DIR = Path("/Users/gnigmat/cernbox/ana/pPb8160")

DATA_DIRS = {
    "pgoing": BASE_DIR / "exp" / "pgoing",
    "Pbgoing": BASE_DIR / "exp" / "Pbgoing",
    "combined": BASE_DIR / "exp",
}

MC_DIRS = {
    "embedding": {
        "pgoing": BASE_DIR / "embedding" / "pgoing",
        "Pbgoing": BASE_DIR / "embedding" / "Pbgoing",
        "combined": BASE_DIR / "embedding",
    },
    "pythia": {
        "pgoing": BASE_DIR / "pythia" / "pgoing",
        "Pbgoing": BASE_DIR / "pythia" / "Pbgoing",
        "combined": BASE_DIR / "pythia",
    },
}
