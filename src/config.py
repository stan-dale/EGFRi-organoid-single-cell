"""
Project paths and dataset registry for the Breault Lab single-cell project.

Usage from any notebook (works regardless of CWD):
    import sys
    from pathlib import Path
    _p = Path(".").resolve()
    while not (_p / "src" / "config.py").exists() and _p != _p.parent:
        _p = _p.parent
    sys.path.insert(0, str(_p))
    from src.config import *
"""

from pathlib import Path

# ── Project root (one level above src/) ──
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# ── Standard directories ──
DATA_DIR = PROJECT_ROOT / "data"
FIGURES_DIR = PROJECT_ROOT / "figures"
UTILITIES_DIR = PROJECT_ROOT / "utilities"
NOTEBOOKS_DIR = PROJECT_ROOT / "notebooks"

# ── Data subdirectories ──
DATA_OBJECTS_DIR = DATA_DIR / "data-objects"
ANALYSIS_FILES_DIR = DATA_DIR / "analysis-files"

RAW_DIR = DATA_OBJECTS_DIR / "raw"
CLUSTERED_DIR = DATA_OBJECTS_DIR / "clustered"
LABELLED_DIR = DATA_OBJECTS_DIR / "labelled"
CELLASSIGN_DIR = DATA_OBJECTS_DIR / "cellassign"
QC_DIR = DATA_OBJECTS_DIR / "qc"

# ── DGE paths ──
DGE_INPUT_DIR = ANALYSIS_FILES_DIR / "dge" / "pseudobulk-csvs"
DGE_OUTPUT_DIR = ANALYSIS_FILES_DIR / "dge" / "pydeseq-output"

# ── GSEA paths ──
GSEA_OUTPUT_DIR = ANALYSIS_FILES_DIR / "gsea" / "gseapy-output"

# ── CellRank paths ──
CELLRANK_DIR = DATA_OBJECTS_DIR / "cellrank"
CELLRANK_FIGURES_DIR = FIGURES_DIR / "cellrank"

# ── Backwards compatibility ──
ANALYSIS_DIR = DATA_DIR

# ── Dataset registry ──
# Maps dataset key -> dict with h5ad paths for different processing stages.
# Not all stages exist for every dataset.
DATASETS = {
    "D2_DZ": {
        "raw": RAW_DIR / "egfDuod_D2_DZ.h5ad",
        "clustered": CLUSTERED_DIR / "clustered_egfDuod_D2_DZ.h5ad",
        "labelled": LABELLED_DIR / "d2_dz_manual_labels.h5ad",
        "condition": "Dz",
    },
    "D2_Lapa": {
        "raw": RAW_DIR / "egfDuod_D2_Lapa_DZ.h5ad",
        "clustered": CLUSTERED_DIR / "clustered_egfDuod_D2_Lapa_DZ.h5ad",
        "labelled": LABELLED_DIR / "d2_lapa_manual_labels.h5ad",
        "condition": "Lapa",
    },
    "D4_DZ": {
        "raw": RAW_DIR / "egfDuod_D4_DZ.h5ad",
        "clustered": CLUSTERED_DIR / "clustered_egfDuod_D4_DZ.h5ad",
        "labelled": LABELLED_DIR /"d4_dz_manual_labels.h5ad",
        "condition": "Dz",
    },
    "D4_Lapa": {
        "raw": RAW_DIR / "egfDuod_D4_Lapa_DZ.h5ad",
        "clustered": CLUSTERED_DIR / "clustered_egfDuod_D4_Lapa_DZ.h5ad",
        "labelled": LABELLED_DIR /"d4_lapa_manual_labels.h5ad",
        "condition": "Lapa",
    },
    "D4_AS": {
        "raw": RAW_DIR / "egfDuod_D4_AS_DZ.h5ad",
        "clustered": CLUSTERED_DIR / "clustered_egfDuod_D4_AS_DZ.h5ad",
        "labelled": LABELLED_DIR /"d4_as_manual_labels.h5ad",
        "condition": "AS",
    },
    "D10_Lapa": {
        "raw": RAW_DIR / "egfDuod_D10_Lapa_DZ.h5ad",
        "clustered": CLUSTERED_DIR / "clustered_egfDuod_D10_Lapa_DZ.h5ad",
        "labelled": CELLASSIGN_DIR / "d10_lapa_predictions.h5ad",
        "condition": "Lapa",
    },
    "G6": {
        "raw": RAW_DIR / "egfDuod_G6_DZ.h5ad",
        "clustered": CLUSTERED_DIR / "clustered_egfDuod_G6_DZ.h5ad",
        "labelled": LABELLED_DIR /"g6_manual_labels.h5ad",
        "condition": "Dz",
    },
}
