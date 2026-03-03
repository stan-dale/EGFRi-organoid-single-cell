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
ANALYSIS_DIR = PROJECT_ROOT / "analysis"
FIGURES_DIR = PROJECT_ROOT / "figures"
UTILITIES_DIR = PROJECT_ROOT / "utilities"
QC_DIR = ANALYSIS_DIR / "qc"
NOTEBOOKS_DIR = PROJECT_ROOT / "notebooks"

# ── Analysis subdirectories ──
MANUAL_LABELLED_DIR = ANALYSIS_DIR / "manual_labelled"
MANUAL_LABELLED_2_DIR = ANALYSIS_DIR / "manual_labelled_2"
CELLASSIGN_DIR = ANALYSIS_DIR / "cellassign_objects"

# ── CellRank paths ──
CELLRANK_DIR = ANALYSIS_DIR / "cellrank"
CELLRANK_FIGURES_DIR = FIGURES_DIR / "cellrank"

# ── DGE paths ──
DGE_INPUT_DIR = DATA_DIR / "dge" / "pseudobulk-csvs"
DGE_OUTPUT_DIR = DATA_DIR / "dge" / "pydeseq-output"

# ── GSEA paths ──
GSEA_OUTPUT_DIR = DATA_DIR / "gsea" / "gseapy-output"

# ── Dataset registry ──
# Maps dataset key -> dict with h5ad paths for different processing stages.
# Not all stages exist for every dataset.
DATASETS = {
    "D2_DZ": {
        "raw": DATA_DIR / "egfDuod_D2_DZ.h5ad",
        "clustered": ANALYSIS_DIR / "clustered_egfDuod_D2_DZ.h5ad",
        "labelled": MANUAL_LABELLED_2_DIR / "d2_dz_manual_labels_2.h5ad",
        "condition": "Dz",
    },
    "D2_Lapa": {
        "raw": DATA_DIR / "egfDuod_D2_Lapa_DZ.h5ad",
        "clustered": ANALYSIS_DIR / "clustered_egfDuod_D2_Lapa_DZ.h5ad",
        "labelled": MANUAL_LABELLED_2_DIR / "d2_lapa_manual_labels_2.h5ad",
        "condition": "Lapa",
    },
    "D4_DZ": {
        "raw": DATA_DIR / "egfDuod_D4_DZ.h5ad",
        "clustered": ANALYSIS_DIR / "clustered_egfDuod_D4_DZ.h5ad",
        "labelled": MANUAL_LABELLED_DIR / "d4_dz_manual_labels.h5ad",
        "condition": "Dz",
    },
    "D4_Lapa": {
        "raw": DATA_DIR / "egfDuod_D4_Lapa_DZ.h5ad",
        "clustered": ANALYSIS_DIR / "clustered_egfDuod_D4_Lapa_DZ.h5ad",
        "labelled": MANUAL_LABELLED_DIR / "d4_lapa_manual_labels.h5ad",
        "condition": "Lapa",
    },
    "D4_AS": {
        "raw": DATA_DIR / "egfDuod_D4_AS_DZ.h5ad",
        "clustered": ANALYSIS_DIR / "clustered_egfDuod_D4_AS_DZ.h5ad",
        "labelled": MANUAL_LABELLED_DIR / "d4_as_manual_labels.h5ad",
        "condition": "AS",
    },
    "D10_Lapa": {
        "raw": DATA_DIR / "egfDuod_D10_Lapa_DZ.h5ad",
        "clustered": ANALYSIS_DIR / "clustered_egfDuod_D10_Lapa_DZ.h5ad",
        "labelled": CELLASSIGN_DIR / "d10_lapa_predictions.h5ad",
        "condition": "Lapa",
    },
    "G6": {
        "raw": DATA_DIR / "egfDuod_G6_DZ.h5ad",
        "clustered": ANALYSIS_DIR / "clustered_egfDuod_G6_DZ.h5ad",
        "labelled": MANUAL_LABELLED_DIR / "g6_manual_labels.h5ad",
        "condition": "Dz",
    },
}
