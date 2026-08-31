"""
Unified colour palette for cell type annotations.

Handles inconsistent naming (e.g. "Inflammed" vs "Inflamed",
"PCs" vs "Progenitor cells") via normalize_celltype_name().
"""


# ── Canonical palette (keys are the *normalised* names) ──

_CANONICAL_PALETTE = {
    "ISCs":                       "#7f7f7f",
    "PCs":                        "#1ce6ff",
    "Secretory PCs":              "#00897b",
    "Proliferating PCs":          "#ff4a46",
    "Inflamed Secretory PCs":     "#355fdb",
    "Enterocytes":                "#ffff00",
    "Enterocyte PCs":             "#ffff00",
    "Goblet cells":               "#b15928",
    "NEUROG3+ progenitor cells":  "#8a7fe0",
    "EECs":                       "#6651d1",
}


# ── Name normalisation ──

_ALIASES = {
    # Spelling variants
    "Inflammed secretory PCs":          "Inflamed Secretory PCs",
    "Inflammed Secretory PCs":          "Inflamed Secretory PCs",
    "Inflamed secretory PCs":           "Inflamed Secretory PCs",
    "Inflammed progenitor cells":       "Inflamed Secretory PCs",
    "Inflamed progenitor cells":        "Inflamed Secretory PCs",
    # Long-form names
    "Secretory progenitor cells":       "Secretory PCs",
    "Proliferating Progenitor cells":   "Proliferating PCs",
    "Proliferating progenitor cells":   "Proliferating PCs",
    "Progenitor cells":                 "PCs",
}


def normalize_celltype_name(name):
    """Map common aliases to canonical cell type names."""
    return _ALIASES.get(name, name)


def get_color(name):
    """Return hex colour for a cell type (handles aliases)."""
    canonical = normalize_celltype_name(name)
    return _CANONICAL_PALETTE.get(canonical, "#808080")


# ── Backwards-compatible flat dict (includes all alias keys) ──

celltype_palette = {}
celltype_palette.update(_CANONICAL_PALETTE)
for alias, canonical in _ALIASES.items():
    celltype_palette[alias] = _CANONICAL_PALETTE[canonical]


# ── EEC subtype palette (D10 Lapa sub-analysis) ──
#
# Frozen from the Set1_9 cycle the sub-analysis originally derived colours from,
# so the granular panels are unchanged. Enterochromaffin is pinned to orange in
# every figure, including the coarse two-way split where position alone had made
# it blue. Early EECs are grey to read as an immature, unresolved state.

eec_palette = {
    "X cells":                "#e41a1c",
    "D cells":                "#377eb8",
    "I cells":                "#4daf4a",
    "K cells":                "#984ea3",
    "Enterochromaffin cells": "#ff7f00",
    "EC cells":               "#ff7f00",
    "Early EECs":             "#999999",
    # Coarse two-way split
    "Peptide EECs":   "#999999",
    # Intermediate labels, reassigned before the final call
    "ARX+ cells":             "#f781bf",
    "NEUROG3+ PCs":           "#a65628",
}
