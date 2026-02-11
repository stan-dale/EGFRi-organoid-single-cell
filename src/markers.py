"""
Marker gene definitions for cell type identification.

Centralises all marker dictionaries used across the project's notebooks.
"""


# ── Main epithelial cell type markers (used by CellAssign & dotplots) ──

cell_type_markers = {
    "ISCs": ["LGR5", "AXIN2", "ASCL2", "OLFM4"],
    "PCs": ["DLL1", "NEUROG3", "CD44"],
    "Proliferating PCs": ["DLL1", "NEUROG3", "CD44", "MKI67"],
    "Secretory PCs": ["NEUROG3", "INSM1", "NEUROD1", "SOX4", "ATOH1"],
    "Enterocytes": ["KRT20", "FABP2", "FABP1", "ALPI"],
    "Goblet cells": ["MUC2", "FCGBP", "GFI1"],
    "EECs": ["CHGA", "CHGB", "NEUROG3", "NEUROD1", "PAX4", "PCSK1"],
}


# ── EEC subtype markers (used in D10_Lapa EEC sub-analysis) ──

EEC_markers = {
    "NEUROG3+ PCs": "NEUROG3",
    "Early EECs": ["HES1", "ASCL2", "ATOH1", "PAX4"],
    "Pre EECs": ["HES6", "ASCL1"],
    "T4 cells": ["PAX6", "PDX1"],
    "X cells": "GHRL",
    "D cells": "SST",
    "I cells": "CCK",
    "K cells": "GIP",
    "Enterochromaffin cells": ["CHGA", "CHGB", "SLC18A1", "LMX1A"],
}


# ── Additional marker sets found across notebooks ──

More_Isc_markers = ["BMI1", "EPHB2", "SOX9", "PROM1"]

Inflammed_Secretory_PC_markers = ["REG3A", "REG1B", "REG1A"]


# ── Cell cycle markers (Tirosh et al.) ──

cell_cycle_markers = {
    "S_genes": [
        "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
        "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP",
        "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
        "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
        "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
        "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
    ],
    "G2M_genes": [
        "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
        "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
        "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
        "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",
        "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
        "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
        "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
        "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
    ],
}


# ── Utility ──

# ── Pathway gene sets (used by GSEA enrichment notebooks) ──

pathway_gene_sets = {
    "EGFR_SIGNALING": [
        "EGFR", "ERBB2", "ERBB3", "ERBB4",
        "GRB2", "SOS1", "SOS2",
        "KRAS", "HRAS", "NRAS",
        "RAF1", "BRAF", "ARAF",
        "MAP2K1", "MAP2K2",
        "MAPK1", "MAPK3",
        "ELK1", "FOS", "JUN",
        "PIK3CA", "PIK3CB", "PIK3CD",
        "AKT1", "AKT2", "AKT3",
        "MTOR", "RPS6KB1",
        "STAT3", "STAT5A", "STAT5B",
    ],
    "INTERFERON_ALPHA_RESPONSE": [
        "STAT1", "STAT2", "IRF9",
        "IRF1", "IRF2", "IRF3", "IRF7",
        "ISG15", "ISG20", "IFIT1", "IFIT2", "IFIT3",
        "IFITM1", "IFITM2", "IFITM3",
        "MX1", "MX2",
        "OAS1", "OAS2", "OAS3",
        "RSAD2", "IFI6", "IFI27", "IFI35", "IFI44",
        "IFIH1", "DDX58",
    ],
    "INTERFERON_GAMMA_RESPONSE": [
        "STAT1", "IRF1", "IRF8",
        "GBP1", "GBP2", "GBP3", "GBP4", "GBP5",
        "CXCL9", "CXCL10", "CXCL11",
        "IDO1", "IDO2",
        "HLA-A", "HLA-B", "HLA-C",
        "HLA-DRA", "HLA-DRB1", "HLA-DQA1",
        "B2M", "TAP1", "TAP2",
        "PSMB8", "PSMB9", "PSMB10",
        "CIITA", "CD274", "PDCD1LG2",
    ],
}


def filter_present_genes(adata, gene_list):
    """
    Return only genes from *gene_list* that exist in *adata.var_names*.

    Also prints which genes are missing (useful for QC).

    Parameters
    ----------
    adata : AnnData
    gene_list : list[str]

    Returns
    -------
    present : list[str]
    missing : list[str]
    """
    varset = set(adata.var_names)
    present = [g for g in gene_list if g in varset]
    missing = [g for g in gene_list if g not in varset]
    if missing:
        print(f"  Missing genes (ignored): {', '.join(missing)}")
    return present, missing


def filter_marker_dict(adata, markers):
    """
    Filter a {cell_type: [genes]} dict to only genes present in *adata*.

    Returns
    -------
    filtered : dict
    all_missing : dict  (cell_type -> list of missing genes)
    """
    varset = set(adata.var_names)
    filtered = {}
    all_missing = {}
    for ct, genes in markers.items():
        if isinstance(genes, str):
            genes = [genes]
        present = [g for g in genes if g in varset]
        miss = [g for g in genes if g not in varset]
        if present:
            filtered[ct] = present
        if miss:
            all_missing[ct] = miss
    if all_missing:
        print("Missing genes (ignored):")
        for ct, lst in all_missing.items():
            print(f"  {ct}: {', '.join(lst)}")
    return filtered, all_missing
