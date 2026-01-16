cell_type_markers = {
    "ISCs": ["LGR5", "AXIN2", "ASCL2"],
    "PCs": ["DLL1", "NEUROG3", "CD44"],
    "Proliferating PCs": ["DLL1", "NEUROG3", "CD44", "MKI67"],
    "Secretory PCs": ["NEUROG3", "INSM1", "NEUROD1", "SOX4", "ATOH1"],
    "Enterocytes": ["KRT20", "FABP2", "FABP1", "ALPI"],
    "Goblet cells": ["MUC2", "FCGBP", "GFI1"],
    "EECs": ["CHGA", "CHGB", "NEUROG3", "NEUROD1", "PAX4", "PCSK1"]
}

cell_cycle_markers = {
    'S_genes': ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS",    "RFC2",     "RPA2",    
 "NASP",     "RAD51AP1", "GMNN",     "WDR76",    "SLBP",     "CCNE2",    "UBR7",     "POLD3",    "MSH2",  "ATAD2",    "RAD51",   "RRM2",     "CDC45",    "CDC6",     "EXO1",     "TIPIN",    "DSCC1",    "BLM",     
 "CASP8AP2", "USP1",     "CLSPN",    "POLA1",    "CHAF1B",   "BRIP1",    "E2F8"],
    'G2M_genes': [
        "HMGB2",   "CDK1",    "NUSAP1",  "UBE2C",   "BIRC5",   "TPX2",    "TOP2A",   "NDC80",   "CKS2",    "NUF2" ,  
"CKS1B"  , "MKI67"  , "TMPO"   , "CENPF"  , "TACC3"   ,"FAM64A" , "SMC4"    ,"CCNB2"   ,"CKAP2L" , "CKAP2"  ,
 "AURKB" ,  "BUB1"  ,  "KIF11" ,  "ANP32E",  "TUBB4B" , "GTSE1" ,  "KIF20B" , "HJURP"  , "CDCA3" ,  "HN1"   , 
 "CDC20" ,  "TTK"   ,  "CDC25C",  "KIF2C" ,  "RANGAP1", "NCAPD2",  "DLGAP5" , "CDCA2"  , "CDCA8" ,  "ECT2"  , 
 "KIF23" ,  "HMMR"  ,  "AURKA" ,  "PSRC1" ,  "ANLN"   , "LBR"   ,  "CKAP5"  , "CENPE"  , "CTCF"  ,  "NEK2"  , 
 "G2E3"  ,  "GAS2L3",  "CBX5"  ,  "CENPA" 
    ]   
}

## TO DO, uninfy the labels for each cell type, eg: "PCs" & "Progenitor cells"
## Should this live in vis.py instead?

celltype_palette = {
    "ISCs": "#7f7f7f",
    "PCs": "#1ce6ff",
    "Secretory PCs": "#355fdb",
    "Secretory progenitor cells": "#355fdb",
    "Proliferating Progenitor cells": "#ff4a46",
    "Proliferating progenitor cells": "#ff4a46",
    "Proliferating PCs": "#ff4a46",
    "Inflammed secretory PCs": "#ffff00",
    "Inflammed Secretory PCs": "#ffff00",
    "Inflammed progenitor cells": "#ffff00",
    "Enterocytes": "#ffff00",
    "Enterocyte PCs": "#ffff00",
    "Goblet cells": "#b15928",
    "NEUROG3+ progenitor cells": "#8a7fe0",
    "EECs": "#6651d1"
}


