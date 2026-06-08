"""
Quick figure: visualise the 'Secretory PC' compartment label that step 2 defined.
  - Global UMAP: all D4_Lapa cells coloured Stem / Absorptive / Secretory PC.
  - Global UMAP: Secretory PC cells highlighted, everything else grey.
  - Subset UMAP: the secretory compartment on its own re-embedding, by subcluster.
"""
import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
FIG = HERE / "figures"; OUT = HERE / "outputs"

full = sc.read_h5ad(OUT / "d4_lapa_reclustered.h5ad", backed="r")
sec = sc.read_h5ad(OUT / "d4_lapa_secretory.h5ad")

# Reconstruct the 3-way lineage label on ALL cells (so we can show what was excluded).
sec_ids = set(sec.obs_names)
lab = np.where([n in sec_ids for n in full.obs_names], "Secretory PC", "other")
import pandas as pd
g = full.obs.copy()
g["compartment"] = lab
# split 'other' into Stem / Absorptive using the same programs (cheap recompute on full)
adf = sc.read_h5ad(OUT / "d4_lapa_reclustered.h5ad")
programs = {
    "Stem":       ["LGR5", "OLFM4", "ASCL2", "AXIN2", "SMOC2", "RNF43", "SOX9", "EPHB2"],
    "Absorptive": ["FABP2", "FABP1", "APOA1", "APOA4", "ALPI", "KRT20", "RBP2", "APOB", "ANPEP", "SI"],
}
for n, gs in programs.items():
    gg = [x for x in gs if x in adf.var_names]
    sc.tl.score_genes(adf, gg, score_name=f"p_{n}", use_raw=False, layer="lognorm")
comp = pd.Series("Secretory PC", index=adf.obs_names, dtype=object)
not_sec = ~adf.obs_names.isin(sec.obs_names)
which = adf.obs.loc[not_sec, ["p_Stem", "p_Absorptive"]].idxmax(axis=1).str.replace("p_", "", regex=False)
comp.loc[not_sec] = which.values
adf.obs["compartment"] = pd.Categorical(comp, categories=["Stem", "Absorptive", "Secretory PC"])

palette = {"Stem": "#7fb069", "Absorptive": "#e09f3e", "Secretory PC": "#9d4edd"}

fig, axes = plt.subplots(1, 3, figsize=(21, 6.5))
sc.pl.umap(adf, color="compartment", ax=axes[0], show=False, size=8,
           palette=palette, title="D4_Lapa compartments (global UMAP)")
# highlight secretory only
adf.obs["is_secretory"] = (adf.obs["compartment"] == "Secretory PC").map(
    {True: "Secretory PC", False: "excluded"}).astype("category")
sc.pl.umap(adf, color="is_secretory", ax=axes[1], show=False, size=8,
           palette={"Secretory PC": "#9d4edd", "excluded": "#d9d9d9"},
           title=f"Secretory PC label  (n={int((adf.obs['compartment']=='Secretory PC').sum())})")
# subset on its own embedding
sc.pl.embedding(sec, "X_umap_sub", color="sub_leiden_r1.0", ax=axes[2], show=False,
                size=12, legend_loc="on data", legend_fontsize=7,
                title="Secretory PC subclusters (subset UMAP)")
fig.tight_layout(); fig.savefig(FIG / "04_secretory_pc_label.png", dpi=140); plt.close(fig)
print("wrote", FIG / "04_secretory_pc_label.png")
