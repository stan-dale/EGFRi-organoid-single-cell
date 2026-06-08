"""
Step 2 — Define the secretory-progenitor compartment in D4_Lapa, then subcluster it.

Rationale (per user):
  - D4 is EARLY: little ATOH1, so a strict secretory-TF gate would be too sparse.
  - Use NEGATIVE SELECTION: secretory compartment = cells that are NOT committed
    stem and NOT committed absorptive. Proliferation is an orthogonal cell-cycle
    state, so proliferating cells are KEPT (flagged), not excluded.
  - Cross-check the negatively-selected pool against a liberal POSITIVE secretory
    score so the definition is transparent / not arbitrary.
  - Re-embed the subset (fresh HVG -> PCA -> UMAP -> Leiden) for real subcluster
    resolution rather than reusing the global embedding.

Inflammation scoring is deferred to step 3.

Inputs : outputs/d4_lapa_reclustered.h5ad   (from step 1; has layers counts+lognorm, X_pca, X_umap)
Outputs: outputs/d4_lapa_secretory.h5ad      (re-embedded secretory subset)
         outputs/02_lineage_class_breakdown.csv
         figures/02_*.png
"""
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

sc.settings.verbosity = 1
HERE = Path(__file__).resolve().parent.parent
FIG = HERE / "figures"; OUT = HERE / "outputs"
FIG.mkdir(exist_ok=True); OUT.mkdir(exist_ok=True)

ad = sc.read_h5ad(OUT / "d4_lapa_reclustered.h5ad")
print("loaded", ad.shape)

# ----------------------------------------------------------------------------
# 1. Lineage programs for negative selection. Liberal secretory set (early D4).
# ----------------------------------------------------------------------------
programs = {
    "Stem":       ["LGR5", "OLFM4", "ASCL2", "AXIN2", "SMOC2", "RNF43", "SOX9", "EPHB2"],
    "Absorptive": ["FABP2", "FABP1", "APOA1", "APOA4", "ALPI", "KRT20", "RBP2", "APOB", "ANPEP", "SI"],
    "Prolif":     ["MKI67", "TOP2A", "CDK1", "PCNA", "UBE2C", "BIRC5", "CCNB1", "CENPF"],
    # liberal secretory: progenitor TFs + early secretory effectors
    "Secretory":  ["ATOH1", "NEUROG3", "NEUROD1", "INSM1", "SOX4", "SPDEF", "DLL1",
                   "DLL3", "HES6", "BCL2", "PROX1", "MUC2", "TFF3", "SPINK4",
                   "CHGA", "CHGB", "FCGBP"],
}
def present(genes): return [g for g in genes if g in ad.var_names]
for name, genes in programs.items():
    g = present(genes)
    sc.tl.score_genes(ad, g, score_name=f"prog_{name}", use_raw=False, layer="lognorm")
    print(f"prog_{name}: {len(g)}/{len(genes)} present -> {g}")

# ----------------------------------------------------------------------------
# 2. z-score programs and classify by negative selection.
#    committed_X = X is the argmax program AND its z-score > ZMIN.
#    Secretory compartment = everything that is NOT committed Stem/Absorptive.
# ----------------------------------------------------------------------------
ZMIN = 0.25
progs = ["Stem", "Absorptive", "Prolif", "Secretory"]
Z = pd.DataFrame({p: ad.obs[f"prog_{p}"] for p in progs}, index=ad.obs_names)
Z = (Z - Z.mean()) / Z.std()
for p in progs:
    ad.obs[f"z_{p}"] = Z[p].values

argmax = Z.idxmax(axis=1)
argmax_pos = Z.max(axis=1) > ZMIN

cls = pd.Series("Secretory", index=ad.obs_names, dtype=object)
cls[(argmax == "Absorptive") & argmax_pos] = "Absorptive"
cls[(argmax == "Stem") & argmax_pos] = "Stem"
# Proliferation is orthogonal: don't pull prolif cells out of the secretory pool,
# but record the call so we can flag them later.
ad.obs["lineage_class"] = pd.Categorical(cls, categories=["Stem", "Absorptive", "Secretory"])
ad.obs["is_proliferating"] = ((argmax == "Prolif") & argmax_pos).map({True: "prolif", False: "non"}).astype("category")

print("\n=== lineage_class breakdown ===")
print(ad.obs["lineage_class"].value_counts())
print("\n=== lineage_class x manual_label ===")
ct = pd.crosstab(ad.obs["lineage_class"], ad.obs["manual_label"])
print(ct.to_string())
ct.to_csv(OUT / "02_lineage_class_breakdown.csv")

# Cross-check: liberal positive secretory score, top tertile, vs the negative-selection set
sec_pos = ad.obs["z_Secretory"] > 0
overlap = pd.crosstab(ad.obs["lineage_class"] == "Secretory", sec_pos,
                      rownames=["neg_select_secretory"], colnames=["pos_score>0"])
print("\n=== negative-selection vs positive-score overlap ===")
print(overlap.to_string())

# ----------------------------------------------------------------------------
# 3. Diagnostic figure on the GLOBAL umap: how the compartment was carved.
# ----------------------------------------------------------------------------
fig, axes = plt.subplots(2, 4, figsize=(22, 10))
for ax, key in zip(axes.flat, ["lineage_class", "is_proliferating", "manual_label",
                               "prog_Stem", "prog_Absorptive", "prog_Secretory",
                               "prog_Prolif"]):
    cmap = None if ad.obs[key].dtype.name == "category" else "viridis"
    sc.pl.umap(ad, color=key, ax=ax, show=False, size=6, title=key,
               legend_fontsize=7, cmap=cmap)
axes.flat[-1].axis("off")
fig.tight_layout(); fig.savefig(FIG / "02_compartment_definition.png", dpi=130); plt.close(fig)

# ----------------------------------------------------------------------------
# 4. Subset to secretory compartment and RE-EMBED from counts.
# ----------------------------------------------------------------------------
sec = ad[ad.obs["lineage_class"] == "Secretory"].copy()
print(f"\nsecretory compartment: {sec.n_obs} cells")

# rebuild a fresh working matrix from raw counts
sec.X = sec.layers["counts"].copy()
sc.pp.normalize_total(sec, target_sum=1e4)
sc.pp.log1p(sec)
sec.layers["lognorm"] = sec.X.copy()
sc.pp.highly_variable_genes(sec, n_top_genes=2000, flavor="seurat")
sec.raw = sec
sec_hvg = sec[:, sec.var["highly_variable"]].copy()
sc.pp.scale(sec_hvg, max_value=10)
sc.tl.pca(sec_hvg, n_comps=30)
sec.obsm["X_pca_sub"] = sec_hvg.obsm["X_pca"]
sec.uns["pca_sub"] = sec_hvg.uns["pca"]
sec.varm = {}

sc.pp.neighbors(sec, n_neighbors=20, n_pcs=30, use_rep="X_pca_sub", key_added="sub")
sc.tl.umap(sec, neighbors_key="sub")
sec.obsm["X_umap_sub"] = sec.obsm["X_umap"].copy()
sec.obsm["X_umap"] = ad[sec.obs_names].obsm["X_umap"]  # keep global umap too

for r in [0.3, 0.5, 0.8, 1.0, 1.5]:
    key = f"sub_leiden_r{r}"
    sc.tl.leiden(sec, resolution=r, key_added=key, flavor="igraph",
                 n_iterations=2, directed=False, neighbors_key="sub")
    print(f"{key} -> {sec.obs[key].nunique()} clusters")

sec.write(OUT / "d4_lapa_secretory.h5ad")
print("wrote", OUT / "d4_lapa_secretory.h5ad", sec.shape)

# ----------------------------------------------------------------------------
# 5. Subcluster overview on the NEW embedding.
# ----------------------------------------------------------------------------
panels = [f"sub_leiden_r{r}" for r in [0.3, 0.5, 0.8, 1.0, 1.5]] + \
         ["manual_label", "is_proliferating", "participant"]
fig, axes = plt.subplots(2, 4, figsize=(22, 10))
for ax, key in zip(axes.flat, panels):
    sc.pl.embedding(sec, "X_umap_sub", color=key, ax=ax, show=False, size=10,
                    title=key, legend_loc="on data" if key.startswith("sub_leiden") else "right margin",
                    legend_fontsize=6)
for ax in axes.flat[len(panels):]: ax.axis("off")
fig.tight_layout(); fig.savefig(FIG / "02_subcluster_overview.png", dpi=130); plt.close(fig)
print("DONE")
