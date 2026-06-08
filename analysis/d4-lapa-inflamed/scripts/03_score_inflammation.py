"""
Step 3 — Score the secretory subclusters for inflammation, triangulating 3 methods,
and ask: DISCRETE inflamed sub-population, or a CONTINUOUS gradient?

Scoring methods (deliberately different assumptions):
  1. score_genes  : Seurat/Tirosh background-subtracted mean        (prefix sg_)
  2. mean z-score : per-gene standardized expression, equal weight   (prefix mz_)
  3. AUCell       : rank-based recovery AUC, magnitude-robust         (prefix auc_)
(decoupler is unimportable in this env due to a numba/GSVA bug, so AUCell is
 implemented directly here.)

Marker sets: see analysis memory / src/markers.py + the shared figures.

Inputs : outputs/d4_lapa_secretory.h5ad   (re-embedded secretory subset, step 2)
Outputs: outputs/03_cluster_scores_<method>.csv
         outputs/03_discreteness_stats.csv
         outputs/03_subcluster_markers.csv
         figures/03_*.png
"""
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import skew, kurtosis
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

sc.settings.verbosity = 1
HERE = Path(__file__).resolve().parent.parent
FIG = HERE / "figures"; OUT = HERE / "outputs"

ad = sc.read_h5ad(OUT / "d4_lapa_secretory.h5ad")
ad.X = ad.layers["lognorm"].copy()
print("loaded", ad.shape)
CLUST = "sub_leiden_r1.0"

# ----------------------------------------------------------------------------
# Marker sets (filtered to genes that are actually expressed here).
# ----------------------------------------------------------------------------
raw_sets = {
    "Inflamed_canonical": ["REG3A", "REG1B", "REG1A"],
    "Inflamed_Gui":       ["REG3A", "REG1B", "REG1A", "SERPINA1", "LYZ", "IGFBP2", "ITLN1", "TFF3"],
    "ISP_core_LND":       ["LCN2", "NOS2", "DUOX2", "DUOXA2"],
    "ISP_Crohn":          ["LCN2", "DUOX2", "DUOXA2", "NOS2", "HLA-DRA", "HLA-DPA1",
                           "CD74", "CIITA", "CASP1", "CXCL1", "DMBT1"],
    "Alarmins":           ["LCN2", "CCL20", "IL1B", "CXCL8", "TNF"],
    "ER_stress":          ["XBP1", "HSPA5", "DDIT3", "ATF4", "ATF6", "ERN1",
                           "EIF2AK3", "CALR", "PDIA4"],
    "IFN_alpha":          ["STAT1", "STAT2", "IRF9", "IRF1", "IRF7", "ISG15", "IFIT1",
                           "IFIT3", "IFITM1", "IFITM3", "MX1", "OAS1", "RSAD2", "IFI6", "IFI27"],
    "IFN_gamma":          ["STAT1", "IRF1", "GBP1", "GBP2", "CXCL9", "CXCL10", "CXCL11",
                           "HLA-A", "HLA-B", "B2M", "TAP1", "PSMB8", "PSMB9", "CIITA"],
}
detected = np.asarray((ad.layers["counts"] > 0).sum(0)).ravel() > 0
det_genes = set(ad.var_names[detected])
sets = {}
for name, genes in raw_sets.items():
    present = [g for g in genes if g in det_genes]
    missing = [g for g in genes if g not in det_genes]
    sets[name] = present
    print(f"{name}: {len(present)}/{len(genes)} expressed; missing: {missing}")

# ----------------------------------------------------------------------------
# Method 1: score_genes
# ----------------------------------------------------------------------------
for name, genes in sets.items():
    sc.tl.score_genes(ad, genes, score_name=f"sg_{name}", use_raw=False)

# ----------------------------------------------------------------------------
# Method 2: mean z-score (per-gene standardized lognorm)
# ----------------------------------------------------------------------------
def gene_matrix(adata, genes):
    idx = [adata.var_names.get_loc(g) for g in genes]
    M = adata.layers["lognorm"][:, idx]
    return M.toarray() if sp.issparse(M) else np.asarray(M)

for name, genes in sets.items():
    M = gene_matrix(ad, genes)
    Z = (M - M.mean(0)) / (M.std(0) + 1e-9)
    ad.obs[f"mz_{name}"] = Z.mean(1)

# ----------------------------------------------------------------------------
# Method 3: AUCell (rank-based recovery AUC), implemented directly.
# ----------------------------------------------------------------------------
def aucell(adata, gene_sets, layer="lognorm", max_rank_frac=0.05, chunk=2500):
    X = adata.layers[layer]
    det = np.asarray((adata.layers["counts"] > 0).sum(0)).ravel() > 0
    gidx = np.where(det)[0]
    pos = {g: i for i, g in enumerate(adata.var_names[gidx])}
    n_genes = len(gidx)
    max_rank = int(np.ceil(max_rank_frac * n_genes))
    set_cols = {n: np.array([pos[g] for g in gs if g in pos]) for n, gs in gene_sets.items()}
    out = {n: np.zeros(adata.n_obs, dtype=np.float32) for n in gene_sets}
    Xd = X[:, gidx]
    for start in range(0, adata.n_obs, chunk):
        end = min(start + chunk, adata.n_obs)
        blk = Xd[start:end]
        blk = blk.toarray() if sp.issparse(blk) else np.asarray(blk)
        order = np.argsort(-blk, axis=1, kind="stable")          # 0 = highest expr
        ranks = np.empty_like(order)
        rows = np.arange(blk.shape[0])[:, None]
        ranks[rows, order] = np.arange(blk.shape[1])[None, :]    # 0-based rank
        for n, cols in set_cols.items():
            if len(cols) == 0:
                continue
            capped = np.minimum(ranks[:, cols], max_rank)
            out[n][start:end] = (max_rank - capped).sum(1) / (max_rank * len(cols))
    return out, max_rank

auc, max_rank = aucell(ad, sets)
print(f"AUCell max_rank = {max_rank}")
for name, v in auc.items():
    ad.obs[f"auc_{name}"] = v

# ----------------------------------------------------------------------------
# Per-cluster means for each method -> heatmaps + CSVs
# ----------------------------------------------------------------------------
methods = {"sg": "score_genes", "mz": "mean_z", "auc": "AUCell"}
order_clusters = sorted(ad.obs[CLUST].unique(), key=int)
for pre, label in methods.items():
    cols = [f"{pre}_{n}" for n in sets]
    tab = ad.obs.groupby(CLUST, observed=True)[cols].mean().loc[order_clusters]
    tab.columns = list(sets.keys())
    tab.to_csv(OUT / f"03_cluster_scores_{label}.csv")

# heatmap: z-score each signature across clusters so methods are visually comparable
fig, axes = plt.subplots(1, 3, figsize=(22, 7))
for ax, (pre, label) in zip(axes, methods.items()):
    tab = pd.read_csv(OUT / f"03_cluster_scores_{label}.csv", index_col=0)
    Z = (tab - tab.mean()) / (tab.std() + 1e-9)
    im = ax.imshow(Z.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(sets))); ax.set_xticklabels(sets.keys(), rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(len(tab))); ax.set_yticklabels(tab.index, fontsize=7)
    ax.set_title(f"{label}  (per-signature z across clusters)"); ax.set_ylabel(CLUST)
    fig.colorbar(im, ax=ax, fraction=0.046)
fig.tight_layout(); fig.savefig(FIG / "03_cluster_score_heatmaps.png", dpi=140); plt.close(fig)

# ----------------------------------------------------------------------------
# Triangulation: do the 3 methods agree on which cells are inflamed?
# ----------------------------------------------------------------------------
key_sig = "ISP_Crohn"
trio = ad.obs[[f"sg_{key_sig}", f"mz_{key_sig}", f"auc_{key_sig}"]].copy()
trio.columns = ["score_genes", "mean_z", "AUCell"]
corr = trio.corr(method="spearman")
print(f"\nSpearman corr of methods on {key_sig}:\n{corr.round(3)}")
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
pairs = [("score_genes", "mean_z"), ("score_genes", "AUCell"), ("mean_z", "AUCell")]
for ax, (a, b) in zip(axes, pairs):
    ax.scatter(trio[a], trio[b], s=3, alpha=0.15, c=ad.obs[f"auc_{key_sig}"], cmap="magma")
    ax.set_xlabel(a); ax.set_ylabel(b)
    ax.set_title(f"{a} vs {b}  (rho={corr.loc[a,b]:.2f})")
fig.suptitle(f"Method triangulation on {key_sig}", y=1.02)
fig.tight_layout(); fig.savefig(FIG / "03_method_triangulation.png", dpi=140, bbox_inches="tight"); plt.close(fig)

# ----------------------------------------------------------------------------
# UMAP (sub) coloured by key inflammation scores, all 3 methods
# ----------------------------------------------------------------------------
show_sigs = ["Inflamed_Gui", "ISP_Crohn", "ISP_core_LND"]
fig, axes = plt.subplots(3, 4, figsize=(22, 16))
for r_i, pre in enumerate(["sg", "mz", "auc"]):
    sc.pl.embedding(ad, "X_umap_sub", color=CLUST, ax=axes[r_i, 0], show=False, size=10,
                    legend_loc="on data", legend_fontsize=6, title=CLUST if r_i == 0 else "")
    for c_i, s in enumerate(show_sigs):
        sc.pl.embedding(ad, "X_umap_sub", color=f"{pre}_{s}", ax=axes[r_i, c_i + 1],
                        show=False, size=10, cmap="magma", title=f"{methods[pre]}: {s}")
fig.tight_layout(); fig.savefig(FIG / "03_umap_inflammation_by_method.png", dpi=130); plt.close(fig)

# ----------------------------------------------------------------------------
# DISCRETE vs GRADIENT diagnostics
# ----------------------------------------------------------------------------
def bimodality_coef(x):
    x = np.asarray(x); n = len(x)
    g = skew(x); k = kurtosis(x, fisher=True)
    denom = k + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    return (g ** 2 + 1) / denom   # >0.555 hints at bimodality/non-normality

rows = []
for s in sets:
    for pre, label in methods.items():
        v = ad.obs[f"{pre}_{s}"].values
        rows.append({"signature": s, "method": label, "BC": bimodality_coef(v),
                     "skew": skew(v), "exc_kurtosis": kurtosis(v, fisher=True)})
disc = pd.DataFrame(rows)
disc.to_csv(OUT / "03_discreteness_stats.csv", index=False)
print("\n=== discreteness (BC>0.555 hints bimodal) ===")
print(disc.pivot(index="signature", columns="method", values="BC").round(3).to_string())

# (a) global histograms of key signature per method  (b) per-cluster violins
fig, axes = plt.subplots(2, 3, figsize=(20, 11))
for c_i, (pre, label) in enumerate(methods.items()):
    v = ad.obs[f"{pre}_{key_sig}"].values
    axes[0, c_i].hist(v, bins=80, color="#444")
    bc = bimodality_coef(v)
    axes[0, c_i].set_title(f"{label}: {key_sig}  (BC={bc:.3f})")
    axes[0, c_i].set_xlabel("score"); axes[0, c_i].set_ylabel("cells")
    # per-cluster violin
    data = [ad.obs.loc[ad.obs[CLUST] == c, f"{pre}_{key_sig}"].values for c in order_clusters]
    axes[1, c_i].violinplot(data, showmeans=True, widths=0.9)
    axes[1, c_i].set_xticks(range(1, len(order_clusters) + 1))
    axes[1, c_i].set_xticklabels(order_clusters, fontsize=7)
    axes[1, c_i].set_title(f"{label}: {key_sig} by {CLUST}")
    axes[1, c_i].set_xlabel("cluster")
fig.tight_layout(); fig.savefig(FIG / "03_discreteness.png", dpi=140); plt.close(fig)

# (c) joint: secretory program vs inflammation — corner (discrete) or continuum?
fig, axes = plt.subplots(1, 3, figsize=(20, 6))
for ax, s in zip(axes, show_sigs):
    sccat = ad.obs["prog_Secretory"].values
    inf = ad.obs[f"auc_{s}"].values
    hb = ax.hexbin(sccat, inf, gridsize=45, cmap="viridis", mincnt=1, bins="log")
    ax.set_xlabel("secretory program (prog_Secretory)"); ax.set_ylabel(f"AUCell {s}")
    ax.set_title(f"secretory vs {s}")
    fig.colorbar(hb, ax=ax, label="log10 cells")
fig.tight_layout(); fig.savefig(FIG / "03_secretory_vs_inflammation.png", dpi=140); plt.close(fig)

# ----------------------------------------------------------------------------
# Characterise subclusters: rank_genes_groups
# ----------------------------------------------------------------------------
sc.tl.rank_genes_groups(ad, CLUST, method="wilcoxon", use_raw=False)
rg = sc.get.rank_genes_groups_df(ad, group=None)
top = (rg.sort_values(["group", "scores"], ascending=[True, False])
         .groupby("group").head(15))
top.to_csv(OUT / "03_subcluster_markers.csv", index=False)
print("\n=== top markers per subcluster ===")
for c in order_clusters:
    genes = top.loc[top["group"] == c, "names"].head(12).tolist()
    print(f"  c{c}: {', '.join(genes)}")

ad.write(OUT / "d4_lapa_secretory.h5ad")
print("\nupdated", OUT / "d4_lapa_secretory.h5ad")
print("DONE")
