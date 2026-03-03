import matplotlib.pyplot as plt
import palettable
import scanpy as sc

from pathlib import Path
import numpy as np
import pandas as pd

from .markers import filter_marker_dict

def plot_cluster_doublets(
    adata,
    cluster,
    cluster_key="leiden",
    score_key="doublet_score",
    cmap="Reds",
    bg_alpha=0.05,
    fg_size=8,
    bg_size=2
):
    """Plot cells on UMAP for a given cluster, colored by doublet score."""
    
    umap = adata.obsm["X_umap"]
    cl = adata.obs[cluster_key].astype(str)
    mask = (cl == str(cluster))
    scores = adata.obs[score_key].values

    fig, ax = plt.subplots(figsize=(5.2, 4.6))

    # background: all other cells in grey
    ax.scatter(
        umap[~mask, 0], umap[~mask, 1],
        s=bg_size, c="lightgrey", alpha=bg_alpha, linewidths=0
    )

    # foreground: this cluster, colored by doublet score (fixed scale 0–1)
    scatt = ax.scatter(
        umap[mask, 0], umap[mask, 1],
        c=scores[mask], cmap=cmap, s=fg_size,
        vmin=0, vmax=1, linewidths=0
    )

    # colorbar (no label, just the scale)
    cbar = plt.colorbar(scatt, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=9)
    cbar.set_label("")   # drop the label text

    # formatting
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(""); ax.set_ylabel("")

    # 🔹 Updated title: bold, larger, with colon
    ax.set_title(
        f"Cluster {cluster}: doublet score",
        fontsize=14, fontweight="bold", pad=12
    )

    plt.tight_layout()
    plt.show()
    return fig, ax






## To-do, amend the absolute paths with an import of Pathlib

def save_doublet_umaps(
    adata,
    object_name: str,
    cluster_key="leiden",
    score_key="doublet_score",
    base_dir="../figures/qc/doublet_score_umaps",
    dpi=300
):
    """
    Generate and save doublet-score UMAPs for each cluster in an AnnData object.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix containing UMAP coordinates.
    object_name : str
        Label used to name the save directory (e.g. 'd2_dz_clustered').
    cluster_key : str
        obs column containing cluster labels.
    score_key : str
        obs column containing doublet scores.
    base_dir : str or Path
        Base directory for saving plots.
    dpi : int
        Resolution for saved images.
    """

    # make save directory
    save_dir = Path(base_dir) / object_name
    save_dir.mkdir(parents=True, exist_ok=True)
    print(f"Saving doublet UMAPs to: {save_dir}")

    # sort clusters naturally
    def _nsort(vals):
        def _key(v):
            s = str(v)
            try:    return (0, int(s))
            except: return (1, s)
        return sorted(map(str, vals), key=_key)

    clusters = _nsort(adata.obs[cluster_key].unique())

    # loop over clusters
    for cl in clusters:
        fig, ax = plot_cluster_doublets(
            adata,
            cluster=cl,
            cluster_key=cluster_key,
            score_key=score_key
        )
        fig.savefig(
            save_dir / f"doublet_score_cluster_{cl}.png",
            dpi=dpi,
            bbox_inches="tight"
        )
        plt.close(fig)  # keep memory usage down







# --- Veres-style panel builder (Scanpy/AnnData) ---

import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.patheffects as PathEffects
from collections import Counter, OrderedDict
from sklearn.neighbors import NearestNeighbors
from palettable.colorbrewer.qualitative import Set1_9, Paired_10, Set2_8

# -------------------- typography to match paper --------------------
def setup_matplotlib_params():
    from matplotlib import rcParams
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Helvetica Neue', 'Arial', 'DejaVu Sans']
    rcParams['axes.titlesize'] = 5
    rcParams['axes.labelsize']  = 5
    rcParams['xtick.labelsize'] = 5
    rcParams['ytick.labelsize'] = 5
    rcParams['axes.linewidth']  = 0.5
    rcParams['xtick.major.width'] = 0.5
    rcParams['ytick.major.width'] = 0.5
    rcParams['pdf.fonttype'] = 42  # editable text in Illustrator

setup_matplotlib_params()

# -------------------- core color registry (as in figure_vis.py) --------------------
def _col(x): return np.array(x)/255.

class CoreColors:
    def __init__(self):
        base = dict(zip(
            ['red','blue','green','purple','orange','yellow','brown','pink','grey'],
            Set1_9.colors
        ))
        for k,v in base.items(): setattr(self, k, _col(v))
        for name,c in zip(['blue','green','red','orange','purple'], Paired_10.colors[::2]):
            setattr(self, 'pale_'+name, _col(c))
        self.teal        = _col(Set2_8.colors[0])
        self.pale_brown  = _col([193,128,108])
        self.light_grey  = _col(Set2_8.colors[-1])
        self.dark_green  = _col([52,117,50])
        self.dark_grey   = np.array([0.5,0.5,0.5])

core_colors = CoreColors()

# -------------------- label params (you can fully customize) --------------------
def make_label_params(category_names, palette=None):
    """
    Build a label->dict(color, short_label) map deterministically from the
    Veres palettes. Edit this to hard-code exact colors if desired.

    Parameters
    ----------
    category_names : list[str]
    palette : dict, optional
        {label: hex_color} to override specific label colors (e.g. celltype_palette).
        Labels not in the dict fall back to the Veres color cycle.
    """
    from matplotlib.colors import to_rgb
    fallback = (
        Set1_9.colors + Paired_10.colors[::2] + Set2_8.colors
    )
    lp = OrderedDict()
    fb_i = 0
    for name in category_names:
        if palette is not None and name in palette:
            col = np.array(to_rgb(palette[name]))
        else:
            col = _col(fallback[fb_i % len(fallback)])
            fb_i += 1
        lp[name] = dict(color=col, short_label=name)
    return lp

# -------------------- identical logic: sort to put homogeneous neighborhoods on top --------------------
def prepare_for_scatter(X, labels, label_params):
    labels = np.asarray(labels)
    nbrs = NearestNeighbors(n_neighbors=5, algorithm='ball_tree').fit(X) # try moving from 10 to 5
    knn = nbrs.kneighbors(X, return_distance=False)
    # fraction of neighbors with the same label; sort so dense/consistent are on top
    same = (labels[knn].T == labels[knn][:,0].T).mean(0)
    order = np.argsort(same)
    # map labels -> RGB
    colors = np.ones((len(labels),3))*0.5
    for k,v in label_params.items():
        colors[labels==k] = v['color']
    return X[order], colors[order]

# -------------------- the panel function (ports the notebook code) --------------------
def plot_veres_panel(
    adata,
    label_key,                 # e.g. "initial_cellassign_prediction_colors"
    stage_text=None,           # e.g. "D10_Lapa"
    ratio_order=None,          # list specifying order in pop-bar (optional)
    label_order=None,          # legend order (optional)
    palette=None,              # dict {label: hex_color} e.g. celltype_palette
    save=None,                 # path to .pdf/.png
    dpi=600
):
    # coords & labels
    X = adata.obsm['X_umap']
    labels = adata.obs[label_key].astype(str).values
    cats = list(pd.unique(labels)) if label_order is None else label_order

    # label params (color + short label)
    lp = make_label_params(cats, palette=palette)

    # ===== panel sizing exactly like their notebook =====
    mm_per_inch = 25.4
    panel_size_in = ((89/2)/mm_per_inch, 25.4/mm_per_inch)  # (W, H)
    cell_pop_bar_h = 0.07
    heights = ((1 - cell_pop_bar_h)*panel_size_in[1], cell_pop_bar_h*panel_size_in[1])
    widths  = (panel_size_in[0] - heights[0], heights[0])  # legend column narrower than scatter

    fig = plt.figure(figsize=panel_size_in, dpi=dpi)
    gs  = gridspec.GridSpec(2, 2, fig, 0,0,1,1,
                            hspace=0, wspace=0,
                            width_ratios=widths, height_ratios=heights)

    # ===== left column: title + legend (with colored dots) =====
    ax = fig.add_subplot(gs[0,0], xticks=[], yticks=[], xlim=[0,1], ylim=[0,1], frameon=False)
    li = 1
    yl = lambda i: 1 - i/10.5
    xl_dot, xl_head, xl_text = 0.22, 0.30, 0.30

    # Title ("Stage X" in Veres; we show stage_text and n_cells)
    if stage_text is None:
        stage_text = ""
    ax.text(xl_head, yl(li), f"{stage_text}", va='center', fontsize=6, fontweight='extra bold', clip_on=False)
    li += 1
    ax.text(xl_head, yl(li), f"({adata.n_obs:,} cells)", va='center', fontsize=5, clip_on=False)
    li += 1.5

    present = [c for c in (label_order or cats) if c in set(labels)]
    for lb in present:
        lb_txt = lp[lb]['short_label']
        ax.scatter(xl_dot, yl(li)+0.008, s=15, c=lp[lb]['color'].reshape(1,-1), clip_on=False)
        for line in str(lb_txt).splitlines():
            ax.text(xl_text, yl(li), line, va='center', fontsize=5, clip_on=False)
            li += 1

    # ===== right column: layered scatter (black -> white -> color), alpha as in notebook =====
    proj, rgb = prepare_for_scatter(X, labels, lp)
    ax = fig.add_subplot(gs[0,1], xticks=[], yticks=[], frameon=False, zorder=-1)
    ax.patch.set_visible(False)

    s_black, s_white, s_type = 4, 2, 1.5   # EXACT dot sizes used in the notebook
    ax.scatter(proj[:,0], proj[:,1], c='k', edgecolor='none', s=s_black, rasterized=True)
    ax.scatter(proj[:,0], proj[:,1], c='w', edgecolor='none', s=s_white, rasterized=True)
    ax.scatter(proj[:,0], proj[:,1], c=rgb, edgecolor='none', s=s_type, alpha=0.7, rasterized=True)

    # ===== bottom row: population ratio bar (with border rectangle and label) =====
    ax = fig.add_subplot(gs[1,0:2], xticks=[], yticks=[], xlim=[0,1], ylim=[0,1], frameon=False)

    counts = pd.Series(Counter(labels))
    fracs  = counts / counts.sum()
    rr = ratio_order or present  # if no ratio_order, follow legend order

    cumul = 0.0
    for lb in rr:
        if lb in fracs:
            w = float(fracs[lb])
            ax.add_patch(patches.Rectangle((0.30 + 0.68*cumul, 0.02), 0.68*w, 1.0,
                                           facecolor=lp[lb]['color'], edgecolor='none', clip_on=False))
            cumul += w

    ax.text(0.175, 0.5, 'Pop. ratios:', va='center', ha='center', fontsize=5, clip_on=False)
    ax.add_patch(patches.Rectangle((0.30, 0.02), 0.68, 1.0, facecolor='none', edgecolor='k', linewidth=0.5, clip_on=False))

    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=dpi, transparent=True, bbox_inches='tight')
    return fig


def plot_veres_panel_highlight(
    adata,
    label_key,
    highlight_label,
    stage_text=None,
    palette=None,              # dict {label: hex_color} e.g. celltype_palette
    save=None,
    dpi=600,
    bg_color=np.array([0.8, 0.8, 0.8]),   # soft grey background
):
    import pandas as pd
    """Single-label highlight using the Veres-style panel structure."""
    # --- setup data ---
    X = adata.obsm['X_umap']
    labels = adata.obs[label_key].astype(str).values
    cats = list(pd.unique(labels))
    lp = make_label_params(cats, palette=palette)

    # --- Veres figure geometry ---
    mm_per_inch = 25.4
    panel_size_in = ((89/2)/mm_per_inch, 25.4/mm_per_inch)
    cell_pop_bar_h = 0.07
    heights = ((1 - cell_pop_bar_h)*panel_size_in[1], cell_pop_bar_h * panel_size_in[1])
    widths = (panel_size_in[0] - heights[0], heights[0])

    fig = plt.figure(figsize=panel_size_in, dpi=dpi)
    gs = gridspec.GridSpec(2, 2, fig, 0, 0, 1, 1,
                           hspace=0, wspace=0,
                           width_ratios=widths, height_ratios=heights)

    # === left column: stage + highlight label ===
    ax = fig.add_subplot(gs[0, 0], xticks=[], yticks=[], xlim=[0, 1], ylim=[0, 1], frameon=False)
    yl = lambda i: 1 - i / 10.5
    xl_dot, xl_head, xl_text = 0.22, 0.30, 0.30
    li = 1
    if stage_text:
        ax.text(xl_head, yl(li), stage_text, va="center", fontsize=6, fontweight="extra bold", clip_on=False)
        li += 1
    ax.text(xl_head, yl(li), f"({adata.n_obs:,} cells)", va="center", fontsize=5, clip_on=False)
    li += 1.5

    # add single label color to legend
    ax.scatter(xl_dot, yl(li) + 0.008, s=15, c=lp[highlight_label]['color'].reshape(1, -1), clip_on=False)
    ax.text(xl_text, yl(li), highlight_label, va="center", fontsize=5, clip_on=False)

    # === right column: UMAP (grey background + color highlight) ===
    # Sort by neighborhood homogeneity (same logic as prepare_for_scatter)
    # but keep labels aligned with sorted coordinates
    labels_arr = np.asarray(labels)
    nbrs = NearestNeighbors(n_neighbors=5, algorithm='ball_tree').fit(X)
    knn = nbrs.kneighbors(X, return_distance=False)
    same = (labels_arr[knn].T == labels_arr[knn[:, 0]]).mean(0)
    order = np.argsort(same)
    proj = X[order]
    sorted_labels = labels_arr[order]

    ax = fig.add_subplot(gs[0, 1], xticks=[], yticks=[], frameon=False)
    ax.patch.set_visible(False)

    # masks on the sorted arrays
    mask_high = sorted_labels == highlight_label
    mask_bg = ~mask_high

    s_black, s_white, s_type = 4, 2, 1.5
    # background in light grey (instead of color)
    ax.scatter(proj[mask_bg, 0], proj[mask_bg, 1], c='k', s=s_black, edgecolor='none', alpha=0.05, rasterized=True)
    ax.scatter(proj[mask_bg, 0], proj[mask_bg, 1], c='w', s=s_white, edgecolor='none', alpha=0.05, rasterized=True)
    ax.scatter(proj[mask_bg, 0], proj[mask_bg, 1], c=[bg_color], s=s_type, edgecolor='none', alpha=0.2, rasterized=True)

    # highlighted label in its assigned color
    highlight_color = lp[highlight_label]['color']
    ax.scatter(proj[mask_high, 0], proj[mask_high, 1],
               c=highlight_color.reshape(1, -1), s=s_type, edgecolor='none', alpha=0.9, rasterized=True)

    # === bottom row: population ratio bar (still full for context) ===
    ax = fig.add_subplot(gs[1, 0:2], xticks=[], yticks=[], xlim=[0, 1], ylim=[0, 1], frameon=False)
    from collections import Counter
    import pandas as pd
    counts = pd.Series(Counter(labels))
    fracs = counts / counts.sum()
    cumul = 0
    for lb, frac in fracs.items():
        color = lp[lb]['color'] if lb == highlight_label else bg_color
        ax.add_patch(patches.Rectangle((0.30 + 0.68 * cumul, 0.02), 0.68 * frac, 1.0,
                                       facecolor=color, edgecolor='none', clip_on=False))
        cumul += frac

    ax.text(0.175, 0.5, "Pop. ratios:", va="center", ha="center", fontsize=5, clip_on=False)
    ax.add_patch(patches.Rectangle((0.30, 0.02), 0.68, 1.0, facecolor='none', edgecolor='k', linewidth=0.5, clip_on=False))

    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=dpi, transparent=True, bbox_inches='tight')
    return fig


# -------------------- marker UMAP panels (from 3_Marker_plotting.ipynb) --------------------

def plot_marker_umaps(adata, markers, ncols=4, cmap="inferno", vmax="p95",
                      save=None, show=True):
    """
    Plot per-gene UMAP panels for a set of markers.

    Parameters
    ----------
    adata : AnnData
    markers : dict or list
        If dict {cell_type: [genes]}, genes are flattened (unique).
        If list, used directly.
    ncols : int
    cmap : str
    vmax : str or float
    save : str or Path, optional
    show : bool
    """
    if isinstance(markers, dict):
        filtered, _ = filter_marker_dict(adata, markers)
        genes = []
        for v in filtered.values():
            genes.extend(v)
        genes = list(dict.fromkeys(genes))  # unique, order-preserved
    else:
        varset = set(adata.var_names)
        genes = [g for g in markers if g in varset]

    if not genes:
        print("No marker genes found in dataset.")
        return

    sc.pl.umap(
        adata, color=genes, ncols=ncols, cmap=cmap,
        vmax=vmax, frameon=False, show=show, save=save,
    )
