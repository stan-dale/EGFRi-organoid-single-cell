# HPC execution

This folder lets you run the D2/D4 transfer-learning pipeline on a SLURM cluster with a GPU. There are exactly two files to know about:

- `discover_env.sh` — run once on the cluster to identify how scANVI is installed.
- `run_pipeline.sbatch` — submits the pipeline (notebooks 0–5 via papermill).

## 1. First-time setup

From the cluster login node, in the project root:

```bash
git fetch origin
git checkout claude/d2-d4-transfer-learning-analysis-GVjst
bash notebooks/d2-d4-integrated/hpc/discover_env.sh
```

The script prints which conda envs (if any) contain `scvi`, whether `module` is available, your current `python` version, GPU visibility, and likely GPU SLURM partitions.

Pick the activation line that matches your cluster and edit `ENV_ACTIVATE` near the top of `run_pipeline.sbatch`. Common choices:

```bash
# (A) conda env on the cluster
ENV_ACTIVATE='conda activate single-cell-env'

# (B) module-based stack
ENV_ACTIVATE='module load python/3.11 cuda/12.1 && source /path/to/single-cell-env/bin/activate'

# (C) pre-packed tarball (per notes/WORKFLOW.md §8)
ENV_ACTIVATE='source ~/envs/single-cell-env/bin/activate && conda-unpack'
```

Also confirm `--partition`, `--gres`, `--time` match your cluster's GPU queue. The defaults assume a partition called `gpu` with a single GPU available via `--gres=gpu:1`.

## 2. Data layout on the cluster

The pipeline reads from `data/data-objects/{clustered,cellassign,labelled}/` and writes back to `data/data-objects/integration/`, `data/data-objects/cellassign/`, and `data/analysis-files/cell_assign_dfs/`. Same paths as on the desktop — the project's `src/config.py` handles them.

Make sure the input h5ads exist on the cluster before submitting:

```bash
ls data/data-objects/cellassign/d10_lapa_predictions.h5ad
ls data/data-objects/clustered/clustered_egfDuod_D{2,4}_{DZ,Lapa}_DZ.h5ad
```

If they're not there, `rsync` from your desktop or the lab's storage.

## 3. Submit

```bash
mkdir -p logs
sbatch notebooks/d2-d4-integrated/hpc/run_pipeline.sbatch
```

Watch progress:

```bash
squeue -u "${USER}"
tail -f logs/egfri-d2d4-scanvi-*.out
```

## 4. What runs

The job uses `papermill` to execute the six notebooks in order, writing executed copies to `notebooks/d2-d4-integrated/_executed/`:

1. `0_Prepare_Slim_h5ads.ipynb` — streams each dataset, writes slim h5ads. Light I/O.
2. `1_Reference_HVGs_and_Training.ipynb` — picks 3000 HVGs from D10, trains SCVI then SCANVI on D10. **GPU-bound.**
3. `2_Map_Queries.ipynb` — for each of D2 and D4: concat Dz + Lapa slim files, `prepare_query_anndata` + `load_query_data` + short fine-tune, predict labels with `Ambiguous_<argmax>` flagging at confidence < 0.5. **GPU-bound.**
4. `3_Sanity_CrossTab_vs_CellAssign.ipynb` — confusion matrices vs the existing per-dataset CellAssign and marker-gene dotplots.
5. `4_Plotting_UMAP_Panels.ipynb` — per-timepoint UMAPs (labels, treatment, confidence) plus rare-class highlights.
6. `5_Population_Ratios.ipynb` — stacked-bar composition by `(Time_point × Treatment)` and per-participant.

Set `RUN_PLOTS=0` in the env (or edit the script) to skip notebooks 3–5 on the cluster and run them locally.

## 5. Outputs to copy back

After the job finishes, these are the files to `rsync` back to your desktop:

```
data/data-objects/cellassign/d2_scanvi_predictions.h5ad
data/data-objects/cellassign/d4_scanvi_predictions.h5ad
data/data-objects/integration/scarches_model/
data/analysis-files/cell_assign_dfs/d{2,4}_scanvi_probabilities.csv
figures/d2-d4-integrated/
notebooks/d2-d4-integrated/_executed/
```

## 6. Troubleshooting

- **`scvi.ImportError` in notebook 1**: `ENV_ACTIVATE` didn't yield an env with scvi-tools. Re-run `discover_env.sh` and try a different activation.
- **`CUDA out of memory`**: drop `n_latent` to 20 in notebook 1 or reduce `batch_size`. With 3000 HVGs, this should not happen on a 16+ GB GPU.
- **Reference recall < 95% in notebook 1**: SCANVI failed to recover the input labels. Inspect `cellassign_confidence` distribution and try setting `REFERENCE_CONFIDENCE_THRESHOLD = 0.85` near the top of notebook 1.
- **Query genes missing**: notebook 2 prints `n_missing_genes`; if > 10% of HVGs are missing for a query, check the slim file's gene set.
