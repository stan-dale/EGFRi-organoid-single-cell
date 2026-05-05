#!/bin/bash
#SBATCH --job-name=egfri-d2d4-scanvi
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Runs notebooks 0-5 sequentially via papermill on a GPU node.
# Edit the FILL-IN variables below for your cluster.
#
# Usage from the repo root:
#     mkdir -p logs
#     sbatch notebooks/d2-d4-integrated/hpc/run_pipeline.sbatch
#
# To run only the heavy notebooks (0/1/2) on GPU and the rest locally,
# set RUN_PLOTS=0 below.

set -euo pipefail

# ─── FILL IN: how to activate the env that has scvi-tools ────────────────
# Pick ONE of the two blocks below.
#
# A) Conda / mamba env name:
#     source ~/.bashrc                       # if needed for `conda` activation
#     conda activate single-cell-env
#
# B) Module load (if your cluster uses Lmod):
#     module load python/3.11 cuda/12.1
#     source /path/to/single-cell-env/bin/activate
#
# Use the helper scripts in the same directory to identify the right
# activation:  notebooks/d2-d4-integrated/hpc/discover_env.sh
ENV_ACTIVATE='conda activate single-cell-env'   # <-- EDIT ME
eval "${ENV_ACTIVATE}"

# ─── FILL IN: project root on the cluster ────────────────────────────────
PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
cd "${PROJECT_ROOT}"

# Optional: choose whether to run the plotting notebooks here too.
RUN_PLOTS=${RUN_PLOTS:-1}

NB_DIR="notebooks/d2-d4-integrated"
OUT_DIR="${NB_DIR}/_executed"
mkdir -p "${OUT_DIR}"

run_nb () {
    local nb="$1"
    local in_path="${NB_DIR}/${nb}"
    local out_path="${OUT_DIR}/${nb%.ipynb}.executed.ipynb"
    echo "===> papermill ${in_path}"
    papermill "${in_path}" "${out_path}" --kernel python3
}

echo "node=$(hostname) gpus=${CUDA_VISIBLE_DEVICES:-none}"
nvidia-smi || true

run_nb 0_Prepare_Slim_h5ads.ipynb
run_nb 1_Reference_HVGs_and_Training.ipynb
run_nb 2_Map_Queries.ipynb

if [[ "${RUN_PLOTS}" -eq 1 ]]; then
    run_nb 3_Sanity_CrossTab_vs_CellAssign.ipynb
    run_nb 4_Plotting_UMAP_Panels.ipynb
    run_nb 5_Population_Ratios.ipynb
fi

echo "done."
