#!/bin/bash
# Step 1 (unified labelling) — sbatch wrapper for 20_integrate_timepoint.py
#
# Subset integrated.h5ad to one timepoint (control + Lapa), scVI batch-correct
# on participant, Leiden + UMAP on the integrated latent. CPU is fine; a GPU
# node cuts the 200-epoch scVI train substantially.
#
# Pick the timepoint with TIMEPOINT=D2 / D4. Submit one job per timepoint:
#     mkdir -p logs
#     sbatch --export=ALL,TIMEPOINT=D2 scripts/hpc-smoke/20_integrate_timepoint.sbatch
#     sbatch --export=ALL,TIMEPOINT=D4 scripts/hpc-smoke/20_integrate_timepoint.sbatch
#
# Smoke first (few epochs) to confirm it runs end-to-end:
#     sbatch --export=ALL,TIMEPOINT=D2,SCVI_EPOCHS=10 scripts/hpc-smoke/20_integrate_timepoint.sbatch

#SBATCH --job-name=unified-integrate
#SBATCH --partition=bch-compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

VENV_PY=$HOME/venvs/single-cell-gpu/bin/python
PROJECT_ROOT=${SLURM_SUBMIT_DIR:-$HOME/single-cell}
SCRIPT=${PROJECT_ROOT}/scripts/hpc-smoke/20_integrate_timepoint.py

cd "${PROJECT_ROOT}"
echo "host=$(hostname)  job=${SLURM_JOB_ID:-?}  partition=${SLURM_JOB_PARTITION:-?}"
echo "cwd=$(pwd)  TIMEPOINT=${TIMEPOINT:-unset}"
echo "python: ${VENV_PY}"

PYTHONUNBUFFERED=1 ${VENV_PY} -u "${SCRIPT}" "${TIMEPOINT:?set TIMEPOINT=D2 or D4}"

echo "=== done ==="
