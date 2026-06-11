#!/bin/bash
# Step 1b — sbatch wrapper for 21_unified_umaps.py (integration QC).
#
# ALWAYS run h5ad-loading scripts on a compute node, never the login node.
#     mkdir -p logs
#     sbatch --export=ALL,TIMEPOINT=D2 scripts/hpc-smoke/21_unified_umaps.sbatch
#     sbatch --export=ALL,TIMEPOINT=D4 scripts/hpc-smoke/21_unified_umaps.sbatch
#
# Light job: load one unified h5ad, print the batch diagnostic, write 2 PNGs.

#SBATCH --job-name=unified-umaps
#SBATCH --partition=bch-compute
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

VENV_PY=$HOME/venvs/single-cell-gpu/bin/python
PROJECT_ROOT=${SLURM_SUBMIT_DIR:-$HOME/single-cell}
SCRIPT=${PROJECT_ROOT}/scripts/hpc-smoke/21_unified_umaps.py

cd "${PROJECT_ROOT}"
echo "host=$(hostname)  job=${SLURM_JOB_ID:-?}  TIMEPOINT=${TIMEPOINT:-unset}"

PYTHONUNBUFFERED=1 ${VENV_PY} -u "${SCRIPT}" "${TIMEPOINT:?set TIMEPOINT=D2 or D4}"

echo "=== done ==="
