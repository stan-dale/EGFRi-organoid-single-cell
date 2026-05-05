#!/bin/bash
# Run on the cluster login node (or in an interactive job) to figure out
# how scANVI / scvi-tools is installed.
#
#     bash notebooks/d2-d4-integrated/hpc/discover_env.sh
#
# The output tells you what to put in ENV_ACTIVATE in run_pipeline.sbatch.

set -u

echo "=== A) Conda / mamba ==="
if command -v conda >/dev/null 2>&1; then
    echo "conda found: $(conda --version)"
    echo "envs containing 'scvi' or 'cell':"
    conda env list 2>/dev/null | awk 'NR>2 && $1 !~ /^#/ {print $1}' | while read -r env; do
        if [[ -n "${env}" ]]; then
            if conda run -n "${env}" python -c "import scvi" 2>/dev/null; then
                ver=$(conda run -n "${env}" python -c "import scvi; print(scvi.__version__)")
                echo "  ${env}    (scvi=${ver})"
            fi
        fi
    done
else
    echo "no conda on PATH"
fi

echo
echo "=== B) Lmod / Environment Modules ==="
if command -v module >/dev/null 2>&1; then
    echo "module command available"
    module avail 2>&1 | grep -iE 'python|conda|scvi|cuda|cell' | head -40 || echo "(no matches)"
else
    echo "no module command on PATH"
fi

echo
echo "=== C) python on PATH right now ==="
python --version 2>&1 || echo "no python"
python -c "import scvi; print('scvi=' + scvi.__version__); from scvi.model import SCANVI; print('SCANVI available')" 2>&1 || echo "scvi NOT importable from current python"

echo
echo "=== D) GPU visibility ==="
nvidia-smi 2>&1 | head -20 || echo "no GPU visible from login node (try inside an srun/sbatch job)"

echo
echo "=== E) Suggested SLURM partitions (from sinfo) ==="
sinfo -h -o '%R %G' 2>&1 | grep -iE 'gpu|a100|v100|h100' || echo "(no obvious GPU partitions; try 'sinfo' alone)"
