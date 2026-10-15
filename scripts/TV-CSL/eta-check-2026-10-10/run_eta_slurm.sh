#!/bin/bash
# TV-CSL with the true eta, m = 2, n = 4000: 20 reps x 3 centerings = 60 tasks (~20 min each).
# Usage (from the project root):
#   sbatch scripts/TV-CSL/eta-check-2026-10-10/run_eta_slurm.sh

#SBATCH -J tvcsl_etacheck
#SBATCH -o scripts/TV-CSL/eta-check-2026-10-10/logs/slurm_%A_%a.out
#SBATCH -e scripts/TV-CSL/eta-check-2026-10-10/logs/slurm_%A_%a.err
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 02:00:00
#SBATCH --array=1-60

project_dir="${PROJECT_DIR:-/homes2/xmeng/TV-CSL}"
SPECS=(cox-linear-censored-only cox-time-varying-prop cox-time-varying-oracle)
t=${SLURM_ARRAY_TASK_ID}
i=$(( (t - 1) / ${#SPECS[@]} + 1 ))
ps=${SPECS[$(( (t - 1) % ${#SPECS[@]} ))]}

module load R/4.3.2
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

cd "${project_dir}"
echo "Task ${t}: i=${i} spec=${ps} on $(hostname)"
Rscript scripts/TV-CSL/eta-check-2026-10-10/eta_rep.R "${i}" "${ps}"
rc=$?
echo "Task ${t}: exit ${rc}"
exit ${rc}
