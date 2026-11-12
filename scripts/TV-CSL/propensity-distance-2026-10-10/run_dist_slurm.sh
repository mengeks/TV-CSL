#!/bin/bash
# Propensity distance check: one task per effect scale m, 20 reps each (~20 s per rep, ~7 min per task).
# Usage (from the project root):
#   sbatch scripts/TV-CSL/propensity-distance-2026-10-10/run_dist_slurm.sh

#SBATCH -J tvcsl_propdist
#SBATCH -o scripts/TV-CSL/propensity-distance-2026-10-10/logs/slurm_%A_%a.out
#SBATCH -e scripts/TV-CSL/propensity-distance-2026-10-10/logs/slurm_%A_%a.err
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 01:00:00
#SBATCH --array=0-4

project_dir="${PROJECT_DIR:-/homes2/xmeng/TV-CSL}"
M=(0 0.5 1 1.5 2)
m=${M[${SLURM_ARRAY_TASK_ID}]}

module load R/4.3.2
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

cd "${project_dir}"
mkdir -p scripts/TV-CSL/propensity-distance-2026-10-10/logs scripts/TV-CSL/propensity-distance-2026-10-10/res
echo "Task ${SLURM_ARRAY_TASK_ID}: m=${m} on $(hostname)"
Rscript scripts/TV-CSL/propensity-distance-2026-10-10/dist_task.R "${m}" 20 4000
rc=$?
echo "Task ${SLURM_ARRAY_TASK_ID}: exit ${rc}"
exit ${rc}
