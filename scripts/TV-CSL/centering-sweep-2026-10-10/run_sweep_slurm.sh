#!/bin/bash
# SLURM array for the centering sweep: one task = one (rep i, effect scale m) job, all three centerings.
# 50 reps x 5 values of m = 250 tasks; task id t -> i = (t-1) %/% 5 + 1, m = M[(t-1) %% 5].
# Resumable: a task whose result file already exists exits immediately, so resubmitting fills gaps.
#
# Usage (from the project root, on the fix-dgp-exact-event-times branch):
#   sbatch scripts/TV-CSL/centering-sweep-2026-10-10/run_sweep_slurm.sh
# Throttle concurrency if needed, e.g. at most 100 at once:
#   sbatch --array=1-250%100 scripts/TV-CSL/centering-sweep-2026-10-10/run_sweep_slurm.sh

#SBATCH -J tvcsl_sweep
#SBATCH -o scripts/TV-CSL/centering-sweep-2026-10-10/logs/slurm_%A_%a.out
#SBATCH -e scripts/TV-CSL/centering-sweep-2026-10-10/logs/slurm_%A_%a.err
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 03:00:00
#SBATCH --array=1-250

project_dir="${PROJECT_DIR:-/homes2/xmeng/TV-CSL}"
D="${project_dir}/scripts/TV-CSL/centering-sweep-2026-10-10"
N=4000
M=(0 0.5 1 1.5 2)

t=${SLURM_ARRAY_TASK_ID}
i=$(( (t - 1) / ${#M[@]} + 1 ))
m=${M[$(( (t - 1) % ${#M[@]} ))]}

module load R/4.3.2
# One BLAS thread per task: multithreaded OpenBLAS oversubscribes the node (~3x slower).
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

cd "${project_dir}"
mkdir -p "${D}/logs" "${D}/res"

out="${D}/res/rep_n${N}_m${m}_i${i}.csv"
if [ -f "${out}" ]; then
  echo "Task ${t}: ${out} exists, skipping."
  exit 0
fi

echo "Task ${t}: i=${i} m=${m} n=${N} on $(hostname)"
Rscript "${D}/sweep_rep.R" "${i}" "${m}" "${N}"
rc=$?
echo "Task ${t}: exit ${rc}"
exit ${rc}
