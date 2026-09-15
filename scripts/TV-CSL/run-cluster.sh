#!/bin/bash
# Generic SLURM array runner for TV-CSL simulations.
#
# Usage:
#   sbatch [--job-name=NAME] run-cluster.sh <relative/path/to/params.json>
#
# Examples:
#   sbatch --job-name=TV_CSL_main   run-cluster.sh scripts/TV-CSL/params-main-methods.json
#   sbatch --job-name=TV_CSL_oracle run-cluster.sh scripts/TV-CSL/params-oracle-only.json
#
# The job name appears in log filenames via %x. Results directory is read from the JSON file.

#SBATCH -J TV_CSL_run
#SBATCH -o scripts/TV-CSL/logs/%x_%A_%a.out
#SBATCH -e scripts/TV-CSL/logs/%x_%A_%a.err
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 24:00:00
#SBATCH --array=1-1000

project_dir="/homes2/xmeng/TV-CSL"
json_file="${project_dir}/${1}"

module load R/4.3.2

mkdir -p "${project_dir}/scripts/TV-CSL/logs"
cd "${project_dir}"

echo "Config:    ${json_file}"
echo "Iteration: ${SLURM_ARRAY_TASK_ID}"

for n in 200 500 1000 2000; do
  for eta_type in "linear" "non-linear"; do
    echo "  n=${n}  eta_type=${eta_type}"
    Rscript -e "
      source('scripts/TV-CSL/TV-CSL-runner.R')
      run_experiment_iteration(
        i         = ${SLURM_ARRAY_TASK_ID},
        json_file = '${json_file}',
        eta_type  = '${eta_type}',
        HTE_type  = 'linear',
        n         = ${n},
        verbose   = 1
      )
    "
  done
done

echo "Finished iteration ${SLURM_ARRAY_TASK_ID}"
