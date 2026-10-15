#!/bin/bash
# Pre-generates all datasets (one short job), then submits the SLURM array running all methods
# (main + oracle) once it succeeds. Writes to scripts/TV-CSL/results/ as specified in params-main-methods.json.
#
# Usage: bash scripts/TV-CSL/run-all.sh
#
# This must be run from the project root directory.

mkdir -p scripts/TV-CSL/logs
pregen=$(sbatch --parsable --job-name=TV_CSL_pregen -c 1 --mem=4G -t 02:00:00 \
  -o scripts/TV-CSL/logs/%x_%j.out -e scripts/TV-CSL/logs/%x_%j.err \
  --wrap="module load R/4.3.2; export OPENBLAS_NUM_THREADS=1; Rscript scripts/TV-CSL/pregenerate-data.R scripts/TV-CSL/params-main-methods.json")
echo "Submitted data pre-generation job ${pregen}"

sbatch --dependency=afterok:${pregen} --job-name=TV_CSL scripts/TV-CSL/run-cluster.sh scripts/TV-CSL/params-main-methods.json
