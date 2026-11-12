#!/bin/bash
# Submits a single SLURM array job running all methods (main + oracle).
# Writes to scripts/TV-CSL/results/ as specified in params-main-methods.json.
#
# Our QOS allows one submitted job at a time, so pre-generate the datasets first, outside SLURM
# (about 15 min; it avoids tasks racing on the shared i + 100 test files):
#   Rscript scripts/TV-CSL/pregenerate-data.R scripts/TV-CSL/params-main-methods.json
#
# Usage: bash scripts/TV-CSL/run-all.sh
#
# This must be run from the project root directory.

sbatch --job-name=TV_CSL scripts/TV-CSL/run-cluster.sh scripts/TV-CSL/params-main-methods.json
