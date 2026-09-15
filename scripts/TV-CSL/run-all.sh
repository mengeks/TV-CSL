#!/bin/bash
# Submits both SLURM array jobs to fully replicate all simulation results.
# Both jobs write to scripts/TV-CSL/results/ (as specified in their JSON params files).
#
# Usage: bash scripts/TV-CSL/run-all.sh
#
# This must be run from the project root directory.

sbatch --job-name=TV_CSL_main   scripts/TV-CSL/run-cluster.sh scripts/TV-CSL/params-main-methods.json
sbatch --job-name=TV_CSL_oracle scripts/TV-CSL/run-cluster.sh scripts/TV-CSL/params-oracle-only.json
