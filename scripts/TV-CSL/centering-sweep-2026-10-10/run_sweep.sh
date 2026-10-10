#!/bin/bash
# Centering sweep: n = 4000, m in {0, 0.5, 1, 1.5, 2}, 50 reps, all three centerings.
# Ordered rep-major so every m has the same number of reps if stopped early.
# Resumable: existing rep files are skipped. Run from project root:
#   setsid nohup bash scripts/TV-CSL/centering-sweep-2026-10-10/run_sweep.sh > scripts/TV-CSL/centering-sweep-2026-10-10/logs/launcher.log 2>&1 < /dev/null &
D=/homes2/xmeng/TV-CSL/scripts/TV-CSL/centering-sweep-2026-10-10
N=${N:-4000}; R=${R:-50}; NPAR=${NPAR:-40}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1   # multithreaded BLAS oversubscribes badly
cd /homes2/xmeng/TV-CSL
for i in $(seq 1 "$R"); do for m in 0 0.5 1 1.5 2; do
  [ -f "$D/res/rep_n${N}_m${m}_i${i}.csv" ] || echo "$i $m"
done; done |
  xargs -P "$NPAR" -L 1 bash -c 'Rscript '"$D"'/sweep_rep.R "$0" "$1" '"$N"' > '"$D"'/logs/rep_n'"$N"'_m"$1"_i"$0".log 2>&1; echo "$(date +%T) finished m=$1 i=$0 (exit $?)"'
echo "$(date) ALL DONE"
