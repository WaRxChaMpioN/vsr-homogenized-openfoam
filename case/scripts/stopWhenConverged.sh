#!/bin/bash
# Ends heat-transfer runs once converged instead of at a fixed endTime: when the last two
# 'energyAudit global' residuals (every 100 iterations) are both below TOL % of the beam and at least
# MIN iterations are done, endTime is set to the next write time (runTimeModifiable re-reads
# controlDict). Allrun.matrix then writes the drift of T over the last write interval (log.drift).
# Run in the run directory of Allrun.matrix; it stops when the matrix has finished:
#   $VSR/case/scripts/stopWhenConverged.sh > log.stopper 2>&1 &     (TOL=0.02 MIN=3000 by default)
TOL=${TOL:-0.02}; MIN=${MIN:-3000}
while :
do
    for c in G*_Re*_*x*
    do
        L=$c/log.chtMultiRegionSimpleFoam; D=$c/system/controlDict
        [ -f "$L" ] && [ -f "$D" ] && [ ! -f "$c/log.energyBalance" ] || continue      # running only
        end=$(awk '/^endTime/{gsub(";",""); print $2}' "$D")
        it=$(grep '^Time = ' "$L" | tail -1 | awk '{print $3}')
        [ -n "$it" ] && [ "$it" -ge "$MIN" ] || continue
        next=$(( (it / 1000 + 1) * 1000 ))
        [ "$next" -lt "$end" ] || continue
        if grep '^energyAudit global' "$L" | tail -2 | sed 's/.*(\(.*\) % of beam)/\1/' \
            | awk -v t="$TOL" '{a = ($1 < 0 ? -$1 : $1); if (a < t) n++} END {exit !(n == 2)}'
        then
            sed -i "s/^endTime .*/endTime         $next;/" "$D"
            echo "$(date '+%F %T') $c: converged at iteration $it, endTime $end -> $next"
        fi
    done
    grep -q 'matrix finished' log.matrix 2>/dev/null && break
    [ -f log.matrix ] && ! pgrep -f 'Allrun.matrix' > /dev/null && break        # matrix stopped early
    sleep 60
done
echo "$(date '+%F %T') matrix finished, stopper done"
