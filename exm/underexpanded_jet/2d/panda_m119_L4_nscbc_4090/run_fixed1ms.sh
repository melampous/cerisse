#!/bin/bash
# Fixed-dt extension: restart from the 10 ms checkpoint and run to 11 ms with
# cns.dt_max=4.0e-7 s (below the historical CFL-dt minimum 4.215e-7 over the
# 6-10 ms window), so dt is constant and probe records are uniformly sampled
# (fs = 2.5 MHz) for spectral analysis. Probes go to a separate file
# probes_fixed.log. Auto-resume on guard trip (max 3), same as production.
D="$(cd "$(dirname "$0")" && pwd)"
cd "$D" || exit 1
rm -f RUN_DONE
MAXR=3
LEG=0
while :; do
  LEG=$((LEG+1))
  CHK=$(ls -d plot/chk????? 2>/dev/null | sort -V | tail -1)
  rm -f dump_and_stop
  echo "fxleg$LEG START restart=$CHK $(date -u +%FT%TZ)"
  setsid nice -n 10 ./main2d.gnu.MPI.CUDA.ex inputs amr.restart=$CHK \
    stop_time=11.0e-3 cns.dt_max=4.0e-7 \
    amr.regrid_int=40 amr.check_int=1000 amr.plot_int=500 amr.data_log=probes_fixed.log \
    > "run4090_fx_leg$LEG.log" 2>&1 < /dev/null &
  SPID=$!
  setsid ./guard_4090.sh < /dev/null > /dev/null 2>&1 &
  setsid ./perf_logger_leg.sh "run4090_fx_leg$LEG.log" 11.0e-3 < /dev/null > /dev/null 2>&1 &
  while kill -0 "$SPID" 2>/dev/null; do sleep 60; done
  sleep 15
  LASTT=$(grep -E "^STEP = " "run4090_fx_leg$LEG.log" | tail -1 | awk '{print $6}')
  if grep -q "finalized" "run4090_fx_leg$LEG.log" && awk -v t="${LASTT:-0}" 'BEGIN{exit !(t>=10.99e-3)}'; then
    echo "FINISHED_11MS t=$LASTT $(date -u +%FT%TZ)"; touch RUN_DONE; exit 0
  fi
  if grep -q "LIMIT HIT" guard.log; then
    if [ "$LEG" -le "$MAXR" ]; then
      echo "GUARD_TRIP_RESUME t=$LASTT fxleg$LEG->fxleg$((LEG+1)) $(date -u +%FT%TZ)"
      sleep 30; continue
    fi
    echo "GUARD_TRIP_RETRIES_EXHAUSTED t=$LASTT $(date -u +%FT%TZ)"; touch RUN_DONE; exit 1
  fi
  echo "DIED_UNEXPECTED t=$LASTT (see run4090_fx_leg$LEG.log) $(date -u +%FT%TZ)"
  touch RUN_DONE; exit 2
done
