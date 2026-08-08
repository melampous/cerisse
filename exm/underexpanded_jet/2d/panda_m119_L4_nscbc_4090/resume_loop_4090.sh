#!/bin/bash
# Auto-resume harness: launch solver from latest chk; on guard trip, resume
# (max 3 resumes). Emits one status line per event (consumed by Monitor).
D="$(cd "$(dirname "$0")" && pwd)"
cd "$D" || exit 1
rm -f RUN_DONE
MAXR=3
LEG=0
while :; do
  LEG=$((LEG+1))
  CHK=$(ls -d plot/chk????? 2>/dev/null | sort -V | tail -1)
  rm -f dump_and_stop
  RST=""
  [ -n "$CHK" ] && RST="amr.restart=$CHK"
  echo "leg$LEG START restart=${CHK:-cold} $(date -u +%FT%TZ)"
  setsid nice -n 10 ./main2d.gnu.MPI.CUDA.ex inputs $RST \
    amr.regrid_int=40 amr.check_int=1000 amr.plot_int=500 amr.data_log=probes.log \
    > "run4090_leg$LEG.log" 2>&1 < /dev/null &
  SPID=$!
  setsid ./guard_4090.sh < /dev/null > /dev/null 2>&1 &
  setsid ./perf_logger_leg.sh "run4090_leg$LEG.log" < /dev/null > /dev/null 2>&1 &
  while kill -0 "$SPID" 2>/dev/null; do sleep 60; done
  sleep 15
  LASTT=$(grep -E "^STEP = " "run4090_leg$LEG.log" | tail -1 | awk '{print $6}')
  if grep -q "finalized" "run4090_leg$LEG.log" && awk -v t="${LASTT:-0}" 'BEGIN{exit !(t>=9.99e-3)}'; then
    echo "FINISHED_10MS t=$LASTT $(date -u +%FT%TZ)"; touch RUN_DONE; exit 0
  fi
  if grep -q "LIMIT HIT" guard.log; then
    if [ "$LEG" -le "$MAXR" ]; then
      echo "GUARD_TRIP_RESUME t=$LASTT leg$LEG->leg$((LEG+1)) $(date -u +%FT%TZ)"
      sleep 30; continue
    fi
    echo "GUARD_TRIP_RETRIES_EXHAUSTED t=$LASTT $(date -u +%FT%TZ)"; touch RUN_DONE; exit 1
  fi
  echo "DIED_UNEXPECTED t=$LASTT (see run4090_leg$LEG.log) $(date -u +%FT%TZ)"
  touch RUN_DONE; exit 2
done
