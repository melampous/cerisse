#!/bin/bash
# 5-min telemetry: step rate, physical progress, ETA, GPU/host resources.
D="$(cd "$(dirname "$0")" && pwd)"
OUT="$D/perf4090.csv"; RLOG="${1:-run4090b.log}"; TSTOP="${2:-10.0e-3}"
[ -s "$OUT" ] || echo "utc,step,t_ms,dt,s_per_step,eta_h,gpu_T,gpu_W,vram_MiB,util,freeRAM_MiB" > "$OUT"
PREV_STEP=0; PREV_EPOCH=$(date +%s)
while pgrep -f "main2d.gnu.MPI.CUDA.e[x]" > /dev/null; do
  sleep 300
  L=$(grep -E "^STEP = " "$D/$RLOG" 2>/dev/null | tail -1)
  STEP=$(echo "$L" | awk '{print $3}')
  T=$(echo "$L" | awk '{print $6}')
  DT=$(echo "$L" | awk '{print $9}')
  NOW=$(date +%s)
  SPS="nan"; ETA="nan"
  if [ -n "$STEP" ] && [ "$STEP" -gt "$PREV_STEP" ] 2>/dev/null; then
    SPS=$(awk -v a=$((NOW-PREV_EPOCH)) -v n=$((STEP-PREV_STEP)) 'BEGIN{printf "%.2f", a/n}')
    ETA=$(awk -v t="$T" -v dt="$DT" -v s="$SPS" -v te="$TSTOP" 'BEGIN{if(dt>0&&s>0)printf "%.1f",(te-t)/dt*s/3600; else print "nan"}')
  fi
  PREV_STEP=${STEP:-$PREV_STEP}; PREV_EPOCH=$NOW
  read -r GT GW GV GU <<< "$(nvidia-smi --query-gpu=temperature.gpu,power.draw,memory.used,utilization.gpu --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ',')"
  FREE=$(awk '/MemAvailable/{print int($2/1024)}' /proc/meminfo)
  echo "$(date -u +%FT%TZ),$STEP,$(awk -v t="$T" 'BEGIN{printf "%.4f", t*1e3}'),$DT,$SPS,$ETA,$GT,$GW,$GV,$GU,$FREE" >> "$OUT"
done
echo "$(date -u +%FT%TZ),SOLVER_EXITED" >> "$OUT"
