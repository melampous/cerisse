#!/usr/bin/env bash

# Read-only health and performance monitor for the active 3D L5 calculation.

set -uo pipefail

D3=/shared/cerisse_reflux/exm/underexpanded_jet/3d/m142_3d_gpu
RUNROOT="$D3/l5_prod_cb_20260715"
OUT="$RUNROOT/speed_monitor.log"
INTERVAL=${INTERVAL:-300}

exec 9>"$RUNROOT/speed_monitor.lock"
flock -n 9 || exit 0

if [[ ! -s "$OUT" ]]; then
  printf '%s\n' \
    '# utc phase step physical_ms dt_us median20_s ms_per_hour gpu_util_pct gpu_mem_mib contexts status' \
    > "$OUT"
fi

previous_step=-1
unchanged_samples=0

while :; do
  now=$(date -u +%FT%TZ)
  process=$(pgrep -af '[m]ain3d\.gnu\.MPI\.CUDA\.ex' | head -1 || true)

  if [[ -z "$process" ]]; then
    if pgrep -f '[r]un_3d_after_m180_l5\.sh' >/dev/null; then
      printf '%s phase_transition - - - - - - - 0 WAIT\n' "$now" >> "$OUT"
      sleep 15
      continue
    fi
    printf '%s complete - - - - - - - 0 DONE\n' "$now" >> "$OUT"
    exit 0
  fi

  phase=$(sed -n 's/.*probes_\([[:alnum:]_-]*\)\.log.*/\1/p' <<<"$process")
  [[ -n "$phase" ]] || phase=unknown
  log_file="$RUNROOT/${phase}.log"

  metrics=$(python3 - "$log_file" <<'PY'
import re
import statistics
import sys

try:
    text = open(sys.argv[1], errors="ignore").read()
except OSError:
    print("- - - - -")
    raise SystemExit

rows = [
    (int(step), float(time), float(dt))
    for step, time, dt in re.findall(
        r"STEP = (\d+) TIME = ([0-9.eE+-]+) DT = ([0-9.eE+-]+)", text
    )
]
timings = [
    float(value)
    for value in re.findall(r"Coarse TimeStep time: ([0-9.eE+-]+)", text)
]
if not rows or not timings:
    print("- - - - -")
    raise SystemExit

step, physical_time, dt = rows[-1]
median = statistics.median(timings[-min(20, len(timings)):])
rate = dt / median * 3600.0 * 1.0e3
print(f"{step} {physical_time*1.0e3:.9f} {dt*1.0e6:.6f} {median:.6f} {rate:.6f}")
PY
  )
  read -r step physical_ms dt_us median20 rate <<<"$metrics"

  gpu_util=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits |
    awk '{sum+=$1; n++} END {if (n) printf "%.1f", sum/n; else print "-"}')
  gpu_mem=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits |
    awk 'NR==1{max=$1} $1>max{max=$1} END {if (NR) print max; else print "-"}')
  contexts=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader,nounits |
    sed '/^[[:space:]]*$/d' | sort -u | wc -l)

  status=OK
  if [[ "$step" != "-" ]]; then
    if [[ "$step" == "$previous_step" ]]; then
      unchanged_samples=$((unchanged_samples + 1))
    else
      unchanged_samples=0
    fi
    previous_step=$step

    if awk -v value="$median20" 'BEGIN {exit !(value > 6.5)}'; then
      status=WARN_SLOW
    elif (( unchanged_samples >= 2 )); then
      status=WARN_NO_PROGRESS
    elif (( contexts != 8 )); then
      status=WARN_GPU_CONTEXTS
    elif [[ "$gpu_mem" != "-" ]] && (( gpu_mem > 65000 )); then
      status=WARN_MEMORY
    fi
  else
    status=STARTING_OR_IO
  fi

  printf '%s %s %s %s %s %s %s %s %s %s %s\n' \
    "$now" "$phase" "$step" "$physical_ms" "$dt_us" "$median20" \
    "$rate" "$gpu_util" "$gpu_mem" "$contexts" "$status" >> "$OUT"

  sleep "$INTERVAL"
done
