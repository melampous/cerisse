#!/bin/bash
# 4090 babysitter: graceful-stop the solver before it can freeze the laptop.
# Triggers: GPU temp >= 87C, VRAM used >= 14000 MiB, or free system RAM < 1.5 GB.
D="$(cd "$(dirname "$0")" && pwd)"
LOG="$D/guard.log"
echo "[guard $(date -u +%FT%TZ)] armed (temp>=87C | vram>=15900MiB x3 consecutive (60s) | freeRAM<1.5G)" > "$LOG"
TRIPPED=0
VHITS=0
while pgrep -f "main2d.gnu.MPI.CUDA.e[x]" > /dev/null; do
  read -r T V <<< "$(nvidia-smi --query-gpu=temperature.gpu,memory.used --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ',')"
  FREE=$(awk '/MemAvailable/{print int($2/1024)}' /proc/meminfo)
  echo "[guard $(date -u +%FT%TZ)] T=${T}C vram=${V}MiB freeRAM=${FREE}MiB" >> "$LOG"
  # vram: debounce single-sample spikes (regrid transients / WSL query glitches)
  if [ "${V:-0}" -ge 15900 ] 2>/dev/null; then VHITS=$((VHITS+1)); else VHITS=0; fi
  if [ "${T:-0}" -ge 87 ] || [ "$VHITS" -ge 3 ] || [ "${FREE:-9999}" -lt 1500 ]; then
    if [ "$TRIPPED" -eq 0 ]; then
      echo "[guard $(date -u +%FT%TZ)] LIMIT HIT -> dump_and_stop" >> "$LOG"
      touch "$D/dump_and_stop"; TRIPPED=$(date +%s)
    elif [ $(( $(date +%s) - TRIPPED )) -gt 180 ]; then
      echo "[guard $(date -u +%FT%TZ)] still alive 180s after stop request -> SIGTERM" >> "$LOG"
      pkill -TERM -f "main2d.gnu.MPI.CUDA.e[x]"; sleep 20
      pkill -KILL -f "main2d.gnu.MPI.CUDA.e[x]"
      break
    fi
  fi
  sleep 20
done
rm -f "$D/dump_and_stop"
echo "[guard $(date -u +%FT%TZ)] solver gone, guard exits (stop file cleaned)" >> "$LOG"
