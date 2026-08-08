#!/bin/bash
# Encode the 19 true SRP CFD states without optical-flow interpolation.
# Usage: encode_srp_tri_paper_movies.sh PAPER_SERIES_DIR

set -euo pipefail

ROOT=${1:?Usage: $0 PAPER_SERIES_DIR}
FPS=${FPS:-5}
mkdir -p "${ROOT}/animations"

for field in mach pressure_ratio schlieren; do
  input_dir=${ROOT}/frames/${field}
  output=${ROOT}/animations/${field}_threeview.mp4
  count=$(find "${input_dir}" -maxdepth 1 -type f -name 'plt*.png' | wc -l)
  if [[ ${count} -ne 19 ]]; then
    echo "Expected 19 ${field} frames, found ${count}" >&2
    exit 1
  fi
  ffmpeg -hide_banner -loglevel warning -y \
    -framerate "${FPS}" -pattern_type glob -i "${input_dir}/plt*.png" \
    -vf "tpad=start_mode=clone:start_duration=0.8:stop_mode=clone:stop_duration=0.8,fps=30,pad=ceil(iw/2)*2:ceil(ih/2)*2,format=yuv420p" \
    -c:v libx264 -preset slow -crf 18 -movflags +faststart \
    -metadata comment="19 true CFD states; fixed scales; no temporal interpolation" \
    "${output}"
  ffprobe -v error -show_entries \
    stream=codec_name,pix_fmt,width,height,r_frame_rate,nb_frames,duration \
    -show_entries format=duration,size -of json "${output}" \
    > "${ROOT}/animations/${field}_threeview.ffprobe.json"
done
