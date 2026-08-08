#!/usr/bin/env bash
set -euo pipefail

case_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$case_root"

build_jobs=${CERISSE_BUILD_JOBS:-4}
make -j"$build_jobs" SHU_EULER_SCHEME=llf-wenoz5 \
  USE_PELEPHYSICS=FALSE USE_EB=FALSE USE_GPIBM=FALSE
executable=./main1d.gnu.MPI.shu_llf-wenoz5.ex

if [[ ! -x "$executable" ]]; then
  echo "expected executable $executable was not produced" >&2
  exit 1
fi

for cells in 200 400 800 6400; do
  output_dir="results/N${cells}"
  if [[ -e "$output_dir" ]]; then
    echo "$output_dir already exists; move it aside before rerunning" >&2
    exit 1
  fi
  mkdir -p "$output_dir"
  "$executable" inputs.thesis \
    amr.n_cell="$cells" \
    amr.max_grid_size=8192 \
    amr.plot_file="${output_dir}/plt" \
    >"${output_dir}/run.log" 2>&1
done

python3 analyze.py results/N6400 results/N200 results/N400 results/N800 \
  --csv results/shuosher_errors.csv \
  --figure results/shuosher_density.png
