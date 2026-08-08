#!/usr/bin/env bash
set -euo pipefail

case_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$case_root"

build_jobs=${CERISSE_BUILD_JOBS:-4}
make -j"$build_jobs" IV_EULER_SCHEME=llf-wenoz5 \
  USE_PELEPHYSICS=FALSE USE_EB=FALSE USE_GPIBM=FALSE
executable=./main2d.gnu.MPI.iv_llf-wenoz5.ex

if [[ ! -x "$executable" ]]; then
  echo "expected executable $executable was not produced" >&2
  exit 1
fi

for cells in 40 80 160; do
  output_dir="results/uniform_N${cells}"
  if [[ -e "$output_dir" ]]; then
    echo "$output_dir already exists; move it aside before rerunning" >&2
    exit 1
  fi
  mkdir -p "$output_dir"
  "$executable" inputs.uniform \
    amr.n_cell="$cells $cells" \
    amr.plot_file="${output_dir}/plt" \
    >"${output_dir}/run.log" 2>&1
done

python3 analyze.py \
  results/uniform_N40 results/uniform_N80 results/uniform_N160 \
  --csv results/uniform_errors.csv
