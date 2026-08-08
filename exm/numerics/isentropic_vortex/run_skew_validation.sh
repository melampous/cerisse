#!/usr/bin/env bash
set -euo pipefail

case_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
solver_root=$(cd "$case_root/../../.." && pwd)
cd "$case_root"

results_root=${CERISSE_IV_SKEW_RESULTS_ROOT:-results/y9000x_skew_validation_20260728}
build_jobs=${CERISSE_BUILD_JOBS:-4}
mpi_ranks=${CERISSE_MPI_RANKS:-4}

if [[ ! "$build_jobs" =~ ^[1-9][0-9]*$ ]]; then
  echo "CERISSE_BUILD_JOBS must be a positive integer" >&2
  exit 2
fi
if [[ ! "$mpi_ranks" =~ ^[1-9][0-9]*$ ]]; then
  echo "CERISSE_MPI_RANKS must be a positive integer" >&2
  exit 2
fi
if [[ -e "$results_root" ]]; then
  echo "$results_root already exists; select a new result root" >&2
  exit 1
fi

mkdir -p "$results_root/source_snapshot/solver" "$results_root/build"
cp GNUmakefile prob.h analyze.py analyze_skew_validation.py \
  inputs.spatial_convergence inputs.uniform run_skew_validation.sh \
  "$results_root/source_snapshot/"
cp "$solver_root/src/rhs/Skew.h" "$results_root/source_snapshot/solver/"
cp "$solver_root/src/rhs/RHS.h" "$results_root/source_snapshot/solver/"
cp "$solver_root/src/tim/compute_rhs.cpp" "$results_root/source_snapshot/solver/"
git status --short >"$results_root/source_snapshot/git_status.txt" || true
git diff -- src/rhs/Skew.h exm/numerics/isentropic_vortex \
  >"$results_root/source_snapshot/worktree.patch" || true
(
  cd "$results_root/source_snapshot"
  find . -type f -print0 | sort -z | xargs -0 sha256sum
) >"$results_root/source_snapshot/SHA256SUMS"

git_revision=$(git rev-parse HEAD 2>/dev/null || printf 'unavailable')
started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
{
  echo "case=isentropic_vortex"
  echo "study=skew_central_and_skew_jst_validation"
  echo "host=$(hostname)"
  echo "git_revision=${git_revision}"
  echo "worktree_clean=false"
  echo "source_snapshot=${results_root}/source_snapshot"
  echo "skew_source_sha256=$(sha256sum "$solver_root/src/rhs/Skew.h" | awk '{print $1}')"
  echo "mpi_ranks=${mpi_ranks}"
  echo "build_jobs=${build_jobs}"
  echo "rk_order=3"
  echo "rk_stages=3"
  echo "started_utc=${started_utc}"
  echo "low_temporal_error_angle_deg=45.0"
  echo "low_temporal_error_final_time=5.75953652404290800e-05"
  echo "low_temporal_error_grid_sequence=40 80 160 320"
  echo "low_temporal_error_step_sequence=64 256 1024 4096"
  echo "low_temporal_error_dt_scaling=dt_proportional_to_dx_squared"
  echo "one_transit_grid_sequence=40 80 160"
  echo "one_transit_final_time=5.759536524042908e-4"
  echo "five_transit_grid=80"
  echo "five_transit_final_time=2.879768262021454e-3"
  echo "vorticity_diagnostic=periodic_fourth_order_centred"
  echo "core_circulation_radius=sqrt(2)_times_vortex_radius"
} >"$results_root/manifest.txt"

build_scheme() {
  local scheme=$1
  local build_log="$results_root/build/${scheme}.log"
  local build_status="$results_root/build/${scheme}.status"
  local build_start
  local build_finish
  local status

  echo "building ${scheme}"
  build_start=$(date +%s)
  set +e
  make -j"$build_jobs" IV_EULER_SCHEME="$scheme" \
    USE_PELEPHYSICS=FALSE USE_EB=FALSE USE_GPIBM=FALSE \
    >"$build_log" 2>&1
  status=$?
  set -e
  build_finish=$(date +%s)
  {
    echo "exit_status=${status}"
    echo "wall_seconds=$((build_finish - build_start))"
    echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  } >"$build_status"
  if [[ "$status" -ne 0 ]]; then
    echo "build failed for ${scheme}; inspect ${build_log}" >&2
    exit "$status"
  fi

  local executable="./main2d.gnu.MPI.iv_${scheme}.ex"
  if [[ ! -x "$executable" ]]; then
    echo "expected executable ${executable} was not produced" >&2
    exit 1
  fi
  {
    echo "scheme=${scheme}"
    echo "executable=${executable}"
    echo "sha256=$(sha256sum "$executable" | awk '{print $1}')"
    echo "size_bytes=$(stat -c %s "$executable")"
  } >"$results_root/build/${scheme}.executable"
}

write_suite_manifest() {
  local suite_root=$1
  local scheme=$2
  local study=$3
  local executable=$4
  local dissipation=$5
  local c2=$6
  local c4=$7

  {
    echo "case=isentropic_vortex"
    echo "study=${study}"
    echo "scheme=${scheme}"
    echo "skew_order=4"
    echo "jst_dissipation=${dissipation}"
    echo "C2skew=${c2}"
    echo "C4skew=${c4}"
    echo "executable=${executable}"
    echo "executable_sha256=$(sha256sum "$executable" | awk '{print $1}')"
    echo "skew_source_sha256=$(sha256sum "$solver_root/src/rhs/Skew.h" | awk '{print $1}')"
    echo "git_revision=${git_revision}"
    echo "host=$(hostname)"
    echo "mpi_ranks=${mpi_ranks}"
    echo "rk_order=3"
    echo "rk_stages=3"
    echo "created_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  } >"$suite_root/manifest.txt"
}

run_case() {
  local output_dir=$1
  shift
  local executable=$1
  shift
  local input_file=$1
  shift
  local run_args=("$@")
  local start_seconds
  local finish_seconds
  local status

  mkdir -p "$output_dir"
  {
    echo '#!/usr/bin/env bash'
    echo 'set -euo pipefail'
    printf 'cd %q\n' "$case_root"
    printf '%q ' mpirun -np "$mpi_ranks" "$executable" "$input_file" \
      "${run_args[@]}"
    printf '\n'
  } >"$output_dir/run_command.sh"
  chmod +x "$output_dir/run_command.sh"

  start_seconds=$(date +%s)
  set +e
  mpirun -np "$mpi_ranks" "$executable" "$input_file" "${run_args[@]}" \
    >"$output_dir/run.log" 2>&1
  status=$?
  set -e
  finish_seconds=$(date +%s)
  {
    echo "exit_status=${status}"
    echo "wall_seconds=$((finish_seconds - start_seconds))"
    echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  } >"$output_dir/run_status.txt"
  if [[ "$status" -ne 0 ]]; then
    echo "run failed in ${output_dir}; inspect run.log" >&2
    exit "$status"
  fi
}

run_low_temporal_error() {
  local scheme=$1
  local dissipation=$2
  local c2=$3
  local c4=$4
  local executable="./main2d.gnu.MPI.iv_${scheme}.ex"
  local suite_root="$results_root/low_temporal_error/${scheme}"
  local cells
  local index
  local output_dir
  local result_dirs=()
  local cell_counts=(40 80 160 320)
  local step_counts=(64 256 1024 4096)
  local time_steps=(
    8.99927581881704387e-07
    2.24981895470426097e-07
    5.62454738676065242e-08
    1.40613684669016310e-08
  )
  local final_time=5.75953652404290800e-05

  mkdir -p "$suite_root"
  write_suite_manifest "$suite_root" "$scheme" \
    low_temporal_error_diagonal "$executable" "$dissipation" "$c2" "$c4"
  {
    echo "flow_angle_deg=45.0"
    echo "final_time=${final_time}"
    echo "grid_sequence=40 80 160 320"
    echo "step_sequence=64 256 1024 4096"
    echo "time_step_scaling=dt_proportional_to_dx_squared"
  } >>"$suite_root/manifest.txt"

  for index in "${!cell_counts[@]}"; do
    cells=${cell_counts[$index]}
    output_dir="$suite_root/N${cells}"
    echo "running ${scheme} low-temporal-error N=${cells}"
    run_case "$output_dir" "$executable" inputs.spatial_convergence \
      "max_step=${step_counts[$index]}" \
      "time_step=${time_steps[$index]}" \
      "amr.n_cell=${cells} ${cells}" \
      "amr.plot_file=${output_dir}/plt" \
      "prob.flow_angle_deg=45.0"
    result_dirs+=("$output_dir")
  done

  MPLCONFIGDIR=/tmp/cerisse-skew-mpl python3 analyze_skew_validation.py \
    "${result_dirs[@]}" --flow-angle-deg 45.0 \
    --csv "$suite_root/errors.csv" \
    >"$suite_root/analysis.log" 2>&1
  cat "$suite_root/analysis.log"
}

run_one_transit() {
  local scheme=skew-jst-o4
  local executable="./main2d.gnu.MPI.iv_${scheme}.ex"
  local suite_root="$results_root/one_transit/${scheme}"
  local result_dirs=()
  local cells
  local output_dir

  mkdir -p "$suite_root"
  write_suite_manifest "$suite_root" "$scheme" one_transit "$executable" \
    true 1.5 0.016
  {
    echo "flow_angle_deg=0.0"
    echo "final_time=5.759536524042908e-4"
    echo "grid_sequence=40 80 160"
    echo "cfl=0.30"
  } >>"$suite_root/manifest.txt"

  for cells in 40 80 160; do
    output_dir="$suite_root/N${cells}"
    echo "running ${scheme} one-transit N=${cells}"
    run_case "$output_dir" "$executable" inputs.uniform \
      "amr.n_cell=${cells} ${cells}" \
      "amr.plot_file=${output_dir}/plt"
    result_dirs+=("$output_dir")
  done

  MPLCONFIGDIR=/tmp/cerisse-skew-mpl python3 analyze_skew_validation.py \
    "${result_dirs[@]}" --csv "$suite_root/errors.csv" \
    >"$suite_root/analysis.log" 2>&1
  cat "$suite_root/analysis.log"
}

run_five_transit() {
  local scheme=skew-jst-o4
  local executable="./main2d.gnu.MPI.iv_${scheme}.ex"
  local suite_root="$results_root/five_transit/${scheme}"
  local output_dir="$suite_root/N80"

  mkdir -p "$suite_root"
  write_suite_manifest "$suite_root" "$scheme" five_transit "$executable" \
    true 1.5 0.016
  {
    echo "flow_angle_deg=0.0"
    echo "final_time=2.879768262021454e-3"
    echo "grid=80"
    echo "cfl=0.30"
  } >>"$suite_root/manifest.txt"

  echo "running ${scheme} five-transit N=80"
  run_case "$output_dir" "$executable" inputs.uniform \
    "stop_time=2.879768262021454e-3" \
    "amr.n_cell=80 80" \
    "amr.plot_file=${output_dir}/plt"
  MPLCONFIGDIR=/tmp/cerisse-skew-mpl python3 analyze_skew_validation.py \
    "$output_dir" --csv "$suite_root/errors.csv" \
    >"$suite_root/analysis.log" 2>&1
  cat "$suite_root/analysis.log"
}

build_scheme skew-central-o4
build_scheme skew-jst-o4
run_low_temporal_error skew-central-o4 false 0.0 0.0
run_low_temporal_error skew-jst-o4 true 1.5 0.016
run_one_transit
run_five_transit

echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)" >>"$results_root/manifest.txt"
echo "completed Skew validation: $results_root"
