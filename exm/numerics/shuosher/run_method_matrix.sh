#!/usr/bin/env bash
set -euo pipefail

# Four-scheme Shu--Osher comparison for an auditable single-node run.
#
# The existing run_study.sh remains the LLF--WENO-Z5 thesis baseline.  This
# script deliberately omits the expensive N=6400 reference and uses N=800 only
# as a same-method resolution-difference reference during optional analysis.

case_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$case_root"

schemes=(
  llf-wenoz5
  llf-teno5
  afd-hllc-wenoz5
  afd-hllc-teno5
)
cells_list=(200 400 800)

build_jobs=${CERISSE_BUILD_JOBS:-4}
python_command=${PYTHON:-python3}
output_root=${SHU_MATRIX_ROOT:-"$case_root/results/y9000x_method_matrix"}
if [[ "$output_root" != /* ]]; then
  output_root="$case_root/$output_root"
fi

# Set, for example, CERISSE_RUN_PREFIX="mpirun -np 1" when a launcher is
# required.  A direct single-process launch is preferable for these small 1-D
# cases and is therefore the default.
run_prefix=()
if [[ -n "${CERISSE_RUN_PREFIX:-}" ]]; then
  read -r -a run_prefix <<<"${CERISSE_RUN_PREFIX}"
fi

if [[ -e "$output_root" ]]; then
  echo "$output_root already exists; select a new SHU_MATRIX_ROOT" >&2
  exit 1
fi
mkdir -p "$output_root"

write_command() {
  local destination=$1
  shift
  {
    printf 'started_utc='
    date -u +%Y-%m-%dT%H:%M:%SZ
    printf 'working_directory=%q\n' "$case_root"
    printf 'command='
    printf '%q ' "$@"
    printf '\n'
  } >"$destination"
}

run_case() {
  local destination=$1
  local executable=$2
  shift 2
  local command=("${run_prefix[@]}" "$executable" "$case_root/inputs.thesis" "$@")

  mkdir -p "$destination"
  write_command "$destination/run_command.txt" "${command[@]}"
  local start_seconds
  start_seconds=$(date +%s)
  if "${command[@]}" >"$destination/run.log" 2>&1; then
    local end_seconds
    end_seconds=$(date +%s)
    {
      printf 'exit_status=0\n'
      printf 'wall_seconds=%d\n' "$((end_seconds - start_seconds))"
      printf 'finished_utc='
      date -u +%Y-%m-%dT%H:%M:%SZ
    } >"$destination/run_status.txt"
  else
    local status=$?
    local end_seconds
    end_seconds=$(date +%s)
    {
      printf 'exit_status=%d\n' "$status"
      printf 'wall_seconds=%d\n' "$((end_seconds - start_seconds))"
      printf 'failed_utc='
      date -u +%Y-%m-%dT%H:%M:%SZ
    } >"$destination/run_status.txt"
    echo "run failed; see $destination/run.log" >&2
    return "$status"
  fi
}

{
  printf 'created_utc='
  date -u +%Y-%m-%dT%H:%M:%SZ
  printf 'hostname=%s\n' "$(hostname)"
  printf 'working_directory=%s\n' "$case_root"
  printf 'output_root=%s\n' "$output_root"
  printf 'build_jobs=%s\n' "$build_jobs"
  printf 'run_prefix=%s\n' "${CERISSE_RUN_PREFIX:-<direct>}"
  printf 'schemes=%s\n' "${schemes[*]}"
  printf 'cells=%s\n' "${cells_list[*]}"
  printf 'git_commit=%s\n' "$(git rev-parse HEAD 2>/dev/null || printf unknown)"
  printf 'uname='
  uname -a
  printf 'compiler='
  "${CXX:-g++}" --version 2>/dev/null | head -n 1 || printf 'unknown\n'
  printf '\ngit_status_begin\n'
  git status --short 2>/dev/null || true
  printf 'git_status_end\n'
} >"$output_root/run_manifest.txt"

cp "$case_root/inputs.thesis" "$output_root/inputs.thesis"
mkdir -p "$output_root/source_snapshot"
cp "$case_root/GNUmakefile" \
  "$case_root/prob.h" \
  "$case_root/inputs.thesis" \
  "$case_root/run_method_matrix.sh" \
  "$case_root/analyze.py" \
  "$case_root/analyze_method_matrix.py" \
  "$output_root/source_snapshot/"

declare -A executables

# Phase 1: build every method and require every short smoke test to pass before
# any production-resolution calculation is started.
for scheme in "${schemes[@]}"; do
  scheme_root="$output_root/$scheme"
  mkdir -p "$scheme_root"
  echo "Building $scheme"
  write_command "$scheme_root/build_command.txt" \
    make -j"$build_jobs" SHU_EULER_SCHEME="$scheme" \
    USE_PELEPHYSICS=FALSE USE_EB=FALSE USE_GPIBM=FALSE
  make -j"$build_jobs" SHU_EULER_SCHEME="$scheme" \
    USE_PELEPHYSICS=FALSE USE_EB=FALSE USE_GPIBM=FALSE \
    >"$scheme_root/build.log" 2>&1

  executable="$case_root/main1d.gnu.MPI.shu_${scheme}.ex"
  if [[ ! -x "$executable" ]]; then
    echo "expected executable $executable was not produced" >&2
    exit 1
  fi
  executables["$scheme"]=$executable
  sha256sum "$executable" >"$scheme_root/executable.sha256"

  method_args=()
  if [[ "$scheme" == afd-hllc-* ]]; then
    method_args=(
      cns.afd_correction=1
      cns.afd_shock_llf=1
      cns.afd_smoothness_threshold=0.08
      cns.afd_teno_cutoff=1.0e-4
    )
  fi

  echo "Smoke testing $scheme"
  run_case "$scheme_root/smoke" "$executable" \
    amr.n_cell=64 \
    max_step=20 \
    stop_time=0.02 \
    amr.plot_files_output=0 \
    amr.checkpoint_files_output=0 \
    cns.nstep_screen_output=20 \
    "${method_args[@]}"
done
touch "$output_root/all_smoke_tests_passed"

# Phase 2: full t=1.8 calculations on the three requested uniform grids.
for scheme in "${schemes[@]}"; do
  executable=${executables["$scheme"]}
  method_args=()
  if [[ "$scheme" == afd-hllc-* ]]; then
    method_args=(
      cns.afd_correction=1
      cns.afd_shock_llf=1
      cns.afd_smoothness_threshold=0.08
      cns.afd_teno_cutoff=1.0e-4
    )
  fi

  for cells in "${cells_list[@]}"; do
    destination="$output_root/$scheme/N${cells}"
    echo "Running $scheme at N=$cells"
    run_case "$destination" "$executable" \
      amr.n_cell="$cells" \
      amr.max_grid_size=8192 \
      amr.plot_files_output=1 \
      amr.plot_file="$destination/plt" \
      amr.checkpoint_files_output=0 \
      cns.nstep_screen_output=200 \
      "${method_args[@]}"
  done
done
touch "$output_root/all_solver_runs_complete"

# Profile analysis needs numpy, yt and matplotlib.  Solver completion is not
# made dependent on optional Python packages because they may not be installed
# on a compute host.  The same command can be run after copying the archive to
# a post-processing machine.
analysis_command=(
  "$python_command"
  "$case_root/analyze_method_matrix.py"
  "$output_root"
  --csv "$output_root/method_matrix_metrics.csv"
  --figure "$output_root/method_matrix_density.png"
)
write_command "$output_root/analysis_command.txt" "${analysis_command[@]}"
if "$python_command" -c 'import matplotlib, numpy, yt' >/dev/null 2>&1; then
  "${analysis_command[@]}" >"$output_root/analysis.log" 2>&1
  touch "$output_root/analysis_complete"
else
  {
    echo "Solver matrix completed, but optional profile analysis was skipped."
    echo "numpy, yt and matplotlib are required."
    printf 'Run: '
    printf '%q ' "${analysis_command[@]}"
    printf '\n'
  } >"$output_root/analysis_skipped.txt"
fi

touch "$output_root/matrix_complete"
echo "Completed Shu--Osher method matrix: $output_root"
