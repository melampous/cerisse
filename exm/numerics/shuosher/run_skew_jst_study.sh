#!/usr/bin/env bash
set -euo pipefail

# Auditable order-4 Skew--JST Shu--Osher study.  This script writes to a new
# archive and does not modify the existing four-scheme result directories.

case_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
solver_root=$(cd "$case_root/../../.." && pwd)
cd "$case_root"

cells_list=(200 400 800)
reference_cells=${SHU_SKEW_REFERENCE_CELLS:-6400}
build_jobs=${CERISSE_BUILD_JOBS:-4}
python_command=${PYTHON:-python3}
skew_c2=${SHU_SKEW_C2:-1.5}
skew_c4=${SHU_SKEW_C4:-0.016}
skew_c2_tag=${skew_c2//./_}
skew_c4_tag=${skew_c4//./_}
output_root=${SHU_SKEW_ROOT:-"$case_root/results/y9000x_skew_jst_c2_${skew_c2_tag}_c4_${skew_c4_tag}_20260728_r1"}
if [[ "$output_root" != /* ]]; then
  output_root="$case_root/$output_root"
fi

run_prefix=()
if [[ -n "${CERISSE_RUN_PREFIX:-}" ]]; then
  read -r -a run_prefix <<<"${CERISSE_RUN_PREFIX}"
fi

if [[ -e "$output_root" ]]; then
  echo "$output_root already exists; select a new SHU_SKEW_ROOT" >&2
  exit 1
fi
mkdir -p "$output_root/source_snapshot"

write_command() {
  local destination=$1
  shift
  {
    printf 'recorded_utc='
    date -u +%Y-%m-%dT%H:%M:%SZ
    printf 'working_directory=%q\n' "$case_root"
    printf 'command='
    printf '%q ' "$@"
    printf '\n'
  } >"$destination"
}

run_logged_command() {
  local log=$1
  local status_file=$2
  shift 2
  local start_seconds
  local status
  start_seconds=$(date +%s)
  if "$@" >"$log" 2>&1; then
    status=0
  else
    status=$?
  fi
  {
    printf 'exit_status=%d\n' "$status"
    printf 'wall_seconds=%d\n' "$(( $(date +%s) - start_seconds ))"
    printf 'finished_utc='
    date -u +%Y-%m-%dT%H:%M:%SZ
  } >"$status_file"
  return "$status"
}

run_case() {
  local destination=$1
  local executable=$2
  shift 2
  local command=("${run_prefix[@]}" "$executable" "$case_root/inputs.thesis" "$@")
  mkdir -p "$destination"
  write_command "$destination/run_command.txt" "${command[@]}"
  if ! run_logged_command \
      "$destination/run.log" "$destination/run_status.txt" "${command[@]}"; then
    echo "run failed; see $destination/run.log" >&2
    return 1
  fi
}

{
  printf 'created_utc='
  date -u +%Y-%m-%dT%H:%M:%SZ
  printf 'hostname=%s\n' "$(hostname)"
  printf 'working_directory=%s\n' "$case_root"
  printf 'solver_root=%s\n' "$solver_root"
  printf 'output_root=%s\n' "$output_root"
  printf 'scheme=skew-jst\n'
  printf 'skew_order=4\n'
  printf 'skew_C2=%s\n' "$skew_c2"
  printf 'skew_C4=%s\n' "$skew_c4"
  printf 'dissipation_enabled=true\n'
  printf 'cells=%s\n' "${cells_list[*]}"
  printf 'numerical_reference_cells=%s\n' "$reference_cells"
  printf 'numerical_reference_is_exact=false\n'
  printf 'formal_order_inferred=false\n'
  printf 'stop_time=1.8\n'
  printf 'cfl=0.30\n'
  printf 'run_prefix=%s\n' "${CERISSE_RUN_PREFIX:-<direct>}"
  printf 'build_jobs=%s\n' "$build_jobs"
  printf 'git_commit=%s\n' "$(git -C "$solver_root" rev-parse HEAD 2>/dev/null || printf unknown)"
  printf 'uname='
  uname -a
  printf 'compiler='
  "${CXX:-g++}" --version 2>/dev/null | head -n 1 || printf 'unknown\n'
} >"$output_root/run_manifest.txt"

git -C "$solver_root" status --short >"$output_root/git_status_short.txt" 2>&1 || true
git -C "$solver_root" diff --binary >"$output_root/git_worktree.patch" 2>&1 || true
git -C "$solver_root" diff --cached --binary \
  >"$output_root/git_index.patch" 2>&1 || true

cp "$case_root/GNUmakefile" \
  "$case_root/Make.package" \
  "$case_root/prob.h" \
  "$case_root/inputs.thesis" \
  "$case_root/run_skew_jst_study.sh" \
  "$case_root/analyze.py" \
  "$case_root/analyze_skew_jst.py" \
  "$output_root/source_snapshot/"
mkdir -p "$output_root/source_snapshot/src_rhs"
cp "$solver_root/src/rhs/Skew.h" "$output_root/source_snapshot/src_rhs/Skew.h"
sha256sum \
  "$output_root/source_snapshot/GNUmakefile" \
  "$output_root/source_snapshot/prob.h" \
  "$output_root/source_snapshot/inputs.thesis" \
  "$output_root/source_snapshot/src_rhs/Skew.h" \
  >"$output_root/source_snapshot/SHA256SUMS"

build_command=(
  make
  -j"$build_jobs"
  SHU_EULER_SCHEME=skew
  SHU_SKEW_C2="$skew_c2"
  SHU_SKEW_C4="$skew_c4"
  USE_PELEPHYSICS=FALSE
  USE_EB=FALSE
  USE_GPIBM=FALSE
)
write_command "$output_root/build_command.txt" "${build_command[@]}"
if ! run_logged_command \
    "$output_root/build.log" "$output_root/build_status.txt" \
    "${build_command[@]}"; then
  echo "Skew--JST build failed; see $output_root/build.log" >&2
  exit 1
fi

executable="$case_root/main1d.gnu.MPI.shu_skew_c2_${skew_c2_tag}_c4_${skew_c4_tag}.ex"
if [[ ! -x "$executable" ]]; then
  echo "expected executable $executable was not produced" >&2
  exit 1
fi
sha256sum "$executable" >"$output_root/executable.sha256"
cp "$executable" "$output_root/"

run_case "$output_root/smoke" "$executable" \
  amr.n_cell=64 \
  max_step=20 \
  stop_time=0.02 \
  amr.plot_files_output=0 \
  amr.checkpoint_files_output=0 \
  cns.nstep_screen_output=20
touch "$output_root/smoke_passed"

for cells in "${cells_list[@]}"; do
  destination="$output_root/N${cells}"
  run_case "$destination" "$executable" \
    amr.n_cell="$cells" \
    amr.max_grid_size=8192 \
    amr.plot_files_output=1 \
    amr.plot_file="$destination/plt" \
    amr.checkpoint_files_output=0 \
    cns.nstep_screen_output=200
done

reference_destination="$output_root/N${reference_cells}_reference"
run_case "$reference_destination" "$executable" \
  amr.n_cell="$reference_cells" \
  amr.max_grid_size=8192 \
  amr.plot_files_output=1 \
  amr.plot_file="$reference_destination/plt" \
  amr.checkpoint_files_output=0 \
  cns.nstep_screen_output=2000
touch "$output_root/all_solver_runs_complete"

analysis_command=(
  "$python_command"
  "$case_root/analyze_skew_jst.py"
  "$output_root"
  --reference "$reference_destination"
  --reference-label "Skew-JST N=${reference_cells} numerical reference"
  --window-lo -1.5
  --window-hi 1.5
  --csv "$output_root/skew_jst_metrics.csv"
  --figure "$output_root/skew_jst_density.png"
  --metadata "$output_root/analysis_metadata.json"
)
write_command "$output_root/analysis_command.txt" "${analysis_command[@]}"
if ! run_logged_command \
    "$output_root/analysis.log" "$output_root/analysis_status.txt" \
    "${analysis_command[@]}"; then
  echo "analysis failed; see $output_root/analysis.log" >&2
  exit 1
fi

touch "$output_root/study_complete"
echo "Completed Skew--JST Shu--Osher study: $output_root"
