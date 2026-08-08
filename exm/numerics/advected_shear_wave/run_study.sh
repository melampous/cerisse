#!/usr/bin/env bash
set -euo pipefail

case_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$case_root/../../.." && pwd)
cd "$case_root"

schemes=(
  llf-wenoz5
  llf-teno5
  afd-hllc-wenoz5
  skew-jst-o4
)
ppw_values=(8 12 16 24 32 48)
repeat_ppw=24
build_jobs=${CERISSE_BUILD_JOBS:-4}
output_root=${SHEAR_WAVE_RESULTS_ROOT:-"$case_root/results/matrix_20260728"}

if [[ "$output_root" != /* ]]; then
  output_root="$case_root/$output_root"
fi
if [[ -e "$output_root" ]]; then
  echo "$output_root already exists; select a new SHEAR_WAVE_RESULTS_ROOT" >&2
  exit 1
fi
mkdir -p "$output_root"

run_prefix=()
if [[ -n "${CERISSE_RUN_PREFIX:-}" ]]; then
  read -r -a run_prefix <<<"${CERISSE_RUN_PREFIX}"
fi

write_command() {
  local destination=$1
  shift
  {
    printf 'working_directory=%q\n' "$case_root"
    printf 'command='
    printf '%q ' "$@"
    printf '\n'
  } >"$destination"
}

write_status() {
  local destination=$1
  local status=$2
  local elapsed=$3
  {
    printf 'exit_status=%d\n' "$status"
    printf 'wall_seconds=%d\n' "$elapsed"
    printf 'finished_utc='
    date -u +%Y-%m-%dT%H:%M:%SZ
  } >"$destination"
}

run_case() {
  local destination=$1
  local executable=$2
  local cells=$3
  local selected_cfl=$4
  shift 4
  local method_args=("$@")
  local command=(
    "${run_prefix[@]}"
    "$executable"
    "$case_root/inputs"
    "amr.n_cell=$cells"
    "cfl=$selected_cfl"
    "amr.plot_file=$destination/plt"
    "${method_args[@]}"
  )

  mkdir -p "$destination"
  write_command "$destination/command.txt" "${command[@]}"
  local started
  started=$(date +%s)
  set +e
  "${command[@]}" >"$destination/run.log" 2>&1
  local status=$?
  set -e
  local finished
  finished=$(date +%s)
  write_status "$destination/status.txt" "$status" "$((finished - started))"
  if [[ "$status" -ne 0 ]]; then
    echo "failed: $destination; inspect run.log" >&2
    return "$status"
  fi
}

{
  printf 'case=periodic_advected_shear_wave\n'
  printf 'created_utc='
  date -u +%Y-%m-%dT%H:%M:%SZ
  printf 'hostname=%s\n' "$(hostname)"
  printf 'repo_root=%s\n' "$repo_root"
  printf 'case_root=%s\n' "$case_root"
  printf 'output_root=%s\n' "$output_root"
  printf 'git_commit=%s\n' "$(git -C "$repo_root" rev-parse HEAD 2>/dev/null || printf unknown)"
  printf 'schemes=%s\n' "${schemes[*]}"
  printf 'ppw=%s\n' "${ppw_values[*]}"
  printf 'primary_cfl=0.05\n'
  printf 'repeat_cfl=0.025\n'
  printf 'repeat_ppw=%d\n' "$repeat_ppw"
  printf 'rho0=1\np0=1\ngamma=1.4\nmach_x=0.5\n'
  printf 'transverse_amplitude_over_c0=1e-5\nmode=1\n'
  printf 'stop_time=1.6903085094570331\n'
  printf 'time_integrator=SSPRK33\n'
  printf 'skew_order=4\nskew_C2=1.5\nskew_C4=0.016\n'
  printf 'afd_correction=1\n'
  printf 'afd_shock_llf=1\n'
  printf 'afd_smoothness_threshold=0.08\n'
  printf 'afd_teno_cutoff=1e-4\n'
  printf 'run_prefix=%s\n' "${CERISSE_RUN_PREFIX:-<direct>}"
  printf 'compiler='
  "${CXX:-g++}" --version 2>/dev/null | head -n 1 || printf 'unknown\n'
  printf 'uname='
  uname -a
  printf '\ngit_status_begin\n'
  git -C "$repo_root" status --short 2>/dev/null || true
  printf 'git_status_end\n'
} >"$output_root/manifest.txt"

mkdir -p "$output_root/source_snapshot/case"
cp GNUmakefile prob.h inputs run_study.sh analyze.py README.md \
  "$output_root/source_snapshot/case/"
mkdir -p "$output_root/source_snapshot/src/rhs"
cp "$repo_root/src/rhs/Skew.h" \
  "$repo_root/src/rhs/Weno.h" \
  "$repo_root/src/rhs/Afd.h" \
  "$repo_root/src/rhs/Riemann.h" \
  "$repo_root/src/rhs/RHS.h" \
  "$output_root/source_snapshot/src/rhs/"
mkdir -p "$output_root/source_snapshot/src/tim"
cp "$repo_root/src/tim/advance.cpp" \
  "$output_root/source_snapshot/src/tim/"
git -C "$repo_root" diff --no-ext-diff -- \
  src/rhs/Skew.h src/rhs/Weno.h src/rhs/Afd.h src/rhs/Riemann.h \
  src/rhs/RHS.h src/tim/advance.cpp \
  >"$output_root/source_snapshot/working_tree.diff"
(
  cd "$output_root/source_snapshot"
  find . -type f -print0 | sort -z | xargs -0 sha256sum
) >"$output_root/source_snapshot.sha256"

declare -A executables
for scheme in "${schemes[@]}"; do
  scheme_root="$output_root/build/$scheme"
  mkdir -p "$scheme_root"
  build_command=(
    make "-j$build_jobs"
    "SHEAR_EULER_SCHEME=$scheme"
    USE_PELEPHYSICS=FALSE
    USE_EB=FALSE
    USE_GPIBM=FALSE
  )
  write_command "$scheme_root/command.txt" "${build_command[@]}"
  started=$(date +%s)
  set +e
  "${build_command[@]}" >"$scheme_root/build.log" 2>&1
  status=$?
  set -e
  finished=$(date +%s)
  write_status "$scheme_root/status.txt" "$status" "$((finished - started))"
  if [[ "$status" -ne 0 ]]; then
    echo "build failed for $scheme" >&2
    exit "$status"
  fi

  executable="$case_root/main1d.gnu.MPI.shear_${scheme}.ex"
  if [[ ! -x "$executable" ]]; then
    echo "missing executable: $executable" >&2
    exit 1
  fi
  executables["$scheme"]=$executable
  sha256sum "$executable" >"$scheme_root/executable.sha256"

  method_args=()
  if [[ "$scheme" == "afd-hllc-wenoz5" ]]; then
    method_args=(
      cns.afd_correction=1
      cns.afd_shock_llf=1
      cns.afd_smoothness_threshold=0.08
      cns.afd_teno_cutoff=1.0e-4
    )
  fi

  smoke_dir="$output_root/runs/$scheme/smoke"
  command=(
    "${run_prefix[@]}"
    "$executable"
    "$case_root/inputs"
    amr.n_cell=8
    max_step=2
    stop_time=0.01
    amr.plot_files_output=0
    amr.checkpoint_files_output=0
    "${method_args[@]}"
  )
  mkdir -p "$smoke_dir"
  write_command "$smoke_dir/command.txt" "${command[@]}"
  started=$(date +%s)
  set +e
  "${command[@]}" >"$smoke_dir/run.log" 2>&1
  status=$?
  set -e
  finished=$(date +%s)
  write_status "$smoke_dir/status.txt" "$status" "$((finished - started))"
  if [[ "$status" -ne 0 ]]; then
    echo "smoke test failed for $scheme" >&2
    exit "$status"
  fi
done
touch "$output_root/all_builds_and_smoke_tests_passed"

for scheme in "${schemes[@]}"; do
  executable=${executables["$scheme"]}
  method_args=()
  if [[ "$scheme" == "afd-hllc-wenoz5" ]]; then
    method_args=(
      cns.afd_correction=1
      cns.afd_shock_llf=1
      cns.afd_smoothness_threshold=0.08
      cns.afd_teno_cutoff=1.0e-4
    )
  fi
  for ppw in "${ppw_values[@]}"; do
    destination="$output_root/runs/$scheme/cfl0p05/N$ppw"
    echo "running $scheme PPW=$ppw CFL=0.05"
    run_case "$destination" "$executable" "$ppw" 0.05 \
      "${method_args[@]}"
  done
  destination="$output_root/runs/$scheme/cfl0p025/N$repeat_ppw"
  echo "running $scheme PPW=$repeat_ppw CFL=0.025"
  run_case "$destination" "$executable" "$repeat_ppw" 0.025 \
    "${method_args[@]}"
done
touch "$output_root/all_solver_runs_complete"

analysis_command=(
  python3 "$case_root/analyze.py" "$output_root"
  --csv "$output_root/metrics.csv"
  --figure "$output_root/dissipation_phase.png"
)
write_command "$output_root/analysis_command.txt" "${analysis_command[@]}"
started=$(date +%s)
set +e
MPLCONFIGDIR=/tmp/cerisse-matplotlib \
  "${analysis_command[@]}" >"$output_root/analysis.log" 2>&1
status=$?
set -e
finished=$(date +%s)
write_status "$output_root/analysis_status.txt" "$status" "$((finished - started))"
if [[ "$status" -ne 0 ]]; then
  echo "analysis failed; inspect $output_root/analysis.log" >&2
  exit "$status"
fi

touch "$output_root/study_complete"
echo "completed: $output_root"
