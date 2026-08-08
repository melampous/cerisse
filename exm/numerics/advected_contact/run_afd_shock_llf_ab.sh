#!/usr/bin/env bash
set -euo pipefail

case_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$case_root/../../.." && pwd)
cd "$case_root"

grids=(64 128 256 512)
modes=(0 1)
build_jobs=${CERISSE_BUILD_JOBS:-4}
output_root=${CONTACT_AB_RESULTS_ROOT:-"$case_root/results/afd_shock_llf_ab_20260729_r1"}
executable="$case_root/main1d.gnu.MPI.contact_afd-hllc-wenoz5.ex"

if [[ "$output_root" != /* ]]; then
  output_root="$case_root/$output_root"
fi
if [[ -e "$output_root" ]]; then
  echo "$output_root already exists; select a new CONTACT_AB_RESULTS_ROOT" >&2
  exit 1
fi
mkdir -p "$output_root"

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

{
  printf 'case=periodic_advected_square_contact\n'
  printf 'created_utc='
  date -u +%Y-%m-%dT%H:%M:%SZ
  printf 'hostname=%s\n' "$(hostname)"
  printf 'git_commit=%s\n' "$(git -C "$repo_root" rev-parse HEAD)"
  printf 'committed_git_tree=%s\n' "$(git -C "$repo_root" rev-parse 'HEAD^{tree}')"
  printf 'grids=%s\n' "${grids[*]}"
  printf 'afd_shock_llf_modes=%s\n' "${modes[*]}"
  printf 'cfl=0.20\n'
  printf 'stop_time=1.6903085094570331\n'
  printf 'time_integrator=SSPRK33\n'
  printf 'scheme=AFD-HLLC-WENO-Z5\n'
  printf 'afd_correction=1\n'
  printf 'afd_smoothness_threshold=0.08\n'
  printf 'afd_shock_pressure_jump=0.01\n'
  printf 'afd_shock_compression=0.001\n'
  printf 'afd_teno_cutoff=1e-4\n'
  printf 'single_executable_for_both_modes=1\n'
  printf 'source_sha256_begin\n'
  sha256sum "$repo_root/src/rhs/Afd.h" \
            "$repo_root/src/rhs/AfdIBM.h" \
            "$repo_root/src/rhs/Riemann.h"
  printf 'source_sha256_end\n'
  printf 'git_status_begin\n'
  git -C "$repo_root" status --short
  printf 'git_status_end\n'
} >"$output_root/manifest.txt"

mkdir -p "$output_root/source_snapshot/case" \
         "$output_root/source_snapshot/src/rhs" \
         "$output_root/source_snapshot/src/tim"
cp GNUmakefile prob.h inputs run_afd_shock_llf_ab.sh analyze_ab.py README.md \
  "$output_root/source_snapshot/case/"
cp "$repo_root/src/rhs/Afd.h" \
   "$repo_root/src/rhs/AfdIBM.h" \
   "$repo_root/src/rhs/Riemann.h" \
   "$repo_root/src/rhs/RHS.h" \
   "$repo_root/src/rhs/Weno.h" \
   "$output_root/source_snapshot/src/rhs/"
cp "$repo_root/src/tim/advance.cpp" "$output_root/source_snapshot/src/tim/"
git -C "$repo_root" diff --no-ext-diff -- \
  src/rhs/Afd.h src/rhs/AfdIBM.h src/rhs/Riemann.h src/rhs/RHS.h \
  src/rhs/Weno.h src/tim/advance.cpp \
  >"$output_root/source_snapshot/working_tree.diff"
(
  cd "$output_root/source_snapshot"
  find . -type f -print0 | sort -z | xargs -0 sha256sum
) >"$output_root/source_snapshot.sha256"

build_command=(
  make "-j$build_jobs"
  USE_PELEPHYSICS=FALSE
  USE_EB=FALSE
  USE_GPIBM=FALSE
)
write_command "$output_root/build_command.txt" "${build_command[@]}"
started=$(date +%s)
set +e
"${build_command[@]}" >"$output_root/build.log" 2>&1
status=$?
set -e
finished=$(date +%s)
write_status "$output_root/build_status.txt" "$status" "$((finished - started))"
if [[ "$status" -ne 0 ]]; then
  echo "build failed; inspect $output_root/build.log" >&2
  exit "$status"
fi
sha256sum "$executable" >"$output_root/executable.sha256"

for mode in "${modes[@]}"; do
  for ncell in "${grids[@]}"; do
    destination="$output_root/runs/shock_llf_${mode}/N${ncell}"
    mkdir -p "$destination"
    command=(
      "$executable"
      "$case_root/inputs"
      "amr.n_cell=$ncell"
      "amr.plot_file=$destination/plt"
      "cns.afd_shock_llf=$mode"
      "cns.afd_shock_pressure_jump=0.01"
      "cns.afd_shock_compression=0.001"
    )
    write_command "$destination/command.txt" "${command[@]}"
    echo "running contact: afd_shock_llf=$mode N=$ncell"
    started=$(date +%s)
    set +e
    "${command[@]}" >"$destination/run.log" 2>&1
    status=$?
    set -e
    finished=$(date +%s)
    write_status "$destination/status.txt" "$status" "$((finished - started))"
    if [[ "$status" -ne 0 ]]; then
      echo "run failed: $destination" >&2
      exit "$status"
    fi
  done
done
touch "$output_root/all_solver_runs_complete"

analysis_command=(python3 "$case_root/analyze_ab.py" "$output_root")
write_command "$output_root/analysis_command.txt" "${analysis_command[@]}"
started=$(date +%s)
set +e
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
