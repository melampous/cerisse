#!/usr/bin/env bash
set -euo pipefail

# Strict same-source and same-binary A/B test of cns.afd_shock_llf.
# Both variants use AFD--HLLC--WENO-Z5.  The only numerical-method argument
# that differs between a paired run is cns.afd_shock_llf.

case_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
solver_root=$(cd "$case_root/../../.." && pwd)
cd "$case_root"

cells_list=(200 400 800)
build_jobs=${CERISSE_BUILD_JOBS:-8}
mpi_ranks=${CERISSE_MPI_RANKS:-4}
build_tag=${SHU_LLF_AB_BUILD_TAG:-shock_llf_ab_20260729_r1}
output_root=${SHU_LLF_AB_ROOT:-"$case_root/results/shock_llf_ab_20260729_r1"}
reference_root=${SHU_REFERENCE_ROOT:-"$case_root/results/y9000x_skew_jst_c2_1_5_c4_0_016_20260728_r3/N6400_reference"}

if [[ "$output_root" != /* ]]; then
  output_root="$case_root/$output_root"
fi
if [[ ! "$build_tag" =~ ^[A-Za-z0-9_]+$ ]]; then
  echo "SHU_LLF_AB_BUILD_TAG must contain only letters, digits and underscores" >&2
  exit 1
fi
if [[ -e "$output_root" ]]; then
  echo "$output_root already exists; select a new SHU_LLF_AB_ROOT" >&2
  exit 1
fi
if [[ ! -d "$reference_root" ]]; then
  echo "diagnostic reference directory does not exist: $reference_root" >&2
  exit 1
fi

mkdir -p "$output_root/build" \
  "$output_root/source_snapshot/case" \
  "$output_root/source_snapshot/solver"

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

hash_build_sources() {
  (
    cd "$solver_root"
    find src \
      -type f \
      \( -name '*.h' -o -name '*.H' -o -name '*.cpp' -o -name '*.cc' \
         -o -name '*.c' -o -name '*.mk' -o -name 'Make.*' \) \
      -print0 \
      | sort -z \
      | xargs -0 sha256sum
    sha256sum \
      exm/numerics/shuosher/GNUmakefile \
      exm/numerics/shuosher/Make.package \
      exm/numerics/shuosher/prob.h \
      exm/numerics/shuosher/inputs.thesis
  )
}

run_case() {
  local mode=$1
  local cells=$2
  local shock_llf=$3
  local destination="$output_root/$mode/N${cells}"
  local executable="$output_root/build/$(basename "$archived_executable")"
  local command=(
    mpirun
    --oversubscribe
    --bind-to
    none
    -np
    "$mpi_ranks"
    "$executable"
    "$case_root/inputs.thesis"
    "amr.n_cell=$cells"
    "amr.max_grid_size=8192"
    "amr.plot_files_output=1"
    "amr.plot_file=$destination/plt"
    "amr.checkpoint_files_output=0"
    "cns.nstep_screen_output=200"
    "cns.afd_correction=1"
    "cns.afd_shock_llf=$shock_llf"
    "cns.afd_smoothness_threshold=0.08"
    "cns.afd_shock_pressure_jump=0.01"
    "cns.afd_shock_compression=0.001"
    "cns.afd_teno_cutoff=1.0e-4"
  )
  mkdir -p "$destination"
  write_command "$destination/run_command.txt" "${command[@]}"
  printf 'cns.afd_shock_llf=%d\n' "$shock_llf" \
    >"$destination/variant_argument.txt"
  if ! run_logged_command \
      "$destination/run.log" "$destination/run_status.txt" \
      "${command[@]}"; then
    echo "run failed; see $destination/run.log" >&2
    return 1
  fi
}

hash_build_sources >"$output_root/source_hashes_before_build.txt"
cp -a "$solver_root/src/." "$output_root/source_snapshot/solver/"
cp "$case_root/GNUmakefile" \
  "$case_root/Make.package" \
  "$case_root/prob.h" \
  "$case_root/inputs.thesis" \
  "$case_root/analyze_shock_llf_ab.py" \
  "$case_root/run_shock_llf_ab.sh" \
  "$output_root/source_snapshot/case/"

git -C "$solver_root" status --short \
  >"$output_root/source_snapshot/git_status_short.txt" 2>&1 || true
git -C "$solver_root" diff --binary \
  >"$output_root/source_snapshot/git_worktree.patch" 2>&1 || true
git -C "$solver_root" diff --cached --binary \
  >"$output_root/source_snapshot/git_index.patch" 2>&1 || true
git -C "$solver_root" submodule status --recursive \
  >"$output_root/source_snapshot/git_submodules.txt" 2>&1 || true

{
  printf 'created_utc='
  date -u +%Y-%m-%dT%H:%M:%SZ
  printf 'case=Shu-Osher shock-entropy-wave interaction\n'
  printf 'study=same-source same-binary cns.afd_shock_llf A/B\n'
  printf 'method=AFD-HLLC-WENO-Z5\n'
  printf 'variants=0 1\n'
  printf 'cells=%s\n' "${cells_list[*]}"
  printf 'stop_time=1.8\n'
  printf 'cfl=0.30\n'
  printf 'time_integrator=SSPRK(3,3)\n'
  printf 'afd_correction=1\n'
  printf 'afd_smoothness_threshold=0.08\n'
  printf 'afd_shock_pressure_jump=0.01\n'
  printf 'afd_shock_compression=0.001\n'
  printf 'afd_teno_cutoff=1.0e-4\n'
  printf 'variant_only_argument=cns.afd_shock_llf\n'
  printf 'mpi_ranks=%s\n' "$mpi_ranks"
  printf 'build_jobs=%s\n' "$build_jobs"
  printf 'build_tag=%s\n' "$build_tag"
  printf 'hostname=%s\n' "$(hostname)"
  printf 'solver_root=%s\n' "$solver_root"
  printf 'output_root=%s\n' "$output_root"
  printf 'diagnostic_reference=%s\n' "$reference_root"
  printf 'diagnostic_reference_is_exact=false\n'
  printf 'git_commit=%s\n' \
    "$(git -C "$solver_root" rev-parse HEAD 2>/dev/null || printf unknown)"
  printf 'git_describe=%s\n' \
    "$(git -C "$solver_root" describe --always --dirty 2>/dev/null || printf unknown)"
  printf 'uname='
  uname -a
  printf 'compiler='
  "${CXX:-g++}" --version 2>/dev/null | head -n 1 || printf 'unknown\n'
} >"$output_root/run_manifest.txt"

cat >"$output_root/common_solver_arguments.txt" <<'EOF'
inputs.thesis
amr.n_cell=<N>
amr.max_grid_size=8192
amr.plot_files_output=1
amr.plot_file=<variant-specific-output>/plt
amr.checkpoint_files_output=0
cns.nstep_screen_output=200
cns.afd_correction=1
cns.afd_smoothness_threshold=0.08
cns.afd_shock_pressure_jump=0.01
cns.afd_shock_compression=0.001
cns.afd_teno_cutoff=1.0e-4
EOF
printf 'The paired numerical inputs differ only in cns.afd_shock_llf.\n' \
  >"$output_root/ab_invariant_statement.txt"

build_command=(
  make
  -j"$build_jobs"
  SHU_EULER_SCHEME=afd-hllc-wenoz5
  SHU_BUILD_TAG="$build_tag"
  USE_PELEPHYSICS=FALSE
  USE_EB=FALSE
  USE_GPIBM=FALSE
)
write_command "$output_root/build/build_command.txt" "${build_command[@]}"
if ! run_logged_command \
    "$output_root/build/build.log" "$output_root/build/build_status.txt" \
    "${build_command[@]}"; then
  echo "build failed; see $output_root/build/build.log" >&2
  exit 1
fi

built_executable="$case_root/main1d.gnu.MPI.shu_afd-hllc-wenoz5_${build_tag}.ex"
if [[ ! -x "$built_executable" ]]; then
  echo "expected executable was not produced: $built_executable" >&2
  exit 1
fi
archived_executable="$output_root/build/$(basename "$built_executable")"
cp "$built_executable" "$archived_executable"
sha256sum "$built_executable" >"$output_root/build/built_executable.sha256"
sha256sum "$archived_executable" >"$output_root/build/archived_executable.sha256"

hash_build_sources >"$output_root/source_hashes_after_build.txt"
if ! cmp -s \
    "$output_root/source_hashes_before_build.txt" \
    "$output_root/source_hashes_after_build.txt"; then
  echo "build sources changed during compilation; refusing to run A/B study" >&2
  diff -u \
    "$output_root/source_hashes_before_build.txt" \
    "$output_root/source_hashes_after_build.txt" \
    >"$output_root/source_hash_change.diff" || true
  exit 1
fi

for shock_llf in 0 1; do
  mode="shock_llf_${shock_llf}"
  smoke_dir="$output_root/$mode/smoke"
  mkdir -p "$smoke_dir"
  smoke_command=(
    mpirun
    --oversubscribe
    --bind-to
    none
    -np
    "$mpi_ranks"
    "$archived_executable"
    "$case_root/inputs.thesis"
    amr.n_cell=64
    max_step=20
    stop_time=0.02
    amr.plot_files_output=0
    amr.checkpoint_files_output=0
    cns.nstep_screen_output=20
    cns.afd_correction=1
    "cns.afd_shock_llf=$shock_llf"
    cns.afd_smoothness_threshold=0.08
    cns.afd_shock_pressure_jump=0.01
    cns.afd_shock_compression=0.001
    cns.afd_teno_cutoff=1.0e-4
  )
  write_command "$smoke_dir/run_command.txt" "${smoke_command[@]}"
  if ! run_logged_command \
      "$smoke_dir/run.log" "$smoke_dir/run_status.txt" \
      "${smoke_command[@]}"; then
    echo "smoke test failed; see $smoke_dir/run.log" >&2
    exit 1
  fi
done
touch "$output_root/all_smoke_tests_passed"

for cells in "${cells_list[@]}"; do
  run_case shock_llf_0 "$cells" 0
  run_case shock_llf_1 "$cells" 1
done
touch "$output_root/all_solver_runs_complete"

hash_build_sources >"$output_root/source_hashes_after_runs.txt"
if ! cmp -s \
    "$output_root/source_hashes_before_build.txt" \
    "$output_root/source_hashes_after_runs.txt"; then
  echo "build sources changed during solver runs; A/B result is not accepted" >&2
  diff -u \
    "$output_root/source_hashes_before_build.txt" \
    "$output_root/source_hashes_after_runs.txt" \
    >"$output_root/source_hash_change.diff" || true
  exit 1
fi

analysis_command=(
  python3
  "$case_root/analyze_shock_llf_ab.py"
  "$output_root"
  --cells
  "${cells_list[@]}"
  --reference
  "$reference_root"
)
write_command "$output_root/analysis_command.txt" "${analysis_command[@]}"
if ! MPLCONFIGDIR=/tmp/cerisse_shuosher_llf_ab_mpl \
    run_logged_command \
      "$output_root/analysis.log" "$output_root/analysis_status.txt" \
      "${analysis_command[@]}"; then
  echo "analysis failed; see $output_root/analysis.log" >&2
  exit 1
fi

find "$output_root/source_snapshot" -type f ! -name SHA256SUMS -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  >"$output_root/source_snapshot/SHA256SUMS"
printf 'finished_utc=' >>"$output_root/run_manifest.txt"
date -u +%Y-%m-%dT%H:%M:%SZ >>"$output_root/run_manifest.txt"
touch "$output_root/complete"
echo "Completed strict Shu--Osher shock_llf A/B study: $output_root"
