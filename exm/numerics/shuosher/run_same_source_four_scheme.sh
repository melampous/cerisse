#!/usr/bin/env bash
set -euo pipefail

# Rebuild and run four Shu--Osher methods from one source state and one input.
# The unique build tag prevents reuse of executables from earlier studies.

case_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
solver_root=$(cd "$case_root/../../.." && pwd)
cd "$case_root"

schemes=(
  llf-wenoz5
  llf-teno5
  afd-hllc-wenoz5
  skew
)
cells_list=(200 400 800)
build_jobs=${CERISSE_BUILD_JOBS:-8}
python_command=${PYTHON:-python3}
skew_c2=1.5
skew_c4=0.016
skew_c2_tag=1_5
skew_c4_tag=0_016
build_tag=${SHU_BUILD_TAG:-same_source_20260728_r1}
output_root=${SHU_SAME_SOURCE_ROOT:-"$case_root/results/y9000x_same_source_four_scheme_20260728_r1"}
if [[ "$output_root" != /* ]]; then
  output_root="$case_root/$output_root"
fi

if [[ ! "$build_tag" =~ ^[A-Za-z0-9_]+$ ]]; then
  echo "SHU_BUILD_TAG must contain only letters, digits and underscores" >&2
  exit 1
fi

run_prefix=()
if [[ -n "${CERISSE_RUN_PREFIX:-}" ]]; then
  read -r -a run_prefix <<<"${CERISSE_RUN_PREFIX}"
fi

if [[ -e "$output_root" ]]; then
  echo "$output_root already exists; select a new SHU_SAME_SOURCE_ROOT" >&2
  exit 1
fi
mkdir -p "$output_root/source_snapshot/src_rhs"

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
  printf 'build_tag=%s\n' "$build_tag"
  printf 'schemes=%s\n' "${schemes[*]}"
  printf 'cells=%s\n' "${cells_list[*]}"
  printf 'stop_time=1.8\n'
  printf 'cfl=0.30\n'
  printf 'rk=SSPRK3\n'
  printf 'skew_order=4\n'
  printf 'skew_C2=%s\n' "$skew_c2"
  printf 'skew_C4=%s\n' "$skew_c4"
  printf 'afd_correction=1\n'
  printf 'afd_shock_llf=1\n'
  printf 'afd_smoothness_threshold=0.08\n'
  printf 'afd_teno_cutoff=1.0e-4\n'
  printf 'run_prefix=%s\n' "${CERISSE_RUN_PREFIX:-<direct>}"
  printf 'build_jobs=%s\n' "$build_jobs"
  printf 'git_commit=%s\n' "$(git -C "$solver_root" rev-parse HEAD 2>/dev/null || printf unknown)"
  printf 'inputs_sha256=%s\n' "$(sha256sum "$case_root/inputs.thesis" | awk '{print $1}')"
  printf 'prob_header_sha256=%s\n' "$(sha256sum "$case_root/prob.h" | awk '{print $1}')"
  printf 'skew_N6400_reference_archive=%s\n' \
    "$case_root/results/y9000x_skew_jst_c2_1_5_c4_0_016_20260728_r3"
  printf 'skew_N6400_used_for_other_schemes=false\n'
  printf 'formal_order_inferred=false\n'
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
  "$case_root/analyze.py" \
  "$case_root/analyze_same_source_four_scheme.py" \
  "$case_root/run_same_source_four_scheme.sh" \
  "$output_root/source_snapshot/"
cp "$solver_root/src/rhs/RHS.h" \
  "$solver_root/src/rhs/Weno.h" \
  "$solver_root/src/rhs/Afd.h" \
  "$solver_root/src/rhs/Riemann.h" \
  "$solver_root/src/rhs/Skew.h" \
  "$output_root/source_snapshot/src_rhs/"
find "$output_root/source_snapshot" -type f ! -name SHA256SUMS -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  >"$output_root/source_snapshot/SHA256SUMS"

declare -A executables

for scheme in "${schemes[@]}"; do
  scheme_root="$output_root/$scheme"
  mkdir -p "$scheme_root"
  echo "Rebuilding $scheme"
  build_command=(
    make
    -j"$build_jobs"
    SHU_EULER_SCHEME="$scheme"
    SHU_SKEW_C2="$skew_c2"
    SHU_SKEW_C4="$skew_c4"
    SHU_BUILD_TAG="$build_tag"
    USE_PELEPHYSICS=FALSE
    USE_EB=FALSE
    USE_GPIBM=FALSE
  )
  write_command "$scheme_root/build_command.txt" "${build_command[@]}"
  if ! run_logged_command \
      "$scheme_root/build.log" "$scheme_root/build_status.txt" \
      "${build_command[@]}"; then
    echo "build failed; see $scheme_root/build.log" >&2
    exit 1
  fi

  if [[ "$scheme" == skew ]]; then
    executable="$case_root/main1d.gnu.MPI.shu_skew_c2_${skew_c2_tag}_c4_${skew_c4_tag}_${build_tag}.ex"
  else
    executable="$case_root/main1d.gnu.MPI.shu_${scheme}_${build_tag}.ex"
  fi
  if [[ ! -x "$executable" ]]; then
    echo "expected executable $executable was not produced" >&2
    exit 1
  fi
  executables["$scheme"]=$executable
  sha256sum "$executable" >"$scheme_root/executable.sha256"
  cp "$executable" "$scheme_root/"
  sha256sum "$scheme_root/$(basename "$executable")" \
    >"$scheme_root/archived_executable.sha256"
done
touch "$output_root/all_builds_complete"

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
  echo "Smoke testing $scheme"
  run_case "$output_root/$scheme/smoke" "$executable" \
    amr.n_cell=64 \
    max_step=20 \
    stop_time=0.02 \
    amr.plot_files_output=0 \
    amr.checkpoint_files_output=0 \
    cns.nstep_screen_output=20 \
    "${method_args[@]}"
done
touch "$output_root/all_smoke_tests_passed"

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

analysis_command=(
  "$python_command"
  "$case_root/analyze_same_source_four_scheme.py"
  "$output_root"
  --window-lo -1.5
  --window-hi 1.5
  --csv "$output_root/same_source_four_scheme_metrics.csv"
  --png "$output_root/same_source_four_scheme_density.png"
  --pdf "$output_root/same_source_four_scheme_density.pdf"
  --metadata "$output_root/analysis_metadata.json"
)
write_command "$output_root/analysis_command.txt" "${analysis_command[@]}"
if ! run_logged_command \
    "$output_root/analysis.log" "$output_root/analysis_status.txt" \
    "${analysis_command[@]}"; then
  echo "analysis failed; see $output_root/analysis.log" >&2
  exit 1
fi

touch "$output_root/matrix_complete"
echo "Completed same-source four-scheme matrix: $output_root"
