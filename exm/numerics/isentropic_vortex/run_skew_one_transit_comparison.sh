#!/usr/bin/env bash
set -euo pipefail

case_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
solver_root=$(cd "$case_root/../../.." && pwd)
cd "$case_root"

results_root=${CERISSE_IV_COMPARISON_RESULTS_ROOT:-results/y9000x_skew_one_transit_comparison_20260728}
build_jobs=${CERISSE_BUILD_JOBS:-4}
mpi_ranks=${CERISSE_MPI_RANKS:-4}
methods=(llf-wenoz5 llf-teno5 afd-hllc-wenoz5 skew-jst-o4)

if [[ ! "$build_jobs" =~ ^[1-9][0-9]*$ ]] ||
   [[ ! "$mpi_ranks" =~ ^[1-9][0-9]*$ ]]; then
  echo "build jobs and MPI ranks must be positive integers" >&2
  exit 2
fi
if [[ -e "$results_root" ]]; then
  echo "$results_root already exists; select a new result root" >&2
  exit 1
fi

mkdir -p "$results_root/source_snapshot/solver" "$results_root/build"
cp GNUmakefile prob.h analyze.py analyze_skew_validation.py \
  plot_skew_one_transit_comparison.py inputs.uniform \
  run_skew_one_transit_comparison.sh "$results_root/source_snapshot/"
cp "$solver_root/src/rhs/Skew.h" "$results_root/source_snapshot/solver/"
cp "$solver_root/src/rhs/Weno.h" "$results_root/source_snapshot/solver/"
cp "$solver_root/src/rhs/Afd.h" "$results_root/source_snapshot/solver/"
git status --short >"$results_root/source_snapshot/git_status.txt" || true
git diff -- src/rhs/Skew.h src/rhs/Weno.h src/rhs/Afd.h \
  exm/numerics/isentropic_vortex \
  >"$results_root/source_snapshot/worktree.patch" || true
(
  cd "$results_root/source_snapshot"
  find . -type f -print0 | sort -z | xargs -0 sha256sum
) >"$results_root/source_snapshot/SHA256SUMS"

git_revision=$(git rev-parse HEAD 2>/dev/null || printf 'unavailable')
{
  echo "case=isentropic_vortex"
  echo "study=same_source_one_transit_four_method_comparison"
  echo "methods=${methods[*]}"
  echo "grid_sequence=40 80 160"
  echo "final_time=5.759536524042908e-4"
  echo "cfl=0.30"
  echo "flow_angle_deg=0.0"
  echo "rk_order=3"
  echo "rk_stages=3"
  echo "host=$(hostname)"
  echo "mpi_ranks=${mpi_ranks}"
  echo "git_revision=${git_revision}"
  echo "source_snapshot=${results_root}/source_snapshot"
  echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} >"$results_root/manifest.txt"

for method in "${methods[@]}"; do
  method_root="$results_root/$method"
  executable="./main2d.gnu.MPI.iv_${method}.ex"
  mkdir -p "$method_root"

  echo "building ${method}"
  build_start=$(date +%s)
  set +e
  make -j"$build_jobs" IV_EULER_SCHEME="$method" \
    USE_PELEPHYSICS=FALSE USE_EB=FALSE USE_GPIBM=FALSE \
    >"$results_root/build/${method}.log" 2>&1
  status=$?
  set -e
  build_finish=$(date +%s)
  {
    echo "exit_status=${status}"
    echo "wall_seconds=$((build_finish - build_start))"
    echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  } >"$results_root/build/${method}.status"
  if [[ "$status" -ne 0 ]]; then
    echo "build failed for ${method}" >&2
    exit "$status"
  fi
  if [[ ! -x "$executable" ]]; then
    echo "expected executable ${executable} was not produced" >&2
    exit 1
  fi

  method_args=()
  if [[ "$method" == "afd-hllc-wenoz5" ]]; then
    method_args=(
      cns.afd_correction=1
      cns.afd_shock_llf=1
      cns.afd_smoothness_threshold=0.08
      cns.afd_teno_cutoff=1.0e-4
    )
  fi

  {
    echo "case=isentropic_vortex"
    echo "study=one_transit"
    echo "scheme=${method}"
    echo "executable=${executable}"
    echo "executable_sha256=$(sha256sum "$executable" | awk '{print $1}')"
    echo "grid_sequence=40 80 160"
    echo "final_time=5.759536524042908e-4"
    echo "cfl=0.30"
    echo "flow_angle_deg=0.0"
    echo "rk_order=3"
    echo "rk_stages=3"
    if [[ "$method" == "skew-jst-o4" ]]; then
      echo "skew_order=4"
      echo "jst_dissipation=true"
      echo "C2skew=1.5"
      echo "C4skew=0.016"
    else
      echo "skew_parameters=not_applicable"
    fi
    if [[ "$method" == "afd-hllc-wenoz5" ]]; then
      echo "cns.afd_correction=1"
      echo "cns.afd_shock_llf=1"
      echo "cns.afd_smoothness_threshold=0.08"
      echo "cns.afd_teno_cutoff=1.0e-4"
    else
      echo "afd_runtime_parameters=not_applicable"
    fi
  } >"$method_root/manifest.txt"

  result_dirs=()
  for cells in 40 80 160; do
    output_dir="$method_root/N${cells}"
    mkdir -p "$output_dir"
    {
      echo '#!/usr/bin/env bash'
      echo 'set -euo pipefail'
      printf 'cd %q\n' "$case_root"
      printf '%q ' mpirun -np "$mpi_ranks" "$executable" inputs.uniform \
        "amr.n_cell=${cells} ${cells}" \
        "amr.plot_file=${output_dir}/plt" \
        "${method_args[@]}"
      printf '\n'
    } >"$output_dir/run_command.sh"
    chmod +x "$output_dir/run_command.sh"

    echo "running ${method} one-transit N=${cells}"
    run_start=$(date +%s)
    set +e
    mpirun -np "$mpi_ranks" "$executable" inputs.uniform \
      "amr.n_cell=${cells} ${cells}" \
      "amr.plot_file=${output_dir}/plt" \
      "${method_args[@]}" \
      >"$output_dir/run.log" 2>&1
    status=$?
    set -e
    run_finish=$(date +%s)
    {
      echo "exit_status=${status}"
      echo "wall_seconds=$((run_finish - run_start))"
      echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    } >"$output_dir/run_status.txt"
    if [[ "$status" -ne 0 ]]; then
      echo "run failed for ${method} N=${cells}" >&2
      exit "$status"
    fi
    result_dirs+=("$output_dir")
  done

  MPLCONFIGDIR=/tmp/cerisse-skew-mpl python3 analyze_skew_validation.py \
    "${result_dirs[@]}" --csv "$method_root/errors.csv" \
    >"$method_root/analysis.log" 2>&1
  cat "$method_root/analysis.log"
done

MPLCONFIGDIR=/tmp/cerisse-skew-mpl python3 plot_skew_one_transit_comparison.py \
  "$results_root" >"$results_root/plot.log" 2>&1
cat "$results_root/plot.log"
echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)" >>"$results_root/manifest.txt"
echo "completed one-transit four-method comparison: $results_root"
