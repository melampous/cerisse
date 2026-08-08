#!/usr/bin/env bash
set -euo pipefail

case_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$case_root"

usage() {
  cat <<'EOF'
Usage:
  ./run_method_matrix.sh SCHEME [SUITE]

SCHEME:
  llf-wenoz5
  llf-teno5
  afd-hllc-wenoz5
  afd-hllc-teno5
  all

SUITE:
  uniform  Run the 40, 80 and 160 uniform-grid cases.
  amr      Run the 40, 80 and 160 base-grid AMR cases.
  all      Run both suites. This is the default.

The four schemes are built as separate executables. Results are written to
results/method_matrix/SCHEME/SUITE and an existing suite is never overwritten.

Optional environment variables:
  CERISSE_BUILD_JOBS  Number of parallel build jobs. Default: 4.
  CERISSE_MPI_RANKS   Run each case with mpirun -np N. Direct execution is
                      used when this variable is unset.
  CERISSE_IV_RESULTS_ROOT
                      Alternative root for the method-matrix results.
EOF
}

all_schemes=(
  llf-wenoz5
  llf-teno5
  afd-hllc-wenoz5
  afd-hllc-teno5
)

requested_scheme=${1:-}
requested_suite=${2:-all}
if [[ -z "$requested_scheme" ]]; then
  usage >&2
  exit 2
fi

case "$requested_suite" in
  uniform|amr|all) ;;
  *)
    echo "unknown suite: $requested_suite" >&2
    usage >&2
    exit 2
    ;;
esac

if [[ "$requested_scheme" == "all" ]]; then
  selected_schemes=("${all_schemes[@]}")
else
  selected_schemes=()
  for candidate in "${all_schemes[@]}"; do
    if [[ "$requested_scheme" == "$candidate" ]]; then
      selected_schemes=("$requested_scheme")
      break
    fi
  done
  if [[ ${#selected_schemes[@]} -eq 0 ]]; then
    echo "unknown scheme: $requested_scheme" >&2
    usage >&2
    exit 2
  fi
fi

if [[ "$requested_suite" == "all" ]]; then
  selected_suites=(uniform amr)
else
  selected_suites=("$requested_suite")
fi

results_root=${CERISSE_IV_RESULTS_ROOT:-results/method_matrix}
cell_counts=(40 80 160)
build_jobs=${CERISSE_BUILD_JOBS:-4}

run_prefix=()
if [[ -n "${CERISSE_MPI_RANKS:-}" ]]; then
  if [[ ! "$CERISSE_MPI_RANKS" =~ ^[1-9][0-9]*$ ]]; then
    echo "CERISSE_MPI_RANKS must be a positive integer" >&2
    exit 2
  fi
  run_prefix=(mpirun -np "$CERISSE_MPI_RANKS")
fi

# Check the complete request before starting so a stale result cannot leave an
# avoidable partial method matrix.
for scheme in "${selected_schemes[@]}"; do
  for suite in "${selected_suites[@]}"; do
    suite_root="${results_root}/${scheme}/${suite}"
    if [[ -e "$suite_root" ]]; then
      echo "$suite_root already exists; move it aside before rerunning" >&2
      exit 1
    fi
  done
done

git_revision=$(git rev-parse HEAD 2>/dev/null || printf 'unavailable')

run_suite() {
  local scheme=$1
  local suite=$2
  local executable=$3
  shift 3
  local afd_args=("$@")

  local input_file
  local analysis_xc0
  if [[ "$suite" == "uniform" ]]; then
    input_file=inputs.uniform
    analysis_xc0=0.05
  else
    input_file=inputs.amr
    analysis_xc0=0.025
  fi

  local suite_root="${results_root}/${scheme}/${suite}"
  mkdir -p "$suite_root"
  cp "$input_file" "${suite_root}/input.archived"
  mkdir -p "${suite_root}/source_snapshot"
  cp GNUmakefile prob.h analyze.py analyze_method_matrix.py \
    plot_flowfield_matrix.py run_method_matrix.sh "$input_file" \
    "${suite_root}/source_snapshot/"

  {
    echo "case=isentropic_vortex"
    echo "scheme=${scheme}"
    echo "suite=${suite}"
    echo "grid_sequence=40 80 160"
    echo "executable=${executable}"
    echo "executable_sha256=$(sha256sum "$executable" | awk '{print $1}')"
    echo "source_input=${input_file}"
    echo "source_input_sha256=$(sha256sum "$input_file" | awk '{print $1}')"
    echo "source_snapshot=${suite_root}/source_snapshot"
    echo "git_revision=${git_revision}"
    echo "host=$(hostname)"
    echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "rk_order=3"
    echo "rk_stages=3"
    echo "cfl=0.30"
    if [[ ${#afd_args[@]} -gt 0 ]]; then
      echo "cns.afd_correction=1"
      echo "cns.afd_shock_llf=1"
      echo "cns.afd_smoothness_threshold=0.08"
      echo "cns.afd_teno_cutoff=1.0e-4"
    else
      echo "afd_runtime_parameters=not_applicable"
    fi
  } >"${suite_root}/manifest.txt"

  local result_dirs=()
  local cells
  for cells in "${cell_counts[@]}"; do
    local output_dir="${suite_root}/N${cells}"
    mkdir -p "$output_dir"
    local run_args=(
      "$executable"
      "$input_file"
      "amr.n_cell=${cells} ${cells}"
      "amr.plot_file=${output_dir}/plt"
    )
    run_args+=("${afd_args[@]}")

    {
      echo '#!/usr/bin/env bash'
      echo 'set -euo pipefail'
      printf 'cd %q\n' "$case_root"
      printf '%q ' "${run_prefix[@]}" "${run_args[@]}"
      printf '\n'
    } >"${output_dir}/run_command.sh"
    chmod +x "${output_dir}/run_command.sh"

    echo "running ${scheme} ${suite} N=${cells}"
    "${run_prefix[@]}" "${run_args[@]}" >"${output_dir}/run.log" 2>&1
    result_dirs+=("$output_dir")
  done

  local analysis_args=(
    python3 analyze.py
    "${result_dirs[@]}"
  )
  if [[ "$suite" == "amr" ]]; then
    analysis_args+=(--xc0 "$analysis_xc0")
  fi
  analysis_args+=(--csv "${suite_root}/errors.csv")
  "${analysis_args[@]}" 2>&1 | tee "${suite_root}/analysis.log"
  echo "completed ${scheme} ${suite}: ${suite_root}/errors.csv"
}

for scheme in "${selected_schemes[@]}"; do
  echo "building ${scheme}"
  make -j"$build_jobs" IV_EULER_SCHEME="$scheme" \
    USE_PELEPHYSICS=FALSE USE_EB=FALSE USE_GPIBM=FALSE

  executable="./main2d.gnu.MPI.iv_${scheme}.ex"
  if [[ ! -x "$executable" ]]; then
    echo "expected executable $executable was not produced" >&2
    exit 1
  fi

  afd_args=()
  case "$scheme" in
    afd-hllc-wenoz5|afd-hllc-teno5)
      afd_args=(
        cns.afd_correction=1
        cns.afd_shock_llf=1
        cns.afd_smoothness_threshold=0.08
        cns.afd_teno_cutoff=1.0e-4
      )
      ;;
  esac

  for suite in "${selected_suites[@]}"; do
    run_suite "$scheme" "$suite" "$executable" "${afd_args[@]}"
  done
done
