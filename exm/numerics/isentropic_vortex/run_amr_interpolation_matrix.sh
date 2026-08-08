#!/usr/bin/env bash
set -euo pipefail

case_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$case_root"

executable=${CERISSE_IV_EXECUTABLE:-./main2d.gnu.MPI.iv_llf-wenoz5.ex}
results_root=${CERISSE_IV_AMR_INTERP_RESULTS_ROOT:-results/amr_interpolation_matrix}
mpi_ranks=${CERISSE_MPI_RANKS:-1}
reflux=${CERISSE_IV_REFLUX:-1}
llf_weno_epsilon_relative=${CERISSE_LLF_WENO_EPSILON_RELATIVE:-}

if [[ ! "$mpi_ranks" =~ ^[1-9][0-9]*$ ]]; then
  echo "CERISSE_MPI_RANKS must be a positive integer" >&2
  exit 2
fi
if [[ "$reflux" != 0 && "$reflux" != 1 ]]; then
  echo "CERISSE_IV_REFLUX must be 0 or 1" >&2
  exit 2
fi
if [[ ! -x "$executable" ]]; then
  echo "missing executable: $executable" >&2
  exit 1
fi
if [[ -e "$results_root" ]]; then
  echo "$results_root already exists; choose a new result root" >&2
  exit 1
fi

read -r -a interpolations <<<"${CERISSE_IV_INTERPOLATIONS:-linear conservative_quartic}"
read -r -a subcycling_modes <<<"${CERISSE_IV_SUBCYCLING_MODES:-Auto None}"
read -r -a cell_counts <<<"${CERISSE_IV_CELL_COUNTS:-40 80 160}"

mkdir -p "$results_root/source_snapshot/solver"
cp inputs.amr prob.h analyze.py run_amr_interpolation_matrix.sh \
  "$results_root/source_snapshot/"
cp ../../../src/CNS.cpp ../../../src/CNS.h ../../../src/set/CNS_setup.cpp \
  ../../../src/tim/advance.cpp ../../../src/tim/compute_rhs.cpp \
  ../../../src/rhs/Weno.h "$results_root/source_snapshot/solver/"

{
  echo "case=isentropic_vortex_amr_interpolation_matrix"
  echo "executable=$executable"
  echo "executable_sha256=$(sha256sum "$executable" | awk '{print $1}')"
  echo "git_revision=$(git rev-parse HEAD 2>/dev/null || printf unavailable)"
  echo "host=$(hostname)"
  echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "grid_sequence=${cell_counts[*]}"
  echo "interpolations=${interpolations[*]}"
  echo "subcycling_modes=${subcycling_modes[*]}"
  echo "reflux=$reflux"
  echo "llf_weno_epsilon_relative=${llf_weno_epsilon_relative:-legacy}"
  echo "mpi_ranks=$mpi_ranks"
} >"$results_root/manifest.txt"

run_prefix=()
if ((mpi_ranks > 1)); then
  run_prefix=(mpirun --oversubscribe -np "$mpi_ranks")
fi

for interpolation in "${interpolations[@]}"; do
  for subcycling in "${subcycling_modes[@]}"; do
    suite_root="$results_root/${interpolation}_${subcycling}"
    mkdir -p "$suite_root"
    result_dirs=()
    for cells in "${cell_counts[@]}"; do
      output_dir="$suite_root/N${cells}"
      mkdir -p "$output_dir"
      run_args=(
        "$executable"
        inputs.amr
        "amr.n_cell=${cells} ${cells}"
        "amr.subcycling_mode=${subcycling}"
        "cns.amr_state_interp=${interpolation}"
        "cns.do_reflux=${reflux}"
        "amr.plot_file=${output_dir}/plt"
      )
      if [[ -n "$llf_weno_epsilon_relative" ]]; then
        run_args+=(
          "cns.llf_weno_epsilon_relative=${llf_weno_epsilon_relative}"
        )
      fi
      {
        echo '#!/usr/bin/env bash'
        echo 'set -euo pipefail'
        printf 'cd %q\n' "$case_root"
        printf '%q ' "${run_prefix[@]}" "${run_args[@]}"
        printf '\n'
      } >"$output_dir/run_command.sh"
      chmod +x "$output_dir/run_command.sh"
      echo "running ${interpolation} ${subcycling} N=${cells}"
      "${run_prefix[@]}" "${run_args[@]}" >"$output_dir/run.log" 2>&1
      result_dirs+=("$output_dir")
    done
    python3 analyze.py "${result_dirs[@]}" --xc0 0.025 \
      --csv "$suite_root/errors.csv" >"$suite_root/analysis.log" 2>&1
  done
done

echo "completed AMR interpolation matrix: $results_root"
