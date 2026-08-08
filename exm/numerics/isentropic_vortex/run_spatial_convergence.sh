#!/usr/bin/env bash
set -euo pipefail

case_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$case_root"

scheme=afd-hllc-wenoz5
input_file=inputs.spatial_convergence
results_root=${CERISSE_IV_SPATIAL_RESULTS_ROOT:-results/afd_hllc_wenoz5_spatial_convergence}
build_jobs=${CERISSE_BUILD_JOBS:-4}
mpi_ranks=${CERISSE_MPI_RANKS:-4}
flow_angle_deg=${CERISSE_IV_FLOW_ANGLE_DEG:-45.0}

if [[ ! "$mpi_ranks" =~ ^[1-9][0-9]*$ ]]; then
  echo "CERISSE_MPI_RANKS must be a positive integer" >&2
  exit 2
fi
if [[ -e "$results_root" ]]; then
  echo "$results_root already exists; select a new result root" >&2
  exit 1
fi

# T = 0.1 L/u_inf. The step counts grow as N^2. Therefore dt scales as
# dx^2 and the global SSPRK(3,3) error scales as dx^6.
cell_counts=(40 80 160 320)
step_counts=(64 256 1024 4096)
time_steps=(
  8.99927581881704387e-07
  2.24981895470426097e-07
  5.62454738676065242e-08
  1.40613684669016310e-08
)
final_time=5.75953652404290800e-05

make -j"$build_jobs" IV_EULER_SCHEME="$scheme" \
  USE_PELEPHYSICS=FALSE USE_EB=FALSE USE_GPIBM=FALSE

executable="./main2d.gnu.MPI.iv_${scheme}.ex"
if [[ ! -x "$executable" ]]; then
  echo "expected executable $executable was not produced" >&2
  exit 1
fi

mkdir -p "$results_root/source_snapshot"
cp GNUmakefile prob.h analyze.py inputs.spatial_convergence \
  run_spatial_convergence.sh "$results_root/source_snapshot/"

git_revision=$(git rev-parse HEAD 2>/dev/null || printf 'unavailable')
{
  echo "case=isentropic_vortex"
  echo "study=low_temporal_error_spatial_convergence"
  echo "scheme=${scheme}"
  echo "grid_sequence=40 80 160 320"
  echo "step_sequence=64 256 1024 4096"
  echo "time_step_scaling=dt_proportional_to_dx_squared"
  echo "expected_ssprk3_temporal_error=O(dx^6)"
  echo "flow_through_fraction=0.1"
  echo "flow_angle_deg=${flow_angle_deg}"
  echo "final_time=${final_time}"
  echo "executable=${executable}"
  echo "executable_sha256=$(sha256sum "$executable" | awk '{print $1}')"
  echo "source_input=${input_file}"
  echo "source_input_sha256=$(sha256sum "$input_file" | awk '{print $1}')"
  echo "git_revision=${git_revision}"
  echo "host=$(hostname)"
  echo "mpi_ranks=${mpi_ranks}"
  echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "rk_order=3"
  echo "rk_stages=3"
  echo "cns.afd_correction=1"
  echo "cns.afd_shock_llf=1"
  echo "cns.afd_smoothness_threshold=0.08"
  echo "cns.afd_teno_cutoff=1.0e-4"
} >"$results_root/manifest.txt"

result_dirs=()
for index in "${!cell_counts[@]}"; do
  cells=${cell_counts[$index]}
  steps=${step_counts[$index]}
  dt=${time_steps[$index]}
  output_dir="$results_root/N${cells}"
  mkdir -p "$output_dir"

  run_args=(
    "$executable"
    "$input_file"
    "max_step=${steps}"
    "time_step=${dt}"
    "amr.n_cell=${cells} ${cells}"
    "amr.plot_file=${output_dir}/plt"
    "prob.flow_angle_deg=${flow_angle_deg}"
    "cns.afd_correction=1"
    "cns.afd_shock_llf=1"
    "cns.afd_smoothness_threshold=0.08"
    "cns.afd_teno_cutoff=1.0e-4"
  )

  {
    echo '#!/usr/bin/env bash'
    echo 'set -euo pipefail'
    printf 'cd %q\n' "$case_root"
    printf '%q ' mpirun -np "$mpi_ranks" "${run_args[@]}"
    printf '\n'
  } >"$output_dir/run_command.sh"
  chmod +x "$output_dir/run_command.sh"

  echo "running ${scheme} spatial convergence N=${cells}, steps=${steps}, dt=${dt}"
  start_seconds=$(date +%s)
  set +e
  mpirun -np "$mpi_ranks" "${run_args[@]}" >"$output_dir/run.log" 2>&1
  status=$?
  set -e
  finish_seconds=$(date +%s)
  {
    echo "exit_status=${status}"
    echo "wall_seconds=$((finish_seconds - start_seconds))"
    echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  } >"$output_dir/run_status.txt"
  if [[ "$status" -ne 0 ]]; then
    echo "N=${cells} failed; inspect $output_dir/run.log" >&2
    exit "$status"
  fi
  result_dirs+=("$output_dir")
done

python3 analyze.py "${result_dirs[@]}" --flow-angle-deg "$flow_angle_deg" \
  --csv "$results_root/errors.csv" 2>&1 | tee "$results_root/analysis.log"

echo "completed ${scheme} low-temporal-error study: $results_root/errors.csv"
