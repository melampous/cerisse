#!/usr/bin/env bash
set -euo pipefail

case_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$case_root"

results_root=${CERISSE_MMS_RESULTS_ROOT:-results/cartesian_mms}
build_jobs=${CERISSE_BUILD_JOBS:-4}
mpi_ranks=${CERISSE_MPI_RANKS:-1}
use_mpi=${CERISSE_USE_MPI:-TRUE}
euler_scheme=${CERISSE_MMS_EULER_SCHEME:-llf-wenoz5}
physics_selection=${CERISSE_MMS_PHYSICS:-both}
viscous_order=${CERISSE_MMS_VISCOUS_ORDER:-2}
transport_model=${CERISSE_MMS_TRANSPORT_MODEL:-constant}
llf_weno_epsilon_relative=${CERISSE_MMS_LLF_WENO_EPSILON_RELATIVE:-}

case "$euler_scheme" in
  llf-wenoz5) executable_scheme_suffix="" ;;
  llf-teno5) executable_scheme_suffix=".llf_teno5" ;;
  llf-teno6) executable_scheme_suffix=".llf_teno6" ;;
  afd-hllc-wenoz5) executable_scheme_suffix=".afd_hllc_wenoz5" ;;
  afd-hllc-teno5) executable_scheme_suffix=".afd_hllc_teno5" ;;
  llf-wenoz5-old) executable_scheme_suffix=".llf_wenoz5_old" ;;
  llf-teno5-old) executable_scheme_suffix=".llf_teno5_old" ;;
  llf-teno6-old) executable_scheme_suffix=".llf_teno6_old" ;;
  *)
    echo "unsupported CERISSE_MMS_EULER_SCHEME=$euler_scheme" >&2
    exit 2
    ;;
esac

case "$viscous_order" in
  2) viscous_order_suffix="" ;;
  4|6) viscous_order_suffix=".visc${viscous_order}" ;;
  *)
    echo "CERISSE_MMS_VISCOUS_ORDER must be one of: 2, 4, 6" >&2
    exit 2
    ;;
esac

case "$transport_model" in
  constant) transport_suffix="" ;;
  sutherland) transport_suffix=".sutherland" ;;
  *)
    echo "CERISSE_MMS_TRANSPORT_MODEL must be one of: constant, sutherland" >&2
    exit 2
    ;;
esac

case "$physics_selection" in
  both) physics_cases=(euler navier-stokes) ;;
  euler) physics_cases=(euler) ;;
  navier-stokes) physics_cases=(navier-stokes) ;;
  *)
    echo "unsupported CERISSE_MMS_PHYSICS=$physics_selection" >&2
    exit 2
    ;;
esac

scheme_run_args=()
if [[ -n "$llf_weno_epsilon_relative" ]]; then
  scheme_run_args+=(
    "cns.llf_weno_epsilon_relative=${llf_weno_epsilon_relative}"
  )
fi
if [[ "$euler_scheme" == afd-hllc-* ]]; then
  scheme_run_args+=(
    "cns.afd_correction=1"
    "cns.afd_shock_llf=1"
    "cns.afd_smoothness_threshold=0.08"
    "cns.afd_teno_cutoff=1.0e-4"
  )
fi
if [[ -e "$results_root" ]]; then
  echo "$results_root already exists; select a new result root" >&2
  exit 1
fi
if [[ ! "$mpi_ranks" =~ ^[1-9][0-9]*$ ]]; then
  echo "CERISSE_MPI_RANKS must be a positive integer" >&2
  exit 2
fi
case "$use_mpi" in
  TRUE|FALSE) ;;
  *)
    echo "CERISSE_USE_MPI must be TRUE or FALSE" >&2
    exit 2
    ;;
esac

grids=(16 32 64 128)
steps=(8 32 128 512)
final_time=1.0e-2

mkdir -p "$results_root/source_snapshot/solver"
cp GNUmakefile inputs prob.h analyze.py run_study.sh \
  "$results_root/source_snapshot/"
cp ../../../src/Make.CNS ../../../src/rhs/Afd.h \
  ../../../src/rhs/AfdIBM.h ../../../src/rhs/AfdIBMDrivers.h \
  ../../../src/rhs/IBMSharedGPFluxUtils.h \
  ../../../src/rhs/Weno.h ../../../src/rhs/Weno_old.h \
  ../../../src/rhs/Riemann.h \
  ../../../src/rhs/RHS.h ../../../src/rhs/viscous.h \
  ../../../src/rhs/diff_ops.H "$results_root/source_snapshot/solver/"

git_revision=$(git rev-parse HEAD 2>/dev/null || printf unavailable)
git_description=$(git describe --always --dirty 2>/dev/null || printf unavailable)
git status --porcelain=v1 >"$results_root/git_status_porcelain.txt"
git diff -- ../../../src/Make.CNS ../../../src/rhs/Afd.h \
  ../../../src/rhs/AfdIBM.h ../../../src/rhs/AfdIBMDrivers.h \
  ../../../src/rhs/IBMSharedGPFluxUtils.h \
  ../../../src/rhs/Weno.h ../../../src/rhs/Weno_old.h \
  ../../../src/rhs/Riemann.h ../../../src/rhs/RHS.h \
  ../../../src/rhs/viscous.h ../../../src/rhs/diff_ops.H \
  >"$results_root/relevant_solver_worktree.patch"
sha256sum "$results_root/source_snapshot/GNUmakefile" \
  "$results_root/source_snapshot/inputs" \
  "$results_root/source_snapshot/prob.h" \
  "$results_root/source_snapshot/analyze.py" \
  "$results_root/source_snapshot/run_study.sh" \
  "$results_root/source_snapshot/solver/Make.CNS" \
  "$results_root/source_snapshot/solver/Afd.h" \
  "$results_root/source_snapshot/solver/AfdIBM.h" \
  "$results_root/source_snapshot/solver/AfdIBMDrivers.h" \
  "$results_root/source_snapshot/solver/IBMSharedGPFluxUtils.h" \
  "$results_root/source_snapshot/solver/Weno.h" \
  "$results_root/source_snapshot/solver/Weno_old.h" \
  "$results_root/source_snapshot/solver/Riemann.h" \
  "$results_root/source_snapshot/solver/RHS.h" \
  "$results_root/source_snapshot/solver/viscous.h" \
  "$results_root/source_snapshot/solver/diff_ops.H" \
  >"$results_root/source_snapshot.sha256"
{
  echo "case=cartesian_time_dependent_mms"
  echo "state_semantics=point_sample"
  echo "physics=${physics_cases[*]}"
  echo "inviscid_scheme=${euler_scheme}"
  echo "viscous_scheme=centred_face_flux"
  echo "viscous_order=${viscous_order}"
  echo "transport_model=${transport_model}"
  if [[ "$transport_model" == sutherland ]]; then
    echo "transport_temperature_reference=1.0"
    echo "transport_viscosity_reference=0.01"
    echo "transport_viscosity_sutherland_constant=0.5"
    echo "transport_conductivity_reference=0.02"
    echo "transport_conductivity_sutherland_constant=0.8"
  fi
  if [[ -n "$llf_weno_epsilon_relative" ]]; then
    echo "llf_weno_epsilon_relative=${llf_weno_epsilon_relative}"
  fi
  echo "grid_sequence=16 32 64 128"
  echo "step_sequence=8 32 128 512"
  echo "time_step_scaling=dt_proportional_to_dx_squared"
  echo "expected_ssprk3_temporal_error=O(dx^6)"
  echo "final_time=${final_time}"
  echo "git_revision=${git_revision}"
  echo "git_description=${git_description}"
  if [[ "$euler_scheme" == afd-hllc-* ]]; then
    echo "afd_correction=1"
    echo "afd_shock_llf=1"
    echo "afd_smoothness_threshold=0.08"
    echo "afd_teno_cutoff=1.0e-4"
  fi
  echo "host=$(hostname)"
  echo "use_mpi=${use_mpi}"
  echo "mpi_ranks=${mpi_ranks}"
  echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} >"$results_root/manifest.txt"

for physics in "${physics_cases[@]}"; do
  make -j"$build_jobs" MMS_PHYSICS="$physics" \
    MMS_EULER_SCHEME="$euler_scheme" \
    MMS_VISCOUS_ORDER="$viscous_order" \
    MMS_TRANSPORT_MODEL="$transport_model" \
    USE_MPI="$use_mpi" USE_PELEPHYSICS=FALSE USE_EB=FALSE USE_GPIBM=FALSE
  suffix=${physics//-/_}
  mpi_suffix=""
  if [[ "$use_mpi" == TRUE ]]; then
    mpi_suffix=".MPI"
  fi
  executable="./main2d.gnu${mpi_suffix}.mms_${suffix}${executable_scheme_suffix}${viscous_order_suffix}${transport_suffix}.ex"
  if [[ ! -x "$executable" ]]; then
    echo "missing executable $executable" >&2
    exit 1
  fi
  sha256sum "$executable" >"$results_root/${suffix}_executable.sha256"

  result_dirs=()
  for index in "${!grids[@]}"; do
    n=${grids[$index]}
    nsteps=${steps[$index]}
    dt=$(python3 -c "print(${final_time}/${nsteps})")
    output_dir="$results_root/$suffix/N$n"
    mkdir -p "$output_dir"
    run_args=(
      "$executable"
      inputs
      "max_step=${nsteps}"
      "stop_time=${final_time}"
      "time_step=${dt}"
      "amr.n_cell=${n} ${n}"
      "amr.plot_int=${nsteps}"
      "amr.plot_file=${output_dir}/plt"
      "${scheme_run_args[@]}"
    )
    launch_command=("${run_args[@]}")
    if [[ "$use_mpi" == TRUE ]]; then
      launch_command=(mpirun --oversubscribe -np "$mpi_ranks" "${run_args[@]}")
    fi
    printf '%q ' "${launch_command[@]}" >"$output_dir/command.txt"
    printf '\n' >>"$output_dir/command.txt"
    echo "running $physics MMS N=$n steps=$nsteps dt=$dt"
    start=$(date +%s)
    set +e
    "${launch_command[@]}" >"$output_dir/run.log" 2>&1
    status=$?
    set -e
    {
      echo "exit_status=$status"
      echo "wall_seconds=$(($(date +%s) - start))"
      echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    } >"$output_dir/status.txt"
    if [[ "$status" -ne 0 ]]; then
      echo "$physics N=$n failed; inspect $output_dir/run.log" >&2
      exit "$status"
    fi
    result_dirs+=("$output_dir")
  done

  python3 analyze.py "${result_dirs[@]}" --physics "$physics" \
    --csv "$results_root/${suffix}_errors.csv" \
    >"$results_root/${suffix}_analysis.log" 2>&1
  cat "$results_root/${suffix}_analysis.log"
done

echo "completed Cartesian MMS study: $results_root"
