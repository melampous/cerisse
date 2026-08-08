#!/usr/bin/env bash
set -euo pipefail

case_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$case_root"

results_root=${CERISSE_RZ_NS_MMS_RESULTS_ROOT:-results/rz_ns_mms}
build_jobs=${CERISSE_BUILD_JOBS:-4}
mpi_ranks=${CERISSE_MPI_RANKS:-1}

if [[ -e "$results_root" ]]; then
  echo "$results_root already exists; select a new result root" >&2
  exit 1
fi
if [[ ! "$mpi_ranks" =~ ^[1-9][0-9]*$ ]]; then
  echo "CERISSE_MPI_RANKS must be a positive integer" >&2
  exit 2
fi

grids=(16 32 64 128)
steps=(8 32 128 512)
final_time=1.0e-2
dt_audit_n=64
dt_audit_steps=256

git_description=$(
  git describe --always --tags --dirty 2>/dev/null || printf unavailable
)
git_revision=$(git rev-parse HEAD 2>/dev/null || printf unavailable)
git_status_snapshot=$(
  git status --short --untracked-files=normal 2>&1 || true
)
mkdir -p "$results_root/source_snapshot"
cp GNUmakefile inputs prob.h analyze.py analyze_operator_difference.py \
  compare_dt_audit.py run_study.sh README.md \
  "$results_root/source_snapshot/"
printf '%s\n' "$git_status_snapshot" >"$results_root/git_status.txt"
sha256sum GNUmakefile inputs prob.h analyze.py \
  analyze_operator_difference.py compare_dt_audit.py run_study.sh README.md \
  >"$results_root/source.sha256"

{
  echo "case=full_axisymmetric_rz_navier_stokes_mms"
  echo "ibm=disabled"
  echo "state_semantics=point_sample"
  echo "inviscid_scheme=LLF-WENO-Z5"
  echo "viscous_scheme=second_order_centred_face_flux"
  echo "viscosity=0.01"
  echo "conductivity=0.02"
  echo "rz_euler_annular_average=0"
  echo "rz_euler_point_flux=0"
  echo "rz_euler_annular_pressure_consistency=0"
  echo "rz_euler_cell_average_deconvolution=0"
  echo "rz_llf_weno_pressure_operator=paired_metric_advection_and_pressure_gradient"
  echo "rz_viscous_annular_average=0"
  echo "rz_visc_hoop=1"
  echo "rz_pressure_split=0"
  echo "rhs_equation=U_t=-div_RZ(Fc)+div_RZ(Fv)+p_over_r-tau_tt_over_r+S"
  echo "mms_source_sign=S=U_t+div_RZ(Fc)-div_RZ(Fv)-p_over_r+tau_tt_over_r"
  echo "order_rk_zero_plot_semantics=assembled_rhs"
  echo "grid_sequence=16 32 64 128"
  echo "step_sequence=8 32 128 512"
  echo "time_step_scaling=dt_proportional_to_dx_squared"
  echo "dt_audit_grid=${dt_audit_n}"
  echo "dt_audit_steps=${dt_audit_steps}"
  echo "final_time=${final_time}"
  echo "git_description=${git_description}"
  echo "git_revision=${git_revision}"
  echo "git_status_snapshot=git_status.txt"
  echo "host=$(hostname)"
  echo "mpi_ranks=${mpi_ranks}"
  echo "mpi_launcher=mpirun --oversubscribe"
  echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} >"$results_root/manifest.txt"

make -j"$build_jobs" RZ_MMS_PHYSICS=navier-stokes \
  USE_PELEPHYSICS=FALSE USE_EB=FALSE USE_GPIBM=FALSE
executable=./main2d.gnu.MPI.mms_rz_ns.ex
if [[ ! -x "$executable" ]]; then
  echo "missing executable $executable" >&2
  exit 1
fi
sha256sum "$executable" >"$results_root/navier_stokes_executable.sha256"

rhs_dirs=()
state_dirs=()
for index in "${!grids[@]}"; do
  n=${grids[$index]}
  nsteps=${steps[$index]}
  dt=$(python3 -c "print(${final_time}/${nsteps})")

  rhs_dir="$results_root/rhs/N$n"
  mkdir -p "$rhs_dir"
  rhs_args=(
    "$executable"
    inputs
    max_step=1
    stop_time=1.0
    time_step=1.0e-5
    cns.order_rk=0
    cns.stages_rk=1
    "amr.n_cell=${n} ${n}"
    amr.plot_int=1
    "amr.plot_file=${rhs_dir}/plt"
  )
  printf '%q ' mpirun --oversubscribe -np "$mpi_ranks" "${rhs_args[@]}" \
    >"$rhs_dir/command.txt"
  printf '\n' >>"$rhs_dir/command.txt"
  echo "running R-Z NS MMS spatial RHS N=$n"
  mpirun --oversubscribe -np "$mpi_ranks" "${rhs_args[@]}" \
    >"$rhs_dir/run.log" 2>&1
  rhs_dirs+=("$rhs_dir")

  state_dir="$results_root/evolve/N$n"
  mkdir -p "$state_dir"
  state_args=(
    "$executable"
    inputs
    "max_step=${nsteps}"
    "stop_time=${final_time}"
    "time_step=${dt}"
    "amr.n_cell=${n} ${n}"
    "amr.plot_int=${nsteps}"
    "amr.plot_file=${state_dir}/plt"
  )
  printf '%q ' mpirun --oversubscribe -np "$mpi_ranks" \
    "${state_args[@]}" \
    >"$state_dir/command.txt"
  printf '\n' >>"$state_dir/command.txt"
  echo "running R-Z NS MMS evolution N=$n steps=$nsteps dt=$dt"
  mpirun --oversubscribe -np "$mpi_ranks" "${state_args[@]}" \
    >"$state_dir/run.log" 2>&1
  state_dirs+=("$state_dir")
done

dt_audit_dt=$(python3 -c "print(${final_time}/${dt_audit_steps})")
dt_audit_dir="$results_root/evolve_dt_half/N${dt_audit_n}"
mkdir -p "$dt_audit_dir"
dt_audit_args=(
  "$executable"
  inputs
  "max_step=${dt_audit_steps}"
  "stop_time=${final_time}"
  "time_step=${dt_audit_dt}"
  "amr.n_cell=${dt_audit_n} ${dt_audit_n}"
  "amr.plot_int=${dt_audit_steps}"
  "amr.plot_file=${dt_audit_dir}/plt"
)
printf '%q ' mpirun --oversubscribe -np "$mpi_ranks" \
  "${dt_audit_args[@]}" >"$dt_audit_dir/command.txt"
printf '\n' >>"$dt_audit_dir/command.txt"
echo "running R-Z NS MMS N=${dt_audit_n} half-dt audit"
mpirun --oversubscribe -np "$mpi_ranks" "${dt_audit_args[@]}" \
  >"$dt_audit_dir/run.log" 2>&1

python3 analyze.py "${rhs_dirs[@]}" --quantity rhs --exact-time 0 \
  --physics navier-stokes \
  --csv "$results_root/rhs_errors.csv" \
  >"$results_root/rhs_analysis.log" 2>&1
cat "$results_root/rhs_analysis.log"

python3 analyze.py "${state_dirs[@]}" --quantity state \
  --physics navier-stokes \
  --csv "$results_root/evolution_errors.csv" \
  >"$results_root/evolution_analysis.log" 2>&1
cat "$results_root/evolution_analysis.log"

python3 analyze.py "$dt_audit_dir" --quantity state \
  --physics navier-stokes \
  --csv "$results_root/evolution_dt_half_error.csv" \
  >"$results_root/evolution_dt_half_analysis.log" 2>&1
cat "$results_root/evolution_dt_half_analysis.log"

python3 compare_dt_audit.py \
  "$results_root/evolution_errors.csv" \
  "$results_root/evolution_dt_half_error.csv" \
  --grid "$dt_audit_n" \
  --csv "$results_root/evolution_dt_audit.csv" \
  >"$results_root/evolution_dt_audit.log" 2>&1
cat "$results_root/evolution_dt_audit.log"

make -j"$build_jobs" RZ_MMS_PHYSICS=euler \
  USE_PELEPHYSICS=FALSE USE_EB=FALSE USE_GPIBM=FALSE
euler_executable=./main2d.gnu.MPI.mms_rz_euler.ex
if [[ ! -x "$euler_executable" ]]; then
  echo "missing executable $euler_executable" >&2
  exit 1
fi
sha256sum "$euler_executable" >"$results_root/euler_executable.sha256"

euler_rhs_dirs=()
for n in "${grids[@]}"; do
  euler_rhs_dir="$results_root/euler_rhs/N$n"
  mkdir -p "$euler_rhs_dir"
  euler_rhs_args=(
    "$euler_executable"
    inputs
    max_step=1
    stop_time=1.0
    time_step=1.0e-5
    cns.order_rk=0
    cns.stages_rk=1
    "amr.n_cell=${n} ${n}"
    amr.plot_int=1
    "amr.plot_file=${euler_rhs_dir}/plt"
  )
  printf '%q ' mpirun --oversubscribe -np "$mpi_ranks" \
    "${euler_rhs_args[@]}" >"$euler_rhs_dir/command.txt"
  printf '\n' >>"$euler_rhs_dir/command.txt"
  echo "running paired R-Z Euler MMS spatial RHS N=$n"
  mpirun --oversubscribe -np "$mpi_ranks" "${euler_rhs_args[@]}" \
    >"$euler_rhs_dir/run.log" 2>&1
  euler_rhs_dirs+=("$euler_rhs_dir")
done

python3 analyze.py "${euler_rhs_dirs[@]}" --quantity rhs --exact-time 0 \
  --physics euler \
  --csv "$results_root/euler_rhs_errors.csv" \
  >"$results_root/euler_rhs_analysis.log" 2>&1
cat "$results_root/euler_rhs_analysis.log"

python3 analyze_operator_difference.py \
  --euler "${euler_rhs_dirs[@]}" \
  --navier-stokes "${rhs_dirs[@]}" \
  --csv "$results_root/viscous_operator_errors.csv" \
  >"$results_root/viscous_operator_analysis.log" 2>&1
cat "$results_root/viscous_operator_analysis.log"

echo "completed_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
  >>"$results_root/manifest.txt"
echo "completed R-Z Navier-Stokes MMS study: $results_root"
