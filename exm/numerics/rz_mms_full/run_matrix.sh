#!/usr/bin/env bash
set -euo pipefail

binary=$(realpath "${1:?usage: run_matrix.sh BINARY INPUTS RESULT_ROOT [SOLVER_OVERRIDE ...]}")
inputs=$(realpath "${2:?usage: run_matrix.sh BINARY INPUTS RESULT_ROOT [SOLVER_OVERRIDE ...]}")
result_root=${3:?usage: run_matrix.sh BINARY INPUTS RESULT_ROOT [SOLVER_OVERRIDE ...]}
global_overrides=("${@:4}")
mpi_root=${MPI_ROOT:-/home/melampous/cerisse_rz_bench_20260720/ompi_user/root/usr}
run_serial=${RUN_SERIAL:-0}

fatal_pattern='MPI_ABORT|amrex::Abort|SIG(SEGV|FPE|ABRT|BUS)|Segmentation fault|Floating point exception|non-finite|nonfinite|(^|[^[:alnum:]_])(nan|[+-]?inf(inity)?)([^[:alnum:]_]|$)'

die()
{
    echo "[RZ-MMS-FULL-MATRIX-FAIL] $*" >&2
    exit 2
}

[[ -x ${binary} && ! -L ${binary} ]] || die "missing executable: ${binary}"
[[ -f ${inputs} && ! -L ${inputs} ]] || die "missing inputs: ${inputs}"
[[ ! -e ${result_root} ]] || die "refusing to reuse result root: ${result_root}"

if [[ ${run_serial} != 0 && ${run_serial} != 1 ]]; then
    die "RUN_SERIAL must be 0 or 1"
fi
if [[ ${run_serial} == 0 ]]; then
    [[ -x ${mpi_root}/bin/orterun ]] || die "missing MPI launcher under ${mpi_root}"
fi

mkdir -p "${result_root}/runs"
sha256sum "${binary}" "${inputs}" >"${result_root}/RUN_INPUTS.sha256"
printf '%q ' "${global_overrides[@]}" >"${result_root}/GLOBAL_OVERRIDES.txt"
printf '\n' >>"${result_root}/GLOBAL_OVERRIDES.txt"
date -u +'%Y-%m-%dT%H:%M:%SZ' >"${result_root}/STARTED_UTC"

export OMP_NUM_THREADS=1
export LC_ALL=C

launcher=()
execution_mode=serial
if [[ ${run_serial} == 0 ]]; then
    export OPAL_PREFIX=${mpi_root}
    export PATH=${mpi_root}/bin:${PATH}
    export LD_LIBRARY_PATH=${mpi_root}/lib/x86_64-linux-gnu:${mpi_root}/lib:${LD_LIBRARY_PATH:-}
    launcher=(
        "${mpi_root}/bin/orterun"
        --mca pml ob1
        --mca btl self,vader
        --mca mca_base_component_show_load_errors 0
        --bind-to none
        -np 1
    )
    execution_mode=mpi_np1
fi

run_case()
{
    local label=$1
    local nr=$2
    local order=$3
    local stages=$4
    local max_steps=$5
    local stop=$6
    local plot_interval=$7
    shift 7
    local stage=${result_root}/runs/${label}
    local -a command=(
        "${launcher[@]}" "${binary}" "${inputs}"
        "amr.n_cell=${nr} ${nr}"
        "cns.order_rk=${order}"
        "cns.stages_rk=${stages}"
        "max_step=${max_steps}"
        "stop_time=${stop}"
        "amr.plot_int=${plot_interval}"
        amr.plot_file=output/plt
        amr.checkpoint_files_output=0
        cns.rz_euler_annular_average=1
        cns.rz_euler_point_flux=1
        cns.rz_euler_annular_pressure_consistency=1
        cns.rz_euler_cell_average_deconvolution=1
        cns.rz_pressure_split=0
        "$@"
        "${global_overrides[@]}"
    )

    mkdir -p "${stage}/output"
    printf '%q ' "${command[@]}" >"${stage}/COMMAND.txt"
    printf '\n' >>"${stage}/COMMAND.txt"
    (
        cd "${stage}"
        /usr/bin/time -v -o time.txt \
            /usr/bin/timeout --signal=TERM --kill-after=10s 600 \
            /usr/bin/taskset --cpu-list 0 \
            "${command[@]}" >console.log 2>&1
    ) || die "solver failed: ${label}"

    grep -Eq 'AMReX .* finalized' "${stage}/console.log" ||
        die "missing AMReX finalization: ${label}"
    if grep -Eqi "${fatal_pattern}" "${stage}/console.log"; then
        die "fatal or nonfinite signature: ${label}"
    fi

    local plotfile
    plotfile=$(find "${stage}/output" -mindepth 1 -maxdepth 1 -type d \
        -name 'plt*' -printf '%p\n' | sort -V | tail -n 1)
    [[ -n ${plotfile} && -f ${plotfile}/Header ]] ||
        die "missing plotfile: ${label}"
    printf '%s\n' "${plotfile}" >"${stage}/FINAL_PLOT.txt"
    (
        cd "${plotfile}"
        find . -type f -print0 | sort -z | xargs -0 sha256sum | sha256sum
    ) >"${stage}/FINAL_PLOT_TREE.sha256"
    date -u +'%Y-%m-%dT%H:%M:%SZ' >"${stage}/EXECUTION_PASS"
}

for nr in 16 32 64 128 256; do
    run_case "rhs_n${nr}" "${nr}" 0 1 1 1.0 1
done
for nr in 16 32 64 128 256; do
    run_case "pressure_rhs_n${nr}" "${nr}" 0 1 1 1.0 1 \
        prob.rho_amp=0 prob.ur_amp=0 prob.uz_amp=0 prob.uz0=0
done
for nr in 16 32 64 128 256; do
    run_case "evolve_n${nr}" "${nr}" 3 3 10000 2.0e-2 100000
done

date -u +'%Y-%m-%dT%H:%M:%SZ' >"${result_root}/COMPLETED_UTC"
printf '%s\n' \
    'purpose=R-Z full Euler time-dependent manufactured-solution matrix' \
    "execution_mode=${execution_mode}" \
    'rhs_case_count=10' \
    'evolution_case_count=5' \
    'analysis_required=1' >"${result_root}/MATRIX_EXECUTION_PASS"

echo '[RZ-MMS-FULL-MATRIX-EXECUTION-PASS] 15 cases completed'
