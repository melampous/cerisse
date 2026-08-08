#!/usr/bin/env bash
set -euo pipefail

binary=$(realpath "${1:?usage: run_matrix.sh BINARY INPUTS RESULT_ROOT}")
inputs=$(realpath "${2:?usage: run_matrix.sh BINARY INPUTS RESULT_ROOT}")
result_root=${3:?usage: run_matrix.sh BINARY INPUTS RESULT_ROOT}
mpi_root=${MPI_ROOT:-/home/melampous/cerisse_rz_bench_20260720/ompi_user/root/usr}
run_serial=${RUN_SERIAL:-0}

fatal_pattern='MPI_ABORT|amrex::Abort|SIG(SEGV|FPE|ABRT|BUS)|Segmentation fault|Floating point exception|non-finite|nonfinite|(^|[^[:alnum:]_])(nan|[+-]?inf(inity)?)([^[:alnum:]_]|$)'

die()
{
    echo "[RZ-EULER-AXIS-MATRIX-FAIL] $*" >&2
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
    local mode=$2
    local nr=$3
    local annular=$4
    local point_flux=$5
    local stage=${result_root}/runs/${label}
    local -a command=(
        "${launcher[@]}" "${binary}" "${inputs}"
        "prob.mode=${mode}"
        "amr.n_cell=${nr} ${nr}"
        "cns.rz_euler_annular_average=${annular}"
        "cns.rz_euler_point_flux=${point_flux}"
        amr.plot_file=output/plt
        amr.checkpoint_files_output=0
    )

    mkdir -p "${stage}/output"
    printf '%q ' "${command[@]}" >"${stage}/COMMAND.txt"
    printf '\n' >>"${stage}/COMMAND.txt"
    (
        cd "${stage}"
        /usr/bin/time -v -o time.txt \
            /usr/bin/timeout --signal=TERM --kill-after=10s 180 \
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

for nr in 16 32 64 128; do
    run_case "m0_legacy_n${nr}" 0 "${nr}" 0 0
done
run_case m0_annular_n16 0 16 1 0
run_case m0_corrected_n16 0 16 1 1
for nr in 16 32 64 128; do
    run_case "m1_corrected_n${nr}" 1 "${nr}" 1 1
done

date -u +'%Y-%m-%dT%H:%M:%SZ' >"${result_root}/COMPLETED_UTC"
printf '%s\n' \
    'purpose=R-Z Euler-axis legacy/correction convergence matrix' \
    "execution_mode=${execution_mode}" \
    'case_count=10' \
    'analysis_required=1' >"${result_root}/MATRIX_EXECUTION_PASS"

echo '[RZ-EULER-AXIS-MATRIX-EXECUTION-PASS] 10 cases completed'
