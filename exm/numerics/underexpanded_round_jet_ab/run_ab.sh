#!/usr/bin/env bash
set -euo pipefail

case_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(git -C "$case_root" rev-parse --show-toplevel)
canonical_case="$repo_root/exm/underexpanded_jet/3d/npr3_lev4_ns"
result_root=${JET_AB_RESULTS_ROOT:-"$case_root/results/sensor_ab_20260729_r1"}
build_jobs=${CERISSE_BUILD_JOBS:-4}
executable="$case_root/main3d.gnu.TPROF.MPI.CUDA.jet_afd-hllc-wenoz5_sensor-ab.ex"

die()
{
    echo "[UNDEREXPANDED-JET-AB-FAIL] $*" >&2
    exit 2
}

if [[ $result_root != /* ]]; then
    result_root="$case_root/$result_root"
fi
[[ ! -e $result_root ]] || die "refusing to reuse result root: $result_root"
[[ -f "$canonical_case/prob.h" ]] || die "missing canonical prob.h"

mkdir -p "$result_root/runs"
date -u +'%Y-%m-%dT%H:%M:%SZ' >"$result_root/STARTED_UTC"
git -C "$repo_root" rev-parse HEAD >"$result_root/GIT_HEAD"
git -C "$repo_root" status --short >"$result_root/GIT_STATUS.txt"
git -C "$repo_root" diff -- \
    src/rhs/Afd.h src/rhs/AfdIBM.h src/rhs/AfdIBMDrivers.h \
    src/rhs/Riemann.h src/rhs/Weno.h \
    >"$result_root/RELEVANT_SOURCE_DIFF.patch"

{
    echo "case=reduced_3d_cartesian_npr3_viscous_round_jet"
    echo "canonical_problem=$canonical_case/prob.h"
    echo "grid=96 64 64"
    echo "domain=0.0 0.06 -0.02 0.02 -0.02 0.02"
    echo "max_level=0"
    echo "cfl=0.30"
    echo "stop_time=7.0e-5"
    echo "time_integrator=SSPRK33"
    echo "scheme=AFD-HLLC-WENO-Z5"
    echo "afd_smoothness_threshold=0.08"
    echo "afd_shock_pressure_jump=0.01"
    echo "afd_shock_compression=0.001"
    echo "paired_variant_only_argument=cns.afd_shock_llf"
    echo "use_gpibm=FALSE"
    echo "use_eb=FALSE"
    echo "use_pelephysics=FALSE"
} >"$result_root/manifest.txt"

make -C "$case_root" -j"$build_jobs" \
    USE_GPIBM=FALSE USE_EB=FALSE USE_PELEPHYSICS=FALSE \
    >"$result_root/build.log" 2>&1 \
    || die "build failed; inspect $result_root/build.log"
[[ -x $executable && ! -L $executable ]] || die "missing executable: $executable"

sha256sum \
    "$executable" \
    "$case_root/GNUmakefile" \
    "$case_root/prob.h" \
    "$case_root/inputs" \
    "$case_root/run_ab.sh" \
    "$case_root/analyze_ab.py" \
    "$canonical_case/prob.h" \
    "$repo_root/src/rhs/Afd.h" \
    "$repo_root/src/rhs/AfdIBM.h" \
    "$repo_root/src/rhs/AfdIBMDrivers.h" \
    "$repo_root/src/rhs/Riemann.h" \
    "$repo_root/src/rhs/Weno.h" \
    >"$result_root/RUN_INPUTS.sha256"
sha256sum "$executable" >"$result_root/executable.sha256"

export OMP_NUM_THREADS=1
export LC_ALL=C
export MPLCONFIGDIR=/tmp/underexpanded_round_jet_ab_mpl

fatal_pattern='MPI_ABORT|amrex::Abort|SIG(SEGV|FPE|ABRT|BUS)|Segmentation fault|Floating point exception|non-finite|nonfinite|(^|[^[:alnum:]_])(nan|[+-]?inf(inity)?)([^[:alnum:]_]|$)'

run_variant()
{
    local mode=$1
    local stage="$result_root/runs/shock_llf_$mode"
    local -a command=(
        mpirun --oversubscribe -np 1
        "$executable" "$case_root/inputs"
        "cns.afd_shock_llf=$mode"
    )

    mkdir -p "$stage/plot"
    printf 'working_directory=%q\n' "$stage" >"$stage/COMMAND.txt"
    printf 'command=' >>"$stage/COMMAND.txt"
    printf '%q ' "${command[@]}" >>"$stage/COMMAND.txt"
    printf '\n' >>"$stage/COMMAND.txt"
    (
        cd "$stage"
        /usr/bin/time -v -o time.txt \
            /usr/bin/timeout --signal=TERM --kill-after=20s 1800 \
            "${command[@]}" >console.log 2>&1
    ) || die "solver failed for afd_shock_llf=$mode"

    grep -Eq 'AMReX .* finalized' "$stage/console.log" \
        || die "missing AMReX finalization for afd_shock_llf=$mode"
    if grep -Eqi "$fatal_pattern" "$stage/console.log"; then
        die "fatal or nonfinite signature for afd_shock_llf=$mode"
    fi
    find "$stage/plot" -mindepth 1 -maxdepth 1 -type d \
        -name 'plt*' -print -quit | grep -q . \
        || die "missing plotfile for afd_shock_llf=$mode"
    date -u +'%Y-%m-%dT%H:%M:%SZ' >"$stage/EXECUTION_PASS"
}

run_variant 0
run_variant 1

python3 "$case_root/analyze_ab.py" "$result_root" \
    >"$result_root/analysis.log" \
    || die "analysis failed; inspect $result_root/analysis.log"

date -u +'%Y-%m-%dT%H:%M:%SZ' >"$result_root/COMPLETED_UTC"
echo "[UNDEREXPANDED-JET-AB-PASS] paired calculations and analysis completed"
