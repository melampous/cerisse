#!/usr/bin/env bash
set -euo pipefail

binary=$(realpath "${1:?usage: run_ab.sh BINARY INPUTS RESULT_ROOT}")
inputs=$(realpath "${2:?usage: run_ab.sh BINARY INPUTS RESULT_ROOT}")
result_root=${3:?usage: run_ab.sh BINARY INPUTS RESULT_ROOT}
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(git -C "${script_dir}" rev-parse --show-toplevel)

die()
{
    echo "[AFD-GRID-SHOCK-FAIL] $*" >&2
    exit 2
}

[[ -x ${binary} && ! -L ${binary} ]] || die "missing executable: ${binary}"
[[ -f ${inputs} && ! -L ${inputs} ]] || die "missing inputs: ${inputs}"
[[ ! -e ${result_root} ]] || die "refusing to reuse result root: ${result_root}"

mkdir -p "${result_root}/runs"
{
    printf 'case=Cartesian stationary Mach-10 grid-aligned shock\n'
    printf 'study=same-source same-binary cns.afd_shock_llf A/B\n'
    printf 'method=AFD-HLLC-WENO-Z5\n'
    printf 'afd_shock_llf_modes=0 1\n'
    printf 'afd_smoothness_threshold=0.08\n'
    printf 'afd_shock_pressure_jump=0.01\n'
    printf 'afd_shock_compression=0.001\n'
    printf 'cfl=0.25\n'
    printf 'stop_time=0.05\n'
    printf 'time_integrator=SSPRK(3,3)\n'
    printf 'hostname=%s\n' "$(hostname)"
    printf 'created_utc='
    date -u +'%Y-%m-%dT%H:%M:%SZ'
} >"${result_root}/RUN_MANIFEST.txt"
sha256sum \
    "${binary}" \
    "${inputs}" \
    "${script_dir}/GNUmakefile" \
    "${script_dir}/prob.h" \
    "${script_dir}/run_ab.sh" \
    "${script_dir}/analyze.py" \
    "${repo_root}/src/rhs/Afd.h" \
    "${repo_root}/src/rhs/AfdIBM.h" \
    "${repo_root}/src/rhs/Riemann.h" \
    "${repo_root}/src/rhs/Weno.h" >"${result_root}/RUN_INPUTS.sha256"
sha256sum \
    "${repo_root}/src/rhs/Afd.h" \
    "${repo_root}/src/rhs/AfdIBM.h" \
    "${repo_root}/src/rhs/Riemann.h" \
    >"${result_root}/FROZEN_CORE_SOURCE.sha256"
git -C "${repo_root}" rev-parse HEAD >"${result_root}/GIT_HEAD"
git -C "${repo_root}" status --short >"${result_root}/GIT_STATUS.txt"
git -C "${repo_root}" diff -- \
    src/rhs/Afd.h src/rhs/Riemann.h src/rhs/Weno.h \
    >"${result_root}/RELEVANT_SOURCE_DIFF.patch"
date -u +'%Y-%m-%dT%H:%M:%SZ' >"${result_root}/STARTED_UTC"

export OMP_NUM_THREADS=1
export LC_ALL=C
export MPLCONFIGDIR=/tmp/afd_gridshock_mpl

fatal_pattern='MPI_ABORT|amrex::Abort|SIG(SEGV|FPE|ABRT|BUS)|Segmentation fault|Floating point exception|non-finite|nonfinite|(^|[^[:alnum:]_])(nan|[+-]?inf(inity)?)([^[:alnum:]_]|$)'

run_case()
{
    local label=$1
    local nx=$2
    local ny=$3
    local llf=$4
    local stop=$5
    local plot_int=$6
    local seed=$7
    local stage=${result_root}/runs/${label}
    local -a command=(
        mpirun --oversubscribe -np 2
        "${binary}" "${inputs}"
        "amr.n_cell=${nx} ${ny}"
        "stop_time=${stop}"
        max_step=100000
        "amr.plot_int=${plot_int}"
        amr.plot_file=plot/plt
        amr.checkpoint_files_output=0
        cns.nstep_screen_output=100
        cns.verbose=0
        "cns.afd_shock_llf=${llf}"
        cns.afd_shock_pressure_jump=0.01
        cns.afd_shock_compression=0.001
        "prob.shock_odd_even_cells=${seed}"
    )

    mkdir -p "${stage}/plot"
    printf '%q ' "${command[@]}" >"${stage}/COMMAND.txt"
    printf '\n' >>"${stage}/COMMAND.txt"
    (
        cd "${stage}"
        /usr/bin/time -v -o time.txt \
            /usr/bin/timeout --signal=TERM --kill-after=10s 900 \
            "${command[@]}" >console.log 2>&1
    ) || die "solver failed: ${label}"
    grep -Eq 'AMReX .* finalized' "${stage}/console.log" ||
        die "missing AMReX finalization: ${label}"
    if grep -Eqi "${fatal_pattern}" "${stage}/console.log"; then
        die "fatal or nonfinite signature: ${label}"
    fi

    mapfile -t plots < <(
        find "${stage}/plot" -mindepth 1 -maxdepth 1 -type d \
            -name 'plt*' -printf '%p\n' | sort -V
    )
    [[ ${#plots[@]} -ge 2 ]] || die "missing initial/final plots: ${label}"
    python3 "${script_dir}/analyze.py" \
        --gamma 1.4 \
        --smoothness-threshold 0.08 \
        --pressure-jump-threshold 0.01 \
        --compression-threshold 0.001 \
        --json "${stage}/metrics_history.json" \
        --csv "${stage}/metrics_history.csv" \
        "${plots[@]}" >"${stage}/analysis.log"
    date -u +'%Y-%m-%dT%H:%M:%SZ' >"${stage}/EXECUTION_PASS"
}

# The 200 by 20 grid is the smoke and coarse stress test.  The 400 by 40
# calculation repeats the same physical problem with square cells at half dx.
run_case n200_planar_llf0 200 20 0 5.0e-2 100 0.0
run_case n200_planar_llf1 200 20 1 5.0e-2 100 0.0
run_case n200_seeded_llf0 200 20 0 5.0e-2 100 0.01
run_case n200_seeded_llf1 200 20 1 5.0e-2 100 0.01
run_case n400_seeded_llf0 400 40 0 5.0e-2 200 0.01
run_case n400_seeded_llf1 400 40 1 5.0e-2 200 0.01

python3 - "${result_root}" <<'PY'
import csv
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
rows = []
for path in sorted(root.glob("runs/*/metrics_history.json")):
    history = json.loads(path.read_text())
    first = history[0]
    final = history[-1]
    rows.append({
        "case": path.parent.name,
        "time_initial": first["time"],
        "time_final": final["time"],
        "shock_p2p_initial_cells": first["shock_peak_to_peak_cells"],
        "shock_p2p_final_cells": final["shock_peak_to_peak_cells"],
        "shock_odd_even_initial_cells":
            first["shock_odd_even_abs_cells"],
        "shock_odd_even_final_cells":
            final["shock_odd_even_abs_cells"],
        "shock_odd_even_growth":
            final["shock_odd_even_abs_cells"] /
            max(first["shock_odd_even_abs_cells"], 1.0e-300),
        "uy_absmax_final": final["uy_absmax_global"],
        "density_min_final": final["density_min_global"],
        "pressure_min_final": final["pressure_min_global"],
        "shock_mean_offset_final_cells":
            final["shock_mean_offset_cells"],
        "sensor_x_fallback_initial_faces":
            first["sensor_x_fallback_faces"],
        "sensor_x_fallback_final_faces":
            final["sensor_x_fallback_faces"],
        "sensor_y_fallback_initial_faces":
            first["sensor_y_fallback_faces"],
        "sensor_y_fallback_final_faces":
            final["sensor_y_fallback_faces"],
        "sensor_y_fallback_initial_rows":
            first["sensor_y_fallback_rows"],
        "sensor_y_fallback_final_rows":
            final["sensor_y_fallback_rows"],
    })

    expected_initial_y_faces = 0
    if "_seeded_" in path.parent.name:
        expected_initial_y_faces = int(first["ny"])
    if first["sensor_y_fallback_faces"] != expected_initial_y_faces:
        raise SystemExit(
            f"{path.parent.name}: expected "
            f"{expected_initial_y_faces} initial y fallback faces, got "
            f"{first['sensor_y_fallback_faces']}"
        )
(root / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")
with (root / "summary.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(json.dumps(rows, indent=2))
PY

date -u +'%Y-%m-%dT%H:%M:%SZ' >"${result_root}/COMPLETED_UTC"
echo "[AFD-GRID-SHOCK-PASS] six A/B cases completed"
