#!/usr/bin/env bash
set -euo pipefail

case_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
solver_root=$(cd "$case_root/../../.." && pwd)
cd "$case_root"

output_root=${SHU_AFD_TENO_ROOT:-"$case_root/results/y9000x_current_afd_hllc_teno5_20260729_r1"}
build_tag=${SHU_AFD_TENO_BUILD_TAG:-current_afd_teno_20260729_r1}
build_jobs=${CERISSE_BUILD_JOBS:-8}
mpi_ranks=${CERISSE_MPI_RANKS:-4}
executable="$case_root/main1d.gnu.MPI.shu_afd-hllc-teno5_${build_tag}.ex"
reference_root=${SHU_REFERENCE_ROOT:-"$case_root/results/y9000x_skew_jst_c2_1_5_c4_0_016_20260728_r3/N6400_reference"}

if [[ -e "$output_root" ]]; then
  echo "$output_root already exists; select a new SHU_AFD_TENO_ROOT" >&2
  exit 1
fi

mkdir -p "$output_root/source_snapshot/src_rhs" "$output_root/build"

cp GNUmakefile Make.package prob.h inputs.thesis analyze.py \
  run_current_afd_teno_check.sh "$output_root/source_snapshot/"
cp "$solver_root/src/rhs/Afd.h" \
  "$solver_root/src/rhs/AfdIBM.h" \
  "$solver_root/src/rhs/AfdIBMDrivers.h" \
  "$solver_root/src/rhs/Weno.h" \
  "$solver_root/src/rhs/Riemann.h" \
  "$solver_root/src/rhs/RHS.h" \
  "$output_root/source_snapshot/src_rhs/"
git -C "$solver_root" status --short \
  >"$output_root/source_snapshot/git_status_short.txt" 2>&1 || true
git -C "$solver_root" diff --binary \
  >"$output_root/source_snapshot/git_worktree.patch" 2>&1 || true

{
  echo "case=Shu-Osher"
  echo "method=AFD-HLLC-TENO5"
  echo "source_state=current_worktree"
  echo "grids=200 400 800"
  echo "stop_time=1.8"
  echo "cfl=0.30"
  echo "time_integrator=SSPRK(3,3)"
  echo "cns.afd_correction=1"
  echo "cns.afd_shock_llf=1"
  echo "cns.afd_smoothness_threshold=0.08"
  echo "cns.afd_teno_cutoff=1.0e-4"
  echo "git_commit=$(git -C "$solver_root" rev-parse HEAD)"
  echo "host=$(hostname)"
  echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} >"$output_root/manifest.txt"

make -j"$build_jobs" \
  SHU_EULER_SCHEME=afd-hllc-teno5 \
  SHU_BUILD_TAG="$build_tag" \
  USE_PELEPHYSICS=FALSE USE_EB=FALSE USE_GPIBM=FALSE \
  >"$output_root/build/build.log" 2>&1

sha256sum "$executable" >"$output_root/executable.sha256"
cp "$executable" "$output_root/"

for cells in 200 400 800; do
  case_dir="$output_root/N${cells}"
  mkdir -p "$case_dir"
  {
    printf '%q ' \
      mpirun --oversubscribe --bind-to none -np "$mpi_ranks" \
      "$executable" "$case_root/inputs.thesis" \
      "amr.n_cell=$cells" \
      "amr.plot_file=$case_dir/plt" \
      "cns.afd_correction=1" \
      "cns.afd_shock_llf=1" \
      "cns.afd_smoothness_threshold=0.08" \
      "cns.afd_teno_cutoff=1.0e-4"
    printf '\n'
  } >"$case_dir/run_command.txt"

  mpirun --oversubscribe --bind-to none -np "$mpi_ranks" \
    "$executable" "$case_root/inputs.thesis" \
    "amr.n_cell=$cells" \
    "amr.plot_file=$case_dir/plt" \
    "cns.afd_correction=1" \
    "cns.afd_shock_llf=1" \
    "cns.afd_smoothness_threshold=0.08" \
    "cns.afd_teno_cutoff=1.0e-4" \
    >"$case_dir/run.log" 2>&1
  echo "exit_status=0" >"$case_dir/run_status.txt"
done

MPLCONFIGDIR=/tmp/cerisse_shu_afd_teno_mpl \
  python3 analyze.py "$reference_root" \
    "$output_root/N200" "$output_root/N400" "$output_root/N800" \
    --csv "$output_root/metrics.csv" \
    --figure "$output_root/density_profiles.pdf" \
    >"$output_root/analysis.log" 2>&1

find "$output_root/source_snapshot" -type f -print0 \
  | sort -z | xargs -0 sha256sum \
  >"$output_root/source_snapshot/SHA256SUMS"
echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)" >>"$output_root/manifest.txt"
touch "$output_root/complete"
