#!/usr/bin/env bash
set -euo pipefail

case_root=/shared/cerisse_test1853_3d_20260807/runs/test1853_single_h100/run165_llfdiag_3way_20260808
build_case=/shared/cerisse_local_0808_regridprov/exm/test1853_regridprov
run_dir="${case_root}/testD_regrid_provenance"
restart="${case_root}/testB_normal_regrid/chk03059"
geometry=/shared/cerisse_test1853_3d_20260807/runs/test1853_single_h100/run165_chain_20260807T115518Z/geometry/test1853_run165_single_master_SI.stl

mkdir -p "${run_dir}"
mkdir -p "${run_dir}/output/run165"
cp "${build_case}/main3d.gnu.TPROF.MPI.CUDA.ex" "${run_dir}/"
cp "${case_root}/inputs.run165" "${run_dir}/"
cp "${case_root}/prob.h" "${run_dir}/"

cd "${run_dir}"
export CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7

set +e
mpirun -np 8 ./main3d.gnu.TPROF.MPI.CUDA.ex inputs.run165 \
  amr.restart="${restart}" \
  ib.filename="${geometry}" \
  amr.regrid_on_restart=1 \
  max_step=3060 \
  stop_time=1.0 \
  dt_max=1.0e-14 \
  cns.stage_positivity=0 \
  cns.llf_fallback_mask_diagnostics=1 \
  cns.regrid_prolongation_diagnostics=1 \
  amr.plot_files_output=0 \
  amr.checkpoint_files_output=0 \
  cns.verbose=0 \
  > run.log 2>&1
status=$?
set -e

printf '%s\n' "${status}" > exit_status.txt
exit "${status}"
