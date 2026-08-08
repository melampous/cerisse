#!/usr/bin/env bash
set -euo pipefail

tree=/shared/cerisse_local_0808_admissible_prolongation
case_dir="${tree}/exm/test1853_admissible_prolongation"
reference_case=/shared/cerisse_local_0808_regridprov/exm/test1853_regridprov
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

test -f "${tree}/src/set/EulerAdmissibleProlongation.H"
test -f "${tree}/src/set/CNS_setup.cpp"
test -f "${tree}/src/CNS.cpp"
test -f "${reference_case}/prob.h"

source /shared/gpu_kit/gpu_env.sh
nvcc_release=$(nvcc --version | sed -n 's/.*release \([0-9][0-9.]*\).*/\1/p' | tail -1)
if [[ "${nvcc_release}" != "12.6" ]]; then
  printf 'ERROR: CUDA 12.6 is required; found %s\n' "${nvcc_release:-unknown}" >&2
  exit 2
fi

mkdir -p "${case_dir}"
cp "${script_dir}/h100_repair_case/GNUmakefile" "${case_dir}/GNUmakefile"
cp "${reference_case}/prob.h" "${case_dir}/prob.h"
cd "${case_dir}"

rm -rf tmp_build_dir
rm -f main3d.gnu.TPROF.MPI.CUDA.ex
make -j24 > build.log 2>&1
test -x main3d.gnu.TPROF.MPI.CUDA.ex

{
  printf 'built_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'host=%s\n' "$(hostname)"
  printf 'cuda_release=%s\n' "${nvcc_release}"
  printf 'cuda_arch=90\n'
  printf 'method=pure_shared_gp_full_cartesian\n'
  printf 'cut_control=disabled\n'
  sha256sum \
    "${tree}/src/set/EulerAdmissibleProlongation.H" \
    "${tree}/src/set/CNS_setup.cpp" \
    "${tree}/src/CNS.cpp" \
    main3d.gnu.TPROF.MPI.CUDA.ex
} > repair_build_manifest.txt

cat repair_build_manifest.txt
