#!/usr/bin/env bash
set -euo pipefail

tree=/shared/cerisse_local_0808_support_hook
source_case=${tree}/exm/test1853_support_hook
case_dir=${tree}/exm/test1853_support_hook_nojet
config_dir=/shared/cerisse_test1853_3d_20260807/analysis/h100_grid_redesign_20260808
exe=main3d.gnu.TPROF.MPI.CUDA.ex

test -f "${tree}/src/CNS.cpp"
test -f "${source_case}/prob.h"
test -f "${config_dir}/GNUmakefile.nojet"
test -f "${config_dir}/inputs.fresh_256x192x192_h100"
grep -Fq 'dispatch=automatic' "${tree}/src/CNS.cpp"

if [[ -e "${case_dir}" ]]; then
  printf 'refusing to reuse no-jet build directory: %s\n' "${case_dir}" >&2
  exit 2
fi

mkdir -p "${case_dir}"
cp "${config_dir}/GNUmakefile.nojet" "${case_dir}/GNUmakefile"
cp "${source_case}/prob.h" "${case_dir}/prob.h"
cp "${config_dir}/inputs.fresh_256x192x192_h100" "${case_dir}/inputs.run165"

source /shared/gpu_kit/gpu_env.sh
cd "${case_dir}"
make -j24 > build.log 2>&1
test -x "${exe}"
test "$(cuobjdump "${exe}" | grep -o 'sm_[0-9]*' | sort -u)" = sm_90
grep -Fq -- '-DTEST1853_JET_ENABLED=0' build.log
if grep -Eq -- '-DTEST1853_JET_ENABLED=1|CNS_IBM_PROBLEM_SUPPORT_CAP' build.log; then
  printf 'wrong jet mode or legacy support-cap macro found in build log\n' >&2
  exit 3
fi

{
  printf 'built_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'host=%s\n' "$(hostname)"
  printf 'cuda_release=%s\n' \
    "$(nvcc --version | awk '/release/ {gsub(/,/,"",$6); print $6}')"
  printf 'cuda_arch=90\n'
  printf 'amr_solver=%s\n' "${tree}"
  printf 'jet_enabled_compile_time=0\n'
  printf 'method=pure_shared_gp_full_cartesian\n'
  printf 'cut_control=disabled\n'
  printf 'support_hook_dispatch=automatic\n'
  printf 'legacy_support_cap_macro=absent\n'
  sha256sum "${tree}/src/CNS.cpp" GNUmakefile prob.h inputs.run165 "${exe}"
} > support_hook_build_manifest.txt

cat support_hook_build_manifest.txt

