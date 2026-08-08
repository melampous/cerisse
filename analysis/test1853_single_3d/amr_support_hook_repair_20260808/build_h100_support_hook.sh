#!/usr/bin/env bash
set -euo pipefail

tree=/shared/cerisse_local_0808_support_hook
case_dir=${tree}/exm/test1853_support_hook
exe=main3d.gnu.TPROF.MPI.CUDA.ex

test -f "${tree}/src/CNS.cpp"
test -f "${case_dir}/GNUmakefile"
test -f "${case_dir}/prob.h"
test -f "${case_dir}/inputs.run165"
grep -Fq 'dispatch=automatic' "${tree}/src/CNS.cpp"
if grep -Eq '^[[:space:]]*DEFINES.*CNS_IBM_PROBLEM_SUPPORT_CAP' \
    "${case_dir}/GNUmakefile"; then
  printf 'legacy support-cap macro must be absent from the acceptance build\n' >&2
  exit 2
fi

source /shared/gpu_kit/gpu_env.sh
cd "${case_dir}"
rm -f "${exe}"
make -j24 > build.log 2>&1
test -x "${exe}"
test "$(cuobjdump "${exe}" | grep -o 'sm_[0-9]*' | sort -u)" = sm_90
if grep -Fq '/shared/cerisse_local_0808_thermodynamic_prolongation' build.log; then
  printf 'stale source tree detected in build log\n' >&2
  exit 3
fi

{
  printf 'built_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'host=%s\n' "$(hostname)"
  printf 'cuda_release=%s\n' "$(nvcc --version | awk '/release/ {gsub(/,/,"",$6); print $6}')"
  printf 'cuda_arch=90\n'
  printf 'amr_solver=%s\n' "${tree}"
  printf 'method=pure_shared_gp_full_cartesian\n'
  printf 'cut_control=disabled\n'
  printf 'support_hook_dispatch=automatic\n'
  printf 'legacy_support_cap_macro=absent\n'
  sha256sum "${tree}/src/CNS.cpp" GNUmakefile prob.h inputs.run165 "${exe}"
} > support_hook_build_manifest.txt

cat support_hook_build_manifest.txt
