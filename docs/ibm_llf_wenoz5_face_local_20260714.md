# Hybrid face-local LLF-WENO-Z5 IBM implementation and validation

Date: 2026-07-14

## 1. Decision

A generic face-local Rusanov/LLF crossing-face path is implemented for the
WENO-Z5/TENO5 flux class.  Smooth-wall MMS subsequently proved that applying
this reflected flux at every offset Cartesian crossing face is inconsistent.
The production path is therefore hybrid.  For a build with
`face_local_euler_slip=true`, its current runtime defaults are

```text
cns.llf_ibm_face_local_crossing = 1
cns.llf_ibm_face_local_activation_threshold = 0.02
cns.llf_ibm_fluid_shell         = 0
```

Smooth crossing faces retain the dense legacy GP/WENO flux.  The sparse
face-local LLF pass overwrites a face only when its fluid-side normalized jump
exceeds the threshold or when it belongs to an annotated 2-D sharp-feature
band.  A zero threshold reproduces the forced-direct ablation.  The
real-fluid-only near-wall candidate mask remains available for experiments
but is default-off because the compression-corner convergence test found a
nonlinear crossing/shell mismatch.

This result establishes a usable hybrid LLF-WENO-Z5 closure for the tested
2-D, Cartesian, stationary, single-species Euler-slip IBM shock problem.  The
smooth-wall study establishes strong global integral-norm convergence but not
uniform second order in the first fluid layer.  NS no-slip compatibility,
moving walls, RZ, multispecies support, AMR/MPI equivalence, and LLF-TENO5
production status remain separate gates.

## 2. Numerical closure

For each Cartesian fluid/solid crossing face, initialization caches one
first-hit surface element, its wall normal, visible fluid support points, and
WLS evaluation weights.  The dense pass first computes the established
GP/WENO crossing flux.  At every RHS evaluation, the sparse pass then:

1. WLS evaluates a fluid extension `W_f` at the Cartesian face centre.
2. The same support polynomial is evaluated at the point reflected through
   the owning-element tangent plane, giving `W_r`.
3. Density and pressure are extended evenly.  Velocity is reflected about the
   wall normal `n`:

   ```text
   u_g = u_r - 2 (u_r . n) n .
   ```

4. Both complete primitive states are converted back to thermodynamically
   consistent conservative states.  The Cartesian-direction physical flux is
   then evaluated on both sides.
5. Computes a normalized fluid-side jump.  In a smooth, non-feature region it
   retains the dense legacy flux.  Otherwise one Rusanov flux is formed:

   ```text
   F_LLF = 0.5 (F_L + F_R) - 0.5 alpha (U_R - U_L),
   alpha = max(|u_d,L| + c_L, |u_d,R| + c_R).
   ```

The STL/polygon normal `n` is used only for the wall reflection.  The Riemann
flux direction `d` remains the target Cartesian face normal.  Mixing these two
directions would be mathematically incorrect.

Only the explicit zero-threshold ablation makes the dense WENO kernel leave a
zero placeholder on every crossing face.  In hybrid mode the sparse record
pass overwrites only activated target faces, exactly once.  Every fluid update
sharing an activated numerical face therefore sees the same conservative
flux.  No direction-local primitive value is written into the shared
primitive MultiFab.

If WLS metadata, weights, positivity, internal energy, or EOS consistency are
invalid, both extension targets fall back together to the adjacent real-fluid
state before velocity reflection.  Missing first-hit metadata or a nonfinite
final flux writes NaN deliberately, causing the existing solver checks to fail
fast instead of silently using a stale or zero wall flux.

## 3. Code organization

- `src/rhs/IbmFluxUtils.h` owns primitive-state validation, stationary-slip
  reflection, WLS pair construction, primitive-to-conservative conversion,
  and the Rusanov flux.
- `src/rhs/Weno.h` owns the WENO/TENO dense pass, runtime switches, sparse
  crossing-face overwrite, and optional real-fluid-shell experiment.
- `src/tim/compute_rhs.cpp` dispatches the boundary-face view through a
  compile-time detected overload; unsupported flux types abort explicitly.
- `src/ibm/ibm_containers.h` and `src/ibm/ibm_solver.h` distinguish LLF from
  HLLC/HLLE in runtime status and report fallback, duplicate-flux, jump, and
  shell counters.
- The inclined-wall, no-solid oblique-shock, and compression-corner cases now
  expose LLF-WENO-Z5 in their build/validation matrices.

Runtime diagnostics are disabled unless `ib.face_local_diag=1`.  The
initialization WLS/visibility cache is reused by every RK stage; BVH visibility
queries are not placed in the RHS hot path.

## 4. Constant-state and memory-isolation gates

The 128 by 128 inclined-wall forced-direct ablation uses a uniform Mach-2
state exactly tangent to a 20-degree immersed wall.  The exact conservative
RHS is zero in every fluid cell.  With face-local LLF forced on:

| Gate | Result |
|---|---:|
| Maximum normalized RHS Linf in wall layers 1--3 | `4.225344e-14` |
| Bulk RHS | exactly zero |
| Active crossing faces | 175 |
| WLS / first-order fallback / fatal | 175 / 0 / 0 |
| Duplicate-face flux inconsistencies | 0 |

Repeating the run with the interior-solid primitive storage poisoned produced
zero absolute and relative plotfile differences in every field.  The new
crossing flux therefore does not depend on the shared solid primitive state.

In the production hybrid smooth-circle MMS, all 48/80/104/152 crossing faces
at `N=32/48/64/96` selected `smooth-legacy`; none selected direct LLF and no
fatal or duplicate-flux inconsistency occurred.  In the `N=80` compression
regression, the final state selected 25 direct LLF faces in the sharp-feature
band and retained 78 smooth legacy faces.  This confirms that selection is
localized rather than globally forced.

## 5. No-solid exact oblique-shock baseline

The no-solid problem uses exact conservative cut-cell averages at `t=0`, exact
time-dependent physical boundary states, a Galilean shock speed of 0.1, RK3,
CFL 0.30, and a fixed large domain `[-10,20] x [-18,18]`.  At `t=2` the
diagnostic region is causally isolated from the boundary/shock intersections.

| h | Grid | downstream p L1 | beta error [deg] | strip-excluded p L1 | overshoot / p1 |
|---:|---:|---:|---:|---:|---:|
| 0.250000 | 120x144 | 0.1322% | +0.16270 | 0.1617% | 2.919% |
| 0.125000 | 240x288 | 0.08238% | +0.03705 | 0.09322% | 3.434% |
| 0.083333 | 360x432 | 0.06467% | -0.00565 | 0.06581% | 2.587% |

At the finest common spacing, the prior LLF-TENO5 no-solid overshoot was
9.672% of upstream pressure; LLF-WENO-Z5 gives 2.587%.  The bulk LLF-WENO-Z5
operator is therefore not the source of the approximately 0.8-p1 IBM corner
overshoot seen with the legacy crossing closure.

These shock norms are not expected to show fifth-order convergence.  A
solution containing a discontinuity is globally first order in L1 under mesh
refinement even when the smooth bulk operator is formally fifth order.

## 6. Forced-direct IBM compression-corner refinement

The historical forced-direct IBM test is the stationary Mach-2, 20-degree
compression corner, integrated to `t=14` with the same fixed polygon,
`iorder=eorder=eorder_surf=2`, `alpha=1`, and one ghost layer.  This ablation
uses face-local crossing LLF on every crossing face and the established
GP-masked WENO-Z5 fluid shell.

| Nx | p2/p1 error | plateau p L1 | flow-angle error [deg] | beta error [deg] | strip-excluded p L1 | overshoot / p1 |
|---:|---:|---:|---:|---:|---:|---:|
| 40  | -1.636% | 1.654% | -0.1629 | -0.9905 | 0.3638% | 4.574% |
| 80  | -0.701% | 0.7908% | -0.0716 | -0.3667 | 0.3725% | 5.149% |
| 120 | -0.323% | 0.5330% | -0.0263 | -0.1710 | 0.3676% | 3.634% |

The plateau-pressure pairwise orders are 1.065 and 0.973.  Approximately first
order is the correct expectation for a shock-containing solution; this test
validates shock location, jump, wall turning, robustness, and removal of the
legacy corner overshoot, not formal smooth-wall order.

At Nx=120, face-local AFD-HLLC-WENO-Z5 gives plateau L1 0.5557%, p2 error
-0.301%, beta error -0.1032 degrees, and overshoot 3.672% of p1.  The new
LLF-WENO-Z5 values are 0.5330%, -0.323%, -0.1710 degrees, and 3.634%.  Within
this forced-direct sharp-corner validation, the LLF path is therefore
comparable to the current best AFD-HLLC-WENO-Z5 pressure result, with a
modestly larger shock-angle error.

After hybridization, the `N=80` regression gives p2 error `-0.841%`, plateau
pressure L1 `0.928%`, beta error `-0.3658 degrees`, and overshoot `3.943%` of
upstream pressure.  Exactly 25 of 103 crossing faces are direct LLF because of
the sharp-feature band.  This retains the corner improvement without applying
the same closure to smooth offset walls.

## 7. N=80 ablation and rejected shell coupling

| Crossing closure | Near-wall WENO candidate policy | plateau p L1 | beta error [deg] | overshoot / p1 |
|---|---|---:|---:|---:|
| legacy GP | legacy GP-masked | 0.9308% | +0.6902 | 80.35% |
| legacy GP | real-fluid-only experiment | 0.8795% | +0.5093 | 74.45% |
| forced-direct face-local LLF | legacy GP-masked | 0.7908% | -0.3667 | 5.149% |
| forced-direct face-local LLF | real-fluid-only experiment | 2.397% | -0.9270 | 4.531% |

At this sharp corner, the face-local crossing closure is the dominant
successful change: for the same LLF-WENO-Z5 bulk scheme it reduces overshoot
from 80.35% to 5.149% and also improves the plateau norm.  The real-fluid-only
shell suppresses the pointwise overshoot when combined with face-local
crossing but biases the post-shock pressure and destroys monotone grid
improvement.  It is therefore not part of the default path.

This is evidence against simply removing every GP value from every near-wall
WENO stencil.  The crossing flux and the neighbouring fluid-fluid closure must
be designed as a matched discrete boundary operator, then tested by smooth
MMS, rather than selected independently from local stencil purity.

## 8. Build and runtime checks

- GCC 2-D GPIBM LLF-WENO-Z5 build: pass.
- NVCC 2-D GPIBM LLF-WENO-Z5, `sm_89`, with
  `--Werror cross-execution-space-call`: pass.
- NVCC no-solid LLF-WENO-Z5 build and all three exact-shock runs: pass.
- CPU/GPU five-step active-path smoke: matching diagnostic counts and values.
- Forced-direct IBM compression runs at Nx=40/80/120: finite through `t=14`.
- Hybrid sharp-feature compression regression at Nx=80: finite through `t=14`.
- Hybrid smooth-circle MMS at Nx=32/48/64/96: CUDA fixed-time and t=0 runs
  complete with zero direct LLF activation, as intended.
- `rho`, `p`, and internal-energy failures: none in the reported runs.

## 9. Current limits and next gate

This work should not be described as a general second-order face-local IBM.
The completed four-grid smooth curved-wall sequence is documented in
`docs/ibm_llf_wenoz5_curved_mms_20260714.md`.  At fixed time, global L1 orders
are 2.24--2.56 and global L2 orders are 1.82--2.41.  However, first-layer L2
orders are only 1.06--1.66, and the exact-state t=0 RHS is non-monotone in that
layer.  A shock refinement sequence and good global L1 cannot replace this
local consistency gate.

Further required work is:

1. derive a wall-location-consistent crossing/shell matched one-sided WENO
   closure and accept it only if the first-layer t=0 RHS and fixed-time MMS
   improve without regressing the compression corner;
2. repeat CPU/GPU/MPI partition equivalence for this LLF path;
3. add a conservative cell-update positivity limiter for more severe states;
4. separately design NS no-slip viscous/thermal face-local coupling; and
5. validate LLF-TENO5 independently instead of inferring it from the common
   class template.

The defensible status is therefore: the legacy LLF-WENO-Z5 sharp-corner defect
has been repaired for the tested fixed 2-D Euler-slip shock configuration, and
the large pressure overshoot is removed.  The hybrid smooth-wall path has
strong global integral accuracy but is not uniformly second order at the
first fluid layer.  The broader physics matrix remains a separate acceptance
gate.
