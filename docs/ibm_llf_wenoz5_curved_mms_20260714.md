# LLF-WENO-Z5 smooth curved-wall IBM MMS

Date: 2026-07-14

> **Superseded result.** This document records the earlier `eorder=2` and
> pre-interpolation-halo study. The final `eorder=3`, partition-independent
> audit is in
> `docs/ibm_llf_wenoz5_curved_mms_e3_final_20260714.md`. In particular, the
> final recommended path keeps the experimental wall-matched crossing and
> shifted same-side shell disabled and demonstrates second-order-or-better
> first-layer accuracy.

## 1. Executive result

The four-grid `N=32/48/64/96` study gives a mixed but useful result.

- At fixed time `t=0.002`, the hybrid LLF-WENO-Z5 IBM is better than second
  order in global L1 for every primitive variable.  Global L2 is order
  2.30--2.41 for density, pressure, and temperature, and about 1.82 for both
  velocity components.
- The first fluid layer is not uniformly second order.  Its four-grid fitted
  L2 orders are 1.35 for density, 1.06--1.07 for velocity, 1.23 for pressure,
  and 1.66 for temperature.  Linf orders are only 0.83--1.42.
- The exact-state `t=0` semi-discrete RHS is non-monotone in the first layer.
  Its fitted L2 orders are negative for mass and energy and about 0.57--0.59
  for momentum.  The largest errors are sparse and move around the circle as
  the Cartesian grid is refined.
- Forcing the reflected face-local LLF closure on every smooth crossing face
  is inconsistent with an immersed wall offset from that Cartesian face.  It
  loses convergence and its error grows catastrophically at `N=64`.  The
  production default is therefore hybrid: smooth curved faces retain the
  legacy GP/WENO crossing flux, while direct LLF is activated only by a
  fluid-side discontinuity or an annotated sharp feature.

Consequently, this result supports strong global integral-norm accuracy, but
it does **not** prove a pointwise or first-layer second-order immersed-wall
operator.  The remaining defect is the matched crossing/shell closure and its
grid-phase sensitivity, not the polynomial reproduction of the WLS solver.

## 2. Exact Euler solution

The domain is `[0,1] x [0,1]`.  The solid is a circle centred at
`(0.5,0.5)` with radius `R=0.2`; the fluid is outside the circle.  A fixed
4096-segment polygon represents the same circle at every fluid resolution.

Let `r` and `e_theta` be polar coordinates about the circle centre.  The
steady exact state is

```text
T       = T0,
u       = V e_theta,
p(r)    = p_w (r/R)^a,
rho(r)  = p(r)/(R_g T0),
a       = V^2/(R_g T0),
```

with `T0=300 K`, `V=50 m/s`, `R_g=287.05 J/(kg K)`, and `gamma=1.4`.
The pressure constant is chosen so that `p(R)=100000 Pa`.  This is an exact
source-free Euler solution because continuity and energy are constant along
circular streamlines and radial momentum satisfies

```text
dp/dr = rho V^2/r.
```

At `r=R`, the velocity is exactly tangent to the stationary wall and the
analytic curvature pressure compatibility supplies the same normal pressure
gradient.  Thus this test exercises a nonzero wall pressure gradient rather
than the simpler, physically incorrect `dp/dn=0` closure for curved slip
flow.  Exact conservative states are imposed on every outer boundary.

## 3. Discretization and protocol

The tested configuration is:

```text
equations                         2-D, single-species ideal-gas Euler
bulk convective flux              characteristic LLF-WENO-Z5
time integration                  four-stage third-order RK
IBM interpolation order           2
IBM ghost extrapolation order     2
IBM surface extrapolation order   2
alpha / alpha_surf                1.0 / 1.0
ghost layers                      1
support visibility mode          2
wall                              stationary isothermal slip
pressure closure                  analytic Euler-slip curvature
face-local crossing switch        on
direct-LLF activation threshold   0.02
real-fluid-only shell experiment  off
```

The fixed-time sequence uses:

| N | h | dt | RK steps | final time |
|---:|---:|---:|---:|---:|
| 32 | 0.0312500 | 1.000000e-5 | 200 | 0.002 |
| 48 | 0.0208333 | 6.666667e-6 | 300 | 0.002 |
| 64 | 0.0156250 | 5.000000e-6 | 400 | 0.002 |
| 96 | 0.0104167 | 3.333333e-6 | 600 | 0.002 |

The initial state is the exact solution evaluated at cell centres.  The time
step scales linearly with `h`, so third-order temporal error is smaller than a
second-order spatial target in this range.

The point-centre initialization is sufficient for this second-order IBM gate,
but it is not an exact finite-volume cell-average protocol.  These curves must
not be used to claim formal fifth-order bulk WENO accuracy; that separate gate
requires quadrature-accurate cell averages (or a clearly finite-difference
interpretation) and a no-solid smooth solution.

For the `t=0` diagnostic, `order_rk=0` and one solver step copy the unscaled
semi-discrete conservative RHS into the output state.  Since the exact source
is zero, the exact target is zero.  This test therefore exposes local spatial
consistency without transient cancellation.

The bands are:

- `all_fluid`: every analytically fluid cell;
- `near_wall_1cell`: the first Cartesian fluid layer adjacent to the solid;
- `interior_fluid`: fluid cells excluding that first layer.

## 4. Why the closure is hybrid

The original face-local experiment reconstructed a fluid state at each
Cartesian crossing face, reflected it through the physical wall tangent, and
formed a direct LLF flux at that Cartesian face.  This construction is
conservative at the numerical face, but it imposes the reflected Riemann
problem at the wrong geometric location whenever the immersed wall is offset
from that face.  On the smooth-circle MMS, forcing this closure on every
crossing face produced fixed-time pressure L2 errors

```text
N=32: 46.70,  N=48: 31.24,  N=64: 660.7,
```

and a fitted order of about `-3.51`.  This is a consistency failure, not a
parameter-tuning issue.

The corrected default first computes the established dense GP/WENO flux.  A
sparse face-local pass overwrites it only when either:

1. the normalized fluid-side density/pressure/velocity change is at least
   `cns.llf_ibm_face_local_activation_threshold` (default `0.02`), or
2. the crossing face lies in the compact band of an annotated 2-D sharp
   feature.

Setting the threshold to zero remains an explicit forced-direct ablation.
For the smooth-circle exact state, all crossing faces correctly remained on
the smooth legacy path:

| N | crossing-face copies | smooth legacy | direct LLF | fatal |
|---:|---:|---:|---:|---:|
| 32 | 48 | 48 | 0 | 0 |
| 48 | 80 | 80 | 0 | 0 |
| 64 | 104 | 104 | 0 | 0 |
| 96 | 152 | 152 | 0 | 0 |

This fact is important: the present smooth MMS validates the production
hybrid LLF-WENO-Z5/IBM combination, but it does not rehabilitate direct
reflected LLF as a smooth-wall closure.

## 5. Fixed-time convergence

### 5.1 Fitted orders

| Norm / band | rho | u | v | p | T |
|---|---:|---:|---:|---:|---:|
| L1, all fluid | 2.51 | 2.24 | 2.26 | 2.52 | 2.56 |
| L2, all fluid | 2.30 | 1.82 | 1.82 | 2.33 | 2.41 |
| Linf, all fluid | 1.05 | 0.83 | 0.84 | 0.96 | 1.42 |
| L1, first layer | 1.38 | 1.11 | 1.12 | 1.30 | 1.64 |
| L2, first layer | 1.35 | 1.06 | 1.07 | 1.23 | 1.66 |
| Linf, first layer | 1.05 | 0.83 | 0.84 | 0.96 | 1.42 |

### 5.2 L2 error values

| Band | Var | N=32 | N=48 | N=64 | N=96 | fitted p |
|---|---|---:|---:|---:|---:|---:|
| all fluid | rho | 3.360e-5 | 1.845e-5 | 7.351e-6 | 2.858e-6 | 2.30 |
| all fluid | u | 1.518e-2 | 7.026e-3 | 2.932e-3 | 2.260e-3 | 1.82 |
| all fluid | v | 1.524e-2 | 7.060e-3 | 2.932e-3 | 2.272e-3 | 1.82 |
| all fluid | p | 4.010 | 2.185 | 8.493e-1 | 3.316e-1 | 2.33 |
| all fluid | T | 3.981e-3 | 2.076e-3 | 8.452e-4 | 2.973e-4 | 2.41 |
| first layer | rho | 6.096e-5 | 3.503e-5 | 1.976e-5 | 1.452e-5 | 1.35 |
| first layer | u | 4.017e-2 | 1.822e-2 | 9.898e-3 | 1.351e-2 | 1.06 |
| first layer | v | 4.049e-2 | 1.866e-2 | 9.898e-3 | 1.363e-2 | 1.07 |
| first layer | p | 5.831 | 2.870 | 1.693 | 1.587 | 1.23 |
| first layer | T | 8.902e-3 | 4.942e-3 | 2.857e-3 | 1.465e-3 | 1.66 |

The `N=64 -> 96` first-layer velocity error increases by about 36%.  Reporting
only the `32/48/64` sequence would therefore give an unjustifiably optimistic
near-wall order near two.

## 6. Exact-state t=0 RHS convergence

### 6.1 Fitted L2 orders

| Band | rho RHS | rhou RHS | rhov RHS | rhoE RHS |
|---|---:|---:|---:|---:|
| all fluid | 0.155 | 1.055 | 1.081 | 0.118 |
| interior fluid | 0.275 | 0.785 | 0.842 | 0.164 |
| first layer | -0.376 | 0.569 | 0.593 | -0.408 |

The first-layer L2 values are:

| Component | N=32 | N=48 | N=64 | N=96 |
|---|---:|---:|---:|---:|
| rho | 9.323e-1 | 6.662e-1 | 7.585e-1 | 1.401 |
| rhou | 2.204e2 | 1.213e2 | 9.704e1 | 1.199e2 |
| rhov | 2.276e2 | 1.267e2 | 9.704e1 | 1.217e2 |
| rhoE | 2.731e5 | 1.947e5 | 2.274e5 | 4.236e5 |

The all-fluid L1 fit is better because the defective set has shrinking
measure: 0.797 for mass, about 1.54 for momentum, and 0.761 for energy.  Linf
is non-convergent for mass and energy and only about 0.34 for momentum.

The largest fixed-time velocity errors also occur in the first layer.  At
`N=64` the maximum speed error is about `0.0233 m/s` at a wall distance of
`0.181 h`; at `N=96` it is about `0.0410 m/s` at `0.100 h`.  Symmetry-related
maxima appear at different angular phases as `N` changes.  This identifies a
Cartesian-grid/circle alignment effect in the local wall closure.

Good fixed-time global L1 convergence can coexist with a poor pointwise RHS:
the bad cells form an `O(h)`-measure shell, signs can cancel during evolution,
and the global norm dilutes sparse maxima.  The fixed-time result is useful,
but it cannot be used as proof of a uniformly second-order boundary operator.

## 7. Sharp-corner regression after hybridization

The smooth-wall safeguard must not restore the original compression-corner
pressure defect.  A Mach-2, 20-degree, `N=80` Euler-slip compression-corner run
to `t=14` used the same `0.02` threshold plus mandatory sharp-feature
activation.  In its final steady diagnostics, 25 of 103 crossing faces were
direct LLF due to the feature band; the remaining 78 were smooth legacy.
There were no first-order WLS fallbacks, fatal states, or duplicate-face flux
inconsistencies.

| Metric | Result |
|---|---:|
| shock-angle error | -0.3658 deg |
| downstream flow-angle error | -0.08295 deg |
| p2/p1 error | -0.841% |
| downstream plateau pressure L1 | 0.928% |
| pressure overshoot / p1 | 3.943% |
| pressure L1 outside `k h` shock strip | 0.443% |

Thus the hybrid switch preserves the direct LLF benefit at the sharp corner
while avoiding its smooth-wall inconsistency.

## 8. Verification gates passed

- NVCC CUDA build for `sm_89`, including cross-execution-space checks: pass.
- Four fixed-time GPU runs and four `t=0` GPU diagnostics: finite and complete.
- Boundary-face diagnostics: zero fatal states and zero duplicate-face flux
  inconsistencies in the reported smooth MMS and compression regression.
- Quadratic WLS reproduction: all 2048 2-D masks, 20128 deterministic/random
  3-D masks, and dual face/reflection targets pass.  Maximum moment errors are
  `6.87e-14`, `1.30e-12`, and `2.49e-13`, respectively.
- Python postprocessor syntax and repository `git diff --check`: pass.

This round did not perform a full multi-rank MPI partition-equivalence run.

## 9. Defensible status and next work

The defensible statement is:

> For this smooth curved Euler-slip problem, LLF-WENO-Z5 with the hybrid IBM
> gives at least second-order global L1 and approximately second-order global
> L2 accuracy, but the first near-wall fluid layer and exact-state RHS are not
> uniformly second order.

The next numerical work should target the boundary operator rather than add
ghost layers:

1. derive a wall-location-consistent crossing flux whose discrete divergence
   matches the pressure-compatibility GP reconstruction;
2. pair it with a same-fluid-side shifted WENO-Z5 stencil for the first
   fluid-fluid shell, rather than independently masking solid candidates;
3. require constrained polynomial moments/EOS consistency at both the wall
   and target Cartesian face;
4. repeat at `N=64/80/96/128` under several sub-cell translations of the
   circle to separate asymptotic order from grid phase; and
5. retain the sharp-feature direct-LLF branch as an explicitly non-smooth
   closure, not as the smooth-wall default.

Increasing `ghost_layers` from one to two does not correct the geometric
offset between the wall and target flux face.  It may provide more stored
states, but without a compatible crossing/shell construction it merely moves
the inconsistency into a wider stencil.

## 10. Reproducibility artifacts

The detailed generated tables are under
`IBM/cases/validation/mms_ibm/`:

- `convergence_rhs_shell_llf_wenoz5_feature_hybrid_t0_32_48_64_96_*.md`
- `convergence_solid_llf_wenoz5_feature_hybrid_tfinal_32_48_64_96_*.md`
- `plot/solid_llf_wenoz5_feature_hybrid_t0_N*/plt00001`
- `plot/solid_llf_wenoz5_feature_hybrid_tfinal_N*/plt*`
- `logs/llf_wenoz5_curved_mms_feature/`

The compression regression is in
`IBM/cases/validation/canonical/expansion_corner_2d/runs/`
`llf_wenoz5_hybrid_feature_threshold002_n080_t14_20260714/`.
