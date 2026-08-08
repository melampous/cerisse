# AFD-HLLC and 2-D sharp-corner tangential reconstruction audit

Date: 2026-07-13

## 1. Scope

This note separates two questions that must not be conflated:

1. Is the production AFD-HLLC Euler operator mathematically consistent and
   compatible with the marker-based IBM?
2. Does a same-edge-biased tangential image-point stencil cure the pressure
   error generated at a singular convex polygon corner?

The first answer is yes within the stated Cartesian, single-species ideal-gas
scope. The second answer is only partially: the prototype reduces several
integral errors, but it neither removes the corner singularity nor raises the
observed convergence rate. It therefore remains disabled by default.

## 2. AFD-HLLC audit

### 2.1 Smooth all-fluid path

For an all-fluid six-cell footprint, `src/rhs/Afd.h` performs:

1. fifth-order point-value TENO5 or WENO-Z5 interpolation of left and right
   states;
2. an HLLC interface solve;
3. the Alternative Finite Difference correction assembled from six
   cell-centred physical fluxes.

The optimal point interpolation has coefficients

```text
[3/128, -5/32, 45/64, 15/32, -5/128]
```

at a face. A symbolic moment check reproduces polynomials of degree zero
through four exactly, so its face-value error starts at order `h^5`.

The two AFD correction brackets have Taylor moments

```text
first bracket  = -(h^2/24) F'' + O(h^6)
second bracket = +(7 h^4/5760) F'''' + O(h^6)
```

which agrees with the implemented coefficients. A fresh rebuild and
full-periodic Euler MMS run with the 2026-07-13 source gives fitted `t=0` RHS
orders 4.91--5.18 for AFD-HLLC-TENO5; the prior WENO-Z5 sequence gives
5.09--5.18 for AFD-HLLC-WENO-Z5. Fixed-time tests with temporal error scaled
to order five give 4.99--5.18. These are smooth no-solid results; HLLC alone is
not fifth order and a shock is not expected to converge at fifth order.

### 2.2 Shock and IBM policy

The current operator applies the full AFD correction only to a valid, smooth,
all-fluid six-cell footprint. Its reduction hierarchy is:

1. invalid or rarefaction-pocket stencil: adjacent-state LLF;
2. nonsmooth all-fluid stencil: high-order characteristic LLF when enabled;
3. wall-intersecting stencil: marker-masked point reconstruction followed by
   HLLC;
4. strong reconstructed-GP/fluid jump, nonphysical reconstructed state, or
   nonfinite HLLC flux: adjacent-state LLF;
5. solid-solid face: zero volume flux.

No deep-solid primitive is intentionally read. The wall-intersecting branch
does not apply the full AFD correction, so the near-wall operator is not a
fifth-order AFD closure. Smooth-sphere IBM MMS previously measured an
approximately second-order wall closure. This distinction is essential:
AFD-HLLC is fifth order in a smooth all-fluid region, not globally fifth order
for a domain containing an immersed boundary.

## 3. Convex-corner verification protocol

The corner case is a two-dimensional inviscid Mach-2 Prandtl-Meyer expansion
through 20 degrees, integrated to `t=7` with SSPRK3 and CFL 0.30. The IBM uses

```text
interp_order=2, extrap_order=2, extrap_order_surf=2
alpha=1, ghost_layers=1
adiabatic slip, shock-aware fluid pressure extrapolation
```

The same fixed polygon and sampling windows are used for every scheme. The
Prandtl-Meyer reference is `M2=2.83059518` and `p2/p1=0.27517765`.

At 120 x 144, all tested AFD configurations remained finite and recovered the
correct downstream turning and thermodynamic state. However, AFD-HLLC did not
improve the local corner pressure relative to the existing LLF path:

| Scheme | angle error [deg] | M2 error [%] | plateau p L1 [%] | plateau p L2 [%] | minimum p error [%] |
|---|---:|---:|---:|---:|---:|
| LLF-TENO5 | -0.0096 | +0.033 | 5.366 | 8.391 | -22.51 |
| LLF-WENO-Z5 | -0.0155 | +0.021 | 4.067 | 6.279 | -15.99 |
| AFD-HLLC-TENO5 | -0.0093 | +0.140 | 5.717 | 8.743 | -27.83 |
| AFD-HLLC-WENO-Z5 | -0.0889 | +0.003 | 5.255 | 8.206 | -24.96 |

Thus HLLC's better contact/shear resolution does not by itself repair an IBM
image-point trace that mixes data across a singular corner. The median
Prandtl-Meyer result validates the bulk physics, while the pressure minimum and
surface trace expose the separate near-corner defect.

## 4. Same-edge-biased tangential WLS prototype

The optional `sharp_feature_tangential_reconstruction` path is currently
implemented only for two-dimensional polygon IBM with quadratic WLS for both
volume and surface image points. It is false when the parameter is absent.

The geometry initialization identifies true polygon vertices from turning
angle and stores graph arc distance to the nearest sharp vertex on every edge.
Within `sharp_feature_tangential_radius_cells * h`, the image point, wall
point, normal, and wall-normal distance are left unchanged. Only the lower
corner of the Cartesian 3 x 3 WLS block is shifted by one neighbouring cell in
the owner-edge direction that points away from the nearest corner. Diagonal
edges choose the best aligned of eight neighbouring blocks. A shifted block is
accepted only if it

- stays inside the owning grown FAB;
- has the required number of real-fluid samples; and
- retains a full-rank quadratic WLS matrix.

Otherwise the original centred block is retained. The source-to-target edge
tangent is reconstructed explicitly, so the direction remains correct when
`interior_is_solid=false` has reversed the local frame. The operation occurs
only during IBM/regrid initialization. It stores the same nine weights and
adds no work to the Runge-Kutta hot path.

This is a directional block shift, not yet a strict per-sample half-space
filter. For an unfavourable sub-cell alignment a translated rectangular block
can still overlap the adjacent corner branch. The defensible name is therefore
"same-edge-biased" reconstruction; strict same-edge support requires an
arbitrary-point WLS mask keyed to the nearest feature vertex.

## 5. Controlled A/B results

Only the tangential WLS switch changes between each pair below. Both use
AFD-HLLC-TENO5 and the same pressure closure.

| Grid | WLS | plateau p L1 [%] | plateau p L2 [%] | area below -5% | surface p L1 [%] | pressure-fallback faces |
|---|---|---:|---:|---:|---:|---:|
| 120 x 144 | centred | 5.717 | 8.743 | 1.4931 | 17.322 | 109 |
| 120 x 144 | edge-biased | 5.293 | 8.837 | 1.2986 | 15.215 | 77 |
| 180 x 216 | centred | 3.809 | 6.424 | 0.9877 | 14.093 | 88 |
| 180 x 216 | edge-biased | 3.564 | 6.304 | 0.9630 | 13.447 | 59 |

The edge-biased block reduces plateau L1 by 7.4% and 6.4% on the two grids. It
reduces surface L1 by 12.2% and 4.6%, and pressure-fallback face count by 29%
and 33%. The 180-grid downstream state is also accurate: angle error is
-0.0347 degrees, Mach error is +0.0044%, and `p2/p1` error is -0.368%.

The limitations are equally important:

- plateau L1 order from 120 to 180 is 1.00 centred and 0.98 edge-biased;
- plateau L2 order is 0.76 centred and 0.83 edge-biased;
- surface L1 order is 0.51 centred and 0.30 edge-biased;
- the pointwise pressure minimum is not consistently improved;
- visible wall-parallel pressure bands remain downstream of the corner.

The prototype changes an error constant but does not demonstrate a higher
order method. A singular inviscid corner also precludes a meaningful claim of
uniform pointwise second order at the vertex itself.

## 6. Build and runtime checks

The final source was checked with:

- GCC 2-D AFD-HLLC-TENO5, tangential reconstruction enabled;
- NVCC `sm_89` 2-D AFD-HLLC-TENO5, enabled, followed by the complete 180 x 216
  run;
- GCC 2-D AFD-HLLC-TENO5, feature disabled;
- GCC/CGAL 2-D AFD-HLLC-TENO5, feature enabled, followed by a two-step finite
  run;
- GCC 3-D sphere IBM AFD-HLLC-TENO5, feature absent/default-off, followed by a
  two-step finite run;
- a fresh GCC 2-D no-solid AFD-HLLC-TENO5 rebuild and periodic `t=0` RHS
  sequence at `N=16/32/64/128` (fitted L2 orders 4.91--5.18);
- `git diff --check`.

The 60 x 72 CPU and GPU active-path initialization selected the same 12 volume
image-point blocks and 100 surface blocks. At 180 x 216 the active CUDA path
selected 12 volume blocks and 32 surface blocks, consistent with an `O(h)`
feature band and the fixed surface discretization.

## 7. Rejected branch-restricted experiments

The proposed stricter support was implemented and tested before any production
decision. Two variants were evaluated:

1. a hard half-space selector that replaced the Cartesian block by the nearest
   nine real-fluid samples on the owning side of the corner; and
2. a fixed-block weighted WLS fit that smoothly reduced the influence of
   samples behind the owning edge while retaining positive weights and exact
   polynomial reproduction.

The hard selector completed the 120 x 144 run, but did not improve the result:
plateau L1/L2 were 5.732%/9.097%, low-pressure area was 1.3125, and surface L1
was 17.794%. At 180 x 216 it deterministically produced negative total energy
at `t=0.1944745909` in a fluid cell adjacent to the downstream ramp. The
strong weighted variant (`minimum weight=0.25`) failed at the same location and
nearly the same time (`t=0.1942748654`). The hard-selector failure was
reproduced on CPU and CUDA; the weighted failure was reproduced on CPU.

An equal-weight control (`minimum weight=1`) crossed the failure time and
finished the `t=0.30` isolation run. This shows that the weighted normal-matrix
implementation and call plumbing were sound; the failure was triggered by the
changed ghost trace, not by an indexing or uninitialised-memory defect.

A weak weighting (`minimum weight=0.90`) remained finite to `t=7`. It reduced
some coarse-grid constants but reversed at the third grid:

| Grid | reconstruction | plateau p L1 [%] | plateau p L2 [%] | area below -5% | surface p L1 [%] |
|---|---|---:|---:|---:|---:|
| 120 x 144 | centred | 5.717 | 8.743 | 1.4931 | 17.322 |
| 120 x 144 | weighted | 5.365 | 8.807 | 1.2917 | 16.764 |
| 180 x 216 | centred | 3.809 | 6.424 | 0.9877 | 14.093 |
| 180 x 216 | weighted | 3.757 | 6.407 | 0.9475 | 13.554 |
| 240 x 288 | centred | 3.066 | 5.110 | 0.8160 | 12.521 |
| 240 x 288 | weighted | 3.427 | 5.180 | 0.8837 | 12.582 |

For plateau L1, the weighted pairwise orders fall from 0.88 to 0.32; centred
orders are 1.00 and 0.75. Surface weighted orders fall from 0.52 to 0.26. The
240-grid A/B comparison used the same NVCC `sm_89` binary family, GPU, time
protocol, geometry, and post-processing windows. The weak weighting therefore
does not pass the convergence gate. The hard and weighted branch paths were
removed rather than retained as unsupported tuning switches.

### 7.1 Rejected directional and dual-projection ghost states

The directional ghost-cell literature was also tested rather than copied
literally. Chi et al.'s direction-local values belong to an incompressible
forcing formulation, so its momentum-only use cannot be copied component by
component into compressible continuity, momentum, and energy equations. A
direction-local *complete* ghost state is not inherently non-conservative,
however. It can define a conservative boundary closure if every Cartesian
face evaluates one thermodynamically consistent state and one complete flux,
and that face flux is shared by every cell update that uses the face. What is
invalid is mixing different directional states between equation components or
evaluating two different fluxes on the same face. The tested LLF-TENO5,
LLF-WENO-Z5, and AFD-HLLC-TENO5 whole-state variants all worsened the corner
metrics; a component-wise momentum-only experiment also produced negative
internal energy and was stopped.

A narrower selector was then implemented. It retained a ghost point only when
both Cartesian grid lines crossed the same geometry but different edge
elements whose normals differed by the configured sharp-feature angle. On the
shifted `120 x 144` geometry this selected exactly one corner GP (two direction
slots). A direction-wise quadratic wall/near/far state passed polynomial and
slip-condition formula tests and remained finite, but still worsened the
LLF-TENO5 plateau and surface errors.

Finally, the compressible dual-image construction of De Vanna et al. was
adapted to stationary Euler slip. Each first-hit edge supplied a symmetric
same-line image state; its velocity was reflected about that edge normal; the
two candidates were combined with inverse-square image-distance weights; and
`p,T,u` were completed through the ideal-gas EOS into one ghost state. This
restored one-state mass/momentum/energy coupling and avoided duplicated flux
kernels, but it did not enforce zero normal velocity on either wall branch
after blending.

The controlled shifted-geometry result was:

| LLF-TENO5 closure | plateau p L1 [%] | plateau p L2 [%] | minimum p error [%] | area below -5% | surface p L1 [%] |
|---|---:|---:|---:|---:|---:|
| normal IBM baseline | 2.030 | 3.124 | -7.837 | 0.6875 | 9.150 |
| strict direction-wise quadratic | 2.377 | 3.661 | -9.001 | 0.9444 | 10.512 |
| strict dual projection | 4.777 | 7.146 | -17.641 | 1.3125 | 16.599 |

The dual-projection run still recovered the bulk PM state (`-0.0056` degree
angle error, `+0.106%` Mach error, and `-0.040%` pressure-ratio error), which
isolates the regression to the corner-generated near-wall trace rather than
the far-field physics. The selector, metadata, runtime kernel, case switches,
and extra halo exchange were removed from production source. The run artifacts
remain under `runs/teno5_directional_quadratic_strict_shift050_070_n120_gpu`
and `runs/teno5_dual_projection_strict_shift050_070_n120_gpu` as negative
evidence.

This rules out two simple endpoints: a component-wise transplant of an
incompressible directional forcing method, and an inverse-distance average of
mutually incompatible slip states. It does not rule out a face-wise,
EOS-consistent directional closure. Any future corner treatment must act on a
wall-face state or surface trace, produce one complete shared flux per face,
and pass the shifted-grid PM gate before being added to the runtime operator.

## 8. Near-wall HLLC robustness follow-up

The strong weighted case was retained only as a diagnostic executable after
its source path had been removed. With the original
`cns.afd_ibm_llf_threshold=0.08`, it fails at `t=0.1942748654`. Forcing LLF on
every nonzero crossing jump (`threshold=0`) crosses that time and completes the
full `t=7` run. A threshold scan to `t=0.30` gives:

| Threshold | Result |
|---:|---|
| 0.080 | negative internal energy at `t=0.1942748654` |
| 0.079 | finite to `t=0.30` |
| 0.050 | finite to the complete `t=7` endpoint |
| 0.075, 0.070, 0.060, 0.040, 0.020, 0 | finite to `t=0.30` or beyond |

The transition identifies an approximately 8% GP-fluid jump that the old
default narrowly classified as smooth. It also shows why unconditional LLF is
not the correct production response: the existing smooth-sphere study measured
only 1.34--1.71 order in the first wall shell for that policy.

The selected default is therefore `0.05`. A fresh three-dimensional
wall-compatible sphere sequence at `N=16/24/32`, fixed `t=0.002`, gives fitted
L2 orders:

| Region | rho | u | v | w | p | T |
|---|---:|---:|---:|---:|---:|---:|
| all fluid | 4.71 | 4.54 | 4.45 | 4.34 | 4.70 | 4.71 |
| first wall cell | 4.37 | 4.00 | 3.87 | 3.45 | 4.37 | 4.38 |

Thus `0.05` catches the sharp-corner stress jump without reproducing the
smooth-wall order loss of all-crossing LLF. It is still a sensor-based
robustness policy, not a formal positivity-preserving update.

The current source was then rebuilt and smoke-tested without a runtime threshold
override: 2-D GCC, 2-D NVCC `sm_89`, and 3-D GCC all compiled and advanced a
finite IBM state. These tests verify build and dispatch coverage; the
`16/24/32` fixed-time sequence above is the accuracy evidence.

## 9. Decision and next work

AFD-HLLC-TENO5/WENO-Z5 can be used with the current Cartesian IBM, subject to
the local order reduction and robustness policy above. AFD-HLLC is not the
solution to the convex-corner pressure defect and should not replace LLF-WENO
solely on the basis of this case. The older translated-block tangential
prototype remains default-off because it improves constants but not order.

The next engineering priority is the near-wall AFD-HLLC positivity mechanism.
Checking reconstructed states and the instantaneous HLLC flux is insufficient
to guarantee positive cell averages after a Runge-Kutta update. A defensible
implementation should form high-order and positivity-preserving low-order
face fluxes, compute a cell admissibility factor for density and internal
energy, and use one conservative face factor shared by both adjacent cells.
The limiter should be active only where required; smooth all-fluid and smooth
IBM MMS cases must report a factor of one and retain their measured order.

Only after that robustness gate should corner reconstruction be revisited. The
appropriate target is a surface-trace or face-state construction tied to the
owner edge, with corner-excluded norms and shrinking `k h` annuli. More ghost
layers or stronger WLS weights do not remove the inviscid vertex singularity
and are not substitutes for a conservative near-wall flux treatment.

The next corner implementation is deliberately split into four gates:

1. **Geometry map.** For every direction and fluid/solid grid face, store the
   first segment hit seen from the fluid cell, global element id, side id,
   intersection point, outward owner-face normal, fluid/ghost distances, and a
   validity flag. A nearest-Euclidean-face result is not an admissible
   substitute. Thin surfaces retain two side ids. This map is rebuilt only
   when geometry or AMR topology changes.
2. **Complete face state.** Construct one ideal-gas primitive state for each
   direction/side. Pressure, temperature, all velocity components, density,
   sound speed, and total energy must be completed through the same EOS. Slip
   uses the first-hit face normal; no equation component may select a different
   owner edge. Invalid, nearly coincident, or positivity-threatening
   reconstruction falls back as a complete state, not component by component.
3. **Flux integration.** LLF-WENO/TENO and AFD-HLLC consume the direction-local
   state only while evaluating that direction's faces. Each face writes one
   mass/momentum/energy flux. All-fluid stencils remain bitwise unchanged, and
   smooth-wall stencils retain the existing measured order. The first
   implementation may use explicit directional scratch storage for auditability;
   a compact GP-index overlay is the production optimisation after validation.
4. **Acceptance.** Unit tests cover first-hit selection, side separation,
   polynomial reproduction, slip reflection, EOS closure, and finite fallback.
   Runtime gates are no-solid bitwise equality, smooth-wall MMS, shifted PM
   grids at `120/180/240`, positivity, and fluid-domain mass/energy plus IBM
   momentum-budget closure. A corner option is retained only if local pressure
   norms improve without degrading the bulk PM state or grid-convergence
   trend.

## 10. Direct-state checked HLLC prototype (2026-07-14)

The next prototype now follows one face-local Riemann path on each Cartesian
fluid/solid crossing face:

1. Find the first hit and owning surface element from the fluid cell to the
   solid cell. The owning-element normal defines wall reflection; the
   Cartesian face direction still defines the Euler flux.
2. Use only visible fluid supports. One WLS polynomial is evaluated at the
   Cartesian face and at that face's reflection through the true wall plane.
3. Reflect the second primitive state about the wall normal, recompute sound
   speed, and rebuild both conservative energies through the same ideal-gas
   EOS. No shared solid-cell primitive participates in this flux.
4. Call HLLC once. Check finite ordered wave speeds, both star densities,
   star internal energies and pressures, and every flux component. Reject an
   inadmissible candidate to Einfeldt-HLLE; if that also fails, retry once from
   the adjacent first-order fluid state. A final failure writes NaN so the
   global finite-state diagnostic stops the run.

The initially proposed condition `rho_f=rho_g` and `p_f=p_g` at the Cartesian
crossing face is not used for an oblique immersed wall. It imposes parity at
the staircase face rather than at the true boundary. A strict one-state
face-target experiment failed near `t=0.887`; treating the reflected point as
the face state failed near `t=0.145`. The distance-aware pair without a
near-wall fluid shell survived longer but failed near `t=1.858`. These are
negative controls, not production choices.

For non-crossing fluid/fluid faces whose nominal AFD stencil intersects the
body, the prototype currently uses marker-safe MUSCL traces: MC when both
one-cell neighbours are fluid and one-sided minmod when the wall blocks one
side. These traces use the same checked HLLC/HLLE path, with adjacent-cell LLF
only as a last fallback. Thus the crossing face never uses masked five-point
TENO/WENO, while the nearby shell is second-order rather than fifth-order.

### 10.1 WLS and runtime audits

The face and reflected-face targets retain independent requested/accepted
ranks and condition numbers. This is required because translating the centred
polynomial basis changes conditioning. In the final GPU smoke test their
maximum condition numbers were approximately `4.23e2` and `1.00e4`,
respectively, while both remained below the configured `1e10` limit.

Runtime first-order selection is split into three causes: invalid static
metadata, invalid/non-partitioning weights, and a non-admissible reconstructed
primitive state. HLLC and HLLE rejection reasons remain separate. The `80 x
96`, `t=7` diagnostic run used WLS and HLLC on all 107 crossing-face copies in
all 2064 RK stages. All 212 near-wall shell copies per stage also used HLLC.
There were no metadata/weight/state fallbacks, no HLLE or LLF selections, no
fatal states, and no numerical or bitwise duplicate-face flux mismatch.

The completed `180 x 216`, `t=7` run gives the high-resolution qualification
that the coarse smoke test cannot provide. Across 6168 RK stages it processed
242 crossing-face copies per stage (1,492,656 total). WLS supplied 1,372,057
states (`91.92%`). Exactly two copies per stage fell back because of static
metadata (`0.83%`); no copy ever failed the weight partition check. A further
108,263 copies (`7.25%`, between 0 and 28 per stage) fell back because the
quadratic WLS evaluation produced a non-admissible primitive state. Every
admissible pair, including each first-order replacement, passed checked HLLC:
there were no HLLC rejections, HLLE selections, retries after the Riemann
solver, fatal states, or duplicate-flux mismatches. All 2,991,480 near-wall
fluid/fluid shell evaluations also selected HLLC without fallback.

This separates two issues that were previously conflated. The checked
HLLC-to-HLLE logic is not causing high-resolution dissipation in this run;
loss of primitive-state admissibility before the Riemann solve is. The
non-monotone `Nx=180` plateau error is consistent with that observation, but
causation still requires a spatial map of the fallback faces. The next
crossing-face experiment must therefore use an admissibility-preserving
polynomial scaling toward a valid fluid state, while retaining the same wall
reflection and EOS completion. It must not reintroduce a ghost/fluid jump
sensor.

The polynomial test covers all 2-D support masks, sampled 3-D masks, and
random non-symmetric face/reflection targets. Maximum moment errors were
`6.87e-14`, `1.95e-12`, and `2.86e-13`; all tests passed.

### 10.2 PM co-refinement result

All runs use AFD-HLLC-TENO5, `iorder=eorder=eorder_surf=2`, `alpha=1`, and one
ghost layer. Errors are relative to the analytic Prandtl-Meyer state.

| Nx | angle error [deg] | Mach error [%] | p2/p1 error [%] | plateau p L1/L2 [%] | plateau minimum [%] | all-field p L1/L2 |
|---:|---:|---:|---:|---:|---:|---:|
| 80  | -0.0044 | +0.069 | -0.015 | 0.310 / 0.526 | -2.017 | 5.599e-3 / 1.254e-2 |
| 120 | -0.0098 | +0.031 | -0.061 | 0.282 / 0.438 | -1.474 | 3.822e-3 / 9.019e-3 |
| 180 | -0.0449 | +0.077 | -0.328 | 0.447 / 0.711 | -3.314 | 2.654e-3 / 6.469e-3 |

Observed orders for `80 -> 120` and `120 -> 180` are `0.94, 0.90` in global
L1 and `0.81, 0.82` in global L2. Excluding a fixed radius `r0=0.5` gives
`0.95, 0.90` in L1 and `0.87, 0.87` in L2. This is consistent with a PM field
containing a vertex, fan head, and fan tail; it is not evidence of formal IBM
order. The non-monotone `Nx=180` plateau error also prevents claiming that the
corner closure is converged.

Relative to Legacy, the new path removes the severe wall-sector defect:
Legacy had about `3.81%` plateau L1 error, `-25.45%` minimum pressure error,
and `19.7%` of plateau cells below theory by more than 5%. The new sequence has
`0.28--0.45%` plateau L1, `-1.47--3.31%` minimum error, and zero such cells.
Global L1/L2 at `Nx=80` are nevertheless slightly worse than Legacy because
fan-head/tail smearing dominates those norms.

### 10.3 Engineering status and remaining gates

GCC, NVCC `sm_89`, AFD-HLLC-TENO5, and AFD-HLLC-WENO-Z5 builds pass. With
runtime diagnostics enabled and disabled, complete `Nx=80` plotfiles are
bitwise identical. Full-stage diagnostics add about 16.4% on this small GPU
case (`39.28 s` versus `33.75 s`) and are default-off.

The instrumented `Nx=180` run completed 2056 steps in `234.10 s`. Its final
plotfile was also compared with the earlier build of the same numerical path.
The files are not bitwise identical, so that comparison is not used as a
bitwise-regression claim. After decoding the plotfiles, however, the maximum
fluid-cell differences were `5.12e-12` in total energy, `3.47e-12` in
x-momentum, and `6.26e-13` in pressure; `sld` and `ghs` were identical, and
all reported PM metrics agreed at their stored precision. This is
roundoff-level numerical equivalence between builds. The same-executable
`Nx=80` diagnostic-on/off comparison remains the valid bitwise test.

A no-solid 3-D AFD-HLLC-TENO5 MPI+CUDA build also passes. Its one-step run
could not be completed locally because both `mpirun -np 1` and OpenMPI
singleton startup stalled before application output under this WSL setup; the
launcher was terminated after a bounded timeout. This is recorded as an
environment-blocked runtime check, not a solver pass or failure.

The current result proves robust stationary, single-species, ideal-gas
Euler-slip operation for this 2-D PM geometry. It does not yet prove:

- formal second-order IBM accuracy; use smooth-wall MMS for that claim;
- a fifth-order near-wall operator; the marker-safe shell is MUSCL;
- a globally positivity-preserving SSPRK update; star-state checks are only a
  necessary local gate;
- MPI-wide single ownership of one crossing flux or fluid/control-volume
  momentum closure;
- moving walls, multispecies states, viscous no-slip walls, or 3-D feature
  edges.

Those five follow-ups have now been executed.  The complete implementation,
negative results and acceptance boundaries are recorded in
`docs/ibm_face_local_12345_20260714.md`.  In summary, the fallback map,
admissibility scaling, positive-weight shifted shell, smooth-wall MMS and
MPI/CUDA five-equation budget are complete.  The direct face-state path gives
an excellent PM pressure trace but fails the smooth-wall MMS convergence gate;
it therefore remains opt-in.  The legacy GP path retains the demonstrated
global second-order result.
