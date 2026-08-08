# Pure GP-IBM Gate 2 protocol: inverse-MOC smooth nozzle

Date: 2026-07-22

## Scope

This gate certifies only the frozen, pure shared-GP, full-Cartesian method:

- two-dimensional Cartesian, level-0, ideal-gas Euler flow;
- fixed smooth Euler-slip walls;
- shared GP states from cell-average-aware BI-CWLS;
- curvature-compatible pressure and one-sided entropy extension;
- full-cell LLF-WENO-Z5 and SSPRK(4,3);
- the full-cell shared-face positivity limiter;
- volume extrapolation order 1 and surface recovery order 2.

Cut-cell, cut-control-volume, aggregate, volume-fraction, area-fraction,
redistribution, and AMReX-EB flow-update paths are excluded. Historical results
from those paths are not admissible evidence for this gate.

## Exact nozzle oracle

The verification problem is the source-free, continuous-curvature planar
inverse-MOC solution in
`IBM/cases/validation/canonical/inverse_moc_nozzle_2d`.
It consists of two separated Prandtl-Meyer simple waves and uniform inlet,
intermediate, and exit states:

\[
M_0=1.2,\qquad M_m=1.6036729033052766,\qquad M_e=2.0.
\]

The two walls are exact streamlines. Their curvature vanishes continuously at
the joins to straight inlet and exit sections. The fixed geometry contains
20,003 vertices per wall. The reference invariants are

\[
\dot m_{\rm ref}=0.6645041151860959,
\qquad
F_{p,x}^{\rm ref}=-0.13764505275517683.
\]

The geometry/oracle metadata reports no construction failures, an inverse
Prandtl-Meyer residual of `2.94e-13`, a mass-flow mismatch of `1.67e-16`, and a
surface/control-volume pressure-force mismatch of `2.37e-10`.

## Frozen numerical matrix

The primary matrix uses one executable and one physical geometry:

| Grid | Phase 0 | Phase (0.37h, 0.23h) |
|---|---|---|
| 256 x 64 | reuse current-source Gate 0 evidence | reuse Gate 0 evidence |
| 384 x 96 | run | run |
| 576 x 144 | run | run |

All cases end at `t=0.5`, use `CFL=0.25`, and retain the same physical geometry.
The translated-grid sequence changes only the Cartesian grid phase.

The executable is
`main2d.gnu.TPROF.MPI.CUDA.inverse_moc_nozzle_llf_wenoz5_i2_e1_s2_a0p6_gl1_vis2_bic_cellavg_sjet_pcurv.src84db3efb00517a42.ex`,
with SHA-256
`cfcb309b66914c8bf2012f23cd59eb1997deffaf5f7454157a2f4d790d5a7020`.

## Boundary-condition audit

The left boundary is an exact manufactured-state `ext_dir` boundary evaluated
from the current SSPRK stage geometry and the analytic MOC field. The right
boundary lies in the uniform Mach-2 region and is a supersonic extrapolation
boundary. The upper and lower Cartesian boundaries are covered by the solid
carrier; the active fluid communicates with the exterior only through the
inlet and outlet.

One finest-grid boundary-location sensitivity case keeps `dx`, `dy`, the wall,
and the MOC solution fixed while adding 12 cells to each straight x-buffer.
The domain therefore contains 600 x 144 cells. This isolates boundary location
from spatial resolution.

## Required observables

The report must include cross-sectional profiles of

\[
\dot m(x),\qquad H_t(x),\qquad
\sigma(x)=\ln p-\gamma\ln\rho,\qquad M(x),
\]

and boundary-intercept profiles of

\[
p_B(s),\qquad u_n(s).
\]

The two principal integral diagnostics are

\[
\epsilon_{\dot m}
=\frac{\max_x\dot m(x)-\min_x\dot m(x)}{\overline{\dot m}},
\qquad
\epsilon_{\rm leak}
=\frac{\left|\int_\Gamma \rho u_n\,dS\right|}{\dot m_{\rm inlet}}.
\]

Physical wall pressure and pressure force come from surface-e2 BI recovery.
Cartesian fluid-solid crossing momentum exchange is retained only as a
numerical diagnostic and is not an acceptance measure for physical traction.

## Acceptance and stop conditions

The gate requires all of the following:

1. Every smooth run has zero positivity-limiter activation and no BI-CWLS
   downgrade from the quadratic reconstruction.
2. Field, near-wall, BI-pressure, impermeability, mass-flow, total-enthalpy,
   entropy, and pressure-force errors converge consistently under refinement.
3. The normalized pressure sensor `h |grad p| / p` decreases with refinement;
   no persistent recompression discontinuity is permitted.
4. Both grid phases approach the same continuum limit. Fine-grid phase spread
   must be no larger than the estimated fine-grid discretization uncertainty.
5. Moving the inlet and outlet within the exact straight buffers changes the
   solution by less than the fine-grid discretization uncertainty.
6. The active-fluid Cartesian mass budget closes to roundoff, while BI
   reconstructed leakage and section mass-flow variation decrease with grid
   refinement.

Any limiter activation, BIC downgrade, persistent recompression wave, distinct
phase limit, excessive boundary-location sensitivity, or nonconvergent mass
flow/BI load is a stop condition. A failure is first assigned to the MOC
oracle, geometry, or physical boundary conditions. It does not authorize a
change to GP-IBM, BI-CWLS, WENO-Z5, or the limiter.

## Permitted claim

A PASS permits only the following statement:

> The frozen pure GP-IBM is a restricted-production candidate for smooth,
> two-dimensional Cartesian, fixed-wall, ideal-gas Euler-slip nozzles.

It does not certify shocked nozzles, AMR, Navier-Stokes wall transport, R-Z,
three-dimensional STL geometry, or exact conservation on the true curved
control volume.
