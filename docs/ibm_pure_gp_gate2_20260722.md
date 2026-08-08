# Pure GP-IBM Gate 2: exact inverse-MOC smooth-nozzle qualification

Date: 2026-07-22

## Decision

\[
\boxed{\text{Gate 2: PASS within the frozen two-dimensional smooth-nozzle scope}}
\]

This result qualifies the frozen method as a restricted-production candidate
for smooth, two-dimensional Cartesian, level-0, fixed-wall, ideal-gas
Euler-slip nozzles. It does not qualify shocked nozzles, AMR, R-Z, viscous
wall transport, general three-dimensional STL geometry, or moving solids.

Gate 2 is an end-to-end engineering qualification against an exact
two-dimensional inverse-method-of-characteristics solution. It is not a
replacement for the formal-order evidence in Gate 1. In particular, this
report does not claim that every Gate-2 norm is uniformly second order.

## Frozen method and provenance

Every run used the pure shared-GP full-Cartesian method:

- one shared state at each real solid-side ghost-cell centre;
- cell-average-aware BI-CWLS with curvature-compatible pressure and a
  one-sided entropy extension;
- one configured solid-side ghost-cell layer;
- full Cartesian cell volumes and face areas;
- LLF-WENO-Z5 and SSPRK(4,3);
- full-cell shared-face positivity limiting;
- volume extrapolation order 1 and BI surface recovery order 2.

The runtime manifest reported
`family=pure_shared_gp_full_cartesian`, `pure_gp_production=1`, and
`cut_control_compiled=0` in every new run. Cut-cell, cut-control, aggregate,
volume/area-fraction, redistribution, and AMReX-EB flow updates were excluded.

The common executable was

```text
main2d.gnu.TPROF.MPI.CUDA.inverse_moc_nozzle_llf_wenoz5_i2_e1_s2_a0p6_gl1_vis2_bic_cellavg_sjet_pcurv.src84db3efb00517a42.ex
```

with SHA-256
`cfcb309b66914c8bf2012f23cd59eb1997deffaf5f7454157a2f4d790d5a7020`.
Its manifest records Git commit
`dc35646d8eb90aaa362a3292a9d702e30c9d5af5`, source-tree SHA-256
`384996cc17fd21e1ae34c8643bb619b8a92bcf90ee103470b97c19a03bbc9dc4`,
and AMReX `25.12-21-gbd922c6216e0`.

The fixed 20,003-vertex walls have SHA-256 values
`90f47995d00597f61740fe31cd0ab9eea7941c9f5acb820f9d8a37001af834ce`
and
`bf2c1666f841279df46e86050e17a4f27213b8843119e9aa7c115251204e7008`.

## Exact inverse-MOC oracle

The source-free steady Euler solution contains two separated smooth
Prandtl-Meyer simple-wave regions and uniform inlet, intermediate, and exit
states:

\[
M_0=1.2,\qquad M_m=1.6036729033052766,\qquad M_e=2.0.
\]

The opposite walls are exact streamlines, their curvature vanishes at each
straight/curved join, and the exit state is uniform and parallel. The exact
integral quantities are

\[
\dot m_{\rm ref}=0.6645041151860959,
\qquad
F_{p,x}^{\rm ref}=-0.13764505275517683.
\]

The oracle construction reports an inverse Prandtl-Meyer residual of
`2.94e-13`, opposite-wall tangency errors below `4.52e-9 rad`, a mass-flow
mismatch of `1.67e-16`, and a surface/control-volume pressure-force mismatch
of `2.37e-10`.

The left boundary uses the exact MOC conservative cell average. The right
boundary is in the uniform Mach-2 supersonic region and uses extrapolation.
The upper and lower Cartesian boundaries lie inside the solid carriers.

## Verification matrix

The main matrix used three grids and two normalized Cartesian phases:

| Grid | Phase |
|---|---|
| `256 x 64` | `(0,0)` and `(0.37h,0.23h)` |
| `384 x 96` | `(0,0)` and `(0.37h,0.23h)` |
| `576 x 144` | `(0,0)` and `(0.37h,0.23h)` |

All runs used `CFL=0.25` and ended at `t=0.5`. The `256 x 64` results are the
current-source Gate-0 runs; the other four were generated specifically for
Gate 2 with the same executable and geometry.

Across the six-run matrix, all 4,864 ghost targets used quadratic BI-CWLS and
none used the linear fallback. The four new runs explicitly record 3,856
Forward-Euler brackets with `mode=HIGH_ORDER` and `theta_min=1`; every
corresponding SSP stage audit passed. The two Gate-0 coarse runs used quiet
logging, but limiter activation/fallback messages are unconditional in the
production path and none occurred. No non-finite or inadmissible active-fluid
state was observed.

## Volume-field convergence

The following values are the maximum over the two grid phases. Errors are
relative L2 norms against conservative-cell-average MOC data.

| Quantity | `256 x 64` | `384 x 96` | `576 x 144` | Adjacent orders |
|---|---:|---:|---:|---:|
| fluid density | `1.618e-4` | `7.738e-5` | `3.901e-5` | `1.82, 1.69` |
| fluid x velocity | `4.592e-4` | `2.703e-4` | `1.626e-4` | `1.31, 1.25` |
| fluid y velocity | `8.444e-4` | `4.077e-4` | `2.194e-4` | `1.80, 1.53` |
| fluid pressure | `2.131e-4` | `1.042e-4` | `5.317e-5` | `1.76, 1.66` |
| fluid temperature | `6.927e-5` | `3.147e-5` | `1.590e-5` | `1.95, 1.68` |
| first-wall-layer pressure | `2.740e-4` | `1.440e-4` | `7.233e-5` | `1.59, 1.70` |

All fields converge to the exact two-dimensional solution, but this sequence
does not exhibit a uniform second-order velocity norm. The first-wall-layer
x- and y-velocity adjacent orders are approximately `0.78--1.02`. This is
reported as an observed property of this finite inverse-MOC sequence, not
hidden by an integrated quantity. Formal approximately second-order smooth
curved-wall complete-state accuracy is supplied independently by Gate 1.

The maximum numerical pressure sensor

\[
h|\nabla p|/p
\]

decreases as `0.05404, 0.03726, 0.02540`, with adjacent orders `0.92` and
`0.95`. The exact MOC cell averages give nearly identical values. This is
consistent with a bounded pressure gradient and provides no evidence of a
grid-persistent recompression shock.

## Section invariants

Section quantities were integrated between the analytic walls with a
16-point Gauss rule and a linear endpoint trace from the active-fluid cell
averages. The table reports the two-phase envelope.

| Diagnostic | `256 x 64` | `384 x 96` | `576 x 144` | Adjacent orders |
|---|---:|---:|---:|---:|
| mass-flow spread, `epsilon_mdot` | `3.892e-4` | `1.922e-4` | `1.098e-4` | `1.74, 1.38` |
| total-enthalpy L2 | `1.594e-4` | `6.858e-5` | `3.502e-5` | `2.08, 1.66` |
| entropy L2 | `1.038e-5` | `4.486e-6` | `1.996e-6` | `2.07, 2.00` |
| area-mean Mach L2 | `8.823e-5` | `3.663e-5` | `1.847e-5` | `2.17, 1.69` |

At the finest grid, the section-mean mass-flow errors are `2.74e-6` and
`2.07e-6` in the two phases. An independent 12-point/constant-endpoint audit
also shows decreasing mass-flow spread, reaching `1.50e-4` and `1.41e-4`.
The difference between these post-processing choices is retained as part of
the finite-grid section-integration uncertainty.

## BI pressure, force, and impermeability

The physical pressure observable is the surface-e2 BI recovery. On the finest
grid, the worst wall values over both phases are

\[
\|C_p-C_{p,\rm ref}\|_{L_1}=2.43\times10^{-5},
\]

\[
\|C_p-C_{p,\rm ref}\|_{L_\infty,\,nonendpoint}
=3.13\times10^{-4}.
\]

The phase-envelope BI `Cp` L1 values are
`1.33e-4, 5.15e-5, 2.43e-5`, with adjacent orders `2.34` and `1.85`.
The worst BI normal-velocity L2 values are
`2.52e-4, 8.86e-5, 3.20e-5`, with adjacent orders `2.58` and `2.51`.

The BI pressure-force relative errors on the finest grid are `8.36e-6` and
`1.01e-5`; their phase spread relative to the exact force is `1.76e-6`.
The integrated force error is cancellation-sensitive and is not monotone in
the last pair, so no formal force order is assigned. The local BI pressure
norms converge, and both the force error and phase spread are substantially
below the predeclared one-percent and 0.25-percent engineering limits.

The required signed leakage diagnostic is

\[
\epsilon_{\rm leak}
=\frac{\left|\int_\Gamma\rho u_n\,dS\right|}{\dot m_{\rm inlet}}.
\]

Because upper- and lower-wall errors cancel, its individual phase sequences
are not monotone. The two-phase envelope decreases as
`7.31e-5, 2.40e-5, 1.62e-5`. The stronger non-cancelling diagnostic

\[
\frac{\int_\Gamma|\rho u_n|\,dS}{\dot m_{\rm inlet}}
\]

decreases as `9.91e-4, 3.21e-4, 1.45e-4`, with adjacent orders `2.78` and
`1.96`. Thus impermeability error decreases rather than approaching a fixed
percentage.

## Conservation semantics

The maximum normalized active-fluid Cartesian mass-budget residual is
`1.67e-13`. This proves consistency of the full-cell Cartesian flux ledger;
it is not machine-precision conservation on the true curved fluid volume.

The Cartesian fluid-solid crossing exchange remains a staircase-domain
diagnostic. Its signed mass exchange and momentum force do not equal the BI
surface integrals at finite grid spacing. In particular, the Cartesian
momentum-exchange force remains several percent from the physical BI pressure
force and is not used as a local traction or a Gate-2 physical-load criterion.
No cut-cell correction, global source, clipping, or a posteriori force
correction was applied.

## Grid-phase and boundary-location sensitivity

The finest-grid surface-force phase spread is `1.76e-6` relative to the exact
force. Fine-grid pressure L2 errors differ by about `1.3%` between phases and
both sequences approach the same exact solution. Fine-grid signed and
absolute leakage differ by `0.6%` and `5.1%`, respectively, while their
phase envelopes decrease under refinement.

For the boundary-location test, 12 cells were added to each straight x-buffer
without changing `dx`, `dy`, the wall, or the exact field. All 1,200 GP targets
remained quadratic and all 1,156 Forward-Euler brackets had `theta=1`.
On the overlapping active-fluid domain, the baseline/extended pressure
difference is

\[
\|\Delta p\|_{L_2}=3.11\times10^{-9},\qquad
\frac{\|\Delta p\|_{L_2}}{\|p\|_{L_2}}=1.31\times10^{-8},
\]

and `Linf=2.86e-7`. The BI pressure-force change is `8.57e-14` relative to
the exact force, and the mass-flow-spread change is `1.85e-11`. These effects
are far below the finest-grid discretization errors.

## Final capability statement

Gate 2 supports the following bounded claim:

> For the frozen LLF-WENO-Z5/SSPRK(4,3) configuration, the pure shared-GP
> full-Cartesian method gives a stable, phase-robust and boundary-location-
> insensitive approximation of a source-free, continuous-curvature,
> two-dimensional inverse-MOC Euler-slip nozzle. BI pressure, pressure force,
> impermeability, mass flow, total enthalpy, entropy and Mach diagnostics
> converge to the exact solution without limiter activation or BI-CWLS
> downgrade.

The next gate is a prescribed-back-pressure shocked nozzle. Gate 3 must
separately verify shock position, Rankine-Hugoniot states, total-pressure loss,
outlet-location sensitivity, wall pressure, thrust, and positivity. Gate 2
does not authorize changes to the frozen IBM algorithm and does not certify
AMR or SRP.

## Evidence

- [Gate-2 protocol](ibm_pure_gp_gate2_protocol_20260722.md)
- [phase-0 field/surface analysis](../IBM/cases/validation/canonical/inverse_moc_nozzle_2d/runs/gate2_pure_20260722/analysis_phase000.json)
- [shifted-phase field/surface analysis](../IBM/cases/validation/canonical/inverse_moc_nozzle_2d/runs/gate2_pure_20260722/analysis_phase037_023.json)
- [phase-0 section/wall profiles](../IBM/cases/validation/canonical/inverse_moc_nozzle_2d/runs/gate2_pure_20260722/profiles_phase000/summary.json)
- [shifted-phase section/wall profiles](../IBM/cases/validation/canonical/inverse_moc_nozzle_2d/runs/gate2_pure_20260722/profiles_phase037_023/summary.json)
- [boundary-location analysis](../IBM/cases/validation/canonical/inverse_moc_nozzle_2d/runs/gate2_pure_20260722/analysis_boundary_extended.json)
- [convergence-envelope figure](../IBM/cases/validation/canonical/inverse_moc_nozzle_2d/runs/gate2_pure_20260722/figures/gate2_convergence_envelopes.png)
- [section-invariant figure](../IBM/cases/validation/canonical/inverse_moc_nozzle_2d/runs/gate2_pure_20260722/figures/gate2_section_invariants.png)
- [BI pressure/impermeability figure](../IBM/cases/validation/canonical/inverse_moc_nozzle_2d/runs/gate2_pure_20260722/figures/gate2_wall_pressure_impermeability.png)
