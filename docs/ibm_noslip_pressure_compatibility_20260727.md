# Pure GP-IBM Stationary No-Slip Pressure Compatibility

Date: 2026-07-27

## Scope

This report records the local RTX 4090 implementation and verification of the
normal-pressure compatibility condition for the frozen pure GP-IBM path. The
method remains:

- one unique shared state at every real solid-side ghost-cell center;
- cell-average-aware BI-CWLS;
- full Cartesian control volumes and face areas;
- LLF-WENO-Z5;
- SSPRK(4,3).

No cut-cell, embedded-boundary, volume-fraction, open-subface, aggregate, or
cut-control update is used.

The verified wall model is limited to a fixed, stationary, isothermal,
no-slip wall in two-dimensional Cartesian, single-fluid Navier-Stokes flow.

## Compatibility Condition

At a stationary no-slip wall, with no wall-normal body force or manufactured
momentum source, the wall velocity and convective acceleration vanish. The
normal momentum equation gives

```text
dp/dn = n_i d(tau_ij)/dx_j,
tau_ij = mu*(du_i/dx_j + du_j/dx_i - 2*theta*delta_ij/3)
         + xi*theta*delta_ij,
theta = div(u).
```

For temperature-dependent transport, the implementation evaluates

```text
div(tau)_i = mu*laplacian(u_i)
           + grad(mu)_j*(du_i/dx_j + du_j/dx_i)
           + (xi + mu/3)*d(theta)/dx_i
           + (d(xi)/dx_i - 2*d(mu)/(3*dx_i))*theta.
```

A quadratic weighted least-squares polynomial is fitted to the current-stage
fluid velocity, `mu(T)`, and `xi(T)`. The resulting `dp/dn` is inserted into
the pressure Neumann functional of the same BI-CWLS extension that produces
the unique shared GP state.

## Support Selection

The established complete first-image `3^D` support is attempted first. If it
is unavailable, the direct BI-CWLS path retries the compatibility fit on the
already accepted visible quadratic support:

1. support indices are taken from `constrained_support_ijk`;
2. the mask is the union of nonzero constrained Dirichlet and Neumann
   functionals;
3. the fit is evaluated at the first image point;
4. rank and factorization must succeed;
5. the normal-matrix condition estimate must not exceed
   `support_visibility_condition_max`.

Failure retains the homogeneous-Neumann pressure extension. The retry does not
alter support visibility, create a second ghost state, or introduce a
face-local boundary flux.

Implementation:

- `src/ibm/ibm_solver_interp.h`: reusable selected-support compatibility fit;
- `src/ibm/ibm_solver.h`: support retry and pressure Neumann coupling.

## Build Provenance

Repository HEAD at test time:

```text
dc35646d8eb90aaa362a3292a9d702e30c9d5af5
```

The worktree was intentionally left uncommitted. The no-slip MMS executable
was:

```text
IBM/cases/validation/mms_ibm/
main2d.gnu.TPROF.MPI.CUDA.mms_ns-viscous-compat_solid_llf-wenoz5_
o4_e2_s2_a0p6_as0p6_gl1_i2_vis2_gp_bic_aframe_
expanded_bic_support_iso_pnsvisc.ex
```

SHA-256:

```text
56d0d23251b83f90b7a5c67cb7b36f770fa2e181f5c8adbda6e3c92986c97e8b
```

Execution used one RTX 4090, one MPI rank, and no concurrent simulations.

## Smooth No-Slip MMS

Four grids, `N=64,96,128,192`, were advanced to `t=0.002`. All GP
reconstructions remained quadratic and no fallback occurred.

### Fixed-Time Solution Error

| Region | rho | u | v | p | T |
|---|---:|---:|---:|---:|---:|
| all active fluid | 3.10 | 3.05 | 3.03 | 3.11 | 3.13 |
| first wall-adjacent layer | 2.72 | 2.77 | 2.65 | 3.19 | 2.84 |
| shifted phase, all fluid | 2.85 | 2.87 | 2.82 | 2.83 | 2.91 |
| shifted phase, first layer | 2.62 | 2.53 | 2.52 | 3.11 | 2.59 |

These observed rates establish second-order-or-better fixed-time consistency
for the complete primitive state in this smooth two-dimensional test. They are
not claims of a third-order boundary method.

### Semi-Discrete `t=0` Residual

| Region | rho | rho-u | rho-v | rho-E |
|---|---:|---:|---:|---:|
| all active fluid | 1.94 | 2.54 | 2.47 | 2.21 |
| first wall-adjacent layer | 1.41 | 2.02 | 1.96 | 1.69 |
| more than three cells from wall | 4.23 | 4.81 | 3.76 | 4.18 |

The wall-adjacent momentum residual is approximately second order. Density and
total-energy residuals do not satisfy a uniform second-order first-layer
criterion, even though their fixed-time solution errors do.

### GP and Surface Reconstruction

One-step GP reconstruction orders were:

| Quantity | observed L2 order |
|---|---:|
| GP density | 3.01 |
| GP u | 2.56 |
| GP v | 2.35 |
| GP pressure | 2.92 |
| GP temperature | 3.03 |
| BI pressure | 2.61 |
| BI heat-flux derivative | 2.13 |
| BI normal viscous traction | 2.13 |
| BI viscous-traction magnitude | 2.03 |

At fixed time, BI pressure remained high order, but derivative-based surface
quantities converged more slowly:

| Quantity | L2 order | Linf order |
|---|---:|---:|
| BI pressure | 3.20 | 3.21 |
| wall heat flux | 1.56 | 1.50 |
| normal viscous traction | 2.19 | 2.04 |
| viscous-traction magnitude | 1.37 | 1.45 |

Consequently, quantitative `q_w`, tangential shear, `C_f`, and heat-transfer
coefficients are not yet certified at second order.

Primary data:

```text
IBM/cases/validation/mms_ibm/diagnostics/
pure_gp_noslip_bic_20260727/
```

## Euler Non-Regression

The no-slip policy is compile-time excluded from Euler-slip builds.
Nevertheless, the frozen pure-GP regression cases were rerun.

- PM expansion, `N=240`: all saved fields are bitwise identical to the frozen
  Gate-0 output. The downstream angle error is `0.000858 deg`, the pressure
  ratio error is `0.0062%`, and no plateau cell lies more than 1% below the
  theoretical pressure.
- Attached compression corner, `N=120`: all saved fields are bitwise
  identical. The shock-angle error is `-0.051865 deg`, and the pressure-ratio
  error is `-0.0033%`.
- Inverse-MOC nozzle, two phases, `N=64`: all volume fields are bitwise
  identical to the frozen Gate-0 outputs. The SSPRK audit recorded 516
  forward-Euler brackets with `theta=1`, and no GP fallback.
- Mach-4 circle, 12-step comparison: initial fields are bitwise identical.
  The current dirty-source build differs from the older frozen executable by
  at most about `6e-7` relative in momentum after 12 steps, while density and
  pressure remain at roundoff-level agreement. This is a small executable
  provenance drift, not a bitwise circle pass. No long circle run was repeated
  to avoid unnecessary local resource use.

Regression outputs:

```text
IBM/cases/validation/canonical/expansion_corner_2d/runs/
noslip_pressure_compat_nonregression_20260727/

IBM/cases/validation/canonical/inverse_moc_nozzle_2d/runs/
noslip_pressure_compat_nonregression_20260727/

IBM/cases/ibm_tests/2d_bvh_gpu/runs/
noslip_pressure_compat_nonregression_20260727/
```

## Qualification Status

### Passed

- pure shared-GP/full-Cartesian method identity;
- stationary isothermal no-slip pressure compatibility in 2-D Cartesian flow;
- quadratic selected-support retry with a condition-number guard;
- fixed-time complete-state convergence of order two or better;
- BI pressure convergence of order two or better;
- no observed SSPRK-stage admissibility failure in the smooth MMS;
- PM, compression-corner, and inverse-MOC non-regression.

### Not Certified

- uniform second-order first-layer semi-discrete density and energy residuals;
- fixed-time second-order wall heat flux and tangential shear magnitude;
- adiabatic or moving no-slip walls;
- body-force or nonzero wall-normal manufactured momentum sources;
- 3-D faceted geometry and feature edges;
- RZ no-slip GP-IBM;
- PelePhysics or reacting/multicomponent transport;
- shock/no-slip interaction;
- AMR viscous/thermal convergence;
- turbulent-wall quantities.

The appropriate current claim is therefore:

```text
The pure GP-IBM stationary isothermal no-slip pressure closure is a verified
two-dimensional smooth-flow candidate. Full quantitative Navier-Stokes wall
heat transfer and skin-friction production capability is not yet certified.
```
