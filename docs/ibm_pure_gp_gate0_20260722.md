# Pure GP-IBM Gate 0 Current-Source Report

Status date: 2026-07-22

## Decision

Gate 0 is **PASS** for the frozen pure GP-IBM method:

\[
\boxed{\text{Gate 0: current-source software non-regression PASS}}
\]

This decision closes the current-source recovery gate. It is a software,
stability, and canonical-physics regression decision. It is not a claim of
absolute Mach-4 circle accuracy, a true-curved-domain conservation proof, or
a completed grid-convergence certification.

The passing evidence consists of:

1. the long-time Mach-4 circle flow and positivity gate;
2. a Mach-2 Prandtl-Meyer expansion corner;
3. a Mach-2 attached compression corner;
4. a source-free continuous-curvature inverse-MOC nozzle at two normalized
   Cartesian grid phases.

No IBM formula was changed while executing these four sub-gates.

## Frozen numerical method

Every sub-gate used the same pure-method family:

- one shared state at each real solid-side ghost-cell centre;
- cell-average-aware BI-constrained weighted least squares (BI-CWLS);
- curvature-compatible normal-pressure closure on smooth Euler-slip walls;
- one-sided fluid entropy extension and ideal-gas EOS closure;
- one configured solid-side ghost-cell layer;
- full Cartesian cell volumes and full Cartesian face areas;
- LLF flux splitting with WENO-Z5 reconstruction;
- conservative shared-face positivity limiting;
- SSPRK(4,3), with each forward-Euler bracket subject to the admissibility
  contract;
- volume extrapolation order one and BI surface recovery order two where
  surface pressure is a gate observable.

The runtime manifest reported
`family=pure_shared_gp_full_cartesian`, `pure_gp_production=1`, and
`cut_control_compiled=0`. The inverse-MOC executable also had no defined
cut-control or aggregate-update symbol. Volume fractions, area fractions,
physical wall-segment replacement fluxes, aggregate states, cut-control
residuals, and redistribution did not enter the RHS.

## Sub-gate results

| Sub-gate | Result | Main quantitative evidence |
|---|---|---|
| Mach-4 circle, `320 x 320`, `D/h=64` | PASS | detached bow shock through `t*=55.56`; `Delta/R=0.5496841171`; final relative plateau span `9.14e-8`; 74 limiter brackets confined to `0.531<=t*<=0.663`; no global/all-low fallback or rejected step |
| Mach-2 PM expansion, `240 x 120` | PASS | flow-angle error `+0.000858 deg`; Mach error below `0.001%`; pressure-ratio error `+0.006%`; no plateau cell more than 1% below theory |
| Mach-2 compression corner, `120 x 120` | PASS | shock-angle error `-0.051865 deg`; flow-angle error `+0.003461 deg`; pressure-ratio error `-0.003%`; all 42 shock-fit columns accepted |
| Inverse-MOC, `256 x 64`, two phases | PASS for Gate 0 | fluid pressure relative L2 `2.11e-4` and `2.13e-4`; BI `Cp` L1 at most `1.33e-4`; BI pressure-force relative error `3.64e-5` and `1.42e-5`; force phase spread `2.22e-5` |

The circle result and its strict limitations are documented separately in
[the circle gate report](ibm_pure_gp_circle_gate_20260722.md).

### Prandtl-Meyer expansion

The numerical downstream state was compared with the analytic
Prandtl-Meyer turning relation:

| Quantity | Theory | Numerical | Error |
|---|---:|---:|---:|
| Flow angle [deg] | -20.000000 | -19.999142 | +0.000858 deg |
| Mach number | 2.830595 | 2.830582 | below 0.001% |
| `p2/p1` | 0.275178 | 0.275195 | +0.006% |
| `rho2/rho1` | 0.397854 | 0.397849 | -0.001% |

The post-fan pressure relative L1/L2 errors were `0.038%/0.095%`. Its
minimum departure from the theoretical plateau was `-0.507%`; no active
fluid cell in the plateau was more than 1% below the theoretical pressure.
No positivity-limiter event was recorded.

Artifacts:

- [machine-readable result](../IBM/cases/validation/canonical/expansion_corner_2d/runs/gate0_current_20260722/pm_N240/result.json)
- [field figure](../IBM/cases/validation/canonical/expansion_corner_2d/runs/gate0_current_20260722/pm_N240/pm_fields.png)
- [execution log](../IBM/cases/validation/canonical/expansion_corner_2d/runs/gate0_current_20260722/pm_N240/run.log)

### Attached compression corner

The weak attached branch of the theta-beta-M relation and the normal-shock
Rankine-Hugoniot relations supplied the reference:

| Quantity | Theory | Numerical | Error |
|---|---:|---:|---:|
| Shock angle [deg] | 53.422941 | 53.371075 | -0.051865 deg |
| Flow angle [deg] | 20.000000 | 20.003461 | +0.003461 deg |
| Downstream Mach | 1.210218 | 1.210884 | +0.055% |
| `p2/p1` | 2.842863 | 2.842769 | -0.003% |
| `rho2/rho1` | 2.042006 | 2.042604 | +0.029% |
| `p02/p01` | 0.892914 | 0.894277 | +0.153% |

No positivity-limiter event was recorded. This is a single-grid canonical
regression; it does not replace a shock-location grid-convergence study.

Artifacts:

- [machine-readable result](../IBM/cases/validation/canonical/expansion_corner_2d/runs/gate0_current_20260722/compression_N120/result.json)
- [shock figure](../IBM/cases/validation/canonical/expansion_corner_2d/runs/gate0_current_20260722/compression_N120/oblique_shock.png)
- [execution log](../IBM/cases/validation/canonical/expansion_corner_2d/runs/gate0_current_20260722/compression_N120/run.log)

## Inverse-MOC two-phase audit

The inverse-MOC field is a source-free two-dimensional Euler solution, not a
quasi-one-dimensional area-Mach reference. Both runs used the same executable,
`256 x 64` grid, physical geometry, boundary model, and final time `t=0.5`.
Only the Cartesian domain origin was shifted by `(0.37 h_x, 0.23 h_y)`.

All 512 ghost points in both phases used the quadratic BI-CWLS path; there was
no linear fallback. The two-step verbose preflight recorded eight of eight
SSPRK forward-Euler brackets as `HIGH_ORDER`, `theta_min=1`, and every
post-stage audit passed. The formal runs recorded no limiter activation.

| Metric | Phase `(0,0)` | Phase `(0.37,0.23)` |
|---|---:|---:|
| Fluid pressure relative L2 | `2.1127e-4` | `2.1309e-4` |
| First wall-neighbour band pressure relative L2 | `2.7397e-4` | `2.5616e-4` |
| Maximum upper/lower BI `Cp` L1 | `1.3307e-4` | `1.1461e-4` |
| Maximum non-endpoint BI `Cp` Linf | `1.2451e-3` | `1.5892e-3` |
| BI surface-e2 pressure-force relative error | `3.6401e-5` | `1.4202e-5` |
| Physical-section mass-flow relative mean error | `1.9258e-5` | `3.5528e-5` |
| Physical-section mass-flow relative spread | `4.4789e-4` | `4.5962e-4` |
| Reconstructed BI signed mass / reference mass flow | `-4.3899e-6` | `-7.3093e-5` |
| Reconstructed BI absolute mass / reference mass flow | `9.9059e-4` | `7.9858e-4` |
| Active-fluid mass-budget residual / reference mass flow | `-8.48e-14` | `-8.15e-14` |

The BI pressure-force phase spread is `2.2199e-5` relative to the exact MOC
force. The physical-section mass flow is evaluated between the analytic walls
with a 12-point Gauss rule; smooth cell averages are interpreted as point
samples with an `O(h^2)` equivalence error. A raw full-cell staircase sum has
a 5.9--6.5% axial spread and is retained in the JSON only as a
geometry-quadrature diagnostic, not as the physical section mass flow.

The reconstructed BI mass figures are finite-grid impermeability diagnostics.
Only one grid was rerun for Gate 0, so their convergence is not certified by
this report. That remains a Gate 1/Gate 2 requirement.

Artifacts:

- [two-phase analysis](../IBM/cases/validation/canonical/inverse_moc_nozzle_2d/runs/gate0_pure_20260722/analysis_two_phase.json)
- [zero-phase log](../IBM/cases/validation/canonical/inverse_moc_nozzle_2d/runs/gate0_pure_20260722/phase000/N0064/run.log)
- [shifted-phase log](../IBM/cases/validation/canonical/inverse_moc_nozzle_2d/runs/gate0_pure_20260722/phase037_023/N0064/run.log)

## Force and conservation interpretation

For this pure full-cell method:

- `surface-e2` BI pressure quadrature is the physical pressure-load
  observation;
- `production_force` and `cv_conservative_force` are the Cartesian
  staircase-domain momentum exchange;
- equality of the latter two proves discrete bookkeeping consistency but not
  equality with the true curved-wall pressure load.

In the two inverse-MOC phases, the Cartesian production/CV forces agree with
each other at roundoff but differ from the MOC pressure force by 7.29% and
7.44%. These values are reported, not hidden, and are excluded from the pure
GP Gate 0 physical-load criterion. True-curved-domain machine conservation is
not claimed.

## Provenance

All runs used Git commit
`dc35646d8eb90aaa362a3292a9d702e30c9d5af5` and tracked dirty-diff SHA-256
`604c6323ebf00368df42bbb9bdcf9da39441070d8af41c8a483709e39c4f5a64`.

The source-content IDs are intentionally case-specific: each hashes all
solver files under `src/` together with that case's compile-time `prob.h`,
`GNUmakefile`, and, where applicable, reference header. Therefore the circle
ID `5b45451c364d120e`, ramp ID `ad0e0c9c747a119b`, and inverse-MOC ID
`84db3efb00517a42` identify three builds of the same working tree, not three
different Git revisions.

| Artifact | SHA-256 |
|---|---|
| PM executable | `77081e00ea2f735416042078d1124666733996f4626b687459a38469a3665252` |
| Compression executable | `9a5b0b11e26cee56ef7137ea3477134eee8d496ede6fa6a699628a94499818ff` |
| Inverse-MOC executable | `cfcb309b66914c8bf2012f23cd59eb1997deffaf5f7454157a2f4d790d5a7020` |
| PM input / geometry / result | `d0e9b7fdcdf47430d87dbf228b14a6ffd0b208497d0b98d41cc5e221808b5ab8` / `7fc7723871ad1dd240d017c9697eed539cf59ed66c8fb3023c44ff29047c51ce` / `f3c13385ebc11907542b869c16a330db4320c77dba43d3504ae8008da92a5b8e` |
| Compression input / geometry / result | `4be4d8c4bb3a19624f224ecb03566921bb4e6c39d43695c7bc2d475aa6bc126e` / `f73412036f0e0d5d3f82a40c238f5c7cc9a3e3acf4d8201374960a0da45bf32c` / `b757e21629d4e7751f1f93c312fd0fbcab2813ebbf6f77f2750c6c903a378679` |
| Inverse-MOC input / upper geometry / lower geometry | `66016f9e58d6774be708f8a7fe98563defb31c069fbf825824489c5226864c0c` / `90f47995d00597f61740fe31cd0ab9eea7941c9f5acb820f9d8a37001af834ce` / `bf2c1666f841279df46e86050e17a4f27213b8843119e9aa7c115251204e7008` |
| Inverse-MOC reference metadata / table | `5568efbee6f00f97675b5d9b521c146babf96fc43be000619180f6aa6644460e` / `0e9b5718cbd28bc37a805147edb1c1d9237e718e1b8fb0769c9c9a32ec23b299` |
| Two-phase analysis JSON | `f38141a79a04593493ce06096d070641d9a12d94ed20166d0c04ed00171bca58` |

## What Gate 0 does not certify

Gate 0 does not establish:

- a converged absolute Mach-4 circle stand-off distance;
- true curved-domain mass or momentum conservation;
- grid-converged BI impermeability;
- a smooth-nozzle production certificate;
- a shocked-nozzle production certificate;
- AMR, R-Z, NS heat transfer/skin friction, general 3-D STL, moving bodies,
  or FSI capability.

The next permitted work is Gate 1 smooth exact-solution co-refinement and
Gate 2 smooth inverse-MOC/nozzle qualification. Gate 0 permits exploratory
smooth-nozzle calculations; it does not by itself authorize quantitative
production use.
