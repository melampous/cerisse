# Cartesian Euler/GP-IBM regression campaign (2026-07-19, audited 2026-07-20)

## Scope

This campaign exercises the frozen two-dimensional Cartesian production
candidate without changing its numerical method:

- calorically perfect, single-component gas with `gamma=1.4`;
- source-free Euler equations;
- characteristic LLF flux splitting and WENO-Z5 reconstruction;
- SSPRK(4,3), with the conservative per-face positivity limiter enabled on
  every Forward-Euler bracket;
- uniform level-0 grids, no AMR or reflux;
- fixed shared-GP IBM for the solid-wall cases;
- BI-constrained degree-two reconstruction, visibility mode 2, expanded
  support, one solid-side ghost-cell layer, `iorder=2`, `eorder=1`, and
  `alpha=0.6`;
- analytic curvature-compatible pressure/entropy extension for the smooth
  circular wall, and the straight-facet zero-curvature (`pzero`) extension for
  the sharp-corner cases.

The source tree is based on Git commit
`dc35646d8eb90aaa362a3292a9d702e30c9d5af5` with the recorded dirty working
tree. CUDA builds target `sm_89` and ran on the local RTX 4090 Laptop GPU.

The tests deliberately separate four questions:

1. Does the limiter preserve the spatial accuracy of a smooth source-free
   exact solution?
2. Does the GP-IBM reproduce the Prandtl-Meyer and oblique-shock relations?
3. Does a sequential two-shock interaction retain the correct intermediate
   states and wave topology?
4. Does the fully discrete no-solid method remain admissible in a canonical
   near-vacuum double-rarefaction problem?

## 1. Source-free stationary isentropic vortex around a circular IBM wall

The exact isentropic vortex is stationary and tangent to a concentric circle,
so it satisfies the homogeneous Euler equations and the exact slip condition
without a manufactured source. The curved-wall pressure condition is

\[
\partial_n p=\rho u_t^2/R.
\]

Five grids, `N=80,120,160,240,320`, were advanced to `t=0.01`. The time step
was reduced approximately as \(h^{4/3}\), so the third-order temporal error is
nominally \(O(h^4)\). Exact finite-volume cell averages were evaluated with
four-point Gaussian quadrature.

### Fixed-time L2 convergence

| Region | rho | u | v | p | T |
|---|---:|---:|---:|---:|---:|
| all active fluid | 2.089 | 2.090 | 2.090 | 2.019 | 2.005 |
| first fluid-cell layer | 2.024 | 2.068 | 2.068 | 2.017 | 1.980 |

The first-layer Linf fitted orders are `1.94, 1.94, 1.94, 1.94, 1.95` for
`rho,u,v,p,T`, respectively. Across the five grids, all 648 SSPRK
Forward-Euler brackets remained high order with `theta=1`; the limiter was
asymptotically inactive in the strongest possible sense. Relative mass drift
decreased from `6.47e-7` to `3.02e-8`, and relative energy drift from
`8.82e-7` to `4.19e-8`.

**Decision: PASS.** The opt-in limiter does not degrade the approximately
second-order smooth curved-wall solution or first-layer accuracy under this
protocol.

Artifacts:

- `IBM/cases/validation/mms_ibm/diagnostics/production_freeze_20260719/source_free_exact_limiter/convergence_L2.md`
- `IBM/cases/validation/mms_ibm/diagnostics/production_freeze_20260719/source_free_exact_limiter/convergence_Linf.md`
- `IBM/cases/validation/mms_ibm/diagnostics/production_freeze_20260719/source_free_exact_limiter/invariants.md`

## 2. Mach-2 Prandtl-Meyer expansion corner

The `240 x 288` calculation was advanced to `t=7`. The exact downstream state
is obtained from the Prandtl-Meyer function and the isentropic relations.

| Quantity | Exact | Numerical | Error |
|---|---:|---:|---:|
| flow angle [deg] | -20.000000 | -20.005117 | -0.005117 deg |
| Mach number | 2.830595 | 2.830869 | +0.010% |
| `p2/p1` | 0.275178 | 0.275073 | -0.038% |
| `rho2/rho1` | 0.397854 | 0.397724 | -0.033% |

The downstream pressure-plateau L2 error is `0.573%`. Its minimum is `3.293%`
below the exact plateau; `6.062%` of plateau cells are more than 1% low, but no
cell is more than 5% low. The global active-fluid pressure minimum is local to
the corner. All 8172 Forward-Euler brackets remained admissible without a
limiter activation.

**Decision: PASS for PM turning and downstream thermodynamics.** The localized
corner pressure defect is retained as a non-oscillatory-accuracy diagnostic;
this single-grid test is not a formal corner convergence proof.

Artifacts:

- `IBM/cases/validation/canonical/expansion_corner_2d/runs/production_freeze_20260719/pm_bic_limiter_N240/result.md`
- `IBM/cases/validation/canonical/expansion_corner_2d/runs/production_freeze_20260719/pm_bic_limiter_N240/pm_fields.png`

## 3. Mach-2 single compression corner

The `120 x 144` calculation was advanced to `t=14`. The weak attached branch
of the theta-beta-M relation and the Rankine-Hugoniot jump conditions provide
the reference.

| Quantity | Exact | Numerical | Error |
|---|---:|---:|---:|
| shock angle beta [deg] | 53.422941 | 53.479231 | +0.056291 deg |
| downstream flow angle [deg] | 20.000000 | 20.011380 | +0.011380 deg |
| downstream Mach number | 1.210218 | 1.210508 | +0.024% |
| `p2/p1` | 2.842863 | 2.842421 | -0.016% |
| `rho2/rho1` | 2.042006 | 2.042011 | below 0.001% |
| `T2/T1` | 1.392191 | 1.392017 | -0.012% |
| `pt2/pt1` | 0.892914 | 0.893407 | +0.055% |

The fitted shock has a normal L1 location error of `0.05653`, or `0.68h`.
Outside a fixed shock strip, the normalized pressure L2 error is `0.596%`.
The run completed 1988 time steps (7952 Forward-Euler brackets). No local,
global-theta, or all-low-order fallback was recorded; every post-SSP audit was
therefore admissible.

**Decision: PASS.** Shock angle, flow turning, all principal jump ratios, and
the attached-shock topology meet the engineering regression gate.

Artifacts:

- `IBM/cases/validation/canonical/expansion_corner_2d/runs/production_freeze_20260719/single_compression_bic_limiter_N120/result.md`
- `IBM/cases/validation/canonical/expansion_corner_2d/runs/production_freeze_20260719/single_compression_bic_limiter_N120/fields.png`

## 4. Sequential double compression

The `120 x 144` wall angles are `0 -> 10 -> 20 deg`, with corners at `x=2`
and `x=4`. The second theta-beta-M solve uses the state behind the first shock
and the incremental second turn.

| Quantity | Exact | Numerical | Error |
|---|---:|---:|---:|
| shock 1 global angle [deg] | 39.313932 | 40.511902 | +1.197970 deg |
| shock 2 global angle [deg] | 59.384042 | 59.531757 | +0.147715 deg |
| state-1 Mach number | 1.640522 | 1.635693 | -0.294% |
| state-1 flow angle [deg] | 10.000000 | 9.986272 | -0.0137 deg |
| state-1 `p/p0` | 1.706579 | 1.705555 | -0.060% |
| state-2 Mach number | 1.284889 | 1.280776 | -0.320% |
| state-2 flow angle [deg] | 20.000000 | 19.986663 | -0.0133 deg |
| state-2 `p/p0` | 2.803191 | 2.801400 | -0.064% |

The fitted incident rays intersect `0.3342` length units, approximately
`4.0h`, from the independent-ray prediction. The first shock has only two
length units in which to establish before the second corner, so its fitted
angle and the extrapolated intersection remain grid-sensitive. All 6824
Forward-Euler brackets remained high order and admissible.

**Decision: CONDITIONAL PASS.** Sequential thermodynamic states, wall turning,
and the two-shock topology pass. The first-ray angle and interaction location
are not yet grid-converged and are not certified as quantitative observables.

Artifacts:

- `IBM/cases/validation/canonical/expansion_corner_2d/runs/production_freeze_20260719/double_compression_bic_limiter_N120/result.md`
- `IBM/cases/validation/canonical/expansion_corner_2d/runs/production_freeze_20260719/double_compression_bic_limiter_N120/fields.png`

## 5. No-solid near-vacuum double rarefaction

The canonical source-free Riemann data are

\[
(\rho,u,p)_L=(1,-2,0.4),\qquad
(\rho,u,p)_R=(1,2,0.4).
\]

The exact non-vacuum star state is

\[
p_*=1.89387342\times10^{-3},\qquad
\rho_*=2.18521182\times10^{-2},\qquad u_*=0.
\]

The loaded polygon is wholly outside the domain; no ghost cell or immersed
boundary face exists in either run. Thus every defect in this subsection is a
bulk LLF-WENO-Z5/SSPRK result, not an IBM wall error.

| Grid | Steps | min rho | min p | min internal-energy density | Limiter events |
|---|---:|---:|---:|---:|---:|
| `512 x 64` | 1056 | 1.13254e-2 | 2.67196e-3 | 6.67990e-3 | 0/4224 brackets |
| `1024 x 128` | 2111 | 9.85077e-3 | 2.19897e-3 | 5.49742e-3 | 0/8444 brackets |

Both runs contain no non-finite values, retain exact transverse uniformity,
and have pressure symmetry defects below `5e-12`. Exact conservative
cell-average integration closes the final finite-domain mass and energy
inventories to approximately `1e-14` on the `512 x 64` run.

The exact finite-volume oracle changes the interpretation of the error. It is
not confined to one or two symmetry-plane cells. Complete-star-region
relative errors are:

| grid | variable | relative L1 | relative L2 | relative Linf |
|---|---|---:|---:|---:|
| `512 x 64` | rho | 21.99% | 27.14% | 48.17% |
| `512 x 64` | p | 46.84% | 47.21% | 60.59% |
| `512 x 64` | T | 90.41% | 100.24% | 172.22% |
| `512 x 64` | sigma | 77.91% | 85.38% | 137.99% |
| `1024 x 128` | rho | 20.55% | 26.58% | 54.92% |
| `1024 x 128` | p | 22.16% | 23.05% | 37.80% |
| `1024 x 128` | T | 54.38% | 68.84% | 157.57% |
| `1024 x 128` | sigma | 52.73% | 65.68% | 138.04% |

At `512 x 64`, density error above 10% occupies physical length `0.07031`
of the `0.10156` complete star region; pressure error above 10% occupies the
whole region. These measures remain `0.06836` and `0.10352` after doubling
the resolution. This is a star-region thermodynamic error, with the largest
error at the symmetry plane, rather than a purely localized center-cell
defect.

The first SSPRK Forward-Euler bracket identifies the startup mechanism. The
center-cell entropy proxy after the numerical update exceeds the exact
finite-volume value by `0.0812463`, while all limiter coefficients remain
one. A controlled run initialized from exact conservative cell averages at
physical solution age `tau=0.01`, then advanced to the same physical age
`0.15`, reduces the `512 x 64` star-region relative L1 errors in
`rho,p,T,sigma` to `0.297%,0.359%,0.459%,0.566%`; its first-bracket entropy
error is only `1.58e-8`.

**Decision: admissibility/conservation PASS; limiter activation N/A;
ab-initio star-region thermodynamic accuracy FAIL.** The t=0 versus exact-
evolved-start comparison identifies the classical bulk receding-flow
overheating mechanism. It does not implicate the IBM or the positivity
limiter and does not justify changing either method.

Artifacts:

- `IBM/cases/validation/canonical/near_vacuum_riemann_2d/runs/production_freeze_20260719/near_vacuum_N512/result.md`
- `IBM/cases/validation/canonical/near_vacuum_riemann_2d/runs/production_freeze_20260719/near_vacuum_N512/profiles.png`
- `IBM/cases/validation/canonical/near_vacuum_riemann_2d/runs/production_freeze_20260719/near_vacuum_N1024/result.md`
- `IBM/cases/validation/canonical/near_vacuum_riemann_2d/runs/production_freeze_20260719/near_vacuum_N1024/profiles.png`
- `IBM/cases/validation/canonical/near_vacuum_riemann_2d/NEAR_VACUUM_OVERHEATING_DIAGNOSIS_20260720.md`
- `IBM/cases/validation/canonical/near_vacuum_riemann_2d/runs/nearvac_diagnostics_20260720/tau001_N512/result.json`

## Provenance

| Executable | SHA-256 |
|---|---|
| source-free circular-wall exact solution | `148f982a1139e03c759731f5f9b827d0e8345042a29851574b334d958d1d670b` |
| PM expansion | `ec7f409e5ad1be9ec71cb1f7e8e411b663edfddd262f6fae0b386f527fece923` |
| single compression | `92fee3e13aef1484519ca87ca05aa2ff13bd7be7cbfc37e6dc69c3bc87e9efec` |
| double compression | `8b38fe5cbab9e5379cc696ba1ab748f53e100a05077e01095425632564223593` |
| no-solid near-vacuum Riemann | `b9aa59d6fb0ab42e475875e88eb74bdeca89c08f59c4b27a8919c501e37f88a5` |

Important geometry hashes are:

- 4096-edge circle: `2fea8e7fe4b0a9911b5c1d640df730a685042816289fda38e9bca6f16f23cfe2`;
- PM ramp: `7fc7723871ad1dd240d017c9697eed539cf59ed66c8fb3023c44ff29047c51ce`;
- single-compression ramp: `f73412036f0e0d5d3f82a40c238f5c7cc9a3e3acf4d8201374960a0da45bf32c`;
- double-compression ramp: `cb169d5131df4176fb40ee76e917551adc4c5ca02194d3f5a965c149fe2a44fa`.

Input hashes are `d0e9b7fdcdf47430d87dbf228b14a6ffd0b208497d0b98d41cc5e221808b5ab8`
(PM), `4be4d8c4bb3a19624f224ecb03566921bb4e6c39d43695c7bc2d475aa6bc126e`
(single compression),
`57159dca7b015412f944bd43eb6952a3e37d3108a039383674be509eb280242b`
(double compression), and
`0f6904ed8251ca3eeacc9e0ec0f8989af83fee03c88d66c869ff693690707a93`
(near vacuum). The corresponding corner analyzers have SHA-256 values
`ae59aa38301ddf4c3605c183e55a994b94d84022e85656876828e4adf60d7741`,
`55b547f1ab2942a6670db9ea198af8075dbdc2d9e21b287bb96f9b285be3712a`,
and `85f42f4b43061655ca0afe49538e78af142dfdf926dd7d3040be736c9e02a3ac`.
The near-vacuum analyzer hash is
`4fb6cd71146de24be023a44e2484ac24ea1a43be7bdb024e1af024bdd288d8dd`.
The exact-evolved-start diagnostic executable has SHA-256
`aff5575d411e75098302f89e7cc948f00325f6ce566d16d99277ee9fbcaf5c18`.

## Overall decision

| Gate | Decision |
|---|---|
| source-free smooth circular-wall accuracy | PASS |
| Prandtl-Meyer expansion | PASS |
| single compression corner | PASS |
| sequential double compression | CONDITIONAL PASS |
| no-solid near-vacuum admissibility | PASS |
| no-solid near-vacuum limiter activation | N/A (`theta=1`) |
| no-solid ab-initio star-region thermodynamic fidelity | FAIL |

**Stage-2 decision: CONDITIONAL GO; proceed to Stage 3 characteristic
boundary conditions.** The campaign demonstrates that the opt-in conservative limiter does not
degrade the smooth curved-wall second-order result and that the frozen GP-IBM
retains accurate PM and oblique-shock states. It does not yet close every gate
for a shocked nozzle: characteristic boundary conditions remain unverified,
the double-shock interaction location requires three-grid convergence before
SRP, and the no-solid ab-initio receding-flow overheating defect remains a
documented bulk limitation. None of these results justifies changing IBM,
BIC, ghost-cell layers, WENO weights, or limiter thresholds.

The current source-content build `75f7ee86428df61f` also reruns the frozen
Mach-4 circle through `t*=55.56` with every final plot variable bitwise
identical, identical final VTP, identical limiter statistics, and unchanged
`Delta/R=0.5496841171`. The provenance report is
`IBM/cases/ibm_tests/2d_bvh_gpu/production_freeze_20260720/circle_N320_limiter_load_src75f7ee86428df61f/REPORT.md`.
