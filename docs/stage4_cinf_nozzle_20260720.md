# Stage 4 C-infinity Cartesian-Nozzle Verification

Date: 2026-07-20

## Decision

\[
\boxed{
\text{Stage 4A all-subsonic integration/grid-phase non-regression: PASS;}
\quad
\text{steady quantitative certification under the present protocol: FAIL}
}
\]

\[
\boxed{
\text{Stage 4B smooth-wall recompression/weak-shock regression: PASS;}
\quad
\text{shock-free verification-fixture designation: INVALID}
}
\]

No IBM geometry, boundary-intercept reconstruction, WENO weight, ghost-cell
layer, or positivity-limiter threshold was changed during this campaign.
Stage 3A and Stage 3B remain frozen at their previously reported status.

## Certified scope used by the test

- two-dimensional Cartesian, uniform level-0 mesh;
- source-free calorically perfect-gas Euler equations;
- fixed smooth analytic wall with a faceted topology carrier;
- shared-GP BI-constrained quadratic reconstruction, one ghost-cell layer;
- cell-average-aware normal momentum and one-sided entropy extension;
- analytic curvature-compatible Euler-slip pressure closure;
- LLF-WENO-Z5 and SSPRK(4,3);
- opt-in conservative shared-face positivity limiter;
- Stage-3A total-condition inlet and characteristic pressure outlet;
- `max_grid_size=256`, `blocking_factor=64`.

The physical wall is

\[
H(x)=H_0-\Delta H\,b((x-x_t)/a),
\qquad
b(\xi)=
\begin{cases}
\exp(1-(1-\xi^2)^{-1}),&|\xi|<1,\\
0,&|\xi|\ge1,
\end{cases}
\]

with exact straight inlet and outlet buffers.  The same wall files are used
on every grid.  The topology-carrier closure edges lie outside the physical
domain; their grid-dependent internal subdivision does not change the in-domain
wall contour.

## All-subsonic matrix

The common-time matrix contains four grids,
`192x64`, `384x128`, `576x192`, and `768x256`, and two fixed normalized
wall-normal phases, 0 and `0.37h`.  All eight calculations reached `t=12`.
Every active state remained admissible, the limiter remained at `theta=1`,
and every shared GP used the quadratic constrained reconstruction.

On the finest axis-aligned grid:

| Diagnostic | Value |
|---|---:|
| minimum density / pressure | 0.883187 / 0.840376 |
| section mass-flow relative range | 2.249e-4 |
| total-enthalpy relative range | 2.745e-5 |
| entropy-proxy range | 9.181e-9 |
| BI pressure symmetry relative Linf | 8.466e-15 |
| wall-normal velocity relative L2 / Linf | 1.308e-5 / 1.411e-4 |
| raw/conservative force difference normalized by \(\int|pn_x|ds\) | 7.616e-4 |

The wall-normal phase sensitivity decreases under refinement:

| \(N_y\) | section-pressure relative L2 | BI-pressure relative L2 |
|---:|---:|---:|
| 64 | 2.060e-5 | 1.595e-4 |
| 128 | 6.738e-6 | 5.780e-5 |
| 192 | 3.189e-6 | 2.261e-5 |
| 256 | 1.406e-6 | 1.391e-5 |

This establishes controlled grid-phase sensitivity and a stable low-amplitude
all-subsonic integration path.  It does not establish a formal nozzle order:
the BI and section-pressure Richardson error vectors are not aligned across successive
triplets, so no GCI is accepted.  The wall-normal trace has fitted zero-limit
orders 1.72 and 1.39 for the two phases, with oscillatory adjacent-grid orders.

The `t=40` stationarity discriminator fails the present steady acceptance
protocol.  The pressure-field change over the final two time units is
`9.93e-4` at `Ny=64` and `1.04e-3` at `Ny=128`, so it does not decrease with
spatial refinement.  The follow-on mode discrimination shows that this is not
a single mean-pressure drift:

- a streamwise round-trip mode follows the quasi-one-dimensional acoustic
  travel time, changes monotonically with boundary-plane spacing, is unchanged
  when CFL is reduced from 0.30 to 0.15, and is strongly dependent on the
  initial condition;
- a symmetric transverse pressure mode remains in the full field after the
  section-mean mode has decayed.  High-cadence probes give frequencies
  `0.819/0.895/0.819` at `x=1.5/3.0/4.5`, compared with local rigid-wall cutoff
  estimates `0.899/0.962/0.899`.

The latter is consistent with a cut-off or trapped duct-acoustic mode, but
that interpretation still requires an independent reference.  It spans the
core and near-wall region and is symmetric to roundoff, so it is not identified
as a wall-localized IBM defect.  At `t=80`, the quasi-one-dimensional-start
axial-mode amplitude has decayed by a factor of about 50, while the final
two-time-unit full-field pressure change remains `8.65e-4`.  The section mean
is therefore not a sufficient stationarity metric.

The strict source-free discrete budget is now closed.  It integrates the
actual finalized numerical flux from all four SSPRK(4,3) Forward-Euler
brackets with weights `dt/6, dt/6, dt/6, dt/2`.  Across the five mode-matrix
runs, the largest absolute mass, momentum, or total-energy closure residual is
`3.93e-12`.  Single-FAB, three-FAB, and two-rank MPI decompositions retain
roundoff-level closure; MPI and serial state integrals differ by no more than
`1.8e-15`.  This rules out RK-flux bookkeeping and the IBM interface reaction
as the source of the pressure unsteadiness.  The older reconstructed-section/
VTP momentum balance remains a separate, non-strict diagnostic.

## Choked, supersonic-exit discriminator

The smooth-wall recompression/weak-shock contour was tested at `Ny=64` and
`Ny=128` through `t=40`.
Both trajectories remain admissible; no limiter event or shared-GP linear
fallback occurs.  The downstream and exit cells remain supersonic.

| Diagnostic | \(N_y=64\) | \(N_y=128\) | observed order |
|---|---:|---:|---:|
| max adjacent \(|\Delta p|/\bar p\) | 8.714e-2 | 7.189e-2 | 0.28 |
| max adjacent \(|\Delta\sigma|\) | 5.030e-4 | 2.974e-4 | 0.76 |
| max \(h|\nabla p|/p\) | 7.758e-2 | 6.907e-2 | 0.17 |
| pressure second-difference sensor | 1.085e-2 | 8.110e-3 | 0.42 |
| section mass-flow relative range | 6.951e-4 | 3.171e-4 | -- |
| total-enthalpy relative range | 1.136e-4 | 5.767e-5 | -- |

These rates are incompatible with a globally smooth supersonic solution.
Nearly constant `h|grad p|/p` and the sharpening X-shaped bulk wave system are
consistent with an `O(h)`-thick weak shock/recompression layer.  A generic
smooth converging-diverging contour is not automatically shock-free: returning
the downstream wall angle to zero generates compression characteristics.

This contour is therefore retained only as a smooth-wall recompression/
weak-shock robustness regression.  Its failure as a shock-free verification
fixture is not a solver failure.  A formal shock-free choked order gate requires a
characteristic-designed contour, such as a verified method-of-characteristics
nozzle, rather than further tuning of the IBM closure.

## Reproducibility

- Cerisse commit: `dc35646d8eb90aaa362a3292a9d702e30c9d5af5` plus recorded dirty changes;
- fixed-time/spatial-matrix source-tree SHA-256:
  `f65e000a045b2f4c94e8ef1baa11acbadd2d1e94db9a09c064b7f28e409d7614`;
- strict-budget/mode-matrix source-tree SHA-256:
  `b59f508bf15a11001c4b4e568135606beb85470f4819605478d73b09876d6325`;
- AMReX: `bd922c6216e0a734f3b1cf0ca73e7d669b90f3ef`;
- upper/lower wall SHA-256: `7028f48fd3b4253c1663fb939771ded5819819f6f331e6d2f90feaac28faf007` /
  `f6279653dfaa2752069c74819fa1a75bf0dc929acf9090c403f65c51b220a964`;
- strict-budget CPU/MPI/CUDA executable SHA-256:
  `9f19b9efbb0e98344cf53154222042bb85be6fddeaf580a842f563f6df747769` /
  `47925fc0256439f623d46ddc505ad41d7467ea447bbe625049858a7c33e673b7` /
  `47c1bcd9d11003e562e4de9e8fdbc04253b8877fefe4d251ffcd5a5069d5bf04`.

The earlier `192x64`, `t=0.05` CPU/CUDA comparison has a maximum absolute field
difference of `1.19904e-14`; IBM markers agree exactly.  The latest strict-
budget CUDA source compiles and links for `sm_89`, but its runtime retest is
blocked before solver initialization by local WSL CUDA error 35 after a host
driver update.  This is recorded as an environment block, not a numerical
pass.  No SU2 or other external CFD solver was used.

Raw evidence is under
`IBM/cases/validation/canonical/cinf_nozzle_2d/runs/`.  The formal matrix,
stationarity, acoustic-mode, and choked reports are respectively:

- `runs/stage4_fixed_t12/subsonic/summary/REPORT.md`;
- `STAGE4_SUBSONIC_STATIONARITY.md`;
- `runs/stage4_low_frequency/summary/REPORT.md`;
- `STAGE4_CHOKED_PILOT.md`.

## Remaining Stage-4 gates

1. Define a steady-reference protocol that either initializes close to the
   target two-dimensional branch or time-averages complete acoustic cycles
   with sampling uncertainty below the spatial-discretization uncertainty.
   Mean-pressure relaxation is not justified solely by the transverse mode.
2. Replace the invalid shock-free contour by a continuous-curvature planar
   inverse-MOC contour.
3. Supply an independently grid-converged, same-geometry two-dimensional
   reference for wall pressure and transverse flow after the contour and
   steady solution branch are fixed.
4. Only after the smooth gates pass, add a prescribed-back-pressure internal
   shock with shock-position and outlet-location sensitivity.

Until these gates close, Stage 4A is an integration/grid-phase non-regression
success but fails steady quantitative certification under the current
protocol.  Stage 4B remains a useful weak-shock regression, not a shock-free
order-verification case.
