# Short nonlinear R-Z axis gate (2026-08-08)

## Scope and reproducibility

This is an all-fluid, Euler, level-0 regression.  It changes no production
flux or update code and does not instantiate IBM, cut-cell, EB, aggregate, or
AMR flow paths.  The four executables were rebuilt after the final
`advance.cpp` diagnostic fail-closed edit; `make -q` returned success for all
four immediately before the runs.

- state: exact stationary normal shock, upstream Mach 10;
- domain: `0 <= r <= 0.25`, `0 <= z <= 1`;
- mesh: `16 x 64`, so `dr=dz=0.015625`;
- R-Z axis symmetry at `r=0`, exact state at `r=0.25`, physical axial states;
- time integration: SSPRK(4,3), CFL 0.25;
- stop: `t=0.02` (67--68 steps, before the `max_step=80` cap);
- output cadence: every 10 steps plus the final state;
- comparison: LLF-WENO-Z5, AFD-HLLC-WENO-Z5, Skew4 central, Skew4-JST;
- each scheme was run with its local first-order LLF fallback off and on.

LLF-WENO-Z5 versus AFD-HLLC-WENO-Z5 is a whole spatial-operator comparison,
not a Riemann-solver-only HLLC-versus-LLF A/B.

Source identity at execution:

| source | SHA-256 |
|---|---|
| `src/rhs/Weno.h` | `47ef843a31fb22ab93ad25f59214706376d52d1feff133adc7482006e2cd68a9` |
| `src/rhs/Afd.h` | `732bf2aec3b4843b043ae03883a3021d784df721f9b0df6cacdf06728c190f89` |
| `src/rhs/Skew.h` | `9d50e7ca27e6671fc19a14b6e243658dc64b925c35bb0910738b34f3f48898c3` |
| `src/tim/compute_rhs.cpp` | `23ef010c74d2d541437d888e258ca0975ff511b9268ca30aa82bd220b5d550a2` |
| `src/tim/advance.cpp` | `0f1043ac64ba525c596a7f6a5e23ee2cc9e0b3cd74d94ea9afa7812081bb33ef` |

The Git base was `dc35646d8eb90aaa362a3292a9d702e30c9d5af5`, with a dirty working
tree as reported in the executable manifest.

## Exact planar-shock preservation

The numbers below are history maxima through `t=0.02`.  Ring 0 is the annular
cell centred at `r=dr/2`; it is not a point sample at `r=0`.

| scheme | fallback | outcome | max radial p spread | max radial rho spread | max all-domain abs(ur/a) | max ring-0 abs(ur/a) | max odd-axis proxy | shock corrugation |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| LLF-WENO-Z5 | off | reached `t=0.02` | 1.19e-13 | 6.02e-14 | 3.25e-14 | 1.46e-14 | 2.18e-14 | 0 cells |
| LLF-WENO-Z5 | on | reached `t=0.02` | 1.25e-13 | 3.86e-14 | 2.86e-14 | 2.86e-14 | 4.01e-14 | 0 cells |
| AFD-HLLC-WENO-Z5 | off | reached `t=0.02` | 4.04e-13 | 3.27e-13 | 1.80e-14 | 6.33e-15 | 9.25e-15 | 0 cells |
| AFD-HLLC-WENO-Z5 | on | reached `t=0.02` | 7.95e-14 | 2.90e-14 | 1.05e-14 | 3.67e-15 | 5.39e-15 | 0 cells |
| Skew4 central | off | failed first step | -- | -- | -- | -- | -- | -- |
| Skew4 central | on | reached `t=0.02` | 1.93e-13 | 5.85e-14 | 4.05e-14 | 3.79e-14 | 4.67e-14 | 0 cells |
| Skew4-JST | off | failed first step | -- | -- | -- | -- | -- | -- |
| Skew4-JST | on | reached `t=0.02` | 1.54e-13 | 4.44e-14 | 2.85e-14 | 2.45e-14 | 3.09e-14 | 0 cells |

The odd-axis proxy is

`(3*u_r(r=dr/2)-u_r(r=3dr/2))/(2*a(r=dr/2))`.

It vanishes for a locally linear odd radial velocity and avoids the incorrect
requirement `u_r(dr/2)=0`.

Both Skew variants without fallback fail after the first update at
`t=3.001258007e-4` because the shock-normal update produces negative internal
energy.  Central Skew reports `rho*e_int=-57.4093` at `(i,j)=(15,31)`;
Skew-JST reports `rho*e_int=-2.52823` at `(15,30)`.  The initial condition and
failed state are identical in every radial ring.  The reported `i=15` is a tie
selected by the reduction, not localization at the outer radial boundary and
not an axis artifact.

Conclusion: every configuration that survives the shock keeps p, rho, and
u_r radially uniform to round-off, and the front remains exactly planar.  This
small gate reproduces no axis-line artifact.  Skew without local LLF fallback
fails earlier for ordinary shock-robustness reasons, before an axis artifact
can form.

## Seeded near-axis corrugation

A second diagnostic deliberately offsets the initial shock by one axial cell
in alternating fashion over the first four radial rings.  This is not an exact
Riemann solution: it intentionally launches transverse waves, so nonzero
`u_r` and parity proxies are physical/numerical responses and cannot alone be
called an axis bug.

Final values at `t=0.02` for the fallback-enabled schemes are:

| scheme | radial p spread | radial rho spread | max all-domain abs(ur/a) | max ring-0 abs(ur/a) | odd-axis proxy | front corrugation |
|---|---:|---:|---:|---:|---:|---:|
| LLF-WENO-Z5 | 0.3348 | 0.1051 | 0.0522 | 0.0232 | 0.0224 | 0 cells |
| AFD-HLLC-WENO-Z5 | 0.3236 | 0.1000 | 0.0874 | 0.0268 | 0.0122 | 0 cells |
| Skew4 central | 0.2318 | 0.0749 | 0.0774 | 0.0328 | 0.0482 | 0 cells |
| Skew4-JST | 0.1835 | 0.0495 | 0.0567 | 0.0171 | 0.0160 | 0 cells |

All four reduce the deliberately seeded one-cell front corrugation to zero by
the final plot.  These damping rates are diagnostic, not an accuracy ranking.
Fallback-off LLF and AFD also reach the stop time, but retain larger final
pressure spreads (0.487 and 0.518) and larger ring-0 radial Mach magnitudes
(0.0540 and 0.0857).  The two fallback-off Skew runs fail at the first step as
in the exact-planar control.

## Files

- exact matrix: `results/nonlinear_axis_matrix_exact_20260808/summary.json`;
- seeded matrix: `results/nonlinear_axis_matrix_seeded05_20260808/summary.json`;
- each subdirectory contains `run.log`, `metrics.json`, and plotfiles;
- figures: `nonlinear_axis_history.png` in each matrix directory;
- analyzer: `analyze_nonlinear_axis_gate.py`;
- reproducible runner: `run_nonlinear_scheme_matrix.sh`.

The full first-four-ring pressure, density, radial-Mach extrema, signed
corrugation extrema, and detected shock position are stored for every output
time in each `metrics.json`.
