# Fixed-STL IBM Pressure Production Contract

This document defines the supported R2 path for pressure on an arbitrary fixed
curve or STL. It applies to the single-component ideal-gas GPIBM solver.

## Recommended Method

Use one ghost layer and second-order image/surface reconstruction:

```cpp
static constexpr int interp_order = 2;
static constexpr int extrap_order = 2;
static constexpr int interp_order_surf = 2;
static constexpr int extrap_order_surf = 2;
static constexpr int ghost_layers = 1;
static constexpr Real alpha = 1.0;
static constexpr Real alpha_surf = 1.0;
static constexpr ibm_pressure_closure_t pressure_closure =
    ibm_pressure_closure_t::shock_aware_fluid_extrapolation;
static constexpr Real pressure_sensor_low = 0.03;
static constexpr Real pressure_sensor_high = 0.10;
```

For Euler use a slip wall. For Navier--Stokes use the wall model required by
the physical problem; generic shock-dominated STL pressure still uses the
shock-aware closure. Analytic curvature compatibility is reserved for a known
analytic surface and is not inferred from STL facets.

## Geometry Contract

The exact same input geometry file must be used at every fluid resolution.
Record its SHA-256, point/facet count, bounds, orientation, closed-manifold
topology, normal closure, and maximum edge length relative to the finest fluid
spacing. The solver may correct global orientation, but it does not improve
surface location or facet normals.

Audit before running:

```bash
python3 tools/audit_ibm_geometry.py body.stl --kind stl \
  --dx DX32 DX48 DX64 DX96 --max-edge-over-dx 0.5 \
  --require-outward-input --json geometry_audit.json
```

For a 2-D curve add `--check-self-intersections`. STL self-intersection still
requires a robust CAD/geometry-kernel check; this tool deliberately does not
claim to provide one.

Resolution labels must state cells per body reference length, not just cells
across the entire domain. For the radius-0.1 circle in a unit box, `D32` means
`160 x 160` domain cells, whereas a `32 x 32` domain has only 6.4 cells across
the diameter and is an under-resolution failure sample.

## Time And Output Contract

Compare either the same physical time or a demonstrated steady/statistically
stationary physical-time window. Equal step counts are invalid when `dt`
changes with grid spacing. The executable writes a final IBM VTP and a final
momentum-budget row when the final step is not an output-interval multiple, so
`stop_time` can be used as the common endpoint.

For a fixed surface sequence require identical VTP points/connectivity:

```bash
python3 tools/integrate_ibm_pressure.py surf_D32.vtp surf_D48.vtp \
  surf_D64.vtp surf_D96.vtp --labels D32 D48 D64 D96 \
  --p-ref P_INF --rho-ref RHO_INF --u-ref U_INF \
  --reference-area A_REF --require-same-mesh --json pressure.json
```

Report the full `Cp` profile, maximum and upstream-cap stagnation pressure,
pressure drag, symmetry lift, nonpositive pressure count, minimum pressure,
image-point quality, shock sensor, and fallback measure fraction. Pressure
below freestream is not automatically a numerical undershoot in an expanding
supersonic flow; it is reported as a diagnostic, not a standalone failure.

## Verified Mach-4 Circle Sequence

The fixed geometry is a 4,096-segment radius-0.1 circle with SHA-256
`53775b9fc2a2b67a354b72ca4fcbb746810656d98a926c5effec123abb424e6e`.
All runs use LLF--WENO-Z5, Euler slip, the recommended IBM settings, one uniform
level, CFL 0.3, and exact `t=0.008`. The force window is `0.007 <= t <= 0.008`.

| Grid | Domain cells | Cp max | Pressure Cd | |Cl| | fallback area | surface/CV force mismatch | drag drift |
|---|---:|---:|---:|---:|---:|---:|---:|
| D32 | 160^2 | 1.77759 | 1.18531 | 2.0e-6 | 2.29% | 5.70% | 0.12% |
| D48 | 240^2 | 1.78378 | 1.24704 | 4.99e-3 | 0.95% | 2.64% | 2.01% |
| D64 | 320^2 | 1.78545 | 1.24111 | 3.29e-3 | 1.44% | 1.77% | 1.20% |
| D96 | 480^2 | 1.78510 | 1.20142 | 4.07e-4 | 1.22% | 0.92% | 0.66% |

For gamma 1.4 and Mach 4, the normal-shock plus isentropic-stagnation
(Rayleigh-pitot) value is `Cp0=1.79179`; D96 is 0.37% low. Modified Newtonian
theory gives the 2-D circular-cylinder estimate `Cd=(2/3)Cp0=1.19453`; the D96
pressure drag is 0.58% high. These are useful supersonic reference checks, not
an exact replacement for a resolved benchmark solution.

The front stagnation pressure is credible at D64/D96 and the D96 integrated
force closes below one percent. The full sequence is not yet a clean steady
grid-convergence result: D48 and D64 retain measurable drag drift and symmetry
breaking, and the complete surface-pressure error is non-monotone. Production
status is therefore:

- front/stagnation pressure: accepted for this circle at D96;
- D96 source-free integrated pressure force: accepted against the conservative
  control-volume force over the stated late window;
- long-time full `Cp` distribution and grid-converged drag: not yet accepted;
- arbitrary STL: requires the double-wedge, 3-D sphere, and audited complex-STL
  matrix before claiming general validation.

## Remaining R2 Matrix

Run Mach-2 double wedge, 3-D supersonic sphere/bow shock, and one fixed audited
high-accuracy complex STL at multiple fluid resolutions. Keep the geometry
hash and facet count invariant. Sharp corners must be reported separately from
smooth faces because a pointwise pressure value at a mathematical corner is
not a smooth convergence target.
