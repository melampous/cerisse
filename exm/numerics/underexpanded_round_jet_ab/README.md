# Reduced 3D underexpanded round-jet sensor A/B test

This harness tests the runtime `cns.afd_shock_llf` switch in a real
three-dimensional underexpanded jet.  It is an isolated numerical regression.
It does not alter the canonical jet case or any production source.

## Physical problem

The local `prob.h` includes
`exm/underexpanded_jet/3d/npr3_lev4_ns/prob.h`.  The following parts are
therefore identical to that canonical source-only case.

- Sonic round-jet exit with stagnation-to-ambient pressure ratio `NPR_0=3`.
- Exit diameter `D_e=10 mm` and a `100 micrometre` tanh lip profile.
- Ambient air at `101325 Pa` and `300 K`.
- Calorically perfect gas with `gamma=1.4`.
- Constant molecular viscosity and thermal conductivity.
- The same initial condition and physical boundary functions.

Only the computational extent and resolution are reduced.  The domain is
`0 <= x <= 6 D_e` and `-2 D_e <= y,z <= 2 D_e`.  The uniform
`96 x 64 x 64` Level-0 mesh gives 16 cells across `D_e`.  The calculation uses
SSPRK(3,3), `CFL=0.30`, and stops at `70 microseconds`.

This resolution is intended for a quick A/B stability and selectivity test.
It is not a grid-converged Mach-disk calculation.  The analysis consequently
reports the strongest centreline compression in a prescribed near-field
window.  That feature may be described as a Mach disk only after the pressure
and Mach profiles confirm that a normal shock has formed.  If the feature is
not established by `70 microseconds`, repeat both variants with
`stop_time=1.2e-4`.

## Controlled comparison

One CUDA executable contains AFD-HLLC-WENO-Z5 and the viscous terms from the
canonical problem.  The paired solver commands differ only in
`cns.afd_shock_llf`.  They run in separate output directories.  The fixed
sensor settings are

```text
cns.afd_smoothness_threshold = 0.08
cns.afd_shock_pressure_jump  = 0.01
cns.afd_shock_compression    = 0.001
```

Both IBM and AMReX EB are disabled.  This is a full-cell Cartesian,
non-immersed calculation.

## Run and analyse

From this directory:

```bash
bash run_ab.sh
```

Set a fresh result directory with `JET_AB_RESULTS_ROOT`.  Set build concurrency
with `CERISSE_BUILD_JOBS`.  The script refuses to reuse an existing result
directory.  It records the executable and relevant source hashes before
running the two variants.

The analyser writes JSON, CSV, Markdown, and
`jet_sensor_ab_diagnostics.png`.  It reports positivity, centreline pressure
and Mach data, the strongest axial compression location and thickness,
cross-plane symmetry, an odd-even pressure measure, and full-field A/B
differences.  It also reports a cell-centred offline proxy for the new joint
pressure and compression gate.  The proxy is diagnostic only.  It is not an
instrumented count of faces selected inside the solver.
