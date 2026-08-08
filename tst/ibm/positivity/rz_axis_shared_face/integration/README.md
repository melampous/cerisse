# Pure-GP R-Z positivity integration smoke

## Qualification status

This fixture currently qualifies only a **single-level (`amr.max_level=0`)**
R-Z limiter update.  It does not qualify multilevel R-Z GP geometry, reflux,
or production physics.  In particular, the raw radial pressure companion is
an auxiliary part of the paired R-Z operator; it is not an independent
conserved flux that can be placed in the ordinary AMReX flux register.

This case compiles the production source-free LLF-WENO-Z5 high-order and
piecewise-constant Rusanov low-order templates in a two-dimensional R-Z
pure shared-GP executable. It evolves complete background cells with complete
background faces. `USE_EB`, cut-control, aggregate, and area-fraction paths
are compile-time forbidden by the local `GNUmakefile`.

The fixed solid is an off-axis circle. `off_axis_circle.dat` supplies only the
BVH topology; `prob.h` supplies one analytic boundary intercept, frame, and
shape operator for the annular BI-CWLS wall closure.

Build and run from this directory:

```bash
make -j8
./main2d.gnu.rz_pure_gp_positivity_smoke.ex inputs
./main2d.gnu.rz_pure_gp_positivity_smoke.ex inputs_force
./main2d.gnu.rz_pure_gp_positivity_smoke.ex inputs_force amr.max_grid_size=16
```

`inputs` is a two-step uniform-state transparency smoke. It must report eight
`HIGH_ORDER theta_min=1` forward-Euler brackets and successful stage and
post-synchronization admissibility audits.

`inputs_force` is a one-step test-only double-rarefaction initial condition.
On the fixed 64x64 grid it must report all four SSPRK43 brackets, one
`LOCAL_SHARED_FACE` result with `theta_min < 1`, a positive number of limited
crossing faces, a roundoff-scale conservation residual, and a PASS verdict.
It must not report `ALL_LOW_ORDER_INADMISSIBLE`, `GLOBAL_THETA`, or
`FINAL_ALL_LOW`.

The `amr.max_grid_size=16` override repeats the forced test with more FABs.
It must retain the same limiting decision (including `theta_min`, to printed
precision) and a roundoff-scale operator-closure residual.  This checks
shared-face synchronization across a different single-level box
decomposition; it is not an AMR or reflux test.

## Expected fail-closed AMR gate

`inputs_amr_noreflux` is intentionally a negative test:

```bash
./main2d.gnu.rz_pure_gp_positivity_smoke.ex inputs_amr_noreflux
```

It must stop before the first RHS evaluation with the existing gate that
RZ-2 annular shared-GP semantics are certified only on a uniform level-0
grid.  Disabling reflux in this input isolates the issue but does not bypass
that geometry/reconstruction gate.  This expected abort is neither a limiter
failure nor an AMR pass.  Do not remove or weaken the annular-GP `max_level=0`
guard to make this fixture run.

Multilevel pressure-companion synchronization and an invariant-domain-
preserving reflux/average-down transaction remain unqualified.  A post-reflux
admissibility audit can detect a bad state but cannot make that synchronization
positivity preserving.

## Operator restrictions

The limiter is a face-flux transaction for source-free Euler.  Problem hooks
that can overwrite the active-cell RHS (including an active-cell
`rhs_nscbc` hook) are forbidden and must fail closed.  The separately guarded
stage-local ghost-cell NSCBC path can pass this particular hook gate because
it does not overwrite the active RHS, but its own coordinate and boundary
qualification gates still apply.  General source terms, viscous terms, and
case-specific RHS corrections are outside this fixture.

The forced case is an integration/branch-execution test. It is not a physical
validation result and does not replace the Mach-4 circle, PM expansion,
compression-corner, curved-wall MMS/BI-pressure, inverse-MOC, or AMR reflux
qualification gates.
