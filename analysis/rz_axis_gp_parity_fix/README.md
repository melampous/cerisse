# Generic point-GP R-Z axis-parity repair

## Scope

This is a pure shared-GP/full-background-cell repair.  It updates only the
negative-radius physical coordinate ghosts after the final generic point-GP
publication.  It does not change the unique GP state, BI-CWLS reconstruction,
thermodynamic extension, active state, full-cell flux, positivity limiter,
surface pressure/force, cell volume, face area, or AMR flux register.  Cut,
EB, cut-control and aggregate flow mechanisms remain disabled.

The `rz_point_target < 0` guard is required.  Annular target calls accumulate
several intermediate point reconstructions and already apply one final parity
refresh in `computeAllGPsRZAnnular`; the generic repair must not alter those
intermediate target states.

## Direct gate

Case:
`tst/ibm/positivity/rz_axis_shared_face/integration/inputs_generic_gp_axis_parity`

The 32 by 32 CPU gate places a pure-GP semicircle across the physical R-Z axis
and checks three negative ghost depths.  Before the repair, the largest stale
post-GP defects were 1.1051e-3 in density, 133.312 Pa in pressure, 53.271 m/s
in axial velocity and 36.842 m/s in radial odd parity.  With the production
repair, all twelve measured defects are exactly zero.  The test now aborts if
any defect exceeds 1e-12.

The accepted manifest is:

```
family=pure_shared_gp_full_cartesian
pure_gp_production=1
cut_control_compiled=0
```

Source and executable hashes at the gate:

```
src/ibm/ibm_solver.h  da4924bb2c294e34504ba5ced951155e99f387c9e6e02e0945042dfe8f457158
generic CPU gate      8daee6a3597354541b66d3d120511f096af4e9e77aa2d692e2974f459980fa27
```

## Annular non-regression

The existing two-step annular shared-GP integration was run with the
pre-repair executable and the current executable.  AMReX `fcompare` with both
absolute and relative tolerance set to zero reported `PLOTFILE AGREE` for
Density, all three momenta, Energy, `sld`, and `ghs`.  This confirms that the
target guard leaves the already-correct annular path bit-identical.

Full R-Z operator results are collected in
`exm/numerics/rz_paired_pressure_axis/RZ_AXIS_DIAGNOSTIC_20260808.md`.
