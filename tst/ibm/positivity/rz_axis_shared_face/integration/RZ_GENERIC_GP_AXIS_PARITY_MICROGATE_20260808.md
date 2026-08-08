# Generic point-GP R-Z axis-parity micro-gate

## Scope

This is a CPU-only, level-0, 32 by 32 diagnostic for the generic
non-annular `computeAllGPs` path.  It retains
`pure_shared_gp_full_cartesian`: full background cells and faces, one shared
solid-side GP state, and no cut cell, EB flow update, aggregate, area fraction,
or partial volume.

The test-local `ProbIB` derives from the production solver and snapshots the
same primitive `MultiFab` immediately before and after
`Base::computeAllGPs`.  It checks the physical-axis pairs

```
negative i  -1  -2  -3
mirror i     0   1   2
```

with even parity for density, pressure, and axial velocity, and odd parity for
radial velocity.  The axis-intersecting semicircle and an axis-regular but
non-wall-compatible initial state make GP publication observable without a
long flow integration.

## Reproduction

From this directory:

```bash
make -j8 TEST_RZ_GP_PARITY_PROBE=TRUE
./main2d.gnu.rz_pure_gp_positivity_smoke.generic_point_gp_axis_parity_probe.ex \
  inputs_generic_gp_axis_parity
./main2d.gnu.rz_pure_gp_positivity_smoke.generic_point_gp_axis_parity_probe.ex \
  inputs_generic_gp_axis_parity prob.apply_post_gp_axis_parity=1
```

The first command now verifies the repaired generic production path.  The
second retains the former test-local A/B hook; after the production repair it
is intentionally redundant and must give the same zero defects.

## Result

There are two positive-radius GP mirrors at every tested depth.  Before GP
publication, every measured parity residual is exactly zero.  After generic
GP publication, the largest GP-pair residuals are:

| pair | density | pressure (Pa) | axial velocity (m/s) | radial velocity (m/s) |
|---|---:|---:|---:|---:|---:|
| `(-1,0)` | 7.657777088e-5 | 9.234057237 | 48.80448275 | 5.907409883 |
| `(-2,1)` | 3.615282023e-4 | 43.60093739 | 40.34369896 | 13.24198895 |
| `(-3,2)` | 1.105096140e-3 | 133.3117427 | 53.27069478 | 36.84154706 |

One direct depth-1 sample at axial index 20 is:

| field | pre negative | pre positive | post negative | post positive |
|---|---:|---:|---:|---:|
| density | 1.176441490 | 1.176441490 | 1.176441490 | 1.176440964 |
| pressure | 101329.9475 | 101329.9475 | 101329.9475 | 101329.8841 |
| axial velocity | 40.00122070 | 40.00122070 | 40.00122070 | -8.803262044 |
| radial velocity | -1.25 | 1.25 | -1.25 | -3.838538713 |

Thus the negative-radius ghost retains the pre-reconstruction mirror while
the positive solid-side GP receives its new shared state.  The final
`FillBoundary(periodicity)` in generic `computeAllGPs` exchanges same-level
FAB data but does not impose a physical boundary condition, so it cannot
repair this mismatch.

The test-local post-GP call to `fillRZAxisPrimitiveParity` reduces every one of
the twelve reported residuals (four fields by three depths) to exactly zero.

## Conclusion and production change

The generic non-annular R-Z GP path had a real stale physical-axis ghost gap.
Production now performs

```cpp
if (rz_point_target < 0) {
  fillRZAxisPrimitiveParity(prims_mf, lev);
}
```

after the final same-level `FillBoundary` in `computeAllGPs`.  Putting the call
inside `computeAllGPs`, rather than at its current callers, covers the RHS,
end-of-step publication, and positivity-limiter reconstruction uniformly.
The target guard leaves annular quadrature passes untouched; their wrapper
already closes all target accumulation with one final parity refresh.  The
helper is a no-op outside R-Z or away from an `r=0` domain boundary.

This change remains pure shared-GP/full-cell IBM: it only mirrors the one
already-published shared GP into coordinate ghost storage.  It does not change
the GP geometry, BI-CWLS state, active-cell update, full-cell volume, full face
area, force definition, or communication semantics, and it introduces no
cut/EB/aggregate mechanism.

The annular R-Z production candidate already performs this final refresh in
`computeAllGPsRZAnnular`, so this specific generic-path defect must not be used
as the explanation for an artifact observed in an annular Run165 case without
separate evidence.

## Post-repair verification

With the production helper active and the test-local post hook disabled, the
three depths and all four parity fields report exactly zero post-GP residual.
The method manifest remains `pure_shared_gp_full_cartesian`, with
`cut_control_compiled=0` and no EB or aggregate flow path.

The guarded annular path was also checked with an exact pre/post executable
A/B.  Its two-step plotfiles agree for every stored variable with absolute and
relative tolerances both set to zero.  The guard therefore repairs the generic
point-GP path without perturbing the already-correct annular target loop.
