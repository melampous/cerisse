# Skew R-Z unweighted-JST validation (2026-08-07)

## Source revisions under test

```text
src/rhs/Skew.h          58a1bbff525bda5df30d7dd38cc78de450e38360961c09a34e26f901a161c8ec
src/tim/compute_rhs.cpp 23ef010c74d2d541437d888e258ca0975ff511b9268ca30aa82bd220b5d550a2
```

The R-Z/Cartesian MMS matrices below correspond to that `58a1bb...` Skew
revision.  While the pure-GP executable was being built, `Skew.h` advanced to

```text
src/rhs/Skew.h          9d50e7ca27e6671fc19a14b6e243658dc64b925c35bb0910738b34f3f48898c3
```

The intervening edit only allows/defaults `skew_axis_collar_faces=0` and makes
the selected axis post-pass follow that choice; it does not change the JST
physical-state dissipation being verified.  The final pure-GP smoke was
rebuilt and rerun against `9d50e7...`, explicitly with collar width zero.

No production source was edited by this validation.  The test-local pure-GP
fixture under `gp_smoke/` is a copy of
`tst/ibm/positivity/rz_axis_shared_face/integration` whose RHS typedef alone
was changed from LLF-WENO-Z5 to Skew4-JST.

## R-Z smooth Euler MMS

The production point-sampled, uniform-Level-0 R-Z MMS was built with
Skew4-JST (`C2=1.5`, `C4=0.016`) and run on N=16, 32, and 64.  Both an
assembled-RHS check and an SSPRK3 evolution to `t=0.01` were run.  The time
step scales as dx squared (8, 32, and 128 steps).

Build command:

```bash
cd exm/numerics/mms_rz_navier_stokes
make -j8 RZ_MMS_PHYSICS=euler RZ_MMS_EULER_SCHEME=skew-jst-o4 \
  USE_MPI=FALSE USE_OMP=FALSE USE_CUDA=FALSE \
  USE_PELEPHYSICS=FALSE USE_EB=FALSE USE_GPIBM=FALSE
```

Executable SHA-256:

```text
db30371cbe9be1ab8fa95f10985caca748f98f9b822a0224e3dad5b0e81e118f
```

Representative RHS command (N is replaced by 16, 32, or 64):

```bash
./main2d.gnu.mms_rz_euler.skew_jst_o4.ex inputs \
  max_step=1 stop_time=1.0 time_step=1.0e-5 \
  cns.order_rk=0 cns.stages_rk=1 "amr.n_cell=N N" \
  amr.plot_int=1 amr.plot_file=RESULT/plt \
  cns.skew_first_order_fallback=1 cns.skew_axis_collar_faces=1
```

RHS cylindrical-volume-weighted L2 errors and rates:

| N | density | rate | radial momentum | rate | axial momentum | rate | energy | rate |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 16 | 1.476346e-3 | - | 1.553210e-3 | - | 1.947422e-3 | - | 5.852980e-3 | - |
| 32 | 1.955734e-4 | 2.916 | 2.162413e-4 | 2.845 | 2.440956e-4 | 2.996 | 7.964932e-4 | 2.877 |
| 64 | 2.479295e-5 | 2.980 | 2.772019e-5 | 2.964 | 3.030157e-5 | 3.010 | 1.026465e-4 | 2.956 |

The first-ring radial-momentum Linf error is `7.837435e-5`,
`5.496271e-6`, and `6.034071e-7`, giving rates 3.834 and 3.187.  Thus the
axis-local error is not hidden by the cylindrical volume weight.

The evolved-state L2 errors exhibit the same asymptotic behavior:

| N | density | rate | radial momentum | rate | axial momentum | rate | energy | rate |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 16 | 1.491353e-5 | - | 1.554685e-5 | - | 1.986334e-5 | - | 5.729773e-5 | - |
| 32 | 1.976638e-6 | 2.916 | 2.165576e-6 | 2.844 | 2.502855e-6 | 2.988 | 7.780603e-6 | 2.881 |
| 64 | 2.501529e-7 | 2.982 | 2.774295e-7 | 2.965 | 3.110963e-7 | 3.008 | 9.991335e-7 | 2.961 |

Fallback/collar transparency was checked by repeating every RHS and evolution
run with `cns.skew_first_order_fallback=0`.  `diff -rq` reports no difference
between ON and OFF plotfiles at all three resolutions; their `Cell_D_00000`
SHA-256 values are identical pairwise.  This is stronger than merely retaining
the same convergence rate: the smooth field did not select the fallback or
axis collar.

Raw tables are in:

- `rz_rhs_on/errors.csv` and `rz_rhs_off/errors.csv`
- `rz_evolve_on/errors.csv` and `rz_evolve_off/errors.csv`

## Cartesian non-regression

The Cartesian Euler MMS was run to `t=0.01` at N=16, 32, and 64 with the
current source and with the pre-fix Skew-fallback binary built at 13:36 on
2026-08-07.  Fallback was explicitly disabled in both binaries.

```text
pre-fix OMP binary SHA-256: 7b853bf91783322eb946a861bbcf12e0bbb310df31afde8922fe5960c030f078
current binary SHA-256:     538932dc7c8ba224b11c1110b6764ffccb4f72cb321097d313163c6ca7ae164b
```

For N=16, 32, and 64, the complete final plotfile trees are byte-for-byte
identical.  The density L2 rates are 3.086 and 3.012; energy L2 rates are
2.988 and 2.916, consistent with the expected approximately third-order
Skew4-JST behavior for the present first-power sensor.

Results are in `cart_current/` and `cart_prefx/`; the analyzed current table is
`cart_current/errors.csv`.

## Pure shared-GP R-Z CPU smoke

The copied fixture compiles with:

```text
USE_GPIBM=TRUE
USE_EB=FALSE
IBM_EXPERIMENTAL_CUT_CONTROL=FALSE
```

Build and run commands:

```bash
cd analysis/skew_axis_diagnosis/jst_unweighted_validation/gp_smoke
make -j8 USE_MPI=FALSE USE_CUDA=FALSE USE_OMP=FALSE \
  USE_GPIBM=TRUE USE_EB=FALSE IBM_EXPERIMENTAL_CUT_CONTROL=FALSE
./main2d.gnu.rz_pure_gp_skew_jst_smoke.ex inputs \
  max_step=1 \
  cns.ibm_positivity_flux_limiter=0 \
  cns.ibm_positivity_flux_limiter_verbose=0 \
  cns.skew_first_order_fallback=1 cns.skew_axis_collar_faces=0
```

The build and one-step SSPRK(4,3) run passed.  The runtime manifest reports:

```text
[Numerics] inviscid_scheme=skew4-jst
[Numerics] rz_axis_pressure_auxiliary=order_matched_parity_reconstruction selector=skew4-jst stencil=marker_aware_shared_gp
[IBM-Method] family=pure_shared_gp_full_cartesian pure_gp_production=1 cut_control_compiled=0 control_volume=full_cartesian_cell face_area=full_cartesian_face solid_extension=shared_ghost_point
```

All 42 ghost points used the quadratic BI-constrained reconstruction.  The
uniform-state plot after one step has density range
`[1.176412770037523, 1.1764127700383231]`, radial/axial momentum extrema below
`2.51e-14`, and energy range
`[253312.49999999977, 253312.5000000003]`.  This is a roundoff-level uniform
transparency result across both all-fluid and shared-GP stencils.

The pure-GP executable SHA-256 is:

```text
ee8622bdedbf6b59308ebb3d89beb3ce46e9dfb51d9c8f191d483de2804648ef
```

## Scope and limits

- The MMS is Euler, smooth, single-level, and non-IBM; it does not qualify
  viscous R-Z, AMR reflux, or discontinuous solutions.
- The pure-GP check is a one-step uniform transparency and template/runtime
  integration smoke with the IBM positivity limiter intentionally disabled.
  It does not qualify Skew's production positivity limiter, shock robustness,
  surface-force accuracy, or the mandatory long IBM regression ladder.
- No cut-cell, cut-control-volume, area/volume-fraction, aggregate, or EB flow
  update is compiled or used by the pure-GP executable.
