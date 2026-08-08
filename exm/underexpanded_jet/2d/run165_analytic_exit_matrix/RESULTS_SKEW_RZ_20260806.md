# Run165 analytic-exit Skew--JST RZ test (2026-08-06)

## Scope and invariant controls

- Two-dimensional axisymmetric RZ Euler equations, no IBM and no solid body.
- Run165 freestream and analytic virtual-exit states are unchanged.
- The axial-high circular aperture imposes the counterflow jet directly.
- Uniform L1 coverage: 64 x 192 = 12,288 finest-level cells, with no AMR
  interface in the active fine solution.
- L1 spacing is 1.06129666667 mm in both directions.  The 12.73556 mm exit
  diameter is resolved by exactly 12 L1 cells.
- Jet conservative-state ramp time: 100 microseconds.
- Skew--JST: fourth-order central base operator, first-power sensor,
  C2 = 1.5 and C4 = 0.016.
- Pure all-fluid test only.  No GP-IBM, cut-cell, EB, aggregate, or
  cut-control path is compiled or used.

The executable was built from Cerisse revision `dc35646d8` with a dirty-tree
manifest because the current solver worktree contains uncommitted research
changes.  The identical locally built OpenMP executable was copied to
Y9000X, so the local and remote runs did not use different solver sources.

## Results

| Test | Result | First failed physical time | Failed L1 cell | Failed location (r,x) | Minimum rho | Minimum rho e_int |
|---|---|---:|---|---|---:|---:|
| Jet, CFL 0.20 | failed | 13.37030308 us | (1,187) | (1.591945,-2.334566) mm | 0.0275849 | -56.1278 J/m3 |
| Jet, CFL 0.10 | failed | 13.36972519 us | (1,187) | (1.591945,-2.334566) mm | 0.0275912 | -46.9173 J/m3 |
| Jet, CFL 0.10, Y9000X replay | failed identically | 13.36972519 us | (1,187) | (1.591945,-2.334566) mm | 0.0275912 | -46.9173 J/m3 |
| Jet-off-equivalent control, CFL 0.20, Y9000X | passed to 120 us | none | none | none | positive | positive |

The jet-off-equivalent control used a `1e30 s` ramp with the same binary, so
the aperture state remained freestream to floating-point accuracy.  It
completed 256 coarse steps to 120 microseconds in 3.216 s of advance time on
one Y9000X CPU thread, or 0.01256 s per coarse step on average.

The last saved positive jet field is step 28 at 13.13604343 microseconds.  It
already contains a local Mach maximum of 24.283, pressure minimum of 21.85 Pa,
and pressure maximum of 24.892 kPa.  The subsequent failure has positive
density but negative internal energy.  It is inside the developing jet, not
on the RZ axis and not in the boundary ghost cells.

## Interpretation

Reducing CFL by two did not materially move the failure in physical time and
did not change its cell.  The baseline RZ/free-stream control is stable.
Therefore this is not an IBM geometry, AMR-interface, or large-time-step
failure.  The fourth-order Skew--JST spatial update lacks a conservative
positivity guarantee for this high-pressure counterflow interaction.  The
current all-fluid run also reports the input `cns.stage_positivity` as unused;
the existing production positivity machinery is not protecting this case.

Skew--JST with these coefficients is not ready for the Run165 jet production
calculation.  Increasing C2 can be tested as a dissipation sensitivity, but it
cannot provide a mathematical positivity guarantee.  The format-independent
repair remains a conservative shared-face low-order LLF/HLLE fallback plus a
positivity-preserving flux limiter applied to every forward-Euler bracket.

## Reproduction

Build:

```bash
make -j8 EULER_SCHEME=skew JET_ENABLED=1 NS_ENABLED=0 LES_ENABLED=0 \
  USE_MPI=FALSE USE_OMP=TRUE \
  EBASE=run165_analytic_exit_skew_jst_o4_rz_omp \
  TMP_BUILD_DIR=tmp_build_run165_skew_jst_o4_rz_omp
```

Run the CFL 0.1 gate:

```bash
OMP_NUM_THREADS=1 ./run165_analytic_exit_skew_jst_o4_rz_omp2d.gnu.TPROF.OMP.ex \
  inputs max_step=2000 stop_time=1.2e-4 cfl=0.1 \
  prob.jet_ramp_time=1.0e-4 amr.plot_files_output=0 \
  amr.checkpoint_files_output=0
```

Y9000X archive directory:

```text
/home/melampous/cerisse_run165_skew_rz_20260806
```
