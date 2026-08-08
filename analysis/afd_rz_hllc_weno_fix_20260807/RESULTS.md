# AFD-HLLC-WENO-Z5 R-Z troubled-face repair and Run165 test

Date: 2026-08-07

## Decision

The Run165 all-fluid R-Z jet case now passes the requested numerical stability
gate with a sensor-gated, shared-face HLLC-MUSCL fallback. The repaired branch
completed 1.8 ms from a clean initial condition without a negative density or
internal-energy state in any logged Heun Forward-Euler candidate or final
stage.

This is a numerical robustness result, not yet a physical validation of the
developed SRP flow. The 1.8 ms run deliberately disabled plot/checkpoint output
and therefore does not establish shock-location, force, or experimental-data
agreement.

## Why the supplied Weno.h patch was not copied literally

The supplied patch repairs nonlinear WENO reconstruction of the radial metric
LLF split flux in `Weno.h`. The failing AFD case uses a different operator:

- AFD performs point-value WENO reconstruction of characteristic primitive
  states and then calls HLLC.
- Its R-Z metric treatment appears in the fixed linear AFD correction to
  `rF`; it does not compute WENO weights from the radial metric split flux.
- The first recorded Run165 failure is on face `(4,171), dir=1`, an axial
  face. A collar acting only on the first positive radial face cannot select
  that face.

The transferable idea is nevertheless valid: detect a genuinely troubled
shared face and replace that face with a conservative, lower-order,
metric-consistent flux. That is the repair used here.

## Failure mechanism and selected repair

With all fallbacks disabled, the case fails at `t=73.8 us`, L1 cell `(4,171)`.
The point-WENO right state at the responsible axial face remains positive in
density and pressure but reconstructs `u_n=-8860.998 m/s` and `p=56.140 Pa`.
HLLC then returns an energy flux of `-1.070123e11`, about three orders of
magnitude larger than the local physical-scale flux, and the RK update creates
negative internal energy.

The selected AFD-WENO default is:

```text
cns.afd_shock_hllc_muscl = 1
cns.afd_shock_llf = 1
cns.afd_hllc_reconstructed_speed_ratio = 0
```

A face switches to HLLC-MUSCL only when the common detector finds all of:

1. a nonsmooth six-cell stencil;
2. a relative pressure jump above `0.01`;
3. dimensionless local compression above `0.001`.

The limited characteristic MUSCL states and HLLC flux are computed once per
shared background-grid face. In R-Z, the fallback returns the ordinary
per-area face flux; the R-Z divergence applies the face-radius metric, and the
axis face remains zero. If the MUSCL candidate is inadmissible or non-finite,
the existing characteristic LLF branch remains available. There is no
cell-wise clipping and no nonconservative state replacement.

The optional reconstructed-speed guard is retained as an experimental A/B
tool but is default-off. By itself it catches the first extreme reconstruction
but is not a positivity guarantee: tested thresholds later failed between
350.8 and 416.6 us.

TENO5 retains `cns.afd_shock_hllc_muscl=0` by default because it was not part
of this validation.

## Run165 stability matrix

All Run165 entries use two-dimensional all-fluid Euler in R-Z coordinates,
Run165 analytic virtual-exit jet data, a 100 us conservative-state jet ramp,
fixed coarse `dt=0.4 us`, Heun RK2, and a uniform L1 finest mesh of `64x192 =
12,288` cells. L1 spacing is `1.0612967 mm`; there is no IBM or solid geometry.

| Configuration | Result |
| --- | --- |
| Original AFD, all fallbacks off | FAIL at 73.8 us, cell `(4,171)` |
| Reconstructed-speed guard, ratio 2.0 | FAIL at 77.6 us |
| Speed guard ratios 1.25 / 1.10 / 1.05 | FAIL at 350.8 / 368.4 / 416.6 us |
| Sensor-gated HLLC-MUSCL, `shock_llf=0` | PASS to 1.8 ms, 4500 coarse steps and 9000 L1 substeps |
| Rebuilt final WENO defaults | PASS to 1.8 ms, 4500 coarse steps and 9000 L1 substeps |
| Rebuilt default checkpoint replay | Known bad face selects HLLC-MUSCL and passes the step |

The 1.8 ms run had zero negative-stage records and a minimum logged internal
energy density of `501.8656557 J/m^3`. It ran in `315.80 s` on 8 OpenMP CPU
threads, or `0.0702 s/coarse-step`, with detailed stage diagnostics enabled.
The exact-default 500 us run also had zero negative-stage records; its minimum
was `586.7274701 J/m^3`. A final exact-default rerun after the HLLC zero-wave
boundary correction reached 1.8 ms in `189.73 s` on one MPI rank. It completed
4500 coarse steps and 9000 L1 substeps without a state-check failure or a
non-finite value. Plot and checkpoint output were disabled for this numerical
stability gate.

Primary logs:

- `exm/underexpanded_jet/2d/run165_analytic_exit_matrix/results/afd_rz_local_fallback_20260807/troubled_muscl_1p8ms_stagecheck.log`
- `exm/underexpanded_jet/2d/run165_analytic_exit_matrix/results/afd_rz_local_fallback_20260807/default_rebuilt_500us_stagecheck.log`
- `exm/underexpanded_jet/2d/run165_analytic_exit_matrix/results/afd_rz_local_fallback_20260807/default_rebuilt_chk184_face_trace.log`
- `exm/underexpanded_jet/2d/run165_analytic_exit_matrix/results/afd_rz_local_fallback_20260807/baseline_pure_80us.log`

## Smooth-solution non-regression

R-Z Euler MMS compared `shock_hllc_muscl=0/1` with the other controls fixed.
For N=16, 32, 64, and 128, every corresponding final plotfile data object was
byte-identical. The N=64 to 128 L2 orders were 4.96999, 4.97034, 5.00639, and
4.98671 for density, radial momentum, axial momentum, and total energy.

Cartesian Euler MMS at N=16 and 32 was likewise byte-identical between the two
branches. Its corresponding L2 orders were 5.1142, 4.9194, 5.0372, and 5.1358.

Artifacts:

- `exm/numerics/mms_rz_navier_stokes/results/afd_speedguard_mms_20260807/latest_muscl/current/RESULTS.md`
- `exm/numerics/mms_cartesian/results/afd_muscl_mms_ab_20260807/`

## Pure shared-GP scope gate

The all-fluid HLLC-MUSCL and reconstructed-speed branches are explicitly
disabled in `AfdIBMDrivers.h`. A current CPU/OMP pure-GP AFD build and 3-step
smoke passed. Its runtime manifest reported:

```text
family=pure_shared_gp_full_cartesian
pure_gp_production=1
cut_control_compiled=0
control_volume=full_cartesian_cell
face_area=full_cartesian_face
solid_extension=shared_ghost_point
```

No cut-cell, cut-control, aggregate, or EB flow path was enabled. During this
gate, a pre-existing C++ scope error around the generic y-face trace variable
was exposed and fixed with braces around the existing IBM/EB `if (!wallflx)`
body. That edit changes no state, flux, or IBM algorithm.

Pure-GP gate artifacts:

- `IBM/cases/ibm_tests/2d_bvh_gpu/results/afd_guard_puregp_smoke_20260807/README.md`
- `IBM/cases/ibm_tests/2d_bvh_gpu/results/afd_guard_puregp_smoke_20260807/smoke_3step_bf32.log`

## Modified production-facing files

- `src/rhs/Afd.h`: WENO default, shared-face selection, optional speed guard,
  validation, and R-Z-safe propagation.
- `src/rhs/AfdIBMDrivers.h`: explicitly disables the new all-fluid branches in
  GP-IBM/EB drivers.
- `src/rhs/Riemann.h`: corrected the HLLC star-state total-energy expression
  to retain the star-density ratio in the pressure-work term. The physical
  left and right branches now include the zero-wave cases, which avoids an
  unnecessary singular star-state evaluation when `S_L=0` or `S_R=0`.
- `docs/input.md`: runtime option and scope documentation.

Current SHA-256 values at validation handoff:

```text
87ec4c1791bc494ac40b8a8dd7beb8d51455aac9805beb9ffece01e076275ecd  src/rhs/Afd.h
5eae587b311ad36aff82fa0f5cc361f4117ef8bf79d749ceb2629d34b12b456e  src/rhs/AfdIBM.h
9988b8d4f7e2af329ed0fcac58b78eb5b611ccd7b593993e61f331887610c256  src/rhs/AfdIBMDrivers.h
f4a78c5c099666132968c7ae7f2044bebe5f3522fa81de92a6edb6be6aa1ba66  src/rhs/Riemann.h
22a95f366c05014d713d41a95efec1a57b6d3765590ef2159982fc59ac2af3cf  docs/input.md
15eb19ef28b04f4782027e27cca52325157d773346239171ea190ef72f9438e3  run165_analytic_exit_afd-hllc-weno_jet1_js1p0_cpu2d.gnu.TPROF.MPI.ex
```

## Remaining gates before physical production use

1. Produce plotfiles and verify shock topology, standoff, mass flow, and force
   histories. Stability alone is not Run165 validation.
2. Profile the default-off speed-guard code on the target GPU; if it increases
   registers or reduces occupancy despite being disabled, isolate it behind a
   compile-time diagnostic switch or remove it.
3. Do not claim TENO5 validation from these AFD-WENO tests.
