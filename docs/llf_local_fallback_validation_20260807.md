# Shared face-local LLF fallback validation

Date: 2026-08-07.

## Implemented scope

The sensor-gated first-order LLF fallback is implemented in the shared
`characteristic_llf_flux_t` driver in `src/rhs/Weno.h`. It therefore applies
to LLF-WENO-Z5, LLF-TENO5, and LLF-TENO6 in Cartesian and R-Z coordinates,
for both all-fluid and pure shared-GP stencils. It does not apply to the AFD
HLLC or skew-symmetric drivers, which have separate robustness mechanisms.

The implementation retains the pure shared-GP method. It writes one complete
background-face flux and does not use cut volumes, open-face fractions,
cut-control residuals, aggregate cells, or EB flow updates. The same finalized
face value is available to both adjacent background cells and to AMR refluxing.
This preserves discrete full-background-grid conservation. It does not imply
exact conservation on the true curved fluid domain.

## Sensor and fallback flux

At a general face, activation requires both

1. an adjacent relative density jump above `0.003`, and
2. a normalized density second difference above `0.08` in at least one of the
   two three-point groups that contains both face-adjacent cells.

The combined gate prevents a smooth but finite density gradient from causing
first-order reconstruction. At the first interior radial face next to the R-Z
axis, the independently validated density-jump gate is retained without the
second-difference condition. The axis face at `r=0` is excluded.
The axis and first interior face are identified by the radial mesh index.
The global R-Z geometry check now requires `geometry.prob_lo[0]` to be
exactly zero, so a roundoff-sized positive radius cannot enter an axis-face
division.

Cartesian and axial faces use the two-state LLF flux. A radial R-Z face uses
the radius-weighted physical-flux average and an unweighted conservative-state
jump. The all-fluid paired-pressure WENO-Z5 path removes radial pressure before
the radius-weighted average, restores its arithmetic face value, and updates
the paired pressure scratch flux consistently.

The resolved settings are reported at startup. Runtime controls are documented
in `docs/input.md`.

## Verification completed on the integrated source

| Test | Coverage | Result |
| --- | --- | --- |
| Reconstruction unit test | WENO-Z5, TENO5, TENO6, CPU | Passed. Measured spatial orders were 5.005, 4.998, and 5.998. |
| CUDA build | LLF-WENO-Z5, Cartesian all-fluid | The final-source `compute_rhs.cpp` and both ordinary and paired-pressure WENO kernels compiled for `sm_89`, and the complete executable linked. A clean one-step GPU smoke passed before the final index-based R-Z axis audit. The later incremental RDC executable was not accepted as a runtime result because it reported a stale AMReX device symbol in `CNS::initData`; a clean rebuild is required before the next CUDA run. The WENO-Z5, TENO5, and TENO6 reconstruction header test also compiled for `sm_89`. |
| Cartesian Euler MMS | LLF-WENO-Z5, all-fluid | At `128^2`, one complete SSPRK stage was bitwise identical with fallback on and off. The final Cell SHA was `06877683...fd61aeb`. |
| Cartesian Euler MMS | LLF-TENO5 and LLF-TENO6, all-fluid | Smooth fallback-on/off comparisons were bitwise identical. |
| R-Z Euler MMS | LLF-WENO-Z5, all-fluid | Eight-step `16^2` fallback-on/off results were bitwise identical. |
| R-Z Euler MMS | LLF-TENO5, all-fluid | The `64^2`, 128-step results were bitwise identical for every conserved variable. Evidence is under `exm/numerics/mms_rz_navier_stokes/results/fallback_teno5_rz_mms_20260807`. |
| R-Z Euler MMS | LLF-TENO6, all-fluid | The `64^2`, 128-step results were bitwise identical. Evidence is under `exm/numerics/mms_rz_navier_stokes/fallback_ab_teno6_20260807`. |
| Cartesian Navier-Stokes MMS smoke | LLF-WENO-Z5 and LLF-TENO5, all-fluid, fourth-order centred viscous flux | Both current-source schemes completed two SSPRK(4,3) steps with finite conserved variables and positive density and total energy. Evidence is under `exm/numerics/mms_cartesian/results/fallback_ns_smoke_20260807`. |
| Shu-Osher | LLF-WENO-Z5, LLF-TENO5, and LLF-TENO6, Cartesian all-fluid | WENO-Z5 and TENO5 completed the standard `t=1.8` run at `N=256`. All three schemes completed short on/off A/B runs and produced different final fields, confirming activation at the shock. Evidence is under `exm/numerics/shuosher/fallback_ab_20260807`. |
| Stationary normal shock | LLF-WENO-Z5 and LLF-TENO5, R-Z all-fluid | Both on and off runs completed without non-finite values. The final fields differed for each scheme, confirming activation in R-Z. The TENO5 evidence is under `exm/numerics/rz_planar_shock/results/fallback_teno5_planar_shock_20260807`. |
| Two-level underexpanded jet smoke | LLF-WENO-Z5, R-Z all-fluid, AMR reflux | On and off runs completed ten coarse steps with two levels and 36 reflux calls. All reported thermodynamic fields remained finite and positive. The final fields differed. Evidence is under `exm/underexpanded_jet/2d/npr3_lev3_euler_axisym/results/fallback_amr_smoke_20260807` and the corresponding `_off_` directory. |
| Smooth curved-wall MMS | LLF-WENO-Z5, Cartesian pure shared-GP | Default on/off results were bitwise identical after two steps. Forcing both thresholds to zero changed the fluid solution while leaving the geometry markers unchanged, which confirms that the GP fallback branch is reachable. |
| Pure shared-GP smoke | LLF-WENO-Z5 and LLF-TENO5, Cartesian | Both completed with `family=pure_shared_gp_full_cartesian` and `cut_control_compiled=0`. |
| Legacy shared-GP underexpanded jet smoke | LLF-TENO5, R-Z | Fixed-final-time on/off runs completed with finite positive states and unchanged geometry markers. The fluid fields differed, which confirms activation in the combined R-Z, GP, and TENO path. This legacy shared-image-point case is integration evidence, not frozen BI-CWLS qualification. Evidence is under `exm/underexpanded_jet/2d/npr3_ibm_solid/results/fallback_universal_teno5_rz_pure_gp_20260807`. |
| Non-WENO build check | AFD-HLLC-WENO-Z5 | Built and ran one step. The optional manifest dispatch did not require a fallback member from non-WENO drivers. |

These checks establish integration, smooth-field dormancy, strong-discontinuity
activation, short-run stability, and AMR execution. They do not establish a
positivity theorem, strict monotonicity, optimal threshold values, or a
conservation-convergence rate. The frozen Mach-4 circle, PM expansion,
attached compression corner, inverse-MOC, and long production gates remain
separate acceptance tests and must be reported independently.

The frozen annular BI-CWLS R-Z configuration currently rejects TENO5 and
TENO6 before time advancement because those schemes do not implement the
required paired annular-pressure capability. Their pure-GP device kernels do
compile, but this existing fail-closed compatibility gap must not be bypassed
to obtain a runtime result.
