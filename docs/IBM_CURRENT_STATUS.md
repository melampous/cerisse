# Current Pure GP-IBM Status

Status date: 2026-07-22. This document is the authoritative short status for
the current working tree. It supersedes every status report that described the
production candidate as a cut-cell, embedded-boundary, cut-control-volume, or
aggregate method.

## Method scope

The only production IBM considered here is a pure sharp-interface ghost-point
/ ghost-cell immersed boundary method:

- one shared state at each real solid-side ghost-cell centre;
- cell-average-aware, BI-constrained weighted least squares (BI-CWLS);
- curvature-compatible normal-pressure closure on smooth Euler-slip walls;
- one-sided fluid entropy extension and ideal-gas EOS closure;
- one solid-side ghost-cell layer in the currently tested configuration;
- full Cartesian cell volumes and full Cartesian face areas;
- characteristic LLF flux splitting with WENO-Z5 reconstruction;
- conservative shared-Cartesian-face positivity limiting;
- SSPRK(4,3), with every forward-Euler bracket audited;
- stage-local characteristic physical-boundary ghost filling;
- second-order BI surface recovery for local pressure, Cp, and pressure force.

The evolved semi-discrete operator is

\[
\frac{d\bar U_i}{dt}
=-\frac{1}{|C_i|}\sum_{f\subset\partial C_i}|A_f|\widehat F_f,
\qquad |C_i|=\prod_d\Delta x_d,
\]

where a stencil crossing the immersed boundary obtains its solid-side data
from the shared GP extension. No volume fraction, area fraction, cut-control
residual, physical wall-segment replacement flux, aggregate state, or
redistribution enters this production operator.

The runtime manifest must print
`family=pure_shared_gp_full_cartesian` and `cut_control_compiled=0`. A pure-GP
run that requests a cut-control input terminates immediately.

## Current-source correction status

The current Git commit is `dc35646d8eb90aaa362a3292a9d702e30c9d5af5` and
the tracked dirty-diff SHA-256 is
`604c6323ebf00368df42bbb9bdcf9da39441070d8af41c8a483709e39c4f5a64`.
Content-addressed source IDs are case-specific because they hash `src/`
together with each case's compile-time files: circle `5b45451c364d120e`,
PM/compression `ad0e0c9c747a119b`, inverse-MOC `84db3efb00517a42`, and
Gate 1 isentropic vortex `74d8125a74751cbc`; the Gate-3 shocked-nozzle
executable uses source ID `c92f3840c8626937`.
They are not Git commits and their difference does not indicate a different
solver core. The CUDA `sm_89` binaries compile without the experimental
cut-control macro; symbol audits find no cut-control update operator in the
pure inverse-MOC or Gate 1 executable.

The SSPRK positivity regression was traced to an incorrect admissibility
definition. A 10 K primitive-conversion safeguard had been treated as part of
the Euler invariant domain. The production limiter now again tests

\[
\rho>\rho_{\min},\qquad
\rho e>p_{\min}/(\gamma-1),
\]

and does not use `CLIP_TEMPERATURE_MIN` to trigger or parameterize conservative
flux limiting. This agrees with the archived source implementation and with
the convex admissible set used by the limiter.

The corrected Mach-4 circle run now passes the full long-time gate at
`t*=55.56`. The detached bow shock has `Delta/R=0.5496841171`, the final
four-point relative plateau span is `9.14e-8`, and the density/pressure
location span is `0.1001h`. The limiter activated in 74 Forward-Euler
brackets only over `0.531<=t*<=0.663`; no later activation, global-theta
fallback, all-low final fallback, or rejected step occurred. Shared-face
conservation remained below `1.76e-15` relative error.

The current and older frozen trajectories are not globally bitwise identical
after limiting activates because the Euler admissibility definition changed.
Their forebody field, bow-shock front, stand-off distance and forebody BI
pressure agree at approximately `1e-12`--`1e-13`; differences are confined to
the unsteady aft wake and rear-half wall pressure. Limiter thresholds must not
be tuned to reproduce an old nonlinear wake trajectory.

## Evidence ledger

The following table separates current-source evidence from frozen pure-GP
evidence. Historical hybrid GP/cut-control results are excluded.

| Capability | Status | Evidence and limitation |
|---|---|---|
| Pure full-Cartesian production path | PASS, current source | compile guard, runtime manifest, binary audit, and fail-closed input test |
| Positivity trigger semantics | PASS, current source | Euler density/internal-energy criterion restored; former premature step-11 limiter activation removed |
| Mach-4 circle short fully discrete evolution | PASS, current source | former step-151 failure crossed; shared-face conservation remains at roundoff |
| Mach-4 circle steady bow-shock/positivity regression | PASS, current source | `t*=55.56`, `Delta/R=0.5496841171`, last-four-point span `9.14e-8`; 74 limiter brackets confined to `t*<=0.663`, no retry/fallback; this binary used surface e1, so surface-e2 Cp is not certified by this run |
| Source-free smooth circular-wall Euler solution | **Gate 1 PASS, current source** | two phases and four asymptotic grids: all-fluid L2 1.96--2.09, first-wall-layer L2 1.90--2.04, BI pressure 2.25--2.27, reconstructed wall-normal mass flux 2.56--2.68; 984/984 FE brackets have `theta=1`, and 1723/1723 GP fits are quadratic; Linf phase oscillation is documented |
| Variable-curvature smooth-wall pressure | CONDITIONAL PASS, frozen pure GP | approximately second-order pressure evidence exists; the worst reported phase/orientation sequence is slightly below 1.9 and the exact current-source matrix remains to be rerun |
| PM expansion | PASS, current source | flow-angle error `+0.000858 deg`, pressure-ratio error `+0.006%`; no limiter event and no plateau cell more than 1% below theory |
| Single compression corner | PASS, current source | shock-angle error `-0.051865 deg`, pressure-ratio error `-0.003%`; no limiter event; all 42 fitted shock columns accepted |
| Sequential double compression | CONDITIONAL PASS, frozen pure GP | two pressure plateaux within `0.064%`; interaction location error about `4h`, without a completed grid-convergence sequence |
| No-solid near-vacuum double rarefaction | admissibility PASS, accuracy FAIL | classical ab-initio receding-flow overheating is a bulk LLF-WENO/finite-volume limitation, not an IBM defect |
| Stage-local characteristic inlet/outlet | operator evidence exists; pure-GP integration rerun pending | the characteristic algebra is usable, but hybrid coupling results cannot certify the pure full-cell executable |
| BI surface pressure and pressure force | PASS for Gate 0, current source | inverse-MOC two-phase maximum `Cp` L1 `1.33e-4`; BI pressure-force errors `3.64e-5/1.42e-5`; relative force phase spread `2.22e-5`; this is one-grid non-regression, not a convergence certificate |
| Smooth internal nozzle | **Gate 2 PASS, current source** | exact continuous-curvature inverse-MOC solution, three grids/two phases: no limiter or BI-CWLS fallback; BI pressure, impermeability, mass flow, total enthalpy, entropy, Mach and pressure force converge; boundary-location effect is below fine-grid uncertainty; local velocity norms are not uniformly second order and remain explicitly documented |
| Shocked nozzle | **Gate 3 CONDITIONAL PASS, current source** | fine-grid engineering physics, RH relations, BI load, admissibility and outlet-location sensitivity pass; original pilot retains an absolute `1e-8` mirror-symmetry FAIL and the outlet-pressure criterion has little margin |
| AMR, R-Z, NS heat flux/skin friction, general 3-D STL, FSI | NOT CERTIFIED | separate future qualification programmes |

## Force and conservation semantics

For pure full-cell GP-IBM:

- local wall pressure and Cp are obtained from BI surface recovery;
- physical integrated pressure force is the BI pressure quadrature;
- Cartesian stencil-crossing momentum exchange is a numerical diagnostic;
- the two forces are not required to be identical at finite resolution;
- true curved-domain machine conservation is not claimed;
- no global source, state clipping, pressure clipping, or a posteriori force
  correction may hide the discrepancy.

The required engineering diagnostics are

\[
\epsilon_{\dot m}=
\frac{\max_x\dot m(x)-\min_x\dot m(x)}{\overline{\dot m}},
\qquad
\epsilon_F=
\frac{|F_{\rm BI}-F_{\rm Cartesian}|}{\max(|F_{\rm BI}|,F_{\rm ref})}.
\]

They must decrease under grid refinement and remain below the uncertainty
budget of the intended calculation. Failure of that gate is a limitation of
the pure method; it is not authorization to reactivate cut-cell geometry.

## Fastest recovery and production route

### Gate 0: restore current-source non-regression

1. **Flow/positivity PASS:** the `N=320` Mach-4 circle reached `t*=55.56`
   with source ID `5b45451c364d120e`; bow-shock existence, symmetry,
   stand-off history, limiter location/time, minimum theta, and shared-face
   conservation are documented in
   `docs/ibm_pure_gp_circle_gate_20260722.md`. This run used surface e1; its
   Cp output is not the production-target surface-e2 gate.
2. **PASS:** the current-source PM expansion and single compression corner
   reproduce the analytic turning, shock-angle, and thermodynamic relations.
3. **PASS:** two current-source inverse-MOC phases reproduce the MOC field,
   BI pressure, BI-integrated pressure force, and physical-section mass flow.
   The complete evidence is in
   `docs/ibm_pure_gp_gate0_20260722.md`.

Any loss of the detached circle bow shock or a material stand-off regression
stops all unrelated development.

### Gate 1: close smooth second-order accuracy

**PASS, current source.** The source-free circular-wall exact solution was
run on `N=64,96,128,192,288` in the axis-aligned phase and a normalized
`(0.37h,0.23h)` phase. Formal fits use `N=96,128,192,288`.

- all-fluid complete-state L2 orders are 1.96--2.09;
- first-wall-layer complete-state L2 orders are 1.90--2.04;
- fixed-physical-band L2 orders are 1.96--2.23;
- BI pressure L2 orders are 2.25--2.27;
- reconstructed BI `rho un` L1/L2/Linf orders are 2.56--2.68;
- all 984 SSPRK forward-Euler brackets have `theta=1`;
- all 1723 GP targets use quadratic BI-CWLS, with no linear fallback.

Some global and axis-phase first-layer Linf fits remain 1.73--1.93 because
the maximum-error cell changes with grid phase; the finest-pair first-layer
orders recover above 2.13. Gate 1 therefore certifies approximately
second-order L1/L2 accuracy and asymptotic Linf recovery, not a uniform strict
Linf theorem. Full evidence and provenance are in
`docs/ibm_pure_gp_gate1_20260722.md`.

This gate certifies the volume GP extension and surface observation operator
separately. Increasing ghost layers or polynomial degree is not permitted
unless this frozen operator fails a diagnosed moment-consistency test.

### Gate 2: qualify a smooth nozzle

**PASS, current source.** The exact continuous-curvature inverse-MOC nozzle
was run at `256 x 64`, `384 x 96`, and `576 x 144` in two normalized Cartesian
grid phases. The exact two-dimensional source-free Euler solution, rather
than a quasi-one-dimensional area-Mach relation, provides the reference.

- all 4,864 GP targets use quadratic BI-CWLS, with no fallback;
- every explicitly logged FE bracket remains high order with `theta=1`, and
  the quiet coarse runs contain no unconditional limiter-event message;
- two-phase-envelope pressure and first-wall-layer pressure errors decrease
  with adjacent orders `1.66--1.76` and `1.59--1.70`;
- BI `Cp` L1 reaches `2.43e-5`, while worst nonendpoint Linf is `3.13e-4`;
- non-cancelling BI mass exchange/reference mass flow decreases from
  `9.91e-4` to `1.45e-4`, with finest-pair order `1.96`;
- section mass-flow spread reaches `1.10e-4`; entropy L2 remains about second
  order and total-enthalpy/Mach errors decrease;
- finest BI pressure-force error is at most `1.01e-5`, with phase spread
  `1.76e-6` relative to the exact force;
- moving the inlet and outlet by 12 unchanged cells changes overlapping-domain
  pressure by only `1.31e-8` in relative L2.

This is an engineering qualification, not a second formal-order proof. The
inverse-MOC sequence has nonuniform local velocity orders, including
approximately `0.78--1.02` in the first wall layer; Gate 1 remains the formal
smooth complete-state order evidence. Full evidence is in
`docs/ibm_pure_gp_gate2_20260722.md`.

### Gate 3: qualify a shocked nozzle

**CONDITIONAL PASS, current source.** A prescribed-back-pressure normal shock
was run in a continuous-curvature divergent nozzle at `384 x 64`, `576 x 96`
and `768 x 128`, with two normalized Cartesian phases. Both fine-grid phases
pass every predeclared physical/admissibility/reconstruction check:

- late shock stations are `6.51148` and `6.51009`, with phase spread
  `1.390e-3 = 0.089h`;
- RH mass/momentum errors are approximately `2.6%/1.5%`, pressure-ratio error
  is at most `0.37%`, downstream-Mach error at most `2.28%`, and
  stagnation-pressure-loss error at most `1.18%`;
- off-shock mass-flow spread is `0.045--0.055%` and absolute reconstructed BI
  leakage is `0.049--0.055%` of inlet flow;
- fine-grid BI pressure-force phase spread is `0.01435%` and its error against
  the quasi-one-dimensional engineering value is below `0.09%`;
- physical-section axial momentum-flux change and BI pressure force agree
  within `0.015--0.046%` across the six runs;
- all 407,380 matrix FE brackets remain high order with `theta=1`, and all
  6,912 GP states use quadratic BI-CWLS.

Moving the outlet from `x=12` to `x=14` at fixed `h` changes the common-time
late mean shock station by `7.89e-5`, BI force by `2.16e-5` relative, and the
common-domain pressure by `4.9609e-4` relative L2. The last value passes the
predeclared `5e-4` gate with little margin.

The original fine pilot retains a strict absolute mirror-error FAIL:
`u_Linf=5.98e-8`, `v_Linf=3.38e-8` versus `1e-8`. The errors are `O(1e-9)`
in L2 and lie in the post-shock core rather than the GP layer; an independent
same-resolution extended-domain run gives velocity mirror Linf below
`2.43e-9`. This prevents an unconditional software-verification PASS but does
not overturn the engineering shocked-nozzle result. Complete evidence is in
`docs/ibm_pure_gp_gate3_20260722.md`.

Passing Gates 0--3 supports only this claim:

> Restricted-production pure GP-IBM for two-dimensional Cartesian, uniform
> level-0, fixed smooth walls, single-component ideal-gas Euler-slip flow,
> LLF-WENO-Z5 and SSPRK(4,3), with BI-integrated pressure loads.

No claim for AMR, R-Z SRP, viscous heat transfer, general 3-D STL, or moving
solids follows from that result.

## Maturity assessment

For the restricted two-dimensional Cartesian Euler-slip target, the core
spatial method is substantially built and has credible pure-GP evidence. The
remaining risk is qualification and integration rather than a demonstrated
need to redesign the IBM. Status is recorded only by gates:

- Gate 0: **PASS**; circle, PM expansion, single compression, and two-phase
  inverse-MOC current-source regressions are closed;
- Gate 1: **PASS**; the current-source two-phase source-free curved-wall
  sequence closes complete-state, BI pressure, impermeability, limiter
  transparency, and no-fallback requirements;
- Gate 2: **PASS**; the current-source exact inverse-MOC three-grid/two-phase
  matrix closes the restricted smooth-nozzle engineering gate;
- Gate 3: **CONDITIONAL PASS**; fine-grid shocked-nozzle physics and
  outlet-location sensitivity pass, while the strict pilot mirror criterion
  and narrow outlet-pressure margin remain explicit qualifications.

This is a restricted-production pure GP-IBM candidate for smooth and
stationary-shock two-dimensional Cartesian Euler-slip nozzles at the
demonstrated resolution. It is not a blanket certification for other
shock topologies or for the broader SRP application.

For the eventual SRP target, the distance is much larger because R-Z metric
conservation/axis compatibility, nozzle-lip feature semantics, AMR/reflux,
viscous/thermal closure, and three-dimensional geometry are still independent
uncertified programmes. The existing two-dimensional result must not be used
to claim that general SRP readiness is more than an early-to-mid-stage effort.

## Superseded evidence

Cut-control geometry, aggregate evolution, true-domain machine conservation,
aggregate positivity, physical-subface restart, and equality of cut-control
and BI forces are hybrid-method results. They may remain archived for research
history, but they are not evidence for the method defined in this document and
must not enter its validation matrix.
