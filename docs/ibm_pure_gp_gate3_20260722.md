# Pure GP-IBM Gate 3: shocked-nozzle qualification

Date: 2026-07-22

## Decision

\[
\boxed{\text{Gate 3: CONDITIONAL PASS within the frozen shocked-nozzle scope}}
\]

No numerical-method change was authorized or made during this gate.

The evidence closes the principal grid/phase, shock-jump, admissibility, BI
pressure-load, outlet-location and full-Cartesian ledger checks. The original
fine-grid pilot nevertheless retains its formal `FAIL` label because an
auxiliary absolute mirror-symmetry threshold was exceeded. The final decision
therefore distinguishes the engineering shocked-nozzle gate from that strict
software-symmetry diagnostic.

This decision authorizes a restricted-production candidate for a smooth-wall
shock-containing nozzle at the demonstrated fine resolution. It is not an
unconditional verification PASS. In particular, the outlet-pressure
sensitivity passes with little margin, the coarse/intermediate local shock
states do not all satisfy the three-percent gate, and the strict pilot mirror
criterion remains failed.

## Frozen method and provenance

Every run uses the pure shared-GP, full-Cartesian method:

- one shared state at each real solid-side ghost-cell centre;
- cell-average-aware BI-CWLS with curvature-compatible pressure and a
  one-sided entropy extension;
- one configured solid-side ghost-cell layer;
- full Cartesian cell volumes and full Cartesian face areas;
- LLF-WENO-Z5 and SSPRK(4,3);
- conservative shared-Cartesian-face positivity limiting;
- stage-local characteristic pressure outflow;
- volume extrapolation order 1 and BI surface recovery order 2.

The runtime manifest reports
`family=pure_shared_gp_full_cartesian`, `pure_gp_production=1`, and
`cut_control_compiled=0`. Cut-cell, cut-control, volume/area-fraction,
aggregate, redistribution and AMReX-EB flow updates are excluded.

The common executable is

```text
main2d.gnu.TPROF.CUDA.gate3_shocked_nozzle_llf_wenoz5_i2_e1_s2_a0p6_gl1_vis2_bic_cellavg_sjet_pcurv.srcc92f3840c8626937.ex
```

with SHA-256
`df6cfb07399cf4e85e46a841855c32c58716ade9a70bf27fa1c71b6f7817bcb1`.
Its manifest records Git commit
`dc35646d8eb90aaa362a3292a9d702e30c9d5af5`, source-tree SHA-256
`384996cc17fd21e1ae34c8643bb619b8a92bcf90ee103470b97c19a03bbc9dc4`,
and AMReX `25.12-21-gbd922c6216e0`.

The two fixed 26,001-point wall carriers have SHA-256 values
`beefc3c8365a37359676eae9e94a2560cdade102a0a6dec7897657aad487158e`
and
`a37c23c547b5ec8011c609d4a73f7da56ff06f61dd8e913354285816fa0f38b8`.

## Verification fixture and reference semantics

The symmetric two-dimensional nozzle has straight inlet/outlet buffers and a
slowly divergent continuous-curvature wall. Its half-height changes from
0.5 to 0.7 through a quintic smootherstep over `2 <= x <= 10`. A complete
Mach-2 state enters at `x=0`; the stage-local characteristic outlet imposes

\[
p_b=0.5226372972150315.
\]

An isentropic/normal-shock/isentropic quasi-one-dimensional construction
targets `x_s=6.5` and supplies the initial conservative cell averages. It also
supplies engineering mass-flow, pressure-force and normal-shock reference
quantities:

\[
M_1=2.253632688,\qquad M_2=0.5400935347,
\]

\[
p_1=0.08599140133,\qquad p_2=0.4951961326,
\]

\[
\dot m_{\rm q1D}=0.4057667889643083,
\qquad F_{p,x}^{\rm q1D}=0.10410171417472933.
\]

This construction is not an exact solution of the two-dimensional Euler
equations. It is not used to assign pointwise two-dimensional discretization
error. Rankine-Hugoniot consistency is evaluated from the numerical pre- and
post-shock states, while shock station and BI pressure load are numerical
observables.

## Pilot and protocol provenance

The first calculation used `768 x 128`, `h=0.015625`, `CFL=0.23`, and
`t=30`. All predeclared physical checks passed:

- late mean shock station `6.5114839311`;
- late range `1.5053e-3` and drift `1.8072e-4` per unit time;
- RH mass, momentum and total-enthalpy residuals
  `2.592%, 1.528%, 0.00805%`;
- pressure-ratio, downstream-Mach and stagnation-pressure-ratio errors
  `0.256%, 2.196%, 1.093%`;
- off-shock mass-flow spread `0.05543%`;
- reconstructed absolute BI wall leakage/inlet flow `0.05470%`;
- positive active-fluid density, internal energy and pressure;
- all 90,536 Forward-Euler brackets at `theta=1`;
- all 1,536 GP states quadratic, with no BI-CWLS downgrade;
- normalized Cartesian mass, energy and momentum ledger residuals no larger
  than `5.04e-10`, `3.01e-10`, and `8.39e-12`.

An auxiliary absolute primitive mirror-error threshold of `1e-8` was added to
the analysis before the final result was inspected, but after the run had
started; it was not part of the pre-launch protocol text. The final axial- and
transverse-velocity mirror Linf values were `5.98e-8` and `3.38e-8`, so the
pilot analysis is conservatively retained as `FAIL` under that added
criterion.

The pilot manifest records protocol SHA-256
`f81a2859d4a0ef0d039c2bf2ef7a17cd44c6d7d326fdf99e12b42dbd6b6c1445`,
whereas the later outlet manifest records
`ea3c6574d7ca5f8c6d9059dabdda74625ea3ae493054512f90e18b8185c470fd`
for the clarified protocol file. The separate amendment preserves the timing
and effect of the symmetry criterion; the current protocol file is therefore
not represented as an untouched pre-launch artifact.

The follow-up audit found that the corresponding L2 errors are only
`1.53e-9` and `9.42e-10`. Their maxima lie in the post-shock core, 16--18
cells downstream of the shock and more than 23 cells from either immersed
wall. Surface-e2 mirror discrepancies remain below `4.24e-9`. This evidence
does not identify an IBM wall-reconstruction defect, but it also does not
permit deletion or retroactive relaxation of the recorded criterion.

## Grid and normalized-phase matrix

The independent diagnostic extension uses three grids and two fixed
dimensionless geometry/grid phases:

| Grid | Phase | late mean `x_s` | RH mass | RH momentum | `M_2` error | absolute BI leakage | BI pressure force |
|---|---|---:|---:|---:|---:|---:|---:|
| `384 x 64` | `(0,0)` | `6.513241` | `5.142%` | `2.930%` | `4.422%` | `0.1428%` | `0.1039040` |
| `384 x 64` | `(0.37h,0.23h)` | `6.512969` | `5.133%` | `2.863%` | `4.529%` | `0.0858%` | `0.1039224` |
| `576 x 96` | `(0,0)` | `6.512262` | `3.461%` | `2.001%` | `2.978%` | `0.0786%` | `0.1039704` |
| `576 x 96` | `(0.37h,0.23h)` | `6.510822` | `3.558%` | `1.991%` | `3.184%` | `0.0622%` | `0.1039928` |
| `768 x 128` | `(0,0)` | `6.511484` | `2.592%` | `1.528%` | `2.196%` | `0.0547%` | `0.1040089` |
| `768 x 128` | `(0.37h,0.23h)` | `6.510094` | `2.603%` | `1.493%` | `2.279%` | `0.0494%` | `0.1040238` |

The coarse and intermediate grids do not satisfy every local three-percent
shock-state threshold. Both fine-grid phases do. RH mass and momentum errors
decrease at approximately first order, as expected for local sampling across
an `O(h)` captured discontinuity. No all-field second-order claim is made for
this shock-containing solution.

The fine-grid shock-station phase spread is `1.3903e-3`, or `0.089h`, and is
comparable with the late temporal range. The phase-spread sequence is not
monotone; it is therefore described as unresolved below the combined
shock-location/temporal uncertainty, not as a formal Richardson limit.

Across the complete six-run matrix, all 407,380 Forward-Euler brackets remain
high order with `theta=1`, and all 6,912 GP states remain quadratic. No
global-theta, all-low-order, rejected-step or BI-CWLS fallback event occurs.

## Mass flow, wall pressure and axial load

The off-shock section mass-flow spread is `0.045%--0.055%` on the two finest
grids. The mean mass-flow error relative to the quasi-one-dimensional
engineering reference decreases to `0.0257%` and `0.0219%` in the two fine
phases.

The non-cancelling reconstructed BI wall leakage decreases in both phases.
For phase zero it is `0.1428%, 0.0786%, 0.0547%`; for the shifted phase it is
`0.0858%, 0.0622%, 0.0494%`. This is BI reconstructed impermeability, not a
claim that each Cartesian crossing face has zero mass flux.

The fine-grid BI surface-e2 pressure-force errors relative to the
quasi-one-dimensional engineering value are `0.0892%` and `0.0748%`; the
two-phase BI-force spread is `0.01435%`. The quasi-one-dimensional value is
not an exact two-dimensional load reference.

An independent physical-section quadrature gives

\[
\int_{A_e}(\rho u^2+p)\,dy-
\int_{A_i}(\rho u^2+p)\,dy .
\]

Its difference from the BI surface-e2 pressure force is `0.015%--0.046%`
across the six runs. The inlet/outlet total-energy-flux difference is
`0.019%--0.057%`. These are finite-time engineering momentum/energy checks;
they are not machine-precision true-curved-domain conservation statements.

The surface-e2 postprocessor and the runtime BI pressure-force budget agree
to `1.8e-8--9.3e-8` relative. Thus the reported physical wall load has one
unambiguous BI pressure source.

## Pure GP conservation semantics

The active full-Cartesian finite-volume ledger closes to reduction tolerance.
It conserves the staircase set of active Cartesian cells, not the exact curved
fluid volume. The Cartesian crossing momentum exchange differs from the BI
physical pressure force by approximately `4.1%` on both fine phases and by as
much as `16.5%` on the coarse phase-zero grid. This quantity remains a
numerical full-cell diagnostic and is not interpreted as local physical
traction or used as a Gate-3 wall-load error.

Similarly, Cartesian crossing mass exchange is not the physical BI wall
leakage. No cut-cell correction, global source, clipping or a posteriori force
correction is used.

## Reflection-symmetry diagnosis

Phase-zero `384 x 64` and `576 x 96` solutions retain primitive mirror Linf
errors of `O(1e-12)`. The `768 x 128` phase-zero run develops an oscillatory
post-shock mode of `O(1e-8)` after approximately `t=24`; the initial field is
symmetric to `O(1e-16)`. The non-monotone grid behavior and the location away
from the immersed wall rule out a demonstrated fixed IBM truncation bias.
They do not prove that the mode is harmless under every decomposition or
longer integration. This is the reason an unconditional Gate-3 verification
PASS is not assigned before the complete review.

## Outlet-location sensitivity

The frozen `896 x 128`, `x_out=14` calculation reached `t=30` with return
code zero. It preserves `h`, wall geometry over the common domain, inlet
state, back pressure, CFL and all numerical methods. The outlet protocol was
fixed shortly after process launch, at approximately
`t=0.15`, but before any outlet-sensitivity result was inspected. It is a
pre-result protocol, not a strict pre-launch protocol; this is an additional
reason for retaining a conditional rather than unconditional rating.

The pre-result thresholds were:

- all fine-grid physical/admissibility/reconstruction checks pass;
- common-time late mean shock-station change no larger than `1.51e-3`;
- BI pressure-force relative change no larger than `2e-4`;
- common-domain pressure relative L2 difference no larger than `5e-4`;
- off-shock mass-flow spread below 1% and absolute BI leakage below 0.5%.

All conditions pass:

| Diagnostic | Result | Limit |
|---|---:|---:|
| common-time late mean shock-station change | `7.8921e-5` | `1.51e-3` |
| BI pressure-force relative change | `2.1570e-5` | `2.0e-4` |
| common-domain pressure relative L2 | `4.9609e-4` | `5.0e-4` |
| off-shock mass-flow spread | `0.08471%` | `1%` |
| absolute reconstructed BI leakage | `0.04947%` | `0.5%` |

The pressure threshold passes by only `3.91e-6`, approximately 0.78% of the
allowed value. This narrow margin is retained as a qualification, not rounded
into a stronger result. The pressure relative L1 is `2.4759e-4`; relative
Linf is `8.426e-3` because the remaining low-amplitude wave field is localized
near the captured shock.

The extended case has late mean shock station `6.5114004`, late range
`1.9186e-3`, and drift `1.9425e-4`. Its RH mass/momentum/total-enthalpy
residuals are `2.638%, 1.602%, 0.00552%`; pressure-ratio, downstream-Mach and
stagnation-pressure-ratio errors are `0.365%, 2.135%, 1.178%`. All remain
inside the predeclared three-percent fine-grid limits.

All 90,536 Forward-Euler brackets remain at `theta=1`; all 1,792 GP states
remain quadratic. The maximum normalized mass, energy and Cartesian momentum
ledger residuals are `6.30e-10`, `4.71e-10`, and `7.80e-12`.

The extended domain also provides an independent same-resolution symmetry
diagnostic with a different x partition. Its axial/transverse velocity mirror
Linf values are `2.43e-9` and `1.34e-9`, and pressure Linf is `9.13e-11`.
Thus the pilot's `O(1e-8)` mode is not a reproducible fixed truncation bias of
the immersed wall. The original pilot criterion remains failed because the
calculation that was tested did exceed it.

## Final assessment

The physical/engineering Gate-3 subgate passes at `768 x 128`, corresponding
to 64 cells across the inlet height and approximately 80 cells across the
shock section. The local RH errors decrease at approximately first order and
only this finest resolution satisfies every three-percent shock-state gate.
Comparable shocked-nozzle production calculations should therefore use at
least this local resolution until a separate case demonstrates otherwise.

The strict software pilot remains `FAIL` under its added `1e-8` absolute
mirror criterion. Since the discrepancy is `O(1e-9)` in L2, is located in the
post-shock core rather than the GP layer, disappears in the extended-domain
same-resolution run, and does not affect the BI pressure symmetry or integral
observables at their error budgets, it is retained as a documented numerical
qualification rather than interpreted as a wall-closure failure.

The resulting bounded capability statement is:

> For the frozen LLF-WENO-Z5/SSPRK(4,3) configuration, the pure shared-GP
> full-Cartesian method is a restricted-production candidate for
> two-dimensional, uniform level-0, fixed smooth-wall, ideal-gas Euler-slip
> nozzles containing a stationary internal shock, provided the demonstrated
> local resolution and the stated mass-flow, wall-leakage, BI-load,
> admissibility and outlet-sensitivity monitors are retained.

## Capability boundary

No Gate-3 evidence changes the certified method or expands its geometry and
physics scope. The result applies only to two-dimensional Cartesian, uniform
level-0, fixed smooth walls,
single-component ideal-gas Euler-slip flow, LLF-WENO-Z5 and SSPRK(4,3), with
BI surface-e2 pressure loads. It does not certify AMR, R-Z, Navier-Stokes wall
shear/heat flux, general three-dimensional STL, moving solids, TENO, AFD or
HLLC.

## Evidence

- [Gate-3 protocol](../IBM/cases/validation/canonical/shocked_nozzle_2d/GATE3_PROTOCOL.md)
- [symmetry amendment](../IBM/cases/validation/canonical/shocked_nozzle_2d/GATE3_PROTOCOL_AMENDMENT_1.md)
- [outlet protocol](../IBM/cases/validation/canonical/shocked_nozzle_2d/GATE3_OUTLET_PROTOCOL.md)
- [pilot analysis](../IBM/cases/validation/canonical/shocked_nozzle_2d/runs/gate3_pilot/analysis/summary.json)
- [pilot symmetry audit](../IBM/cases/validation/canonical/shocked_nozzle_2d/runs/gate3_pilot/analysis/symmetry_audit.json)
- [matrix analysis](../IBM/cases/validation/canonical/shocked_nozzle_2d/runs/gate3_symmetry_matrix/matrix_analysis.json)
- [matrix diagnostics](../IBM/cases/validation/canonical/shocked_nozzle_2d/runs/gate3_symmetry_matrix/gate3_grid_phase_diagnostics.png)
- [outlet-sensitivity analysis](../IBM/cases/validation/canonical/shocked_nozzle_2d/runs/gate3_outlet_extended/analysis/outlet_sensitivity.json)
- [outlet-sensitivity figure](../IBM/cases/validation/canonical/shocked_nozzle_2d/runs/gate3_outlet_extended/analysis/gate3_outlet_sensitivity.png)
