# R-Z axis micro-diagnostics, 2026-08-08

## Scope

These are small, level-0 diagnostics for the current coordinate-aware R-Z
operators.  The all-fluid tests are 16 by 16 unless noted otherwise.  They do
not compile IBM, EB, cut-cell, cut-control or aggregate flow updates.  A
separate pure shared-GP/full-cell micro-gate is cited below for the GP parity
check.

The uniform state is

\[
  \rho=1,\qquad p=1,\qquad u_z=1,\qquad u_r=0.
\]

The R-Z radial index is `i`; `ring_i=0` is centred at
\(r=\Delta r/2\), not at the mathematical axis.

## Uniform free-stream results

The default-off `cns.rz_stage_rhs_diagnostics=1` instrument captures the
actual accepted face flux and paired-pressure companion.  For the current
source-free Euler micro-gate it reports

\[
 R_{m_r}=R_{\mathrm{metric},F-P}+R_{\partial p/\partial r}+R_z.
\]

The current paired operator does **not** evaluate a separate Euler
\(D_r(rp)/r\) and pointwise \(p/r\) source.  Consequently those two old-form
terms must not be presented as the active discretization.

| scheme | one RHS, worst conserved residual | SSPRK(4,3), all four RHS calls | split closure |
|---|---:|---:|---:|
| LLF-WENO-Z5, q=2 | \(1.13\times10^{-13}\) | \(1.13\times10^{-13}\) | \(6.31\times10^{-30}\) or zero |
| AFD-HLLC-WENO-Z5 | zero in the initial RHS | \(1.42\times10^{-14}\) | zero |
| Skew4 central | \(2.84\times10^{-14}\) | \(4.26\times10^{-14}\) | zero |
| Skew4 + JST | \(2.84\times10^{-14}\) | \(4.26\times10^{-14}\) | zero |

All four schemes pass the \(5\times10^{-13}\) gate.  No scheme produces a
\(10^{-6}\) to \(10^{-4}\) axis residual.

Reproducible stage logs and JSON summaries are under:

- `results/freestream_llfweno_q2_n16_ssprk43_postgate/`
- `results/freestream_afdhllc_wenoz5_n16_ssprk43/`
- `results/freestream_skew4_central_final_ssprk43_n16/`
- `results/freestream_skew4_jst_final_ssprk43_n16/`

`analyze_stage_rhs_log.py` parses the first-four-ring stage records and fails
if the state residual or split closure exceeds the requested tolerance.

## Smooth pressure moment

For \(p=2+0.25r^2\), the exact radial-momentum RHS is \(-0.5r\).

| scheme | axis error, Nr=16 | Nr=32 | Nr=64 |
|---|---:|---:|---:|
| LLF-WENO-Z5 q=2 | \(2.477\times10^{-10}\) | \(1.968\times10^{-12}\) | \(1.421\times10^{-14}\) |
| AFD-HLLC-WENO-Z5 | \(3.553\times10^{-15}\) | \(7.105\times10^{-15}\) | \(1.421\times10^{-14}\) |
| Skew4 central | \(3.553\times10^{-15}\) | \(7.105\times10^{-15}\) | \(1.421\times10^{-14}\) |
| Skew4 + JST | \(3.553\times10^{-15}\) | \(7.105\times10^{-15}\) | \(1.421\times10^{-14}\) |

Thus the basic axis pressure moment is correct; there is no observed
second-order mixed-form residual.  Skew-JST does add an energy damping term in
this non-uniform pressure field.  Its first-ring magnitude for Nr=16/32/64 is
\(1.2766\times10^{-4}\), \(1.5958\times10^{-5}\), and
\(1.9947\times10^{-6}\), an eightfold reduction per refinement (the expected
third-order background-JST behaviour).  Skew central has no such term.

## Axis parity and vorticity interpretation

With the analytic regular field \(u_r=2r, u_z=1\), the physical-axis ghosts
satisfy, exactly,

- even: \(\rho(-r)=\rho(r)\), \(p(-r)=p(r)\), \(u_z(-r)=u_z(r)\);
- odd: \(u_r(-r)=-u_r(r)\).

For example, the three negative ghost values are
`-0.0078125,-0.0234375,-0.0390625`, paired with the positive values
`+0.0078125,+0.0234375,+0.0390625`.  The active first-ring value is nonzero,
so no active-cell \(u_r=0\) clamp is present.  The same irrotational field has
stored `magvort=0` in rings 0--3.  A plotted first-ring vorticity is a value at
\(r=\Delta r/2\), and must not be interpreted as an exact \(r=0\) sample.

Evidence:

- `results/parity_irrotational_ur2r_n16/parity_trace.log`
- `results/parity_irrotational_ur2r_n16/axis_regularity.json`

## Viscous flux / hoop-source split

The independent `../rz_viscous_axis` gate uses \(u_r=ar\), \(a=0.2\),
\(\mu=0.01\).  A hoop-on/off A/B isolates the actual metric viscous-face
divergence and \(-\tau_{\theta\theta}/r\) source.

| ring | viscous face | hoop source | total radial-momentum RHS |
|---:|---:|---:|---:|
| 0 | \(+4.26667\times10^{-2}\) | \(-4.26667\times10^{-2}\) | \(+1.39\times10^{-17}\) |
| 1 | \(+1.42222\times10^{-2}\) | \(-1.42222\times10^{-2}\) | \(-2.43\times10^{-17}\) |
| 2 | \(+8.53333\times10^{-3}\) | \(-8.53333\times10^{-3}\) | \(+1.91\times10^{-17}\) |
| 3 | \(+6.09524\times10^{-3}\) | \(-6.09524\times10^{-3}\) | \(-1.91\times10^{-17}\) |

The point-state viscous balance is therefore correct cell by cell.  Feeding an
annular cell average to the old point operator instead leaves a first-ring
error of \(-3.7926\times10^{-2}\); the annular-aware operator restores
roundoff balance.  Current Run165 inputs use point semantics, so this mixed
annular/point failure is a guard case, not the active defect.

## Confirmed and repaired generic pure-GP parity defect

The generic non-annular `computeAllGPs()` published a changed positive-radius
GP but left the corresponding negative-radius physical ghost at its pre-GP
value.  The dedicated pure shared-GP/full-cell micro-gate measured
post-publication parity defects as large as 133.3 Pa in pressure, 53.3 m/s in
axial velocity and 36.8 m/s in radial odd parity.

The production repair is deliberately narrow.  Immediately after the final
same-level `FillBoundary`, the generic point-GP call reapplies the existing
coordinate parity helper when `rz_point_target < 0`.  It changes only
negative-radius physical coordinate ghosts; it does not change the unique GP,
active state, BI-CWLS reconstruction, thermodynamic extension, full-cell face
flux, positivity limiter, surface recovery, volume, area or force.  The
annular target loop is excluded because its wrapper already performs one final
parity refresh after all target states have been accumulated.

With the production repair active, the same micro-gate reports exactly zero
for all 12 post-GP parity defects (three ghost depths times rho, pressure,
axial velocity and radial-velocity odd parity).  The runtime manifest remains
`pure_shared_gp_full_cartesian`, with cut, EB and aggregate mechanisms all
disabled.  A two-step annular shared-GP A/B was also compared against the
pre-repair executable: every plot variable agrees bit for bit at zero absolute
and relative tolerance.

This defect is relevant to the current Run165 short-domain inputs because they
set `cns.rz_gp_annular_bic=0`.  The annular GP path already performed the final
parity refresh and did not have this particular gap.

Full evidence and reproduction commands:

- `../../../tst/ibm/positivity/rz_axis_shared_face/integration/` relative to
  the repository root;
- `tst/ibm/positivity/rz_axis_shared_face/integration/` from the repository
  root, file `RZ_GENERIC_GP_AXIS_PARITY_MICROGATE_20260808.md`.

## Nonlinear strong-gradient gate

A separate M=10 stationary planar-shock gate used 16 by 64 level-0 cells,
SSPRK(4,3), CFL 0.25, and the same R-Z axis treatment.  LLF-WENO-Z5 and
AFD-HLLC-WENO-Z5 both reached t=0.02 with fallback on or off.  Skew4 central
and Skew4-JST reached t=0.02 only with their complete-face LLF fallback on;
with fallback off, both failed during the first advance with negative internal
energy.  The failed state was radially uniform, so the reported radial argmin
was a reduction tie rather than an axis defect.

Every successful unperturbed run retained radial pressure/density uniformity
at about 8e-14--4e-13, |u_r|/a at about 1e-14--4e-14, and zero shock-front
corrugation.  A one-cell first-four-ring perturbation was damped rather than
amplified by every fallback-enabled scheme.  This gate therefore identifies
Skew's complete-face fallback as necessary shock protection, but it does not
show an R-Z operator spontaneously generating an axis stripe.

Evidence: `../rz_planar_shock/NONLINEAR_AXIS_GATE_20260808.md`.

## Current diagnosis

The tests rule out the following as the basic all-fluid defect: uniform-flow
preservation, physical-axis ghost parity, the paired pressure moment, and the
point-state viscous/hoop balance.  They found and repaired one independent
generic point-GP parity publication bug.  They also show that Skew needs its
local complete-face LLF fallback at a strong shock; this is a shock-robustness
requirement, not evidence of a free-stream or axis-pressure imbalance.  The
remaining Run165 field question is how much of the observed first-ring
structure survives the GP parity repair and spatial refinement once the
physical flow has reached a comparable developed time.
