# R--Z First-Axis-Ring Consistency Audit

## Follow-up clarification (2026-07-31)

The pressure-flux/geometric-source conclusion in this report remains valid:
the paired radial pressure flux and geometric source are not the source of the
observed first-axis-ring defect.

A later source-free conical-flow audit exposed a separate data-semantics error
in the axial flux of radial momentum. The annular conservative average had
been replaced by its recovered cell-centre value before axial WENO
reconstruction. For an axis-regular odd field \(m_r=a r\), the first-ring
centre value and annular average differ by the fixed factor \(3/4\), leaving an
\(O(h)\) residual in the axial divergence of radial momentum. The opt-in
production correction retains centre-state recovery for radial point fluxes
and BI-CWLS, but constructs axial flux support from the stored annular
conservative averages.

The corrected four-grid, two-phase matrix gives fitted first-axis-ring
radial-momentum residual orders of \(2.0299\) and \(2.0296\). The complete
follow-up evidence and intended-use qualification are recorded in
`../../../IBM/cases/validation/canonical/RZ_PURE_GP_QUALIFICATION_20260731.md`.
Consequently, the original statement below that no *pressure/source* change
was justified remains correct; it must not be read as excluding the later,
independent axial-flux semantics correction.

## Scope

This validation-only audit tests the frozen full-cell R--Z finite-volume
operator used by the pure GP-IBM path. It does not modify the production
pressure flux, geometric source, WENO reconstruction, BI-CWLS extension, or
SSPRK(4,3) update. No embedded-boundary, cut-cell, cut-control-volume, or
aggregate update is used.

The audit addresses the previously observed approximately first-order radial
momentum residual in the first axis-adjacent ring of the RZ-2b wall--axis
junction MMS.

## Pressure-flux/source decomposition

For the radial momentum pressure contribution, four validation routes were
evaluated:

\[
\begin{aligned}
A&=\text{analytic face pressure}+\text{analytic cell source integral},\\
B&=\text{analytic face pressure}+\text{production geometric source},\\
C&=\text{production face pressure}+\text{analytic cell source integral},\\
D&=\text{production face pressure}+\text{production geometric source}.
\end{aligned}
\]

The pressure profiles were constant, axis-regular quadratic, and axis-regular
quartic functions of radius. Grids \(N=16,32,64,128\) were used.

### Results

- Constant pressure: all four routes are exactly balanced to roundoff.
- Quadratic pressure: all four routes reproduce the radial pressure balance
  to roundoff.
- Quartic pressure: the first-ring \(L_\infty\) errors converge at order
  \(3.00\) for routes B, C, and D.
- The complete production pair D decreases from
  \(3.4516\times10^{-6}\) at \(N=16\) to
  \(6.7552\times10^{-9}\) at \(N=128\).

These results reject the hypothesis that the current production
pressure-flux/geometric-source pair is intrinsically first order in the first
axis-adjacent ring.

Raw data and figures are in
`abcd_matrix_20260730/summary.json`,
`abcd_matrix_20260730/summary.csv`, and
`abcd_matrix_20260730/abcd_first_ring_convergence.{png,pdf}`.

## Current-source periodic R--Z MMS

The current executable was independently rerun with periodic axial boundaries
and the complete axis-regular, non-zero-\(u_r\) manufactured solution. The
production annular-average, point-flux, pressure-consistency, and centre-state
recovery options were explicitly enabled; all historical pressure-split,
well-balanced, axis-dissipation, and shock-only experimental options were
disabled.

For radial momentum, the adjacent-grid orders in the first ring are

\[
2.73,\qquad 2.81,\qquad 2.40,
\]

and those over the first three rings are

\[
2.84,\qquad 2.97,\qquad 2.99.
\]

The four-grid fitted first-ring orders are \(2.66\) for radial momentum and
\(1.96\) for total energy. The current periodic-axis production operator is
therefore second-order consistent or better in the first axis-adjacent ring.

Data and provenance are in
`../rz_mms_full/axis_consistency_current_20260730/`.

## RZ-2b spatial provenance

The archived RZ-2b wall--axis junction matrix was reanalysed without rerunning
or changing the numerical method.

For both geometry phases and all four grids, the maximum first-ring radial
momentum residual occurs at index

\[
(i_r,i_z)=(0,N_z-1),
\]

the upper axial physical-boundary cell. Its fitted first-ring \(L_2\) order is
\(1.014\) for phase 0 and \(1.017\) for the shifted phase.

In contrast, the wall--axis patch radial-momentum \(L_2\) orders are
\(2.096\) and \(2.215\), while the wall--axis first-wall-adjacent orders are
\(1.874\) and \(2.040\). The approximately first-order aggregate is therefore
caused by the axial physical-boundary closure included in the first-ring norm,
not by the axis pressure/source pair or the wall--axis GP reconstruction.

The machine-readable provenance is stored in
`../../../IBM/cases/validation/rz_wall_axis_mms/axis_region_provenance_20260730.json`.

## Decision

\[
\boxed{
\text{R--Z first-axis-ring production operator: PASS}
}
\]

\[
\boxed{
\text{RZ-2b wall--axis GP closure: PASS}
}
\]

No production numerical change is justified by this audit. The physical
axial-boundary closure remains a separate local verification limitation and
must not be attributed to the GP-IBM or the R--Z geometric source.
