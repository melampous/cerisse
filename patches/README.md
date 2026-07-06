# Isolated RZ-NS audit patches (session 2026-06-17)

These two patches isolate THIS session's audit fixes from the contaminated working
tree (which mixes prior-session + other-AI uncommitted changes). Both verified with
`git apply --check` against a pristine HEAD worktree (exit 0, individually and together).

## rz_audit_standalone.patch  — STANDALONE, mergeable now
Applies AND compiles on HEAD (strict subset of an already-built+run config; the
trait/idiom are unused-but-valid without the guard). Contents:
 - Rusanov energy-flux bug fix: cp->cv in prims2fluxes (Thermodynamics.h, Thermodynamics2.h).
   rhoet was rho*(cp*T+ekin)=rhoE+p, so (rhoet+P)*u double-counted P. (Rusanov path only;
   WENO path prims2flux already correct.)
 - IBM ghost P/T positivity floor 0.1 -> 1e-3 (ibm_solver_interp.h): true positivity-rescue;
   stops clamping legitimately-low expansion-fan pressures while still catching zero-crossing.
   Regression: srp_70cone NS (CLIP_MINTEMP=FALSE) 400 steps clean, min rho=4e-4, min T=42K.
 - IBM RZ axis face-skip: comment only (ibm_solver.h). Verified skip is correct for the SRP
   center-nozzle geometry (skipped faces are the contour axis-closure, not real surface);
   caveat documented for solid-on-axis-apex geometries.
 - rz_pressure_split capability traits (weno_t=true; riemann/rusanov/skew/centraldif/no_euler
   =false) + detection idiom rz_psplit_capable<> (defaults false for any flux, incl keep_euler_t).
   Inert without the guard below; included as the safe infrastructure.

## rz_pressure_split_guard.patch  — MERGE WITH the other-AI route-A code
One hunk: the runtime guard in compute_rhs.cpp that Aborts if cns.rz_pressure_split is set
with a non-WENO flux. It references `rz_pressure_split` (route-A's variable, currently
uncommitted) so it only COMPILES once route-A lands. Apply together with route-A.

## NOT included (still experimental / not done)
 - route-A pressure split itself (cns.rz_pressure_split) and route-B near-axis dissipation
   (cns.rz_axis_diss): other-AI / experimental, default-off, not final fixes.
 - A2-a well-balanced pressure: direction decided (use WENO face pressure p^f), NOT implemented.
 - K1 RZ viscous hoop terms: not done.
