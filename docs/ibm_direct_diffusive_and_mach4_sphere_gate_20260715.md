# Direct IBM Diffusive Closure and Mach-4 Sphere Gate

Date: 2026-07-15

## Scope

This report closes two validation priorities for the fixed, Cartesian,
single-species ideal-gas IBM path:

1. diagnose the strong-shock `invalid-state` fallback in the direct isothermal
   heat/viscous crossing flux;
2. run a long three-dimensional Mach-4 sphere case with strict stage
   positivity, strict cubic crossing closure, surface loads, Cartesian
   interface fluxes, and an outer-control-volume momentum budget.

The two verdicts are deliberately separate.  A closure can be mathematically
correct and robust while an extremely coarse sphere grid remains inadequate
for local pressure or integrated-load production.

## Direct diffusive fallback

### Failure mechanism

Let the IBM normal point into the fluid and let

\[
  s_f=(\mathbf{x}_f-\mathbf{x}_w)\mathbin{\cdot}\mathbf{n}
\]

be the signed position of a Cartesian crossing-face centre relative to the
wall.  The constrained polynomial is built from admissible real-fluid support
states and the isothermal/no-slip wall constraints.  When `s_f < 0`, its value
at the Cartesian face is a solid-side mathematical extension used to evaluate
a wall-location-consistent gradient.  It is not a thermodynamic state advanced
by the finite-volume method.

The old gate required the extended temperature to exceed a physical
temperature floor on both sides of the wall.  Behind a Mach-4 bow shock, a
cold 300 K wall and hot fluid support can produce a negative backward cubic
extension even though every real-fluid state, transport coefficient, gradient,
stress, and final conservative flux is finite.  The code then replaced a valid
cubic derivative by the lower-order `invalid-state` fallback.

### Correction

The temperature floor is now enforced only for `s_f >= 0`, where the target is
on the physical fluid side.  The following checks remain mandatory on both
sides:

- finite, positive temperature at every real-fluid support;
- finite velocity supports and polynomial weights;
- non-negative reconstructed conductivity and shear viscosity;
- finite velocity and temperature gradients, stress, viscous work, heat flux,
  and final conservative flux.

Transport is reconstructed as `I_f[mu(T)]`, `I_f[lambda(T)]`, and `I_f[xi(T)]`
from the real-fluid support fields.  The code never evaluates a constitutive
law at a negative solid-side extension temperature.

### Diagnostic evidence

| Test | Old path | Corrected path |
| --- | ---: | ---: |
| Mach-4 sphere targets | 312 | 312 |
| Cubic targets | 288 | 312 |
| `face-temperature-floor` invalid targets | 24 | 0 |
| Solid-side targets | 120 | 120 |
| Corrected GPU run | not applicable | 40 steps, all RK stages cubic |
| Corrected CPU run | not applicable | compile, link, one step, clean exit |

The 24 old failures were all solid-side targets.  In the corrected run, a
solid-side extension may lie below zero or below the local floor, while the
minimum real-fluid support temperature remains positive.  The final per-face
CSV records path, reason, signed distance, target temperature/floor, support
temperature, and transport coefficients.  `direct_diffusive_strict_cubic=1`
turns any quadratic, linear, legacy, or invalid selection into a runtime abort.

A smooth 32-cubed full-NS sphere MMS regression reached `t=0.002` with 744/744
cubic targets and no fallback.  The polynomial reproduction test also passes
exhaustive two-dimensional masks and 20,000 sampled three-dimensional cases;
the maximum constrained cubic value and gradient moment errors are of order
`1e-10`.

With the current executable, a one-step run with runtime diagnostics enabled
is bitwise identical to the same run with diagnostics disabled.  A rebuilt
current executable is not bitwise identical to the archived pre-change binary
after 200 GPU steps, despite taking the same closure-order branch.  Most MMS L2
errors are smaller in the rebuilt run, while the near-wall pressure Linf error
changes from about 62.7 Pa to 117.7 Pa on this single 32-cubed grid.  This is
therefore a branch-safety regression and a single-grid accuracy check, not a
claim of bitwise equivalence with an old compiler artifact.  Formal accuracy
continues to rest on the existing multi-grid MMS sequences.

This closes the strong-shock fallback defect.  It does not expand the direct
diffusive production contract beyond a stationary, constant-isothermal,
no-slip Cartesian wall.

## Mach-4 sphere protocol

The local validation uses:

- Mach 4, Reynolds number 10,000, constant transport;
- WENO-Z5 convection and `viscous_t`;
- stationary 300 K isothermal no-slip sphere, diameter 0.2;
- a 112 x 80 x 80 uniform grid, giving `D/h=8`;
- RK3, CFL 0.30, 2,000 steps, approximately 60 convective times;
- direct wall-location-consistent inviscid and diffusive crossing fluxes;
- stage-by-stage fluid-cell positivity and strict cubic diffusive gates;
- surface VTP output every 25 steps and momentum accounting every 10 steps.

The final analysis uses steps 1,000--2,000.  Its independent quantities are:

- pressure and viscous traction reconstructed on the physical STL surface;
- the actual combined conservative flux through every Cartesian fluid/solid
  marker face;
- the outer-domain conservative flux and fluid momentum time derivative;
- a centreline pressure/density-gradient estimate of bow-shock standoff.

The exact stagnation reference is obtained from a normal shock followed by
isentropic deceleration:

\[
  C_{p,0}=1.7917929466, \qquad p_{0,2}/p_\infty=21.0680810.
\]

Modified Newtonian pressure is plotted only as an approximate forebody shape;
it is not treated as an exact surface-pressure solution.

For the whole sphere, modified Newtonian theory gives the pressure-only value
`Cd,p=Cp0/2=0.8958965`; it is not a total-drag theory.  The Singh et al.
general sphere-drag correlation gives `Cd=0.9716921` for `M=4`,
`Re_D=10000`, `gamma=1.4`, and `Tw/Tinf=1`.  That correlation is informed by
broad ballistic-range data but assumes temperature-dependent air viscosity,
whereas this case uses constant transport, so it is reported as context and
not as an exact acceptance value.  Its reported relative L2 error over the
complete experimental comparison set is 7.6%, so close agreement at this one
state is only a total-drag sanity check.  Turner et al. (JCP 500, 112748, 2024)
document low-Reynolds-number transonic/supersonic sphere validation; they do
not supply an exact solution for this high-Reynolds-number bridge case.

## Final sphere results

The completed-run values are recorded in
`IBM/cases/validation/canonical/sphere_M4_Re10000/diagnostics/direct_pressure_D8/analysis/summary.json`.

The run reached step 2,000 at `t=0.00807396`, or `t U_inf/D=56.07`, and
finalized normally in 2,453 s on the local RTX 4090.  The statistical window is
steps 1,000--2,000 (`t U_inf/D=26.10--56.07`).

| Quantity | Result | Interpretation |
| --- | ---: | --- |
| Exact normal-shock stagnation `Cp` | 1.791793 | reference |
| 10-degree cap area-mean `Cp` | 1.794917 | +0.174%, but cancellation dominated |
| 10-degree cap min / max / standard deviation | 1.59608 / 2.07768 / 0.18699 | local pressure fail |
| Nearest stagnation-bin `Cp` | 1.60194 | -10.60% |
| Fitted stagnation `Cp` from first 20 degrees | 1.81499 | +1.295% |
| Forebody azimuthal `Cp` scatter / exact stagnation `Cp` | 7.463% | symmetry/grid-imprint fail |
| Surface pressure `Cd` | 1.191766 | statistically stationary |
| Surface total `Cd` | 1.192598 | statistically stationary |
| Cartesian conservative-interface `Cd` | 1.105618 | statistically stationary |
| Surface/interface drag difference / surface drag | 7.293% | 1% gate fail |
| Mean vector surface/conservative mismatch / CV force | 8.018% | 1% gate fail |
| Mean transverse surface force / drag | 1.813% | coarse-grid symmetry error |
| VTP/solver surface-pressure integration difference | 0.00333% | pass |
| Surface/CV momentum residual, mean / p95 | 3.856% / 3.954% | 1% gate fail |
| Discrete interface momentum residual, maximum | 6.16e-12 | pass |
| Discrete mass / energy residual, maximum | 2.91e-12 / 1.20e-11 | pass |
| Net Cartesian-interface mass flux, maximum | 96.98% of `rho_inf U_inf A_ref` | impermeability fail |
| Measured pressure-density shock consensus | 0.2868 D (2.295 cells) | underresolved |
| Reference shock standoff resolution | 0.0875489 D = 0.700 cell | cannot validate position |
| Final fluid minimum density / pressure / temperature | 0.6410 / 73.76 kPa / 214.8 K | positive and finite |
| Direct diffusive targets | 312/312 cubic, all reason `none` | strict gate pass |
| Surface area error / normal closure | -0.1195% / 2.10e-16 | geometry audit pass |

The surface pressure-drag relative standard deviation is 0.075% and its linear
drift over the final window is 0.237%.  The persistent surface/interface force
difference is therefore not attributable to an unsteady averaging window.

## Interpretation

### Passed

- The 2,000-step strong-shock run completes without a negative or non-finite
  fluid density, pressure, internal energy, or temperature.
- Every selected direct heat/viscous crossing face remains cubic; no corrected
  face returns to a legacy or invalid path.
- Discrete mass, total-energy, and Cartesian interface momentum accounting is
  at roundoff, proving that the combined flux used by the fluid update is
  conservative on this single-level grid.
- The independent VTP integration and solver surface-pressure accounting agree
  closely, validating their area, normal, sign, and pressure conventions.

### Not passed

- `D/h=8` does not resolve the Ambrosio-Wortman/Billig estimate
  `delta/D=0.0875489`: the physical distance is only 0.700 cell.  A numerical
  gradient peak on this grid is a
  robustness indicator, not a quantitative shock-position result.
- Local forebody pressure has persistent Cartesian/STL imprint.  A stagnation
  cap average close to theory is produced partly by cancellation between local
  undershoot and overshoot and must not be reported as pointwise accuracy.
- Physical-surface traction and Cartesian conservative interface force are not
  the same discrete load.  Their mismatch exceeds the 1% production target,
  even though each accounting path is internally consistent.
- The near-roundoff interface/CV budget does not imply a no-penetration wall:
  the signed net Cartesian-interface mass flux is approximately one reference
  freestream mass flux.  The interface reaction is therefore a discrete
  accounting result, not yet a physically qualified body load.
- The transverse resultant is not negligible on the coarse nominally symmetric
  grid.
- At Reynolds number 10,000 this grid cannot resolve the boundary layer.
  Skin friction and heat flux are checked only for finiteness and closure-path
  robustness; they are not accuracy results.

The historical sphere data in the case directory are not a direct A/B
comparison.  They use an adiabatic wall, four AMR levels, finest
`h/D approximately 0.007`, and a much longer averaging interval.  The direct
run is isothermal and has `h/D=0.125`.  A claim that direct IBM improves or
degrades the legacy method requires matched wall physics, geometry, AMR grids,
initial condition, and time window.

## Required next gate

The next production-pressure test is a matched direct-versus-legacy AMR
sequence with at least 3--5 cells across the expected bow-shock standoff and a
substantially finer near-wall pressure grid.  It must report:

1. pointwise and azimuthally averaged `Cp(theta)`, including cap spread rather
   than only cap mean;
2. shock standoff with pressure and density estimators and cell-count
   uncertainty;
3. STL surface force, conservative interface force, and outer-CV momentum
   balance as separate quantities;
4. transverse-force convergence;
5. strict cubic/fallback maps and fluid positivity extrema.

If the surface/interface load mismatch does not fall below 1% under
refinement, the next algorithmic task is a conservative IBM load closure that
uses the same discrete crossing flux for fluid momentum exchange and body
force, with physical-surface pressure retained as a separate diagnostic.

That accounting closure is now implemented, but the D/h8 mass-flux audit adds
a stricter prerequisite: the crossing flux must first satisfy and converge
under a global impermeability gate.  Otherwise the algorithmic task is a true
wall-location-consistent no-penetration flux/cut-control-volume correction,
not merely reporting the equal-and-opposite crossing reaction.

Follow-up status: the composite-AMR Cartesian crossing-flux load and its
restart-safe reflux lifecycle are implemented and documented in
`docs/ibm_amr_composite_crossing_force_20260715.md`.  The statistically
stationary `D/h=32` and `D/h=48` surface/interface convergence results remain
in progress.
