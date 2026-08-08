# Composite AMR Cartesian Crossing-Flux IBM Load

Date: 2026-07-15

## Scope

This implementation defines a discretely conservative Cartesian crossing-flux
reaction for a fixed IBM body in Cartesian, single-species ideal-gas
calculations.  It does not alter the fluid flux, rescale the reconstructed
surface traction, or implement FSI.  Physical load qualification additionally
requires impermeability, stationarity, positivity, and grid convergence.

The physical-surface and conservative-interface loads answer different
questions:

- the STL surface integral diagnoses the accuracy of reconstructed pressure,
  heat flux, and viscous traction at the physical wall;
- the Cartesian crossing-flux integral is the equal-and-opposite momentum
  exchange associated with the numerical flux that actually updates the
  fluid control volumes.

## Discrete definition

For Runge-Kutta stage `s`, level `l`, and Cartesian fluid/solid crossing face
`f`, let `F[l,f,s]` be the combined Euler plus diffusive conservative flux and
`A[l,f]` its Cartesian face area.  The time-integrated interface exchange is

\[
 I_{\Gamma}
 = \sum_s b_s \Delta t
   \sum_l \sum_{f\in\Gamma_l^{\rm owner}}
   \boldsymbol{F}_{l,f,s} A_{l,f}.
\]

The SSPRK stage weights are the final quadrature weights, not the forward
Euler substep sizes.  For SSPRK(3,3), for example, they are
`(1/6, 1/6, 2/3)`.

For a source-free fixed-body problem, the composite fluid balance is

\[
 \frac{\boldsymbol{P}_f^{n+1}-\boldsymbol{P}_f^n
       -\Delta\boldsymbol{P}_{\rm regrid}}{\Delta t}
 + \boldsymbol{\Phi}_{\rm outer}
 + \boldsymbol{F}_{\Gamma}=0,
 \qquad
 \boldsymbol{F}_{\Gamma}=I_{\Gamma}/\Delta t.
\]

The same accounting is applied to mass and total energy as a diagnostic.  An
isothermal wall can exchange energy through heat conduction; a stationary
impermeable wall should drive the net mass exchange toward zero with
refinement.

## Composite AMR ownership

Every active level retains its actual combined face flux when
`ib.momentum_budget=1`.  A cached `makeFineMask` marks coarse cells covered by
the next level.  A crossing flux is counted only from an uncovered fluid-cell
owner, so no coarse/fine representation is counted twice.  Every level adds
its time-integrated contribution to the level-0 budget.

The implementation enforces these topology conditions:

1. `cns.do_reflux=1` for AMR budgets;
2. the fluid and solid cell centers adjacent to a coarse crossing face have
   the same coarse/fine ownership;
3. fine grids remain strictly inside the physical domain because the outer
   boundary accumulator currently uses level-0 physical faces;
4. AMR remap changes in composite mass, momentum, and energy are recorded and
   removed from the time derivative;
5. every boundary face has one valid first hit and one accounting owner.

If any condition fails, the run aborts instead of emitting a nominal load.

## Restart correction

AMReX restores a level through the default constructor, so the parameterized
constructor does not define that level's `FluxRegister`.  Reflux had
historically been disabled, which hid this lifecycle defect.  The first AMR
restart step therefore segfaulted in `FluxRegister::Reflux` after conservative
reflux was enabled.

`post_restart` now rebuilds every nonzero level's register from its restored
BoxArray, DistributionMapping, refinement ratio, and conservative-component
count.  `advance` also aborts explicitly if a required register is absent.
The transient IBM budget is initialized from the accepted synchronized state
before the first restart advance swaps StateData time levels.

The ordinary AMR regrid path was audited separately: AMReX constructs every
new or remade level through the parameterized `CNS` constructor before calling
`init(old)`, so the register is defined against the new BoxArray and
DistributionMapping.  The lifecycle defect was specific to restart's default
constructor path.

## Geometry-classification correction

The 3-D BVH inside/outside test previously rejected an entire ray whenever it
encountered any nearly parallel candidate triangle.  Axis rays through a
sphere therefore fell through to one symmetric oblique ray, which could
double-count a shared edge.  At `D/h=48` this produced one isolated false
solid cell and six crossing faces with no first hit.

Parallel triangles are now skipped because they do not cross the ray.  A ray
is retried only when a proper forward intersection lies near a triangle edge
or vertex.  Three axis and two asymmetric oblique directions provide a
deterministic fallback.

The `D/h=48` audit changed as follows:

| Quantity | Before | After |
| --- | ---: | ---: |
| GP targets | 17499 | 17496 |
| crossing copies | 10830 | 10824 |
| missing first hit | 6 | 0 |
| invalid WLS reconstruction | 6 | 0 |
| coarse-fine missing supports | 0 | 0 |

The six removed crossing copies all surrounded the same false solid cell.

## Verification

The ratio-3 hierarchy has `D/h=48`.  The M=4 sphere correlation used by this
case gives `delta/D=0.087548896`, or 4.202 finest cells across the theoretical
bow-shock stand-off.

### Physical-reference hierarchy

No exact full-Navier--Stokes solution exists for this `M=4`, `Re_D=10000`
cold-isothermal sphere.  Three references therefore have distinct roles:

1. The normal-shock plus isentropic-stagnation relation gives the exact
   calorically-perfect-gas stagnation reference `Cp0=1.7917929466`.
2. The Ambrosio--Wortman/Billig relation gives `delta/D=0.087548896` for a
   continuum sphere.  It is an empirical shock-shape reference, not an exact
   viscous stand-off distance.
3. Modified Newtonian theory gives the windward pressure-only estimate
   `Cd,p=Cp0/2=0.8958964733`.  It must not be compared as if it were the total
   viscous drag.  The general Singh et al. sphere-drag correlation, evaluated
   at `M=4`, `Re_D=10000`, `gamma=1.4`, and `Tw/Tinf=1`, instead gives the
   experimental-data-informed total-drag reference `Cd=0.9716920610`.

The Singh correlation assumes a temperature-dependent air viscosity while
this matched direct/legacy experiment intentionally uses constant transport.
It is therefore context for total drag, not a one-percent acceptance target.
Singh et al. report a 7.6% relative L2 error over the complete experimental
database used to assess the correlation; agreement with its single M=4 value
must be treated as an engineering sanity check rather than a high-fidelity
Navier--Stokes validation datum.
Turner, Seo, and Mittal (JCP 500, 112748, 2024) validate their method on
low-Reynolds-number transonic and supersonic spheres; that paper does not
provide an exact reference for the present `M=4`, `Re_D=10000` case.

### Two-step AMR conservation

| Path | MPI ranks | max momentum residual | max mass residual | max energy residual |
| --- | ---: | ---: | ---: | ---: |
| direct | 4 | 6.706e-10 | 1.175e-11 | 1.616e-09 |
| legacy isothermal | 4 | 9.296e-11 | 1.322e-12 | 5.925e-10 |
| direct | 1 | 7.154e-10 | 6.990e-12 | 1.572e-09 |
| legacy isothermal | 1 | 1.025e-10 | 7.504e-13 | 5.765e-10 |

The direct Cartesian interface drag differs between one and four ranks by at
most `3.84e-16` relative in these steps.  Surface drag differs by about
`1e-14`; the finite-difference CV force differs by at most `3.04e-10` because
the time derivative amplifies reduction-roundoff in total momentum.

### Restart

A four-rank direct run was checkpointed after step 1 and restarted for step 2.
The uninterrupted and restarted step-2 records are identical for time,
composite momentum, all three surface-force components, all three interface
force components, all three CV-force components, and the discrete residual.

### Uniform hierarchy

The same direct path on `max_level=0` completed with maximum relative
residuals `1.422e-9` (momentum), `3.584e-11` (mass), and `7.039e-9` (energy).
This checks the no-mask branch independently of the composite AMR branch.

### Crossing-mode mass-flux ablation

All seven no-slip crossing modes were restarted from the same statistically
stationary direct `D/h=8` checkpoint and advanced for two CPU steps.  The table
is an instantaneous flux-path diagnostic, not a stability or convergence
result.  Mismatch is normalized by the Cartesian interface force and mass flux
by `rho_inf U_inf A_ref`.

| Mode | net mass flux | surface/interface drag mismatch |
| ---: | ---: | ---: |
| 0 | 187.10% | 20.23% |
| 1 | 271.00% | 36.54% |
| 2 | 263.59% | 1.68% |
| 3 | 254.06% | 41.40% |
| 4 | 107.35% | 31.56% |
| 5 | 8.72% | 21.31% |
| 6 | 97.03% | 7.79% |

Mode 5 removes the one-sided two-state LLF pair and reduces leakage by more
than an order of magnitude relative to mode 6, showing that the LLF density
jump term is the dominant source.  Its remaining leakage still fails the 1%
gate, and its force mismatch is worse.  No existing mode satisfies both
requirements, so parameter selection alone cannot close the problem.

### Compile coverage

- GCC, MPI, direct and legacy-isothermal ratio-3 builds: pass.
- NVCC `sm_89`, direct and legacy-isothermal ratio-3 builds including the
  restart lifecycle patch: compile and link pass.

## Production gate

Report these quantities separately:

1. `interface-CV`: discrete conservation error.  The run gate requires every
   sample to be valid and each relative mass, momentum, and energy residual to
   remain below `1e-7`.
2. net `interface_mass_flux/(rho_inf U_inf A_ref)`: fixed-wall global
   impermeability error.  It must remain below 1 percent and decrease with
   refinement; cancellation-free local leakage maps should also be inspected.
3. load stationarity and fluid-state admissibility: drag drift must remain
   below 1 percent over a common physical-time window and all fluid density,
   pressure, and internal energy values must remain positive and finite.
4. `surface-interface`: physical-wall traction reconstruction error.  Surface
   traction qualifies as the integrated production load only below 1 percent.
5. surface `Cp(theta)`, heat flux, and skin friction: local physical-wall
   diagnostics; the Cartesian interface load does not replace them.

If `surface-interface` remains above 1 percent under refinement while
`interface-CV`, impermeability, stationarity, and positivity all pass, use the
Cartesian crossing flux for the integrated fixed-body force and report the
surface discrepancy as unresolved wall-load reconstruction error.  If net
mass leakage does not converge below the threshold, the crossing reaction is
only a discrete accounting quantity and the local wall flux requires a
conservative no-penetration or cut/control-volume correction.

## Remaining limitations

- The long direct/legacy `D/h=32` sequence and the `D/h=48` statistically
  stationary sequence are still required before assigning a production
  uncertainty to the surface pressure and total drag.
- A stationary `D/h=8` direct sphere closes its discrete interface/CV mass,
  momentum, and energy residuals near `1e-11`, but has maximum signed net mass
  leakage `0.970 rho_inf U_inf A_ref`.  Co-refinement must show that this
  quantity vanishes before the crossing reaction can be called a physical
  production load.
- The current closure provides an integrated fixed-body load.  FSI would also
  require a conservative spatial force and torque distribution.
- RZ geometry is rejected because its swept IBM surface measure and metric
  source accounting are not implemented in this budget.
- Fine grids touching a physical boundary are rejected until composite
  fine-level outer-boundary flux accumulation is implemented.
- Budget mode retains all-level face fluxes and performs reductions at every
  RK stage.  It is opt-in; its memory and runtime overhead should be profiled
  before enabling it in large production calculations.

## Physical-reference sources

- F. S. Billig, *Shock-Wave Shapes around Spherical- and
  Cylindrical-Nosed Bodies*, Journal of Spacecraft and Rockets 4 (1967),
  <https://doi.org/10.2514/3.28969>.
- N. Singh et al., *A General Drag Coefficient for Flow over a Sphere*,
  <https://arxiv.org/abs/2012.04813>.
- J. M. Turner, J. H. Seo, and R. Mittal, *A high-order sharp-interface
  immersed boundary solver for high-speed flows*, JCP 500 (2024),
  <https://doi.org/10.1016/j.jcp.2023.112748>.
