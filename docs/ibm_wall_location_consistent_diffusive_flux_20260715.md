# Wall-location-consistent IBM heat and viscous flux closure

Date: 2026-07-15

## Scope

This work targets the derivative amplification observed after coupled
Navier--Stokes evolution near a stationary immersed wall.  The supported
production contract is currently:

- Cartesian geometry;
- one ideal-gas species and five conservative equations;
- stationary, constant-isothermal, no-slip wall;
- `ghost_layers=1`, visibility-filtered support, and quadratic IBM support
  blocks;
- `viscous_t` or `diffusiveheat_t` with the direct closure capability enabled.

It does not yet cover RZ, moving FSI walls, adiabatic walls, species diffusion,
LES closures, or PelePhysics transport.

## Failure mechanism

Reconstructing an immersed-wall primitive value with an `O(h^2)` error is not
enough for a second-order conservative RHS.  A first-layer cell applies

\[
  R_i^{(d)}=-\frac{F_{i+1/2}^{(d)}-F_{i-1/2}^{(d)}}{h_d}.
\]

If either the crossing-face or the neighbouring fluid-face heat/stress flux
has only `O(h^2)` error, division by `h_d` can produce an `O(h)` first-layer
residual.  The same mechanism amplifies small evolved errors in temperature
and tangential velocity because Fourier heat flux and viscous stress
differentiate those fields before the conservative divergence is taken.

The required local contract is therefore an `O(h^3)` face-gradient error (or
better) on both faces that enter the first-layer divergence.

## Discrete closure

### Geometry-time reconstruction

For every actual Cartesian solid/fluid crossing, initialization records the
raw first hit, owning element, oriented wall normal, visible fluid supports,
and a separate wall closure point.  The raw first hit remains the topology and
accounting key.  A problem may override only the closure point and frame for
formal analytic-geometry verification.

A constrained cubic least-squares polynomial is formed from real-fluid
support values.  Relative to wall point `x_w`, its basis contains the scaled
normal linear mode and all Cartesian quadratic and cubic monomials.  This
enforces, at `x_w`,

\[
  q(x_w)=q_w, \qquad \nabla_t q(x_w)=0.
\]

The explicit fallback ladder is cubic, quadratic, normal-linear, then legacy.
Rank, condition number, finite weights, and an L1 amplification bound are all
checked before a level is accepted.

### Finite-volume target consistency

The support values are point-indexed primitive samples, while the existing
fourth-order viscous shell uses the finite-volume coefficients in
`src/rhs/diff_ops.H`.  The direct target functionals therefore reproduce the
same modified-equation values for every fitted cubic:

\[
  I_f q = q_f-\frac{h_d^2}{24}q_{dd},
  \qquad
  D_{f,g}q = q_g\vert_f-\frac{h_d^2}{24}q_{ddg}.
\]

This applies to face values, the normal derivative, and every tangential
derivative.  It avoids subtracting a pointwise crossing flux from a
finite-volume neighbouring flux and then amplifying that mismatch by `1/h`.

### Runtime conservative flux

At each RK stage, the cached geometry-only functionals are applied to current
temperature and velocity.  A unique crossing-face flux is then assembled as

\[
  F_{\rho}=F_{\rho}^{I},
  \quad
  F_{\rho u_j}=F_{\rho u_j}^{I}-\tau_{jd},
  \quad
  F_{\rho E}=F_{\rho E}^{I}-u_j\tau_{jd}-\lambda T_{,d}.
\]

The stress is

\[
  \tau_{ij}=\mu\left(u_{i,j}+u_{j,i}
  -\frac{2}{3}\delta_{ij}\nabla\cdot u\right)
  +\xi\delta_{ij}\nabla\cdot u.
\]

Face velocity, all velocity derivatives, temperature derivative, stress,
heat flux, and viscous work use the same target face and polynomial.

A Cartesian crossing-face centre is not necessarily in the physical fluid.
When it lies inside the solid, the constrained value is only a polynomial
extension used to obtain a wall-location-consistent derivative.  In
particular, a cold wall behind a strong shock can give a negative solid-side
temperature extension without implying a negative thermodynamic state.
Positivity floors are therefore enforced on fluid-side target values and on
every real-fluid support state, but not on a solid-side extension value.
Both sides still require finite weights, gradients, stresses, and final
conservative fluxes; reconstructed conductivity and viscosity must also be
non-negative.  This distinction prevents a mathematically valid derivative
from being replaced by a lower-order legacy flux.

For temperature-dependent transport, the code reconstructs the coefficient
field itself:

\[
  \mu_f=I_f[\mu(T)],\quad
  \lambda_f=I_f[\lambda(T)],\quad
  \xi_f=I_f[\xi(T)].
\]

Using `mu(I_f[T])` would omit nonlinear `O(h^2)` terms such as
`mu''(T)|grad(T)|^2`, which can again become first order after divergence.

The inviscid face flux is captured before the dense diffusive operator runs.
After that operator, only the selected crossing face is overwritten.  The
combined conservative flux remains in the ordinary divergence and AMR reflux
path.  Non-crossing near-wall fluid/fluid faces use a common real-fluid cubic
finite-volume shell; no wall state is written into the shared primitive
MultiFab.

## Verification

All convergence fits below use `N=32,48,64,96`, an exact analytic circular
closure point and normal, and the unchanged faceted polygon for topology and
visibility.

| Test | Quantity/band | fitted order |
|---|---|---:|
| heat-only, `t=0` RHS | energy, first fluid layer | 2.192 |
| heat-only, `t=0` RHS | energy, all fluid | 2.729 |
| heat-only, fixed `t=0.002` | temperature, first fluid layer | 2.19 |
| heat-only, fixed `t=0.002` | temperature, all fluid | 2.73 |
| viscous-only, `t=0` RHS | energy, first fluid layer | 2.008 |
| viscous-only, `t=0` RHS | energy, all fluid | 1.998 |
| viscous-only, `t=0` RHS | momentum, first fluid layer | roundoff |
| viscous-only, fixed `t=0.002` | temperature, first fluid layer | 2.01 |
| viscous-only, fixed `t=0.002` | pressure, all fluid | 2.00 |
| full NS WENO-Z5, fixed `t=0.002` | temperature, first fluid layer | 2.17 |
| full NS WENO-Z5, fixed `t=0.002` | velocity, first fluid layer | 2.37--2.38 |
| full NS WENO-Z5, fixed `t=0.002` | all-fluid primitive variables | 2.56--2.78 |

In the viscous-only evolved test, first-layer velocity error remains
approximately `1e-13`; it is not converted into a growing wall-shear error.
This is the direct regression for the tangential-velocity amplification
failure.  The evolved heat and full-NS tests provide the corresponding
temperature regression.

The `t=0` full-NS test retains approximately fifth order in the bulk.  Its
first-layer momentum fit is `1.85--1.97`, while the fixed-time coupled solution
is above second order.  This identifies the remaining local defect as the
convective/pressure IBM shell, not the heat/stress closure.

Additional gates passed:

- constrained polynomial reproduction in 2-D exhaustive masks and sampled
  3-D masks;
- maximum cubic face value/gradient moment errors
  `5.67e-11 / 5.56e-10`;
- CPU GCC 2-D and 3-D full-NS compilation;
- CUDA `sm_89` compilation and an actual RTX 4090 run;
- strict GPU run: `104/104` crossing faces cubic, with zero quadratic,
  linear, legacy, or invalid fallbacks;
- one-rank versus two-rank 64-square run: all five conservative RHS arrays
  bitwise identical;
- non-IBM full-NS build;
- default production run with no direct-flux diagnostic output or diagnostic
  atomic increments.

Primary result files:

- `IBM/cases/validation/mms_ibm/convergence_rhs_shell_heat_target_fv_exact_wall_v1_32_48_64_96.md`
- `IBM/cases/validation/mms_ibm/convergence_fixed_heat_target_fv_exact_wall_v1_32_48_64_96.md`
- `IBM/cases/validation/mms_ibm/convergence_rhs_shell_viscous_target_fv_exact_wall_v7_32_48_64_96.md`
- `IBM/cases/validation/mms_ibm/convergence_fixed_viscous_target_fv_exact_wall_v7_32_48_64_96.md`
- `IBM/cases/validation/mms_ibm/convergence_rhs_shell_ns_wenoz5_target_fv_exact_wall_v1_32_48_64_96.md`
- `IBM/cases/validation/mms_ibm/convergence_fixed_ns_wenoz5_target_fv_exact_wall_v1_32_48_64_96.md`

## Runtime gates

The following optional `ib.*` controls are available:

- `direct_diffusive_diag=1`: report the path count after each RHS call;
- `direct_diffusive_strict_order=1`: abort on linear, legacy, or invalid
  crossing closure;
- `direct_diffusive_strict_cubic=1`: additionally abort on a quadratic
  crossing closure;
- `direct_diffusive_runtime_csv=prefix`: write the final per-face path,
  failure reason, signed wall-normal target distance, target temperature,
  local floor, minimum support temperature, conductivity, and viscosity to
  `prefix_direct_diffusive_runtime_L<level>_R<rank>.csv`;
- `direct_face_admissibility_floor_fraction=f`: set the local primitive-state
  admissibility floor used by the wall-consistent inviscid and diffusive
  crossing closures. The default is `f=0.5`, and the accepted range is
  `0<f<=1`. Pressure, temperature, and (where present) `gamma-1` are compared
  with `f` times the minimum admissible value in the local fluid support. A
  non-admissible inviscid candidate is convexly scaled toward its admissible
  local reference state. A non-admissible physical diffusive target takes the
  explicit `invalid-state` fallback; a finite solid-side polynomial extension
  is not itself treated as a thermodynamic state. This is a face-state
  safeguard, not a conservative-state clip or a proof of positivity for the
  complete update.

Counter reset and atomic path counting are compiled into the feature but do
not execute unless a diagnostic or strict gate is enabled.

The following optional `cns.*` validation controls audit the time integrator:

- `stage_positivity=1`: after every completed Runge--Kutta stage, audit every
  fluid IBM cell for finite conservative components, positive density, and
  positive internal-energy density. Constant-gamma ideal-gas builds also
  audit pressure;
- `stage_positivity_verbose=1`: print the global minimum density,
  internal-energy density, and pressure for every stage. It requires
  `stage_positivity=1`;
- stage auditing never clips the state and is therefore incompatible with
  `soft_positivity=1`. A failure reports the scheme, stage, time step, global
  minima, and representative bad cells before aborting.

## Interpretation and remaining work

The direct heat/viscous subsystem now meets the formal second-order gate for a
smooth, accurately located wall and remains second order after coupled time
evolution.  This does not prove that a fixed low-resolution STL has second-order
geometry.  A fixed geometric location error `delta` can enter a second
derivative as `delta/h^2`; formal fluid refinement must therefore use an
analytic surface or geometry whose location and normal errors are below the
fluid truncation error.

The next technical priorities are:

1. add an adiabatic constrained basis and moving-wall work consistently;
2. consume real PelePhysics coefficient FABs before advertising that path;
3. add species enthalpy/diffusion fluxes for multicomponent flow;
4. validate the same gate in 3-D with an analytic sphere and with a fixed,
   sufficiently accurate STL;
5. quantify and cache the runtime-selected real-fluid cubic shell stencil in
   3-D, where its current marker search is more expensive than the crossing
   flux itself;
6. finish the independent convective/pressure near-wall shell gate, which is
   now the limiting part of the full NS `t=0` first-layer residual.
