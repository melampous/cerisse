# IBM Pressure Closures

This document defines the pressure policies available to the built-in static
GPIBM wall models. The implementation is in `src/ibm/ibm_containers.h` and
`src/ibm/ibm_walltypes.h`.

Scope: single-fluid static IBM. The historical default remains
`zero_gradient`; existing parameter structs and custom wall models are not
silently changed.

## Policies

```cpp
enum class ibm_pressure_closure_t : int {
    zero_gradient,
    fluid_extrapolation,
    prescribed_gradient,
    shock_aware_fluid_extrapolation,
    euler_slip_analytic_curvature,
    navier_stokes_noslip_viscous_compatibility
};
```

For surface-to-fluid image-point distances `0 < x1 < x2`, the smooth
two-point extrapolate is

```text
p_s = (x2*p1 - x1*p2)/(x2-x1).
```

The prescribed-gradient reconstruction inverts the nonuniform one-sided
Lagrange derivative through the surface, IP1, and IP2. Given `g=dp/dn`,

```text
g  = cs*p_s + c1*p1 + c2*p2,
cs = -(x1+x2)/(x1*x2),
c1 = x2/[x1*(x2-x1)],
c2 = -x1/[x2*(x2-x1)].
```

All higher-order reads are protected by `if constexpr` and `n_valid`. With the
default policy and `extrap_order=1`, the implementation retains the historical
nearest-image-point behavior.

## Shock-Aware Extrapolation

The Stage-2 sensor removes the reversible acoustic part of an IP1-to-IP2
compression:

```text
dp    = p1-p2
drho  = rho1-rho2
cbar2 = 0.5*(c1^2+c2^2)
S     = clip(2*max(dp-cbar2*drho,0)/(p1+p2),0,2).
```

It is zero unless both `dp` and `drho` are positive. A cubic smoothstep blends
the unlimited extrapolate to IP1 between `pressure_sensor_low` and
`pressure_sensor_high`; the validated defaults are `0.03` and `0.10`.
Nonfinite or nonpositive candidates also fall back deterministically.

This is a local protection for a compression that crosses IP1--IP2. It is not
an entropy-stability proof and cannot detect a shock that does not intersect
that segment. The thresholds are engineering parameters validated for the
single-component ideal-gas path.

Set `ib.surf_pressure_diag=1` to print trigger statistics. Surface VTP output
contains `PressureShockSensor` and `PressureFallbackFraction`. The reporting
path is optional because it traverses managed surface arrays on the host.

## Two-Dimensional Sharp-Feature Regularization

A polygon corner is not a smooth boundary point: its normal is discontinuous,
so a nominally second-order normal polynomial through IP1 and IP2 has no local
smoothness argument at the vertex. For opt-in two-dimensional sharp polygons,
the shock-aware pressure closure can therefore impose a compact geometric
fallback floor without changing the owner edge or averaging its normal.

For polygon vertex `i`, the turning angle is

```text
theta_i = atan2(abs(t_(i-1) cross t_i), t_(i-1) dot t_i).
```

Vertices above `sharp_feature_angle_deg` seed a shortest-path calculation on
the cyclic edge graph. If `d` is the resulting arc distance from a surface or
ghost projection to the nearest sharp vertex and
`R = sharp_feature_radius_cells * min(dx)`, the fallback floor is

```text
eta = clip(1-d/R, 0, 1),
w_f = eta^2 * (3-2*eta),
w   = max(w_f, w_shock),
p_s = (1-w)*p_unlimited + w*p_IP1.
```

The band is `O(h)`: it contracts in physical units under grid refinement and
is an exact no-op outside `R`. Constant pressure is preserved exactly. The
feature does not modify slip velocity, temperature, species, ghost
extrapolation, or the convective flux; a broader all-variable fallback was
tested and rejected because it increased the post-fan pressure error.

Selection is compile-time and currently requires a built-in shock-aware wall:

```cpp
struct ibmparm_t {
    static constexpr ibm_pressure_closure_t pressure_closure =
        ibm_pressure_closure_t::shock_aware_fluid_extrapolation;
    static constexpr bool sharp_feature_pressure_limiter = true;
    static constexpr Real sharp_feature_angle_deg = Real(10.0);
    static constexpr Real sharp_feature_radius_cells = Real(10.0);
};
```

Choose the angle threshold below the smallest intentional corner and above
the per-vertex turn of any tessellated smooth contour. The radius is a
regularization width in cells, not a physical rounding radius. This policy is
not enabled by default and is not a replacement for a resolved rounded lip.
When the optional members are omitted, the dormant angle and radius defaults
are `30 deg` and `3 h`; production cases should set them explicitly when the
feature is enabled. Both BVH and CGAL geometry backends use the same
topological annotation and a scale-relative degenerate-edge tolerance.

## Analytic Euler-Slip Curvature

For a fixed impermeable inviscid wall, define the IBM normal from solid into
fluid and

```text
B_ab = t_a dot (grad n) t_b.
```

The instantaneous normal Euler compatibility condition is

```text
dp/dn = rho_s * u_t^T B u_t.
```

Consequently, `zero_gradient` is not a high-order-compatible extension on a
smooth curved slip wall when the tangential velocity and curvature are both
nonzero. It remains exact for a plane wall, and at a stagnation point where
`u_t=0`, but those special cases do not justify using it over an entire curved
body.

Under this convention, an exterior convex circle or sphere has principal
curvature `+1/R`; reversing the normal reverses both `B` and `dp/dn`. The flow
may be unsteady, but the wall must be fixed. A moving wall, normal body force,
normal manufactured source, or viscous no-slip wall requires additional terms.

The problem parameter type must provide a device-safe callback:

```cpp
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ibm_shape_operator_t
wall_shape_operator(
    const Array1D<Real, 0, AMREX_SPACEDIM-1>& xyz,
    const Array1D<Real, 0, AMREX_SPACEDIM-1>& normal,
    const Array1D<Real, 0, AMREX_SPACEDIM-1>& tangent1,
    const Array1D<Real, 0, AMREX_SPACEDIM-1>& tangent2);
```

The returned `k11,k12,k22` are the symmetric components of `B` in the supplied
tangent basis; two-dimensional problems use only `k11`. Return `valid=0`
outside the callback's declared smooth region.

For a smooth analytic or spline-defined two-dimensional wall, the topology
carrier and the closure geometry are deliberately separated. The polygon
continues to provide inside/outside classification and an owning element, but
the problem type may provide the complete closure geometry through three
device-safe callbacks:

```cpp
smooth_wall_boundary_intercept(raw_point, boundary_intercept);
smooth_wall_local_frame(boundary_intercept, frame);
wall_shape_operator(boundary_intercept, normal, tangent1, tangent2);
```

The boundary intercept, frame, and shape operator must describe the same
smooth wall point. The contract is restricted to fixed two-dimensional smooth
walls and must not be applied across corners or other feature edges. It has
been verified for axis-aligned, translated-grid-phase, and rotated analytic
ellipses. It is not a curvature-recovery certificate for a piecewise-linear
polygon or a faceted STL surface. The BI-constrained visibility functionals
currently require the native BVH backend; selecting this reconstruction with
the CGAL backend is rejected at compile time.

The sign follows directly from the implemented convention. For the exterior
of a convex circle, the solid-to-fluid normal is `n=e_r`, so
`t dot (grad n) t=+1/R`. The normal component of the steady Euler momentum
equation then gives `dp/dn=+rho*u_t^2/R`. A case callback that uses the opposite
normal must reverse the shape operator as well.

The BI-constrained GP reconstruction and the current-stage fluid wall trace
are deliberately separate functionals. The primary constrained functional
alone selects the support, polynomial degree, condition number, and weight-L1
acceptance. The trace supplies only `rho_s` and `u_t` to the compatibility
datum, so changing the pressure policy cannot silently change the support used
for temperature or velocity.

The high-order curvature candidate is accepted only when the complete ghost
primitive state is finite and has admissible pressure and temperature. If it
is rejected on a non-smooth startup or shock stencil, the complete state is
reduced together: pressure, temperature, and tangential velocity are copied
from the nearest visible active-fluid support actually used by the accepted
functional, and the wall-normal velocity is reflected for the stationary slip
wall. This avoids combining a lower-order pressure with high-order temperature
and velocity. A homogeneous-gradient BIC state, the fluid trace, and finally
the established EOS floor are progressively more exceptional fallbacks; all
paths are exposed by validation diagnostics. The smooth exact-solution test
uses the high-order path for every GP.

A validation-only forced-rejection test executes this complete-state reduction
through the production branch. In an assertion-enabled `N=64`
nonconstant-entropy ellipse run, all 68 ghost cells select fallback code 1 and
satisfy exact assertions that pressure, temperature, and tangential velocity
come from the same nearest admissible fluid support while the normal velocity
is reflected. This verifies branch execution and primitive-state consistency;
it is not a long-time accuracy test of the reduced state.

Selection is compile-time:

```cpp
struct ibmparm_t {
    static constexpr int interp_order = 2;
    static constexpr int extrap_order = 1;
    static constexpr int interp_order_surf = 2;
    static constexpr int extrap_order_surf = 1;
    static constexpr int ghost_layers = 1;
    static constexpr int support_visibility_mode = 2;
    static constexpr bool expanded_bic_support = true;
    static constexpr ibm_shared_gp_reconstruction_t shared_gp_reconstruction =
        ibm_shared_gp_reconstruction_t::bi_constrained;
    static constexpr ibm_pressure_closure_t pressure_closure =
        ibm_pressure_closure_t::euler_slip_analytic_curvature;
    // wall_shape_operator(...) must also be defined.
};
```

`extrap_order=1` is intentional here. The BI-constrained polynomial evaluates
the real ghost-cell centre directly; it does not perform the historical
image-point-to-GP extrapolation. Higher surface reconstruction order is a
separate diagnostic choice and is not required by this volume closure.

The curvature policy is rejected at compile time for built-in no-slip walls.
It does not infer curvature from STL facets and is not combined with the shock
fallback. Use `shock_aware_fluid_extrapolation` for shock-dominated generic
STL calculations.

This ghost-state admissibility logic is distinct from admissibility of evolved
active-fluid cell averages. A strict SSPRK-stage audit of the Mach-4 circle
finds negative internal energy in a first-layer active-fluid cell for both this
curvature policy and the matched `zero_gradient` policy. This is an observed
cell-average failure, not merely the absence of a positivity theorem. It is not
unique to the curvature closure, but the available comparison does not isolate
the no-solid bulk operator from the IBM near-wall coupling. No positivity
limiter is enabled by selecting the pressure closure.

## Entropy-Based Complete-State Euler Extension

For the single-component calorically perfect gas, the smooth-wall Euler-slip
path reconstructs the entropy proxy

```text
sigma = log(p) - gamma*log(rho)
```

from the visible fluid support. Its auxiliary normal extension is homogeneous,
while pressure uses the normal-momentum compatibility condition above. The
ghost-cell density and temperature are then recovered from `(p,sigma)` through
the ideal-gas EOS, and total energy is recomputed from the resulting primitive
state. Pressure, density, temperature, velocity, and total energy therefore
belong to one thermodynamically consistent state.

Only wall impermeability is a physical Euler boundary condition. The
homogeneous normal extension of `sigma` is an auxiliary smooth-extension
choice; it is not a viscous adiabatic-wall condition and does not imply
`dT/dn=0`. The nonconstant-entropy ellipse verification has a wall-normal
temperature gradient as large as 0.15185 and exercises this distinction.

The entropy extension is currently restricted to the built-in stationary
Euler-slip wall, a single ideal gas with constant gamma, BI-constrained shared
ghost cells, and the smooth-wall geometry contract. If the high-order pressure
or entropy state is nonfinite, nonpositive, rank-deficient, ill-conditioned,
or otherwise inadmissible, the complete primitive state is reduced together;
independent variable-by-variable fallback is not permitted.
The current order test has entropy varying tangentially but satisfies
`d(sigma)/dn=0`; a solution class with nonzero wall-normal entropy gradient is
not covered by that certificate.

## Navier-Stokes No-Slip Compatibility

At a fixed, stationary no-slip wall with no wall-normal body force or source,
the wall velocity and convective acceleration vanish. The normal momentum
equation therefore requires

```text
dp/dn = n_i d(tau_ij)/dx_j,
tau_ij = mu*(du_i/dx_j + du_j/dx_i - 2*theta*delta_ij/3)
         + xi*theta*delta_ij,
theta = div(u).
```

The implementation fits Cartesian velocity, `mu(T)`, and `xi(T)` around the
first image point with a quadratic least-squares basis. It evaluates the full
variable-property expression

```text
div(tau)_i = mu*laplacian(u_i)
           + grad(mu)_j*(du_i/dx_j + du_j/dx_i)
           + (xi + mu/3)*d(theta)/dx_i
           + (d(xi)/dx_i - 2*d(mu)/(3*dx_i))*theta.
```

The estimate is sampled one image-point distance from the wall. Its truncation
error is `O(h)` because it contains second velocity derivatives; multiplication
by the `O(h)` wall-reconstruction coefficient still gives an `O(h^2)` pressure
value. This is the accuracy target of the closure, not a claim that
`n dot div(tau)` itself is second-order pointwise.

Selection is compile-time and requires all four IBM reconstruction orders to
be at least two:

```cpp
struct ibmparm_t {
    static constexpr int interp_order = 2;
    static constexpr int extrap_order = 2;
    static constexpr int interp_order_surf = 2;
    static constexpr int extrap_order_surf = 2;
    static constexpr int ghost_layers = 1;
    static constexpr ibm_pressure_closure_t pressure_closure =
        ibm_pressure_closure_t::navier_stokes_noslip_viscous_compatibility;
};
```

The built-in slip walls reject this policy at compile time. It currently
supports fixed Cartesian 2-D/3-D, single-fluid transport. RZ and PelePhysics
return an invalid compatibility estimate rather than using inconsistent
coefficients. The historical complete fluid `3^D` candidate block remains the
first estimator. In the direct BI-CWLS shared-GP path, an unavailable complete
block is retried on the already accepted, visibility-filtered quadratic
BI-CWLS support. The retry uses the same support selected for the constrained
Dirichlet/Neumann functionals and rejects rank-deficient or ill-conditioned
normal matrices using `support_visibility_condition_max`; it does not search a
second stencil or create a face-local ghost state.

If neither quadratic estimator is admissible, the BI-CWLS pressure extension
retains its explicit homogeneous-Neumann fallback. Nonfinite temperature or
transport coefficients, insufficient rank, excessive condition number, or a
failed factorization all take this deterministic fallback. No pressure
clipping is applied.

Moving walls, wall acceleration, body forces, and manufactured momentum sources
add terms to the compatibility equation and are outside this policy. A source
MMS is valid only when its wall-normal momentum source vanishes, as in the
provided compatibility MMS.

The current production-path verification is documented in
`docs/ibm_noslip_pressure_compatibility_20260727.md`. In the two-dimensional
stationary isothermal no-slip MMS, the fixed-time complete state and BI
pressure are second order or better on both tested grid phases. This is not yet
a production certificate for wall heat flux or skin friction: the fixed-time
surface `qwall` and `tau_mag` fits remain about 1.56 and 1.37, respectively.
Adiabatic or moving walls, RZ, three-dimensional faceted geometry, reacting
transport, shock/no-slip interaction, and AMR remain outside the verified
scope.

## Verification Snapshot

The latest formal gate is a globally smooth, nonconstant-entropy and
nonconstant-temperature Euler MMS around a variable-curvature ellipse. The
same field is run with and without the solid over `N=64,96,144,216`, including
a `(0.37h,0.23h)` grid-phase shift and a 27-degree wall rotation.

- The no-solid fixed-time L2 orders are 3.09--3.75, so the source and bulk
  discretization do not explain the wall-adjacent rate.
- Physical-edge-length-weighted pressure at the actual closure BI converges at
  L2 orders 1.892, 1.973, and 1.928 for the axis, shifted, and rotated series;
  the maximum-error phase/orientation envelope has order 1.957.
- Fixed-time first-fluid-layer pressure converges at order 2.44--2.75, while
  velocity components converge at order 1.84--2.08.
- The strict first-layer semi-discrete gate is not passed: momentum residual
  orders are 1.67--1.85, and density/energy orders are 1.64--1.98.
- Fresh CPU and CUDA sm_89 fixed-time fields agree to at most `7.99e-15`
  absolute and `7.2e-16` relative L2.

The detailed equations, norms, phase envelope, fallback audit, and provenance
are in
`IBM/cases/validation/mms_ibm/diagnostics/variable_curvature_ellipse_entropy/VARIABLE_CURVATURE_ELLIPSE_ENTROPY_MMS_REPORT_20260718.md`.

The earlier concentric-vortex test below remains supporting evidence for the
curvature-pressure mechanism, but it is not the current formal order gate.

The earlier two-dimensional curvature gate is a source-free, stationary
isentropic vortex with a concentric radius-0.2 circular slip wall. It is an
exact Euler solution, not an MMS. The exact wall condition is
`dp/dn=rho*u_t^2/R`.

- The former `pzero` path gives about order 0.56 for global pressure and 0.51
  for first-layer pressure, while its first-layer energy residual does not
  converge. This directly verifies the curved-wall incompatibility described
  above.
- With the numerical curvature trace, expanded visible BIC support, and no
  exact ghost-state replacement, the five-grid `N=80,120,160,240,320` t=0
  first-layer residual orders are 2.02 for density, 1.89 for both momentum
  components, and 2.02 for total energy.
- At fixed `t=0.01`, first-layer solution-error orders are 2.66 for density,
  2.42 for both velocity components, 2.65 for pressure, and 2.33 for
  temperature. A normal SSPRK physical step uses the high-order curvature path
  at all 88 GPs on the `N=80` grid.
- CPU and CUDA results from the same final source differ by at most
  `7.99e-15` over all saved conservative and primitive fields.
- A Mach-4 circular-cylinder regression at `D/h=64`, fixed 4,096-segment
  geometry, and `t*=55.56` remains stable and symmetric. It gives
  `Delta/R=0.5498128788`, a 17.59-cell stand-off, and a final-window relative
  span of `1.95e-8`. Billig and Hornung values are empirical/fitted references,
  not exact errors for this computation.

The shock regression also exposes a separate limitation. At CFL 0.30 the
curvature configuration produces `rho*e<0` at SSPRK43 stage 4, which is the
accepted `U^(n+1)` state; the matched `pzero` configuration fails at stage 1,
which feeds the next RHS. Before the end-of-step check was corrected, these
states could evade detection because density and total energy remained
positive, after which primitive conversion applied an internal-energy floor.
The current code rejects `rho*e<=0` unconditionally. The curvature fallback is
therefore not a proof of positivity preservation, and shock-containing nozzle
use remains uncertified.

Earlier component and surface-reconstruction checks remain useful supporting
evidence:

- Algebraic reconstruction, normal reversal, and Stage-2 sensor checks pass.
- Independent centered residual checks recover orders 2, 4, and 6.
- The 3-D manufactured source has wall-normal momentum below `1.5e-11`.
- A fixed 20,480-facet frozen sphere sequence (`N=24,32,48`) gives surface
  pressure L2 orders 1.04 (zero gradient), 1.77 (fluid extrapolation), and
  2.55 (analytic curvature).
- A 4,096-segment 2-D circle evolved to fixed `t=0.002` (`N=32,48,64`) gives
  pressure order 2.34 over all fluid cells and 2.13 on the surface.
- The same 3-D GPU frozen reconstruction agrees with CPU to
  `1.02e-9 Pa` in surface-pressure Linf and passes finite/positive VTP audits.
- A nonzero-gradient no-slip NS MMS on a fixed 4,096-segment circle gives fitted
  L2 orders 2.32 for pressure, 1.87 for `dTdn`, and 1.89 for shear over
  `N=32,48,64,96`. The normal viscous traction gives all-grid L2 order 1.82
  and a last-pair order 2.27; its Cartesian stress identity agrees with an
  independent symbolic evaluation to `3.6e-15`. All surface stencils are
  complete and no pressure fallback is used.
- A fixed 20,480-facet sphere over `N=16,24,32,48` gives pressure order 2.83,
  but `dTdn` order 1.80 and shear L2 order 1.65 because its facet-normal error
  reaches the fine-grid floor. Repeating the entire sequence with one fixed
  81,920-facet STL gives pressure order 2.84, `dTdn` order 1.96, and shear L2
  order 1.92 (L-infinity 1.96). No grid-dependent geometry is used in either
  sequence; the A/B result isolates the input-geometry accuracy requirement.
- For an `M=2`, 20-degree sharp Prandtl-Meyer expansion, the geometric pressure
  fallback reduced the `240x288` near-corner surface-pressure minimum from
  `0.0125` to `0.1497` and the ramp pressure L1 error from `16.05%` to `7.50%`.
  Over `120x144`, `180x216`, and `240x288`, surface L1 converges at about first
  order, while post-fan pressure L1 converges at order 1.4--1.5 and the
  under-pressure area contracts. Bulk angle, Mach, pressure, and density stay
  within 0.1% of the exact Prandtl-Meyer state on the finest grid. Pointwise
  corner `Linf` is not claimed because the exact boundary normal is singular.

A fixed triangulated sphere is not an exact analytic sphere. If its facet
location and normal errors are not reduced with grid spacing, they impose an
error floor on an evolved solution even though frozen reconstruction is
second-order. The solver does not repair STL geometry; production geometry
must be audited to the required location, normal, closure, and scale tolerance.

The generic fixed-STL production and validation contract is documented in
`docs/ibm_fixed_stl_pressure.md`. Increasing `ghost_layers` does not supply
pressure compatibility; that requires a consistent wall closure.
