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
    euler_slip_analytic_curvature
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
outside the callback's declared smooth region. Invalid geometry or state data
falls back to fluid-side pressure extrapolation.

Selection is compile-time:

```cpp
struct ibmparm_t {
    static constexpr int interp_order = 2;
    static constexpr int extrap_order = 2;
    static constexpr int interp_order_surf = 2;
    static constexpr int extrap_order_surf = 2;
    static constexpr int ghost_layers = 1;
    static constexpr ibm_pressure_closure_t pressure_closure =
        ibm_pressure_closure_t::euler_slip_analytic_curvature;
    // wall_shape_operator(...) must also be defined.
};
```

The curvature policy is rejected at compile time for built-in no-slip walls.
It does not infer curvature from STL facets and is not combined with the shock
fallback. Use `shock_aware_fluid_extrapolation` for shock-dominated generic
STL calculations.

## Verification Snapshot

The analytic-curvature MMS uses a radius-0.2 circle/sphere and an isothermal
circumferential Euler field satisfying the compatibility condition exactly.

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

A fixed triangulated sphere is not an exact analytic sphere. If its facet
location and normal errors are not reduced with grid spacing, they impose an
error floor on an evolved solution even though frozen reconstruction is
second-order. The solver does not repair STL geometry; production geometry
must be audited to the required location, normal, closure, and scale tolerance.

The next pressure-accuracy item is the no-slip Navier--Stokes compatibility
term involving `n dot div(tau)`, followed by surface-force/control-volume
momentum-budget closure. Increasing `ghost_layers` does not supply either one.
