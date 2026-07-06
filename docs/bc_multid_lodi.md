# Multidimensional LODI / NSCBC Boundary Architecture

## Problem

The current `bc_nscbc_farfield` helper is a local, algebraic ghost-cell
closure.  It uses only the nearest interior state and the outward normal.  This
is enough for a normally outgoing acoustic wave, but it is not enough for a
multidimensional oblique wave.

For a right boundary and a small acoustic perturbation,

```text
delta u_n = delta p / (rho c)                 normal wave
delta u_n = cos(theta) delta p / (rho c)      oblique wave
delta u_t = sin(theta) delta p / (rho c)
```

The one-dimensional incoming invariant is not zero for an oblique outgoing
wave:

```text
delta J_- = delta u_n - delta p/(rho c)
          = (cos(theta) - 1) delta p/(rho c).
```

Thus a boundary that blindly sets `J_-` to the far-field value will reflect or
distort oblique acoustic content.

## LODI Model

At a locally planar boundary, write the inviscid Euler operator in normal and
tangential coordinates,

```text
partial_t Q + A_n partial_n Q + A_t partial_tan Q = 0.
```

For a 2D right boundary linearized about a uniform state, the acoustic normal
characteristic variables are

```text
W+ = p' + rho c u_n'
W- = p' - rho c u_n'.
```

Their local one-dimensional inviscid equations become

```text
partial_t W+ + (u_n + c) partial_n W+ = -rho c^2 partial_tan u_t'
partial_t W- + (u_n - c) partial_n W- = -rho c^2 partial_tan u_t'.
```

The right-hand side is the transverse LODI term.  Dropping it recovers a
one-dimensional characteristic boundary.  For oblique waves that term is the
piece that balances the apparently incoming normal invariant `W-`.

For nonlinear 2D/3D use, the boundary should work with local primitive
variables and characteristic amplitudes,

```text
L+ = (u_n + c) (partial_n p + rho c partial_n u_n)
L- = (u_n - c) (partial_n p - rho c partial_n u_n)
Ls = u_n (partial_n rho - partial_n p / c^2)
Lt = u_n partial_n u_t
```

plus the transverse source terms from each tangential direction in
`A_t partial_t Q`.  Incoming amplitudes are then set by a policy, for example
pure radiation, weak pressure relaxation, or sponge/buffer damping.

For the 2D right-boundary validation case with primitive variables
`q=(rho,u,v,p)` and mean flow mostly in `+x`, the primitive equations imply

```text
p_t   + (L+ + L-)/2
      + v p_y + rho c^2 v_y = 0

u_t   + (L+ - L-)/(2 rho c)
      + v u_y = 0

v_t   + Lt
      + v v_y + p_y/rho = 0

rho_t + Ls + (L+ + L-)/(2 c^2)
      + v rho_y + rho v_y = 0
```

For this 2D notation, the incoming acoustic compatibility source is therefore

```text
T- = v p_y + rho c^2 v_y - rho c v u_y .
```

The generic 3D helper applies the same form summed over tangential directions:

```text
T- = sum_t [ u_t p_t + rho c^2 partial_t u_t - rho c u_t partial_t u_n ] .
```

In a pointwise boundary-equation implementation, setting `L-=-T-` corresponds
to freezing the incoming characteristic variable.  In a ghost-cell
reconstruction this is too strong for oblique outgoing sound waves.  For a
locally outgoing acoustic branch with normal direction cosine `n_x`, the exact
normal derivative contribution is smaller by

```text
beta = (1 - M_n) n_x / (1 + n_x)
```

for the uniform-flow validation setup.  The current native implementation uses
that idea locally: it estimates `n_x` from the tangential pressure/velocity
gradient relation, gates the correction off when the signal looks vortical
rather than acoustic, and uses

```text
L- = sigma (p - p_inf) - beta T- .
```

This is a characteristic-amplitude ghost-fill LODI approximation.  It is the
right theoretical level for a conservative finite-volume code that fills ghost
cells before the Riemann solve, but it is still not a full Poinsot-Lele
boundary-equation update in which the boundary-cell RHS is replaced by the LODI
relations.

## Required Information

A true multidimensional LODI / NSCBC implementation needs more than
`bcnormal()` currently receives:

1. The boundary face normal and side.
2. The first interior state.
3. Tangential neighbor states on the boundary plane.
4. Tangential derivatives of primitive variables.
5. A policy for incoming characteristic amplitudes:
   - pure non-reflecting,
   - weak relaxation to far-field,
   - sponge-layer damping,
   - user-prescribed inflow/outflow data.

## Code Architecture

The existing pointwise API is kept:

```cpp
bcnormal(x, dratio, s_int, s_refl, s_ext, idir, sgn, time, geom, cls, pparm)
```

A new optional API is introduced for cases that need multidimensional boundary
data:

```cpp
bool bcnormal_lodi(iv, state, x, dratio,
                   s_int, s_refl, s_ext,
                   idir, sgn, time, geom, cls, pparm)
```

`bcnormal_lodi` returns `true` if it handled the ghost cell.  If it is absent or
returns `false`, the code falls back to `bcnormal`.  This keeps all existing
cases source-compatible.

This hook is now wired in `src/set/bcs.cpp`.  It is deliberately optional:
existing `prob.h` files do not need to change, and cases without
`bcnormal_lodi` take the compile-time fallback with no LODI work.

The reusable production-oriented helper lives in `src/set/bc_types.h`:

```cpp
GlobalBC::nscbc_lodi_bc_parm_t bp;
bp.rho_inf = rho_inf;
bp.u_inf = u_inf;
bp.v_inf = v_inf;
bp.w_inf = w_inf;
bp.p_inf = p_inf;
bp.pressure_relax = sigma_l_ref_over_c;
bp.transverse_relax = 1.0;  // projected transverse term, default
bp.l_ref = l_ref;

return GlobalBC::bc_nscbc_lodi_farfield(
    iv, state, s_ext, idir, sgn, geomdata, &closures, bp, Yinf);
```

The helper is dimension-generic for Cartesian 2D/3D boundaries.  It reads the
first interior cell, the inward neighbor, and available tangential neighbors;
computes `L+`, `L-`, entropy, and tangential outgoing gradients; applies the
summed transverse term to the incoming acoustic amplitude; and reconstructs a
conservative ghost state.  Supersonic cases and inflow cases intentionally
delegate to the fixed characteristic far-field helper.

## Implemented Native Modes

The native validation case now has two opt-in multidimensional modes:

```text
prob.bc_mode = 2   legacy plane-wave acoustic extrapolation
prob.bc_mode = 3   generic characteristic-amplitude LODI ghost-fill helper
```

The legacy `bc_mode=2` model is deliberately limited to a matched acoustic
test:

- right x-boundary only,
- calorically perfect single-species gas,
- small-amplitude acoustic packet,
- y-periodic tangential direction,
- local tangential derivative from first interior cells.

It estimates the outgoing acoustic direction from the local pressure and
tangential velocity perturbations, then extrapolates the acoustic primitive
state along the inferred wave-normal direction:

```text
sin(theta) ~= rho_inf c_inf delta u_t / delta p
partial_x q ~= cot(theta) partial_y q
q_g = q_i + distance_to_ghost partial_x q
```

Degenerate points fall back to the fixed characteristic far-field.  This is not
production NSCBC.  It remains only as a comparison curve.

Because it infers a single local acoustic angle from `delta u_t / delta p`, the
prototype is expected to work best for clean single-angle acoustic content.  It
is not designed to decompose mixed acoustic/vortical/entropy/shock content.

The new `bc_mode=3` calls `bc_nscbc_lodi_farfield` from `src/set/bc_types.h`.
It computes local amplitudes `L+`, `L-`, entropy, and tangential modes from
normal primitive gradients.  It extrapolates outgoing acoustic, entropy, and
tangential information from the interior, prescribes the incoming acoustic
amplitude with pressure relaxation plus a projected transverse term, and
reconstructs ghost primitive gradients.

The production helper also contains local mixed-flow limiters:

- strong pressure/compression sensor: a shock or strong compression crossing
  the boundary is copied out by low-constraint extrapolation;
- pressure-poor velocity packet sensor: wake/vorticity packets with large
  velocity content and weak acoustic coherence are copied out instead of being
  forced through the acoustic incoming amplitude;
- optional convective extrapolation: a limited primitive-gradient extrapolation
  for pressure-poor velocity packets is available for experiments, but it is
  disabled by default because the native wake packet test was worse than the
  conservative copy-out fallback;
- entropy/contact mode: pure density/entropy content is not copied out by the
  wake sensor; it remains an outgoing entropy characteristic;
- pressure-relaxation weighting: low-coherence signals receive only weak
  far-field pressure anchoring, avoiding unnecessary conversion into sound.

The exposed native parameters are:

```text
prob.lodi_relax                         pressure-relaxation strength, default 0
prob.lodi_transverse                    projected transverse multiplier, default 1
prob.lodi_shock_pressure_jump           shock pressure-jump sensor, default 1e-2
prob.lodi_shock_compression             compression sensor, default 1e-3
prob.lodi_acoustic_coherence            acoustic coherence threshold, default 0.35
prob.lodi_convective_velocity_jump      velocity-packet floor, default 1e-5
prob.lodi_convective_pressure_fraction  pressure-poor packet threshold, default 0.25
prob.lodi_pressure_relax_min_weight     minimum pressure-relax weight, default 0.10
prob.lodi_use_shock_sensor              enable shock fallback, default 1
prob.lodi_use_convective_sensor         enable wake/vorticity fallback, default 1
prob.lodi_use_convective_extrapolation  experimental limited extrapolation, default 0
prob.lodi_use_pressure_relax_sensor     enable relaxation weighting, default 1
```

The native validation case also exposes an RHS sponge/buffer source through
`prob.h`.  It is off by default and is activated by:

```text
prob.sponge_strength  damping rate, default 0
prob.sponge_width     buffer thickness, default 0.20
prob.sponge_power     polynomial ramp exponent, default 2
prob.sponge_x_hi      enable high-x buffer, default 1
prob.sponge_x_lo      enable low-x buffer, default 0
prob.sponge_y_lo/hi   enable transverse buffers, default 0
```

The source term is conservative-state relaxation toward the far-field state,

```text
dU/dt += -sigma(x) (U - U_inf).
```

For production cylinder/SRP cases this should be ported into the relevant
`prob.h` or a reusable source functor, then tuned by side-domain sensitivity.

## Native Validation Snapshot

All results below use the compiled Cerisse path, comparing a short domain
against a long-domain reference at the same grid spacing.

```text
clean acoustic angle sweep, quick N=48/96:
case        BC        N=48       N=96
15 deg      lodi      1.835e-04  4.703e-04
15 deg      char      9.761e-04  1.825e-03
15 deg      foextrap  6.931e-04  1.631e-03
30 deg      lodi      1.553e-03  3.399e-03
30 deg      char      2.232e-03  2.027e-02
30 deg      foextrap  3.139e-03  1.828e-02
45 deg      lodi      9.884e-03  1.430e-02
45 deg      char      3.780e-02  1.014e-01
45 deg      foextrap  3.397e-02  6.015e-02
60 deg      lodi      2.959e-02  4.190e-02
60 deg      char      1.047e-01  1.913e-01
60 deg      foextrap  7.575e-02  9.018e-02
```

The over-strong unprojected `L-=-T-` variant was also tested.  It was worse than
the projected form, especially at 30 and 45 degrees.  That is the concrete bug
fixed by the current formula.

```text
component / robustness quick N=48/96:
case          lodi N=96   char N=96   foextrap N=96   best note
acoustic30    4.623e-03  2.027e-02  1.828e-02      lodi best
entropy       1.308e-02  6.763e-01  5.943e-02      lodi best; equals raw LODI
vortex        3.305e-01  1.637e-01  2.677e-01      not solved by LODI
wake          1.027e-01  2.594e-01  9.316e-02      lodi improves raw slightly; foextrap lower
weak shock    3.341e-03  1.394e-02  1.392e-02      lodi best in weak-compression test
```

After adding the sponge source and testing the wake/pressure-pulse physical
region (`x <= 0.75`, sponge begins at `x=0.80`, quick N=48/96), the refined
result was:

```text
case            BC                N=96 metric   note
wake            lodi              4.990e-02    conservative copy fallback
wake            lodi_convective   6.560e-02    worse; remains experimental
wake            lodi_sponge       3.809e-02    better than lodi
wake            foextrap          3.915e-02    close to lodi_sponge
wake            foextrap_sponge   3.288e-02    best in this simple wake packet
pressure pulse  lodi              5.924e-02    baseline LODI
pressure pulse  lodi_sponge       4.176e-02    better
pressure pulse  foextrap          6.153e-02    baseline extrapolation
pressure pulse  foextrap_sponge   4.116e-02    best in this simple pressure pulse
```

This confirms two practical points.  First, a naive linear convective ghost
extrapolation is not an improvement for the current wake packet, so it is not a
production default.  Second, a weak buffer layer is the more promising path for
wake and low-frequency pressure feedback, but it must be placed outside the
measurement region and tuned with the actual flow case.

After moving `bc_mode=3` to the generic helper, a smoke run on the native
`acoustic30` case gave:

```text
case        BC        N=48       N=96
30 deg      lodi      2.897e-03  4.623e-03
30 deg      char      2.232e-03  2.027e-02
```

The generic helper remains much better than the fixed characteristic far-field
for the refined oblique acoustic test, but it is slightly less tuned than the
previous case-local projected prototype.  That is expected: the new code avoids
using the known validation wave angle and is meant to be reusable in 3D.

RHS-level NSCBC hook snapshot:

```text
implementation:
  compute_rhs calls optional rhs_nscbc() after FV flux divergence and EB
  redistribution, before problem source terms.
  bc_mode=4 in exm/numerics/bc_native applies a boundary-cell characteristic
  amplitude correction, not a full overwrite of all primitive modes.
  The existing FV RHS is converted to primitive time derivatives, the current
  incoming acoustic amplitude L- is inferred with full transverse terms, and
  only L- is corrected toward the pressure-relax/projected-transverse target.
  Outgoing acoustic, entropy, tangential, shock/wake-sensor fallback content is
  left as close as possible to the compiled FV/WENO path.
```

The all-mode RHS overwrite prototype was rejected: on `acoustic30` it gave
`1.432e-02` at N=96, and on `entropy` it gave `1.961e-02`, both worse than the
ghost LODI baseline.  The accepted characteristic-correction form gives:

```text
case        BC          N=48       N=96       note
acoustic30  rhs_lodi    2.731e-03  5.439e-03 close to ghost LODI
acoustic30  lodi        2.897e-03  4.623e-03 ghost-fill baseline
entropy     rhs_lodi    3.927e-02  1.308e-02 identical to ghost LODI
entropy     lodi        3.927e-02  1.308e-02 ghost-fill baseline
```

The projected transverse target is required.  Setting the RHS target to the
full unprojected transverse correction over-forces the oblique acoustic case
(`acoustic30`, N=96: `2.653e-02`).  Setting no transverse target is stable but
less accurate (`1.347e-02` at N=96).  The projected form is therefore the
current RHS-level default.

Mixed-flow quick results using the whole short-domain metric (`x in [0,1]`):

```text
case            BC                N=96 metric   note
wake            rhs_lodi          1.194e-01    decreasing but not converged
wake            rhs_lodi_sponge   1.328e-01    sponge setting not better here
wake            lodi_sponge       1.295e-01    similar residual floor
wake            foextrap_sponge   1.311e-01    similar residual floor
pressure pulse  rhs_lodi          5.427e-02    best in this whole-domain quick metric
pressure pulse  rhs_lodi_sponge   7.723e-02    current sponge too strong
pressure pulse  lodi_sponge       6.755e-02    current sponge too strong
pressure pulse  foextrap_sponge   6.615e-02    current sponge too strong
weak shock      rhs_lodi          3.337e-03    decreasing, best of this set
weak shock      rhs_lodi_sponge   2.834e-02    sponge smears/reflection metric
weak shock      lodi_sponge       2.834e-02    same as RHS+sponge
weak shock      foextrap_sponge   3.569e-02    larger residual
```

These results do not overturn the production strategy: RHS-level acoustic
control is now usable for simple acoustic/entropy/weak-compression tests, but
wake-dominated lateral boundaries still require side-domain and buffer-layer
sensitivity.  Sponge strength and placement must be tuned against the actual
cylinder/SRP flow; the quick pressure-pulse result shows that a poorly placed
or too-strong sponge can be worse than no sponge.

The wake and vortex results are the main caution for cylinder lateral
boundaries: if the side boundary is dominated by slow vortical/wake content,
foextrap, a buffer layer, or a larger side domain can still be lower risk even
though the acoustic LODI is better for clean outgoing sound.  The LODI helper
now avoids obvious acoustic over-forcing of wake packets, but it does not
magically convect all rotational content out of a short lateral domain.

Generated reports:

```text
temp/bc_lodi_amp_angle_projected_quick/native_bc_validation_report.pdf
temp/bc_lodi_amp_components_projected_quick/native_bc_validation_report.pdf
temp/bc_lodi_amp_shock_projected_quick/native_bc_validation_report.pdf
temp/bc_lodi_physics_completeness_quick4/native_bc_validation_report.pdf
temp/bc_lodi_sponge_convective_quick2/native_bc_validation_report.pdf
temp/bc_lodi_regression_quick/native_bc_validation_report.pdf
temp/bc_rhs_nscbc_projected_target_quick/native_bc_validation_report.pdf
temp/bc_rhs_nscbc_mixed_quick/native_bc_validation_report.pdf
temp/cylinder_side_domain_sensitivity_manifest/cylinder_side_domain_sensitivity.md
```

## Production Path

For production use, this should remain opt-in until restart-based cylinder/SRP
side-domain sensitivity checks have been run.  The current code is a real
2D/3D characteristic-amplitude ghost-fill implementation in the reusable BC
layer, but it is not a universal non-reflecting boundary and it is not a
replacement for a sponge/buffer layer in wake-dominated side boundaries.

If the cylinder lateral boundaries remain dominated by shocks, vortical wake,
and entropy content, the current evidence points to LODI or foextrap plus a
weak sponge/buffer and side-domain sensitivity, not LODI alone.  A future
production BC code should be a separate opt-in mode after 2D/3D acoustic and
flow-case validation, not a silent replacement for current code 7.

Current implementation status:

- `code 7` in input files still means fixed characteristic far-field through
  the existing `bcnormal` path.
- Multidimensional LODI is activated only by a case-level `bcnormal_lodi`
  implementation that calls `bc_nscbc_lodi_farfield`.
- RHS-level NSCBC is activated only by a case-level `rhs_nscbc` implementation
  that calls `rhs_nscbc_lodi_farfield`; it is currently an inviscid open-face
  characteristic correction, not a full viscous/thermal/species NSCBC closure.
- 2D native validation and 3D native compilation have passed.
- Viscous far-field terms are not explicitly characteristic-decomposed.  This is
  consistent with the common NSCBC practice of treating open boundaries through
  inviscid characteristics while keeping the boundary far enough that viscous
  terms are small, but near-wake production runs still need domain/sponge
  validation.
