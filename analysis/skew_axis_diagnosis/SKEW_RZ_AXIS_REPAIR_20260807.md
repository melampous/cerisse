# Skew symmetric R--Z axis diagnosis and repair

## Conclusion

The prominent near-axis line in the Run165 analytic-exit calculation was not
caused by IBM geometry, an AMR coarse--fine interface, or an incorrect
symmetry ghost fill.  Its primary deterministic source was the R--Z JST
dissipation in `src/rhs/Skew.h`: the old radial path applied the dissipative
stencil to the metric-weighted state `r U`.  The R--Z divergence already
multiplies an ordinary face flux by the face radius, so the dissipative state
must be the physical conservative state `U`.

The repaired path:

1. applies radial JST dissipation to `U`, not `r U`, away from the axis;
2. adds no physical JST increment at the zero-area axis face, while retaining
   the high-order central auxiliary metric flux stored there;
3. when the first interior radial face is actually replaced by local LLF,
   replaces its axis companion in the same pass using a parity-consistent
   metric auxiliary flux and `P_axis = p_1`;
4. uses the same paired pressure in the complete flux and the ordinary
   pressure-gradient channel;
5. keeps the optional fixed-width axis collar disabled by default
   (`cns.skew_axis_collar_faces=0`).

The final source revisions assessed here are:

```text
src/rhs/Skew.h          9d50e7ca27e6671fc19a14b6e243658dc64b925c35bb0910738b34f3f48898c3
src/tim/compute_rhs.cpp 23ef010c74d2d541437d888e258ca0975ff511b9268ca30aa82bd220b5d550a2
src/tim/advance.cpp     9f7a288bfe00f6b32fb2d6ba2e1d4338784d4a72ba46beaa74816b106ba314e3
```

## Why the old metric JST creates a line

Let a radial numerical face flux contain JST dissipation

\[
  d_f = \epsilon_f \sum_\ell c_\ell U_{f+\ell}.
\]

For the fourth-order shock term,

\[
  c=\frac{1}{6}(-1,-5,7,-1),\qquad
  \sum c_\ell=0,\qquad \sum c_\ell\xi_\ell=1.
\]

If `r U` is inserted instead of `U`, even a constant physical state produces

\[
  d_f^{old}=\epsilon_f\frac{\Delta r}{r_f}U.
\]

The contribution cancels only in the accidental special case that
`epsilon_f` is identical on adjacent faces.  Near a shock the JST sensor, and
therefore `epsilon_f`, changes sharply; the R--Z divergence then converts the
uncancelled term into a grid-attached pseudo-source.  Its `1/r_f` scaling makes
it most visible in the first few radial rings.

The algebra regression in
`analysis/skew_axis_diagnosis/jst_physical_state_algebra_gate.py` gives, for a
constant `U` and deliberately varying face coefficients,

```text
old metric-state JST RHS: [1000, -166.6666667, 200]
new physical-state JST RHS: roundoff (about 1e-13)
```

## Selected first-face/axis closure

Denote the first two ring centres by

\[
  r_1=\frac{h}{2},\qquad r_2=\frac{3h}{2},
\]

and the first interior face by `r_f=h`.  The metric value corresponding to the
two-state physical LLF face is

\[
  G_1=\frac12(r_1F_1+r_2F_2)
      -\frac12 r_f\alpha(U_2-U_1).
\]

When this face is selected, the non-radial-momentum axis auxiliary value is

\[
  G_0=\frac12\left[(r_1-r_f)F_1+(r_2-r_f)F_2\right].
\]

It satisfies the exact identity

\[
  G_1-G_0=r_f\left[
      \frac12(F_1+F_2)-\frac12\alpha(U_2-U_1)\right],
\]

so the first-ring metric difference is the ordinary physical two-state LLF
operator rather than a mixture of a low-order high face and a high-order axis
face.

For radial momentum, pressure is removed from the metric advective channel,
the axis pressure companion is set to `p_1`, and the odd-parity leading LLF
jump is changed from `m_2-m_1` to `m_2-3m_1`.  This removes the `O(1)` defect
for a regular axis state `m_r=a r+O(r^3)`.  This is a nonlinear rescue for a
selected face; it is not an invariant-domain proof for the full Skew scheme.

## What was ruled out

- **IBM:** the reproducing Run165 case has no immersed geometry.
- **AMR interface:** its active field is a uniformly covered L1 grid; there is
  no coarse--fine interface in the first rings.
- **axis ghost fill:** R--Z MMS with AMReX symmetry fill and exact analytic
  odd/even extension produced byte-identical plotfiles.
- **lack of background dissipation alone:** with `C2=C4=0`, local LLF prevents
  the immediate positivity abort but leaves a visibly dirtier field.  JST is
  still needed; its radial metric implementation had to be corrected.
- **a wider fixed low-order collar:** K=1, 2, and 4 moved the strongest stripe
  to the first high-order face instead of removing it.  The production default
  is therefore K=0.

## Verification completed

### Smooth accuracy

Skew4-JST Euler MMS at N=16, 32, and 64 retained approximately third-order
global convergence, as expected for the present first-power JST sensor.  The
first-ring radial-momentum error converged at rates 3.83 and 3.19.  Fallback
on/off plotfiles were byte-for-byte identical in the smooth runs.  Cartesian
MMS plotfiles were also byte-for-byte unchanged by the R--Z-only repair.

Details: `analysis/skew_axis_diagnosis/jst_unweighted_validation/README.md`.

### Run165 nonlinear gate

The analytic-exit, jet-on, no-IBM Run165 case used a uniform L1 field with
64 x 192 finest cells, `dr=dz=1.0612967 mm`, fixed coarse step `0.4 us`,
Skew4-JST `C2=1.5`, `C4=0.016`, local LLF enabled, and fixed collar K=0.

| physical time | coarse steps | wall time | wall/coarse step | result |
|---:|---:|---:|---:|---|
| 500 us | 1250 | 24.47 s | 0.01957 s | pass |
| 1.5 ms | 3750 | 72.80 s | 0.01941 s | pass |

There was no negative-pressure, negative-internal-energy, non-finite-state, or
positivity abort.  The residual near-axis extrema did not grow or migrate
between 500 us and 1.5 ms.

A synchronous 200 us L1/L2 comparison then halved `dr` from 1.0613 mm to
0.53065 mm.  The dominant pressure local-extremum p95 residual fell from
7.28% to 2.61%, density from 4.90% to 2.49%, temperature from 5.80% to 3.17%,
radial velocity from 9.57% to 2.96%, and Mach from 2.67% to 1.47%.  The
pressure residual's radial width above 2% shrank from 5.31 mm to 3.18 mm, and
the Mach width shrank from 1.06 mm to zero.  At the physical radius of the L1
main line, the five residuals decreased by factors 5.19, 4.03, 4.09, 3.42,
and 3.95, respectively.  The pressure even-parity p95 defect fell from 11.84%
to 2.59%; the `u_r/c` odd-parity p95 defect fell from 2.42% to 0.984%.

An L1 control run reduced the finest time step by a factor four while keeping
the same spatial grid.  Its five main residuals changed by no more than 1.57%,
which rules out the smaller L2 time step as the explanation.  The remaining
line is therefore a spatially converging coarse-grid shock/gradient sampling
error, not the old nonconvergent metric-JST pseudo-source.  This evidence also
rejects adding a permanent low-order axis collar.

Refinement evidence:
`analysis/skew_axis_diagnosis/refinement_audit/README.md`.

Main nonlinear report:
`analysis/skew_axis_diagnosis/jst_physical_jump_final/README.md`.

### Pure shared-GP integration

A two-dimensional R--Z pure shared-GP CPU build and one-step uniform-flow
smoke passed with the IBM positivity limiter disabled.  The manifest reported

```text
family=pure_shared_gp_full_cartesian
cut_control_compiled=0
control_volume=full_cartesian_cell
face_area=full_cartesian_face
solid_extension=shared_ghost_point
```

All 42 ghost points used quadratic BI-CWLS.  The uniform-state variation after
one step was at roundoff.  The repair reads only the existing full-cell/shared
GP states and marker usability; it introduces no cut volume, area fraction,
wall subface, aggregate, redistribution, or EB flow update.

## Current qualification boundary

- The corrected all-fluid R--Z Skew-JST path has passed smooth accuracy and a
  1.5 ms strong-jet stability gate.
- The pure shared-GP result is only a build/interface/uniform-transparency
  smoke.  It is not the mandatory Mach-4 circle, PM expansion,
  compression-corner, curved-wall MMS, and inverse-MOC production ladder.
- The production conservative shared-face positivity limiter is currently
  certified only for source-free LLF-WENO-Z5; enabling it with Skew still
  fails closed.  Skew therefore has no strict stage-wise positivity or
  invariant-domain guarantee yet.
- Dynamic R--Z AMR reflux with the paired pressure channel is not certified by
  the uniform-L1 Run165 calculation.
- The same source completed a CUDA compile and link, but GPU runtime was not
  exercised because the current OS environment denied GPU/NVML access.
- With `rz_gp_annular_bic=1`, all-fluid JST currently dissipates the stored
  annular conservative state while the reduced shared-GP JST path reconstructs
  a point conservative state from primitives.  Those semantics coincide in
  ordinary point-state mode but not in annular-BIC mode.  Quantitative annular
  pure-GP Skew use is therefore blocked until that variable choice is unified
  and verified by the mandatory GP regression ladder.

Consequently this repair makes Skew a substantially cleaner and stable R--Z
alternative for all-fluid tests.  It must not yet replace the frozen
LLF-WENO-Z5 pure-GP production method without the remaining production gates.
