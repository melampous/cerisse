# Conservative VTP load and production shifted shell

## Scope

This change closes two production ambiguities in the fixed-body Cartesian
GPIBM path:

1. VTP wall fields and the conservative Cartesian interface load use different
   discrete operators. An old mode-6 checkpoint continued with mode 8 showed
   an apparent 15--18 percent difference because the operator-change transient,
   instantaneous WLS reconstruction, and an RK-averaged load window were mixed.
2. High-order same-fluid near-wall faces had only been exercised with the
   validation-only `ib.noslip_inviscid_shifted_shell=1` path.

It does not force the physical-wall reconstruction to equal the Cartesian
finite-volume momentum exchange. Those remain independent verification
quantities.

## Conservative load mapped to VTP

For every unique active fluid/solid Cartesian face, the body reaction is

\[
  \boldsymbol F_{b,f}
  =-s_f\,\boldsymbol F_{\rho\boldsymbol u,f}\,A_f,
\]

where `s_f` is the fluid-side orientation and the face flux is the actual
inviscid-plus-diffusive flux used by the Runge--Kutta update. The solver
accumulates the RK-stage impulse, excludes AMR coarse faces covered by a fine
level, maps each face to its cached first-hit surface element, and divides by
the exact surface-output time window.

The VTP fields are:

- `ConservativeElementForce`: element contribution to the integrated body
  force;
- `ConservativeTraction`: that force divided by physical surface-element area;
- `ConservativeCrossingFaces`: number/weight of mapped Cartesian crossings;
- `ConservativeLoadValid`: validity flag;
- `ConservativeLoadWindowStart/End/Dt`: the exact averaging interval.

The total reaction is now separated into
`ConservativeInviscidElementForce` and
`ConservativeDiffusiveElementForce`. The output additionally contains
`LoadConsistentPressure`, `LoadConsistentTau*`, and
`LoadConsistentViscousTraction`. They are minimum surface-L2 corrections of the
raw WLS fields constrained to reproduce the separated conservative forces. A
primal/dual active set enforces positive projected pressure. Raw `Pressure`,
`Tau*`, `Cf`, and `dTdn` remain unchanged.

`Pressure`, `Tau*`, `Cf`, and `dTdn` retain their physical-wall WLS
reconstruction meaning. Their surface integral is useful for local wall
physics, but `production_force_*` and the sum of `ConservativeElementForce`
are the canonical conservative total load.

### AMR/MPI regression

A 40-step, two-rank continuation of the Mach-4, Re=10000 sphere from the old
mode-6 checkpoint produced:

| Quantity | x | y | z |
| --- | ---: | ---: | ---: |
| VTP conservative force | 32554.266744064607 | 457.871524840837 | 189.325068244176 |
| time-window production force | 32554.266744064527 | 457.871524840837 | 189.325068244175 |
| reconstructed pressure force | 37498.827881598394 | 65.375626156486 | -10.711698718969 |

The first two rows agree to roundoff (streamwise relative difference about
`2.5e-15`). The reconstructed pressure difference is deliberately visible.
All 5120 surface elements passed the VTP audit, and the final discrete momentum
budget residual was approximately `6.7e-8` relative.

A subsequent two-rank, two-level direct-mode continuation over ten coarse
steps isolated the current extraction difference on one identical output
window:

| Quantity | x | y | z |
| --- | ---: | ---: | ---: |
| conservative total force | 31261.483496908797 | 312.772476525376 | 349.108777398876 |
| conservative inviscid force | 31132.905154742886 | 312.361876884002 | 348.874742072968 |
| conservative diffusive force | 128.578342165943 | 0.410599641368 | 0.234035325912 |
| raw WLS pressure force | 32766.849699681108 | 353.419699529120 | 270.079242887962 |

The raw streamwise pressure-force difference from the separated conservative
inviscid target is 5.25 percent, not 18 percent. (A separate instantaneous
step-1500 total-load comparison gave about 3.2 percent and must not be mixed
with this RK-window decomposition.) The projected pressure and viscous forces
recover their respective conservative targets, and their sum closes with a
relative residual of `3.8e-15`. The pressure correction is 3.08 percent in the
absolute-pressure surface L2 norm, with no pressure-floor activation. This
does not establish pointwise Cp accuracy; it makes the discrepancy explicit,
separated, and usable without corrupting the raw wall reconstruction.

## Production same-fluid shifted shell

The production path first forms the established reconstructed-GP/masked
LLF-WENO-Z5 or LLF-TENO5 face flux. It then searches for a complete five-cell
block entirely on the connected fluid side. A shifted high-order candidate is
accepted only if both directional translations have positive optimal weights.
The two conservative fluxes are combined with

\[
  \boldsymbol F_f=(1-\theta_f)\boldsymbol F_{\mathrm{GP/masked}}
                  +\theta_f\boldsymbol F_{\mathrm{shifted}},
  \qquad 0\leq\theta_f\leq1,
\]

where `theta` is a cubic-smoothstep function of normalized fluid-side density
and pressure jumps and second differences. Strong compression drives
`theta` to zero. If a stable complete block is unavailable, the GP/masked
flux is retained exactly.

Reduced-candidate extrapolation and extreme signed-optimal-weight translations
were rejected: both generated a growing first-shell mode and eventually
negative density or internal energy in fixed-time full-NS tests. Therefore the
current implementation improves the second same-fluid face when the geometry
permits, while keeping the first face robust. It does not claim a uniformly
fifth-order boundary closure.

## Verification

- Static polynomial, marker-safety, and blend-transition tests pass.
- LLF-WENO-Z5 `t=0` curved-wall MMS: first-face flux L1/L2 orders are about
  2.7--3.0; the second same-fluid face is mostly 4.6--5.6 order.
- LLF-WENO-Z5 `t=0` first-cell momentum-residual orders are 1.85--1.94. This is
  the remaining local semi-discrete limitation.
- Fixed-time full-NS `N=32/48/64/96`: all-fluid L2 orders are 3.02 for density,
  2.72/2.70 for velocity, 2.57 for pressure, and 3.67 for temperature; every
  run reached `t=0.002` without a positivity failure.
- LLF-TENO5 builds and reaches the fixed terminal time in the smoke regression.
- Mach-4, Re=10000 sphere AMR continuation passes the production strict-
  positivity gate, reflux, VTP validity, and conservative-load closure.

## Production interpretation

Use `ib.noslip_inviscid_shifted_shell=1` for the qualified no-slip
LLF-WENO-Z5/TENO5 IBM configurations. The global default remains zero to avoid
silently changing legacy cases. A production case should record the switch,
crossing closure mode, stage-positivity gate used during qualification, and the
conservative load-window metadata.

The implementation establishes a conservative total-force output and a stable
adaptive near-wall convective shell. The remaining accuracy target is the
first-cell semi-discrete momentum operator, not the already closed VTP force
bookkeeping.
