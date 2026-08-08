# IBM AMR metadata, impermeability, and load closure

Date: 2026-07-16

## Scope

This change addresses three fixed-body Cartesian GPIBM invariants for the
single-species ideal-gas stationary no-slip path:

1. one geometry/reconstruction record per Cartesian crossing face and AMR
   level, independent of FAB overlap and MPI decomposition;
2. a shock/shear-aware no-penetration mass and inviscid-energy closure that
   retains the smooth-wall spatial-accuracy contract; and
3. one canonical integrated production load taken from the same conservative
   crossing flux that updates the fluid.

Moving walls, RZ coordinates, multiple species, and FSI force distribution are
outside this result.

## AMR metadata ownership

Boundary-face construction now emits a record only from the FAB whose valid
cell box owns the adjacent fluid cell. A solid-side ghost copy at a FAB
boundary can no longer construct a second record with a different local WLS
support block. The existing composite AMR uncovered-cell mask remains a
separate coarse/fine accounting rule.

On the Mach-4 sphere L1 hierarchy, both one and four MPI ranks report:

- 1248 active crossing records and 1248 accounting owners;
- zero duplicate copies and zero duplicate accounting-owner groups;
- zero first-hit or reconstruction inconsistencies; and
- zero coarse-fine missing GP or surface supports.

## Mode-8 no-penetration closure

In smooth flow, the mass and inviscid-energy flux scalars are reconstructed
with the cached hard-zero wall-constrained cubic functional. Pressure and
momentum retain the translated characteristic LLF-WENO/TENO construction.
The smooth path therefore remains an analytic-extension ghost-cell method; it
does not replace a curved wall by a zero-flux staircase.

Two fluid-side sensors identify a nonsmooth or unresolved wall band:

- normalized support pressure range, with default limits 0.03 and 0.10; and
- support velocity span divided by minimum sound speed, with default limits
  0.20 and 0.50.

A cubic smoothstep blends momentum toward the robust face-state/one-sided LLF
closure. The same activation blends stationary-wall advective mass and energy
to zero. If the cubic metadata path is unavailable, the already assembled
legacy momentum flux is retained while mass and inviscid energy are set to
zero. Diffusive heat and viscous work remain owned by the direct diffusive
closure.

The rejected alternative was to force zero flux on every smooth Cartesian
crossing. It gave 100 percent relative mass and energy flux defects and only
about first-order convergence because a Cartesian crossing face is not the
physical curved wall.

## Accuracy gates

For the 2-D smooth circular full-NS MMS on N=64/96/128/192, the minimum fitted
orders over L1/L2/Linf at the crossing face are:

| component | minimum fitted order |
| --- | ---: |
| mass | 5.093 |
| x momentum | 4.215 |
| y momentum | 3.749 |
| total energy | 3.994 |

The earlier validation-only, unconditional shifted-shell experiment gave
high initial face-flux orders but was not a stable production boundary
operator. The production hybrid now rejects extreme signed-weight
translations and reduced-candidate extrapolation. It accepts a complete
shifted fifth-order block only when both translations have positive optimal
weights and the connected fluid-side density/pressure sensor is smooth. If
that block is unavailable, as is normally the case at the first same-fluid
face, the established reconstructed-GP/masked flux is retained. This is an
adaptive conservative shell closure, not a claim that every first-shell face
uses a fifth-order stencil.

For the resulting fixed-time 2-D full-NS `N=32/48/64/96` sequence, all-fluid
L2 fitted orders are 3.02 (density), 2.72/2.70 (velocity), 2.57 (pressure),
and 3.67 (temperature). The first-fluid-cell fits are 3.10, 2.65/2.69, 2.80,
and 3.72, respectively. These exceed the second-order production requirement
for this smooth MMS; they are verification margins and not a claim that an
immersed boundary is globally third or fourth order. At `t=0`, the fitted
first-layer momentum-residual orders remain about 1.85--1.94, so the local
semi-discrete pressure/convective shell is still the limiting operator.

The same-fluid-face and first-layer figures in the preceding paragraph use
`ib.noslip_inviscid_shifted_shell=1`. The hybrid same-fluid-side shell is now
the recommended production path for LLF-WENO-Z5/TENO5; value 0 retains the established
legacy GP/masked-WENO shell for regression. The mode-8 crossing-face orders do not
depend on this switch.

## Impermeability and load evidence

Restarting the D/h=8 Mach-4 sphere from the same previous mode-6 checkpoint
gives an initial normalized net interface mass flux of about 5.8e-6 with mode
8, compared with about 0.95 for mode 6. Of 312 crossing faces in this test,
311 are fully activated by the nonsmooth-band closure and one is blended; 274
use the cubic WLS path and 38 use the conservative fallback state. This is a
local regression and not a new statistically stationary sphere validation; a
fresh mode-8 production run is still required before quoting drag uncertainty.

The momentum-budget CSV now writes `production_force_*` from the composite
stage-weighted Cartesian crossing flux and records source code 1. The analyzer
rejects another source code or any difference between production and interface
forces. Reconstructed STL surface force remains separate so a traction error
cannot be hidden by relabelling it as conservative load.

For the same one-step restart, the all-conserved discrete budget has relative
residuals of 2.9e-10 for mass, 1.2e-9 for momentum, and 2.8e-9 for total
energy. The reconstructed surface force and conservative interface force
differ by about 18 percent on this non-equilibrated transition state. This
difference is exposed as a diagnostic and is not used to alter the production
load. A later 40-step, two-rank AMR continuation maps the actual stage-weighted
crossing reaction into each VTP: summing `ConservativeElementForce` gives
`(32554.266744064607, 457.871524840837, 189.325068244176)`, while the matching
budget-window `production_force` is
`(32554.266744064527, 457.871524840837, 189.325068244175)`. Their streamwise
relative difference is about `2.5e-15`. The independent reconstructed pressure
force is `(37498.827881598394, 65.375626156486, -10.711698718969)`, and remains
a wall-traction diagnostic rather than the conservative total force.

The complete 3-D build and one-step restart also pass with CUDA `sm_89`. The
GPU run has 312 unique crossing records, no duplicate or inconsistent
metadata, no positivity failure, and all production-load identity checks pass.
Its relative all-conserved residuals are 5.2e-12 for mass, 1.2e-11 for
momentum, and 6.2e-12 for total energy. Against the CPU run, every per-face
path decision is identical and the maximum relative difference among the five
conservative crossing-flux components is 6.8e-11. On the smooth N=64 circle
MMS with the same validation shell configuration, the CPU/GPU actual-flux
difference is at most 2.5e-14 relative to the component scale.

A separate 3-D CPU executable built with `IBM_VALIDATION=FALSE` also compiles,
links, and advances the same restart by one complete four-stage Runge-Kutta
step with strict state positivity enabled. Mode 8 therefore has no runtime or
type dependency on validation-only audit storage.

## Remaining gates

- fresh statistically stationary mode-8 Mach-4 sphere refinement, including
  surface/interface load convergence;
- production qualification (or rejection) of the same-side shifted WENO/TENO
  shell currently used by the first-layer accuracy gate;
- full ratio-3 AMR budget and restart reproducibility with mode 8; and
- strong-shock long-time positivity monitoring beyond the short local smoke
  tests performed here.
