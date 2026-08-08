# Cerisse IBM Scope Lock

## Non-negotiable research scope

All current and future production IBM work in this repository targets a pure,
sharp-interface ghost-point / ghost-cell immersed boundary method (GP-IBM).
This scope is a hard project requirement. Do not replace it, silently augment
it, or describe it as an embedded-boundary, cut-cell, cut-control-volume, or
aggregate finite-volume method.

Both Cartesian and axisymmetric R--Z formulations are required production
targets. The generic R--Z metric, geometric-source, axis-parity, and GP-IBM
paths are not disposable experimental branches. Preserve them when cleaning or
reorganising the source. In R--Z, the pure-GP requirement still means evolution
on complete background-grid cells without cut-cell geometry; it does not mean
replacing the coordinate-aware R--Z operator with a Cartesian operator.

Production source code must remain independent of individual calculations.
Do not add numerical branches, compile-time macros, hard-coded indices, or
diagnostic ownership rules for Run165, a particular nozzle, circle, MMS, or any
other named case. Case setup and case-specific instrumentation belong under
`exm/` or `tst/` and must call generic production interfaces. Run165 is one
downstream validation case, not the design basis for the R--Z or GP-IBM source
and not a routine regression unless the user explicitly requests it.

The frozen production candidate is:

- one unique shared ghost-point state at each real solid-side ghost-cell center;
- cell-average-aware, boundary-intercept-constrained weighted least squares
  (BI-CWLS) for the ghost-state extension;
- curvature-compatible normal-pressure closure for smooth Euler-slip walls;
- one-sided fluid entropy extension and thermodynamically consistent state
  construction;
- full background-grid LLF-WENO-Z5 flux divergence, using the Cartesian form
  in Cartesian coordinates and the metric/source-consistent form in R--Z;
- a conservative shared-background-face positivity limiter;
- SSPRK(4,3), with every forward-Euler bracket checked independently;
- stage-local characteristic physical-boundary ghost filling (NSCBC/CBC).

Every active fluid cell must be advanced with its full Cartesian volume and
full Cartesian face areas. A fluid-solid stencil crossing obtains its solid-side
state exclusively from the shared GP extension. The standard production RHS is

    dU_i/dt = -(1 / |C_i|) sum_f |A_f| Fhat_f,

with |C_i| = product(dx) and full Cartesian |A_f|. The R--Z formulation uses
the corresponding full background-cell metric balance and geometric source.
Neither coordinate system may introduce cut volume or open-face fractions.

## Prohibited production-path mechanisms

The following mechanisms must never enter the pure GP-IBM production RHS,
state update, positivity update, force definition, checkpoint state, or AMR
flux register:

- cut-cell or embedded-boundary flow updates;
- cut-control-volume residuals;
- fluid volume fractions or cut-cell volume fractions;
- open-face, area-fraction, or subface-fraction fluxes;
- physical wall-segment fluxes used to replace the Cartesian RHS;
- aggregate control volumes, aggregate topology, gather/scatter evolution, or
  state/flux redistribution;
- aggregate or cut-control positivity limiting;
- AMReX EB as a flow discretization or geometric conservation correction.

Experimental code implementing any item above may remain only in a clearly
isolated experimental location. It must be default-off, must require an
explicit experimental opt-in, and must abort immediately when requested by a
pure-GP production case. Results from such a method are hybrid GP/cut-cell
results and must never be cited as verification or validation of pure GP-IBM.

## Pure GP-IBM force and conservation semantics

For pure GP-IBM:

- local wall pressure and Cp come from the boundary-intercept surface recovery
  (surface reconstruction order may differ from the volume ghost extension);
- physical integrated pressure force comes from quadrature of the recovered BI
  pressure over the immersed surface;
- Cartesian fluid-solid stencil-crossing momentum exchange is a discrete
  numerical diagnostic, not pointwise physical wall traction;
- disagreement between BI-integrated force and Cartesian momentum exchange is
  a reported discretization uncertainty, not a reason to introduce cut-cell
  geometry into the production method;
- do not enforce zero mass flux independently on each Cartesian crossing face;
- do not hide leakage or force disagreement with global sources, state
  clipping, pressure clipping, or a posteriori force correction.

Pure full-cell GP-IBM is not claimed to conserve exactly on the true curved
fluid control volume. Do not claim machine-precision true-domain conservation
or exact equality of Cartesian and BI forces for this method.

## Required safeguards and regression gates

Before changing any IBM production code:

1. State explicitly how the change remains within pure shared-GP IBM.
2. Identify whether it affects GP geometry, BI-CWLS, thermodynamic extension,
   Cartesian flux reconstruction, the full-cell positivity limiter, surface
   recovery, or only diagnostics.
3. Confirm that cut-control, cut-cell, aggregate, and EB flow paths remain
   disabled and unreferenced by the intended executable.
4. Preserve a runtime method manifest that unambiguously reports
   `pure_shared_gp_full_cartesian` for production runs.
5. Fail closed if a pure-GP run requests any prohibited mechanism.

Minimum non-regression sequence after an IBM production change:

1. current pure-GP build and configuration manifest;
2. short field-by-field regression against the frozen pure-GP Mach-4 circle;
3. one PM expansion and one attached compression-corner regression;
4. smooth curved-wall exact-solution/MMS checks, including first-wall-adjacent
   cells and BI pressure;
5. inverse-MOC phase checks using the fluid field, BI pressure, BI-integrated
   pressure force, and mass-flow variation. Do not require Cartesian/CV force
   equality as a pure-GP acceptance criterion.

Any loss of the detached bow shock or a material circle standoff regression is
a stop condition. Diagnose and restore that baseline before unrelated work.

## Change-control rule

Do not delete, migrate, or reactivate existing experimental cut-control code
without explicit user approval. Do not use it to fix a pure GP-IBM validation
failure. New algorithmic work outside the frozen pure-GP scope requires the
user to change this file explicitly first.
