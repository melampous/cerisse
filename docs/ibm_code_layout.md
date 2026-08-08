# Pure GP-IBM code layout

This document describes the frozen production method, not every IBM-related
prototype present in the working tree. The code-history baseline is Salvador's
unmodified `upstream/fsi` commit
`7b50e3cc9666c867f8f3c301f9502ab7aeef4f90`.

## Production contract

The production method is `pure_shared_gp_full_cartesian`:

- one unique reconstructed state at each real solid-side ghost-cell centre;
- cell-average-aware, BI-constrained WLS (BI-CWLS);
- curvature-compatible normal-pressure closure for smooth Euler-slip walls;
- one-sided fluid entropy extension and EOS-consistent state construction;
- full Cartesian cell volumes and full Cartesian face areas;
- characteristic LLF-WENO-Z5 fluxes;
- shared-Cartesian-face positivity limiting;
- SSPRK(4,3), with every Forward-Euler bracket checked;
- stage-local characteristic physical-boundary ghost filling.

No cut-cell volume, open-face fraction, physical wall-segment residual,
aggregate state, or redistribution term belongs to this production update.

## Pure-GP production modules

- `src/ibm/ibm_method_config.h`
  - prints the unambiguous runtime method manifest;
  - rejects cut-control input keys in a pure-GP executable.
- `src/ibm/ibm_bvh*.h`, `ibm_backend_*.h`
  - closest-point, first-intersection, inside/outside, and visibility geometry.
- `src/ibm/ibm_containers.h`
  - GP, BI, support, surface, and reconstruction storage.
- `src/ibm/ibm_solver.h`
  - level orchestration, marker construction, GP initialization, stage-local
    ghost-state reconstruction, and surface recovery.
- `src/ibm/ibm_solver_interp.h`
  - cell-average moments, BI-CWLS, constrained thermodynamic extension, and
    admissible fallback.
- `src/ibm/ibm_rz_moments.h`
  - annular-average moment functionals used by the R-Z extension.
- `src/ibm/ibm_positivity_limiter.h`
  - full-cell, shared-face positivity limiting and stage admissibility.
- `src/ibm/ibm_solver_amr_support.h`
  - support-halo and coarse-fine coverage checks. This is an audit/safety
    dependency; it does not change the full-cell update.
- `src/ibm/ibm_solver_io.h`
  - geometry input and BI surface output.
- `src/ibm/ibm_solver_surface_load.h`
  - BI pressure/traction integration and Cartesian crossing-momentum
    diagnostics. The latter is not physical pointwise wall traction.
- `src/rhs/Weno.h`, `src/rhs/IBMSharedGPFluxUtils.h`
  - the certified path uses the unique shared-GP state and the historical
    solid-safe LLF-WENO-Z5 reconstruction;
  - common marker and Cartesian-stencil helpers used by that path live in
    `IBMSharedGPFluxUtils.h`.
- `src/tim/advance.cpp`, `src/tim/compute_rhs.cpp`
  - SSPRK orchestration, stage-local GP refresh, full Cartesian flux
    divergence, and final limited-face-flux accounting.
- `src/set/nscbc.h`, `src/set/nscbc_ghost.cpp`
  - stage-local characteristic physical-boundary states.

## Production extensions, not part of the frozen Euler-slip certification

These modules may be useful, but their results must not be attributed to the
frozen Gate 0--3 method unless the selected case explicitly validates them:

- AFD-HLLC IBM support in `Afd.h`, `AfdIBM.h`, and `AfdIBMDrivers.h`;
- wall-location-consistent inviscid and diffusive direct-flux closures;
- no-slip viscous/thermal near-wall reconstruction;
- overlap communication;
- AMR/reflux execution beyond the audited support checks;
- general three-dimensional STL feature handling.

## Historical IBM prototypes requiring isolation

The following default-off paths are not the frozen smooth-wall shared-GP
method:

- face-local Euler-slip states;
- shifted same-side fluid shells;
- wall-matched or Cartesian conservative crossing-face replacement;
- patch-conservative AFD crossing modes;
- direct-crossing ablations used by earlier MMS studies.

Their implementation is now separated from the shared-GP helpers:

- `IBMFaceLocalExperimentalUtils.h` contains direct/face-local inviscid
  crossing utilities;
- `IBMDirectDiffusiveExperimentalUtils.h` contains direct viscous/thermal
  crossing utilities;
- `WenoIBMSameSideExperimental.h` contains shifted same-side and face-local
  WENO/TENO class fragments;
- `WenoIBMDirectWallExperimental.h` contains direct stationary-wall class
  fragments.

`IbmFluxUtils.h` and `WenoIBMReconstruction.h` remain compatibility umbrellas
only. No production translation unit should include the broad
`IbmFluxUtils.h` umbrella. Enabling a prototype must be explicit in the
executable manifest.

`Weno.h` still includes the small `WenoIBMReconstruction.h` member-fragment
umbrella, so these default-off experimental members remain visible to the
compiler. The present cleanup separates their implementation and method
identity; it does not yet redesign `weno_t` to remove them from the production
compile graph.

## Experimental hybrid code

`src/ibm/ibm_cut_control_geometry.h` and all aggregate/cut-control update
branches implement a hybrid GP/cut-control method. They are not pure GP-IBM.
They may compile only with
`CERISSE_ENABLE_EXPERIMENTAL_IBM_CUT_CONTROL`, are default-off, and must not
enter a pure-GP RHS, positivity update, checkpoint, force definition, or flux
register.

## Validation-only code

`IBM_VALIDATION=TRUE` defines `CNS_IBM_VALIDATION` and adds
`tst/ibm/validation/include/` to the include path. Files there may:

- compare against analytic states or exact annular averages;
- audit visibility, rank, conditioning, and stencil use;
- record GP, face-flux, R-Z axis, or SSPRK-stage CSV data;
- poison/perturb a state for dependency tests;
- compute A/B/C operator decompositions.

They must not select a production state or flux. Production binaries do not
add that directory to the include path and must not allocate its managed
arrays or execute its scans. Guarded call sites may remain in `src/`, but the
test implementation belongs under `tst/`.

## Force and conservation semantics

For pure GP-IBM:

- local wall pressure and `Cp` come from BI surface recovery;
- physical integrated pressure force comes from BI quadrature;
- Cartesian fluid-solid stencil-crossing momentum exchange is a discrete
  full-cell-domain diagnostic;
- disagreement between those forces is reported as discretization
  uncertainty, not corrected with cut-cell geometry;
- zero mass flux is not imposed independently on each Cartesian crossing face.

## Current refactoring boundary

The mixed utility and WENO class-fragment headers have been split and
validation implementations have moved under `tst/`. The retired
boundary-face cache, projected-load correction, and momentum-budget CSV paths
have been removed. The frozen shared-GP Euler RHS obtains solid-side stencil
states only from the shared GP field.

Further safe refactoring must preserve the frozen shared-GP state, BI-CWLS
coefficients, thermodynamic extension, WENO face flux, SSPRK stages, and
surface-e2 observation. Any numerical change requires the pure-GP regression
sequence defined in `AGENTS.md`.
