# GP-IBM code reorganization audit

Date: 2026-07-23

## 1. Authoritative baseline

The comparison baseline is Salvador's original FSI branch:

```text
remote:  upstream = git@github.com:salvadornm/cerisse.git
ref:     upstream/fsi
commit:  7b50e3cc9666c867f8f3c301f9502ab7aeef4f90
subject: rename ibm_wallmodel.h back to ibm_walltypes.h
```

The local clean `HEAD` is not the baseline. Results produced earlier by
comparing only against local `HEAD` are superseded by this audit.

At audit time:

```text
branch:  ibm_fsi
HEAD:    dc35646d8eb90aaa362a3292a9d702e30c9d5af5
distance from upstream/fsi: 18 commits
```

The working tree is intentionally dirty and contains user work. This audit
does not revert or overwrite unrelated changes.

## 2. Size of the change

Approximate source-line counts:

| Scope | Salvador baseline | Current tree | Increase |
|---|---:|---:|---:|
| `src/` C/C++ sources and headers | 21,973 | 62,209 | 40,236 |
| `src/ibm/` sources and headers | 7,598 | 29,178 | 21,580 |

Tracked source differences relative to `upstream/fsi` are approximately:

```text
43 files, +25,685 / -2,044 lines
```

The current tree also has 16,876 lines in untracked files below `src/`.
After the header split, the largest untracked modules are:

| File | Lines | Current classification |
|---|---:|---|
| `ibm_cut_control_geometry.h` | 5,162 | hybrid experiment |
| `AfdIBM.h` | 1,790 | uncertified extension |
| `ibm_solver_boundary_faces.h` | 1,622 | mixed direct-flux/diagnostic metadata |
| `ibm_positivity_limiter.h` | 1,449 | pure-GP production |
| `WenoIBMSameSideExperimental.h` | 1,169 | experimental WENO/TENO fragments |
| `ibm_solver_surface_load.h` | 1,010 | mixed physical observation/diagnostic |
| `WenoIBMDirectWallExperimental.h` | 967 | experimental direct-wall fragments |

`WenoIBMReconstruction.h` is now a 9-line compatibility umbrella rather than
the former 2,128-line mixed implementation. Validation implementation headers
are no longer below `src/`; they are stored under
`tst/ibm/validation/include/`.

The main tracked growth is concentrated in:

| File | Added | Removed |
|---|---:|---:|
| `ibm_solver.h` | 5,847 | 772 |
| `compute_rhs.cpp` | 2,905 | 18 |
| `ibm_solver_interp.h` | 2,835 | 173 |
| `advance.cpp` | 1,887 | 131 |
| `ibm_containers.h` | 1,868 | 51 |
| `Weno.h` | 1,816 | 131 |
| `CNS.cpp` | 1,762 | 146 |
| `IbmFluxUtils.h` | 1,392 | 0 |

This growth is too large to review as one IBM change. It represents several
different numerical methods and substantial diagnostic infrastructure.

## 3. What was added after Salvador's baseline

### 3.1 Pure shared-GP functionality to retain

These features form the frozen production candidate or its R-Z extension:

- unique shared ghost state at each real solid-side ghost-cell centre;
- cell-average-aware BI-CWLS;
- expanded visible support with rank and conditioning checks;
- curvature-compatible pressure closure for smooth Euler-slip walls;
- one-sided fluid entropy extension and EOS-consistent conservative state;
- surface-e2 BI pressure and integrated pressure-force observation;
- full-cell, shared-Cartesian-face positivity limiting;
- SSPRK stage rollback/retry and admissibility accounting;
- stage-local characteristic boundary ghost filling;
- AMR image-point support coverage auditing;
- R-Z annular moments and centre-state recovery;
- BVH robustness, geometry provenance, and runtime method manifest.

These algorithms must not be changed by a code-layout cleanup.

### 3.2 Useful extensions to retain separately

These are not part of the frozen two-dimensional Euler-slip qualification:

- AFD-HLLC with WENO-Z5/TENO reconstruction;
- no-slip inviscid, viscous, and thermal direct-flux closures;
- LES/viscous compatibility fixes;
- overlap communication;
- three-dimensional and AMR support beyond the currently certified scope.

They should be kept under explicit method/capability labels. Their tests must
not be used as evidence for the frozen LLF-WENO-Z5 shared-GP method.

### 3.3 Validation and diagnostic infrastructure

The following belongs under `tst/ibm/validation/include` or in test tools:

- GP visibility/rank/condition CSV output;
- poison and finite-perturbation dependency tests;
- WENO split-flux, coefficient, and A/B/C operator audits;
- R-Z axis-RHS and shock-face CSV decomposition;
- analytic GP/component oracles;
- full-domain SSPRK-stage observer scans;
- one-shot diagnostic flux dumps.

These observers are scientifically valuable and should not be deleted merely
because they are large. They must, however, be absent from production binaries.

### 3.4 Historical numerical prototypes

The following paths were explored but are not the frozen shared-GP method:

- face-local Euler-slip state construction;
- shifted same-side WENO/TENO shells;
- wall-matched and Cartesian conservative crossing-face replacement;
- patch-conservative AFD modes;
- several direct-crossing validation ablations.

Repository input search shows nonzero face-local/shifted settings primarily in
validation cases. Production Gate 0--3 configurations explicitly disable them.
They are candidates for migration to an `experimental` or `extensions`
directory after their common utilities are split from the shared-GP path.

### 3.5 Hybrid cut-control experiment

The cut-control database, open subfaces, aggregate states, redistribution, and
aggregate positivity limiter implement a hybrid GP/cut-control method. They
are not GP-IBM production functionality.

The current build guard is correct:

```text
IBM_EXPERIMENTAL_CUT_CONTROL=FALSE
```

removes the database type and most call sites from a pure build. However, the
guarded branches remain interleaved through `CNS.cpp`, `advance.cpp`,
`compute_rhs.cpp`, `ibm_solver.h`, and the positivity limiter. This makes the
pure path hard to review. Moving these branches behind one experimental
adapter is recommended, but deletion or migration requires explicit approval.

## 4. Refactoring decisions

### 4.1 Completed low-risk moves

These changes do not alter a reconstructed state, numerical flux, or RK
coefficient:

1. Moved the full-cell shared-face positivity implementation from a
   validation-named header into `src/ibm/ibm_positivity_limiter.h`.
2. Kept SSPRK stage observation in
   `tst/ibm/validation/include/ibm_stage_positivity.h`.
3. Guarded `first_fe_bracket_flux_csv` with `CNS_IBM_VALIDATION`.
4. Moved the R-Z shock-face diagnostic descriptor to
   `tst/ibm/validation/include/rz_shock_face_diagnostic.h`.
5. Compiled R-Z axis/shock CSV capture only in validation builds.
6. Replaced the stale code-layout document with the frozen pure-GP contract.
7. Verified both sides of the compile-time boundary:
   - the RZ wall-axis validation executable compiled and linked with
     `CNS_IBM_VALIDATION`;
   - a separate CPU circle executable compiled and linked with
     `IBM_VALIDATION=FALSE`, without `CNS_IBM_VALIDATION`.

The production executable was

```text
build id:
sm89.release.dp.llf_wenoz5.e2.i2.a1p0.g1.slip.uniform.shared.gp_ip.
pshock_l0p03_h0p10.prod.src6156a0cd6676f051

SHA-256:
cf75845429329be9c59d927d6c7f6881b16b87c6a7f3e505fb9a33f53ab66cc5
```

Its binary contains none of the `first_fe_bracket`, `rz_axis_rhs`, or
`rz_shock_face` diagnostic strings. A one-step CPU smoke test completed and
printed:

```text
family=pure_shared_gp_full_cartesian
pure_gp_production=1
cut_control_compiled=0
validation=0
control_volume=full_cartesian_cell
face_area=full_cartesian_face
solid_extension=shared_ghost_point
```

This is a compile/startup isolation check, not a replacement for the numerical
Gate 0 regression.

### 4.2 Completed mixed-header split

The low-risk dependency split is complete:

```text
IbmFluxUtils.h
  -> IBMSharedGPFluxUtils.h
  -> IBMFaceLocalExperimentalUtils.h
  -> IBMDirectDiffusiveExperimentalUtils.h

WenoIBMReconstruction.h
  -> WenoIBMSameSideExperimental.h
  -> WenoIBMDirectWallExperimental.h
```

The certified shared-GP LLF-WENO-Z5 implementation remains in `Weno.h`;
`IBMSharedGPFluxUtils.h` contains only marker, stencil, and Cartesian helper
functions shared with that path. The other four implementation headers are
explicitly labelled experimental. The old two filenames remain as small
compatibility umbrellas so the class layout and include order did not have to
change in the same patch.

Definitions were moved mechanically without changing expressions, template
parameters, or call order. Production and validation CPU builds and the
CUDA-sm_89 production build link successfully.

This split establishes source ownership and prevents experimental methods from
being mistaken for the certified shared-GP implementation. It is not yet a
complete compile-graph isolation: `Weno.h` still includes the 9-line
`WenoIBMReconstruction.h` compatibility member-fragment umbrella, so the
default-off experimental member functions are parsed as part of `weno_t`.
Removing that dependency would require a separate class/API refactor and is
intentionally deferred until after this cleanup is committed.

### 4.3 Lazy boundary-face metadata

`buildBoundaryFaceAudit()` is no longer called for every non-CGAL GP setup.
The owning-face cache is invalidated with GP geometry and constructed by
`ensureBoundaryFaceData()` only when a consumer requests it. Consumers are:

- compile-time face-local/direct-flux extensions;
- wall-location-consistent inviscid or diffusive closures;
- strict boundary invariant checks;
- validation visibility/flux-defect reports;
- Cartesian crossing-momentum surface diagnostics.

The frozen shared-GP Euler update does not request this data. A production
one-step run with surface output disabled has no
`IBM::buildBoundaryFaceAudit` profiler entry. Enabling the surface crossing
diagnostic builds it once and reuses it across subsequent RHS calls.

The cache retains the original strict ownership checks whenever it is built.
Physical BI pressure recovery does not depend on it.

A stronger regrid test exposed one lifecycle defect after the initial split.
The first surface-load accumulation following a regrid requested the cache
from inside the per-FAB flux loop. The cache builder starts a level-wide
`MFIter`, so this produced a nested-`MFIter` assertion. The fix is deliberately
structural:

- `compute_rhs.cpp` ensures requested level-wide metadata before entering the
  per-FAB loop;
- the per-FAB surface-load function now asserts that the cache is ready;
- the level-wide post-limiter function remains responsible for its own
  pre-loop ensure.

No GP state, WENO flux, R-Z source, or SSPRK expression changed. Level-0 fields
and all non-rank VTP arrays are bitwise identical before and after this fix.

### 4.4 Validation implementation location

Validation implementation headers now live in:

```text
tst/ibm/validation/include/
```

`IBM_VALIDATION=TRUE` adds that directory to the build. Production builds do
not add it. The production tree retains only guarded integration call sites;
analytic oracles, poison/perturbation hooks, CSV observers, and R-Z diagnostic
descriptors are not stored under `src/`.

### 4.5 Do not merge the core files back together

`ibm_solver.h` is already about 7,000 lines. The extracted class-body fragments
are directionally correct. The target is smaller cohesive modules, not a
single monolithic header:

- shared-GP geometry and BI-CWLS;
- thermodynamic wall extension;
- R-Z semantics;
- physical BI surface observation;
- optional direct-flux extensions;
- validation observers.

### 4.6 Deletion policy

No numerical prototype is deleted in this audit. Deletion is justified only
after:

1. its compile-time and runtime callers are listed;
2. no retained case requires it;
3. its last scientifically useful result is archived;
4. shared helpers are extracted;
5. the user explicitly approves deletion.

Generated binaries, plotfiles, VTP/CSV histories, and images are a separate
storage cleanup and must not be mixed with source refactoring.

## 5. Performance findings

Potential avoidable production overheads are:

- unconditional sparse crossing-face cache construction;
- production compilation of validation-only R-Z CSV code before this audit;
- large mixed headers that increase NVCC compile time and template memory;
- optional Cartesian crossing-load accumulation being coupled too closely to
  physical BI surface output;
- repeated runtime option parsing inside method-specific flux templates.

The hot numerical kernels were not rewritten during this audit. Performance
changes require measurement, not line-count assumptions.

## 6. Cleanup regression results

The refactoring was tested with the frozen pure-GP method manifest. No test
enabled cut-control, aggregate, open-subface, volume-fraction, or EB flow
updates.

| Gate | Configuration | Result |
|---|---|---|
| build | CPU production, CPU validation, CUDA `sm_89` production, CUDA `sm_89` MMS validation | PASS |
| lazy cache | current circle executable, cache unrequested versus explicitly requested by the read-only audit | all plotfile fields bitwise identical |
| Mach-4 circle | `N=320`, shared GP, LLF-WENO-Z5, 12-step field comparison | all fields bitwise identical between lazy-cache modes |
| Mach-4 circle long run | `N=320`, `D/h=64`, `t=0.008` (`t*=55.56`) | detached bow shock retained; `Delta/R=0.5496841170717145` |
| PM expansion | `N=240`, `M=2`, `20 deg`, frozen Gate 0 protocol | all 12 fields bitwise identical to the frozen plotfile |
| compression corner | `N=120`, `M=2`, `20 deg`, frozen Gate 0 protocol | all 12 fields bitwise identical to the frozen plotfile |
| smooth curved-wall MMS | `N=96,128,192`, axis phase and `(0.37h,0.23h)` phase | all six terminal plotfiles bitwise identical to frozen Gate 1 |
| inverse-MOC | two geometry phases, `256x64`, fixed post-cleanup executable | all principal flow metrics unchanged; maximum of 413 numeric-field differences `1.43e-8`, surface-force difference `2.96e-9` |
| RZ-0 | annular BI-CWLS moments, CPU single/multi-FAB, MPI-2, CUDA | polynomial and constrained moment tests PASS |
| RZ-1 | no-solid R-Z operator with explicit centre-state recovery | frozen pressure-balance and convergence values reproduced |
| RZ-2a | two phases, four grids, RHS/fixed-time plus surface recovery | summaries reproduce the frozen PASS |
| RZ-2b | wall-axis junction matrix | cleanup result reproduces the frozen summary; known first-axis-ring and pole-pressure limits remain |
| cache lifecycle | level-0, MPI-2, AMR regrid, checkpoint/restart | PASS; build counts `1`, `1`, `4`, `1`, respectively |
| Salvador clean reproduction | detached worktree at `upstream/fsi`, production and validation builds | PASS; clean binary and current binary produce bitwise-identical short fields and VTP output |

The regenerated three-grid MMS fits retain approximately second-order global
and wall-adjacent convergence. In the axis phase, the fitted `L2` orders are
`1.99/2.00/2.00/1.98/2.00` for
`rho/u/v/p/T`; the first-wall-adjacent velocity and pressure orders are
`2.00` and `1.84`. In the shifted phase, the corresponding global orders are
`2.00/2.00/2.00/2.01/2.01`, with wall-adjacent velocity and pressure orders
of approximately `2.01` and `1.88`. These are short three-grid cleanup
regressions; the authoritative five-grid formal-order statement remains the
frozen Gate 1 report.

The PM and compression-corner postprocessors also reproduce the frozen
physical diagnostics:

- PM downstream turning-angle error: `0.000858 deg`;
- PM pressure-ratio error: `0.006%`;
- compression shock-angle error: `-0.051865 deg`;
- compression pressure-ratio error: `-0.003%`.

The long-run circle stand-off differs from the frozen
`0.5496841170717` value by approximately `2.6e-14` in relative terms. The
density- and pressure-gradient shock locators, fitted front shape, and shock
thickness are unchanged. A whole-domain terminal plotfile comparison is not
bitwise because the aft unsteady wake has decorrelated over more than 15,000
steps; this is reported explicitly and is not used to weaken the detached
bow-shock stop criterion.

The production run with `ib.plot_surf=0` has no
`IBM::buildBoundaryFaceAudit` profiler entry. Validation and Cartesian
crossing-load consumers still construct the cache exactly once.

The formal cache-lifecycle result is
`tst/ibm/regression/results/code_cleanup_20260723/cache_lifecycle/formal_r2/summary.json`.
For the production build:

- three level-0 surface-output windows use one cache construction;
- `np=1` and `np=2` fields and VTP arrays are bitwise identical;
- a two-level run with regridding at every coarse step performs four expected
  level constructions and retains 48 crossing faces in both windows;
- a checkpoint restart rebuilds once and matches the continuous plotfile
  bitwise.

Validation executables intentionally request boundary-face data eagerly for
visibility reporting. Their equivalent regrid run therefore records six
constructions, not the production lazy-cache count of four. This is a
validation lifecycle distinction, not a production regression.

The R-Z cleanup regression must be interpreted carefully. RZ-2a remains a
PASS. RZ-2b is an exact non-regression against its frozen result, but its
overall summary remains `pass=false` because the previously known first-axis
ring and phase-0 pole-pressure orders are below the formal acceptance
threshold. Consequently:

```text
RZ-2b GP-specific cleanup non-regression: PASS
complete R-Z first-axis-ring consistency: CONDITIONAL / pending A-B-C-D audit
```

RZ-1 also has a provenance issue: the reported 2.26--2.40 axis-order sequence
is reproduced only when
`cns.rz_euler_cell_average_deconvolution=1` is supplied explicitly. The
default matrix driver does not currently encode that requirement. The
numerical operator was not changed during this cleanup; the runner and report
must be made unambiguous before the next R-Z numerical study.

### 6.1 Clean-worktree reproduction

An independent detached worktree was created from:

```text
7b50e3cc9666c867f8f3c301f9502ab7aeef4f90
```

Only the current `src/` patch, required untracked pure-GP source headers,
test-only headers, and a minimal circle fixture were overlaid. The
experimental `ibm_cut_control_geometry.h` was deliberately omitted. The
production and validation executables both compiled without it, proving that
the pure build has no implicit dependency on that hybrid header.

```text
source patch SHA-256:
97b654ed728809c74fe8cbfb04b915cbbded326b592bd1c801601ddc67377ad5

AMReX:
bd922c6216e0a734f3b1cf0ca73e7d669b90f3ef (clean)

clean production binary SHA-256:
18dc701c67948787c47a46b3bf3f4ee846cff3ae837a21f48b9d90b5f441aace

clean validation binary SHA-256:
f1a1b75f3667d495c2a66f74cbe1aa7518cad4876a1649cd943811597c3414cc
```

The clean production binary reproduces the current-tree production binary
bitwise for all short-run fields and all VTP arrays except the intentionally
excluded MPI `Rank` field. Production and validation binaries are also
bitwise identical when validation observers are disabled. The generated
summary is stored at
`tst/ibm/regression/results/code_cleanup_20260723/clean_reproduction/summary.json`.

This proves that a self-contained source bundle can be reconstructed from the
Salvador baseline. It does not replace version-control integration: required
new production headers and `tst/ibm` validation sources are still untracked in
the current worktree and must be added in the eventual cleanup commit.

## 7. Required regression sequence

After any production-source split:

1. print `pure_shared_gp_full_cartesian` and confirm cut-control is absent;
2. compare a short Mach-4 circle run field by field;
3. run one PM expansion and one attached compression corner;
4. run the smooth curved-wall first-layer/BI-pressure check;
5. run inverse-MOC phase checks using BI pressure, BI force, and mass flow.

Detached bow-shock loss or a material standoff regression stops the cleanup.

## 8. Current recommendation

Freeze both the numerical method and the present file split. The cleanup
gates now have the following status:

```text
code structure:                         PASS
Cartesian non-regression:               PASS
R-Z cleanup non-regression:             PASS
boundary-face cache lifecycle:          PASS
Salvador clean-worktree reproduction:   PASS
Git integration of new source/tests:    PENDING
```

Do not continue moving headers before the cleanup is committed. Do not delete
the cut-control or face-local prototypes without explicit approval, and do
not use them to repair a pure shared-GP validation failure. The next numerical
task, in a separate change, is the frozen first-axis-ring pressure
flux/geometric-source A/B/C/D audit.
