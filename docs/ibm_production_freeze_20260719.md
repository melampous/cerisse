# GP-IBM production-route freeze (2026-07-19)

## Purpose and change control

The project is frozen around the intended application rather than around a
claim that every Cerisse configuration is validated. Until this document is
superseded, IBM changes are limited to:

1. defects that block an already defined acceptance gate;
2. boundary conditions, geometric features, and R--Z discretization required
   by the nozzle/SRP application;
3. verification, diagnostics, conservation budgets, and provenance.

New ghost-state closures, WENO/TENO variants, Riemann solvers, interpolation
orders, and geometry algorithms are outside the frozen route unless an
existing gate proves that they are necessary.

## Phase-0 historical baseline

The pre-integration reference is a 2-D Cartesian, level-0, fixed-circle,
single-species ideal-gas Euler-slip calculation with LLF-WENO-Z5 and
SSPRK(4,3). It uses shared GP states, BI-constrained reconstruction, expanded
BIC support, the finite-volume normal-momentum correction, one-sided entropy
jet, entropy wall extension, and analytic-curvature pressure compatibility.
The production positivity limiter and all validation oracles are disabled.

Frozen provenance:

| Item | Value |
| --- | --- |
| Git HEAD | `dc35646d8eb90aaa362a3292a9d702e30c9d5af5` (dirty worktree recorded separately) |
| Tracked-diff SHA-256 | `c8a324eaf3959286a463a769aff1bcf36336a88d648cbeb5426c019039212168` |
| Solver source-content ID | `5d0decac5024c791` |
| AMReX commit | `bd922c6216e0a734f3b1cf0ca73e7d669b90f3ef` |
| `prob.h` SHA-256 | `af9f913c56ccf7c0d6466a63ce5eac66c2b66a06a0898816f849f7a091288d5e` |
| GNUmakefile SHA-256 | `32026d44485801c94ec94fb3655e5e1f7ffc1829d0830912c4badb8c3a2b8ad0` |
| inputs SHA-256 | `81a058287a7db070cb4f717ff73c2709f1ea2f6e8023d79ca6005542863fb9c9` |
| 4096-edge circle SHA-256 | `53775b9fc2a2b67a354b72ca4fcbb746810656d98a926c5effec123abb424e6e` |
| CUDA binary SHA-256 | `87d69175102d80b0d83c16b4fb1e102318b3b571197ea965ba8b80a22778ba70` |
| CUDA target | CUDA 12.0.140, `sm_89`, NVIDIA RTX 4090 Laptop GPU |

The reference run uses `320 x 320`, `max_level=0`, `max_grid_size=256`,
`blocking_factor=64`, and 150 steps. Its plotfile is:

`IBM/cases/ibm_tests/2d_bvh_gpu/runs/ibm_production_freeze_20260719/preintegration_default_off_N320/plt00150`

Reference payload hashes are:

| Plotfile object | SHA-256 |
| --- | --- |
| `Header` | `e0a934faa4672373cf66617acca91e8c4d187678315252daafeaf0ff8dd7473a` |
| `Level_0/Cell_H` | `519136651b6617c4225d0c223cf6bf81036057c0f47d7b39445cb554ebc9f589` |
| `Level_0/Cell_D_00000` | `09df2e28230b93a04a76873fc756607ef3be2e68e1b7a7ebacae0fffb25f0c8e` |

### Baseline supersession after deterministic closest-facet ownership

The plotfile above is retained as a historical, pre-fix reference. It is not
the current bitwise reference. A Phase-1 CPU/CUDA audit found one ambiguous
ghost cell, with global Cartesian index `(137,137)`, whose boundary intercept
lies at the common vertex of polygon elements 2559 and 2560. The two squared
closest-point distances agree to roundoff. The previous BVH search selected
the owning surface element according to traversal order and floating-point
roundoff: different backends could therefore assign a different normal and a
different constrained ghost state to the same ghost-cell centre.

The closest-point query now applies a roundoff-scaled distance tie and selects
the lower primitive ID on a tie. BVH pruning uses the same tolerance. This is
an intentional geometry-semantics correction, not an effect of the positivity
limiter. On the fixed 4096-element circle, 179 of 180 ghost-cell owners remain
unchanged and the one vertex tie changes from element 2560 to element 2559.
That single change is sufficient to separate the later nonlinear shock
trajectory, so comparison with the pre-fix plotfile cannot satisfy a bitwise
criterion.

Consequently:

1. the pre-fix reference remains available only for historical regression;
2. the deterministic closest-facet implementation defines the new Phase-0
   geometry semantics;
3. the limiter-off bitwise gate compares builds and executions using those
   identical deterministic geometry semantics;
4. changing the tie tolerance or primitive-ID rule requires an explicit new
   geometry baseline and cannot be hidden inside an IBM algorithm comparison.

The case build directory now includes the complete source-content ID for both
audit and production builds. This closes a separate stale-object risk in which
production object files could previously be reused across source changes.

## Phase-1 Cartesian positivity contract

For each SSPRK(4,3) Forward-Euler bracket, define the high- and low-order
numerical face fluxes from the identical stage state:

\[
U_i^H=U_i^{(s)}-\frac{\Delta t_{FE}}{V_i}
\sum_f A_f s_{if}F_f^H,
\qquad
U_i^L=U_i^{(s)}-\frac{\Delta t_{FE}}{V_i}
\sum_f A_f s_{if}F_f^L.
\]

The low-order flux is piecewise-constant Rusanov with
\(\alpha_f=\max(|u_{n,L}|+c_L,|u_{n,R}|+c_R)\), using admissible fluid and
GP-IBM states. The actual all-low-order cell update, not only its nominal CFL,
must satisfy

\[
\rho>\epsilon_\rho,
\qquad
\rho e=\rho E-\frac{|\rho\mathbf u|^2}{2\rho}>\epsilon_e.
\]

The antidiffusive contribution of face \(f\) to cell \(i\) is

\[
P_{if}=-\frac{\Delta t_{FE}A_f}{V_i}s_{if}
       (F_f^H-F_f^L).
\]

For a Cartesian cell with \(N_f=2d\), each cell-face bound is chosen so that

\[
U_i^L+N_f\theta_{if}P_{if}\in\mathcal G_\epsilon.
\]

An internal face uses one shared coefficient

\[
\theta_f=\min(\theta_{if},\theta_{jf}),
\]

and one shared limited flux

\[
F_f^{lim}=F_f^L+\theta_f(F_f^H-F_f^L).
\]

The final update is therefore the convex average

\[
U_i^{lim}=\frac{1}{N_f}\sum_f
\left(U_i^L+N_f\theta_fP_{if}\right)
\in\mathcal G_\epsilon.
\]

The four SSPRK(4,3) Forward-Euler brackets use \(\Delta t_{FE}=\Delta t/2\).
In particular, bracket 3 is limited as

\[
U_{FE}^{(3)}=U^{(2)}+\frac{\Delta t}{2}L(U^{(2)}),
\qquad
U^{(3)}=\frac23U^n+\frac13U_{FE}^{(3)}.
\]

Limiting the already combined stage with an effective `dt/6` is not the SSP
construction and is prohibited.

## Fail-closed hierarchy

Every fallback is reconstructed from the same
\(U^{(s)},F^H,F^L,U^L\); fallbacks are never applied cumulatively:

1. retain the high-order bracket if it is admissible;
2. apply the local shared-face convex limiter;
3. if the final audit fails, reconstruct a global-\(\theta\) face-flux blend;
4. if that fails, reconstruct the all-low-order update;
5. if the all-low-order update is inadmissible, abort the step without state
   clipping.

AMReX currently advances its coarse time after `CNS::advance()` returns, so
returning a smaller `dt` from inside `advance()` is not a valid rejection
mechanism. Transactional whole-step rollback/retry must restore the state,
stage time metadata, GP cache, and accumulated flux data before automatic time
step reduction can be enabled.

## Current opt-in support boundary

The phase-1 implementation must terminate for unsupported configurations. Its
initial supported set is:

- 2-D Cartesian coordinates;
- one uniform level, no reflux;
- source-free, inviscid, single-species constant-gamma Euler equations;
- fixed stationary smooth Euler-slip GP-IBM wall;
- LLF-WENO-Z5 and SSPRK(4,3);
- non-periodic boundaries, NSCBC disabled;
- conservative IBM load and surface output using the final limited
  crossing-face flux.

The conservative load window uses the finalized face flux from every SSPRK
bracket with weights `dt/6`, `dt/6`, `dt/6`, and `dt/2`. The same limited
crossing flux therefore enters the fluid update and the IBM momentum-exchange
load. Reconstructed BI pressure remains a separate local wall diagnostic.

## Phase-1 acceptance gates

1. Limiter disabled: the deterministic-geometry, limiter-off 150-step circle
   plotfile is bitwise identical across a clean rebuild of the same backend.
   The historical pre-fix plotfile is explicitly excluded from this gate.
2. Limiter enabled: the known Mach-4 circle failure is crossed without
   clipping; all four FE brackets and all SSP convex combinations are
   admissible.
3. Every limiter event records the all-low-order CFL, minimum \(\theta_f\),
   limited face count, crossing-face count, and conservation residual.
4. A multi-FAB and MPI run gives one shared flux and one shared \(\theta_f\)
   per face, with conservation residual at roundoff scale.
5. CPU and CUDA runs give the same limiter decisions and quantitatively
   equivalent states.
6. Sourced smooth-wall MMS retains its established convergence with the
   source-free limiter disabled. A source-free smooth exact solution with the
   limiter enabled has \(\theta_f=1\) in smooth regions, or activation vanishes
   under refinement. Source-term positivity is a separate future contract.
7. Mach-4 circle, PM expansion, single/double compression corners, and a
   no-solid near-vacuum Riemann problem pass before the method is used for a
   shocked nozzle.

Passing these gates permits a limited claim for 2-D uniform Cartesian
Euler-slip nozzle calculations. It does not certify R--Z, AMR, NS, moving
bodies, TENO/HLLC, or general 3-D STL geometries.

### Verification status at source-content ID `c68c0fa9bdd7c4e8`

All results in this subsection use `320 x 320`, `max_level=0`,
`max_grid_size=256`, `blocking_factor=64`, and the fixed 4096-element circle
whose SHA-256 is listed above.

| Gate | Result | Status |
| --- | --- | --- |
| Fixed geometry actually loaded | Run log reports `4096 -> 4096` elements | PASS |
| Four SSPRK(4,3) FE brackets | Production limiter is called on every bracket with `dt_FE=dt/2` | PASS (code and runtime audit) |
| Long-time CUDA admissibility | 15,436 steps and 61,744 FE brackets; 80 early local-limiter events; no clipping, global-theta, all-low, or inadmissible accepted state | PASS |
| Shared-face conservation | Worst reported relative residual is `1.138e-15` | PASS |
| CPU one-rank vs two-rank | Plotfiles agree exactly; event histories agree | PASS |
| CPU vs CUDA | Conservative variables differ by approximately `2.2e-4` to `3.5e-4` at step 200; limiter event histories differ | CONDITIONAL |
| Final limited flux used by IBM load | SSPRK-weighted crossing flux feeds both fluid update and conservative force; one-/two-rank load metadata agree | PASS |
| Wall/conservative drag closure | Final `Cd=1.24615/1.25162`, a `0.4375%` difference; load-consistent projection residual below `7e-15` | PASS for this grid/window |
| Historical pre-fix bitwise comparison | One intentionally changed owning element causes nonlinear trajectory separation | SUPERSEDED |
| Deterministic limiter-off clean-rebuild reference | Must be generated after the Phase-1 source is frozen | OPEN |
| Full Mach-4 circle at `t*=55.56` | `Delta/R=0.5496841171`; final four-sample span `1.056e-7`; density/pressure estimator span `0.1001h` | PASS |
| Automatic failed-step retry | All-low inadmissibility currently aborts without clipping; transactional rollback/retry is not implemented | OPEN |

The full circle report and machine-readable records are in
`IBM/cases/ibm_tests/2d_bvh_gpu/production_freeze_20260719/circle_N320_limiter_load`.
The current source-content build `75f7ee86428df61f` was rerun to the same
final time on 2026-07-20. AMReX `fcompare` reports zero error in every plot
variable at the final step, the final VTP is byte-identical, and all limiter
statistics and the stand-off estimate are unchanged. This closes the current-
source provenance concern; details are in
`IBM/cases/ibm_tests/2d_bvh_gpu/production_freeze_20260720/circle_N320_limiter_load_src75f7ee86428df61f/REPORT.md`.
The forebody stand-off is stationary, but late wake loads still drift by about
0.5 percent over the final five output windows; force statistics therefore
remain a separate convergence gate.

The CPU/CUDA result is quantitative agreement, not bitwise agreement. Before
the first limiter event, conservative-field differences are already of order
`2e-5` to `4e-5`; local admissibility decisions subsequently amplify this
difference. This does not invalidate the convex limiting construction, but it
does require the long-time circle gate and a documented cross-backend
tolerance before production certification.

### Source-free and canonical-wave gate update (audited 2026-07-20)

The frozen limiter and IBM method were subsequently exercised without any
algorithm change:

| Gate | Result | Status |
|---|---|---|
| Source-free stationary isentropic vortex with circular IBM, N=80/120/160/240/320 | all-fluid primitive L2 orders 2.01--2.09; first-layer orders 1.98--2.07; all 648 FE brackets retain `theta=1` | PASS |
| Mach-2 PM expansion | flow-angle error `-0.0051 deg`, pressure-ratio error `-0.038%`; no limiter event | PASS |
| Mach-2 single compression | shock-angle error `+0.0563 deg`, pressure-ratio error `-0.016%`; no limiter event | PASS |
| Sequential double compression | both pressure plateaus within `0.064%`; interaction location differs by about `4.0h` | CONDITIONAL PASS |
| No-solid near-vacuum double rarefaction | all 12,668 FE brackets admissible and conservative with `theta=1`; exact FV averages show broad star-region thermodynamic error for t=0 discontinuous startup; an exact-evolved `tau=0.01` startup reduces 512-grid star-region relative L1 errors in rho/p/T/sigma to 0.30/0.36/0.46/0.57% | admissibility PASS; limiter activation N/A; ab-initio star-region thermodynamic accuracy FAIL |

The near-vacuum defect is present with no solid or ghost cell and therefore
must not be attributed to the GP-IBM. It also shows that an admissibility
guarantee alone does not guarantee entropy accuracy or monotone primitive
profiles. The first FE bracket and exact-evolved-start control identify the
classical receding-flow overheating mechanism in the bulk finite-volume
update. The complete data and provenance are in
`docs/ibm_cartesian_euler_regressions_20260719.md`.

The formal Stage-2 decision is **CONDITIONAL GO**. The frozen method may enter
Stage 3 characteristic-boundary certification. Double-compression
interaction-location convergence remains required before SRP, and the bulk
near-vacuum limitation remains documented, but neither result authorizes a
change to IBM, BIC, ghost-cell layers, WENO weights, or limiter thresholds.

## Frozen downstream sequence

After phase 1, the only approved sequence is:

1. characteristic inlet, outlet, and far-field boundary conditions;
2. globally smooth curved-wall MMS and smooth/shocked Cartesian nozzle gates;
3. nozzle lip, shoulder, curvature-break, owning-feature, connected-fluid-side,
   and two-sided lip semantics;
4. independent R--Z metric-conservative discretization, axis parity/source
   balance, and R--Z positivity certification;
5. NASA CDV three-condition nozzle validation;
6. nozzle-only underexpanded jet and jet-off sphere-cone isolation cases;
7. NASA Test 1853 R--Z mean-flow and axisymmetric-mode statistical validation;
8. AMR/reflux/regrid positivity and conservation certification if production
   SRP uses AMR.

NASA CDV and NASA Test 1853 are distinct validation databases. Quasi-1-D
area--Mach relations are physics-consistency references for a finite-angle
2-D nozzle, not exact 2-D Euler solutions. R--Z results certify only the
axisymmetric mean and azimuthal mode \(m=0\).

## Required downstream contracts

### Characteristic boundary conditions

Characteristic inlet, back-pressure outlet, and non-reflecting far-field
conditions are an independent acceptance gate. Supersonic inflow prescribes
the complete incoming state. Subsonic inlet and outlet boundaries prescribe
only the incoming characteristic information, including the required total
pressure, total temperature, or back pressure. A zero-normal-gradient
boundary must not be used to tune an IBM or nozzle result when an acoustic
characteristic enters the domain.

### Nozzle reference hierarchy

Nozzle comparisons have three distinct purposes and must not be conflated:

1. a globally smooth curved-wall MMS establishes formal spatial order;
2. quasi-one-dimensional area--Mach, mass-flow, total-enthalpy, entropy, and
   Rankine--Hugoniot relations establish physics consistency;
3. a grid-converged, matched body-fitted calculation establishes the reference
   for the actual two-dimensional wall pressure and internal-shock location.

The quasi-one-dimensional solution is not an exact solution of a finite-angle
two-dimensional Euler nozzle and cannot be used alone to tune the IBM.

### Feature and nozzle-lip semantics

Expanded BIC support is permitted only within the owning smooth geometric
patch and the same connected fluid side. A sharp lip is two-sided: its two
fluid sectors have distinct boundary intercepts, normals, and admissible
support sets. Stencils must not cross a lip, shoulder, curvature break, thin
wall, or disconnected fluid region. At a curvature discontinuity the method
must use a documented lower-order admissible fallback rather than blend two
smooth quadratic reconstructions.

### R--Z metric, source, and positivity contract

Cartesian positivity evidence does not apply to R--Z. The axisymmetric Euler
update is defined in metric-conservative form,

\[
\partial_t(rU)+\partial_z(rF_z)+\partial_r(rF_r)
=(0,0,p,0)^T,
\]

with physical annular control-volume measures and radial face areas. The axis
face has zero area; no `1/r` expression is evaluated at the axis. If stored
states are cylindrical physical volume averages, radial WENO reconstruction
uses the corresponding cylindrical volume moments rather than Cartesian
moments.

For every SSPRK Forward-Euler bracket, the complete low-order R--Z update is

\[
U_i^{L,RZ}=U_i^{(s)}+\Delta t_{FE}
\left[L_i^L(U^{(s)})+S_i^{L,RZ}(U^{(s)})\right].
\]

It must be admissible before antidiffusive flux corrections are applied. If
the high- and low-order source discretizations differ, then
\(\Delta t_{FE}(S^{H,RZ}-S^{L,RZ})\) is an additional correction and must be
limited or handled by a separately verified positivity-preserving source
step. No unbounded source or axis correction may be applied after limiting.

R--Z has two separate gates:

1. density and internal-energy admissibility after every actual FE bracket;
2. metric and axis compatibility, including uniform-pressure equilibrium,
   parity, radial pressure-flux/source balance, and integral \(2\pi r\)
   conservation.

Passing either gate does not imply passing the other. A three-dimensional
spot check is required before an R--Z result is used to infer axisymmetric
three-dimensional physics.

### NASA validation cases and SRP loads

NASA CDV is the inviscid axisymmetric converging--diverging nozzle validation
with \(p_e/p_t=0.89,0.75,0.16\) and quasi-one-dimensional steady reference
data. NASA Test 1853 is the 5-inch, 70-degree sphere-cone with a 4:1 area-ratio
nozzle at \(M_\infty=2.4,3.5,4.6\), with surface-pressure, schlieren, and
internal-state measurements. They are separate databases and retain these
names in all manifests and reports.

Test 1853 requires at least three grids, two time steps, and two independent
statistical windows. Reports include mean and RMS pressure coefficients,
axisymmetric-mode spectra and peak frequencies, bow-shock and Mach-disk
locations, and axial force. R--Z spectra certify only \(m=0\).

The SRP force audit reports, with one sign and normal convention,

\[
F_{conservative},\qquad
\int_{wall}p n_z\,dA,\qquad
F_{nozzle},\qquad
R_{CV,momentum}.
\]

The final limited crossing-face flux must feed both the fluid update and the
conservative IBM force. The previously observed discrepancy between
reconstructed wall-pressure force and conservative interface force is not
hidden by selecting one diagnostic. If AMR is used, final limited fluxes with
their actual RK weights enter the flux register; admissibility and conservation
are re-audited after reflux, regrid, and coarse--fine correction.

### Production claim boundaries

Passing Cartesian phases 1--4 permits only the claim: fixed, smooth-wall,
two-dimensional Cartesian, single-species ideal-gas Euler-slip nozzle solver
for the certified flux, reconstruction, time integrator, and boundary
conditions. Passing the feature, R--Z, NASA, load, and any required AMR gates
permits only the additional claim: axisymmetric, centreline single-nozzle,
zero-angle-of-attack Euler SRP solver. Neither claim certifies NS heat transfer,
skin friction, moving bodies, arbitrary three-dimensional STL geometry, or
non-axisymmetric SRP modes.
