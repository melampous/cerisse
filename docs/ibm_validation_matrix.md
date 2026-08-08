# Fixed-Body IBM Verification and Validation Matrix

Status date: 2026-07-22.  Scope: single-component perfect gas, fixed bodies,
Cartesian and AMR meshes.  FSI and moving rigid bodies are outside this
contract.

No one benchmark validates the complete IBM.  A manufactured solution proves
order against a known field but not physical fidelity; a cylinder benchmark
tests physical loads but cannot isolate every term; a shock case does not test
wall heat flux.  Production acceptance therefore requires the matrix below.

## Production Configuration

There is no blanket production certificate for all IBM operators. The only
current production candidate is the pure sharp-interface ghost-point method
defined in `docs/IBM_CURRENT_STATUS.md`: shared-GP cell-average-aware BI-CWLS,
curvature-compatible pressure and one-sided entropy extension, full Cartesian
cell/face finite-volume update, LLF-WENO-Z5, one configured solid-side ghost
layer, conservative full-cell positivity limiting, SSPRK(4,3), stage-local
characteristic boundaries, and surface-e2 BI pressure recovery. It is
restricted to two-dimensional Cartesian, uniform level-0, fixed-wall,
ideal-gas Euler-slip flow. Face-local states, direct crossing closures,
shifted shells, TENO, AFD, cut-control, aggregate, and EB evolution are not
part of this candidate.

Within that scope, one stage-local Giles/GC-NSCBC type-2/type-3 outflow is also
verified when the tangential direction is periodic and no open cut subface
touches a periodic boundary. This does not certify a general external
farfield, multiple open boundaries or physical boundary corners.

Cut-control, aggregate, physical-subface, and EB rows retained later in this
document are historical research evidence only. They must not be used to
certify the pure GP-IBM candidate or enter its production executable.

The older two-stage `image-point WLS -> boundary condition -> GP
extrapolation` path with `interp_order=2` is explicitly unsupported for strong
shock-wall cases: it fails to establish the Mach-4 circle bow shock.  Smooth
MMS results obtained with that path remain component/formal-order evidence,
not production evidence.

For smooth stationary no-slip NS walls, the direct thermal/viscous operators
and `navier_stokes_noslip_viscous_compatibility` retain their separate MMS
status below.  They cannot be combined into a production certificate until
the complete coupled operator passes the circle hard gate and its own NS wall
heat-flux/shear gates.

The identical high-accuracy geometry file must be used on every fluid grid.
The solver is not allowed to refine or alter that geometry during a formal
fluid-grid sequence.

The 2026-07-18 curved-wall campaigns add an important qualification. A
homogeneous normal-pressure extension is not high-order compatible with a
smooth stationary Euler-slip wall carrying nonzero tangential velocity. The
numerical two-dimensional curvature closure evaluates
`dp/dn=rho_B*u_t^T*B*u_t` from the current RK-stage fluid trace. A smooth-wall
geometry callback supplies one consistent BI, local frame and shape operator,
while the polygon remains the topology carrier. The current thermodynamic
extension reconstructs `sigma=log(p)-gamma*log(rho)` and recovers density and
temperature from pressure, entropy and the ideal-gas EOS; it does not impose
an isothermal or adiabatic viscous condition on Euler.

The latest nonconstant-entropy/temperature ellipse MMS gives a BI
surface-pressure phase/orientation envelope order of 1.957 and fixed-time
first-layer pressure orders 2.44--2.75. Its original primitive-point BIC path
had only 1.67--1.85 order in first-layer semi-discrete momentum residuals.
The 2026-07-19 cell-average-aware normal-momentum BIC plus one-sided fluid
entropy jet raises the complete first-layer conservative L2 residual fits to
1.905--2.020 on both homogeneous- and nonzero-normal-entropy analytic-nozzle
MMS fields. These options are enabled in the frozen pure-GP candidate and
have now passed the current-source Gate 0 strong-shock/canonical-regression
matrix. Their multi-grid smooth-nozzle qualification is now closed by the
current-source exact inverse-MOC Gate 2 sequence. This is not a feature-edge,
narrow-gap, AMR, general 3-D STL, NS, or general strong-shock certificate.

## Mandatory Circle Hard Gate

Every change that can affect IBM geometry, ghost states, near-wall inviscid or
diffusive fluxes, reconstruction, AMR transfer, or wall loads must first pass
the frozen Mach-4 `N=320` circle regression.  Run

```text
python3 tools/check_ibm_circle_gate.py <circle-history-directory>
```

The gate requires a detached stagnation-line bow shock, `t*=t U_inf/D >= 40`,
a stationary final stand-off history, and density/pressure shock locations
within 0.5 cell. `Delta/R` must be compared with the frozen reference for the
same pressure closure and support policy. The current expanded-support,
numerical-curvature, positivity-protected pure-GP reference is
`0.549684117072`. Older compact-BIC and pre-limiter values remain historical
diagnostics and are not interchangeable baselines.
A final density/pressure/schlieren image and an active-fluid conservative-state
audit are mandatory. The unrestricted high-order update has produced negative
internal energy in a first-layer active-fluid cell and is not admissible. The
opt-in source-free Cartesian limiter must therefore be enabled for a strong-
shock time-integration certificate. On gate failure, all unrelated IBM work
stops until the circle regression is repaired.

This is a software-regression gate, not the final physical reference.  Billig
and Hornung values are diagnostics; physical certification requires a matched
body-fitted or independently converged Euler solution.

## Coverage Matrix

| Contract | Required case | Primary observables | Acceptance |
|---|---|---|---|
| Fluid solver without IBM | periodic no-solid Euler/NS MMS | RHS residual and fixed-time `rho,u,v,w,p,T` norms | expected scheme order after temporal error is negligible |
| Ghost/image reconstruction | smooth fixed circle/sphere MMS | GP values and one-, two-, three-cell fluid shells | errors decrease; finest shell is at least second-order asymptotic within sampling uncertainty |
| Surface pressure | prescribed-gradient, Euler-curvature, and NS-viscous-compatible MMS | pointwise pressure L2/Linf and integrated pressure force | approximately second order on a fixed surface; no silent fallback |
| Viscous wall traction | nonzero-shear NS MMS | `tau_normal`, `tau_mag`, viscous force | approximately second order; correct sign and tensor projection |
| Wall heat flux | nonzero thermal-Dirichlet MMS and Couette heat flux | `dTdn`, `qwall` L2/Linf | approximately second order and correct conductivity/sign |
| Low-speed physical NS | Schaefer-Turek DFG 2D-1/2 | `Cd`, `Cl`, `St`, pressure difference, recirculation | joint zero-grid/zero-Mach limit inside primary-source ranges |
| Supersonic smooth pressure | Mach-4 circle and 3-D sphere | full `Cp`, stagnation `Cp`, bow-shock position, drag, symmetry lift | monotone grid study plus reference/theory agreement |
| Corners and shocks | Mach-2 double wedge | facewise `Cp`, shock angles/locations, corner fallback | smooth-face convergence; corners treated as integrated, not pointwise, targets |
| Generic geometry | one fixed audited complex STL | geometry hash, pressure map, loads, limiter/fallback fractions | invariant STL, stable normals/ownership, converged loads |
| Pure-GP conservation and load uncertainty | source-free uniform physical cases | active-Cartesian-domain budget; BI pressure force; Cartesian crossing exchange reported separately | active-domain budget closes to reduction tolerance; BI/crossing discrepancy and reconstructed impermeability decrease under refinement and remain below the application uncertainty budget; finite-grid equality is not required |
| Implementation parity | CPU, CUDA, MPI, MPI-CUDA where supported | budgets and surface arrays | same physical result to floating-point reduction tolerance |

## DFG Physical Benchmark

`IBM/cases/validation/canonical/dfg_cylinder` implements the original
Schaefer-Turek channel and cylinder geometry as a compressible low-Mach ideal
gas.  The [primary benchmark paper](https://www.mathematik.tu-dortmund.de/lsiii/cms/papers/SchaeferTurek1996.pdf)
defines the geometry, force normalization, and accepted intervals.

Strict 2D-2 targets are

```text
Cd_max = 3.22 ... 3.24
Cl_max = 0.99 ... 1.01
St     = 0.295 ... 0.305
Delta_p at half period = 2.46 ... 2.50
```

The finite-Mach density solver is not the incompressible benchmark.  Run the
complete Cartesian matrix `Ma={0.20,0.10,0.05}` and
`N_D={20,40,60,80}`.  The primary extrapolation is

```text
q(Ma,h) = q_inf + a Ma^2 + b_Ma (h/D)^p,
```

with an independent grid-error coefficient at every Mach.  This prevents the
LLF acoustic dissipation coefficient, which can grow as Mach decreases, from
being incorrectly forced into an analytic `Ma^2 h^p` coefficient.  Report the
fixed `p=2` limit, a free common observed order, the result after dropping the
highest Mach, and the result after dropping the coarsest grid.  A raw `Ma=0.1`
result is diagnostic only.

At least six settled shedding cycles are required.  Cycle variation, surface
quality, and the full momentum budget are hard audits.  Repeat the finest
`Ma=0.1` run with first-order extrapolation at the downstream boundary; outlet
sensitivity must be smaller than the grid-extrapolation uncertainty.
Also repeat that finest run at `CFL=0.35,0.175,0.0875`.  The original benchmark
requires three time discretizations on the finest mesh.  Demonstrate the RK3
time trend or show that the two finer results differ by less than half a
reference-interval width before attributing the remaining error to IBM space
discretization.

## Current Evidence

The following results are evidence, not a blanket production certificate.

> **2026-07-17 provenance correction.** Results obtained with direct crossing
> modes, face-local states, shifted stencils, or conservative-interface
> closures certify only those explicitly selected experimental operators. They
> do not certify the standard shared-GP plus masked LLF-WENO/TENO operator.
> The production candidate is currently restricted to shared-GP
> LLF-WENO-Z5. LLF-TENO5 is unsupported for IBM shock-wall interaction until
> it passes an independent circle/corner and smooth-MMS campaign.

| Test | Result | Interpretation |
|---|---|---|
| Current pure-GP Gate 0, source IDs `5b45451c364d120e` / `ad0e0c9c747a119b` / `84db3efb00517a42` | Mach-4 circle reaches `t*=55.56` with stable detached shock and conservative early-transient limiting; current-source PM angle/pressure errors `+0.000858 deg/+0.006%`; compression shock-angle/pressure errors `-0.051865 deg/-0.003%`; inverse-MOC two-phase BI force errors `3.64e-5/1.42e-5` with relative phase spread `2.22e-5` | **Gate 0 PASS** for current-source pure-GP software non-regression; this is not absolute circle accuracy, true-curved-domain conservation, or smooth-nozzle grid-convergence certification; see `docs/ibm_pure_gp_gate0_20260722.md` |
| Current pure-GP Gate 1, source ID `74d8125a74751cbc`, source-free isentropic vortex with concentric circle, `N=96/128/192/288`, axis and `(0.37h,0.23h)` phases | all-fluid L2 orders 1.96--2.09; first-wall-layer L2 1.90--2.04; fixed physical near-wall-band L2 1.96--2.23; BI pressure 2.25--2.27; reconstructed `rho un` L1/L2/Linf 2.56--2.68; 984/984 FE brackets use `theta=1`; 1723/1723 GP states use quadratic BI-CWLS | **Gate 1 PASS** for approximate-second-order L1/L2 smooth curved-wall accuracy, BI surface pressure, impermeability convergence, limiter transparency, and two-phase consistency; global/axis-phase Linf fits of 1.73--1.93 remain documented grid-phase/max-location variation, so a uniform strict-second-order Linf claim is excluded; see `docs/ibm_pure_gp_gate1_20260722.md` |
| Current pure-GP Gate 2, source ID `84db3efb00517a42`, exact continuous-curvature inverse-MOC nozzle, `256x64/384x96/576x144`, two phases | 4864/4864 GP states use quadratic BI-CWLS; no limiter event; phase-envelope fluid/first-layer pressure errors decrease with adjacent orders 1.66--1.76/1.59--1.70; finest BI `Cp` L1 `2.43e-5`, nonendpoint Linf `3.13e-4`; absolute BI mass exchange/reference flow decreases `9.91e-4 -> 3.21e-4 -> 1.45e-4`; mass-flow spread reaches `1.10e-4`; entropy remains about second order; finest BI force error at most `1.01e-5` and phase spread `1.76e-6`; inlet/outlet-location pressure effect is `1.31e-8` relative L2 | **Gate 2 PASS** for restricted smooth two-dimensional Cartesian Euler-slip nozzle engineering use; this is not a second formal-order proof, local velocity orders are nonuniform, Cartesian crossing exchange is not physical BI traction, and shocked nozzle/AMR/R-Z/NS/3-D remain uncertified; see `docs/ibm_pure_gp_gate2_20260722.md` |
| Current pure-GP Gate 3, source ID `c92f3840c8626937`, smooth divergent shocked nozzle, `384x64/576x96/768x128`, two phases, plus fine-grid `x_out=12/14` comparison | both fine phases pass the 3% shock-state/RH gate; RH mass/momentum residuals are about 2.6%/1.5%; shock phase spread is `0.089h`; off-shock mass-flow spread is at most `0.055%`; absolute BI leakage is at most `0.055%`; fine BI-force phase spread is `0.01435%`; six-run matrix has 407380/407380 high-order FE brackets and 6912/6912 quadratic GP states; outlet move changes shock mean by `7.89e-5`, BI force by `2.16e-5`, and common pressure by `4.9609e-4` relative L2 | **Gate 3 CONDITIONAL PASS** for the demonstrated fine-resolution stationary internal-shock nozzle; the original pilot retains `u/v` mirror-Linf `5.98e-8/3.38e-8` versus the added `1e-8` criterion, coarse/intermediate local shock errors exceed 3%, and outlet-pressure sensitivity passes with little margin; see `docs/ibm_pure_gp_gate3_20260722.md` |
| Analytic-nozzle cell-average A/B/C operator audit with cell-average normal momentum and one-sided entropy jet, N=96/144/216/324 | nonzero-normal-entropy first-layer C L2 orders: `rhou=2.007`, `rhov=1.905`, `rhoE=2.018`, `rho=2.020`; second layer 1.979--2.022; homogeneous-entropy sequence gives first layer 1.905--2.020; every GP uses the joint P2 path | restricted smooth 2-D Euler-slip semi-discrete second-order gate passes; final-pair transverse momentum remains grid-phase sensitive, and this is not fixed-time or shock certification |
| Same candidate, disabled-branch and CPU/CUDA N=96 audits | disabled cell-operator CSV is bitwise identical to the frozen baseline; CPU/CUDA GP primitive states differ by at most `1.49e-14`, flux divergence by `1.32e-13` | default-off behavior and production ghost/flux execution-space parity pass; the generated symbolic MMS source has a separate `9.04e-4` CUDA cancellation discrepancy and is excluded from this parity claim |
| Nonconstant-entropy/temperature variable-curvature ellipse MMS, no solid, N=64/96/144/216 | fixed-time all-fluid primitive orders 3.09--3.75; t=0 conservative residual orders 3.10--3.66 | validates the globally smooth source, exact outer boundary and bulk LLF-WENO-Z5 protocol independently of IBM |
| Same field with analytic smooth-wall BI/frame/shape operator, axis/shifted/rotated sequences | BI surface-pressure L2 orders 1.892/1.973/1.928; phase/orientation envelope 1.957; fixed-time first-layer pressure 2.44--2.75 and velocity 1.84--2.08 | curvature-pressure plus entropy-based complete-state extension passes the restricted pressure gate; first-layer velocity remains conditional |
| Same ellipse, t=0 semi-discrete residual | all-fluid orders 2.14--2.51; first-layer density/energy 1.64--1.98 and momentum 1.67--1.85 | strict local second-order gate fails; grid phase is influential but is not the sole source |
| Same ellipse, fresh CPU/CUDA sm_89 N=64 comparison | relative L2 difference 3.5e-16--7.2e-16; maximum absolute difference 7.99e-15 | single-level execution-space equivalence passes for the latest entropy path |
| Assertion-enabled forced complete-state fallback, ellipse N=64 | 68/68 GPs select fallback code 1 and all joint-state assertions pass | pressure, temperature and velocity reduce together; branch coverage only, not fallback accuracy evidence |
| Source-free stationary isentropic vortex, no solid, N=80/120/160/240/320 | fixed-time global primitive rates 3.66--4.24 | high-order smooth bulk verification under this RK/boundary protocol; not a standalone fifth-order claim |
| Source-free stationary isentropic vortex with concentric circular IBM and production limiter, N=80/120/160/240/320 | fixed-time all-fluid L2 orders `rho/u/v/p/T=2.09/2.09/2.09/2.02/2.01`; first-layer L2 orders `2.02/2.07/2.07/2.02/1.98`; all 648 FE brackets use `theta=1` | PASS: the limiter is inactive and does not degrade the restricted smooth curved-wall second-order solution |
| Mach-2 PM expansion and single compression corners with production limiter, current pure-GP source | PM flow angle/Mach/pressure errors `+0.000858 deg/<0.001%/+0.006%`; compression shock-angle/flow-angle/pressure errors `-0.051865 deg/+0.003461 deg/-0.003%`; no limiter event | PASS for turning relations and downstream thermodynamic states; no PM plateau cell is more than 1% below theory |
| Mach-2 sequential `0 -> 10 -> 20 deg` double compression with production limiter | state-1/state-2 pressure errors `-0.060%/-0.064%`; flow-angle errors about `-0.013 deg`; first/second fitted shock-angle errors `+1.20/+0.148 deg`; interaction-location error about `4.0h`; no limiter event | CONDITIONAL PASS: sequential states and topology pass, but the short first ray and interaction location are not grid-converged |
| No-solid symmetric near-vacuum double rarefaction, t=0 discontinuous startup, 512x64 and 1024x128 | all 12,668 FE brackets remain admissible with `theta=1`; exact transverse uniformity and symmetry; exact FV star-region relative L1 errors on the 512 grid are rho/p/T/sigma = 22.0/46.8/90.4/77.9%, and the anomalous physical width does not collapse on the 1024 grid | admissibility/conservation PASS, limiter activation N/A, ab-initio star-region thermodynamic accuracy FAIL; this is a bulk receding-flow overheating defect, not an IBM wall failure |
| Same exact solution initialized from exact conservative averages at physical age `tau=0.01`, 512x64 | reaches the same physical age 0.15 with all 3944 FE brackets at `theta=1`; star-region relative L1 errors in rho/p/T/sigma = 0.297/0.359/0.459/0.566%; first-bracket sigma error changes from +0.08125 for t=0 startup to -1.58e-8 | causal diagnostic PASS: removing the initial discontinuity removes almost all persistent overheating without changing LLF-WENO-Z5, SSPRK43 or the limiter |
| Same exact solution with a concentric circular IBM and homogeneous `dp/dn=0` | global/first-layer pressure rates 0.56/0.51; first-layer energy residual remains `O(1)` | `pzero` is not high-order compatible for this curved Euler-slip state with nonzero tangential velocity |
| Earlier circle oracle with analytic compatible `dp/dn=rho*u_t^2/R` | ordinary path: global pressure 2.14, first-layer pressure 1.88, surface pressure 2.05, first-layer velocity 1.29 | historical mechanism-isolation result; the numerical-trace entropy path below supersedes it as the current restricted candidate |
| Same earlier oracle, replacing only degree-one constrained GP states by exact values | first-layer velocity/pressure orders 2.63/2.73; momentum/energy RHS orders 1.99/1.95; `N=320` velocity and momentum-RHS errors reduced about 100/377 times | causally identifies reduced ghost-state reconstruction in that test; validation-only and not a production no-fallback method |
| Same earlier oracle with expanded visible full-rank BIC support | all five grids use degree-two reconstruction; first-layer velocity/pressure orders 2.46/2.62; first-layer conserved RHS orders 1.94--2.01; shifted phase and CPU/GPU short gates pass | established the support-expansion mechanism later exercised by the nonconstant-entropy ellipse; still restricted to two-dimensional smooth walls |
| Same exact solution with numerical `rho_B,u_t,B` curvature compatibility and complete-state admissibility | first-layer t=0 orders: density 2.02, momentum 1.89, energy 2.02; fixed-time first-layer orders: velocity 2.42, pressure 2.65; 88/88 GPs use the high-order path; CPU/CUDA maximum field difference `7.99e-15` | the 2-D smooth fixed-wall Euler-slip pressure closure passes its exact-solution and implementation-parity gates |
| Earlier validation-only forced complete-state reduction, constant-entropy circle `N=80` | 88/88 GPs select fallback code 1; assertions verify one support supplies `p,T,u_t` and its normal velocity is reflected | historical temperature-based branch coverage; the current entropy-based branch is covered separately by the `N=64` ellipse test above |
| No-solid 3-D WENO thermal NS MMS, N=32/48/64/96 | L2 pair orders are about 4.2-5.0 for primitive fields | bulk solver retains high order in this smooth case |
| Standard shared-GP LLF-WENO-Z5 2-D full-NS circle MMS, N=64/96/128/192, all experimental crossing/shell modes off | fixed-time all-fluid L1/L2/Linf orders are approximately 2.00-2.15; first-fluid-cell-layer L2 orders are 2.09 for `u`, 2.22 for `v`, 2.57 for `rho/T`, and 2.86 for `p` | this is the first full-chain second-order solution evidence for the standard shared-GP operator; it applies only to this smooth LLF-WENO-Z5 case |
| Same standard shared-GP case, t=0 semi-discrete residual | first-layer L2 orders are 1.83 `rho`, 1.92 `rhou`, 1.99 `rhov`, and 1.83 `rhoE`; bulk L2 is about 3.9-4.2 | the IBM layer limits the operator and remains grid-phase sensitive; do not claim uniformly second-order local truncation error in every norm |
| Independent production/validation builds of the standard shared-GP case, N=64 | the t=0 RHS plotfile and the 400-step `t=0.002` plotfile are recursively byte-identical | `CNS_IBM_VALIDATION` instrumentation is numerically inert for this campaign; the reported convergence belongs to the production operator |
| Legacy shared-GP LLF-WENO-Z5 Mach-4 circle, N=320, all experimental crossing/shell modes off | `Delta/R=0.5159359` at `t=0.002`, but the settled value at `t=0.008` is `0.5712879` | the archived fixed-time result is reproducible, but its apparent agreement with Billig was transient and is not a physical certificate |
| BI-centered constrained shared-GP LLF-WENO-Z5 Mach-4 circle, N=320 | clean detached shock; settled `Delta/R=0.5643135`; density/pressure location span `0.1427` cell; steady legacy difference `-1.22%` | passes the frozen software circle gate; matched body-fitted physical accuracy remains open |
| Fixed-geometry `D/h=64` compact versus expanded-support BIC circle, `t*=55.56` | final `Delta/R=0.5618378/0.5605581`; expanded-support difference `-0.228%`; last-four-sample range `2.78e-6`; no negative state | expanded support passes the strong-shock non-regression gate; this does not measure absolute physical accuracy or validate curvature pressure compatibility |
| Fixed-geometry `D/h=64` expanded-support numerical-curvature BIC circle, `t*=55.56` | stable symmetric bow shock; `Delta/R=0.5498128788`; stand-off 17.59 cells; density/pressure estimator span 0.101 cell | spatial strong-shock non-regression passes, but time-integration production certification does not; 12--30 of 180 GPs use local complete-state reduction during the 100--500-step audit; Billig/Hornung differences are not exact error estimates |
| Strict RK-stage admissibility audit for the same circle | curvature: step 132, SSPRK43 stage 4, first-layer active cell `(187,143)` has `rho*e=-2.33918`, `p=-0.935671`; matched `pzero` fails at stage 1 near the same physical time; CFL 0.15 still gives two invalid first-layer active cells | actual cell-average admissibility failure in the coupled near-wall shared-GP/LLF-WENO-Z5/SSPRK update; stage 4 is the accepted `U^(n+1)`, and stage 1 feeds the next RHS; `pzero` excludes the new curvature closure as the unique cause but does not prove a no-solid bulk origin |
| Frozen step-151 all-low-order and flux-limiting audit, source `b7102d469cdae3b9` | piecewise-constant shared-face Rusanov update is admissible at multidimensional CFL 0.22995; global-theta oracle passes through step 200; validation-only local convex limiter passes 65/65 triggered stages while limiting 2--36 of 199112 faces per event; maximum active-fluid conservation residual is `2.44e-15` relative | the high/low flux pair can restore fully discrete admissibility and a conservative local shared-face construction is feasible; production high-order update still fails without the oracle, and AMR/reflux/long-time accuracy are not certified |
| Production-path shared-face limiter and conservative IBM load, source `c68c0fa9bdd7c4e8`, Mach-4 circle N=320 to `t*=55.56` | 15,436 steps; 80 early local-limiter events; zero global/all-low fallback; minimum theta 0.25993; maximum relative conservation residual `1.138e-15`; final `Delta/R=0.5496841171`, four-sample span `1.056e-7`; final wall/conservative `Cd=1.24615/1.25162` | PASS for opt-in, uniform level-0, source-free Cartesian Euler; finalized limited crossing flux feeds both the fluid update and conservative IBM load; automatic time-step retry, characteristic BC, AMR, R--Z and nozzle gates remain open |
| Current source `75f7ee86428df61f`, same N=320 frozen circle protocol to `t*=55.56` | AMReX `fcompare` is exactly zero for all 11 plot variables at steps 1000/6000/8000/10000/15436; final VTP SHA-256 is identical; 61,744 bracket modes, limiter counts, minimum theta, conservation residual and `Delta/R=0.5496841171` all reproduce | PASS: the default-off first-FE audit added for near-vacuum diagnosis is numerically inert and the current-source circle provenance is closed |
| Stage-3A characteristic inlet/outlet matrix | oblique uniform-flow RHS below `3.5e-13`; total-condition duct mass-flow error `-1.85e-6`; supersonic outlet is bitwise independent of back pressure; normal acoustic incoming/outgoing ratio `O(1e-5)` | PASS for straight-buffer internal-flow inlet/outlet use; this is not a general external far-field certificate |
| Stage-3B Giles/GC-NSCBC local vorticity oracle, N=48/96/192/384 | pressure-conversion error converges at about 2.06--2.25; CPU/CUDA difference `2.44e-15`; one/two-rank decomposition is bitwise identical | operator-level PASS; physical corners, multiple open boundaries and boundary-time Fourier reflection remain system-level gates |
| Stage-4A C-infinity all-subsonic nozzle, four grids/two phases plus CFL/initialization/buffer/acoustic matrix | all runs admissible with no limiter activation or GP linear fallback; strict finalized-flux SSPRK mass/momentum/energy residual is at most `3.93e-12`; the axial fitted frequency follows domain acoustic travel time and is CFL-insensitive; high-cadence probes identify a symmetric `f=0.819--0.895` transverse mode | integration/grid-phase and discrete-budget PASS; steady quantitative protocol FAIL because a transverse acoustic field mode leaves the final-two-time-unit pressure change at `8.65e-4` even after the axial mode decays; no IBM change is authorized |
| Stage-4B smooth-wall recompression/weak-shock nozzle, N64/128 to `t=40` | no limiter/fallback; max adjacent pressure jump `0.0871/0.0719`, `h|grad p|/p=0.0776/0.0691`, and pressure second-difference sensor `0.0109/0.00811`; downstream X-wave sharpens under refinement | robustness PASS; shock-free fixture designation INVALID, not a solver failure; retain as a weak-shock regression without changing IBM, WENO or limiter parameters |
| Historical inverse-MOC full-cell Cartesian exchange audit, volume-e1/surface-e2, four grids and two phases | BI pressure recovery and reconstructed wall-normal velocity pass, while Cartesian production/CV exchange differs from the BI physical load by 1.38--12.26% and signed Cartesian IBM mass exchange is non-monotone | retained as evidence of the known stair-step/full-cell true-domain conservation limitation; it is not a failure of the current pure-GP BI pressure/force definition and is superseded as a Gate-2 acceptance test by the current-source row above |
| Same inverse-MOC state, exact solid-side continuation through the frozen Cartesian crossing operator | finest-grid exact-state load phase spread 7.803%; maximum production-BIC increment 0.133% and masking increment 0.0084% | causal localization PASS: surface recovery, BIC and WENO masking are not the dominant load error; the missing term is true curved-wall control-volume geometry |
| Validation-only cut-control geometry and true-domain residual oracles | faceted polygons close geometrically, but fixed-facet MOC tangency is `7.19e-5` and fails as an exact oracle; one Hermite provider gives geometric closure `4.45e-15`, tangency `1.34e-8`, and at least 16.7x lower energy residual | common smooth geometry semantics PASS; production conservative residual and mandatory small-cell stabilization are not implemented, so the solver remains NO-GO for phase-robust nonzero IBM loads |
| Validation-only inverse-MOC active-cell-rooted agglomeration, four grids/two phases | all 4634 center-in-solid physical fluid fragments are assigned once to an existing center-in-fluid active root; no orphan/cycle; maximum 2 members and 1 parent edge; minimum aggregate volume `0.5004328 Vcell`; maximum aggregate closure `2.14e-13`; fluid-volume error `1.45e-14` | geometry/topology PASS for this 2-D smooth channel; production state semantics, redistribution, positivity, MPI/FAB ownership and finalized conservative RHS remain unimplemented |
| Validation-only aggregate/Cartesian cell-average transfer, four grids/two phases | bidirectional P2 reproduction error at most `5.88e-15`; maximum condition number `23.24`; maximum weight L1 norm `2.78`; radius at most 2 cells; direct identification of the two averages has an approximately first-order near-wall defect | formula-level state-semantics PASS; production storage, nonlinear admissibility, MPI/FAB ownership and time advancement remain unimplemented |
| Validation-only aggregate pressure residual, four grids/two phases | one shared pressure flux per open face segment plus physical-wall pressure traction; two-point quadrature RHS L2 at most `2.33e-12`; reference Linf at most `7.86e-11`; shared interval mismatch `2.22e-16`; one-point wall evaluation retains an approximately first-order residual | stationary Euler momentum operator PASS at formula level; nonlinear mass/energy flux, SSPRK advancement and production integration remain open |
| Validation-only full-Euler aggregate residual, straight oblique channel, four grids/two phases | physical aggregate averages plus P2 reconstruction, shared two-point open-face Rusanov flux and true-wall pressure traction; minimum finest-pair component order `2.0447`; exact-reference L2 at most `2.51e-12`; zero inadmissible states | smooth nonlinear Euler operator PASS at formula level; production RHS/time advancement, curved-wall shocks and distributed ownership remain open |
| Production cut-control inverse-MOC matrix with one stage-local BIC wall-pressure trace, Ny=32/48/64, two phases | force-error envelopes `0.07027/0.04972/0.01255%`; finest two-phase spread `0.00025%`; physical wall mass flux exactly zero; maximum semi-discrete/RK/true-domain relative closure `3.94e-16/1.73e-14/1.73e-14`; one-/two-rank, multi-FAB and CUDA `sm_89` checks pass | PASS for the two-dimensional, uniform level-0, fixed smooth Euler-slip hybrid GP--cut-control operator; local pressure accuracy, AMR, NS, R--Z, moving walls and general STL are not implied |
| Aggregate physical-subface limiter smooth transparency, inverse-MOC Ny=32/48/64, two phases | 2328/2328 SSPRK43 FE brackets remain high order with `theta=1`; limiter-on/off budgets and wall forces differ only by `O(1e-14)`; forced `theta=0.37` branch passes CPU, two-rank/four-FAB and CUDA | implementation and smooth transparency PASS; the forced coefficient is branch coverage rather than a physical activation test |
| Exact inclined-wall reflected shock, natural limiter activation, N=64/96/128, two phases | natural local events occur on every activation grid; minimum `theta=0.2277--0.2355`; limited-entity fraction decreases from about 3.66% to 1.74% and remains in the `O(h)` shock band; no clipping/global/all-low fallback; true-domain closure at most `1.27e-14`; shock location is within about `0.22--0.27h` | PASS for naturally activated spatially nonuniform aggregate limiting and strong-shock conservation in the restricted level-0 scope |
| Transactional all-low rejection/retry, same inclined-wall case | original full interval reaches an inadmissible all-low FE bracket, rolls back, and is recomputed as two half-step SSPRK43 advances; CPU one-FAB and two-rank/multi-FAB conserved states are bitwise identical to the a-priori half-step reference; CUDA maximum conserved-state difference `1.14e-13`; wall-pressure VTP/force and additive budgets agree within `4.69e-13`; exhausted retry terminates without clipping | PASS for transactional state/time/physical-flux/wall-load/budget rollback in the same restricted configuration; this is not general adaptive stepping, NSCBC, AMR, R--Z, NS or moving-wall certification |
| Stage-local Giles/GC-NSCBC plus cut-control limiting and retry, N=64 | 16 identity-path FE brackets retain `theta=1`; forced local `theta=0.37` closes conservation within `1.35e-15`; CFL sweeps exercise all-low, BI-CWLS wall-pressure and high-order cut-residual rejection; direct versus rejected/subdivided states are bitwise identical separately on CPU and CUDA; CPU/CUDA difference is at most `6.57e-14`; true-domain closure is at most `1.24e-14` | PASS for the restricted one-outflow, periodic-tangent, level-0 joint update; persistent NSCBC, sonic/reverse-flow behavior, multiple open boundaries, corners and general farfield reflection remain uncertified |
| Conservative cut-control checkpoint/restart | checkpoint metadata records aggregate-root state semantics, cut-control manifest hash, dimension and state layout; continuous versus CPU restart final state and post-restart budget rows are bitwise identical; CPU checkpoint to CUDA differs by at most `9.51e-14` | PASS for metadata-gated level-0 restart; legacy checkpoints and mismatched geometry/state semantics fail closed; AMR/regrid restart is not implied |
| `pzero` disabled-branch regression, old source `a10cd62fda9911da` versus current `a950002d8fa5ecea` | same 20-step Mach-4 circle protocol; AMReX `fcompare` reports zero absolute and relative error for every conservative/primitive field and IBM marker | selecting no curvature closure leaves the historical `pzero` numerical path bitwise unchanged |
| Same BI-centered circle, N=192/320/512 | steady `Delta/R=0.5783158/0.5643135/0.5557441`; shock remains detached and symmetric | monotone grid trend; three-point diagnostic order is about 0.83, but no body-fitted reference or rigorous asymptotic uncertainty yet |
| Original two-stage quadratic image-point WLS, i2/e1/alpha=0.6, Mach-4 circle N=320 | no detached bow shock through `t=0.00537` (`t*=37.3`) | strong-shock failure is in the GP-BI-IP extension before WENO; this operator is unsupported for production |
| Shared-GP LLF-TENO5 Mach-4 circle, N=192/320 | the unmodified path has no resolved detached shock; an experimental WENO-Z boundary band is grid-sensitive and rejected | LLF-TENO5 is not certified for IBM shock-wall production |
| Fixed 4096-face 3-D NS-viscous-compatible MMS, N=32/48/64/96 | surface pressure 2.32, `dTdn` 1.87, `tau_mag` 1.89 in L2; corresponding Linf values about 2.46, 2.00, 2.00 | strongest direct evidence for the current fixed-geometry second-order wall-flux path |
| Earlier WENO/TENO e2 direct-flux N=24/32/48 sequences | `dTdn` and `tau_mag` about 1.43-1.49 | coarse-grid/pre-asymptotic warning; must remain visible in reports |
| Fixed Mach-4 4096-segment circle, D32/48/64/96 | D96 stagnation `Cp` 0.37% low; D96 surface/CV force mismatch 0.92% | stagnation pressure and finest integrated force accepted for that case; full long-time `Cp` convergence still pending |
| DFG inlet-compatible CPU/CUDA 20-step check | VTP pressure difference below 2.6e-13; budget closure difference below 1.2e-9 absolute | CPU/CUDA implementation parity passed |
| DFG 500-step CUDA gate, N_D=20, Ma=0.1 | finite positive surface state, all 4096 fits order 2, no pressure fallback, endpoint momentum residual 0.65% | startup/stability gate passed; not a settled physical result |
| DFG 2D-2 TENO5, N_D=20, Ma=0.1, t=8 | six-cycle `Cd_max=3.17287`, `Cl_max=0.74136`, `St=0.30719`, `Delta_p=2.39108`; mean surface/CV vector-force mismatch 1.30% | periodic/surface/geometry checks pass, but the direct 1% force-closure gate fails; finite-Mach coarse values are not incompressible acceptance values |
| DFG 2D-2 scheme compatibility restart, N_D=20, Ma=0.1, t=8 | CUDA fixed-IBM restarts ran with LLF-TENO5, LLF-WENO-Z5, AFD-HLLC-TENO5, AFD-HLLC-WENO-Z5, MUSCL-HLLC, and Central4. The conservative interface/CV mismatch stayed at `3.5e-12` to `5.1e-12` mean relative level. Surface/interface mismatch means were WENO-Z5 1.20%, TENO5 1.53%, AFD-HLLC-TENO5 2.55%, Central4 2.79%, AFD-HLLC-WENO-Z5 3.20%, MUSCL-HLLC 3.63%. | these schemes are IBM/GPU/viscous compatible in the finite-volume budget sense; WENO-Z5 is the best short-restart surface-load path in this low-Mach case; long-time Cd/Cl/St ranking still requires settled-cycle runs |
| Two-rank MPI and AMR regrid/reflux smoke tests | build/run and topology/remap audits pass | infrastructure works; physical AMR force convergence remains pending |
| IBM AMR image-support gate, 2-D CPU/MPI, CUDA `sm_89`, and 3-D CPU/MPI | automatic refinement gives zero missing nonzero supports; a deliberately narrow patch reports 176 affected GP targets and aborts before RHS | coarse-fine image support is now fail-closed; physical AMR solution/load convergence remains a separate gate |
| Direct surface, heat-only, and viscous-only NS MMS, N=64/96/128/192 | `dTdn` 2.107, `tau_normal` 2.128, and `tau_mag` 2.094 in L2 at fixed time | independent Fourier and viscous IBM operators pass the second-order gate |
| Fully coupled NS MMS, N=64/96/128/192 at t=0.002 | pressure 3.247 and `tau_normal` 2.242, but `dTdn` 1.624 and `tau_mag` 1.440 in L2 | do not claim fully coupled second-order wall heat flux or shear; wall-location-consistent direct diffusive flux remains open |

The R4 surface data use a fixed 4096-element geometry and the production
`e2/s2/alpha=1/ghost_layers=1` path.  Co-refined-geometry MMS files are useful
diagnostics but are not accepted as proof for the fixed-STL production
contract.

## Hard Reporting Rules

1. Record physical time, CFL, RK scheme, flux/reconstruction, grid spacing per
   body diameter, geometry SHA-256 and facet count, IBM orders, alpha, ghost
   layers, and pressure closure.
2. Report L1, L2, and Linf, including the first fluid-cell shell.  A good bulk
   norm cannot hide a low-order wall shell.
3. Report pointwise surface distributions and integrated loads.  Agreement of
   one scalar drag coefficient is insufficient.
4. Separate pressure and viscous force and state the IBM normal/sign
   convention.  Compare surface force with the conservative control-volume
   force.
5. Keep failed and pre-asymptotic grids in the record.  They may be excluded
   from a fit only with a stated, reproducible criterion.
6. Do not apply incompressible reference intervals to a single finite-Mach
   run, and do not use equal step counts as equal physical time.
7. A smoke test establishes implementation viability only.  It is never a
   physical or convergence validation.
8. A fixed-time circle value is not a steady shock stand-off result.  Report
   `Delta/R(t)`, `t*=t U_inf/D`, the plateau criterion, density/pressure shock
   agreement, and the final field image.
9. A failed circle hard gate blocks all downstream IBM validation and
   production work.  It cannot be waived by a good MMS order or an integrated
   force result.

## Remaining Gates

1. Restrict the smooth BI/frame/shape-operator and entropy-extension path to
   fixed smooth 2-D analytic/spline walls until owning-feature,
   connected-fluid-side, narrow-gap, AMR coarse-fine, and MPI/FAB-decomposition
   gates are explicitly authorized and passed.
2. Gate 3 conditionally closes the stationary internal-shock nozzle sequence
   at the demonstrated fine resolution. Preserve the original pilot mirror
   failure, the narrow outlet-pressure margin, and the minimum local-resolution
   requirement in every production manifest. A materially different shock
   topology, geometry, boundary layout or numerical method requires its own
   qualification rather than inheriting this result.
3. Record active-fluid admissibility at every RK stage in the nozzle gates.
   The final limited flux is integrated into the production uniform-level
   Cartesian budget and conservative IBM load paths. A shock-containing AMR
   nozzle remains blocked until that finalized flux also enters FluxRegister/
   reflux and post-reflux admissibility is audited. Per-cell clipping or
   relying on primitive/EOS floors is not acceptable.
4. After each isolated change, pass the mode-matched frozen Mach-4 circle
   software gate.
   The absolute stand-off error remains open until a separately approved,
   matched, independently converged Euler reference is available.
5. Complete the DFG 2D-1 steady and 2D-2 Mach/grid matrices and outlet
   sensitivity run.
6. Complete Mach-2 double-wedge and 3-D sphere grid sequences with fixed
   geometry and force closure.
7. Add one high-quality fixed complex STL and retain its geometry audit and
   hash in every result.
8. Repeat the accepted uniform physical cases with AMR and require both
   solution/load agreement and reflux-aware momentum closure. The image-point
   support audit is necessary but does not establish physical AMR convergence.
9. Add a physical nonzero-wall-heat-flux dataset before claiming hypersonic
   aeroheating validation; MMS and Couette verify the operator but not the
   complete physical model.
10. Raise the fully coupled fixed-time NS `dTdn` and `tau_mag` sequences to
   second order before certifying production wall heat flux and skin friction.
