# Test 1853 single-nozzle 3-D failure audit

Date: 2026-08-07; deterministic provenance audit updated 2026-08-08

## Direct answer

The failed cells did not miss a second-order WENO fallback. No second-order
fallback exists in this flux path. The production choice is either
characteristic LLF-WENO-Z5 or a two-state, first-order LLF/Rusanov face flux.
The face sensor is not a cell-admissibility test, so fallback on one or more
faces does not guarantee positive internal energy after the six-face 3-D
update.

The AMR regrid path is a separate issue. It uses AMReX limited conservative
linear interpolation. That interpolation is second order in smooth regions,
but it does not enforce the nonlinear Euler condition

```text
rho > 0 and rho*e = rhoE - |rho*u|^2/(2*rho) > 0.
```

## Frozen-run evidence

Run directory:

```text
/shared/cerisse_test1853_3d_20260807/runs/test1853_single_h100/
run165_chain_20260807T115518Z/stage2_jet_l1
```

The failed executable reported:

```text
inviscid_scheme=characteristic_llf_weno_z5
local_first_order_llf_fallback=enabled
density_jump=0.003 density_curvature=0.08
3-D positivity limiter unavailable/off
AMR state interpolation = limited conservative linear
```

The input contained `cns.stage_positivity=1`, but the exact source tree used
for the run, `/shared/cerisse_local_0807c`, has no reader or implementation for
that keyword. It was therefore a silent no-op. The actual conservative
shared-face limiter was explicitly disabled with
`cns.ibm_positivity_flux_limiter=0`.

At coarse step 3059, time `6.845347345e-4 s`, regridding changed L1 from the
previous layout to 257 grids and 7,326,208 cells. The first following L1
substep failed at `6.845905585e-4 s` with eight pure-fluid cells,
`min(rho*e)=-931.6492316`, and no NaN/Inf. Every reported bad cell had both IBM
markers equal to zero.

An exact grid-log containment audit gives:

| L1 cell | present after step-3051 regrid | present after step-3059 regrid |
| --- | --- | --- |
| `(186,120,257)` | no | yes |
| `(184,124,257)` | no | yes |
| `(188,120,259)` | no | yes |
| `(190,120,261)` | no | yes |
| `(184,227,257)` | no | yes |
| `(188,231,259)` | no | yes |
| `(188,120,92)` | no | yes |
| `(188,231,92)` | no | yes |

Thus all eight failures occurred in newly created fine cells. This proved that
the failure chain passed through coarse-to-fine interpolation. Before the H100
replay below, it did not yet distinguish these two possibilities:

1. the interpolated state is already inadmissible before flux evaluation;
2. the interpolated state is admissible, but its first WENO SSPRK forward-Euler
   bracket is inadmissible.

The old executable did not record a stage-input audit or the face fallback
mask, so either stronger causal claim would exceed the stored evidence.

## Diagnostic added

The current source has a default-off mode:

```text
cns.llf_fallback_mask_diagnostics=1
```

With the production positivity limiter disabled, this mode preserves the
historical update and ordinary reflux transaction. It performs these checks:

1. audit each SSPRK43 forward-Euler input immediately after FillPatch;
2. if the input is admissible but the high-order bracket fails, record the
   exact fallback decision on all six faces of the worst cell;
3. compare the production WENO/mixed flux with an all-face first-order
   Rusanov baseline and reconstruct `dt*div(F_high-F_low)`;
4. abort without changing state or flux.

The diagnostic changes neither GP geometry nor BI-CWLS, thermodynamic ghost
extension, surface recovery, AMR state data, or the accepted full-cell flux.
It retains `pure_shared_gp_full_cartesian`; cut-cell, cut-control, aggregate,
and EB flow mechanisms remain disabled.

## Deterministic H100 replay

The diagnostic executable was rebuilt on the H100 node with CUDA 12.6 and
native `sm_90` code. Its deployed source hashes matched the locally verified
source. All replay runs kept the pure shared-GP, full-background-cell WENO
operator, reflux enabled, and the production positivity limiter disabled.
The ignored historical input was explicitly overridden with
`cns.stage_positivity=0` so that the current fail-closed parameter validation
did not mistake it for an implemented limiter.

Four tests were run from the frozen state:

| test | regrid and time-step control | result |
| --- | --- | --- |
| A | restart `chk03000`, `amr.regrid_int=-1 -1` | Passed through coarse step 3061 with no inadmissible stage input or FE bracket. |
| B | restart `chk03000`, original `amr.regrid_int=8 8`, normal CFL step | Reproduced the step-3059 regrid and failed before the first L1 flux evaluation. |
| C | restart the test-B `chk03059`, force the same regrid, global `dt_max=1e-14` | Failed at the same pre-flux input audit with identical counts and extrema. |
| D | repeat C with retained old-fine coverage and coarse-source stencils | Classified all 12 inadmissible cells as new valid cells prolonged from L0; no bad cell came from old-fine overlap. |

The target regrid in the original run, test B, and test C had 257 boxes and
7,326,208 L1 cells. After removing only the `DistributionMapping` rank suffix
from each grid-log row, all three target BoxArrays had the same SHA-256:

```text
0f7e288d6b9450e86f8d63b580c30ab55ede9f11d4f187d5fb1ce183815822d4
```

Test B used coarse/fine time steps `1.116478628e-7/5.582393138e-8 s`.
Test C used `1e-14/5e-15 s`. Both stopped at exactly the same audit:

```text
[LLF-FALLBACK-STAGE-INPUT] scheme=ssprk43 FE_bracket=1/4 level=1
t=0.00068453473453651255 verdict=INADMISSIBLE_BEFORE_FLUX_EVALUATION
[IBM-Stage-Admissibility-Oracle] FE-input active=5278364 nonfinite=0
invalid_rho=0 invalid_rhoe=12 min_rho=0.0028075071815553462
min_rhoe=-984.86910627266479
```

This closes the earlier evidence gap. The newly prolonged fine state is
already outside the Euler admissible set before WENO reconstruction, local
LLF fallback selection, flux divergence, or any SSPRK update. Consequently,
the six-face fallback mask and WENO/Rusanov flux difference are not defined
for this failure: the diagnostic correctly aborts before evaluating either
flux.

The responsible code path is `CNS::init(AmrLevel&) -> FillPatch` with the
StateDescriptor mapper `lincc_interp`. AMReX limited conservative linear
interpolation limits reconstructed conservative components, but it does not
enforce the coupled nonlinear condition
`rhoE - |rho*u|^2/(2*rho) > 0`. Near this strong shock, the prolonged density
remains positive while the independently reconstructed energy and momentum
produce negative internal energy.

### Test D: exact regrid provenance

Test D enabled the additional default-off, read-only mode:

```text
cns.regrid_prolongation_diagnostics=1
```

Immediately before `FillPatch`, it retained (1) the intersection of the old
and new fine valid regions and (2) the exact coarse state stencil supplied to
the interpolater. It did not change the interpolater, state, GP extension,
flux, RHS, time step, or RK update. The deployed H100 executable had SHA-256
`5e449184a210ac9e35b4f7dfa4005c7e73fd3970974a92c45c9bdfe1d10ae4b5`.

The target L1 layout contained 7,326,208 valid cells. Of these, 7,263,744
overlapped the old L1 layout and 62,464 were newly constructed from L0. The
pre-flux audit reported:

```text
bad_total=12
bad_new_from_coarse=12
bad_old_fine_overlap=0
bad_unknown_origin=0
```

The 12 cells had 12 distinct, admissible coarse parents. Parent
`rho*e` ranged from `2.869e3` to `7.239e3 J/m^3`; the prolonged child range
was `-984.869` to `-12.479 J/m^3`. Every parent had exactly one inadmissible
child among its eight logged children. Thus no bad state was copied from L0
and no old L1 state was corrupted: the limited-linear child construction
itself crossed the Euler admissible-set boundary.

All 12 cells were pure fluid and lay at `x=0.01556--0.02231 m`,
`sqrt(y^2+z^2)=0.10517--0.11120 m`. This is far outside both the central
nozzle and the `0.0635 m` maximum body radius and is consistent with the
strong external bow-shock band. The points form four Cartesian azimuthal
clusters, and every bad point is the low-x child of its parent. This is a
multidimensional coarse-to-fine interpolation failure at an oblique curved
shock. It is separate from the Mach-disk odd-even deformation investigation.

For each parent, the arithmetic mean of the eight reconstructed children
recovered the parent conservative state to a maximum absolute residual of
`1.819e-12`. An offline admissibility calculation then evaluated

```text
U_child(theta) = U_parent + theta * (U_child_linear - U_parent).
```

The largest common parent-wise `theta` that keeps every child admissible was
`0.846596--0.998071` over the 12 affected parents. A common `theta` leaves the
zero-sum child deviations unchanged and therefore preserves the parent
average. The evidence supports a local conservative prolongation limiter,
not cellwise clipping and not a global first-order fallback.

Reproducible local records:

```text
analysis/test1853_single_3d/regrid_provenance_20260808/
  testD_regrid_provenance/run.log
  regrid_bad_cells.csv
  regrid_provenance_summary.json
  parse_regrid_provenance.py
  run_h100_replay.sh
```

The first three launcher attempts did not execute a CFD step: they stopped,
respectively, at an unsupported CUDA-aware OpenMPI request, a missing output
directory, and a relative STL path invalidated by the isolated run directory.
The recorded Test D result is the subsequent run using default host-staged
MPI, the required output directory, and the audited absolute STL path.

The appropriate production repair is a generic, conservative,
Euler-admissibility-preserving prolongation. For each parent refinement block,
scale all child conservative deviations from the admissible parent state by a
single `theta in [0,1]`, chosen to make every child admissible. Because the
limited-linear child deviations sum to zero, this retains the parent coarse
average; `theta=0` is the fail-closed piecewise-constant limit. This is not
independent state clipping and does not introduce cut cells, face fractions,
aggregate cells, or EB flow geometry.

## Verification status

- Complete 3-D CPU MPI rebuild: passed.
- Benign one-step diagnostic-only run: passed.
- Deliberately excessive-CFL run: correctly stopped at the failing FE bracket
  and printed the exact six-face masks and WENO/Rusanov flux differences.
- CUDA 12.6 native `sm_90` H100 build: passed.
- No-regrid frozen-checkpoint control through step 3061: passed.
- Original-regrid normal-step replay: deterministically reproduced the
  pre-flux inadmissible state.
- Same-regrid `1e-14 s` replay: reproduced identical pre-flux admissibility
  statistics.
- Same-regrid provenance replay: proved `12/12` failures are in newly
  prolonged L1 valid cells and `0/12` are in retained old-fine cells.
- Offline parent-block audit: verified child-average conservation to
  `1.819e-12` and bounded the required local common `theta` to
  `0.846596--0.998071`.

## Repair candidate and local regression

The generic production candidate now uses an
`EulerAdmissibleConservativeLinear` StateDescriptor mapper. It evaluates the
existing AMReX limited conservative linear children first. For each complete
parent refinement block, it then finds the largest common `theta` that gives
every child positive density and internal-energy density above
`max(p_min/(gamma-1), rho e_min)`. Here `e_min=0` in pressure-floor builds and
`e_min=c_v T_min` in temperature-floor builds. It applies

```text
U_child = U_parent + theta * (U_child_linear - U_parent)
```

to every child and every conservative component. A coarse parent that is
already inadmissible emits a failing state rather than being clipped. The
legacy and conservative-quartic mappers are rejected by pure-GP and EB
executables. The pure-GP method manifest remains
`pure_shared_gp_full_cartesian`; GP geometry, BI-CWLS, thermodynamic extension,
full-background-face WENO/LLF, SSPRK, reflux, and surface recovery are
unchanged. No cut-cell, cut-control, open-face, volume-fraction, aggregate, or
EB flow mechanism is referenced by the intended executable.

The direct CPU regression is stored in
`tst/amr/euler_admissible_prolongation/`. Its 3-D Cartesian branch embeds the
exact worst Test-D parent and six-neighbour stencil:

```text
legacy rho*e             = -984.86910627266479
protected minimum rho*e  =  5.4331030696630478e-07
theta                    =  0.84659599064979041
child-average residual   =  4.4408920985006262e-16
```

Its 2-D R--Z branch independently forces the same nonlinear failure in an
annular parent block and verifies metric conservation:

```text
legacy minimum rho*e       = -2.7222222222222232
protected minimum rho*e    =  2.5400000547293189e-08
theta                      =  0.65079136932344317
annular-average residual   =  0
```

Both branches also require bitwise transparency for uniform and nonuniform
smooth admissible states.

The strengthened temperature-floor branch of the frozen 3-D Cartesian test
uses `e_min=7000 J/kg` and gives:

```text
protected minimum e       = 7000.0000126397645 J/kg
theta                     = 0.80364958353778715
child-average residual    = 4.4408920985006262e-16
```

The corresponding R--Z thermodynamic-floor branch gives minimum
`e=1.0000000002000005` for a unit floor and retains exactly zero annular
average residual.
The complete 2-D R--Z pure-GP executable rebuild and its two-step level-zero
integration smoke passed with the new mapper and retained the required method
manifest. That smoke does not exercise coarse-to-fine interpolation; the
direct annular regression supplies that local coverage. A negative
configuration test with `cns.amr_state_interp=legacy_linear` aborted during
pure-GP setup as required.

A complete 3-D CPU MPI all-fluid executable was also rebuilt and run on the
small `pulse_amr_test` hierarchy. With `amr.regrid_int=1 1`, the log confirmed
regrids at initialization and after coarse step 1. Level 0 contained 1,024
cells and Level 1 contained 3,328 cells; two coarse steps and four fine
substeps completed without a state failure. This exercises the production
StateDescriptor, `FillCoarsePatch`, and regrid `FillPatch` wiring in 3-D. It
does not replace the pending 3-D pure-GP CUDA checkpoint replay.

## H100 Test E and Test F results

The first repair candidate was built natively with CUDA 12.6 for `sm_90` and
replayed from the frozen `chk03059` checkpoint on eight H100 GPUs. The build
and runtime both retained `pure_shared_gp_full_cartesian`; cut-control remained
disabled. The executable SHA-256 was
`04965bd1d27da24ab025a16628a6840bba5fd7be4cb211daf6f559cbdf9efd97`.
The complete local record is under
`analysis/test1853_single_3d/regrid_provenance_20260808/testE_admissible_prolongation/`.

Test E removed the original failure mode. There was no
`INADMISSIBLE_BEFORE_FLUX_EVALUATION` report, level 0 completed step 3060, and
no regridded cell had negative internal-energy density before flux evaluation.
The first finite L1 update then stopped the strict post-step audit with 18
cells below the compiled `T_min=10 K` closure floor:

```text
nonfinite                         = 0
strict-floor failures             = 18
minimum rho                       = 0.002807343791 kg/m^3
minimum rho*e                     = 8.619831296 J/m^3
floor at the minimum density      = 20.14982304 J/m^3
reported failing specific e range = 238.153075--3999.260767 J/kg
```

All 12 cells from Test D occur in this 18-cell set; the other six are immediate
neighbours in the same bow-shock interpolation band. Thus the first mapper did
not leave negative energy, but its `rho*e>0` boundary placed the repaired
children arbitrarily close to zero temperature. Primitive conversion then had
to apply its 10 K safeguard, and the first finite update remained below that
configured floor. Weakening `strict_positivity` would hide this state mismatch
and is not an acceptable repair.

The candidate has therefore been strengthened to use the density-dependent
thermodynamic floor described above. Local 3-D Cartesian, 2-D R--Z annular,
full 3-D two-level forced-regrid, and pure-GP R--Z level-zero regressions pass.
The full 3-D test regridded at initialization and after coarse step 1 and
completed two coarse plus four fine advances. The pure-GP R--Z smoke retained
the required method manifest and reported `cut_control_compiled=0`.

Test F completed that acceptance gate on eight H100 GPUs. The final executable
was built natively with CUDA 12.6 for `sm_90` from
`/shared/cerisse_local_0808_thermodynamic_prolongation`; its SHA-256 was
`90a733b431a9a27782862877ee29f3aa8025910c8571032ff0c4cc077f131a5d`.
The build manifest retained `pure_shared_gp_full_cartesian` and reported
cut-control disabled. Startup reported the configured prolongation floors

```text
rho  = 1e-19 kg/m^3
rho*e = 2.5e-8 J/m^3
e    = 7177.540244 J/kg
```

The forced restart regrid retained 7,263,744 old L1 cells, constructed 62,464
new L1 cells from coarse data, and constructed a new 22,436,864-cell L2. With
the normal checkpoint CFL path, it then completed coarse steps 3060 and 3061:

```text
L0 advances = 2,  dt = 2.232957255e-8 and 2.232986985e-8 s
L1 advances = 4
L2 advances = 8
exit status = 0
```

`cns.strict_positivity=1` remained active for every accepted finite update.
There was no `INADMISSIBLE_BEFORE_FLUX_EVALUATION`, no post-step closure-floor
failure, no non-finite state report, and no CUDA or signal failure. This closes
the deterministic `chk03059` regrid acceptance gate: the conservative mapper
now prevents both the original negative-`rho*e` children and the low-temperature
states exposed by Test E.

The first Test-F compile was deliberately rejected before calculation because
the copied `GNUmakefile` still pointed to the old admissible-prolongation source
tree. The corrected build now verifies the exact `AMR_SOLVER` path, cleans the
object directory, and fails if the build log references the stale tree. No CFD
step used the rejected executable.

The complete local record is under
`analysis/test1853_single_3d/regrid_provenance_20260808/testF_thermodynamic_prolongation/`.
This is a two-coarse-step deterministic replay, not a long-duration stability
claim. Any continuation beyond step 3061 remains a separate calculation.
