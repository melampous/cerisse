# Pure GP-IBM Mach-4 Circle Long-Run Gate

Status date: 2026-07-22

## Scope

This is a nonlinear software-regression gate for the frozen pure GP-IBM
operator. It is not a validation against an exact two-dimensional Euler
circle solution. Billig and Hornung values are reported only as literature
diagnostics.

The tested method is:

- shared solid-side ghost-point extension;
- cell-average-aware BI-CWLS with expanded smooth-wall support;
- one-sided entropy extension and curvature-compatible pressure closure;
- one ghost-cell layer;
- full Cartesian cell volumes and face areas;
- LLF-WENO-Z5;
- SSPRK(4,3) with conservative shared-face positivity limiting;
- uniform level-0 grid and fixed Euler-slip circle.

The runtime manifest reported
`family=pure_shared_gp_full_cartesian`, `pure_gp_production=1`, and
`cut_control_compiled=0`. A defined-symbol audit of the executable found no
cut-control, aggregate, volume-fraction, or open-subface symbol.

## Provenance

| Item | Value |
|---|---|
| Git commit | `dc35646d8eb90aaa362a3292a9d702e30c9d5af5` |
| Tracked dirty diff SHA-256 | `604c6323ebf00368df42bbb9bdcf9da39441070d8af41c8a483709e39c4f5a64` |
| Source content ID | `5b45451c364d120e` |
| Executable SHA-256 | `1c45397fbbbd694c0d80b3b56ba105b7c411994426c1921cc93f0db340c74bb0` |
| Inputs / geometry SHA-256 | `81a058287a7db070cb4f717ff73c2709f1ea2f6e8023d79ca6005542863fb9c9` / `53775b9fc2a2b67a354b72ca4fcbb746810656d98a926c5effec123abb424e6e` |
| GPU | NVIDIA GeForce RTX 4090 Laptop GPU, `sm_89` |
| Geometry | `circle_r0p1_4096.dat`, radius `R=0.1`, 4096 segments |
| Domain | `[-0.5,0.5]^2` |
| Grid | `320 x 320`, `D/h=64`, `max_grid_size=256`, `blocking_factor=64` |
| Flow | ideal-gas Euler, `M_inf=4`, `gamma=1.4` |
| Volume / surface orders | `iorder=2, eorder=1` for both paths |
| CFL | 0.3 |
| End time | `t=0.008`, `t*=t U_inf/D=55.560026` |
| Steps / wall time | 15,287 / 820.34 s |
| Run directory | `IBM/cases/ibm_tests/2d_bvh_gpu/production_freeze_20260722/circle_N320_pure_gp_src5b45451c364d120e` |

The read-only limiter-location audit is disabled by default. Audit-off and
audit-on runs were bitwise identical through 12 and 1000 steps, respectively,
so the diagnostic does not alter the update.

`5b45451c364d120e` is a build-system source-content ID, not a Git commit or a
plain dirty-diff digest. It hashes all files under `src/`, including untracked
files, together with the circle `prob.h` and `GNUmakefile`. Full definitions
and artifact digests are in `PROVENANCE.md` and `SHA256SUMS` in the run
directory.

## Bow-shock result

The detached bow shock remained present and stable through the prescribed
end time.

| Metric | Result |
|---|---:|
| Final stand-off distance, `Delta/R` | 0.5496841170717 |
| Stand-off distance in cells | 17.58989175 |
| Density/pressure locator span | 0.10011365 cells |
| Last-four-point relative plateau span | `9.1391e-8` |
| Difference from frozen source-75 stand-off | `6.46e-14` relative |
| Difference from Billig correlation | +6.357% |
| Difference from Hornung fit | +5.846% |

The Billig and Hornung differences are not interpreted as numerical error,
because neither is an exact matched Euler reference for this discrete
problem.

## Positivity-limiter audit

The local shared-face limiter activated in 74 SSPRK Forward-Euler brackets.
Every associated all-low-order update was admissible, and every limited
update passed the post-stage admissibility audit.

| Metric | Result |
|---|---:|
| Activation-time interval | `7.6512354e-5 <= t <= 9.5449388e-5` |
| Convective-time interval | `0.531379 <= t* <= 0.662896` |
| Minimum face coefficient | `theta_min=0.2182016715` |
| Limited-face occurrences | 1,295 |
| Crossing-face occurrences | 92 |
| Maximum relative shared-face conservation residual | `1.7540e-15` |
| Activations after `t*=0.662896` | 0 |
| Global-theta / final all-low fallback / rejected step | 0 / 0 / 0 |

The occurrence counts are accumulated over FE brackets and are not counts of
unique geometric faces. The location envelope was:

- x-normal faces: `x in [-0.103125,0.1]`,
  `y in [-0.0578125,0.0578125]`;
- y-normal faces: `x in [0.0921875,0.1171875]`,
  `y in [-0.04375,0.04375]`.

The strongest correction occurred at approximately `(0.0953,0.0313)`, on the
leeward side of the circle. Thus the limiter was confined to the initial
wall/shock/wake establishment transient, with a substantial aft-side
component; it was not a persistent limiter acting along the steady wall or
bow shock.

## Frozen-baseline comparison

The previous source-75 run used the earlier admissibility semantics and
activated the limiter in 80 FE brackets. The current run activated it in 74,
so the complete wake field is not expected to be bitwise identical after the
first activation.

The physically relevant forebody comparison is nevertheless essentially at
roundoff:

| Region/quantity | Relative difference |
|---|---:|
| Forebody density `L2` | `2.40e-13` |
| Forebody pressure `L2` | `3.10e-13` |
| Stagnation shock-layer pressure `L2` | `1.98e-13` |
| Shock-front RMS displacement | `1.25e-12` cells |
| BI pressure, forebody `L2` | `9.45e-14` |

Differences are concentrated in the unsteady aft wake and rear-half wall
pressure. They should not be hidden by a whole-domain bitwise claim. The
current global symmetry defects are smaller than in the frozen run: pressure
even-symmetry error decreased from 0.0255 to 0.00852, and density from 0.0278
to 0.00941.

The e1 BI-integrated axial pressure force changed from 282777.78 to
282844.92, or +0.0237%. This is a finite-grid software diagnostic, not an
exact force reference and not a revalidation of the production-target
surface-e2 observation operator.

## Decision

The Mach-4 circle flow/positivity sub-gate is **PASS** for source
`5b45451c364d120e`:

- the detached bow shock exists;
- its stand-off distance and front shape reproduce the frozen pure-GP result;
- the final stand-off history is stationary;
- positivity limiting is conservative and restricted to the early transient;
- no cut-control or aggregate operator is present.

This does not close Gate 0 by itself. The circle binary used surface e1, so
its Cp output does not close the separate surface-e2 observation gate. The
frozen PM expansion, single compression corner, and inverse-MOC phase checks
remain required. No IBM formula change is justified by this circle result.

The run did not write a checkpoint. Consequently, surface e2 cannot be added
by a legitimate zero-step restart from this run, and the final plotfile must
not be presented as a checkpoint. Repeating 15,287 steps only to change an
observation operator is not justified. The remaining surface-e2 gate will be
closed by the already required inverse-MOC two-phase runs.

## Artifacts

- `circle_gate.json`: machine-readable gate decision;
- `plt*_standoff.json`: subcell shock-location histories;
- `standoff_history.png` and `standoff_history.csv`: current/frozen history;
- `run.log`: complete manifest, limiter and profiler record;
- `PROVENANCE.md`, `RUN_COMMAND.txt`, and `SHA256SUMS`: reproducibility data;
- `plt15287`: final field;
- `circle_r0p1_4096_15287.vtp`: final BI surface fields.
