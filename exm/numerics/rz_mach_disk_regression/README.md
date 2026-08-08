# R-Z Mach-disk and shock regression

## Current decision

The 2026-07-23 controlled Mach-10 matrix passed 15/15 checks, but the real SRP
checkpoint v2 A/B failed four frozen gates. An exact-time baseline decomposition
also showed that the annular/point-flux control amplifies the first-ring axis
artifact and v2 amplifies it further. V2 is rejected for production use and
remains off by default.

The requested four-way control split and per-face WENO/LLF audit are complete.
They identify annular recovery as the largest immediate near-axis RHS
perturbation. A separate additive H-correction now retains the full WENO flux
and adds only smooth transverse Rusanov dissipation. Its planar calibration and
real one-stage flux identity checks pass at strength `0.125`. The exact-time
250-step real-SRP pair formally passes every frozen pre-result gate, but does
not demonstrate a material repair. The apparent 9.16% reduction is only
`0.00451` L4 cell (`0.705 micrometre`) because the baseline first-eight-cell
front range was already near zero. The full fields and profiles are visually
indistinguishable. The engineering decision is therefore `NO-GO`: do not use
this exact `H=0.125` formulation for a pre-defect long restart. It remains
experimental, default-off, and unpromoted. See
[H_ENGINEERING_DECISION.md](H_ENGINEERING_DECISION.md) for the decision and
[RESULTS_20260723.md](RESULTS_20260723.md) for metrics, plots, provenance, and
limitations.

A later exact-time audit of the original NPR15 `plt05098` field corrected the
target definition. The apparent 20-cell axis cliff was produced by a
max-density-gradient detector switching from a post-disk subsonic compression
at `z/D=2.47557` in the first three radial rings to the actual Mach disk at
`z/D=2.39648` from the fourth ring onward. The physical Mach disk is the first
compression with upstream Mach at least 1.5 and downstream Mach at most 1.5.
Under that criterion its `r/D<=0.08` core range is `0.768` L4 cell, and H
changes it by only `-0.039` cell. The old 37-cell/32.8-cell values are retained
only as detector-error provenance and must not be cited as Mach-disk geometry.

The remaining visible cap is a genuine post-disk recirculation/stagnation
compression: on the axis, axial velocity changes from about `+145 m/s` to
`-89 m/s` across it. It is not the Mach disk. The next discriminating work is
a grid study and full-azimuth 3-D audit of this coherent axisymmetric
recirculation, not another axis pressure-source or H-strength adjustment. See
[NPR15_ROOT_CAUSE_20260723.md](NPR15_ROOT_CAUSE_20260723.md) for the complete
identity audit and corrected acceptance metrics.

This regression has two deliberately separate layers.

1. `../rz_planar_shock/run_matrix.sh` is the controlled Mach-10 mechanism
   test. It checks both flow directions, exact planar preservation, a seeded
   two-cell shock displacement, three radial resolutions, positivity, and the
   directed-filter bypass.
2. `run_real_ab.sh` restarts the developed SRP `chk40000` field and compares
   the same annular point-flux shock-gate control against that control plus the
   directed transverse v2 filter. All other command-line settings are shared.

The real-geometry continuation is a diagnostic regression of a legacy
low-storage-RK2 checkpoint. It requires the runtime IBM manifest
`pure_shared_gp_full_cartesian`, but it is not a qualification result for the
later frozen SSPRK(4,3) production method. No cut-cell, cut-control,
aggregate, redistribution, or AMReX-EB flow path is enabled.

Typical Y9000X execution:

```bash
MPI_ROOT=/home/melampous/cerisse_rz_bench_20260720/ompi_user/root/usr \
NP_PER_CASE=8 CONTROL_CPUSET=0-7 CANDIDATE_CPUSET=8-15 \
./run_real_ab.sh \
  /path/to/main2d.gnu.MPI.ex \
  /path/to/inputs_aws_steady \
  /path/to/chk40000 \
  /path/to/results/real_ab_250
```

The runner refuses to reuse an output directory, records the exact commands
and SHA256 inputs, checks the pure-GP runtime contract, and requires a clean
AMReX finalization for both continuations. Multi-branch runners also read the
final physical time from each plot Header and refuse to write
`EXECUTION_PASS` if the branches differ by more than `1e-12 s`; equal step
counts are not treated as equal physical times under adaptive CFL stepping.

`run_control_decomposition.sh` keeps the frozen A/B runner unchanged and
separates the current control into `baseline`, `annular_only`,
`annular_point`, and `annular_point_gate`. With `RUN_FACE_DIAG=1`, the first
RHS evaluation on L4 also writes one CSV row per owned radial face. Each row
contains the point-flux/gate masks, shock sensor, incident-flow gate, LLF
wavespeed, actual WENO flux, alternate first-order LLF flux, and their
componentwise difference. The diagnostic is default-off.

The frozen-checkpoint one-stage audit contains 16,256 owned L4 radial faces.
The directed H-stencil shock sensor triggers on 510 faces. The point-flux gate
would block 526 faces, including all 510 shock-triggered faces; this directly
supersedes the unsupported historical claim that only 17 of 96 columns
triggered for this exact checkpoint. Both masks reach all 128 recorded radial
columns, but the directed
mask occupies 8 axial rows while the point gate occupies 11. The 16 gate-only
faces include the separate strong-jump cluster near `z=-123 mm`, which the
incident-flow direction check excludes. In the first eight radial cells, the
radial-momentum RHS increments have RMS values:

| isolated transition | RMS RHS increment | absolute maximum |
|---|---:|---:|
| annular minus baseline | 1.2031e4 | 8.1215e4 |
| point minus annular | 4.5092e3 | 3.8174e4 |
| gate minus point | 1.0117e3 | 2.1452e4 |

Thus annular midpoint recovery is the largest immediate near-axis perturbation
in this snapshot, point flux is second, and the gate is smaller but nonzero.
This is a one-RHS-stage causal decomposition, not by itself a claim about the
fully evolved root cause. `analyze_face_audit.py` archives the masks,
WENO/LLF flux differences, and all three annular-divergence increments.

`run_h_ab.sh` is the independent additive H-correction A/B. It explicitly
disables annular recovery, point flux, its gate, and the complete-face LLF
filter in both branches. The candidate retains the full WENO flux and adds

```text
-0.5 * weight_H * alpha_r * (U_R - U_L)
```

only on a directed shock-triggered radial face. `RUN_FACE_DIAG=1` additionally
records `h_correction_active`, `h_correction_weight`, the componentwise H flux
increment, and the final selected flux. `analyze_h_face_audit.py` verifies that
the base WENO flux is unchanged and that the selected flux equals WENO plus
the recorded correction. The controlled Mach-10 formal matrix selects
`H_STRENGTH=0.125`; stronger values `0.25` and `0.5` are rejected for
high-resolution vorticity amplification.

For a single decomposition branch, set `ONLY_VARIANT` to `baseline`,
`control`, or `candidate`. `STOP_TIME` can pin an independently computed
branch to the exact physical endpoint of the pair:

```bash
ONLY_VARIANT=baseline STOP_TIME=4.30332634858483e-3 \
BASELINE_CPUSET=0-7 ./run_real_ab.sh BINARY INPUTS CHECKPOINT RESULT_ROOT
```

The result is evaluated with `analyze_real_ab.py` at L4 around the Mach disk
and L3 along the bow shock. It records the subcell pressure-gradient front,
center-to-annulus indentation, first-8/first-16 front range, radial velocity,
azimuthal vorticity, positivity, and bow position/gradient. `field_ab.png`
uses identical color limits for the two branches.

`verify_real_ab.py` contains the predeclared acceptance gates. In particular,
the candidate may not amplify either axis artifact by more than 5%, must damp
at least one by 10%, and must reduce center indentation by 10% unless its
absolute magnitude is already at most 0.25 mm. It may not materially worsen
the first-eight- or first-sixteen-cell front range and must preserve the
bow-shock position and gradient. A failed gate is retained as a failed
diagnostic; thresholds are not retuned after inspecting the result. Passing
these frozen relative gates is necessary but was not sufficient here: the
baseline defects were too small for the relative reductions to establish a
physically meaningful repair.
