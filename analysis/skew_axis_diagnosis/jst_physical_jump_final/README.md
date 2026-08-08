# Run165 RZ Skew-JST axis fallback audit

## Scope

This directory contains read-only analysis and plotfiles for the analytic-exit,
jet-on Run165 RZ case without IBM geometry.  No production source was edited by
this audit.

The final K=0 runs used:

- `src/rhs/Skew.h` SHA-256
  `9d50e7ca27e6671fc19a14b6e243658dc64b925c35bb0910738b34f3f48898c3`;
- `src/tim/advance.cpp` SHA-256
  `9f7a288bfe00f6b32fb2d6ba2e1d4338784d4a72ba46beaa74816b106ba314e3`;
- executable `run165_skew_jst_k0stable2d.gnu.TPROF.OMP.ex` SHA-256
  `aa5ce66b20d453e908151ae43b9b099d0f3c58564ab5c52172e85764ad551492`;
- input `inputs_compare_fixed_dt` SHA-256
  `25321386f5bf1b028b852085beb2f636ec27d9243e3021b851d46efbd80fc74b`.

## Build and run

From `exm/underexpanded_jet/2d/run165_analytic_exit_matrix`:

```bash
make -j8 EULER_SCHEME=skew SKEW_C2=1.5 SKEW_C4=0.016 \
  USE_MPI=FALSE USE_OMP=TRUE \
  EBASE=run165_skew_jst_k0stable \
  TMP_BUILD_DIR=tmp_build_run165_skew_jst_k0stable
```

The 500 us gate used `max_step=1250`; the 1.5 ms extension used
`max_step=3750`.  Both used the following overrides:

```text
cns.stage_positivity=0
cns.skew_first_order_fallback=1
cns.skew_axis_collar_faces=0
amr.checkpoint_files_output=0
```

The fixed coarse step is 0.4 us.  The full-domain uniform L1 field has
64 x 192 = 12,288 finest-level cells and `dr=dz=1.06129666667 mm`.

## Stability and cost

Both gates completed normally with strict positivity enabled and no NaN,
non-finite-state, or positivity abort:

| gate | coarse steps | physical time | wall time | wall/coarse-step |
|---|---:|---:|---:|---:|
| K=0 | 1250 | 500 us | 24.4653 s | 0.01957 s |
| K=0 | 3750 | 1500 us | 72.7972 s | 0.01941 s |

The run-time method manifest reported Skew-JST-O4, C2=1.5, C4=0.016,
fallback enabled, `rz_axis_collar_faces=0`, paired metric advection/pressure
gradient, Euler equations, jet on, and no IBM/LES.

## Main result

K=0 removes the deterministic hard-collar transition seen with K=1, 2, and 4.
In the near-nozzle interval `-60 <= z <= -20 mm`, the worst p95 radial
local-extremum residual is:

| setting at 500 us | pressure | density | temperature | radial velocity | Mach | dominant ring trend |
|---|---:|---:|---:|---:|---:|---|
| fallback off | 0.08213 | 0.04784 | 0.06384 | 0.09896 | 0.02853 | ring 1/2 |
| fallback on, K=0 | **0.07210** | **0.04683** | **0.05492** | **0.09299** | **0.02407** | ring 1/2 |
| fallback on, K=1 | 0.09730 | 0.05698 | 0.05997 | 0.16619 | 0.03005 | ring 2/3 |
| fallback on, K=2 | 0.09287 | 0.05402 | 0.04568 | 0.13224 | 0.02470 | ring 3 |
| fallback on, K=4 | 0.08204 | 0.05498 | 0.03222 | 0.08886 | 0.01777 | ring 5 |

Thus widening a static collar only moves the strongest stripe outward.  K=0
retains the generic troubled-face fallback and invokes the paired axis rescue
only when the actual first interior face is selected.

The K=0 near-nozzle residual is stationary between 500 us and 1.5 ms to the
reported precision: pressure p95 0.0720965 -> 0.0720975 and Mach p95
0.0240683 -> 0.0240986.  It neither grows nor migrates to a fixed outer ring.

This is a substantial improvement, but not an axis-accuracy proof.  Over
`-80 <= z <= -5 mm`, the K=0 maximum parity defects at 500 us are 12.18% in
pressure, 9.00% in density, 11.73% in temperature, 0.769% in axial velocity,
and 2.46% in radial velocity normalized by sound speed.  K=0 is therefore the
least intrusive and most stable option in this Run165 gate, while a separate
smooth-axis accuracy/regression gate is still required before calling the
axis treatment fully validated.

## Evidence

- `fallback_on_k0_500us/plt01250`: final 500 us plotfile.
- `fallback_on_k0_1500us/plt03750`: final 1.5 ms plotfile.
- `fallback_on_k0_500us/figures/mach_pressure_500p00us.png`.
- `fallback_on_k0_1500us/figures/mach_pressure_1500p00us.png`.
- `collar_sweep_audit/collar_sweep_heatmap_500us.png`: common-scale
  off/K0/K1/K2/K4 comparison.
- `collar_sweep_audit/near_nozzle_minus60_to_minus20mm.csv`: collar-sweep
  metrics.
- `k0_long_audit/k0_heatmap_500_1500us.png`: common-scale K=0 time comparison.
- `k0_long_audit/axis_parity_metrics_near_exit.csv`: parity metrics.

The fallback-off and K=1/2/4 reference plotfiles were generated with the
immediately preceding Skew source SHA (`58a1bb...`).  They are valid for the
hard-collar migration diagnosis, but an exact bitwise same-source comparison
would require rerunning those reference settings with SHA `9d50e7...`.
