# L5-lip-only configuration for the x_fc profile sweep

Everything here is derived from `m142_2d_L5` on AWS
(`/shared/cerisse_reflux/exm/underexpanded_jet/2d/m142_2d_L5`). Domain, base
grid, boundary conditions, scheme, CFL and RK settings are unchanged, so the
results remain directly comparable with the existing 2D and 3D data.

## What changed and why

| item | production | here | reason |
|---|---|---|---|
| `amr.max_level` | 6 (24.8 µm) | **5 (49.6 µm)** | inlet-flux bias at 49.6 µm is ≤ 3.4e-5 D_e |
| L5 source | sensor + shear corridor to 6 D_e | **forced, lip band + axis strip only** | sensor-driven L5/L6 was 86% of the cost |
| L6 | sensor to 4 D_e | **removed** | carries no x_fc signal |
| `amr.n_error_buf` | 8 8 8 8 8 8 | 4 4 4 4 2 | forced zones need no halo |
| `amr.blocking_factor` | 32 | 32 32 32 16 16 16 | less padding in the thin band |
| `amr.grid_eff` | 0.75 | 0.85 | tighter boxes |
| `amr.regrid_int` | 20 | 40 | L5 is static |
| `amr.plot_int` | 50 | −1 | `main.cpp:96` writes a final plotfile anyway; periodic I/O was 5–32% of wall time |
| `cns.record_probe` | 1 (32 probes, every step) | 0 | unused by this study |

## Cost

Measured production cost from `m142_gpu_prod/phaseB.log`:

| level | dx | cells | ×2^l | cell-updates | share |
|---|---|---|---|---|---|
| 0 | 1587.5 µm | 98,304 | 1 | 0.10 M | 0.04% |
| 1 | 793.8 µm | 98,304 | 2 | 0.20 M | 0.08% |
| 2 | 396.9 µm | 183,296 | 4 | 0.73 M | 0.31% |
| 3 | 198.4 µm | 569,344 | 8 | 4.55 M | 1.95% |
| 4 | 99.2 µm | 1,642,496 | 16 | 26.3 M | 11.3% |
| 5 | 49.6 µm | 2,649,088 | 32 | 84.8 M | 36.4% |
| 6 | 24.8 µm | 1,821,696 | 64 | 116.6 M | 50.0% |
| | | | | **233.2 M** | |

This configuration: **15.3 M cell-updates per coarse step, a 15× reduction.**

## Run protocol

```
phase 1   ./main2d.gnu.CUDA.ex inputs_L5lip cns.record_stats=0 stop_time=4.0e-4
phase 2   ./main2d.gnu.CUDA.ex inputs_L5lip cns.record_stats=1 \
              amr.restart=./plot/chk????? stop_time=8.0e-4
```

`CNS.cpp:1515` zeroes both the stats arrays and `time_stat_level` on restart,
so the phase-2 mean is a clean window over [0.4, 0.8] ms with no startup
transient. One GPU per case — at ~1.1 M total cells, multi-GPU is pure
overhead.

## Must be validated before the sweep is trusted

1. **Averaging window** (zero compute): recompute x_fc from the existing
   `m142_2d_L6stat` statistics over different sub-windows. Sets `stop_time`.
2. **Grid bias** (one case): run a = 254 µm, s = 3 on this configuration and
   compare with the existing full-grid value x_fc = 0.883. Accept if
   |Δx_fc| < 0.002 D_e, i.e. under 4% of the 0.052 D_e s-ladder signal.

## Honest limitation

49.6 µm resolves the prescribed inlet flux essentially exactly, but not the
Kelvin–Helmholtz roll-up (δ_ω spans 4–14 cells across the matrix). That biases
the *absolute* x_fc, not the *slope*: along the s-ladder δ_ω is held fixed at
340 µm, so the under-resolution is common-mode and cancels in the differences
that the study measures. Absolute x_fc from this configuration should not be
compared with the L6 production values or the experiment until check 2 above
has quantified the offset.
