# Skew RZ repair validation

Date: 2026-08-07

## Frozen source and executable

- `src/rhs/Skew.h` SHA-256:
  `f8f97b3e412feb8ea6b7c349acb08b30a994e1a66658eaab7405b47e272d5b3b`
- `run165_phase2_skew_o4_jst2d.gnu.TPROF.OMP.ex` SHA-256:
  `49369de9e8b24dc1aa84e5cb923048189221359dc989f7770377a73eadc55f2e`

The case uses the all-fluid RZ Run165 analytic-exit configuration in
`inputs_compare_fixed_dt`. The fixed coarse-grid time step is
`4.0e-7 s`. The inviscid operator is fourth-order Skew with JST
coefficients `C2=1.5` and `C4=0.016`.

## Results

| Local LLF fallback | Coarse steps | Final physical time | Exit status |
|---|---:|---:|---:|
| enabled | 200 | `80 us` | 0 |
| disabled | 200 | `80 us` | 0 |
| enabled | 3750 | `1.5 ms` | 0 |

All three runs used `cns.strict_positivity=1`. No state check failed and no
NaN or infinity was reported. The final paired-pressure RZ discretisation
therefore crosses the former early axis failure even with the local fallback
disabled. The fallback remains enabled by default as an additional
troubled-face robustness measure.

The 1.5 ms command was:

```text
OMP_NUM_THREADS=8 ./run165_phase2_skew_o4_jst2d.gnu.TPROF.OMP.ex \
  inputs_compare_fixed_dt max_step=3750 \
  cns.skew_first_order_fallback=1 \
  cns.nstep_screen_output=500 amr.v=0 \
  amr.plot_files_output=0 amr.checkpoint_files_output=0
```

This is an engineering stability regression. It is not a proof of strict
positivity or invariant-domain preservation.
