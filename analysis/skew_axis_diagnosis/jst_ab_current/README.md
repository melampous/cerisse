# Run165 Skew central/JST orthogonal A/B

Date: 2026-08-07

## Scope and invariant setup

This is the all-fluid R--Z Run165 analytic-exit case: jet enabled, Euler,
no IBM, no LES, uniform full-domain L1, fixed `dt=4.0e-7 s`.  Both central
branches use the same fourth-order Skew operator compiled with `C2=0,C4=0`;
the sole runtime A/B variable is `cns.skew_first_order_fallback=0/1`.

The L1 plotfile header confirms a `64 x 192` finest grid with
`dr=dz=1.06129666667 mm` and no AMR interface.  The old input value
`cns.stage_positivity=1` is now rejected because it was historically ignored,
so every run explicitly sets it to zero.  `cns.strict_positivity=1` remains
enabled; no state clipping or IBM path is used.

Frozen production-source/build hashes:

```text
Skew.h       5de943653e1cbac15a9a49ede5deb7c7f470fe736073f68e29ed8a3fdc88b7d0
advance.cpp  9f7a288bfe00f6b32fb2d6ba2e1d4338784d4a72ba46beaa74816b106ba314e3
central exe  ec466e7530a8310293f7387e4a8e0845019eff8b8f0d46346430327bee6c46b2
JST exe      cf601ba06c57f45d59a5a8e4a210ea54f6072f4576c972363be7fdb5a90733ac
inputs       25321386f5bf1b028b852085beb2f636ec27d9243e3021b851d46efbd80fc74b
```

Production source was not edited by this audit.  Only this analysis script,
CSV, figures, logs, and plotfiles were generated.

## Build and run commands

```text
make -j8 EULER_SCHEME=skew SKEW_C2=0 SKEW_C4=0 \
  USE_MPI=FALSE USE_OMP=TRUE \
  EBASE=run165_skew_central_o4_jst_ab_current_v2 \
  TMP_BUILD_DIR=tmp_build_run165_skew_central_o4_jst_ab_current_v2
```

The two 80-us commands differed only in the value shown as `FALLBACK`:

```text
OMP_NUM_THREADS=8 timeout 120 \
  ./run165_skew_central_o4_jst_ab_current_v22d.gnu.TPROF.OMP.ex \
  inputs_compare_fixed_dt max_step=200 cns.stage_positivity=0 \
  cns.skew_first_order_fallback=FALLBACK cns.nstep_screen_output=25 \
  amr.v=0 amr.plot_files_output=1 amr.plot_int=200 \
  amr.plot_file=RESULT_DIRECTORY/plt amr.checkpoint_files_output=0
```

The matched current-source JST reference used the same command and
`run165_skew_axis_parity_final2d.gnu.TPROF.OMP.ex`, which reports
`C2=1.5,C4=0.016` and fallback enabled.

## Gate result

| Operator | Fallback | Result | Final time | Wall time |
|---|---:|---|---:|---:|
| pure central O4 | 0 | failed at step 11 | 2.4 us | about 1.5 s including abort |
| pure central O4 | 1 | passed 200 steps | 80 us | 3.587 s |
| Skew-JST O4 | 1 | passed 200 steps | 80 us | 3.885 s |

Because both pure-central branches did not pass the 80-us gate, neither was
extended to 500 us.  This follows the pre-declared gate and avoids treating a
short surviving branch as a developed-flow result.

Fallback-off first loses internal energy at finest-level cell `(3,191)`,
`r=3.714538 mm`, `z=1.910620 mm`, inside the jet aperture and one cell below
the exit boundary, not on the axis and not at an AMR interface.  The state is

```text
rho=0.02985603785 kg/m3
momentum=(-5.05005e-5, 23.59734731, 0) kg/(m2 s)
E=9139.129902 J/m3
rho*e_int=-186.1999134 J/m3
```

The surrounding top-row rings 0--5 also have negative internal energy in the
failed candidate.  Thus, without JST and without local LLF replacement, the
central update creates kinetic energy larger than total energy next to the
time-ramped jet boundary almost immediately.  Enabling the local fallback is
sufficient to cross this early positivity failure.

## Matched 80-us first-eight-ring result

The table gives the worst radial centre ring for the 95th-percentile absolute
local-extremum residual over `-80 <= z <= -5 mm`.

| Field | central + fallback | JST + axis-pair + fallback | ratio central/JST |
|---|---:|---:|---:|
| pressure | 15.62% (ring 2) | 10.93% (ring 1) | 1.43 |
| density | 9.74% (ring 3) | 9.76% (ring 1) | 1.00 |
| temperature | 10.71% (ring 1) | 4.78% (ring 2) | 2.24 |
| Mach | 13.81% (ring 1) | 4.03% (ring 3) | 3.43 |

At the same time, pure central has `Mach_max=4.7486` and
`p_min=522.79 Pa`, whereas the JST reference has `Mach_max=4.6059` and
`p_min=531.02 Pa` (freestream values are 4.6 and 534.58 Pa).  The inferred
axis Mach-one location differs by 1.285 mm: `-30.354 mm` central versus
`-29.069 mm` JST.  The central profiles and contour lines show visibly larger
ring-to-ring and dispersive oscillation despite remaining positive to 80 us.

The previously generated 500-us files provide context, not a matched-source
or matched-time A/B.  Their worst residuals are:

| Older 500-us reference | pressure | density | temperature | Mach |
|---|---:|---:|---:|---:|
| JST baseline | 6.53% | 5.84% | 8.84% | 4.30% |
| first-face axis-pair | 13.80% | 10.39% | 3.98% | 2.04% |

The 80-us matched comparison is therefore the causal evidence: the selected
LLF fallback supplies enough local positivity robustness for pure central,
but it does not replace the background/shock-selective JST damping needed to
control the early SRP wave system.  A 500-us pure-central field was not
qualified by this gate.

## Evidence

- `central_fallback_off_80us/run.log`: exact failed state and stencil.
- `central_fallback_on_80us/plt00200`: surviving pure-central field at 80 us.
- `jst_axis_pair_fallback_on_80us/plt00200`: matched current JST field.
- `radial_filament_metrics_80us.csv`: all fields and centre rings 1--6.
- `first8_profiles_80us.png`: first eight radial profiles.
- `radial_filament_heatmaps_80us.png`: fixed-ring local-extremum map.
- Each successful case also has `figures/mach_pressure_80p00us.png` and
  `figures/axis_profiles_80p00us.png`.

