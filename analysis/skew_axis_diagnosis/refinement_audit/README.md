# Skew RZ K=0 near-axis refinement audit

## Conclusion

The thin line remaining in the uniform-L1 K=0 result is predominantly a
coarse-grid, axis-attached discretization error.  It is not a resolved physical
curvature at a fixed radius.  Halving `dr=dz` makes the line much weaker, and
the thermodynamic/Mach peak that remains closest to it stays at centre-ring 1:
its physical radius therefore halves from 1.592 to 0.796 mm.  At the original
L1 peak radius, the pressure, density, temperature, radial-velocity and Mach
local-extremum residuals fall by factors 5.19, 4.03, 4.09, 3.42 and 3.95,
respectively.

This two-level test is strong evidence about the visible line, but it is not a
formal asymptotic-order proof for the complete discontinuous flow.  Smooth
physical radial curvature remains in the L2 solution, especially around
`r=3--6 mm`, and the parity defects are reduced rather than identically zero.

## Reproducible setup

Both principal snapshots use the same analytic-exit, jet-on Run165 RZ problem,
Skew4-JST (`C2=1.5`, `C4=0.016`), generic local fallback enabled, static axis
collar disabled (`cns.skew_axis_collar_faces=0`), no IBM geometry, and the same
physical time `t=200 us` (agreement better than `3e-18 s`).

| case | target level | finest field | `dr=dz` | coarse steps / `dt` | finest `dt` |
|---|---:|---:|---:|---:|---:|
| uniform L1 | 1 | 64 x 192 = 12,288 cells | 1.0612967 mm | 500 / 0.4 us | 0.2 us |
| uniform L2 | 2 | 128 x 384 = 49,152 cells | 0.5306483 mm | 1,000 / 0.2 us | 0.05 us |
| L1 time control | 1 | 64 x 192 = 12,288 cells | 1.0612967 mm | 2,000 / 0.1 us | 0.05 us |

All three runs completed without a non-finite-state, strict-positivity or NaN
abort.  The L2 run took 235.05 s on this workstation; it advanced 49,152 L2
cells on each of 4,000 fine substeps.  Timing should not be used as a scaling
comparison because another short job overlapped part of this run.

Executable SHA-256:
`aa5ce66b20d453e908151ae43b9b099d0f3c58564ab5c52172e85764ad551492`.
Input SHA-256:
`25321386f5bf1b028b852085beb2f636ec27d9243e3021b851d46efbd80fc74b`.

The abandoned 500-us L2 attempt was stopped at 141.85 us and has no terminal
plotfile; it is not used anywhere in this audit.

## Line-position and thickness test

The diagnostic is the normalized radial local-extremum residual

```text
|q_j - (q_{j-1}+q_{j+1})/2| / max(|q_{j-1}|,|q_j|,|q_{j+1}|,floor)
```

over `-60 <= z <= -20 mm`.  The table reports the centre ring with the largest
axial p95 residual.  Width is the contiguous physical radial width whose p95
exceeds 2%.

| field | L1 peak `r` / p95 / width | L2 peak `r` / p95 / width |
|---|---:|---:|
| pressure | 1.592 mm / 7.277% / 5.306 mm | 3.980 mm / 2.606% / 3.184 mm |
| density | 1.592 mm / 4.899% / 5.306 mm | 0.796 mm / 2.491% / 0.531 mm |
| temperature | 1.592 mm / 5.799% / 3.184 mm | 0.796 mm / 3.166% / 0.531 mm |
| radial velocity | 2.653 mm / 9.571% / 5.306 mm | 2.388 mm / 2.958% / 5.306 mm |
| Mach | 1.592 mm / 2.670% / 1.061 mm | 0.796 mm / 1.468% / 0 mm |

The L2 pressure maximum at 3.980 mm is a different, much weaker part of the
smooth curved field.  At the *same physical radius* as each L1 maximum, the L2
p95 residual is smaller by the factors quoted in the conclusion.  Thus the L1
stripe does not converge to a finite-radius feature.  Mach has no L2 ring whose
p95 residual exceeds 2%.

## Axis even/odd defects

Even scalar parity is tested by a fourth-order `q(r)=a+b r^2+c r^4`
extrapolation to the first centre.  Odd radial-velocity parity is tested by an
even extrapolation of `u_r/r`; its residual is normalized by local sound speed.
The interval is `-80 <= z <= -5 mm`.

| field | L1 max / p95 | L2 max / p95 |
|---|---:|---:|
| pressure, even | 12.233% / 11.840% | 3.966% / 2.589% |
| density, even | 9.173% / 8.667% | 6.520% / 5.245% |
| temperature, even | 11.452% / 11.326% | 8.746% / 6.532% |
| axial velocity, even | 9.447% / 7.472% | 9.522% / 6.404% |
| radial velocity / sound, odd | 3.024% / 2.417% | 1.974% / 0.984% |

Pressure p95 improves by 4.57x and radial-velocity odd p95 by 2.46x.  The
axial-velocity maximum alone is essentially unchanged, although its p95 and RMS
both decrease; this is consistent with a localized discontinuity crossing the
axis diagnostic and prevents claiming that every near-axis error is eliminated.

## Time-step control

To separate spatial refinement from the smaller L2 time step, L1 was repeated
to the same 200 us with its finest time step reduced from 0.2 to 0.05 us, exactly
matching L2.  The dominant p95 residual changed by only:

| pressure | density | temperature | radial velocity | Mach |
|---:|---:|---:|---:|---:|
| -1.57% | -0.48% | -0.70% | -0.29% | +0.38% |

The dominant ring was unchanged for every field.  Therefore the 3.4--5.2x
same-radius improvement under L2 is a spatial-resolution effect, not merely the
smaller time step.

## Recommendation on an axis taper

Do **not** add a static collar or ad-hoc radial taper to the production
operator on the basis of this result.  Earlier fixed collars moved the stripe
to the collar edge, while K=0 has no artificial transition and the present
error converges strongly with `h`.  For Run165, refine the near-axis/nozzle
region to about 0.53 mm or finer if this visual error matters.  A conservative,
face-shared sensor blend could be studied separately if a coarse L1 production
run is unavoidable, but it needs its own constant-state, MMS, shock and
conservation gates and is not justified as the default fix here.

## Evidence

- `mach_pressure_L1_L2_200us_common_scale.png`: synchronous Mach and pressure
  fields on identical physical/color scales.
- `near_axis_residual_L1_L2_200us.png`: first 8/16-ring residual maps at the
  same physical-radius extent.
- `radial_residual_vs_physical_radius_200us.png`: p95 residual versus physical
  radius.
- `first8_16_ring_metrics.csv`: every ring used in the physical-radius audit.
- `filament_resolution_metrics.csv`: peak position, amplitude and width.
- `axis_parity_resolution_metrics.csv`: even/odd max, p95 and RMS defects.
- `l1_time_step_control.csv`: matched-finest-time-step control.
- `analyze_refinement.py`: reproducible read-only post-processor.
- `uniform_l1_200us/plt00500`, `uniform_l2_200us/plt01000`, and
  `uniform_l1_matched_fine_dt_200us/plt02000`: source plotfiles.

No production source file was modified by this audit.
