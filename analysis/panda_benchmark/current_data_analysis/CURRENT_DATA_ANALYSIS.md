# Current Mj=1.42 Data Analysis

Date: 2026-07-16

## Scope and provenance

This audit uses only the local NPZ files in `analysis/panda_benchmark`. No AWS
instance or simulation was accessed. The calculations can be reproduced with:

```bash
python3 analysis/panda_benchmark/analyze_current_data.py
```

The primary axial experiment is the Panda `m142den` survey at exactly
`Mj=1.42`. Its 37 phase bins are averaged with equal weight here, so the result
is a phase-cycle mean. The radial comparison uses Panda figure 5a at
`Mj=1.44`, because that is the available 30-station radial dataset. Simulation
density is normalized by `rho_j=1.6413 kg/m3`.

The records are not equivalent in duration:

| Dataset | Statistics interval | Length |
|---|---:|---:|
| TENO5 2D | 7.200-8.089 ms | 0.889 ms |
| WENO-Z5 short block | 7.240-8.089 ms | 0.849 ms |
| WENO-Z5 2D comparison mean | 6.500-10.053 ms | 3.553 ms |
| WENO-Z5 2D final mean | 6.500-15.747 ms | 9.247 ms |
| WENO-Z5 3D mean | 6.500-10.160 ms | 3.660 ms |

At 5.4 kHz, the TENO record contains only about 4.8 candidate screech periods.
Its start differs from the reconstructed WENO block by 0.040 ms, about 0.22 of
one 5.4 kHz period. The short comparison is therefore close in duration but is
not an exact scheme-only experiment.

## Main findings

1. The first shock cell is robust. Scheme and dimensionality changes have much
   smaller effects there than they have after the first closure.
2. The present TENO record does not show a systematic advantage over WENO that
   exceeds finite-window variability. It is too short to establish scheme
   equivalence or a stable downstream shock spacing.
3. The long-time 2D mean loses the second-cell amplitude and no longer resolves
   the third compression maximum. This is a real validation deficiency of the
   mean field, not something that should be dismissed as a bad averaging
   procedure.
4. The 3D calculation restores cells 2-4 and matches their amplitude much
   better, but its compression maxima drift progressively upstream and its
   fourth cell decays too slowly.
5. The persistent 2D near-axis density deficit has an entropy-like mean-field
   signature. It is strongly reduced in 3D downstream, but its numerical and
   physical causes are not yet isolated.

## TENO5 versus WENO-Z5

The direct short-record differences are smaller than WENO's own block-to-block
variation:

| Metric | TENO-WENO | Neighboring WENO blocks |
|---|---:|---:|
| Axial RMS, `0.9 <= x/D <= 6` | 0.0314 | 0.0736 |
| Median radial-profile RMS, 8 stations | 0.0085 | 0.0241 |

At every radial station, the TENO-WENO difference is below the change between
the first two adjacent WENO blocks. This supports the limited statement:

> No systematic TENO5 improvement larger than the present sampling uncertainty
> is detected.

It does not support either "TENO and WENO are equivalent" or "scheme effects
are negligible." The apparent second compression maximum is at `x/D=2.300`
for TENO and `2.425` for the reconstructed WENO block, while adjacent WENO
blocks alone span approximately `2.27-2.49`. The short-record extrema are phase
sensitive.

The WENO cumulative mean over 3.553 ms still differs from the 9.247 ms final
mean by an axial RMS of 0.0241 over `0.6 <= x/D <= 5`. Thus the 3.553 ms 2D
mean is useful but not fully converged at the precision needed to distinguish
WENO from TENO.

## Averaging-window interpretation

As the WENO average grows, the second compression maximum remains near
`x/D=2.40-2.44`, still about `0.10-0.14D` upstream of the same-Mach experiment.
Its preceding trough rises from about 0.61 in a short block to 0.74 in the final
mean, and the third compression maximum becomes unresolved.

This should not be called a "long-window artifact." The validation target is a
statistically converged mean, and the experimental profiles are also averaged.
The result instead shows that the simulated 2D mean contains much less
downstream shock-cell contrast than the experiment. Excessive axisymmetric
shock motion is a plausible explanation, but proving it requires an
instantaneous shock-position distribution or phase-conditioned analysis.

## Two-dimensional versus three-dimensional results

Compression maxima and preceding peak-to-trough amplitudes are:

| Cell | Experiment x/D | 2D x/D | 3D x/D | Experiment amp. | 2D/exp. | 3D/exp. |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.175 | 1.120 | 1.185 | 0.931 | 1.00 | 1.13 |
| 2 | 2.535 | 2.425 | 2.475 | 0.804 | 0.53 | 1.02 |
| 3 | 3.730 | 3.650 | 3.585 | 0.675 | 0.16 | 1.01 |
| 4 | 4.850 | unresolved | 4.580 | 0.407 | unresolved | 1.36 |

The 3D peak-position error evolves from `+0.01D` in cell 1 to `-0.06D`,
`-0.145D`, and `-0.27D` in cells 2-4. The physical interpretation is therefore
not "3D fixes 2D." It is:

> Three-dimensional freedom prevents the rapid collapse of the downstream
> mean shock cells, but the current 3D configuration shortens the downstream
> spacing and underpredicts shock-cell decay.

Raw axis RMS against the same-Mach experiment is 0.183 for 2D and 0.193 for
3D over `0.6 <= x/D <= 5`. This scalar does not rank the physics correctly:
the smooth 2D profile can score better after its shock cells disappear, while
the 3D profile is penalized twice for an oscillatory phase error.

The 30 radial stations tell the same qualified story:

| Region | 2D mean RMS | 3D mean RMS |
|---|---:|---:|
| First cell, `x/D <= 1.25` | 0.0469 | 0.0406 |
| Closure transition, `1.4-2.3` | 0.1073 | 0.0957 |
| Cells 2-4, `2.45-4.0` | 0.1251 | 0.1490 |
| Downstream, `x/D > 4` | 0.0668 | 0.1123 |

These values are nominal-station errors and remain sensitive to axial phase.
They must be reported with peak positions and amplitudes, not used alone to
claim that 2D or 3D is superior.

The 3D comparison also contains two extraction limitations. Its axial curve is
a four-cell near-axis proxy because a cell-centred Cartesian mesh has no cell
exactly on the axis. Its radial curves average only four Cartesian arms, not a
complete circular shell. In addition, native AMR spacing varies downstream:
the plotted repeated sampling interval must not be described as uniform
`D/256` resolution.

## Near-axis density deficit

For the 2D final mean, the normalized deficit
`rho(r/D=0.05)-rho(axis)` reaches 0.082 after the first closure and has a median
of 0.023 for `x/D >= 3.05`. The corresponding 3D downstream median is only
0.0007.

In 2D, the axis pressure is generally slightly higher than at `r/D=0.05`, while
the ratio-of-means temperature proxy `pbar/(R rhobar)` is 7-17 K higher after
the first closure. The mean-field entropy proxy

```text
ln(pbar_axis/pbar_0.05) - gamma ln(rhobar_axis/rhobar_0.05)
```

is positive and reaches 0.093 after the first closure. The data therefore
support the descriptive term "axis-centred entropy-like low-density core."

They do not yet prove the proposed mechanism. Axisymmetric shock focusing and
shock-capturing dissipation are plausible, but RZ-axis discretization and
nonlinear averaging remain alternatives. Also, `pbar/(R rhobar)` is not the
same quantity as directly averaged temperature. A causal claim requires
instantaneous entropy or total-pressure fields, an axis-resolution study, and a
matched RZ/Cartesian comparison.

## Corrections to earlier interpretations

- The nozzle-exit shoulder near `x/D=0.4` is not counted as a shock-cell peak.
- A short-block peak spacing is not a converged mean spacing.
- A longer average does not itself introduce a numerical bias; it exposes the
  simulated mean after shock motion has been averaged.
- Flat RMS versus disc radius only shows weak spatial sampling sensitivity near
  the axis. It does not demonstrate temporal convergence.
- The 3D four-arm extraction is not a true azimuthal average.
- The Mj=1.42 topology is a conical/intercepting-shock closure, not a clearly
  formed normal Mach disk.

## Defensible paper-level statement

The current data support one compact physical story:

> The first-cell gas-dynamic core is robust to reconstruction and
> dimensionality, whereas the downstream wave train is controlled by
> unsteady shock motion and three-dimensional mixing. Axisymmetric WENO and
> TENO give indistinguishable first-cell validation within the present sampling
> uncertainty, but the 2D time mean rapidly loses downstream shock-cell
> contrast. Three-dimensional freedom restores the downstream cells, while the
> current fixed-grid L5 calculation still predicts an upstream phase drift and
> insufficient decay. It is therefore suitable for topology and mode
> exploration, but not yet for quantitative screech amplitude or far-field
> acoustics.

## Highest-value next analyses

1. Track instantaneous compression-peak positions and build their probability
   distributions. This distinguishes true spacing error from mean smearing.
2. If scheme sensitivity remains a paper claim, extend TENO to at least the
   same 3.5 ms interval and start WENO/TENO statistics at exactly the same
   checkpoint time.
3. Compare the L5 3D result with a short-L6 or G2 shear-layer refinement case.
   Dimensionality and resolution are currently confounded.
4. Replace the four-arm radial extraction with full-shell azimuthal averages
   and report `m=0`, `m=1`, and Cartesian `m=4` content separately.
5. Extract instantaneous entropy and total pressure near the RZ axis and repeat
   at a second axis resolution before assigning a cause to the low-density core.

All numerical tables used above are stored as CSV files in this directory, and
the complete machine-readable summary is `summary_metrics.json`.
