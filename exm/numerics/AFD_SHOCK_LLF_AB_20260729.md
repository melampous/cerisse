# AFD shock-region LLF fallback A/B study

> Historical baseline. This report describes the broad pre-gate sensor.
> The implemented pressure-and-compression joint sensor and its repeated
> five-case assessment are documented in
> `AFD_JOINT_SHOCK_SENSOR_20260729.md`.

## Purpose

This study compares `cns.afd_shock_llf=0` and
`cns.afd_shock_llf=1` in the current AFD--HLLC--WENO-Z5
implementation. Each A/B pair uses one executable and the same grid, CFL
number, time integrator, final time, initial condition, boundary conditions,
AFD correction, reconstruction, and smoothness threshold. The paired command
changes only `cns.afd_shock_llf` and the required output path.

When the option is zero, an admissible non-rarefaction stencil proceeds to the
AFD--HLLC branch. When the option is one, a stencil classified as non-smooth
uses characteristic high-order LLF reconstruction. The independent
cell-centred LLF safeguard for invalid or rarefaction stencils is active in
both variants.

The solver sources used by all four studies had the following SHA256 values.

```text
Afd.h      5d461b134bed5e81f833b07b1b9fe6fb2c4267f79b37728801979bea8e817af9
AfdIBM.h   b1b3318561f19485847fa229432f25912166e2c71b5d09d1337
Riemann.h  d8a453fb9a886943453b6aae60368a5b6ea43ce7a7da172d987ae69593e3119d
```

The source hashes remained unchanged after all runs. No production solver or
IBM source was changed for this study.

## Test matrix

| Test | Resolution | CFL | Final time | Main purpose |
|---|---:|---:|---:|---|
| Shu--Osher | 200, 400, 800 cells | 0.30 | 1.8 | Shock and entropy-wave resolution |
| Advected square contact | 64, 128, 256, 512 cells | 0.20 | One transit | Contact smearing and conservation |
| Advected shear wave | 8, 12, 16, 24, 32 points per wavelength | 0.05 | One transit | Shear-wave dissipation and sensor selectivity |
| Stationary Mach 10 shock | 200 by 20 and 400 by 40 cells | 0.25 | 0.05 | Grid-aligned odd-even shock stability |

All runs use SSPRK(3,3), a uniform Cartesian Level 0 grid, the AFD correction,
and WENO-Z5 point reconstruction.

## Results

### Shu--Osher

Both variants completed on all three grids with positive density and pressure.
The measured leading-shock position was \(x=2.4\) for both settings. At 800
cells, the equivalent shock thickness was 2.300 cells with the fallback off
and 2.282 cells with it on. The fallback therefore did not move or broaden
the leading shock in this test.

In the post-shock entropy-wave interval, the fallback reduced density total
variation from 8.6643 to 8.2138 at 800 cells. It reduced the peak-to-peak
density range from 1.6812 to 1.5726. This is about five to six per cent
additional damping. Error relative to a shared 6400-cell Skew--JST diagnostic
profile decreased when the fallback was enabled. That profile is not an exact
solution, so this comparison is not an accuracy proof.

Detailed result:
[`shuosher/results/shock_llf_ab_20260729_r2/SUMMARY.md`](shuosher/results/shock_llf_ab_20260729_r2/SUMMARY.md).

### Advected contact

The initial discontinuity activated the density part of the sensor on ten
faces. Pressure and velocity were constant. The final offline snapshot did
not activate the sensor. This confirms that the current non-smoothness test
does not distinguish a contact from a compression shock.

The practical difference was small in this case. Across 64 to 512 cells, the
ratio of the density \(L_1\) error with the fallback on to that with it off was
0.9971 to 0.9986. The one-percent mixed-cell count was identical in each pair.
Mass, momentum, and energy conservation errors remained below
\(4.5\times10^{-13}\).

Detailed result:
[`advected_contact/RESULTS_AFD_SHOCK_LLF_AB_20260729.md`](advected_contact/RESULTS_AFD_SHOCK_LLF_AB_20260729.md).

### Advected shear wave

The small-amplitude mode-one control did not activate the sensor. The final
conservative states from the two switch values were bitwise identical at
every tested resolution.

The finite-amplitude mode-four stress test exposed a selectivity problem. At
eight points per wavelength, the transverse-velocity part of the sensor
activated all 32 faces. The retained amplitude decreased from 0.9135 with the
fallback off to 0.6945 with it on. The retained perturbation kinetic energy
decreased from 0.8370 to 0.4852. The equivalent damping increased by a factor
of about 4.03. The LLF branch also reduced a large phase error at this
severely under-resolved resolution. This result is therefore a dissipation and
dispersion trade-off, not a resolved-wave accuracy comparison.

At 12 or more points per wavelength, the sensor did not activate and every A/B
pair was bitwise identical.

Detailed result:
[`advected_shear_wave/RESULTS_AFD_SHOCK_LLF_AB_20260729.md`](advected_shear_wave/RESULTS_AFD_SHOCK_LLF_AB_20260729.md).

### Mach 10 grid-aligned shock

The unperturbed 200 by 20 controls retained a planar front for both settings.
The seeded cases began from the same alternating subcell shock displacement.

| Grid | Final front distortion, fallback off | Final front distortion, fallback on | Reduction factor |
|---:|---:|---:|---:|
| 200 by 20 | 0.1791 cells | 0.02465 cells | 7.26 |
| 400 by 40 | 0.9840 cells | 0.03932 cells | 25.0 |

Without the fallback, the physical front distortion increased under
refinement. With the fallback, it decreased. At 400 by 40 cells, the
row-to-row pressure range was 51.15 with the fallback off and 2.168 with it
on. The transverse velocity remained below \(4\times10^{-17}\). The observed
mode is therefore row-wise odd-even shock-front decoupling rather than
transverse-velocity growth.

Detailed result:
[`afd_grid_aligned_shock_ab/RESULTS.md`](afd_grid_aligned_shock_ab/RESULTS.md).

## Decision

The present evidence does not support setting `cns.afd_shock_llf=0` as the
default for strong-shock calculations. Pure AFD--HLLC strongly amplifies the
seeded grid-aligned Mach 10 shock perturbation, and the amplification becomes
larger under refinement. The characteristic high-order LLF fallback suppresses
that mode without moving or broadening the Shu--Osher leading shock.

The current sensor is broader than a shock sensor. It can activate on a pure
contact and on a severely under-resolved finite-amplitude shear wave. The
fallback should therefore remain enabled for the present strong-shock
calculations, but the next implementation should gate it with compression and
pressure information. The revised gate must retain the Mach 10 stability
result while avoiding activation on the contact and shear tests.

Before changing a thesis production configuration, the revised method should
also be compared in a representative underexpanded-jet or Mach-disk case. The
four tests here isolate numerical behaviour. They do not replace a
target-flow sensitivity study.
