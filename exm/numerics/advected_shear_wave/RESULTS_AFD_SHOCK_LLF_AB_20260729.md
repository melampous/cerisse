# AFD shock-LLF switch: advected-shear-wave results

> Historical baseline for the broad pre-gate sensor. The authoritative
> joint-sensor rerun is
> `results/afd_joint_sensor_20260729_r2` and is summarised in
> `../AFD_JOINT_SHOCK_SENSOR_20260729.md`.

The authoritative archive is
`results/afd_shock_llf_ab_20260729_r2`.  All twenty calculations and the
analysis completed successfully.  Its `comparison.csv` is identical to the
earlier successful `r1` result.  The `r2` archive is retained because its
analysis source snapshot also includes the initial-sensor diagnostic.

Every pair uses the same AFD HLLC WENO-Z5 executable, grid, CFL 0.05,
SSPRK(3,3), one-transit stop time and all other method settings.  The only
paired-run difference is `cns.afd_shock_llf`.

An automated command audit checked all ten pairs.  After removing the
required output-path token and the switch token, every paired command is
identical.  All twenty run status files report exit status zero.  The common
executable SHA256 is
`32471763b75ffc8ab21ea385ad4e93f64b5d84708b0835db20636875903af1ff`.

## Linear no-trigger control

The original shear wave has amplitude `1e-5 c0` and Fourier mode one.  The
sensor does not trigger at the initial or final state for 8, 12, 16, 24 or
32 points per wavelength.  The final conservative states from switch 0 and
switch 1 are bitwise identical on all five grids.

At 8 points per wavelength the amplitude and perturbation kinetic-energy
ratios are 0.978617 and 0.957692.  At 32 points per wavelength they are
0.999979 and 0.999957.

## Finite-amplitude short-wave stress test

The second exact solution uses amplitude `0.5 c0` and Fourier mode four.

| PPW | Sensor faces, initial/final | Amplitude, switch 0 | Amplitude, switch 1 | KE, switch 0 | KE, switch 1 |
|---:|---:|---:|---:|---:|---:|
| 8  | 32 / 32 | 0.913469 | 0.694468 | 0.837013 | 0.485248 |
| 12 | 0 / 0   | 0.988343 | 0.988343 | 0.977224 | 0.977224 |
| 16 | 0 / 0   | 0.997252 | 0.997252 | 0.994583 | 0.994583 |
| 24 | 0 / 0   | 0.999638 | 0.999638 | 0.999284 | 0.999284 |
| 32 | 0 / 0   | 0.999914 | 0.999914 | 0.999829 | 0.999829 |

At 8 points per wavelength the transverse-velocity part of the sensor flags
all 32 faces.  Enabling shock LLF reduces the retained amplitude by 24 per
cent relative to switch 0 and reduces the retained perturbation kinetic
energy by 42 per cent.  Its equivalent damping coefficient is 3.415e-4,
compared with 8.477e-5 for switch 0.

The trade-off is not purely dissipative.  The phase error is -0.667 rad for
switch 0 and 0.0100 rad for switch 1 at this severely under-resolved
resolution.  The LLF branch suppresses the large dispersive error but
increases the pressure L-infinity error from 5.95e-3 to 1.88e-2.  Neither
coarse result should be treated as resolved.

At 12 points per wavelength and above the sensor does not trigger.  The
paired conservative states are bitwise identical.  The abrupt change between
8 and 12 points per wavelength demonstrates that the current velocity-based
sensor can classify an under-resolved pure shear wave as a shock.  A
compression or pressure gate is therefore required if the fallback is to be
used as a shock-only stabilization.

Core source SHA256 values:

```text
Afd.h      5d461b134bed5e81f833b07b1b9fe6fb2c4267f79b37728801979bea8e817af9
AfdIBM.h   b1b3318561f19485847fa229432f25912166e2b9315b7ef00e2c71b5d09d1337
Riemann.h  d8a453fb9a886943453b6aae60368a5b6ea43ce7a7da172d987ae69593e3119d
```
