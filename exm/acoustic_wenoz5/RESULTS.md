# WENO-Z5 acoustic validation results

All tests use WENO-Z5, RK3, CFL = 0.30 and an acoustic perturbation small
enough to remain in the linear regime.  The periodic-wave and AMR tests use
`rho0 = p0 = 1`, `gamma = 1.4` and `amplitude = 1e-5`.

## Periodic plane wave

The domain contains two wavelengths.  Each case was advanced for 20 acoustic
periods and compared with the exact translated wave using its Fourier mode.

| PPW | amplitude loss | phase error | relative L2 error |
|---:|---:|---:|---:|
| 16 | 2.933% | +0.1636 deg | 2.947e-2 |
| 24 | 0.520% | +0.00823 deg | 5.251e-3 |

The 24-PPW result is clean.  The 16-PPW result is usable for identifying a
dominant frequency, but its accumulated damping after 20 periods is already
visible.

## Oblique outgoing acoustic disturbance

The native two-dimensional boundary test uses WENO-Z5 and a Mach-0.3 mean
flow.  A short domain is compared with a longer domain at the same grid
spacing.  The quoted metric is the square root of error energy divided by the
initial disturbance energy, not a literal reflected-wave coefficient.

At N = 96, for the clean 30-degree Gaussian-windowed acoustic wave train:

| outlet treatment | normalized boundary error |
|---|---:|
| characteristic | 2.027% |
| first-order extrapolation | 1.828% |
| LODI | 0.462% |

For the non-carrier oblique Gaussian blob, the corresponding values are
9.187%, 3.974% and 3.591%.  About 47% of that blob remains in the reference
domain at the sample time, so this second result is dominated by its slow and
non-radiating content and should not be interpreted as a reflection
coefficient.

The N = 48 to N = 96 wave-train metric is not monotone.  LODI is best at the
higher resolution, but this test does not establish asymptotic convergence of
the outlet treatment.

## AMR interface crossing

A Gaussian-windowed carrier with 8 wavelengths per unit length crosses both
edges of a fixed level-1 patch.  The base grid has 16 PPW, the patch has 32
PPW, and the comparison solution is a globally uniform 32-PPW grid.  The
right- and left-running acoustic characteristics are
`p+ = (p' + rho0*c0*u')/2` and `p- = (p' - rho0*c0*u')/2`.

| reflux | AMR-only reflected L2 ratio | transmitted carrier amplitude loss | carrier phase error | transmitted waveform L2 error |
|---|---:|---:|---:|---:|
| off | 4.55e-7 | 1.13% | -0.353 deg | 9.90% |
| on | 5.00e-7 | 1.27% | +0.055 deg | 9.90% |

The raw left-running content is about `3e-6` in both AMR and uniform-fine
runs; subtracting the fine-grid reference isolates an AMR-interface component
of only about `5e-7`.  Thus no meaningful spurious reflection is detected in
this case.  The important AMR effect is instead the dispersion/dissipation
change along the mixed 16/32-PPW path: roughly 1.3% carrier-amplitude loss and
9.9% full-waveform L2 difference relative to an everywhere-fine grid.  Reflux
has little effect on reflection, although it reduces the carrier phase error.

## Practical conclusion

For a screech calculation, 24 PPW is a sensible minimum target for the
frequency-bearing acoustic mode with this setup.  Sixteen PPW can still locate
the spectral peak, but it should not be trusted for long propagation paths or
sound-pressure amplitude.  The AMR interface itself looks acoustically clean;
the resolution change across it remains dispersive, so microphones and the
dominant feedback path should preferably stay on a uniform refinement level.
