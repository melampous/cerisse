# Periodic advected shear wave

This small Euler test measures the dissipation and dispersion of four inviscid
methods without a shock or a physical viscosity.  It does not modify any
production solver source.

The exact state is

```text
rho = 1
p   = 1
ux  = 0.5 c0
uy  = 1e-5 c0 sin(k x)
k   = 2 pi / L
c0  = sqrt(1.4)
L   = 1
```

All calculations use a periodic, uniform, Level-0 Cartesian grid and
SSPRK(3,3).  The primary matrix uses CFL 0.05 and 8, 12, 16, 24, 32 and 48
points per wavelength.  Each method is repeated at 24 points per wavelength
with CFL 0.025.

The analysis reports the Fourier amplitude ratio, transverse perturbation
kinetic-energy ratio, phase-speed error and

```text
nu_num = -log(A(T)/A(0)) / (k^2 T).
```

`nu_num` is a wavenumber-specific equivalent damping coefficient.  It is not
a constant material viscosity and must not be extrapolated to shocks or
general nonlinear flow.

For the order-four Skew-JST method, uniform density and pressure make the
second-order JST sensor zero in exact arithmetic.  Its fourth-difference term
has the semi-discrete decay rate

```text
sigma4 = 16 C4 lambda sin(theta/2)^4 / dx,
lambda = abs(ux) + c0,
theta  = k dx.
```

The corresponding analytic amplitude ratio is `exp(-sigma4 T)`.  The study
compares this curve with the solver result.  The remaining difference includes
SSPRK time error, roundoff and any finite sensor response.

The analysis writes both a 300 dpi PNG and a vector PDF.  The PDF uses
embedded TrueType fonts.

Run the complete study with

```bash
./run_study.sh
```

Set `SHEAR_WAVE_RESULTS_ROOT` to select a different new output directory.
The script refuses to overwrite an existing archive.

## AFD shock-LLF switch comparison

Run

```bash
./run_afd_shock_llf_ab.sh
```

to compare `cns.afd_shock_llf=0` and `1` with one freshly checked
AFD HLLC WENO-Z5 executable.  The linear mode uses the original amplitude of
`1e-5 c0`.  It is a no-trigger control.  A second exact shear solution uses
amplitude `0.5 c0` and Fourier mode four.  This finite-amplitude short wave
checks whether the pressure and compression gates prevent transverse
velocity curvature alone from switching the flux to LLF.  Paired runs use
the same grid, CFL, stop time and all other method settings.  The analysis
reports the broad curvature sensor, pressure gate, compression gate, final
intersection and maximum conservative-state difference between final
plotfiles.  The three thresholds are 0.08, 0.01 and 0.001.
