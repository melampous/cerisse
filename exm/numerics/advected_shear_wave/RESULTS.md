# Advected shear-wave results

The completed auditable archive is
`results/matrix_20260728_r2`.  The earlier `matrix_20260728` directory records
a sandbox MPI socket failure and is not a numerical result.

All four builds, all four smoke tests, all 28 full calculations and the
Fourier analysis completed successfully.  The primary matrix used CFL 0.05.
The values below are measured after one exact transit.

| Method | PPW | A(T)/A(0) | KE(T)/KE(0) | Phase-speed error | nu_num |
|---|---:|---:|---:|---:|---:|
| LLF WENO-Z5 | 8 | 0.90350840 | 0.81632913 | -1.9551e-3 | 1.5206e-3 |
| LLF WENO-Z5 | 24 | 0.99959705 | 0.99919426 | -2.4073e-6 | 6.0397e-6 |
| LLF WENO-Z5 | 48 | 0.99998784 | 0.99997568 | -3.6135e-8 | 1.8225e-7 |
| LLF TENO5 | 8 | 0.92274290 | 0.85145445 | -1.4868e-3 | 1.2049e-3 |
| LLF TENO5 | 24 | 0.99962026 | 0.99924066 | -2.2693e-6 | 5.6918e-6 |
| LLF TENO5 | 48 | 0.99998798 | 0.99997595 | -3.5814e-8 | 1.8020e-7 |
| AFD HLLC WENO-Z5 | 8 | 0.97861738 | 0.95769241 | -1.6877e-3 | 3.2391e-4 |
| AFD HLLC WENO-Z5 | 24 | 0.99991051 | 0.99982102 | -2.2760e-6 | 1.3411e-6 |
| AFD HLLC WENO-Z5 | 48 | 0.99999718 | 0.99999436 | -3.5815e-8 | 4.2273e-8 |
| Skew JST | 8 | 0.87654387 | 0.76832915 | -1.1785e-2 | 1.9746e-3 |
| Skew JST | 24 | 0.99466415 | 0.98935677 | -1.5531e-4 | 8.0175e-5 |
| Skew JST | 48 | 0.99932570 | 0.99865186 | -9.7667e-6 | 1.0108e-5 |

At 24 points per wavelength, halving the CFL from 0.05 to 0.025 changed the
amplitude ratio by about `1.9e-8` for every method.  The largest change in the
normalised phase speed was `1.12e-11`.  The comparison is therefore controlled
by spatial error at the reported precision.

For Skew-JST, the density and pressure sensors are zero in the exact solution.
The measured damping was compared with

```text
A(T)/A(0) =
  exp[-16 C4 lambda T sin(k dx/2)^4 / dx].
```

The maximum relative difference over the six grids was `5.59e-7`.  This small
difference contains SSPRK(3,3) time error, roundoff and any finite sensor
response.  The result confirms the implemented order-four JST damping
operator and the use of `lambda=abs(ux)+c0` in this all-fluid branch.

The Skew-JST phase error also agrees with the modified wave number of the
fourth-order central derivative,

```text
k_star dx = [8 sin(theta) - sin(2 theta)] / 6.
```

The largest absolute difference between the measured and analytic normalised
phase-speed errors was `6.65e-9`.

For the selected `C4=0.016`, Skew-JST damps this smooth mode more strongly than
the other three methods.  At 24 points per wavelength, its equivalent
coefficient is about 13 times the LLF WENO-Z5 value and 60 times the AFD HLLC
WENO-Z5 value.  This statement applies only to this Fourier mode, grid, CFL and
small-amplitude smooth state.

The reported

```text
nu_num = -log(A(T)/A(0)) / (k^2 T)
```

is a wavenumber-specific equivalent damping coefficient.  It is not a
material viscosity.  It depends on the grid, wave number, time integrator,
local state and nonlinear sensor response.  The present result does not
measure shock thickness, contact-wave damping, multidimensional coupling or
AMR-interface error.
