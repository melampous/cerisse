# Skew and Skew--JST verification, 28 July 2026

## Scope and provenance

This assessment adds only case-level compile selectors, run drivers, and
post-processing. It does not change `src/rhs/Skew.h` or another production
source file. The calculation uses the current worktree version of `Skew.h`.
Its SHA-256 is

```text
d08aee655c50d3b8ada706b250e07e572a4b8bfded450afe0bea60085cea2c09
```

The successful run is stored in
`results/y9000x_skew_validation_20260728_r2`. The host name is `Melampous`.
All twelve calculation exit codes are zero. All final densities and pressures
are positive. The source snapshot includes the active case files, `Skew.h`,
the worktree patch, build logs, run logs, exact replay commands, and hashes.

The two Skew executables are:

| Executable | Dissipation | \(C_2\) | \(C_4\) | SHA-256 |
|---|---:|---:|---:|---|
| `skew-central-o4` | off | 0 | 0 | `ed2b0788c2d115e33035cd75af804245949c775c596d803a75741dbe6f43a08e` |
| `skew-jst-o4` | on | 1.5 | 0.016 | `758a63581af61f463806d8ddf1b289e0b9a2cda5ba7c10df3d7688bd4cdce242` |

The first sandboxed launch is retained in
`results/y9000x_skew_validation_20260728`. It stopped before the first time
step because PMIx could not create a local listener. The executables in that
attempt have the same hashes as the successful run.

## Low-temporal-error diagonal vortex

The vortex travels through 0.1 domain lengths at a mean-flow angle of 45
degrees. The grids are \(40^2\), \(80^2\), \(160^2\), and \(320^2\).
The step counts are 64, 256, 1024, and 4096. Thus
\(\Delta t\) is proportional to \(\Delta x^2\). The final time is
\(5.759536524042908\times10^{-5}\) s on every grid.

### Undamped central operator

| \(N\) | Density \(L_2\) | Observed order | \(K^\prime/K^\prime_e\) | \(Z/Z_e\) | Advance time, s |
|---:|---:|---:|---:|---:|---:|
| 40  | \(3.4043\times10^{-5}\) | --   | 1.00004661 | 1.00022128 | 0.0226 |
| 80  | \(2.6721\times10^{-6}\) | 3.67 | 1.00000073 | 1.00001197 | 0.252 |
| 160 | \(3.0989\times10^{-7}\) | 3.11 | 0.99999957 | 1.00000042 | 3.29 |
| 320 | \(6.9132\times10^{-8}\) | 2.16 | 0.99999987 | 1.00000001 | 89.3 |

The finest-grid orders are 2.23, 2.16, and 2.47 for the density
\(L_1\), \(L_2\), and \(L_\infty\) errors. The pressure orders are 2.22,
2.10, and 2.22. The velocity orders lie between 2.51 and 2.59. The result
does not support a fourth-order accuracy claim for the current order-four
coefficient branch.

### Complete Skew--JST method

| \(N\) | Density \(L_2\) | Observed order | \(K^\prime/K^\prime_e\) | \(Z/Z_e\) | Advance time, s |
|---:|---:|---:|---:|---:|---:|
| 40  | \(3.9372\times10^{-5}\) | --   | 0.92160869 | 0.86158396 | 0.0294 |
| 80  | \(5.6516\times10^{-6}\) | 2.80 | 0.98756207 | 0.97585721 | 0.394 |
| 160 | \(7.5323\times10^{-7}\) | 2.91 | 0.99836501 | 0.99676597 | 4.46 |
| 320 | \(1.0833\times10^{-7}\) | 2.80 | 0.99979386 | 0.99959086 | 131.4 |

The velocity \(L_2\) orders reach 2.96 on the finest pair. The pressure
\(L_2\) order is 2.85. These data support an approximately third-order
observed rate for the complete method in this test. They do not support a
fourth-order claim.

## One-transit and five-transit dissipation

| Grid and duration | Density \(L_2\) | \(K^\prime/K^\prime_e\) | Peak \(\omega/\omega_e\) | \(Z/Z_e\) | Core circulation ratio |
|---|---:|---:|---:|---:|---:|
| \(40^2\), one transit | \(1.5970\times10^{-4}\) | 0.60339 | 0.51536 | 0.41224 | 0.69116 |
| \(80^2\), one transit | \(4.0141\times10^{-5}\) | 0.90054 | 0.83856 | 0.81891 | 0.94909 |
| \(160^2\), one transit | \(6.2309\times10^{-6}\) | 0.98508 | 0.96946 | 0.97084 | 0.99433 |
| \(80^2\), five transits | \(1.2817\times10^{-4}\) | 0.67640 | 0.54925 | 0.48896 | 0.76779 |

The equivalent numerical viscosity is evaluated from

\[
\nu_{\mathrm{eq},K}
=-\frac{R_v^2}{4t}
\log\left(\frac{K^\prime}{K^\prime_e}\right).
\]

This is a Gaussian-vortex, scale-specific dissipation indicator. It is not a
constant material property or a universal value for the scheme.

| Grid and duration | \(\nu_{\mathrm{eq},K}\), \(\mathrm{m^2\,s^{-1}}\) |
|---|---:|
| \(40^2\), one transit | \(5.4822\times10^{-3}\) |
| \(80^2\), one transit | \(1.1369\times10^{-3}\) |
| \(160^2\), one transit | \(1.6309\times10^{-4}\) |
| \(80^2\), five transits | \(8.4852\times10^{-4}\) |

The change between one and five transits confirms that this value depends on
the evolving resolved spectrum and the averaging interval.

## Same-source four-method comparison

The separate directory
`results/y9000x_skew_one_transit_comparison_20260728` contains a complete
same-source rerun. It uses the current `prob.h`, `inputs.uniform`, SSPRK(3,3),
and four MPI ranks for every method.

| Method, \(N=160\) | Density \(L_2\) | \(K^\prime/K^\prime_e\) | \(Z/Z_e\) | \(\nu_{\mathrm{eq},K}\), \(\mathrm{m^2\,s^{-1}}\) | Advance time, s |
|---|---:|---:|---:|---:|---:|
| LLF WENO-Z5 | \(1.7742\times10^{-6}\) | 0.999164 | 0.998038 | \(9.0739\times10^{-6}\) | 21.7 |
| LLF TENO5 | \(1.3033\times10^{-6}\) | 0.999235 | 0.998104 | \(8.3069\times10^{-6}\) | 23.0 |
| AFD HLLC WENO-Z5 | \(1.1405\times10^{-7}\) | 0.999840 | 0.999629 | \(1.7372\times10^{-6}\) | 18.7 |
| Skew--JST | \(6.2309\times10^{-6}\) | 0.985083 | 0.970837 | \(1.6309\times10^{-4}\) | 6.56 |

Skew--JST is faster in this small CPU test. It is also more dissipative at
the same grid size. The timings are indicative. They are not a formal
performance benchmark. The large one-transit AFD order between \(80^2\) and
\(160^2\) is pre-asymptotic and must not be described as a sixth-order method.

The comparison figure is
`one_transit_four_method_comparison.pdf` in the comparison result directory.
