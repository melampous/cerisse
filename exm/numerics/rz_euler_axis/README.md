# R-Z Euler-axis regression

This inviscid manufactured case checks the radial LLF-WENO operator with
smooth regular velocity fields. It fills either point samples or cylindrical
control-volume averages and writes the spatial RHS after one stage.

- `mode=0`: `rho=const`, `p=const`, `u_r=a*r`.
- `mode=1`: `rho=const`, `p=const`, `u_r=a*r+b*r^3`.

The analytic mass and radial-momentum right-hand sides are

```text
R_rho = -rho * (2*a + 4*b*r^2)
R_rhou = -rho * (3*a^2*r + 10*a*b*r^3 + 7*b^2*r^5).
```

For cell-average tests, the powers of `r` on the right-hand side are replaced
by their annular averages. This specifically exercises radial advection, which
is absent from the planar-shock and pressure-only R-Z regressions.

## Defects isolated by this case

For a regular odd field `u_r=a*r`, an R-Z finite-volume cell stores the annular
average

```text
<u_r>_i = a * (r_i + dr^2 / (12*r_i)).
```

Passing that value directly to a Cartesian midpoint WENO stencil is an O(dr)
error in the first ring and produces an O(1) continuity-RHS error at the axis.
The error therefore does not converge even though the parity ghost cells are
correct. The first correction recovers the regular midpoint radial velocity
before radial reconstruction.

There is a second, independent mismatch. Standard finite-difference WENO
reconstructs a *numerical flux*, not the physical point flux at the face. For a
quadratic radial-momentum flux its smooth optimal formula contains a constant
`-dr^2/12` face error. Adjacent Cartesian face differences cancel that term,
but an R-Z divergence first multiplies each face by a different radius, so it
does not cancel. The second correction uses half-cell point-interpolation
polynomials for radial WENO fluxes.

Enable both corrections with

```text
cns.rz_euler_annular_average = 1
cns.rz_euler_point_flux = 1
```

Both switches are off by default while SRP checkpoint A/B validation is in
progress. Point-flux reconstruction requires annular-average correction. The
corrections are radial-only and are not applied to an IBM stencil containing
reflected solid ghost states, because those states are already point values.

For the annular midpoint and radial point-flux treatment, `mode=0`, `Nr=16` gives
density axis `Linf=8.5321e-14` and radial-momentum axis
`Linf=1.4256e-12`.  For `mode=1`, the density global error is
`3.28517e-4, 8.21348e-5, 2.05338e-5, 5.13344e-6` for
`Nr=16,32,64,128` (order `2.000`), while radial momentum is
`3.37141e-4, 8.68877e-5, 2.20487e-5, 5.55337e-6` (orders
`1.956, 1.979, 1.989`).  These are the smooth-field Gate-B values for source
SHA-256 `611ba288940fdc19775d34bc0c3cac29b29f2880badac803d74edca550019bf6`.
Strong-shock behavior is covered separately by the planar-shock and Run165
gates.

## Reference results

For annular-average input and `mode=0`, the legacy density axis error is
`5.799e-2` at every tested resolution (`Nr=16..128`). The first correction
reduces density error to roundoff but leaves a first-order radial-momentum
error. With both corrections, both equations are exact to roundoff at `Nr=16`:

| radial treatment | density axis Linf | momentum axis Linf | `R_rhou/r` |
|---|---:|---:|---:|
| legacy | `5.799e-2` | `3.500e-2` | `0.1460` |
| annular midpoint only | `1.1e-13` | `4.167e-4` | `0.1466` |
| midpoint + point flux | `8.5e-14` | `1.4e-12` | `0.1600` |

For `mode=1`, both corrections give stable second-order global convergence:

| Nr | density global Linf | momentum global Linf | momentum axis Linf |
|---:|---:|---:|---:|
| 16 | `3.285e-4` | `3.371e-4` | `5.367e-6` |
| 32 | `8.213e-5` | `8.689e-5` | `2.304e-6` |
| 64 | `2.053e-5` | `2.205e-5` | `6.954e-7` |
| 128 | `5.133e-6` | `5.553e-6` | `1.887e-7` |

An axis-spanning Mach-10 planar shock remains radially uniform to roundoff.
A one-cell radial shock-position perturbation decays to `0.0448` cell by
`t=0.02`; this guards against replacing the smooth-axis fix with a carbuncle
instability.

The 2026-07-22 Y9000X matrix is archived at
`output/y9000x_matrix_20260722_r1`. After the default-off radial
conservative-state deconvolution code was added, the case was rebuilt and the
10-case serial matrix was repeated at
`output/local_default_off_regression_20260722_r5`. After removing only result
paths from the JSON reports, their complete numeric payloads are identical.
This is the non-regression evidence that the new opt-in path does not alter the
legacy/default or the existing annular-plus-point-flux configurations.

Build the case with `make -j`, run one step, then analyze its plotfile with:

```bash
python3 analyze.py output/plt00001 --mode 0
```
