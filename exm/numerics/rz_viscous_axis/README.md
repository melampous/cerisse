# R-Z viscous-axis regression

This diffusion-only case stores the spatial RHS directly and checks four
smooth fields with constant viscosity:

- `mode=0`, `u_r=a*r`: radial viscous RHS is exactly zero.
- `mode=1`, `u_r=a*r+b*r^3`: radial RHS is `(32/3)*mu*b*r`.
- `mode=2`, `u_z=c*r^2`: axial RHS is `4*mu*c`.
- `mode=3`, `u_z=d*(z-z0)`: radial viscous RHS is exactly zero.

Modes 0 and 3 are direct cancellation tests between the metric face-flux
divergence and `-tau_theta_theta/r`.

`prob.cell_average=1` is the default and initializes cylindrical control-volume
averages. Set `prob.cell_average=0` and pass `--point-values` to `analyze.py`
for the point-sample control matrix.

## 2026-07-19 results

Two implementation defects were isolated:

1. `compute_rhs` added the Euler `+p/r` source even for `no_euler_t`. A pure
   diffusion test therefore produced `max(abs(RHS_r))=320`, exactly `p/r` in
   the first cell. `no_euler_t` now disables that source explicitly.
2. Face stresses used effective molecular plus SGS viscosity (and Pele bulk
   viscosity), while the hoop source used only `cls->visc(T)` and `xi=0`.
   `viscous_t` now rebuilds the same cell transport coefficients for the hoop
   source. The Smagorinsky control in `../rz_viscous_axis_les` falls from
   `9.85e-5` to roundoff at Nr=16.

After those fixes, point-sample controls give:

| field | result |
|---|---:|
| `u_r=a*r` | radial RHS <= 1.54e-15 through Nr=64 |
| `u_r=a*r+b*r^3` | axis error 1.67e-4, 8.33e-5, 4.17e-5 |
| `u_z=c*r^2` | axial error <= 9.48e-16 through Nr=64 |
| `u_z=d*(z-z0)` | both momentum errors <= 8.05e-16 through Nr=64 |

## Annular-average axis correction

The finite-volume-consistent cylindrical average exposed a third defect. For
`u_r=a*r`, the first-cell average is `(2/3)*a*dr`, whereas the point value at
the geometric midpoint is `(1/2)*a*dr`. The legacy viscous derivative and
`u_r/r` formulas consumed the former as the latter. The resulting axis radial
RHS was `3.793e-2`, `7.585e-2`, and `1.517e-1` on Nr=16/32/64, growing as
`1/dr`.

Set `cns.rz_viscous_annular_average=1` to use the corrected operator. It writes
the regular radial velocity as `u_r=r*w`, recovers `w` from the annular
centroid, and evaluates both face stresses and the hoop source with
`du_r/dx_j = delta_rj*w + r*dw/dx_j`. The same representation is therefore
used on both sides of the radial-stress/hoop-stress cancellation.

The cylindrical-average matrix with the correction enabled gives:

| field | Nr=16 | Nr=32 | Nr=64 |
|---|---:|---:|---:|
| `u_r=a*r`, radial Linf | 2.78e-16 | 7.21e-16 | 4.26e-15 |
| `u_r=a*r+b*r^3`, axis error | 2.70e-5 | 1.35e-5 | 6.75e-6 |
| `u_z=c*r^2`, axial error | 2.60e-17 | 3.33e-16 | 1.35e-15 |
| `u_z=d*(z-z0)`, max momentum error | 5.55e-17 | 1.94e-16 | 8.05e-16 |

The option defaults to off while the same-checkpoint SRP A/B test is evaluated.
It currently supports the second-order non-LES viscous operator; requesting it
with LES aborts instead of silently mixing inconsistent strain definitions.
