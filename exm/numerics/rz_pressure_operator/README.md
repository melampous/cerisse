# R-Z pressure-operator regression

This case stores the spatial RHS directly (`cns.order_rk=0`) and compares the
radial momentum component with `-dp/dr`. By default, both the initial pressure
and the exact RHS use the cylindrical `r dr dz` control-volume average.
Set `prob.cell_average=0` and pass `--point-values` to `analyze.py` for the
point-sample control matrix.

Profiles:

- `prob.profile=0`: smooth quadratic radial pressure.
- `prob.profile=1`: a smooth pressure jump with quadratic shock curvature.

The matrix is intended to compare the default `p_i/r_i` source against
`rz_wb_pressure_ncell`, `rz_pressure_quad_ncell`, and `rz_pressure_split` on
successively refined uniform R-Z grids.

## 2026-07-19 results

Quadratic pressure, cylindrical volume-average axis `Linf` error:

| treatment | Nr=16 | Nr=32 | Nr=64 | global Nr=64 |
|---|---:|---:|---:|---:|
| default | 2.37e-14 | 9.47e-15 | 1.89e-14 | 9.28e-14 |
| `wbp2` | 1.042e-2 | 5.208e-3 | 2.604e-3 | 2.604e-3 |
| `pquad2` | 3.472e-3 | 1.736e-3 | 8.681e-4 | 8.681e-4 |
| pressure split | 5.208e-3 | 2.604e-3 | 1.302e-3 | 4.340e-3 |

Smooth curved pressure layer (`width=0.04`), axis `Linf` error:

| treatment | Nr=16 | Nr=32 | Nr=64 |
|---|---:|---:|---:|
| default | 2.433e-4 | 3.235e-5 | 4.069e-6 |
| `wbp2` | 4.979e-2 | 2.573e-2 | 1.298e-2 |
| `pquad2` | 1.660e-2 | 8.576e-3 | 4.327e-3 |
| pressure split | 2.656e-2 | 1.301e-2 | 6.501e-3 |

For a curved layer kept one cell thick as the mesh is refined, the default
axis error decreases from `8.45e-5` to `2.05e-5`. The `wbp2` and `pquad2`
axis errors remain approximately constant at `2.54e-2` and `8.47e-3`;
pressure split develops a growing global error.

The accepted result is therefore the default pressure flux/source pair. The
three optional treatments are experimental failures and must not be enabled in
the SRP production case. A planar normal shock is covered separately by
`../rz_planar_shock` and is preserved to roundoff when no radial AMR interface
crosses the shock.
