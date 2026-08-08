# R-Z LES hoop-stress coefficient regression

This build reuses `../rz_viscous_axis/prob.h` with Smagorinsky viscosity
enabled. Run it with the parent `inputs` file and `prob.cell_average=0`.

For `prob.mode=0`, `u_r=a*r`, the exact molecular plus SGS radial viscous RHS
is zero. A nonzero result isolates a mismatch between the effective viscosity
used by the face stresses and the viscosity used by the cell-centred
`-tau_theta_theta/r` source.

## 2026-07-19 result

Before the coefficient repair, axis errors for Nr=16/32/64 were
`9.853e-5`, `4.927e-5`, and `2.463e-5`. After `viscous_t` was changed to use
the same molecular/Pele/SGS coefficients in the hoop source, the errors are
`1.39e-17`, `2.78e-17`, and `5.55e-17`.

This build tests `viscous_t`. The older `viscousLES_t` implementation still
recomputes the hoop source with molecular `mu` and `xi=0`; its combined
LES/ATF coefficient path remains unsupported in R-Z until it receives the same
refactor.
