# Full axisymmetric R-Z Navier--Stokes MMS

This is a standalone, non-IBM verification case. It exercises the production
point-sampled R-Z configuration. It does not use the annular-average Euler or
viscous options.

The default build is the full Navier--Stokes operator. An Euler-only build
with the same manufactured fields is also available through
`RZ_MMS_PHYSICS=euler`. Comparing the two spatial residual errors isolates the
viscous truncation error from the inviscid R-Z error. The paired comparison is
implemented in `analyze_operator_difference.py`.

The inviscid operator defaults to `RZ_MMS_EULER_SCHEME=llf-wenoz5`.
`llf-teno5`, `llf-teno6`, `afd-hllc-wenoz5`, and `afd-hllc-teno5` select the
corresponding production WENO/TENO operator. `skew-central-o4`,
`skew-jst-o4`, and `skew-central-o6` select the metric-consistent R-Z Skew
variants while retaining the same point-sampled manufactured solution. Every
build suffix and the test-level runtime manifest identify the selected scheme.

The smooth manufactured fields have the required parity at the axis. Density,
temperature, pressure, axial velocity, and total energy are even in the radial
coordinate. Radial velocity and radial momentum are odd. The lower radial
boundary therefore uses the solver's symmetry condition. The axial direction
is periodic. The outer radial ghost cells use the exact manufactured state.

The manufactured source matches the solver equation

```text
dU/dt = -div_RZ(Fc) + div_RZ(Fv)
        + [0, p/r - tau_theta_theta/r, 0, 0] + S.
```

It is therefore constructed with the sign convention

```text
S = dU_exact/dt + div_RZ(Fc) - div_RZ(Fv)
    - [0, p/r - tau_theta_theta/r, 0, 0].
```

This convention follows the implemented operator. The face array stores
`Fc-Fv`. `compute_rhs.cpp` applies its negative metric divergence, adds
`p/r`, calls `rz_geometric_source` to subtract
`tau_theta_theta/r`, and then adds the problem source. The
`order_rk=0` branch in `advance.cpp` copies this assembled RHS into the new
state. Its plotfile therefore contains `dU/dt`, not an evolved conserved
state.

Constant viscosity and conductivity isolate the differential operator. The
fields activate radial and axial convection, the pressure metric source,
normal and shear stresses, the hoop stress, viscous work, and radial and axial
Fourier heat conduction.

The radial momentum source is evaluated in the regular combined form

```text
(Fc_r - p - Fv_r + tau_theta_theta)/r.
```

This avoids subtracting two separately singular pressure or normal-stress
terms near the axis.

The default run performs two checks on 16, 32, 64, and 128 square meshes.

1. A one-step `order_rk=0` check compares the assembled spatial RHS with the
   exact time derivative at time zero.
2. An SSPRK3 evolution reaches time 0.01. The time step scales with the square
   of the mesh spacing so that time error is below the second-order viscous
   truncation error.
3. A second N=64 evolution halves the baseline time step at the same final
   time.
4. A paired Euler-only RHS matrix separates the viscous operator error from
   the inviscid R-Z error.

Both checks report cylindrical-volume-weighted L1 and L2 errors. They also
report global Linf errors, interior norms that exclude the outer three radial
rings, and separate Linf errors in the first and first three radial rings. The
interior norms prevent the exact outer-boundary closure from masking the
interior operator rate. The local axis norms are required because an axis
error can be hidden by the small physical volume of the first ring.

The expected acceptance result is second-order convergence in the
cylindrical-volume-weighted global L1 and L2 norms of the full
Navier--Stokes RHS and evolved state. The default point-valued R-Z
discretisation is not uniformly second order at the axis. In particular, the
first-ring Linf error of radial momentum is expected to approach first order.
This local result is reported as a method limitation. It is not hidden by the
small volume of the first annular ring. Global Linf and first-ring rates must
therefore not be presented as a uniform second-order result.

Run with

```bash
CERISSE_RZ_NS_MMS_RESULTS_ROOT=results/name ./run_study.sh
```

The launcher defaults to one MPI rank and records
`mpirun --oversubscribe` in every command file. The result archive contains
the source snapshot and hashes, both executable hashes, the Git description
and status snapshot, all explicit R-Z switches, raw logs, plotfiles, norm
tables, the paired viscous-operator analysis, and the N=64 time-step audit.

The build explicitly disables IBM, embedded boundaries, PelePhysics, LES, and
all optional R-Z shock or near-axis filters.

## Fixed AMR interface variant

The problem also provides an optional static Level-1 radial refinement band.
It is disabled by default, so the established uniform study is unchanged.
Enable it with

```text
prob.static_radial_refinement=1
prob.refine_radial_lo=0.25
prob.refine_radial_hi=0.75
```

The band spans the complete periodic axial direction. Its two radial edges
therefore isolate the R-Z coarse-fine treatment. Use a refinement ratio of two
and a blocking factor of four so that both requested edges remain aligned on
the grid sequence used by this case. The paired radial-pressure operator is
automatic for the all-fluid R-Z LLF-WENO-Z5 configuration.

Use `cns.amr_state_interp=linear` for this R-Z test. The production mapper
starts from AMReX's metric-aware linear conservative interpolation. If a child
leaves the configured closure-admissible set, including an active temperature
floor, it scales every conserved-variable deviation in that parent block by
one common coefficient, preserving the annular-volume weighted parent average.
The optional `conservative_quartic` interpolator uses
Cartesian equal-volume coefficients and is therefore not a valid basis for an
R-Z high-order AMR claim. With the linear interpolator, the uniform portions
of the spatial RHS retain the fifth-order LLF-WENO rate. The coarse-fine
interface has a lower order because the ghost-state interpolation is lower
order. This limitation must be reported separately from the R-Z face-flux
result.

`analyze_amr.py` evaluates the analytic solution on the active AMR leaf cells.
It reports cylindrical-volume-weighted global errors, errors in a configurable
band around both coarse-fine interfaces, errors within each level away from
the interfaces, first-axis-ring errors, and composite integral errors. For
example,

```bash
python3 analyze_amr.py results/amr/N16 results/amr/N32 \
  results/amr/N64 results/amr/N128 \
  --quantity state --physics euler \
  --refine-radial-lo 0.25 --refine-radial-hi 0.75 \
  --csv results/amr/errors.csv
```
