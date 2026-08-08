# Current all-fluid R-Z pressure-axis gate

This is an all-fluid, two-dimensional R-Z regression for LLF-WENO/TENO,
AFD-HLLC-WENO/TENO and Skew4 central/JST. The variants exercise their current
paired metric-advection/ordinary-pressure operator on identical fields. It
does not compile GP-IBM, EB, cut-cell, cut-control or aggregate flow paths.

With `cns.order_rk=0`, `plt00001/Xmom` is the spatial radial-momentum RHS.
The case therefore tests the pressure operator directly, before time
integration can mix pressure, advection and kinetic-energy errors.

The point-value profiles and exact RHS values are

```text
constant:   p=p0                                  RHS=0
quadratic:  p=p0+A*r^2                            RHS=-2*A*r
quartic:    p=p0+A*r^2+B*r^4                      RHS=-2*A*r-4*B*r^3
curved:     p=p0+A*tanh((z-z0-k*r^2)/width)       RHS=2*A*k*r/width*sech^2(...)
```

Select the reconstruction at build and run time. Each build has a unique
executable name; the runner also checks a scheme-specific test manifest, so a
stale or mismatched binary is rejected. For example:

```bash
make -j4 RZ_AXIS_SCHEME=llf-wenoz5
./main2d.gnu.llf_wenoz5.ex inputs \
  amr.n_cell="64 256" \
  prob.profile=1 \
  amr.plot_file=results/curved_n64/plt
python3 analyze.py results/curved_n64/plt00001 \
  --profile 1 \
  --json results/curved_n64/metrics.json \
  --figure results/curved_n64/rhs_error.png
```

The two AFD-HLLC selectors are `afd-hllc-wenoz5` and `afd-hllc-teno5`; their
non-MPI executables are respectively `main2d.gnu.afd_hllc_wenoz5.ex` and
`main2d.gnu.afd_hllc_teno5.ex`.

The Skew selectors are `skew4-central` and `skew4-jst`.  The latter uses the
current Run165 coefficients `C2=1.5`, `C4=0.016`, and `sensor_power=1`.

`inputs_freestream` is the small 16 by 16 uniform R-Z gate requested for
stage-by-stage diagnosis.  Enable the default-off instrumentation with:

```bash
./<scheme-executable> inputs_freestream \
  cns.order_rk=3 cns.stages_rk=4 \
  cns.rz_stage_rhs_diagnostics=1 \
  cns.rz_stage_rhs_diagnostics_rings=4 \
  amr.plot_files_output=0 amr.checkpoint_files_output=0
```

The diagnostic reports the first four radial rings for every SSPRK(4,3) RHS
call and splits radial momentum into metric non-pressure flux, ordinary
pressure gradient, axial flux and closure.  See
`RZ_AXIS_DIAGNOSTIC_20260808.md` for the frozen results and interpretation.

For an LLF build, the run log must additionally contain

```text
rz_llf_weno_pressure_operator=paired_metric_advection_and_pressure_gradient_when_RZ
```

Run the complete four-grid matrix with:

```bash
./run_matrix.sh --scheme llf-wenoz5
./run_matrix.sh --scheme llf-teno5
./run_matrix.sh --scheme llf-teno6
./run_matrix.sh --scheme afd-hllc-wenoz5
./run_matrix.sh --scheme afd-hllc-teno5
```

This writes one log, JSON metric file, and diagnostic figure per case under a
new UTC-stamped result directory, plus a combined `summary.json`. The script
refuses to overwrite results and fails if the production scheme manifest is
mismatched, if an LLF build lacks its paired-pressure manifest, or if any
accuracy gate fails.

Do not enable `cns.rz_euler_annular_pressure_consistency` or
`cns.rz_euler_cell_average_deconvolution`; both are intentionally
incompatible with this paired operator.  Do not use the removed historical
`cns.rz_companion_pressure_flux` option.

The static case uses analytic user boundary data in both axial directions.
`prob.cell_average` is fail-closed at zero because the active all-fluid WENO
operator has POINT semantics, and `prob.gamma` is fail-closed at 1.4 to match
the compiled calorically-perfect-gas closure.  In the order-zero plots the
header time is the nominal driver time; the stored variables are the RHS
evaluated from the initial state, not a physical state advanced to that time.

The first matrix should use `Nr=16,32,64,128` for constant, quadratic and
quartic pressure, and `(Nr,Nz)=(16,64),(32,128),(64,256),(128,512)` for the
curved layer.  Acceptance requires finite output, a roundoff-scale constant
pressure RHS, and decreasing axis/first-four-ring errors for every smooth
nonconstant profile.  A dynamic radially uniform Mach-3/Mach-10 shock gate is
provided separately by `../rz_planar_shock` and must be rebuilt from the same
source revision before its results are combined with this gate.

## Constant-pressure density-contact gate

`inputs_contact` advances an exact stationary Euler contact for 100 SSPRK(4,3)
steps. Its default is a cylindrical density jump four cells from the R-Z axis;
set `prob.contact_location=0.015625` on the 64-cell mesh to put the jump at the
first interior radial face. Set `prob.contact_orientation=1` and
`prob.contact_location=0.5` for the orthogonal, radially uniform axial contact.
Both states have `p=2`, `u_r=u_z=0`, and density ratio 4:1. LLF may diffuse
density, but it must not generate pressure or velocity.

```bash
./main2d.gnu.llf_wenoz5.ex inputs_contact \
  amr.plot_file=results/contact_radial4/plt
python3 analyze_contact.py results/contact_radial4/plt00100 \
  --json results/contact_radial4/metrics.json
```

The gate remains on point-state semantics (`prob.cell_average=0`) and requires
the runtime manifest
`rz_llf_weno_pressure_operator=paired_metric_advection_and_pressure_gradient_when_RZ`.
Run the first-face, fourth-face and axial-plane matrix together with
`bash run_contact_matrix.sh --scheme llf-wenoz5` (or either TENO selector); it
refuses to overwrite an existing result tree
and returns nonzero unless pressure pollution and radial Mach remain below
`1e-10`. Density monotonicity is reported but deliberately is not an acceptance
condition for this pressure-balance gate.
