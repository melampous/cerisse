# AFD pressure-compression joint shock sensor

## Purpose and implementation

The Cartesian non-IBM AFD operator previously allowed any variable in its
broad nonsmoothness test to select the characteristic high-order LLF branch.
This included density contacts and under-resolved tangential velocity waves.
The revised selection requires three conditions on the same six-point face
stencil:

1. the existing broad test classifies the stencil as nonsmooth;
2. the maximum adjacent relative pressure jump exceeds 0.01;
3. the maximum dimensionless compression exceeds 0.001.

For a face stencil, the two new measures are

\[
S_p=\max_m
\frac{|p_{m+1}-p_m|}
     {\max(p_m,p_{m+1},\epsilon)}
\]

and

\[
S_c=\max_m
\frac{\max(-\nabla\mathbin{\cdot}\boldsymbol{u}_m,0)\Delta_{\min}}
     {\max_{\mathcal N_m}c}.
\]

The selected high-order flux is therefore

```text
characteristic high-order LLF
    iff shock_llf && nonsmooth && S_p > 0.01 && S_c > 0.001
```

The divergence uses centred physical derivatives in every Cartesian
direction. The scale is the minimum cell spacing and the maximum local sound
speed over the derivative stencil. Pressure and compression are required on
the same face stencil, but their maxima are not required at the same cell.

The first-order cell-centred LLF safeguards are separate. They remain active
for invalid stencils, rarefaction pockets, inadmissible reconstructed states,
and failed HLLC fluxes. The result therefore shows that a contact or
tangential velocity wave cannot activate the sensor-controlled high-order LLF
by itself. It does not claim that no LLF safeguard can ever be used.
The broad nonsmoothness test also continues to control the AFD flux
correction. A contact or tangential velocity wave can therefore suppress that
correction even when it does not select the high-order LLF branch.

The new runtime options are

```text
cns.afd_shock_pressure_jump = 0.01
cns.afd_shock_compression   = 0.001
```

Both values must be finite and lie in `[0,1]`. The implementation used by all
authoritative reruns has

```text
src/rhs/Afd.h
6ebc88cedd833e251ac88b0adf514fe3bca224981983b98ea006738f46494fa7
```

The change is confined to the bulk Cartesian AFD operator. The legacy IBM
wrapper retains its previous marker-aware selection. `AfdIBM.h`,
`AfdIBMDrivers.h`, and `Riemann.h` remained byte-identical to the archived
pre-change snapshot. No cut-cell, cut-control, aggregate, or AMReX EB flow
mechanism was introduced or used by these tests.

## Controlled A/B matrix

Every pair uses one executable and changes only `cns.afd_shock_llf` between
zero and one. The reconstruction, grid, CFL number, time integrator, initial
condition, boundary conditions, and new thresholds are fixed within each
pair.

| Test | Resolution | Purpose |
|---|---:|---|
| Advected square contact | 64, 128, 256, 512 cells | Exclude density-only activation |
| Advected shear wave | 8, 12, 16, 24, 32 points per wavelength | Exclude tangential-velocity-only activation |
| Shu--Osher | 200, 400, 800 cells | Retain shock and entropy-wave behaviour |
| Stationary Mach 10 shock | 200 by 20 and 400 by 40 cells | Retain grid-aligned shock stability |
| Underexpanded round jet | 96 by 64 by 64 cells | Practical near-field and main-compression sensitivity |

## Results

### Advected contact

The initial broad sensor marks ten faces on every grid. All ten marks come
from density. The pressure gate marks zero faces and the compression gate
marks zero faces. The final joint fallback count is also zero.

The paired solution metrics are identical at every resolution. At 512 cells,
both variants give a density L1 error of 0.0091745996. Direct paired
conservative-state comparison gives a maximum difference of zero. Relative
mass, momentum, and energy conservation errors remain below
\(4.47\times10^{-13}\).

Authoritative archive:
`advected_contact/results/afd_joint_sensor_20260729_r2`.

### Advected shear wave

The finite-amplitude mode-four test at eight points per wavelength is the
strongest selectivity case. The broad sensor marks all 32 faces because of
tangential velocity curvature. The pressure, compression, and final joint
fallback counts are all zero.

The paired conservative states are identical at every tested resolution. At
eight points per wavelength, both variants retain an amplitude ratio of
0.91346895 and a perturbation kinetic-energy ratio of 0.83701335. This removes
the false LLF damping seen with the historical broad sensor.

Authoritative archive:
`advected_shear_wave/results/afd_joint_sensor_20260729_r2`.

### Shu--Osher

All six calculations remain finite and positive. The leading shock is at
\(x=2.4\) for both switch values on all three grids. At 800 cells, its
equivalent thickness is 2.30025 cells with the LLF branch disabled and
2.28197 cells with it enabled. The offline final snapshot contains six
joint-sensor candidate faces in the leading-shock window. The entropy-wave
window contains none.

At 800 cells, the density total variation in the entropy-wave interval is
8.66430 with the LLF branch disabled and 8.21378 with it enabled. The density
L1 difference from the shared fine-grid diagnostic reference decreases from
0.0125199 to 0.00728458. This reference is not an exact solution, so the
comparison is diagnostic rather than an order or accuracy proof.

Authoritative archive:
`shuosher/results/joint_sensor_20260729_r2`.

### Mach 10 grid-aligned shock

The unperturbed controls retain zero transverse front distortion for both
switch values. In the offline initial-state analysis, the seeded condition
satisfies the transverse joint sensor on all 20 rows of the 200 by 20 mesh and
all 40 rows of the 400 by 40 mesh.

| Grid | Final front distortion, LLF off | Final front distortion, LLF on | Suppression |
|---:|---:|---:|---:|
| 200 by 20 | 0.179084 cells | 0.0246506 cells | 7.26 times |
| 400 by 40 | 0.983977 cells | 0.0393209 cells | 25.0 times |

All density and pressure minima remain positive. The result retains the
intended suppression of grid-aligned odd-even shock-front decoupling.

Authoritative archive:
`afd_grid_aligned_shock_ab/results/afd_joint_sensor_20260729_r2`.

### Underexpanded round jet

The practical A/B test uses the canonical NPR 3 round-jet physical model in a
reduced three-dimensional Level 0 domain. The mesh has 16 cells across the
exit diameter. Both variants reach \(70\,\mu\mathrm{s}\) with finite positive
states. The normalized global cross-plane symmetry errors are below
\(3\times10^{-12}\), and the near-field odd-even measures are below
\(5\times10^{-16}\).

| Metric | LLF off | LLF on |
|---|---:|---:|
| Strongest centreline compression, \(x/D_e\) | 0.96875 | 0.96875 |
| Compression half-maximum width, \(D_e\) | 0.125 | 0.125 |
| Upstream centreline Mach number | 2.10131 | 1.89165 |
| Downstream centreline Mach number | 0.936756 | 0.893365 |
| Minimum pressure, Pa | 34241.8 | 37118.4 |
| Maximum Mach number | 2.20606 | 2.03401 |

The main compression location and measured width are unchanged, while the
coarse-grid amplitudes are sensitive to the LLF branch. The full-field
pressure L1 difference is \(8.07\times10^{-4}p_a\), and its maximum local
difference is \(0.116p_a\). This is useful stability and sensitivity evidence.
It is not a grid-converged Mach-disk validation.

Authoritative archive:
`underexpanded_round_jet_ab/results/joint_sensor_20260729_r2`.

## Decision

The revised sensor passes the intended regression matrix. Density
discontinuities and tangential velocity curvature no longer activate the
sensor-controlled high-order LLF in the two selectivity tests. Strong
compression still satisfies the joint sensor in Shu--Osher and Mach 10. The
grid-aligned shock stabilization is retained. The reduced jet calculation
shows a stable main compression location but non-negligible coarse-grid
amplitude sensitivity.

The current evidence supports retaining `cns.afd_shock_llf=1` with the default
pressure and compression thresholds for Cartesian AFD strong-shock
calculations. A production jet study should still include its normal grid and
AMR sensitivity checks. Highly anisotropic meshes also require a separate
sensor sensitivity study because the current compression scale uses
\(\Delta_{\min}\).
