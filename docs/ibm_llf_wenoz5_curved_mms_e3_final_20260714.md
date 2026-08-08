# LLF-WENO-Z5 curved-wall IBM MMS: e3 final audit

Date: 2026-07-14

## 1. Final result

The validated smooth-wall configuration is the established masked
GP/LLF-WENO-Z5 crossing and near-wall path with the IBM communication-halo
repair. The experimental wall-matched crossing and shifted same-side shell
must remain disabled.

```text
interp_order                  = 2
extrap_order                  = 3
interp_order_surf             = 2
extrap_order_surf             = 2
alpha / alpha_surf            = 1.0 / 1.0
ghost_layers                  = 1
support_visibility_mode       = 2
llf_ibm_face_local_crossing   = 0
llf_ibm_fluid_shell           = 0
llf_ibm_shifted_shell         = 0
llf_ibm_wall_matched_crossing = 0
```

For the source-free circular Euler vortex, the five-grid `t=0` RHS fit over
`N=128/160/192/224/256` gives the following first-fluid-layer L2 orders:

| Conservative component | fitted order |
|---|---:|
| rho | 2.175 |
| rho u | 1.964 |
| rho v | 1.964 |
| rho E | 2.176 |

At fixed time `t=0.002`, using RK(4,3) and `dt proportional to h`, the
`N=64/96/128/192` first-layer primitive L2 fits are:

| Primitive | fitted order |
|---|---:|
| rho | 2.89 |
| u | 2.63 |
| v | 2.63 |
| p | 2.91 |
| T | 2.83 |

The corresponding whole-fluid fits are `3.19, 3.00, 3.00, 3.23, 3.17`.
These numbers establish a stable second-order-or-better smooth-wall result for
this verification problem. They do not establish fifth-order global IBM
accuracy.

## 2. Exact problem and geometry protocol

The exact solution is the steady circular Euler vortex outside a circle of
radius `0.2`, centred at `(0.5,0.5)`:

```text
T       = 300 K,
u       = 50 e_theta m/s,
p(r)    = p_w (r/R)^a,
rho(r)  = p(r)/(R_g T),
a       = 50^2/(R_g T),
p_w     = 100000 Pa.
```

It satisfies `dp/dr = rho V^2/r`, so the wall uses the nonzero Euler-slip
curvature pressure compatibility rather than `dp/dn=0`. Exact states are
imposed at the outer boundary.

A fixed 65536-edge polygon is used at every fluid resolution. The volume-GP
normal is overridden by the analytic circle normal. This is intentional: it
isolates the fluid/IBM discretisation from faceted-geometry error. Therefore
the reported order proves the numerical closure in a geometry-exact limit; an
arbitrary fixed STL still needs a separate normal/location/curvature audit.

## 3. Partition-dependence defect and repair

Before this repair, `cls_t::NGHOST=3` was incorrectly used both for the WENO
stencil and for IBM image-point interpolation. With `eorder=3, alpha=1`, a
volume image-point support can reach about nine Cartesian cells from the
owning FAB. Records outside the three-cell halo were deactivated at FAB
boundaries. A `192^2` four-FAB run consequently differed from a one-FAB run in
28--32 cells near `i/j=127/128`; first-layer momentum RHS L2 increased from
`0.476` to `23.59`.

The repair separates the two concepts:

- the conservative RHS state retains the normal WENO halo;
- primitive and marker scratch data use geometry-derived volume and surface
  interpolation halos;
- valid primitive cells are exchanged before GP and surface reconstruction;
- periodic directions exchange the complete IBM halo;
- non-periodic physical boundaries remain limited to the region actually
  populated by FillPatch/BC data;
- the optional one-cell sharp-feature WLS translation is included in the
  halo upper bound.

After repair, all 648 volume image records and all 131072 surface image
records are active for `N=192`. One-FAB and four-FAB GPU RHS arrays are
bitwise identical. One-rank and two-rank CPU arrays are also bitwise identical
for `rho`, both momentum components, and `rhoE`.

## 4. Why the new crossing/shell prototype was rejected

The wall-matched crossing uses one reconstructed GP plus four same-side fluid
cells. The shifted shell extrapolates complete five-cell real-fluid WENO
blocks to near-wall fluid-fluid faces. Both improve or preserve the `t=0`
order, but the first near-wall cell then receives two extrapolatory high-order
fluxes without a proven stable boundary closure.

The fixed-time ablation found:

| Path | Failure/result |
|---|---|
| wall-matched crossing + shifted shell, `N=64` | negative internal energy at `t=0.000570` |
| shifted shell only, `N=64` | negative internal energy at `t=0.001575` |
| wall-matched crossing only, `N=192` | negative internal energy at `t=0.001497` |
| validated legacy crossing/shell, `N=192` | completed to `t=0.002` |

A fluid-side jump fallback was also tested. It preserved the initial RHS
bitwise and prevented the `N=192` crash, but did not suppress the unstable
mode: at `t=0.002` the whole-fluid pressure L2 error was about `405 Pa`, versus
`0.00423 Pa` for the validated legacy path. That fallback was removed rather
than accepted on a positivity-only criterion.

The experimental crossing/shell options remain default-off research paths.
They are not part of the recommended configuration.

## 5. Accuracy tables

### 5.1 Semi-discrete RHS, L2

Five-grid fitted orders:

| Band | rho | rho u | rho v | rho E |
|---|---:|---:|---:|---:|
| all fluid | 2.680 | 2.458 | 2.457 | 2.681 |
| interior excluding first layer | 2.238 | 2.348 | 2.346 | 2.226 |
| first fluid layer | 2.175 | 1.964 | 1.964 | 2.176 |

The first-layer pairwise rates remain grid-phase sensitive, especially over
`N=192 -> 224`, but the five-grid regression is approximately second order.

### 5.2 Fixed time, primitive L2

Four-grid fitted orders:

| Band | rho | u | v | p | T |
|---|---:|---:|---:|---:|---:|
| all fluid | 3.19 | 3.00 | 3.00 | 3.23 | 3.17 |
| first fluid layer | 2.89 | 2.63 | 2.63 | 2.91 | 2.83 |
| second layer | 2.50 | 2.77 | 2.77 | 2.75 | 2.25 |
| third layer | 2.34 | 2.76 | 2.76 | 2.66 | 2.57 |

At `N=192`, the whole-fluid pressure L2 error is `4.23e-3 Pa` and the first
layer pressure L2 error is `2.23e-2 Pa`.

## 6. Verification gates

- CUDA/NVCC `sm_89` build with cross-execution-space errors enabled: pass.
- CPU GCC build: pass.
- Final-source rebuild after removal of the rejected smoothness fallback: pass
  on both CPU and CUDA.
- Final-source `N=192` four-FAB GPU RHS smoke: all 648 volume images and all
  131072 surface images active; both `plt00000` and `plt00001` are bitwise
  identical to the converged baseline (matching SHA-256 data and metadata).
- Two-rank MPI CPU run: pass.
- One-FAB/four-FAB GPU equality: bitwise pass.
- One-rank/two-rank CPU equality: bitwise pass.
- 2-D exhaustive and sampled 3-D WLS polynomial reproduction: pass.
- Shifted WENO/TENO reconstruction unit test: pass as algebra, but its IBM
  runtime path remains rejected by the stability test above.
- Every reported volume/surface image record active: pass.
- `git diff --check`: pass.

## 7. Remaining limits

1. This result is Euler-slip and geometry-exact. The thermal/viscous NS gates
   remain separate.
2. A generic STL must demonstrate surface location, normal, and any required
   curvature accuracy independently.
3. The enlarged halo is fully validated across same-level FAB and MPI
   ownership. At AMR coarse-fine boundaries, IBM-tagged fine grids must retain
   an interpolation-width refinement buffer; a dedicated coarse-fine support
   availability audit is still required.
4. The current result supports stable second-order IBM accuracy, not a claim
   that the complete solid-domain calculation is globally fifth order.

## 8. Reproducibility artifacts

Generated tables are in `IBM/cases/validation/mms_ibm/`:

- `convergence_rhs_shell_llf_wenoz5_baseline_final_128_160_192_224_256_L2.md`
- `convergence_tfinal_llf_wenoz5_baseline_correctedhalo_64_96_128_192_L2.md`
- `error_tfinal_llf_wenoz5_baseline_correctedhalo_N192.md`

Plotfiles use the prefixes:

- `plot/solid_llf_wenoz5_baseline_final_G65536_e3_aframe_rhs_N*`
- `plot/solid_llf_wenoz5_baseline_correctedhalo_N*_t002`
