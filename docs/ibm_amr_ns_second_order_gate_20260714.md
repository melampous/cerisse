# IBM AMR Support and NS Wall-Flux Gate

Status date: 2026-07-14. Scope: fixed IBM, single-species ideal gas, Cartesian
meshes. The tested reconstruction is `iorder=2`, `extrap_order=3`,
`extrap_order_surf=3`, `ghost_layers=1`, and
`alpha=alpha_surf=0.6`.

## AMR coarse-fine support

Image-point values are reconstructed from same-level primitive states. A
coarse-fine ghost value inside a nominal quadratic WLS stencil changes the
moment conditions and is not an order-preserving fallback. The implementation
therefore has three fail-closed layers:

1. `ib.amr_support_buffer=1` tags a solid-fluid interface band. Its radius is
   derived from the next level's volume and surface interpolation halos and the
   refinement ratio. The band is formed by a radius-one transition mask and
   separable dilation, reducing setup work from an `O(r^D)` neighborhood scan
   per cell to `O(D r)` passes.
2. `ib.amr_support_audit=1` checks every nonzero cached GP and surface weight
   against the level's global valid `BoxArray`. Periodic supports are wrapped;
   physical-boundary supports are classified separately.
3. `ib.amr_support_strict=1` aborts during IBM initialization, before RHS, if
   any interior support exists only in a coarse-fine ghost region.

Runtime interpolation skips exact zero weights, so the audited and consumed
read sets are identical.

### AMR results

The 2-D positive test used a 64 by 64 coarse circle grid, refinement ratio two,
and a fine level covering 39.0625% of the domain. A coarse RK3 step and both
fine substeps completed on two CPU/MPI ranks and on CUDA `sm_89`.

| Target | Level | Targets | Active nonzero supports | Missing |
| --- | ---: | ---: | ---: | ---: |
| GP | 0 | 216 | 1,928 | 0 |
| GP | 1 | 432 | 3,856 | 0 |
| Surface | 1 | 12,288 | 110,592 | 0 |

The 3-D two-rank sphere test exercised the 27-point quadratic WLS path:

| Target | Level | Targets | Active nonzero supports | Missing |
| --- | ---: | ---: | ---: | ---: |
| GP | 0 | 1,176 | 31,736 | 0 |
| GP | 1 | 4,824 | 130,200 | 0 |
| Surface | 1 | 15,360 | 414,720 | 0 |

With the automatic buffer disabled, a deliberately narrow fine patch produced
176 affected GP targets and 1,248 missing nonzero supports. CPU/MPI and CUDA
both aborted during IBM initialization. This is the required negative control.

## Independent NS second-order gates

The tests separate semi-discrete reconstruction from time-evolved operators:

1. one `dt=1e-10` RK3 step tests direct surface reconstruction without
   interpreting RHS-output mode as a physical state;
2. heat-only fixed-time MMS isolates `diffusiveheat_t`;
3. viscous-only fixed-time MMS isolates `viscous_t`;
4. fully coupled NS is retained as a stronger stress test.

### Certified results

| Gate | Quantity | L2 order | Linf order |
| --- | --- | ---: | ---: |
| Direct surface | `dTdn` | 2.049 | 2.040 |
| Direct surface | normal stress | 2.077 | 1.842 |
| Direct surface | shear magnitude | 2.023 | 2.051 |
| Heat-only, fixed time | `dTdn` / `qwall` | 2.107 | 1.996 |
| Viscous-only, fixed time | normal stress | 2.128 | 2.014 |
| Viscous-only, fixed time | shear magnitude | 2.094 | 1.975 |

These sequences certify the independent surface differentiation, Fourier
heat-conduction, and viscous-stress paths at approximately second order. The
fully coupled first-fluid-cell state also passes: pressure is 2.962 order and
temperature is 2.333 order in L2.

### Coupled-NS release gap

The N=64/96/128/192 sequence at fixed `t=0.002` gives:

| Quantity | L2 order | Linf order | Status |
| --- | ---: | ---: | --- |
| Pressure | 3.247 | 3.191 | pass |
| `dTdn` / `qwall` | 1.624 | 1.606 | open |
| Normal stress | 2.242 | 2.171 | pass |
| Shear magnitude | 1.440 | 1.587 | open |

The fine sequence rules out a simple coarse-grid transient. It localizes the
remaining defect to differentiation of the coupled, evolved temperature and
tangential-velocity state; it does not indicate a general Cartesian viscous
tensor, Fourier-law, normal, or direct WLS reproduction failure.

The next numerical prototype should be a wall-location-consistent direct
thermal/viscous flux closure. It must form one constrained wall polynomial,
one wall/transport state, and one conservative wall flux for adjacent updates.
The merge gate remains the same fixed-time N=64/96/128/192 sequence. Until
`dTdn` and `tau_mag` pass it, the code must not be described as fully
second-order for coupled NS wall heat flux and skin friction.
