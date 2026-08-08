---
icon: laptop-arrow-down
cover: >-
  https://images.unsplash.com/photo-1577401239170-897942555fb3?crop=entropy&cs=srgb&fm=jpg&ixid=M3wxOTcwMjR8MHwxfHNlYXJjaHw2fHxwcm9ibGVtfGVufDB8fHx8MTczMDk3NjczN3ww&ixlib=rb-4.0.3&q=85
coverY: 0
---

# Input

This page provides detailed info on the `input` file options. These arguments are passed to AMReX, not all arguents are reuired, see examples

## AMReX options

The AMReX options covers control of

* the problem domain definition
* time-stepping
* gridding and load balancing
* output files
* checkpoint and restarting

The reference is available on [AMReX's documentation](https://amrex-codes.github.io/amrex/docs_html/Inputs_Chapter.html).

{% hint style="info" %}
In the tables below, _DIM_ means the number of dimensions, _Int_ means integer, and _Bool_ means boolean value (0 for False and 1 for True). If the option has no default value, a value must be given by the user.
{% endhint %}

### Problem definition and time-stepping

<table><thead><tr><th width="228">Option</th><th>Type</th><th align="center">Default</th><th>Description</th></tr></thead><tbody><tr><td><strong>max_step</strong></td><td>Int</td><td align="center"></td><td>Maximum number of time steps to take</td></tr><tr><td><strong>stop_time</strong></td><td>Real</td><td align="center"></td><td>Maximum time to reach</td></tr><tr><td><strong>time_step</strong></td><td>Real</td><td align="center"></td><td>dt (base level), higher level time step is based on number of subcycles</td></tr><tr><td><strong>cfl</strong></td><td>Real</td><td align="center"></td><td>CFL (incompatible option with time_step)</td></tr><tr><td><strong>geometry.is_periodic</strong></td><td>DIM * Int</td><td align="center">0 0 0</td><td>1 for true, 0 for false (one value for each coordinate direction)</td></tr><tr><td><strong>geometry.coord_sys</strong></td><td>Int</td><td align="center">0</td><td>0 = Cartesian; 1 = Cylindrical; 2 = Spherical (only support Cartesian)</td></tr><tr><td><strong>geometry.prob_lo</strong></td><td>DIM * Real</td><td align="center">0 0 0</td><td>Low corner of physical domain (physical not index space)</td></tr><tr><td><strong>geometry.prob_hi</strong></td><td>DIM * Real</td><td align="center"></td><td>High corner of physical domain (physical not index space)</td></tr><tr><td><strong>geometry.prob_extent</strong></td><td>DIM * Real</td><td align="center"></td><td>Extent of physical domain, choose between this or <code>prob_hi</code></td></tr><tr><td><strong>amr.n_cell</strong></td><td>DIM * Int</td><td align="center"></td><td>Number of cells at level 0 in each coordinate direction</td></tr></tbody></table>

For cylindrical R-Z calculations, `geometry.prob_lo[0]` must be exactly zero
because the lower radial face is the symmetry axis.

### Gridding and load balancing

<table><thead><tr><th width="252">Option</th><th width="160">Type</th><th width="134" align="center">Default</th><th>Description</th></tr></thead><tbody><tr><td><strong>amr.max_level</strong></td><td>Int</td><td align="center">0</td><td>Maximum level of refinement allowed<br>(0 when single-level)</td></tr><tr><td><strong>amr.ref_ratio</strong></td><td>Int *(nlev-1)</td><td align="center"></td><td>Refeniment ratio per level. If the number of <code>ref_ratio</code> is less than the number of levels - 1, the last entry will be automatically propagated</td></tr><tr><td><strong>amr.regrid_int</strong></td><td>Int</td><td align="center">-1</td><td>How often to regrid (in number of steps). No regridding will occur if set to &#x3C; 0</td></tr><tr><td><strong>amr.max_grid_size</strong></td><td>Int</td><td align="center">32</td><td>Maximum number of cells in each grid in all directions</td></tr><tr><td><strong>amr.blocking_factor</strong></td><td>Int</td><td align="center">8</td><td>Each grid must be divisible by blocking_factor in all directions (must be 1 or power of 2)</td></tr><tr><td><strong>amr.refine_grid_layou</strong>t</td><td>Bool</td><td align="center">1</td><td>Split grids in half until the number of grids is no less than the number of procs</td></tr><tr><td><strong>amr.n_error_buf</strong></td><td>DIM * Int</td><td align="center">1 1 1</td><td>Buffer in added around tagged cells</td></tr><tr><td><strong>amr.grid_eff</strong></td><td>Real</td><td align="center">0.7</td><td>Target value of the percentage of tagged cells in the grids</td></tr><tr><td><strong>amr.loadbalance_level0_int</strong></td><td>Int</td><td align="center">2</td><td>How often to do load balance (in number of steps). For single level (i.e., amr.max_level=0) only</td></tr><tr><td><strong>amr.loadbalance_with_workestimates</strong></td><td>Bool</td><td align="center">0</td><td>For multi-level runs, load balance is done during regrid and thus the load balance interval is controlled by <code>regrid_int</code></td></tr><tr><td><strong>amr.loadbalance_max_fac</strong></td><td>Real</td><td align="center">1.5</td><td>This controls the change in the maximum number of boxes that can be assigned to an MPI rank in load balancing</td></tr></tbody></table>

#### Conservative state prolongation

| Option | Type | Default | Description |
| --- | --- | ---: | --- |
| **cns.amr_state_interp** | String | `linear` | `linear` and `euler_admissible_linear` select the closure-admissibility-preserving limited conservative mapper. It first evaluates AMReX limited-linear interpolation, then applies one common `theta` to every conserved-variable deviation in a parent refinement block if any child violates positive density, the pressure-derived internal-energy-density floor, or the configured temperature-derived specific-internal-energy floor. `legacy_linear` and `conservative_quartic` are diagnostic smooth-flow choices and are rejected by pure-GP and EB builds. |

The common parent-block coefficient preserves the Cartesian child average and,
for R--Z, the annular-volume-weighted child average. It is not independent
cell clipping: an inadmissible accepted parent state fails closed instead of
being replaced. When `CLIP_MINTEMP=TRUE`, prolongation enforces the same
`c_v T_min` lower bound used by conservative-to-primitive conversion, so
regridding cannot manufacture a state that immediately requires hidden
temperature clipping. Startup prints the selected mapper and active floors.

### Outputs and Restarting

<table><thead><tr><th width="247">Option</th><th>Type</th><th width="133" align="center">Default</th><th>Description</th></tr></thead><tbody><tr><td><strong>amr.plot_files_output</strong></td><td>Bool</td><td align="center">1</td><td>Output plotfile or not (redundent because one can set plot_int = -1 to disable output)</td></tr><tr><td><strong>amr.plot_file</strong></td><td>String</td><td align="center">./plot/plt</td><td>Prefix of plotfile output</td></tr><tr><td><strong>amr.plot_int</strong></td><td>Int</td><td align="center">-1</td><td>Frequency of plotfile output; if -1 then no plotfiles will be written</td></tr><tr><td><strong>amr.derive_plot_vars</strong></td><td>Strings</td><td align="center">NONE</td><td>List of derived variables to plot; can use "ALL" or "NONE" to select all or none of the variables. See the full list of derived variables available in Cerisse below.</td></tr><tr><td><strong>amr.checkpoint_files_output</strong></td><td>Bool</td><td align="center">1</td><td>Same as plot_files_output, but for checkpoint files</td></tr><tr><td><strong>amr.check_file</strong></td><td>String</td><td align="center">chk</td><td>Prefix of checkpoint file output</td></tr><tr><td><strong>amr.check_int</strong></td><td>Int</td><td align="center">-1</td><td>Frequency of checkpoint file output; if -1 then no plotfiles will be written</td></tr><tr><td><strong>amr.restart</strong></td><td>String</td><td align="center"></td><td>If present, then the name of checkpoint file to restart from</td></tr><tr><td><strong>amr.plotfile_on_restar</strong>t</td><td>Bool</td><td align="center">0</td><td>Write a plotfile when immediately after restart or not</td></tr></tbody></table>

#### Restart recovery from a stop-time-clipped time step

| Option | Type | Default | Description |
| --- | --- | ---: | --- |
| **cns.restart_first_dt_from_cfl** | Bool | 0 | With CFL-based time stepping, arm a one-shot exception when a checkpoint is read. The first normal `computeNewDt` uses the current CFL estimate without applying the checkpoint time step's legacy 1.1 growth cap. |

This option is intended for explicitly segmented calculations whose final step
was shortened to land on `stop_time`. It must be `0` or `1`, and value `1`
requires `cfl` rather than `time_step`. The current CFL estimate, the absolute
`dt_max` cap, AMR subcycling constraints, and the new segment's `stop_time`
remain active. A post-regrid estimate is always capped by the pre-regrid time
step and does not consume the one-shot exception; ordinary later steps retain
the 1.1 growth limit. Fresh calculations do not arm the exception.

When enabled, startup reports the option, checkpoint reconstruction reports
that the exception is armed, and its single use reports the old coarse time
step, the CFL candidate, and the finally selected time step. Keep the default
for legacy restart trajectories, and record this opt-in in run provenance
because it changes the first continued time step.

### GPU-related parameters

<table><thead><tr><th width="252">Option</th><th>Type</th><th width="163" align="center">Default</th><th>Description</th></tr></thead><tbody><tr><td><strong>amrex.the_arena_init_size</strong></td><td>Int</td><td align="center">3/4 of total device memory</td><td>GPU device memory allocated to The_Arena (in bytes)</td></tr></tbody></table>

### Other parameters

<table><thead><tr><th width="244">Option</th><th>Type</th><th align="center">Default</th><th>Description</th></tr></thead><tbody><tr><td><strong>amrex.omp_threads</strong></td><td>String or Int</td><td align="center"><code>system</code></td><td><code>nosmt</code>: avoid using threads for virtual cores (aka Hyperthreading or SMT), as is default in OpenMP; <code>system</code>: use the environment variable <code>OMP_NUM_THREADS</code>. For Integer values, <code>OMP_NUM_THREADS</code> is ignored.</td></tr><tr><td><strong>amrex.fpe_trap_invalid</strong></td><td>Bool</td><td align="center">0</td><td>Produce error when invalid floating-point arthematic is detected. Helpful for debugging</td></tr></tbody></table>

### Boundary conditions

```ini
# 0 = Interior                               3 = Symmetry
# 1 = Inflow / UserBC                        4 = SlipWall =3
# 2 = Outflow (First Order Extrapolation)    5 = NoSlipWall (adiabatic)
cns.lo_bc = 1 5 0
cns.hi_bc = 2 5 0
```

In the above example, the bc in _x_ would be inflow (at lower boundary) and outflow (at upper boundary) while no specific boudnary will be defined in _z_

| Option         | Type       | Default | Description                           |
| -------------- | ---------- | ------- | ------------------------------------- |
| **cns.lo\_bc** | DIM \* Int |         | BC flags at lower boundaries in x,y,z |
| **cns.hi\_bc** | DIM \* Int |         | BC flags at upper boundaries in x,y,z |

If option "0" is selected, the corresponding `geometry.is_periodic` must also be set.

If option "1" is selected, the `bcnormal` function in `prob.H` will be activated (see boundary conditions).

### Characteristic LLF local fallback

The characteristic LLF driver provides one sensor-gated, face-local,
first-order fallback for LLF-WENO-Z5, LLF-TENO5, and LLF-TENO6. It is used in
Cartesian and R-Z coordinates by both all-fluid and pure shared-GP
calculations. The fallback replaces one shared background-face flux. It does
not introduce cut cells, open-face fractions, or wall-face fluxes.

The general sensor requires both an adjacent relative density jump and a
normalized density second difference. This prevents an ordinary smooth
gradient from activating the first-order flux on practical verification
grids. The first interior R-Z radial face retains the density-jump gate alone
because this is the independently verified near-axis contact safeguard. The
Cartesian, general R-Z, TENO, and shared-GP extensions use the common
construction but require the regression checks for each solver configuration.
The axis face at `r=0` is excluded because it stores an auxiliary
radius-weighted flux.

| Option | Type | Default | Description |
| --- | --- | ---: | --- |
| **cns.llf_first_order_fallback** | Bool | 1 | Enable the shared face-local first-order LLF fallback. |
| **cns.llf_fallback_density_jump** | Real | 0.003 | Lower adjacent relative-density-jump threshold. |
| **cns.llf_fallback_density_curvature** | Real | 0.08 | Lower normalized density second-difference threshold away from the first interior R-Z radial face. |

The R-Z radial fallback is formed directly in the radius-weighted flux. The
same final face value is used by both adjacent full background cells and by
AMR refluxing. This mechanism is a robustness safeguard. It is not a proof of
positivity or monotonicity.

### All-fluid AFD shock fallback

These controls apply to the Cartesian and R-Z all-fluid AFD-HLLC-WENO-Z5 and
AFD-HLLC-TENO5 operators. Pure shared-GP IBM and experimental EB drivers keep
the all-fluid HLLC-MUSCL and reconstructed-speed branches disabled. The broad
smoothness test still controls the AFD flux correction. The troubled-face
fallback requires all three
conditions: a nonsmooth stencil, a pressure jump, and local compression.
Compression is measured by the negative velocity divergence and is scaled
with the minimum cell spacing and the local sound speed. A density jump or
tangential-velocity curvature cannot activate this sensor-controlled
high-order LLF fallback by itself. The existing first-order cell-centred LLF
safeguards remain available for invalid stencils, rarefaction pockets,
inadmissible reconstructed states, and failed HLLC fluxes.

| Option | Type | Default | Description |
| --- | --- | ---: | --- |
| **cns.afd_correction** | Bool | 1 | Apply the AFD flux correction on stencils accepted by the broad smoothness test. |
| **cns.afd_shock_llf** | Bool | 1 | Allow the characteristic high-order LLF fallback on stencils that satisfy the joint sensor. When HLLC-MUSCL is enabled, LLF is used only if its MUSCL candidate is invalid. |
| **cns.afd_shock_hllc_muscl** | Bool | 1 for WENO-Z5; 0 for TENO5 | Replace a sensor-identified shared troubled face by the existing limited characteristic MUSCL reconstruction and HLLC flux. This is the default AFD-HLLC-WENO-Z5 shock fallback; TENO5 retains its previous default until separately validated. |
| **cns.afd_hllc_reconstructed_speed_ratio** | Real | 0 | Experimental a posteriori all-fluid fallback for finite but extreme point-WENO/HLLC states. If the reconstructed total signal speed exceeds this multiple of the six-cell physical envelope, replace that shared face with HLLC-MUSCL (and use two-point LLF only if the MUSCL candidate is invalid). Set to `0` to disable. This optional fallback catches the first observed extreme reconstruction but is not a positivity guarantee and is insufficient by itself for the Run165 long test. The pure shared-GP IBM driver explicitly leaves it disabled. |
| **cns.afd_smoothness_threshold** | Real | 0.08 | Upper normalized second-difference limit used by the broad density, pressure, and velocity smoothness test. |
| **cns.afd_shock_pressure_jump** | Real | 0.01 | Lower limit for the maximum adjacent relative pressure jump on the six-point face stencil. |
| **cns.afd_shock_compression** | Real | 0.001 | Lower limit for the maximum dimensionless negative velocity divergence on the four interior stencil points. |

The two shock thresholds must lie in `[0,1]`. The frozen pure shared-GP IBM
path retains its established marker-aware fallback and does not use these two
bulk sensor thresholds.

### Axisymmetric LLF--WENO pressure discretisation

The two-dimensional, all-fluid R--Z LLF--WENO-Z5 operator evaluates the
radial advective metric divergence and pressure gradient as a fixed matched
pair,

```text
-(1/r) d(r rho ur^2)/dr - dp/dr.
```

The pressure contribution reuses the nonlinear WENO weights obtained from the
complete metric LLF reconstruction. The complete Euler flux remains in the
standard face array, while the paired pressure flux is stage-local scratch
data. This discretisation is automatic and has no runtime switch. It does not
alter the GP-IBM, AMReX EB, AFD, TENO, or skew-symmetric operators.

Smooth single-level R--Z MMS tests retain fifth-order accuracy. On refined
R--Z hierarchies, the present coarse--fine accuracy is limited by the
R--Z-aware linear candidate inside the Euler-admissible parent-block mapper.
Its common limiter coefficient preserves AMReX's annular-volume-weighted
parent average. The Cartesian conservative-quartic interpolator is not
metric-aware in R--Z and must not be used to claim high-order R--Z AMR
accuracy.

### Characteristic outflow (experimental)

The Stage 3B Giles/GC-NSCBC path is opt-in and default-off. Its currently
verified scope is one subsonic outflow on a two-dimensional, Cartesian,
uniform level-0 grid with a periodic tangential direction, single-species
ideal-gas Euler, LLF-WENO-Z5, and SSPRK(4,3). Unsupported R--Z, refined-level,
multiple-boundary, physical-corner, and reverse-flow configurations terminate
instead of silently reverting to extrapolation. See
`docs/stage3_characteristic_boundary_20260720.md` for the equations, evidence,
and remaining gates.

| Option | Type | Default | Description |
| --- | --- | ---: | --- |
| **cns.nscbc_lo / cns.nscbc_hi** | DIM \* Int | all 0 | Per-side characteristic model: `0` disabled, `1` relaxed inflow, `2` outflow without mean-pressure relaxation, `3` pressure-relaxed outflow. Stage-local GC mode currently requires exactly one type `2` or `3` side and rejects type `1`. |
| **cns.nscbc_outflow_transverse_model** | Int | 0 | `0` retains the historical projected-NSCBC model; `1` selects the Giles second-order two-dimensional unsteady outflow closure. |
| **cns.nscbc_ghost_update_model** | Int | 0 | `0` retains the historical persistent ghost-state ODE; `1` constructs three stage-local ghost layers using GC-NSCBC. Mode `1` requires the Giles transverse model. |
| **cns.nscbc_order** | Int | 2 | One-sided normal derivative order (`1` or `2`). The certified Giles/GC configuration uses `2`. |
| **cns.nscbc_use_transverse** | Bool | 1 | Enable transverse derivatives. Giles mode requires `1`. |
| **cns.nscbc_Ptarget** | Real | 1 | Mean target pressure used only by type `3`; must be specified explicitly for that type. |
| **cns.nscbc_sigma** | Real | 0 | Non-negative mean-pressure relaxation coefficient used only by type `3`. |
| **cns.nscbc_Lchar** | Real | 1 | Positive characteristic relaxation length used by type `3`. |
| **cns.nscbc_Mmax** | Real | 0 | Reference Mach number for type-3 pressure relaxation; require `0 <= Mmax < 1`. |
| **cns.nscbc_min_rho / min_T / min_p** | Real | `1e-12` | Admissibility thresholds. Violation aborts the stage; these values are not clipping floors. |

### Immersed-boundary AMR support

The fixed-body ghost-point IBM reconstructs image points from same-level
primitive states. When AMR is enabled, the complete nonzero interpolation
support must therefore remain inside the fine level; silently substituting a
coarse-fine ghost value changes the configured WLS moment conditions and its
formal order.

| Option | Type | Default | Description |
| --- | --- | ---: | --- |
| **ib.amr_support_buffer** | Bool | 1 | Add an IBM-driven refinement band whose width is derived from the next level's volume and surface image-point halos. |
| **ib.amr_support_audit** | Bool | 1 | Audit every nonzero cached GP and surface interpolation weight against the global same-level valid `BoxArray`. Periodic indices are wrapped and physical-boundary supports are reported separately. |
| **ib.amr_support_strict** | Bool | 1 | Abort during IBM initialization, before the first RHS evaluation, if an interior support crosses a coarse-fine boundary. |

Keep all three options enabled in production. Setting
`ib.amr_support_buffer=0` is intended only for the negative validation fixture;
setting `ib.amr_support_strict=0` converts a correctness failure into a warning
and does not make the resulting reconstruction order-valid.

### IBM positivity controls

| Option | Type | Default | Description |
| --- | --- | ---: | --- |
| **cns.ibm_positivity_flux_limiter** | Bool | 0 | Opt-in pure-GP shared-background-face convex limiter for every Forward-Euler bracket of SSPRK(4,3). The established scope is fixed Euler-slip GP-IBM, single-species ideal gas, source-free LLF-WENO-Z5, full background cells/faces, and non-periodic boundaries. Cartesian stages blend the complete flux. The R-Z implementation candidate uses the paired operator and applies one shared face `theta` to both the complete flux and its radial pressure companion. Its uniform, forced-limiter, and multi-FAB integration smokes are qualified only on a single level (`amr.max_level=0`); multilevel annular GP and reflux remain fail-closed and unqualified, so this is not yet a full R-Z production qualification. |
| **cns.ibm_positivity_flux_limiter_verbose** | Bool | 0 | Print the high-order, limited-flux, and post-SSP admissibility verdict for every SSPRK(4,3) bracket. |
| **cns.ibm_positivity_retry** | Bool | 0 | Reserved. The pure full-cell transaction path currently terminates if this option is enabled. |
| **cns.strict_positivity** | Bool | 0 | Strengthen the accepted-state audit from `rho>0` and `rho e>0` to margins above the configured density and specific-internal-energy floors. This is an aborting check. |
| **cns.soft_positivity** | Bool | 0 | Replace a finite inadmissible active-fluid cell with a floor state. This changes the conservative solution and is forbidden when the conservative flux limiter is enabled. It is not an accepted verification or production setting. |

The limiter constructs high- and piecewise-constant Rusanov low-order fluxes
on the same full Cartesian faces. Each internal face has one final shared
limited flux. Local limiting falls back to a global flux blend and then the
all-low-order update; failure of the all-low-order update terminates without
cell-wise clipping. Problem/physics source terms, active-cell RHS overrides,
AMR/reflux, moving bodies, Navier--Stokes and chemistry remain outside the
certified limiter scope.  The coordinate-required R--Z metric/pressure pair is
part of the flux operator, not an optional problem source.

The R--Z radial pressure companion is stage-local auxiliary data in the
paired metric-advection/ordinary-pressure-gradient operator.  The same
limiting `theta` is applied to the complete face flux and this companion, but
the companion is not registered as an independent conservative AMR flux.
Current annular shared-GP R--Z runs fail closed when `amr.max_level>0`, even
if reflux is disabled.  Multilevel reflux/average-down does not yet provide an
invariant-domain-preserving transaction and is not qualified by a post-sync
admissibility audit alone.

Enabling the limiter also forbids problem or physics source terms and any
problem hook (such as an active-cell `rhs_nscbc` hook) that may overwrite the
assembled active-cell RHS.  Such a hook terminates rather than bypassing the
high/low shared-face transaction.  Only the separately guarded stage-local
ghost-cell NSCBC coupling can pass this hook gate because it does not modify
the active-cell RHS; its own coordinate and boundary qualification gates still
apply.

For pure GP-IBM, physical wall pressure and integrated pressure force are
obtained from BI surface recovery and surface quadrature. Cartesian
fluid--solid stencil-crossing momentum exchange is not a pointwise physical
wall traction and is not projected back onto the BI pressure.

### Geometry EB options

In the **input** file, users should specify the geometry of the embedded boundary with `eb2.geom_type`, then supply the required parameters in the format of `eb2.{geom_param}`.

| eb2.geom\_type      | additional parameters required                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| ------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `all_regular`       | no EB, no additional parameters needed                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| `plane`             | `plane_point` - a point where the plane intersects, `plane_normal` - the normal vector of the plane that points into the solid                                                                                                                                                                                                                                                                                                                                                                             |
| `sphere`            | `sphere_center`, `sphere_radius`, `sphere_has_fluid_inside` - bool value fluid inside or outside                                                                                                                                                                                                                                                                                                                                                                                                           |
| `cylinder`          | `cylinder_center`, `cylinder_radius`, `cylinder_height`, `cylinder_direction` - (0,1,2) for (x,y,z), and `cylinder_has_fluid_inside`                                                                                                                                                                                                                                                                                                                                                                       |
| `box`               | `box_lo` and `box_hi` - lower and upper corners of the box, and `box_has_fluid_inside`                                                                                                                                                                                                                                                                                                                                                                                                                     |
| `stl`               | `stl_file` - the STL file name, `stl_scale` - the scaling factor in all directions, and `stl_center` - center of object in relation to the cooridinate system in the file, and `stl_reverse_normal` - essentially stl\_has\_fluid\_inside                                                                                                                                                                                                                                                                  |
| `triangles`         | <p><code>num_tri</code> - number of triangles, up to 5 (change the value in <code>custom_geometry.cpp</code> if needed),<br>for each triangle, <code>{i}</code> from 0 to num_tri-1, <code>tri_{i}_point_0</code>, <code>tri_{i}_point_1</code>, <code>tri_{i}_point_2</code> - three points that define the triangle, give the points in anti-clockwise direction to set solid inside of the triangle. The z-coordinate isn't really needed because the triangle will be extruded in the z-direction.</p> |
| `combustor`         | `far_wall_loc`, `ramp_plane1_point`, `ramp_plane2_point`, `ramp_plane2_normal`, `ramp_plane3_point`, `pipe_lo`, `pipe_hi`                                                                                                                                                                                                                                                                                                                                                                                  |
| `converging-nozzle` | `d_inlet`, `l_inlet`, `d_exit`, `l_nozzle`                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |

NOTE: You can only choose one geometry type and one geometry. If you want to use multiple geometries, you need to define your own geometry

Below are some examples:

```ini
eb2.geom_type = all_regular

eb2.geom_type = cylinder
eb2.cylinder_direction = 2
eb2.cylinder_radius = 0.25
eb2.cylinder_center = 1.0 2.0 0.0
eb2.cylinder_has_fluid_inside = 0

eb2.geom_type = box
eb2.box_lo = -1.0 -1.0 0.0
eb2.box_hi =  1.0  2.0  0.0
eb2.box_has_fluid_inside = 0

# This gives the same geometry as the box above
eb2.geom_type = triangles 
triangles.num_tri = 2
triangles.tri_0_point_0 = -1.0  1.0 0.0
triangles.tri_0_point_1 = -1.0 -1.0 0.0
triangles.tri_0_point_2 =  1.0 -1.0 0.0
triangles.tri_1_point_0 = -1.0  1.0 0.0
triangles.tri_1_point_1 =  1.0 -1.0 0.0
triangles.tri_1_point_2 =  1.0  1.0 0.0
```
