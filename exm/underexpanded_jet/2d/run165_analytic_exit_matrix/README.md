# Run165 RZ analytic-exit mechanism block

This case is the non-IBM half of the axis-mode root-cause matrix.  It is not a
replacement for the full Run165 forebody calculation and must not be compared
to the IBM block by absolute Mach-disk location or force.

The axial-high boundary contains an exact quasi-1D Run165 virtual-exit disk;
the rest of that boundary is extrapolating.  The external state, gas model,
LLF/WENO-Z5 operator, staged jet scaling, and RZ switches match the IBM case.

The base spacing is `r_exit/3`.  Consequently the L1/L2/L3 exit disks contain
exactly 6/12/24 complete annuli.  `prob.target_level=1,2,3` refines the complete
small domain uniformly, giving about 12k/49k/197k leaf cells and no AMR
interface.  Within each level, legacy and annular-only branches must restart
the same jet-off checkpoint and keep `rz_euler_point_flux=0`.

This block can test whether an analytic, face-based jet source creates the
same first-ring `2 dx` entropy mode.  Because the IBM and analytic blocks do
not share geometry or checkpoints, their cross-block comparison is a
mechanism fingerprint, not a formal nozzle main effect.

For the all-fluid RZ Skew--JST comparison, build with
`EULER_SCHEME=skew`.  The default Skew parameters are fourth order,
`C2=1.5`, `C4=0.016`, with the original first-power discontinuity sensor.
This changes only the Euler face-flux operator; geometry, boundary states,
jet ramp, AMR coverage, and Run165 thermodynamic states remain unchanged.

The legacy second-order characteristic MUSCL reconstruction with the HLLC
Riemann flux is selected with `EULER_SCHEME=hllc-muscl`.  It uses the same
RZ finite-volume divergence and the same case controls as the WENO and
Skew--JST variants.

The fifth-order point-state WENO-Z reconstruction, HLLC interface flux, and
alternative finite-difference correction are selected with
`EULER_SCHEME=afd-hllc-weno`.  The runtime AFD shock-to-LLF fallback remains
available through the standard `cns.afd_*` controls.

The face-level MUSCL/AFD diagnostics used by the pointwise failure report are
compiled only when `FACE_TRACE=1` is passed to `make`; normal and production
builds contain no active trace selection.  The target face is selected at run
time with `cns.afd_face_trace_*` or `cns.muscl_face_trace_*`.
